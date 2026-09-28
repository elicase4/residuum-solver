#include <algorithm>
#include <vector>

#include <gtest/gtest.h>

#include "core/Types.hpp"

#include "fem/dof/DOFOrdering.hpp"

#include "linalg/solver/base/SolverReport.hpp"
#include "linalg/types/Vector.hpp"
#include "linalg/types/backend/CPU.hpp"

#include "solver/config/NonlinearSolverConfig.hpp"
#include "solver/nonlinear/Newton.hpp"
#include "solver/stage/NonlinearCapableStage.hpp"

#include "utils/logging/core/NullLogger.hpp"
#include "utils/logging/nonlinear/Logger.hpp"

using namespace residuum;

namespace {

	// Minimal NonlinearCapableStage: not backed by any real FEM problem -- heat conduction has
	// no NonlinearTangentForm-conforming physics term yet (conductivity is always linear in T),
	// so there's nothing real to assemble a Newton step against (see #42). This mock exercises
	// NewtonRunner's OWN convergence-loop mechanics (iteration counting, tolerance checks, the
	// solveLinearStep-fails path) in isolation from that blocker.
	struct MockNonlinearStage {

		std::vector<Real> residualPerAssemble; // residualPerAssemble[k] is returned after the (k+1)-th assemble() call
		bool failSolveLinearStep = false;

		Index assembleCalls = 0;
		Index solveLinearStepCalls = 0;

		void initialize() {}
		void assemble() { ++assembleCalls; }
		bool solve() { return true; } // unused by NewtonRunner -- part of the Stage concept only
		void finalize() {}

		Real residualNorm() const {
			const Index idx = std::min<Index>(assembleCalls - 1, static_cast<Index>(residualPerAssemble.size()) - 1);
			return residualPerAssemble[idx];
		}

		// single-DOF mock, so the per-DOF split is just the aggregate residual
		linalg::types::Vector<Real, linalg::types::backend::CPU> residual() const {
			linalg::types::Vector<Real, linalg::types::backend::CPU> r(1);
			r.data()[0] = residualNorm();
			return r;
		}

		bool solveLinearStep() {
			++solveLinearStepCalls;
			return !failSolveLinearStep;
		}

	}; // struct MockNonlinearStage

	using VectorType = linalg::types::Vector<Real, linalg::types::backend::CPU>;
	using NewtonRunnerT = solver::nonlinear::NewtonRunner<MockNonlinearStage, VectorType>;

	solver::config::NonlinearSolverConfig makeConfig(Index maxIterations = 20) {
		solver::config::NonlinearSolverConfig cfg;
		cfg.absoluteTolerance = 1e-10;
		cfg.relativeTolerance = 1e-8;
		cfg.maxIterations = maxIterations;
		return cfg;
	}

	utils::logging::nonlinear::Logger makeQuietLogger() {
		return utils::logging::nonlinear::Logger(utils::logging::NullLogger{});
	}

} // namespace

static_assert(solver::stage::NonlinearCapableStage<MockNonlinearStage>);

TEST(NewtonRunnerTest, ConvergesWhenResidualDropsBelowAbsoluteTolerance) {

	MockNonlinearStage stage;
	stage.residualPerAssemble = {1.0, 1e-2, 1e-11}; // iter2's residual is below absoluteTolerance

	NewtonRunnerT runner(stage, makeConfig(), makeQuietLogger(), {"x"}, 1, fem::dof::DOFOrdering::Interleaved);

	linalg::solver::SolverReport<VectorType> report;
	const bool converged = runner.solve(report);

	EXPECT_TRUE(converged);
	EXPECT_TRUE(report.converged);
	EXPECT_EQ(report.iterations, 2);
	EXPECT_DOUBLE_EQ(report.finalResidual, 1e-11);

	EXPECT_EQ(stage.assembleCalls, 3); // once per iteration, including the converged one
	EXPECT_EQ(stage.solveLinearStepCalls, 2); // NOT called on the iteration that converges

}

TEST(NewtonRunnerTest, StopsImmediatelyWhenSolveLinearStepFails) {

	MockNonlinearStage stage;
	stage.residualPerAssemble = {1.0, 1.0, 1.0}; // never converges on its own
	stage.failSolveLinearStep = true;

	NewtonRunnerT runner(stage, makeConfig(), makeQuietLogger(), {"x"}, 1, fem::dof::DOFOrdering::Interleaved);

	linalg::solver::SolverReport<VectorType> report;
	const bool converged = runner.solve(report);

	EXPECT_FALSE(converged);
	EXPECT_FALSE(report.converged);
	EXPECT_EQ(stage.assembleCalls, 1);
	EXPECT_EQ(stage.solveLinearStepCalls, 1); // failed on the very first attempt, loop exits immediately

}

TEST(NewtonRunnerTest, ReportsNotConvergedAfterMaxIterations) {

	const Index maxIterations = 4;

	MockNonlinearStage stage;
	stage.residualPerAssemble.assign(maxIterations, 0.5); // stays flat -- relative residual never drops

	NewtonRunnerT runner(stage, makeConfig(maxIterations), makeQuietLogger(), {"x"}, 1, fem::dof::DOFOrdering::Interleaved);

	linalg::solver::SolverReport<VectorType> report;
	const bool converged = runner.solve(report);

	EXPECT_FALSE(converged);
	EXPECT_FALSE(report.converged);
	EXPECT_EQ(report.iterations, maxIterations);

	EXPECT_EQ(stage.assembleCalls, maxIterations);
	EXPECT_EQ(stage.solveLinearStepCalls, maxIterations); // every iteration attempted a correction

}

TEST(NewtonRunnerTest, FirstIterationRelativeResidualIsAlwaysOneRegardlessOfAbsoluteScale) {

	// a residual that starts far above absoluteTolerance must NOT be reported converged on
	// iteration 0 just because relative residual is trivially 1.0 there -- only absoluteTolerance
	// can converge on the very first assemble.
	MockNonlinearStage stage;
	stage.residualPerAssemble = {5.0, 5.0}; // large, constant -- would only "converge" via a
	                                         // buggy relative check that misreads iter 0

	NewtonRunnerT runner(stage, makeConfig(1), makeQuietLogger(), {"x"}, 1, fem::dof::DOFOrdering::Interleaved);

	linalg::solver::SolverReport<VectorType> report;
	const bool converged = runner.solve(report);

	EXPECT_FALSE(converged);

}
