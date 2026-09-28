#include <gtest/gtest.h>

#include "solver/linear/LinearSolverFactory.hpp"

#include "fem/dof/DOFOrdering.hpp"

#include "linalg/types/Vector.hpp"
#include "linalg/types/CSRMatrix.hpp"
#include "linalg/operator/CSROperator.hpp"

using namespace residuum;

namespace {

	using Backend = linalg::types::backend::CPU;
	using Vec = linalg::types::Vector<Real, Backend>;
	using Mat = linalg::types::CSRMatrix<Real, Backend>;

	Mat make2x2SPDMatrix() {

		// A = [[4,1], [1,3]]
		Mat A(2,2);
		A.resize(4);

		A.rowPtr()[0] = 0;
		A.rowPtr()[1] = 2;
		A.rowPtr()[2] = 4;

		A.colIdx()[0] = 0; A.data()[0] = 4.0;
		A.colIdx()[1] = 1; A.data()[1] = 1.0;
		A.colIdx()[2] = 0; A.data()[2] = 1.0;
		A.colIdx()[3] = 1; A.data()[3] = 3.0;

		return A;

	}

	Mat make2x2NonsymmetricMatrix() {

		// A = [[4,1], [3,3]] -- nonsymmetric, exercises GMRES specifically (CG requires SPD)
		Mat A(2,2);
		A.resize(4);

		A.rowPtr()[0] = 0;
		A.rowPtr()[1] = 2;
		A.rowPtr()[2] = 4;

		A.colIdx()[0] = 0; A.data()[0] = 4.0;
		A.colIdx()[1] = 1; A.data()[1] = 1.0;
		A.colIdx()[2] = 0; A.data()[2] = 3.0;
		A.colIdx()[3] = 1; A.data()[3] = 3.0;

		return A;

	}

} // namespace

TEST(LinearSolverFactory, CGConfigProducesConvergedSolution) {

	Mat A = make2x2SPDMatrix();
	linalg::op::CSROperator<Mat> op(A);

	Vec b(2), x(2);
	b.data()[0] = 1.0;
	b.data()[1] = 2.0;
	x.zero();

	solver::config::LinearSolverConfig cfg;
	cfg.type = solver::config::LinearSolverConfig::Type::CG;
	cfg.tolerance = 1e-12;
	cfg.maxIterations = 1000;

	solver::config::LinearLoggerConfig loggerCfg;
	loggerCfg.type = solver::config::LoggerConfig::Type::None;

	auto runner = solver::linear::makeLinearSolverRunner<decltype(op), Vec>(op, op.size(), cfg, loggerCfg, "Test", {"x"}, op.size(), fem::dof::DOFOrdering::Interleaved);

	linalg::solver::SolverReport<Vec> report;
	bool converged = runner->solve(b, x, report);

	const Real tol = 1e-9;
	EXPECT_TRUE(converged);
	EXPECT_NEAR(x.data()[0], 1.0/11.0, tol);
	EXPECT_NEAR(x.data()[1], 7.0/11.0, tol);

}

TEST(LinearSolverFactory, GMRESConfigProducesConvergedSolutionOnNonsymmetricSystem) {

	Mat A = make2x2NonsymmetricMatrix();
	linalg::op::CSROperator<Mat> op(A);

	Vec b(2), x(2);
	b.data()[0] = 1.0;
	b.data()[1] = 2.0;
	x.zero();

	solver::config::LinearSolverConfig cfg;
	cfg.type = solver::config::LinearSolverConfig::Type::GMRES;
	cfg.tolerance = 1e-12;
	cfg.maxIterations = 1000;
	cfg.params = solver::config::GMRESParams{5};

	solver::config::LinearLoggerConfig loggerCfg;
	loggerCfg.type = solver::config::LoggerConfig::Type::None;

	auto runner = solver::linear::makeLinearSolverRunner<decltype(op), Vec>(op, op.size(), cfg, loggerCfg, "Test", {"x"}, op.size(), fem::dof::DOFOrdering::Interleaved);

	linalg::solver::SolverReport<Vec> report;
	bool converged = runner->solve(b, x, report);

	const Real tol = 1e-9;
	EXPECT_TRUE(converged);
	EXPECT_NEAR(x.data()[0], 1.0/9.0, tol);
	EXPECT_NEAR(x.data()[1], 5.0/9.0, tol);

}

TEST(LinearSolverFactory, UnimplementedSolverTypesThrow) {

	Mat A = make2x2SPDMatrix();
	linalg::op::CSROperator<Mat> op(A);

	for (auto type : {solver::config::LinearSolverConfig::Type::BiCGSTAB, solver::config::LinearSolverConfig::Type::LU}) {

		solver::config::LinearSolverConfig cfg;
		cfg.type = type;

		solver::config::LinearLoggerConfig loggerCfg;
		loggerCfg.type = solver::config::LoggerConfig::Type::None;

		EXPECT_THROW((solver::linear::makeLinearSolverRunner<decltype(op), Vec>(op, op.size(), cfg, loggerCfg, "Test", {"x"}, op.size(), fem::dof::DOFOrdering::Interleaved)), std::runtime_error);

	}

}

TEST(LinearSolverFactory, NonIdentityPreconditionerThrows) {

	Mat A = make2x2SPDMatrix();
	linalg::op::CSROperator<Mat> op(A);

	solver::config::LinearSolverConfig cfg;
	cfg.type = solver::config::LinearSolverConfig::Type::CG;
	cfg.preconditioner.type = static_cast<solver::config::PreconditionerConfig::Type>(-1); // not Identity

	solver::config::LinearLoggerConfig loggerCfg;
	loggerCfg.type = solver::config::LoggerConfig::Type::None;

	EXPECT_THROW((solver::linear::makeLinearSolverRunner<decltype(op), Vec>(op, op.size(), cfg, loggerCfg, "Test", {"x"}, op.size(), fem::dof::DOFOrdering::Interleaved)), std::runtime_error);

}
