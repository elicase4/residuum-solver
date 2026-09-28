#include <gtest/gtest.h>

#include <type_traits>

#include "linalg/types/Vector.hpp"
#include "linalg/types/CSRMatrix.hpp"
#include "linalg/operator/CSROperator.hpp"
#include "linalg/solver/iterative/cg/CGRunner.hpp"
#include "linalg/solver/preconditioner/Identity.hpp"
#include "utils/logging/core/NullLogger.hpp"

using namespace residuum;

TEST(CGRunner, Solve2x2SPD) {

	// declare types
	using Backend = linalg::types::backend::CPU;
	using Vec = linalg::types::Vector<Real, Backend>;
	using Mat = linalg::types::CSRMatrix<Real, Backend>;

	// build csr matrix
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

	linalg::op::CSROperator<Mat> op(A);

	// setup linear system Ax = b, with b = [[1], [2]]
	Vec b(2);
	Vec x(2);

	b.data()[0] = 1.0;
	b.data()[1] = 2.0;

	x.zero();

	// declare the runner directly (linalg-level, config-agnostic -- takes cg::Config, not
	// solver::config::LinearSolverConfig; the latter is exercised in
	// tests/unit/solver/linear/LinearSolverFactory.cpp instead)
	using Precond = linalg::solver::preconditioner::Identity<Vec>;
	using Logger = utils::logging::NullLogger;

	linalg::solver::iterative::cg::Config<Vec> cfg{1e-12};

	linalg::solver::iterative::cg::CGRunner<decltype(op), Vec, Precond, Logger> runner(op, op.size(), cfg, Logger{});

	linalg::solver::SolverReport<Vec> report;
	bool converged = runner.solve(b, x, report);

	// check solution x = [[1/11], [7/11]]
	const Real tol = 1e-10;
	EXPECT_TRUE(converged);
	EXPECT_NEAR(x.data()[0], 1.0/11.0, tol);
	EXPECT_NEAR(x.data()[1], 7.0/11.0, tol);

}

TEST(CGRunner, SatisfiesLinearSolverRunnerInterface) {

	using Backend = linalg::types::backend::CPU;
	using Vec = linalg::types::Vector<Real, Backend>;
	using Mat = linalg::types::CSRMatrix<Real, Backend>;
	using Precond = linalg::solver::preconditioner::Identity<Vec>;
	using Logger = utils::logging::NullLogger;
	using Op = linalg::op::CSROperator<Mat>;

	using RunnerT = linalg::solver::iterative::cg::CGRunner<Op, Vec, Precond, Logger>;

	// CGRunner should be usable polymorphically through the base interface, since that's the
	// whole point of the *Runner pattern -- callers (e.g. HeatProblem) only ever hold a
	// linalg::solver::LinearSolverRunner<Vec>*.
	static_assert(std::is_base_of_v<linalg::solver::LinearSolverRunner<Vec>, RunnerT>);

}
