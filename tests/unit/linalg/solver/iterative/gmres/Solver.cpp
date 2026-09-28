#include <gtest/gtest.h>

#include "linalg/types/Vector.hpp"
#include "linalg/solver/iterative/gmres/Solver.hpp"
#include "linalg/operator/CSROperator.hpp"
#include "linalg/solver/preconditioner/Identity.hpp"
#include "utils/logging/core/NullLogger.hpp"
#include "utils/logging/linear/ConsoleLogger.hpp"

using namespace residuum;

TEST(GMRESSolver, SolveNonsymmetric2x2) {

	// declare types
	using Backend = linalg::types::backend::CPU;
	using Vec = linalg::types::Vector<Real, Backend>;
	using Mat = linalg::types::CSRMatrix<Real, Backend>;

	// build csr matrix
	// A = [[4,1], [3,3]] -- nonsymmetric, CG would not be valid here
	Mat A(2,2);
	A.resize(4);

	A.rowPtr()[0] = 0;
	A.rowPtr()[1] = 2;
	A.rowPtr()[2] = 4;

	A.colIdx()[0] = 0; A.data()[0] = 4.0;
	A.colIdx()[1] = 1; A.data()[1] = 1.0;
	A.colIdx()[2] = 0; A.data()[2] = 3.0;
	A.colIdx()[3] = 1; A.data()[3] = 3.0;

	linalg::op::CSROperator<Mat> op(A);

	// setup linear system Ax = b, with b = [[1], [2]]
	Vec b(2);
	Vec x(2);

	b.data()[0] = 1.0;
	b.data()[1] = 2.0;

	x.zero();

	// setup solver workspace & report
	linalg::solver::iterative::gmres::Workspace<Vec> W(2, 5);
	linalg::solver::SolverReport<Vec> report;

	// setup preconditioner & logger
	linalg::solver::preconditioner::Identity<Vec> M;
	utils::logging::linear::ConsoleLogger logger("Linear System", "GMRES", "Identity", {"x"});

	// setup solver config
	const Real solverTol = 1e-12;
	linalg::solver::iterative::gmres::Config<Vec> cfg{solverTol};
	cfg.krylovDim = 5;

	// declare solver
	linalg::solver::iterative::gmres::Solver<decltype(op), Vec, decltype(M), decltype(logger)> solver(cfg);

	// solve system
	bool converged = solver.solve(report, logger, W, M, op, b, x);

	// check solution x = [[1/9], [5/9]] -- A^-1 b for A = [[4,1],[3,3]], b = [1,2]
	const Real tol = 1e-10;
	EXPECT_TRUE(converged);
	EXPECT_NEAR(x.data()[0], 1.0/9.0, tol);
	EXPECT_NEAR(x.data()[1], 5.0/9.0, tol);

}

TEST(GMRESSolver, SolveSPD2x2MatchesKnownSolution) {

	// GMRES should also work fine on a symmetric system, just without exploiting the symmetry
	using Backend = linalg::types::backend::CPU;
	using Vec = linalg::types::Vector<Real, Backend>;
	using Mat = linalg::types::CSRMatrix<Real, Backend>;

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

	Vec b(2);
	Vec x(2);

	b.data()[0] = 1.0;
	b.data()[1] = 2.0;

	x.zero();

	linalg::solver::iterative::gmres::Workspace<Vec> W(2, 5);
	linalg::solver::SolverReport<Vec> report;

	linalg::solver::preconditioner::Identity<Vec> M;
	utils::logging::NullLogger logger;

	const Real solverTol = 1e-12;
	linalg::solver::iterative::gmres::Config<Vec> cfg{solverTol};
	cfg.krylovDim = 5;

	linalg::solver::iterative::gmres::Solver<decltype(op), Vec, decltype(M), decltype(logger)> solver(cfg);

	bool converged = solver.solve(report, logger, W, M, op, b, x);

	const Real tol = 1e-10;
	EXPECT_TRUE(converged);
	EXPECT_NEAR(x.data()[0], 1.0/11.0, tol);
	EXPECT_NEAR(x.data()[1], 7.0/11.0, tol);

}

TEST(GMRESSolver, RestartsWhenKrylovDimSmallerThanSystem) {

	// krylovDim = 1 forces at least one restart cycle to reach the tolerance on a 2x2 system
	using Backend = linalg::types::backend::CPU;
	using Vec = linalg::types::Vector<Real, Backend>;
	using Mat = linalg::types::CSRMatrix<Real, Backend>;

	Mat A(2,2);
	A.resize(4);

	A.rowPtr()[0] = 0;
	A.rowPtr()[1] = 2;
	A.rowPtr()[2] = 4;

	A.colIdx()[0] = 0; A.data()[0] = 4.0;
	A.colIdx()[1] = 1; A.data()[1] = 1.0;
	A.colIdx()[2] = 0; A.data()[2] = 3.0;
	A.colIdx()[3] = 1; A.data()[3] = 3.0;

	linalg::op::CSROperator<Mat> op(A);

	Vec b(2);
	Vec x(2);

	b.data()[0] = 1.0;
	b.data()[1] = 2.0;

	x.zero();

	linalg::solver::iterative::gmres::Workspace<Vec> W(2, 1);
	linalg::solver::SolverReport<Vec> report;

	linalg::solver::preconditioner::Identity<Vec> M;
	utils::logging::NullLogger logger;

	const Real solverTol = 1e-10;
	linalg::solver::iterative::gmres::Config<Vec> cfg{solverTol};
	cfg.krylovDim = 1;
	cfg.maxIters = 1000;

	linalg::solver::iterative::gmres::Solver<decltype(op), Vec, decltype(M), decltype(logger)> solver(cfg);

	bool converged = solver.solve(report, logger, W, M, op, b, x);

	const Real tol = 1e-8;
	EXPECT_TRUE(converged);
	EXPECT_GT(report.iterations, 1); // must have taken more than one restart cycle's worth of steps
	EXPECT_NEAR(x.data()[0], 1.0/9.0, tol);
	EXPECT_NEAR(x.data()[1], 5.0/9.0, tol);

}
