#include <gtest/gtest.h>

#include "linalg/types/Vector.hpp"
#include "linalg/solver/iterative/cg/Solver.hpp"
#include "linalg/operator/CSROperator.hpp"
#include "linalg/solver/preconditioner/Identity.hpp"
#include "utils/logging/core/NullLogger.hpp"
#include "utils/logging/linear/ConsoleLogger.hpp"

using namespace residuum;

TEST(CGSolver, Solve2x2SPD) {

	// declare types
	using Backend = linalg::types::backend::CPU;
	using Vec = linalg::types::Vector<Real, Backend>;
	using Mat = linalg::types::CSRMatrix<Real, Backend>;

	// build csr matrix
	// A = [[4,1], [1,3]]
	Mat A(2,2);
	A.resize(4);

	// fill A
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

	// setup solver workspace & report
	linalg::solver::iterative::cg::Workspace<Vec> W(2);
	linalg::solver::SolverReport<Vec> report;

	// setup preconditioner & logger
	linalg::solver::preconditioner::Identity<Vec> M;
	utils::logging::linear::ConsoleLogger logger("Linear System", "CG", "Identity", {"x"});

	// setup solver config
	const Real solverTol = 1e-12;
	linalg::solver::iterative::cg::Config<Vec> cfg{solverTol};

	// declare solver
	linalg::solver::iterative::cg::Solver<decltype(op), Vec, decltype(M), decltype(logger)> solver(cfg);

	// solver system
	solver.solve(report, logger, W, M, op, b, x);

	// check solution x = [[1/11], [7/11]]
	const Real tol = 1e-10;
	EXPECT_NEAR(x.data()[0], 1.0/11.0, tol);
	EXPECT_NEAR(x.data()[1], 7.0/11.0, tol);

}
