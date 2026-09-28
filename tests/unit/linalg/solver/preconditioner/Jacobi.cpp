#include <cmath>

#include <gtest/gtest.h>

#include "linalg/operator/CSROperator.hpp"
#include "linalg/solver/iterative/cg/Solver.hpp"
#include "linalg/solver/preconditioner/Identity.hpp"
#include "linalg/solver/preconditioner/Jacobi.hpp"
#include "linalg/types/CSRMatrix.hpp"
#include "linalg/types/Vector.hpp"
#include "utils/logging/core/NullLogger.hpp"

using namespace residuum;

namespace {

	using Backend = linalg::types::backend::CPU;
	using Vec = linalg::types::Vector<Real, Backend>;
	using Mat = linalg::types::CSRMatrix<Real, Backend>;

	// A = [[4,1,0], [1,5,1], [0,1,6]] -- symmetric tridiagonal, diagonal entries 4,5,6
	Mat make3x3Tridiagonal() {

		Mat A(3,3);
		A.resize(7);

		A.rowPtr()[0] = 0;
		A.rowPtr()[1] = 2;
		A.rowPtr()[2] = 5;
		A.rowPtr()[3] = 7;

		A.colIdx()[0] = 0; A.data()[0] = 4.0;
		A.colIdx()[1] = 1; A.data()[1] = 1.0;
		A.colIdx()[2] = 0; A.data()[2] = 1.0;
		A.colIdx()[3] = 1; A.data()[3] = 5.0;
		A.colIdx()[4] = 2; A.data()[4] = 1.0;
		A.colIdx()[5] = 1; A.data()[5] = 1.0;
		A.colIdx()[6] = 2; A.data()[6] = 6.0;

		return A;

	}

	// symmetric tridiagonal SPD matrix whose diagonal spans several orders of magnitude
	// (diag[i] = 10^(6*i/(n-1)), weak constant off-diagonal coupling), so Jacobi
	// preconditioning collapses its condition number while CG on the raw system struggles
	Mat makeBadlyScaledTridiagonal(Index n) {

		Mat A(n, n);
		A.resize(3*n - 2);

		const Real offDiag = 0.05;

		Index* rowPtr = A.rowPtr();
		Index* colIdx = A.colIdx();
		Real* data = A.data();

		Index p = 0;
		for (Index i = 0; i < n; ++i) {

			rowPtr[i] = p;

			const Real diagVal = std::pow(10.0, 6.0 * static_cast<Real>(i) / static_cast<Real>(n - 1));

			if (i > 0) { colIdx[p] = i - 1; data[p] = offDiag; ++p; }
			colIdx[p] = i; data[p] = diagVal; ++p;
			if (i + 1 < n) { colIdx[p] = i + 1; data[p] = offDiag; ++p; }

		}
		rowPtr[n] = p;

		return A;

	}

} // namespace

TEST(JacobiPreconditioner, UpdateComputesInverseDiagonalFromCSROperator) {

	Mat A = make3x3Tridiagonal();
	linalg::op::CSROperator<Mat> op(A);

	linalg::solver::preconditioner::Jacobi<Vec> M;
	M.update(op);

	Vec ones(3);
	ones.data()[0] = 1.0;
	ones.data()[1] = 1.0;
	ones.data()[2] = 1.0;

	Vec z(3);
	M.apply(ones, z);

	const Real tol = 1e-12;
	EXPECT_NEAR(z.data()[0], 1.0/4.0, tol);
	EXPECT_NEAR(z.data()[1], 1.0/5.0, tol);
	EXPECT_NEAR(z.data()[2], 1.0/6.0, tol);

}

TEST(JacobiPreconditioner, ReducesCGIterationsOnBadlyScaledSystem) {

	const Index n = 50;
	Mat A = makeBadlyScaledTridiagonal(n);
	linalg::op::CSROperator<Mat> op(A);

	Vec b(n);
	for (Index i = 0; i < n; ++i) b.data()[i] = 1.0;

	const Real solverTol = 1e-8;
	const Index maxIters = 5000;
	linalg::solver::iterative::cg::Config<Vec> cfg{solverTol};
	cfg.maxIters = maxIters;

	utils::logging::NullLogger logger;

	// solve with Identity preconditioning
	Vec xIdentity(n);
	xIdentity.zero();
	linalg::solver::iterative::cg::Workspace<Vec> WIdentity(n);
	linalg::solver::SolverReport<Vec> reportIdentity;
	linalg::solver::preconditioner::Identity<Vec> MIdentity;
	linalg::solver::iterative::cg::Solver<decltype(op), Vec, decltype(MIdentity), decltype(logger)> solverIdentity(cfg);
	solverIdentity.solve(reportIdentity, logger, WIdentity, MIdentity, op, b, xIdentity);

	// solve with Jacobi preconditioning
	Vec xJacobi(n);
	xJacobi.zero();
	linalg::solver::iterative::cg::Workspace<Vec> WJacobi(n);
	linalg::solver::SolverReport<Vec> reportJacobi;
	linalg::solver::preconditioner::Jacobi<Vec> MJacobi;
	linalg::solver::iterative::cg::Solver<decltype(op), Vec, decltype(MJacobi), decltype(logger)> solverJacobi(cfg);
	solverJacobi.solve(reportJacobi, logger, WJacobi, MJacobi, op, b, xJacobi);

	EXPECT_TRUE(reportJacobi.converged);
	EXPECT_LT(reportJacobi.iterations, reportIdentity.iterations);

}
