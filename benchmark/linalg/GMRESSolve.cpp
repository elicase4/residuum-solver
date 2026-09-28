#include <cmath>
#include <iostream>
#include <vector>

#include "benchmark/core/CSVReporter.hpp"
#include "benchmark/core/Runner.hpp"

#include "linalg/operator/CSROperator.hpp"
#include "linalg/solver/iterative/gmres/Solver.hpp"
#include "linalg/solver/preconditioner/Identity.hpp"
#include "linalg/solver/preconditioner/Jacobi.hpp"
#include "linalg/types/CSRMatrix.hpp"
#include "linalg/types/Vector.hpp"
#include "linalg/types/backend/CPU.hpp"
#include "utils/logging/core/NullLogger.hpp"

using namespace residuum;

namespace {

	// nonsymmetric tridiagonal with a badly-scaled diagonal (10^(6*i/(n-1))) and
	// asymmetric off-diagonal coupling -- CG is not valid here, this is exactly the
	// shape of system GMRES was built for (see the ubend_nonlinear Newton tangent case)
	template<typename MatT>
	MatT makeBadlyScaledNonsymmetricTridiagonal(Index n) {

		MatT A(n, n);
		A.resize(3*n - 2);

		const Real lowerOffDiag = 0.15;
		const Real upperOffDiag = 0.05;

		Index* rowPtr = A.rowPtr();
		Index* colIdx = A.colIdx();
		Real* data = A.data();

		Index p = 0;
		for (Index i = 0; i < n; ++i) {

			rowPtr[i] = p;

			const Real diagVal = std::pow(10.0, 6.0 * static_cast<Real>(i) / static_cast<Real>(n - 1));

			if (i > 0) { colIdx[p] = i - 1; data[p] = lowerOffDiag; ++p; }
			colIdx[p] = i; data[p] = diagVal; ++p;
			if (i + 1 < n) { colIdx[p] = i + 1; data[p] = upperOffDiag; ++p; }

		}
		rowPtr[n] = p;

		return A;

	}

	template<typename PreconditionerT, typename OperatorT, typename VecT>
	benchmark::core::Result runGMRES(const std::string& backendName, const std::string& preconditionerName, const OperatorT& op, Index n, Index krylovDim, Index warmupIters, Index measureIters) {

		VecT b(n);
		for (Index i = 0; i < n; ++i) b.data()[i] = 1.0;

		VecT x(n);
		linalg::solver::iterative::gmres::Workspace<VecT> W(n, krylovDim);
		linalg::solver::SolverReport<VecT> report;
		PreconditionerT M;
		utils::logging::NullLogger logger;

		linalg::solver::iterative::gmres::Config<VecT> cfg{Real(1e-8)};
		cfg.krylovDim = krylovDim;
		cfg.maxIters = 1000;

		linalg::solver::iterative::gmres::Solver<OperatorT, VecT, PreconditionerT, decltype(logger)> solver(cfg);

		Index lastIterations = 0;

		auto result = benchmark::core::run("GMRESSolve", "", {{"backend", backendName}, {"n", std::to_string(n)}, {"preconditioner", preconditionerName}, {"krylov_dim", std::to_string(krylovDim)}}, warmupIters, measureIters, [&]() {
			x.zero();
			solver.solve(report, logger, W, M, op, b, x);
			lastIterations = report.iterations;
		});

		if (!report.converged) {
			std::cerr << "GMRESSolve: warning -- " << preconditionerName << " (m=" << krylovDim << ") did not converge for n=" << n << std::endl;
		}

		// amortized per-Arnoldi-step cost: operator + preconditioner apply, plus
		// orthogonalization against the average-sized basis within a restart cycle
		const Real avgOrthoCost = Real(4) * (static_cast<Real>(krylovDim) / Real(2) + Real(1)) * static_cast<Real>(n);
		result.estimatedFlops = static_cast<Real>(lastIterations) * (static_cast<Real>(op.flopsPerApply()) + static_cast<Real>(M.flopsPerApply()) + avgOrthoCost);

		return result;

	}

	// see CGSolve.cpp's runSuite for why this is backend-templated even though only
	// CPU is instantiated below
	template<typename BackendT>
	void runSuite(benchmark::core::CSVReporter& reporter, const std::string& backendName) {

		using Vec = linalg::types::Vector<Real, BackendT>;
		using Mat = linalg::types::CSRMatrix<Real, BackendT>;

		const std::vector<Index> sizes = {200, 800};
		const std::vector<Index> krylovDims = {10, 25, 50};
		const Index warmupIters = 1;
		const Index measureIters = 5;

		for (Index n : sizes) {

			Mat A = makeBadlyScaledNonsymmetricTridiagonal<Mat>(n);
			linalg::op::CSROperator<Mat> op(A);

			for (Index m : krylovDims) {

				auto identityResult = runGMRES<linalg::solver::preconditioner::Identity<Vec>, decltype(op), Vec>(backendName, "identity", op, n, m, warmupIters, measureIters);
				reporter.append(identityResult);

				auto jacobiResult = runGMRES<linalg::solver::preconditioner::Jacobi<Vec>, decltype(op), Vec>(backendName, "jacobi", op, n, m, warmupIters, measureIters);
				reporter.append(jacobiResult);

				std::cout << "backend=" << backendName << " n=" << n << " m=" << m << " identity: " << identityResult.meanSeconds << "s, jacobi: " << jacobiResult.meanSeconds << "s" << std::endl;

			}

		}

	}

} // namespace

int main() {

	benchmark::core::CSVReporter reporter(std::string(BENCHMARK_RESULTS_PATH) + "/GMRESSolve.csv");

	runSuite<linalg::types::backend::CPU>(reporter, "cpu");
	// runSuite<linalg::types::backend::CUDA>(reporter, "cuda"); // once the CUDA backend lands

	return 0;

}
