#include <cmath>
#include <iostream>
#include <vector>

#include "benchmark/core/CSVReporter.hpp"
#include "benchmark/core/Runner.hpp"

#include "linalg/operator/CSROperator.hpp"
#include "linalg/solver/iterative/cg/Solver.hpp"
#include "linalg/solver/preconditioner/Identity.hpp"
#include "linalg/solver/preconditioner/Jacobi.hpp"
#include "linalg/types/CSRMatrix.hpp"
#include "linalg/types/Vector.hpp"
#include "linalg/types/backend/CPU.hpp"
#include "utils/logging/core/NullLogger.hpp"

using namespace residuum;

namespace {

	// symmetric tridiagonal SPD matrix whose diagonal spans several orders of magnitude
	// (diag[i] = 10^(6*i/(n-1)), weak constant off-diagonal coupling) -- badly conditioned
	// for plain CG, but collapses to a well-conditioned system once Jacobi-preconditioned
	template<typename MatT>
	MatT makeBadlyScaledTridiagonal(Index n) {

		MatT A(n, n);
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

	template<typename PreconditionerT, typename OperatorT, typename VecT>
	benchmark::core::Result runCG(const std::string& backendName, const std::string& preconditionerName, const OperatorT& op, Index n, Index warmupIters, Index measureIters) {

		VecT b(n);
		for (Index i = 0; i < n; ++i) b.data()[i] = 1.0;

		VecT x(n);
		linalg::solver::iterative::cg::Workspace<VecT> W(n);
		linalg::solver::SolverReport<VecT> report;
		PreconditionerT M;
		utils::logging::NullLogger logger;

		linalg::solver::iterative::cg::Config<VecT> cfg{Real(1e-8)};
		cfg.maxIters = 5000;

		linalg::solver::iterative::cg::Solver<OperatorT, VecT, PreconditionerT, decltype(logger)> solver(cfg);

		Index lastIterations = 0;

		auto result = benchmark::core::run("CGSolve", "", {{"backend", backendName}, {"n", std::to_string(n)}, {"preconditioner", preconditionerName}}, warmupIters, measureIters, [&]() {
			x.zero();
			solver.solve(report, logger, W, M, op, b, x);
			lastIterations = report.iterations;
		});

		if (!report.converged) {
			std::cerr << "CGSolve: warning -- " << preconditionerName << " did not converge for n=" << n << std::endl;
		}

		result.estimatedFlops = static_cast<Real>(lastIterations) * (static_cast<Real>(op.flopsPerApply()) + static_cast<Real>(M.flopsPerApply()) + Real(10) * static_cast<Real>(n));

		return result;

	}

	// runs the full CGSolve sweep for one backen
	template<typename BackendT>
	void runSuite(benchmark::core::CSVReporter& reporter, const std::string& backendName) {

		using Vec = linalg::types::Vector<Real, BackendT>;
		using Mat = linalg::types::CSRMatrix<Real, BackendT>;

		const std::vector<Index> sizes = {50, 200, 800, 3200};
		const Index warmupIters = 2;
		const Index measureIters = 10;

		for (Index n : sizes) {

			Mat A = makeBadlyScaledTridiagonal<Mat>(n);
			linalg::op::CSROperator<Mat> op(A);

			auto identityResult = runCG<linalg::solver::preconditioner::Identity<Vec>, decltype(op), Vec>(backendName, "identity", op, n, warmupIters, measureIters);
			reporter.append(identityResult);

			auto jacobiResult = runCG<linalg::solver::preconditioner::Jacobi<Vec>, decltype(op), Vec>(backendName, "jacobi", op, n, warmupIters, measureIters);
			reporter.append(jacobiResult);

			std::cout << "backend=" << backendName << " n=" << n << " identity: " << identityResult.meanSeconds << "s, jacobi: " << jacobiResult.meanSeconds << "s" << std::endl;

		}

	}

} // namespace

int main() {

	benchmark::core::CSVReporter reporter(std::string(BENCHMARK_RESULTS_PATH) + "/CGSolve.csv");

	runSuite<linalg::types::backend::CPU>(reporter, "cpu");
	// runSuite<linalg::types::backend::CUDA>(reporter, "cuda"); // once the CUDA backend lands

	return 0;

}
