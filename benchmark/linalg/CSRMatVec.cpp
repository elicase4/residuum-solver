#include <iostream>
#include <vector>

#include "benchmark/core/CSVReporter.hpp"
#include "benchmark/core/Runner.hpp"

#include "linalg/operator/CSROperator.hpp"
#include "linalg/types/CSRMatrix.hpp"
#include "linalg/types/Vector.hpp"
#include "linalg/types/backend/CPU.hpp"

using namespace residuum;

namespace {

	// plain symmetric tridiagonal (1D discrete Laplacian-like), diag=2, offdiag=-1 --
	// conditioning is irrelevant here, only nnz/row structure matters for matvec cost
	template<typename MatT>
	MatT makeTridiagonal(Index n) {

		MatT A(n, n);
		A.resize(3*n - 2);

		Index* rowPtr = A.rowPtr();
		Index* colIdx = A.colIdx();
		Real* data = A.data();

		Index p = 0;
		for (Index i = 0; i < n; ++i) {

			rowPtr[i] = p;

			if (i > 0) { colIdx[p] = i - 1; data[p] = -1.0; ++p; }
			colIdx[p] = i; data[p] = 2.0; ++p;
			if (i + 1 < n) { colIdx[p] = i + 1; data[p] = -1.0; ++p; }

		}
		rowPtr[n] = p;

		return A;

	}

	// see CGSolve.cpp's runSuite for why this is backend-templated even though only
	// CPU is instantiated below
	template<typename BackendT>
	void runSuite(benchmark::core::CSVReporter& reporter, const std::string& backendName) {

		using Vec = linalg::types::Vector<Real, BackendT>;
		using Mat = linalg::types::CSRMatrix<Real, BackendT>;

		const std::vector<Index> sizes = {1000, 10000, 100000, 1000000};
		const Index warmupIters = 10;
		const Index measureIters = 100;

		for (Index n : sizes) {

			Mat A = makeTridiagonal<Mat>(n);
			linalg::op::CSROperator<Mat> op(A);

			Vec x(n), y(n);
			for (Index i = 0; i < n; ++i) x.data()[i] = 1.0;

			auto result = benchmark::core::run("CSRMatVec", "", {{"backend", backendName}, {"n", std::to_string(n)}}, warmupIters, measureIters, [&]() {
				op.apply(x, y);
			}, static_cast<Real>(op.flopsPerApply()));

			reporter.append(result);

			std::cout << "backend=" << backendName << " n=" << n << " mean=" << result.meanSeconds << "s gflops/s=" << result.gflopsPerSecond().value_or(0.0) << std::endl;

		}

	}

} // namespace

int main() {

	benchmark::core::CSVReporter reporter(std::string(BENCHMARK_RESULTS_PATH) + "/CSRMatVec.csv");

	runSuite<linalg::types::backend::CPU>(reporter, "cpu");
	// runSuite<linalg::types::backend::CUDA>(reporter, "cuda"); // once the CUDA backend lands

	return 0;

}
