#include <iostream>
#include <vector>

#include "benchmark/core/CSVReporter.hpp"
#include "benchmark/core/Runner.hpp"

#include "linalg/operator/CSROperator.hpp"
#include "linalg/solver/preconditioner/Identity.hpp"
#include "linalg/solver/preconditioner/Jacobi.hpp"
#include "linalg/types/CSRMatrix.hpp"
#include "linalg/types/Vector.hpp"
#include "linalg/types/backend/CPU.hpp"

using namespace residuum;

namespace {

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

	// isolates update()/apply() cost -- this is what a preconditioner adds on top of
	// whatever it saves in iteration count, so it should be small next to CGSolve/GMRESSolve
	template<typename PreconditionerT, typename OperatorT, typename VecT>
	void runOverhead(benchmark::core::CSVReporter& reporter, const std::string& backendName, const std::string& preconditionerName, const OperatorT& op, VecT& r, VecT& z, Index n, Index warmupIters, Index measureIters) {

		PreconditionerT M;

		auto updateResult = benchmark::core::run("JacobiOverhead", "", {{"backend", backendName}, {"n", std::to_string(n)}, {"preconditioner", preconditionerName}, {"operation", "update"}}, warmupIters, measureIters, [&]() {
			M.update(op);
		});
		reporter.append(updateResult);

		M.update(op); // apply() needs a valid state; measured separately from update() above

		auto applyResult = benchmark::core::run("JacobiOverhead", "", {{"backend", backendName}, {"n", std::to_string(n)}, {"preconditioner", preconditionerName}, {"operation", "apply"}}, warmupIters, measureIters, [&]() {
			M.apply(r, z);
		}, static_cast<Real>(M.flopsPerApply()));
		reporter.append(applyResult);

		std::cout << "backend=" << backendName << " n=" << n << " " << preconditionerName << ": update=" << updateResult.meanSeconds << "s apply=" << applyResult.meanSeconds << "s" << std::endl;

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

			Vec r(n), z(n);
			for (Index i = 0; i < n; ++i) r.data()[i] = 1.0;

			runOverhead<linalg::solver::preconditioner::Identity<Vec>, decltype(op), Vec>(reporter, backendName, "identity", op, r, z, n, warmupIters, measureIters);
			runOverhead<linalg::solver::preconditioner::Jacobi<Vec>, decltype(op), Vec>(reporter, backendName, "jacobi", op, r, z, n, warmupIters, measureIters);

		}

	}

} // namespace

int main() {

	benchmark::core::CSVReporter reporter(std::string(BENCHMARK_RESULTS_PATH) + "/JacobiOverhead.csv");

	runSuite<linalg::types::backend::CPU>(reporter, "cpu");
	// runSuite<linalg::types::backend::CUDA>(reporter, "cuda"); // once the CUDA backend lands

	return 0;

}
