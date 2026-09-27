#ifndef RESIDUUM_BENCHMARK_CORE_RESULT_HPP
#define RESIDUUM_BENCHMARK_CORE_RESULT_HPP

#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "core/Types.hpp"

namespace residuum {
	namespace benchmark {
		namespace core {

			struct Result {

				std::string kernelName;

				// which equation module this result is specific to (e.g. "heateq"), empty
				// for equation-agnostic kernels (pure linalg/mesh/fem-layer benchmarks)
				std::string equation;

				// kernel-supplied axis name/value pairs, e.g. {"n","1000"}, {"preconditioner","jacobi"}
				std::vector<std::pair<std::string, std::string>> axisValues;

				Index warmupIters = 0;
				Index measureIters = 0;

				Real minSeconds = 0.0;
				Real meanSeconds = 0.0;
				Real medianSeconds = 0.0;

				// total FLOPs for one measured call, if the kernel can estimate it
				// (e.g. iterations * operator.flopsPerApply())
				std::optional<Real> estimatedFlops;

				std::optional<Real> gflopsPerSecond() const {
					if (!estimatedFlops || meanSeconds <= 0.0) return std::nullopt;
					return (*estimatedFlops) / meanSeconds / 1e9;
				}

			}; // struct Result

		} // namespace core
	} // namespace benchmark
} // namespace residuum

#endif
