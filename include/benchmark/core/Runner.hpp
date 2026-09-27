#ifndef RESIDUUM_BENCHMARK_CORE_RUNNER_HPP
#define RESIDUUM_BENCHMARK_CORE_RUNNER_HPP

#include <algorithm>
#include <numeric>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "benchmark/core/Result.hpp"
#include "benchmark/core/Timer.hpp"

#include "core/Types.hpp"

namespace residuum {
	namespace benchmark {
		namespace core {

			// runs callable() warmupIters times (discarded), then measureIters times,
			// timing each measured call individually; callable takes no arguments
			template<typename CallableT>
			Result run(const std::string& kernelName, const std::string& equation, const std::vector<std::pair<std::string, std::string>>& axisValues, Index warmupIters, Index measureIters, CallableT&& callable, std::optional<Real> estimatedFlopsPerCall = std::nullopt) {

				for (Index i = 0; i < warmupIters; ++i) callable();

				Timer timer;
				std::vector<Real> samples(measureIters);

				for (Index i = 0; i < measureIters; ++i) {
					timer.start();
					callable();
					samples[i] = timer.stop();
				}

				Result result;
				result.kernelName = kernelName;
				result.equation = equation;
				result.axisValues = axisValues;
				result.warmupIters = warmupIters;
				result.measureIters = measureIters;

				result.minSeconds = *std::min_element(samples.begin(), samples.end());
				result.meanSeconds = std::accumulate(samples.begin(), samples.end(), Real(0)) / static_cast<Real>(measureIters);

				std::vector<Real> sorted = samples;
				std::sort(sorted.begin(), sorted.end());
				result.medianSeconds = sorted[sorted.size() / 2];

				result.estimatedFlops = estimatedFlopsPerCall;

				return result;

			}

		} // namespace core
	} // namespace benchmark
} // namespace residuum

#endif
