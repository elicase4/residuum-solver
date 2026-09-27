#ifndef RESIDUUM_BENCHMARK_CORE_TIMER_HPP
#define RESIDUUM_BENCHMARK_CORE_TIMER_HPP

#include <chrono>

namespace residuum {
	namespace benchmark {
		namespace core {

			class Timer {
			public:

				void start() { start_ = Clock::now(); }

				// seconds elapsed since the last start()
				double stop() {
					return std::chrono::duration<double>(Clock::now() - start_).count();
				}

			private:
				using Clock = std::chrono::steady_clock;
				Clock::time_point start_;

			}; // class Timer

		} // namespace core
	} // namespace benchmark
} // namespace residuum

#endif
