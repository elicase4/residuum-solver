#ifndef RESIDUUM_UTILS_LOGGING_CORE_NULLLOGGER_HPP
#define RESIDUUM_UTILS_LOGGING_CORE_NULLLOGGER_HPP

#include <vector>

#include "core/Types.hpp"

namespace residuum {
	namespace utils {
		namespace logging {

			struct NullLogger {

				template<typename DataT>
				inline void log(Index, const std::vector<DataT>&, DataT = DataT(0)) const {}

				// timestepper::ConsoleLogger's log() shape: step, time, dt, attempts, residualNorm
				inline void log(Index, Real, Real, Index, Real) const {}

				// timestepper::ConsoleLogger's startStep() shape: step, time
				inline void startStep(Index, Real) const {}

				// nonlinear::ConsoleLogger's log() shape: iter, residualRelPerDOF
				inline void log(Index, const std::vector<Real>&) const {}

				template<typename DataT>
				inline std::vector<DataT> computePerDOFNorms(const DataT*, Index) const { return {}; }

				inline void reset() const {}

				inline void summary(bool) const {}

				template<typename Args>
				inline void event(Args&&) const {}

				template<typename Args>
				inline void warn(Args&&) const {}

				template<typename Args>
				inline void error(Args&&) const {}

			}; // struct NullLogger

		} // namespace logging
	} // namespace utils
} // namespace residuum

#endif
