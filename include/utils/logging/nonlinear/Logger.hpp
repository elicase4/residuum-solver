#ifndef RESIDUUM_UTILS_LOGGING_NONLINEAR_LOGGER_HPP
#define RESIDUUM_UTILS_LOGGING_NONLINEAR_LOGGER_HPP

#include <utility>
#include <variant>
#include <vector>

#include "core/Types.hpp"
#include "utils/logging/core/NullLogger.hpp"
#include "utils/logging/nonlinear/ConsoleLogger.hpp"

namespace residuum {
	namespace utils {
		namespace logging {
			namespace nonlinear {

				class Logger {
				public:

					using Variant = std::variant<ConsoleLogger, NullLogger>;

					explicit Logger(Variant impl) : impl_(std::move(impl)) {}

					void log(Index iter, const std::vector<Real>& residualRelPerDOF) const {
						std::visit([&](const auto& l){ l.log(iter, residualRelPerDOF); }, impl_);
					}

					void reset() const {
						std::visit([&](const auto& l){ l.reset(); }, impl_);
					}

					void summary(bool converged) const {
						std::visit([&](const auto& l){ l.summary(converged); }, impl_);
					}

				private:
					Variant impl_;

				}; // class Logger

			} // namespace nonlinear
		} // namespace logging
	} // namespace utils
} // namespace residuum

#endif
