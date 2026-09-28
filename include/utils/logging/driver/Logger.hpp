#ifndef RESIDUUM_UTILS_LOGGING_DRIVER_LOGGER_HPP
#define RESIDUUM_UTILS_LOGGING_DRIVER_LOGGER_HPP

#include <utility>
#include <variant>

#include "utils/logging/core/NullLogger.hpp"
#include "utils/logging/driver/ConsoleLogger.hpp"

namespace residuum {
	namespace utils {
		namespace logging {
			namespace driver {

				class Logger {
				public:

					using Variant = std::variant<ConsoleLogger, NullLogger>;

					explicit Logger(Variant impl) : impl_(std::move(impl)) {}

					template<typename Args>
					void event(Args&& msg) const {
						std::visit([&](const auto& l){ l.event(std::forward<Args>(msg)); }, impl_);
					}

					template<typename Args>
					void warn(Args&& msg) const {
						std::visit([&](const auto& l){ l.warn(std::forward<Args>(msg)); }, impl_);
					}

					template<typename Args>
					void error(Args&& msg) const {
						std::visit([&](const auto& l){ l.error(std::forward<Args>(msg)); }, impl_);
					}

				private:
					Variant impl_;

				}; // class Logger

			} // namespace driver
		} // namespace logging
	} // namespace utils
} // namespace residuum

#endif
