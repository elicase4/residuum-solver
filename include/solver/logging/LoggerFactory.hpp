#ifndef RESIDUUM_SOLVER_LOGGING_LOGGERFACTORY_HPP
#define RESIDUUM_SOLVER_LOGGING_LOGGERFACTORY_HPP

#include <filesystem>
#include <string>
#include <utility>
#include <vector>

#include "core/Types.hpp"

#include "fem/dof/DOFOrdering.hpp"

#include "solver/config/LoggingConfig.hpp"
#include "solver/config/TimeStepperConfig.hpp"

#include "utils/logging/core/NullLogger.hpp"
#include "utils/logging/driver/Logger.hpp"
#include "utils/logging/linear/Logger.hpp"
#include "utils/logging/nonlinear/Logger.hpp"
#include "utils/logging/timestepper/Logger.hpp"

namespace residuum {
	namespace solver {
		namespace logging {

			inline void createParentDirectory(const std::string& filePath) {

				if (filePath.empty()) return;

				const std::filesystem::path parent = std::filesystem::path(filePath).parent_path();
				if (!parent.empty()) std::filesystem::create_directories(parent);

			}

			// every logger config in LoggingConfig may name a text/CSV file; call before any make*Logger
			inline void ensureLoggingDirectories(const config::LoggingConfig& cfg) {

				createParentDirectory(cfg.driver.textFile);

				createParentDirectory(cfg.linear.textFile);
				createParentDirectory(cfg.linear.csvFile);

				createParentDirectory(cfg.timestepper.textFile);
				createParentDirectory(cfg.timestepper.csvFile);

				createParentDirectory(cfg.nonlinear.textFile);
				createParentDirectory(cfg.nonlinear.csvFile);

			}

			inline utils::logging::driver::Logger makeDriverLogger(const config::LoggerConfig& cfg, std::string tag) {

				const bool consoleEnabled = (cfg.type == config::LoggerConfig::Type::Console);

				if (consoleEnabled || !cfg.textFile.empty()) {
					return utils::logging::driver::Logger(utils::logging::driver::ConsoleLogger{std::move(tag), consoleEnabled, cfg.textFile});
				}

				return utils::logging::driver::Logger(utils::logging::NullLogger{});

			}

			inline utils::logging::linear::Logger makeLinearLogger(const config::LinearLoggerConfig& loggerCfg, const std::string& equationName, const std::string& solverName, const std::string& preconditionerName, const std::vector<std::string>& dofNames, const std::vector<std::pair<std::string, std::string>>& extraParams, Index freeDOFsPerField, fem::dof::DOFOrdering ordering) {

				using LoggerT = utils::logging::linear::Logger;

				const bool consoleEnabled = (loggerCfg.type == config::LoggerConfig::Type::Console);
				const bool anyOutput = consoleEnabled || !loggerCfg.textFile.empty() || !loggerCfg.csvFile.empty();

				if (!anyOutput) return LoggerT(utils::logging::NullLogger{});

				return LoggerT(utils::logging::linear::ConsoleLogger(equationName, solverName, preconditionerName, dofNames, extraParams, freeDOFsPerField, ordering, loggerCfg.interval, consoleEnabled, loggerCfg.textFile, loggerCfg.csvFile));

			}

			inline utils::logging::nonlinear::Logger makeNonlinearLogger(const config::NonlinearLoggerConfig& loggerCfg, const std::string& equationLabel, const std::string& solverName, const std::vector<std::string>& dofNames) {

				using LoggerT = utils::logging::nonlinear::Logger;

				const bool consoleEnabled = (loggerCfg.type == config::LoggerConfig::Type::Console);
				const bool anyOutput = consoleEnabled || !loggerCfg.textFile.empty() || !loggerCfg.csvFile.empty();

				if (!anyOutput) return LoggerT(utils::logging::NullLogger{});

				return LoggerT(utils::logging::nonlinear::ConsoleLogger(equationLabel, solverName, dofNames, loggerCfg.interval, consoleEnabled, loggerCfg.textFile, loggerCfg.csvFile));

			}

			inline utils::logging::timestepper::Logger makeTimestepperLogger(const config::TimeStepperConfig& cfg, const config::TimeStepperLoggerConfig& loggerCfg, const std::string& equationLabel, const std::string& timestepperName) {

				using LoggerT = utils::logging::timestepper::Logger;

				const bool consoleEnabled = (loggerCfg.type == config::LoggerConfig::Type::Console);
				const bool anyOutput = consoleEnabled || !loggerCfg.textFile.empty() || !loggerCfg.csvFile.empty();

				if (!anyOutput) return LoggerT(utils::logging::NullLogger{});

				return LoggerT(utils::logging::timestepper::ConsoleLogger(equationLabel, timestepperName, cfg.t0, cfg.tf, loggerCfg.interval, consoleEnabled, loggerCfg.textFile, loggerCfg.csvFile));

			}

		} // namespace logging
	} // namespace solver
} // namespace residuum

#endif
