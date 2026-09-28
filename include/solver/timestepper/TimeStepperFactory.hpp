#ifndef RESIDUUM_SOLVER_TIMESTEPPER_TIMESTEPPERFACTORY_HPP
#define RESIDUUM_SOLVER_TIMESTEPPER_TIMESTEPPERFACTORY_HPP

#include <memory>
#include <stdexcept>
#include <string>

#include "solver/config/LoggingConfig.hpp"
#include "solver/config/TimeStepperConfig.hpp"
#include "solver/logging/LoggerFactory.hpp"
#include "solver/problem/TransientCapableProblem.hpp"
#include "solver/stage/BackwardEulerStage.hpp"
#include "solver/timestepper/BackwardEuler.hpp"
#include "solver/timestepper/StepSizePolicyFactory.hpp"
#include "solver/timestepper/TimeStepper.hpp"
#include "solver/timestepper/TimeStepperRunner.hpp"

#include "utils/logging/core/NullLogger.hpp"
#include "utils/logging/timestepper/Logger.hpp"

namespace residuum {
	namespace solver {
		namespace timestepper {

			template<typename StageT>
			std::unique_ptr<TimeStepperRunner> makeTimeStepperRunner(StageT& stage, const config::TimeStepperConfig& cfg, const config::TimeStepperLoggerConfig& loggerCfg, const std::string& equationLabel) {

				switch (cfg.type) {

					case config::TimeStepperConfig::Type::BackwardEuler:
						static_assert(TimeStepper<BackwardEulerRunner<StageT>>, "BackwardEulerRunner<StageT> must satisfy the TimeStepper concept, including declaring TemporalOrder Order");
						return std::make_unique<BackwardEulerRunner<StageT>>(stage, cfg, makeStepSizePolicy(cfg.stepSize, equationLabel), logging::makeTimestepperLogger(cfg, loggerCfg, equationLabel, "Backward Euler"));

					case config::TimeStepperConfig::Type::ForwardEuler:
						throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: ForwardEuler not yet implemented");

					case config::TimeStepperConfig::Type::GeneralizedAlpha:
						throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: GeneralizedAlpha not yet implemented");

					case config::TimeStepperConfig::Type::RK4:
						throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: RK4 not yet implemented");

				}

				throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: unknown timestepper type");

			}

			template<problem::TransientCapableProblem ProblemT, typename VisitorT>
			auto dispatchTransientStage(ProblemT& problem, const config::TimeStepperConfig& cfg, const config::TimeStepperLoggerConfig& loggerCfg, const std::string& equationLabel, VisitorT&& visitor) {

				switch (cfg.type) {

					case config::TimeStepperConfig::Type::BackwardEuler: {
						using StageT = stage::BackwardEulerStage<ProblemT>;
						StageT stage(problem);
						auto runner = makeTimeStepperRunner<StageT>(stage, cfg, loggerCfg, equationLabel);
						return visitor(stage, *runner);
					}

					case config::TimeStepperConfig::Type::ForwardEuler:
						throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: ForwardEuler not yet implemented");

					case config::TimeStepperConfig::Type::GeneralizedAlpha:
						throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: GeneralizedAlpha not yet implemented");

					case config::TimeStepperConfig::Type::RK4:
						throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: RK4 not yet implemented");

				}

				throw std::runtime_error("TimeStepperFactory[" + equationLabel + "]: unknown timestepper type");

			}

		} // namespace timestepper
	} // namespace solver
} // namespace residuum

#endif
