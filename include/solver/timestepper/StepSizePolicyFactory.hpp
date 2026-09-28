#ifndef RESIDUUM_SOLVER_TIMESTEPPER_STEPSIZEPOLICYFACTORY_HPP
#define RESIDUUM_SOLVER_TIMESTEPPER_STEPSIZEPOLICYFACTORY_HPP

#include <memory>
#include <stdexcept>
#include <string>

#include "solver/config/TimeStepSizeConfig.hpp"
#include "solver/timestepper/ConstantStepSizePolicy.hpp"
#include "solver/timestepper/TimeStepSizePolicy.hpp"

namespace residuum {
	namespace solver {
		namespace timestepper {

			inline std::unique_ptr<TimeStepSizePolicy> makeStepSizePolicy(const config::TimeStepSizeConfig& cfg, const std::string& equationLabel) {

				switch (cfg.mode) {

					case config::TimeStepSizeConfig::Mode::Constant:
						return std::make_unique<ConstantStepSizePolicy>(cfg.dt);

					case config::TimeStepSizeConfig::Mode::Adaptive:
						throw std::runtime_error("StepSizePolicyFactory[" + equationLabel + "]: adaptive step size not yet implemented");

					case config::TimeStepSizeConfig::Mode::PseudoTransient:
						throw std::runtime_error("StepSizePolicyFactory[" + equationLabel + "]: pseudo-transient step size not yet implemented");

				}

				throw std::runtime_error("StepSizePolicyFactory[" + equationLabel + "]: unknown step size mode");

			}

		} // namespace timestepper
	} // namespace solver
} // namespace residuum

#endif
