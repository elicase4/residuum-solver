#ifndef RESIDUUM_SOLVER_CONFIG_TIMESTEPSIZECONFIG_HPP
#define RESIDUUM_SOLVER_CONFIG_TIMESTEPSIZECONFIG_HPP

#include <optional>

#include "core/Types.hpp"

namespace residuum {
	namespace solver {
		namespace config {

			struct TimeStepSizeConfig {

				enum class Mode { Constant, Adaptive, PseudoTransient };

				Mode mode = Mode::Constant;

				// Constant.
				Real dt = 1e-3;

				// Adaptive-only
				std::optional<Real> dtMin;
				std::optional<Real> dtMax;
				std::optional<Real> growthFactor;
				std::optional<Real> shrinkFactor;
				std::optional<Real> errorTolerance;

				// PseudoTransient-only
				std::optional<Real> convergenceTolerance;

			}; // struct TimeStepSizeConfig

		} // namespace config
	} // namespace solver
} // namespace residuum

#endif
