#ifndef RESIDUUM_SOLVER_CONFIG_MONITORTERMCONFIG_HPP
#define RESIDUUM_SOLVER_CONFIG_MONITORTERMCONFIG_HPP

#include <vector>

#include "core/Types.hpp"

namespace residuum {
	namespace solver {
		namespace config {

			struct MonitorTermConfig {

				Int boundary;
				Real coefficient = 1.0;

				// per-component weights for a weighted/projection term (e.g. lift); empty for an
				// ordinary term, which uses coefficient uniformly instead
				std::vector<Real> weights;

			}; // struct MonitorTermConfig

		} // namespace config
	} // namespace solver
} // namespace residuum

#endif
