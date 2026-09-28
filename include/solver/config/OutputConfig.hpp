#ifndef RESIDUUM_SOLVER_CONFIG_OUTPUTCONFIG_HPP
#define RESIDUUM_SOLVER_CONFIG_OUTPUTCONFIG_HPP

#include <string>

#include "core/Types.hpp"

namespace residuum {
	namespace solver {
		namespace config {

			struct OutputConfig {

				enum class Format { VTK, VTU };

				std::string directory = "output";

				Format format = Format::VTK;

				Index writeFrequency = 1;

				std::string prefix = "solution";

			}; // struct OutputConfig

		} // namespace config
	} // namespace solver
} // namespace residuum

#endif
