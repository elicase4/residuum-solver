#ifndef RESIDUUM_APPLICATION_HEATEQ_CONFIG_SCALARMATERIALPROPERTYCONFIG_HPP
#define RESIDUUM_APPLICATION_HEATEQ_CONFIG_SCALARMATERIALPROPERTYCONFIG_HPP

#include <string>

#include "core/Types.hpp"

namespace residuum {
	namespace application {
		namespace heateq {
			namespace config {

				struct ScalarMaterialPropertyConfig {

					Real value;
					std::string unit;

				}; // struct ScalarMaterialPropertyConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace residuum

#endif
