#ifndef RESIDUUM_APPLICATION_HEATEQ_PARSER_BOUNDARYCONDITIONCONFIGPARSER_HPP
#define RESIDUUM_APPLICATION_HEATEQ_PARSER_BOUNDARYCONDITIONCONFIGPARSER_HPP

#include <string>
#include <stdexcept>
#include <yaml-cpp/yaml.h>

#include "application/heateq/config/BoundaryConditionConfig.hpp"

namespace residuum {
	namespace application {
		namespace heateq {
			namespace parser {

				class BoundaryConditionConfigParser {
				public:

					static config::BoundaryConditionConfig::Type parseBoundaryConditionType(const std::string& str);

					static config::BoundaryConditionConfig::Form parseBoundaryConditionForm(const std::string& str);
					
					static config::BoundaryConditionConfig parse(const YAML::Node& node);
					
				}; // class BoundaryConditionConfigParser

			} // namespace parser
		} // namespace heateq
	} // namespace application
} // namespace residuum

#endif
