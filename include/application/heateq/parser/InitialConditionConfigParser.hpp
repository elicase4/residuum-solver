#ifndef RESIDUUM_APPLICATION_HEATEQ_PARSER_INITIALCONDITIONCONFIGPARSER_HPP
#define RESIDUUM_APPLICATION_HEATEQ_PARSER_INITIALCONDITIONCONFIGPARSER_HPP

#include <string>
#include <stdexcept>
#include <yaml-cpp/yaml.h>

#include "application/heateq/config/InitialConditionConfig.hpp"

namespace residuum {
	namespace application {
		namespace heateq {
			namespace parser {

				class InitialConditionConfigParser {
				public:

					static config::InitialConditionConfig parse(const YAML::Node& node);

				}; // class InitialConditionConfigParser

			} // namespace parser
		} // namespace heateq
	} // namespace application
} // namespace residuum

#endif
