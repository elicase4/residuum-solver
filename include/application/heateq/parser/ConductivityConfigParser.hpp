#ifndef RESIDUUM_APPLICATION_HEATEQ_PARSER_CONDUCTIVITYCONFIGPARSER_HPP
#define RESIDUUM_APPLICATION_HEATEQ_PARSER_CONDUCTIVITYCONFIGPARSER_HPP

#include <string>
#include <stdexcept>
#include <vector>
#include <yaml-cpp/yaml.h>

#include "application/heateq/config/ConductivityConfig.hpp"

namespace residuum {
	namespace application {
		namespace heateq {
			namespace parser {

				class ConductivityConfigParser {
				public:

					static config::ConductivityConfig::Type parseConductivityType(const std::string& str);

					static config::ConductivityConfig parse(const YAML::Node& node);
				
				}; // class ConductivityConfigParser

			} // namespace parser
		} // namespace heateq
	} // namespace application
} // namespace residuum

#endif
