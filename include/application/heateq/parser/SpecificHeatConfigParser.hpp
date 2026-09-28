#ifndef RESIDUUM_APPLICATION_HEATEQ_PARSER_SPECIFICHEATCONFIGPARSER_HPP
#define RESIDUUM_APPLICATION_HEATEQ_PARSER_SPECIFICHEATCONFIGPARSER_HPP

#include <string>
#include <stdexcept>
#include <yaml-cpp/yaml.h>

#include "application/heateq/config/SpecificHeatConfig.hpp"

namespace residuum {
	namespace application {
		namespace heateq {
			namespace parser {

				class SpecificHeatConfigParser {
				public:

					static config::SpecificHeatConfig::Type parseSpecificHeatType(const std::string& str);

					static config::SpecificHeatConfig parse(const YAML::Node& node);

				}; // class SpecificHeatConfigParser

			} // namespace parser
		} // namespace heateq
	} // namespace application
} // namespace residuum

#endif
