#ifndef RESIDUUM_APPLICATION_HEATEQ_PARSER_MONITORCONFIGPARSER_HPP
#define RESIDUUM_APPLICATION_HEATEQ_PARSER_MONITORCONFIGPARSER_HPP

#include <string>
#include <yaml-cpp/yaml.h>

#include "application/heateq/config/MonitorConfig.hpp"

namespace residuum {
	namespace application {
		namespace heateq {
			namespace parser {

				class MonitorConfigParser {
				public:

					static config::MonitorConfig::Quantity parseQuantity(const std::string& str);

					static config::MonitorConfig parse(const YAML::Node& node);

				}; // class MonitorConfigParser

			} // namespace parser
		} // namespace heateq
	} // namespace application
} // namespace residuum

#endif
