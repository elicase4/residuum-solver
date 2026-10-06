#ifndef RESIDUUM_APPLICATION_HEATEQ_PARSER_OUTPUTFIELDCONFIGPARSER_HPP
#define RESIDUUM_APPLICATION_HEATEQ_PARSER_OUTPUTFIELDCONFIGPARSER_HPP

#include <string>
#include <yaml-cpp/yaml.h>

#include "application/heateq/config/OutputFieldConfig.hpp"

namespace residuum {
	namespace application {
		namespace heateq {
			namespace parser {

				class OutputFieldConfigParser {
				public:

					static config::OutputFieldQuantity parseQuantity(const std::string& str);

					static config::OutputFieldConfig parse(const YAML::Node& node);

				}; // class OutputFieldConfigParser

			} // namespace parser
		} // namespace heateq
	} // namespace application
} // namespace residuum

#endif
