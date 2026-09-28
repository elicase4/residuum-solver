#ifndef RESIDUUM_APPLICATION_HEATEQ_PARSER_HEATCONFIGPARSER_HPP
#define RESIDUUM_APPLICATION_HEATEQ_PARSER_HEATCONFIGPARSER_HPP

#include <string>

#include "application/heateq/config/HeatConfig.hpp"

namespace residuum {
	namespace application {
		namespace heateq {
			namespace parser {

			class HeatConfigParser {
			public:

				static config::HeatConfig read(const std::string& filename);

			}; // HeatConfigParser

			} // namespace parser
		} // namespace heateq
	} // namespace application
} // namespace residuum

#endif
