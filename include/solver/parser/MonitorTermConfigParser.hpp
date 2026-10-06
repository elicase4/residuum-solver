#ifndef RESIDUUM_SOLVER_PARSER_MONITORTERMCONFIGPARSER_HPP
#define RESIDUUM_SOLVER_PARSER_MONITORTERMCONFIGPARSER_HPP

#include <yaml-cpp/yaml.h>

#include "solver/config/MonitorTermConfig.hpp"

namespace residuum {
	namespace solver {
		namespace parser {

			class MonitorTermConfigParser {
			public:

				static config::MonitorTermConfig parse(const YAML::Node& node);

			}; // class MonitorTermConfigParser

		} // namespace parser
	} // namespace solver
} // namespace residuum

#endif
