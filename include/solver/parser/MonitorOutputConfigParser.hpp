#ifndef RESIDUUM_SOLVER_PARSER_MONITOROUTPUTCONFIGPARSER_HPP
#define RESIDUUM_SOLVER_PARSER_MONITOROUTPUTCONFIGPARSER_HPP

#include <yaml-cpp/yaml.h>

#include "solver/config/MonitorOutputConfig.hpp"

namespace residuum {
	namespace solver {
		namespace parser {

			class MonitorOutputConfigParser {
			public:

				static config::MonitorOutputConfig parse(const YAML::Node& node);

			}; // class MonitorOutputConfigParser

		} // namespace parser
	} // namespace solver
} // namespace residuum

#endif
