#ifndef RESIDUUM_SOLVER_PARSER_OUTPUTCONFIGPARSER_HPP
#define RESIDUUM_SOLVER_PARSER_OUTPUTCONFIGPARSER_HPP

#include <string>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "solver/config/OutputConfig.hpp"

namespace residuum {
	namespace solver {
		namespace parser {
	
			class OutputConfigParser {
			public:

				static config::OutputConfig parse(const YAML::Node& node);

			}; // class OutputConfigParser

		} // namespace parser
	} // namespace solver
} // namespace residuum

#endif
