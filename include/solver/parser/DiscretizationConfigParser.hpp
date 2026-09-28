#ifndef RESIDUUM_SOLVER_PARSER_DISCRETIZATIONCONFIGPARSER_HPP
#define RESIDUUM_SOLVER_PARSER_DISCRETIZATIONCONFIGPARSER_HPP

#include <string>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "solver/config/DiscretizationConfig.hpp"

namespace residuum {
	namespace solver {
		namespace parser {
	
			class DiscretizationConfigParser {
			public:

				static config::DiscretizationConfig parse(const YAML::Node& node);

			}; // class DiscretizationConfigParser

		} // namespace parser
	} // namespace solver
} // namespace residuum

#endif
