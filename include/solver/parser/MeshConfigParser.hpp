#ifndef RESIDUUM_SOLVER_PARSER_MESHCONFIGPARSER_HPP
#define RESIDUUM_SOLVER_PARSER_MESHCONFIGPARSER_HPP

#include <string>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "solver/config/MeshConfig.hpp"

namespace residuum {
	namespace solver {
		namespace parser {
	
			class MeshConfigParser {
			public:

				static config::MeshConfig parse(const YAML::Node& node);

			}; // class MeshConfigParser

		} // namespace parser
	} // namespace solver
} // namespace residuum

#endif
