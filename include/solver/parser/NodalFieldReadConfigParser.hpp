#ifndef RESIDUUM_SOLVER_PARSER_NODALFIELDREADCONFIGPARSER_HPP
#define RESIDUUM_SOLVER_PARSER_NODALFIELDREADCONFIGPARSER_HPP

#include <string>

#include <yaml-cpp/yaml.h>

#include "solver/config/NodalFieldReadConfig.hpp"

namespace residuum {
	namespace solver {
		namespace parser {

			class NodalFieldReadConfigParser {
			public:

				static config::NodalFieldReadConfig::Mode parseMode(const std::string& str);

				static config::NodalFieldReadConfig parse(const YAML::Node& node);

			}; // class NodalFieldReadConfigParser

		} // namespace parser
	} // namespace solver
} // namespace residuum

#endif
