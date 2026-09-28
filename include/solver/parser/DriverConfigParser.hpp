#ifndef RESIDUUM_SOLVER_PARSER_DRIVERCONFIGPARSER_HPP
#define RESIDUUM_SOLVER_PARSER_DRIVERCONFIGPARSER_HPP

#include <string>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "solver/config/DriverConfig.hpp"

namespace residuum {
	namespace solver {
		namespace parser {
	
			class DriverConfigParser {
			public:

				static config::DriverConfig::Type parseDriverType(const std::string& str);

				static config::DriverConfig parse(const YAML::Node& node);

			}; // class DriverConfigParser

		} // namespace parser
	} // namespace solver
} // namespace residuum

#endif
