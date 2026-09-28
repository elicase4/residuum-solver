#ifndef RESIDUUM_SOLVER_PARSER_TIMESTEPSIZECONFIGPARSER_HPP
#define RESIDUUM_SOLVER_PARSER_TIMESTEPSIZECONFIGPARSER_HPP

#include <string>

#include <yaml-cpp/yaml.h>

#include "solver/config/TimeStepSizeConfig.hpp"

namespace residuum {
	namespace solver {
		namespace parser {

			class TimeStepSizeConfigParser {
			public:

				static config::TimeStepSizeConfig::Mode parseMode(const std::string& str);

				static config::TimeStepSizeConfig parse(const YAML::Node& node);

			}; // class TimeStepSizeConfigParser

		} // namespace parser
	} // namespace solver
} // namespace residuum

#endif
