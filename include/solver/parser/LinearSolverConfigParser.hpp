#ifndef RESIDUUM_SOLVER_PARSER_LINEARSOLVERCONFIGPARSER_HPP
#define RESIDUUM_SOLVER_PARSER_LINEARSOLVERCONFIGPARSER_HPP

#include <string>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "solver/config/LinearSolverConfig.hpp"

namespace residuum {
	namespace solver {
		namespace parser {

			class LinearSolverConfigParser {
			public:

				static config::LinearSolverConfig::Type parseLinearSolverType(const std::string& str);

				static config::LinearSolverConfig parse(const YAML::Node& node);

			}; // class LinearSolverConfigParser

		} // namespace parser
	} // namespace solver
} // namespace residuum

#endif
