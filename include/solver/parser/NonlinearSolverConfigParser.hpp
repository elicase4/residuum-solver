#ifndef RESIDUUM_SOLVER_PARSER_NONLINEARSOLVERCONFIGPARSER_HPP
#define RESIDUUM_SOLVER_PARSER_NONLINEARSOLVERCONFIGPARSER_HPP

#include <string>
#include <stdexcept>

#include <yaml-cpp/yaml.h>

#include "solver/config/NonlinearSolverConfig.hpp"

namespace residuum {
	namespace solver {
		namespace parser {

			class NonlinearSolverConfigParser {
			public:

				static config::NonlinearSolverConfig::Type parseNonlinearSolverType(const std::string& str);

				static config::NonlinearSolverConfig parse(const YAML::Node& node);

			}; // class NonlinearSolverConfigParser

		} // namespace parser
	} // namespace solver
} // namespace residuum

#endif
