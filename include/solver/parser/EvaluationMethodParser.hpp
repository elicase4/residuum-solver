#ifndef RESIDUUM_SOLVER_PARSER_EVALUATIONMETHODPARSER_HPP
#define RESIDUUM_SOLVER_PARSER_EVALUATIONMETHODPARSER_HPP

#include <string>

#include "solver/config/EvaluationMethod.hpp"

namespace residuum {
	namespace solver {
		namespace parser {

			class EvaluationMethodParser {
			public:

				static config::EvaluationMethod parse(const std::string& str);

			}; // class EvaluationMethodParser

		} // namespace parser
	} // namespace solver
} // namespace residuum

#endif
