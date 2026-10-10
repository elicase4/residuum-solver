#ifndef RESIDUUM_SOLVER_PARSER_REDUCTIONPARSER_HPP
#define RESIDUUM_SOLVER_PARSER_REDUCTIONPARSER_HPP

#include <string>

#include "fem/quantity/Reduction.hpp"

namespace residuum {
	namespace solver {
		namespace parser {

			class ReductionParser {
			public:

				static fem::quantity::Reduction parse(const std::string& str);

			}; // class ReductionParser

		} // namespace parser
	} // namespace solver
} // namespace residuum

#endif
