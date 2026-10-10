#include "solver/parser/ReductionParser.hpp"

#include <stdexcept>

residuum::fem::quantity::Reduction residuum::solver::parser::ReductionParser::parse(const std::string& str) {

	if (str == "integral") {
		return residuum::fem::quantity::Reduction::Integral;
	} else if (str == "average") {
		return residuum::fem::quantity::Reduction::Average;
	} else {
		throw std::runtime_error("ReductionParser: unknown reduction '" + str + "'. Valid options: integral, average");
	}

}
