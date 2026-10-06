#include "solver/parser/EvaluationMethodParser.hpp"

#include <stdexcept>

residuum::solver::config::EvaluationMethod residuum::solver::parser::EvaluationMethodParser::parse(const std::string& str) {

	if (str == "strong") {
		return residuum::solver::config::EvaluationMethod::Strong;
	} else if (str == "weak") {
		return residuum::solver::config::EvaluationMethod::Weak;
	} else {
		throw std::runtime_error("EvaluationMethodParser: unknown method '" + str + "'. Valid options: strong, weak");
	}

}
