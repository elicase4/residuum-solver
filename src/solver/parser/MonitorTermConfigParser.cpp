#include "solver/parser/MonitorTermConfigParser.hpp"

#include "io/YAMLReader.hpp"

residuum::solver::config::MonitorTermConfig residuum::solver::parser::MonitorTermConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	residuum::solver::config::MonitorTermConfig term;

	term.boundary = YAMLReader::required<Int>(node, "boundary");
	term.coefficient = YAMLReader::optional<Real>(node, "coefficient", Real(1.0));

	const YAML::Node& weightsNode = node["weights"];
	if (weightsNode) {
		for (const auto& w : weightsNode) {
			term.weights.push_back(w.as<Real>());
		}
	}

	return term;

}
