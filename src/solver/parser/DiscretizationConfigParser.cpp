#include "solver/parser/DiscretizationConfigParser.hpp"
#include "io/YAMLReader.hpp"

namespace {

	residuum::fem::dof::DOFOrdering parseDOFOrdering(const std::string& str) {

		if (str == "interleaved") {
			return residuum::fem::dof::DOFOrdering::Interleaved;
		}

		if (str == "block") {
			return residuum::fem::dof::DOFOrdering::Block;
		}

		throw std::runtime_error("Unknown dof_ordering type: " + str);

	}

} // namespace

residuum::solver::config::DiscretizationConfig residuum::solver::parser::DiscretizationConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	config::DiscretizationConfig cfg;

	const YAML::Node& quadrature = node["quadrature"];
	if (!quadrature) {
		throw std::runtime_error("DiscretizationConfigParser: missing required 'quadrature' section");
	}
	cfg.quadrature.xi = YAMLReader::optional<Index>(quadrature, "xi", 2);
	cfg.quadrature.eta = YAMLReader::optional<Index>(quadrature, "eta", 2);
	cfg.quadrature.zeta = YAMLReader::optional<Index>(quadrature, "zeta", 2);

	const YAML::Node& dofOrdering = node["dof_ordering"];
	cfg.dofOrdering = dofOrdering ? parseDOFOrdering(YAMLReader::required<std::string>(dofOrdering, "type")) : fem::dof::DOFOrdering::Interleaved;

	return cfg;

}
