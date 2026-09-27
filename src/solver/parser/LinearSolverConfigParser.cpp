#include "solver/parser/LinearSolverConfigParser.hpp"
#include "io/YAMLReader.hpp"

namespace {

	residuum::solver::config::LinearSolverConfig::OperatorType parseOperatorType(const std::string& str) {

		if (str == "csr")
			return residuum::solver::config::LinearSolverConfig::OperatorType::CSR;
		if (str == "fem")
			return residuum::solver::config::LinearSolverConfig::OperatorType::FEM;

		throw std::runtime_error("Unknown linear solver operator type: " + str + ". Valid options: csr, fem");

	}

	residuum::solver::config::PreconditionerConfig::Type parsePreconditionerType(const std::string& str) {

		if (str == "identity")
			return residuum::solver::config::PreconditionerConfig::Type::Identity;
		if (str == "jacobi")
			return residuum::solver::config::PreconditionerConfig::Type::Jacobi;

		throw std::runtime_error("Unknown preconditioner type: " + str + ". Valid options: identity, jacobi");

	}

} // namespace

residuum::solver::config::LinearSolverConfig::Type residuum::solver::parser::LinearSolverConfigParser::parseLinearSolverType(const std::string& str) {

	if (str == "cg")
		return residuum::solver::config::LinearSolverConfig::Type::CG;
	if (str == "gmres")
		return residuum::solver::config::LinearSolverConfig::Type::GMRES;
	if (str == "bicgstab")
		return residuum::solver::config::LinearSolverConfig::Type::BiCGSTAB;
	if (str == "lu")
		return residuum::solver::config::LinearSolverConfig::Type::LU;

	throw std::runtime_error("Unknown linear solver type: " + str + ". Valid options: cg, gmres, bicgstab, lu");

}

residuum::solver::config::LinearSolverConfig residuum::solver::parser::LinearSolverConfigParser::parse(const YAML::Node& node) {

	using io::YAMLReader;

	residuum::solver::config::LinearSolverConfig cfg;

	cfg.type = parseLinearSolverType(YAMLReader::required<std::string>(node, "type"));
	cfg.tolerance = YAMLReader::optional<Real>(node, "tolerance", 1e-10);
	cfg.maxIterations = YAMLReader::optional<Index>(node, "max_iterations", 1000);

	switch (cfg.type) {

		case config::LinearSolverConfig::Type::CG:
			cfg.params = config::CGParams{};
			break;

		case config::LinearSolverConfig::Type::GMRES: {
			config::GMRESParams p;
			p.krylovDim = YAMLReader::optional<Index>(node, "krylov_dim", 50);
			cfg.params = p;
			break;
		}

		case config::LinearSolverConfig::Type::BiCGSTAB:
			cfg.params = config::BiCGSTABParams{};
			break;

		case config::LinearSolverConfig::Type::LU:
			cfg.params = config::LUParams{};
			break;

	}

	const YAML::Node& op = node["operator"];
	cfg.operatorType = op ? parseOperatorType(YAMLReader::required<std::string>(op, "type")) : config::LinearSolverConfig::OperatorType::CSR;

	const YAML::Node& precond = node["preconditioner"];
	cfg.preconditioner.type = precond ? parsePreconditionerType(YAMLReader::required<std::string>(precond, "type")) : config::PreconditionerConfig::Type::Identity;

	switch (cfg.preconditioner.type) {

		case config::PreconditionerConfig::Type::Identity:
			cfg.preconditioner.params = config::IdentityPreconditionerParams{};
			break;

		case config::PreconditionerConfig::Type::Jacobi:
			cfg.preconditioner.params = config::JacobiPreconditionerParams{};
			break;

	}

	return cfg;

}
