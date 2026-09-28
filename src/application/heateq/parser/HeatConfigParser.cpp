#include "application/heateq/parser/HeatConfigParser.hpp"
#include "application/heateq/parser/BoundaryConditionConfigParser.hpp"
#include "application/heateq/parser/ConductivityConfigParser.hpp"
#include "application/heateq/parser/InitialConditionConfigParser.hpp"
#include "application/heateq/parser/MonitorConfigParser.hpp"
#include "application/heateq/parser/ScalarMaterialPropertyConfigParser.hpp"
#include "application/heateq/parser/SourceConfigParser.hpp"
#include "application/heateq/parser/SpecificHeatConfigParser.hpp"

#include "solver/parser/BackendConfigParser.hpp"
#include "solver/parser/DiscretizationConfigParser.hpp"
#include "solver/parser/MeshConfigParser.hpp"
#include "solver/parser/SolverConfigParser.hpp"
#include "solver/parser/OutputConfigParser.hpp"
#include "solver/parser/LoggingConfigParser.hpp"

#include "io/YAMLReader.hpp"

residuum::application::heateq::config::HeatConfig residuum::application::heateq::parser::HeatConfigParser::read(const std::string& filename) {

	using io::YAMLReader;

	const YAML::Node root = YAMLReader::loadFile(filename);

	residuum::application::heateq::config::HeatConfig cfg;

	const YAML::Node& mesh = root["mesh"];
	if (!mesh) {
		throw std::runtime_error("HeatConfigReader: missing required 'mesh' section in " + filename);
	}
	cfg.mesh = residuum::solver::parser::MeshConfigParser::parse(mesh);

	const YAML::Node& disc = root["discretization"];
	if (!disc) {
		throw std::runtime_error("HeatConfigReader: missing required 'discretization' section in " + filename);
	}
	cfg.discretization = residuum::solver::parser::DiscretizationConfigParser::parse(disc);

	const YAML::Node& solver = root["solver"];
	if (!solver) {
		throw std::runtime_error("HeatConfigParser: missing required 'solver' section in " + filename);
	}
	cfg.solver = residuum::solver::parser::SolverConfigParser::parse(solver);
	
	// optional; absent means no field-snapshot output at all
	const YAML::Node& output = root["output"];
	if (output) {
		cfg.output = residuum::solver::parser::OutputConfigParser::parse(output);
	}

	cfg.logging = residuum::solver::parser::LoggingConfigParser::parse(root["logging"]);

	const YAML::Node& initialCondition = root["initial_condition"];
	if (!initialCondition) {
		throw std::runtime_error("HeatConfigReader: missing required 'initial_condition' section in " + filename);
	}
	cfg.initialCondition = residuum::application::heateq::parser::InitialConditionConfigParser::parse(initialCondition);

	const YAML::Node& bcs = root["boundary_conditions"];
	if (!bcs) {
		throw std::runtime_error("HeatConfigReader: missing required 'boundary_conditions' section in " + filename);
	}
	for (auto& bc : bcs) {
		cfg.boundaryConditions.push_back(residuum::application::heateq::parser::BoundaryConditionConfigParser::parse(bc));
	}

	const YAML::Node& materials = root["materials"];
	if (!materials) {
		throw std::runtime_error("HeatConfigReader: missing required 'materials' section in " + filename);
	}

	const YAML::Node& cond = materials["conductivity"];
	if (!cond) {
		throw std::runtime_error("HeatConfigReader: missing required 'materials.conductivity' section in " + filename);
	}
	cfg.conductivity = residuum::application::heateq::parser::ConductivityConfigParser::parse(cond);

	// density/specific_heat are used for a transient driver;
	const YAML::Node& density = materials["density"];
	if (density) {
		cfg.density = residuum::application::heateq::parser::ScalarMaterialPropertyConfigParser::parse(density);
	}

	const YAML::Node& specificHeat = materials["specific_heat"];
	if (specificHeat) {
		cfg.specificHeat = residuum::application::heateq::parser::SpecificHeatConfigParser::parse(specificHeat);
	}

	const YAML::Node& physics = root["physics"];
	if (!physics) {
		throw std::runtime_error("HeatConfigReader: missing required 'physics' section in " + filename);
	}

	const YAML::Node& models = physics["models"];
	if (!models) {
		throw std::runtime_error("HeatConfigReader: missing required 'physics.models' section in " + filename);
	}

	const YAML::Node& src = models["source"];
	if (!src) {
		throw std::runtime_error("HeatConfigReader: missing required 'physics.models.source' section in " + filename);
	}
	cfg.source = residuum::application::heateq::parser::SourceConfigParser::parse(src);

	cfg.backend = residuum::solver::parser::BackendConfigParser::parse(root);

	const YAML::Node& monitors = root["monitors"];
	if (monitors) {
		for (auto& monitor : monitors) {
			cfg.monitors.push_back(residuum::application::heateq::parser::MonitorConfigParser::parse(monitor));
		}
	}

	return cfg;

}
