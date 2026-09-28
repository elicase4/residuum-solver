#include "application/mesh/MeshConfigParser.hpp"

#include "io/YAMLReader.hpp"
#include "solver/parser/LoggingConfigParser.hpp"

residuum::application::mesh::MeshConfig::Type residuum::application::mesh::MeshConfigParser::parseMeshType(const std::string& str) {

	if (str == "block2d") {
		return residuum::application::mesh::MeshConfig::Type::Block2D;
	}

	if (str == "block3d") {
		return residuum::application::mesh::MeshConfig::Type::Block3D;
	}

	if (str == "gmsh") {
		return residuum::application::mesh::MeshConfig::Type::Gmsh;
	}

	throw std::runtime_error("Unknown mesh type: " + str);

}

residuum::application::mesh::MeshConfig residuum::application::mesh::MeshConfigParser::read(const std::string& filename) {

	using io::YAMLReader;

	const YAML::Node root = YAMLReader::loadFile(filename);

	residuum::application::mesh::MeshConfig cfg;

	const YAML::Node& meshNode = root["mesh"];
	if (!meshNode) {
		throw std::runtime_error("MeshConfigReader: missing required 'mesh' section in " + filename);
	}

	cfg.type = residuum::application::mesh::MeshConfigParser::parseMeshType(YAMLReader::required<std::string>(meshNode, "type"));

	if (cfg.type == residuum::application::mesh::MeshConfig::Type::Block2D) {

		cfg.block2D.nx   = YAMLReader::required<Index>(meshNode, "nx");
		cfg.block2D.ny   = YAMLReader::required<Index>(meshNode, "ny");
		cfg.block2D.Px   = YAMLReader::required<Index>(meshNode, "px");
		cfg.block2D.Py   = YAMLReader::required<Index>(meshNode, "py");
		cfg.block2D.xmin = YAMLReader::required<Real>(meshNode, "xmin");
		cfg.block2D.xmax = YAMLReader::required<Real>(meshNode, "xmax");
		cfg.block2D.ymin = YAMLReader::required<Real>(meshNode, "ymin");
		cfg.block2D.ymax = YAMLReader::required<Real>(meshNode, "ymax");

	} else if (cfg.type == MeshConfig::Type::Block3D) {

		cfg.block3D.nx   = YAMLReader::required<Index>(meshNode, "nx");
		cfg.block3D.ny   = YAMLReader::required<Index>(meshNode, "ny");
		cfg.block3D.nz   = YAMLReader::required<Index>(meshNode, "nz");
		cfg.block3D.Px   = YAMLReader::required<Index>(meshNode, "px");
		cfg.block3D.Py   = YAMLReader::required<Index>(meshNode, "py");
		cfg.block3D.Pz   = YAMLReader::required<Index>(meshNode, "pz");
		cfg.block3D.xmin = YAMLReader::required<Real>(meshNode, "xmin");
		cfg.block3D.xmax = YAMLReader::required<Real>(meshNode, "xmax");
		cfg.block3D.ymin = YAMLReader::required<Real>(meshNode, "ymin");
		cfg.block3D.ymax = YAMLReader::required<Real>(meshNode, "ymax");
		cfg.block3D.zmin = YAMLReader::required<Real>(meshNode, "zmin");
		cfg.block3D.zmax = YAMLReader::required<Real>(meshNode, "zmax");

	} else if (cfg.type == residuum::application::mesh::MeshConfig::Type::Gmsh) {

		cfg.inputFile = YAMLReader::required<std::string>(meshNode, "file");
	
	}

	const YAML::Node& outNode = root["output"];
	if (!outNode) {
		throw std::runtime_error("MeshConfigReader: missing required 'output' section in " + filename);
	}
	
	cfg.outputFile = YAMLReader::required<std::string>(outNode, "output_file");

	// 'logging:' is optional; absent means console logging with no file mirror
	cfg.logging = residuum::solver::parser::LoggingConfigParser::parse(root["logging"]);

	return cfg;

}
