#include "io/YAMLReader.hpp"

YAML::Node residuum::io::YAMLReader::loadFile(const std::string& filename) {
	return YAML::LoadFile(filename);
}
