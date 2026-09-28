#ifndef RESIDUUM_IO_YAMLREADER_HPP
#define RESIDUUM_IO_YAMLREADER_HPP

#include <string>
#include <yaml-cpp/yaml.h>

namespace residuum {
	namespace io {

		class YAMLReader {
		public:
			
			static YAML::Node loadFile(const std::string& filename);

			template<typename T>
			static T required(const YAML::Node& node, const std::string& key);

			template<typename T>
			static T optional(const YAML::Node& node, const std::string& key, const T& defaultValue);

		}; // class YAMLReader

	} // namespace io
} // namespace residuum

#include "YAMLReader.tpp"

#endif
