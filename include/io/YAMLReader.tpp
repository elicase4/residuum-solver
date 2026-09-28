namespace residuum::io {

	template<typename T>
	T YAMLReader::required(const YAML::Node& node, const std::string& key) {

		if (!node[key]) {
			throw std::runtime_error("Missing YAML key: " + key);
		}

		return node[key].as<T>();
	}

	template<typename T>
	T YAMLReader::optional(const YAML::Node& node, const std::string& key, const T& defaultValue) {

		if (!node[key]) {
			return defaultValue;
		}

		return node[key].as<T>();
	}

} // namespace residuum::io
