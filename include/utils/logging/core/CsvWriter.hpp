#ifndef RESIDUUM_UTILS_LOGGING_CORE_CSVWRITER_HPP
#define RESIDUUM_UTILS_LOGGING_CORE_CSVWRITER_HPP

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/Types.hpp"

namespace residuum {
	namespace utils {
		namespace logging {

			class CsvWriter {
			public:

				CsvWriter(const std::string& path, std::vector<std::string> columnNames) : numColumns_(columnNames.size()) {

					if (path.empty()) return;

					file_.open(path);
					if (!file_.is_open()) {
						throw std::runtime_error("CsvWriter: could not open '" + path + "' for writing");
					}

					for (Index i = 0; i < columnNames.size(); ++i) {
						file_ << columnNames[i];
						if (i + 1 < columnNames.size()) file_ << ",";
					}
					file_ << "\n";
					file_.flush();

				}

				bool enabled() const { return file_.is_open(); }

				void writeRow(const std::vector<Real>& values) const {

					if (!enabled()) return;

					if (values.size() != numColumns_) {
						throw std::runtime_error("CsvWriter::writeRow: expected " + std::to_string(numColumns_) + " columns, got " + std::to_string(values.size()));
					}

					for (Index i = 0; i < values.size(); ++i) {
						file_ << values[i];
						if (i + 1 < values.size()) file_ << ",";
					}
					file_ << "\n";
					file_.flush();

				}

			private:

				mutable std::ofstream file_;
				std::size_t numColumns_;

			}; // class CsvWriter

		} // namespace logging
	} // namespace utils
} // namespace residuum

#endif
