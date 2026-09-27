#ifndef RESIDUUM_BENCHMARK_CORE_CSVREPORTER_HPP
#define RESIDUUM_BENCHMARK_CORE_CSVREPORTER_HPP

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "benchmark/core/Result.hpp"

namespace residuum {
	namespace benchmark {
		namespace core {

			// appends one row per result to a plain CSV file (creating it with a header
			// on first write), one column per axis plus timing/FLOPs columns -- kept
			// simple on purpose so it opens directly in ParaView's spreadsheet/plot views
			class CSVReporter {
			public:

				explicit CSVReporter(const std::string& path) : path_(path) {}

				// all results appended to the same file must share the same axis names
				// (in the same order), since a CSV has one fixed column set
				void append(const Result& result) const {

					namespace fs = std::filesystem;

					const bool needsHeader = !fs::exists(path_) || fs::is_empty(path_);

					if (!needsHeader) checkAxesMatchExistingHeader(result);

					std::ofstream out(path_, std::ios::app);
					if (!out) throw std::runtime_error("CSVReporter: could not open " + path_ + " for writing");

					if (needsHeader) writeHeader(out, result);
					writeRow(out, result);

				}

			private:

				std::string path_;

				static void writeHeader(std::ofstream& out, const Result& result) {

					out << "kernel,equation";
					for (const auto& [axisName, axisValue] : result.axisValues) {
						(void)axisValue;
						out << "," << axisName;
					}
					out << ",warmup_iters,measure_iters,min_seconds,mean_seconds,median_seconds,estimated_flops,gflops_per_second\n";

				}

				static void writeRow(std::ofstream& out, const Result& result) {

					out << result.kernelName << "," << result.equation;
					for (const auto& [axisName, axisValue] : result.axisValues) {
						(void)axisName;
						out << "," << axisValue;
					}

					out << "," << result.warmupIters << "," << result.measureIters << "," << result.minSeconds << "," << result.meanSeconds << "," << result.medianSeconds;

					if (result.estimatedFlops) out << "," << *result.estimatedFlops << "," << *result.gflopsPerSecond();
					else out << ",,";

					out << "\n";

				}

				void checkAxesMatchExistingHeader(const Result& result) const {

					std::ifstream in(path_);
					std::string headerLine;
					std::getline(in, headerLine);

					std::vector<std::string> existingColumns;
					std::stringstream ss(headerLine);
					std::string column;
					while (std::getline(ss, column, ',')) existingColumns.push_back(column);

					// columns are: kernel, equation, <axes...>, warmup_iters, ...
					const Index numFixedTrailingColumns = 7;
					if (existingColumns.size() < 2 + numFixedTrailingColumns) return;

					const Index numExistingAxes = existingColumns.size() - 2 - numFixedTrailingColumns;

					if (numExistingAxes != result.axisValues.size()) {
						throw std::runtime_error("CSVReporter: " + path_ + " already has " + std::to_string(numExistingAxes) + " axis columns, but this result has " + std::to_string(result.axisValues.size()) + " -- use a separate file per kernel family");
					}

					for (Index i = 0; i < numExistingAxes; ++i) {
						if (existingColumns[2 + i] != result.axisValues[i].first) {
							throw std::runtime_error("CSVReporter: " + path_ + " axis column '" + existingColumns[2 + i] + "' does not match result's axis '" + result.axisValues[i].first + "'");
						}
					}

				}

			}; // class CSVReporter

		} // namespace core
	} // namespace benchmark
} // namespace residuum

#endif
