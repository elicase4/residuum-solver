#ifndef RESIDUUM_UTILS_LOGGING_NONLINEAR_CONSOLELOGGER_HPP
#define RESIDUUM_UTILS_LOGGING_NONLINEAR_CONSOLELOGGER_HPP

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/Types.hpp"

#include "utils/logging/core/CsvWriter.hpp"
#include "utils/logging/core/TeeStreamBuf.hpp"

namespace residuum {
	namespace utils {
		namespace logging {
			namespace nonlinear {

				struct ConsoleLogger {

					// config
					std::string equationName;
					std::string solverName;
					std::vector<std::string> dofNames;
					Index interval;

					// consoleEnabled=false with an empty textFilePath prints nothing; use NullLogger directly for that case instead
					explicit ConsoleLogger(std::string equationNameIn, std::string solverNameIn, std::vector<std::string> dofNamesIn, Index reportInterval = 1, bool consoleEnabled = true, const std::string& textFilePath = "", const std::string& csvFilePath = "") :
						equationName(std::move(equationNameIn)), solverName(std::move(solverNameIn)), dofNames(std::move(dofNamesIn)), interval(reportInterval) {

						// see solver::ConsoleLogger's identical note: teeBuf_/textFile_/out_ are
						// heap-allocated so a move of this struct (temporary -> variant) doesn't
						// dangle the pointer out_ holds into teeBuf_.
						teeBuf_ = std::make_unique<TeeStreamBuf>();

						if (consoleEnabled) teeBuf_->addTarget(std::cout.rdbuf());

						if (!textFilePath.empty()) {
							textFile_ = std::make_unique<std::ofstream>(textFilePath, std::ios::app);
							if (!textFile_->is_open()) {
								throw std::runtime_error("nonlinear::ConsoleLogger: could not open '" + textFilePath + "' for writing");
							}
							teeBuf_->addTarget(textFile_->rdbuf());
						}

						out_ = std::make_unique<std::ostream>(teeBuf_.get());

						std::vector<std::string> csvColumns = {"outer_tick", "iter"};
						for (const auto& name : dofNames) csvColumns.push_back("rel[" + name + "]");
						csvColumns.push_back("elapsed_s");
						csv_ = std::make_unique<CsvWriter>(csvFilePath, std::move(csvColumns));

						startTime_ = std::chrono::steady_clock::now();

					}

					// marks the start of a fresh nonlinear solve (e.g. a new timestep) so rows can be
					// correlated with the outer context that produced them, and so the banner/column
					// header reprint to mark where this solve begins in the console/text output
					void reset() const {
						printHeader_ = true;
						++outerTick_;
					}

					void log(Index iter, const std::vector<Real>& residualRelPerDOF) const {

						const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime_).count();

						if (csv_->enabled()) {
							std::vector<Real> row = {static_cast<Real>(outerTick_), static_cast<Real>(iter)};
							for (const auto& v : residualRelPerDOF) row.push_back(v);
							row.push_back(static_cast<Real>(elapsed));
							csv_->writeRow(row);
						}

						lastIter_ = iter;
						lastRelPerDOF_ = residualRelPerDOF;

						// console/text output throttle; the CSV row above is always written regardless
						if (interval == 0) return;
						if ((iter > 0) && ((iter % interval) != 0)) return;

						if (printHeader_) {
							printBanner();
							printColumnHeader();
							printHeader_ = false;
						}

						*out_ << tag() << std::left << std::setw(6) << iter << "  ";
						*out_ << std::scientific << std::setprecision(4);
						for (const auto& v : residualRelPerDOF) *out_ << std::setw(14) << v << "  ";
						*out_ << std::fixed << std::setprecision(3) << elapsed << "s\n";

					}

					void summary(bool converged) const {

						const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime_).count();

						*out_ << tag() << (converged ? "converged" : "did not converge") << ": " << lastIter_ << " iters, ";

						for (Index i = 0; i < lastRelPerDOF_.size(); ++i) {
							*out_ << "rel[" << dofNames[i] << "] " << std::scientific << std::setprecision(4) << lastRelPerDOF_[i];
							if (i + 1 < lastRelPerDOF_.size()) *out_ << ", ";
						}

						*out_ << ", " << std::fixed << std::setprecision(3) << elapsed << "s\n";

					}

				private:

					std::string tag() const { return "    [NL:" + solverName + "] "; }

					void printBanner() const {

						const int width = 60;

						*out_ << "\n";
						*out_ << tag() << std::string(width, '=') << "\n";

					}

					void printColumnHeader() const {

						*out_ << tag() << std::left << std::setw(6) << "Iter" << "  ";
						for (const auto& name : dofNames) *out_ << std::setw(14) << ("Rel[" + name + "]") << "  ";
						*out_ << "Elapsed" << "\n";
						*out_ << tag() << std::string(60, '-') << "\n";

					}

					std::unique_ptr<TeeStreamBuf> teeBuf_;
					std::unique_ptr<std::ofstream> textFile_;
					std::unique_ptr<std::ostream> out_;
					std::unique_ptr<CsvWriter> csv_;

					mutable bool printHeader_ = true;
					mutable Index lastIter_ = 0;
					mutable std::vector<Real> lastRelPerDOF_;
					mutable Index outerTick_ = 0;
					std::chrono::steady_clock::time_point startTime_;

				}; // struct ConsoleLogger

			} // namespace nonlinear
		} // namespace logging
	} // namespace utils
} // namespace residuum

#endif
