#ifndef RESIDUUM_UTILS_LOGGING_TIMESTEPPER_CONSOLELOGGER_HPP
#define RESIDUUM_UTILS_LOGGING_TIMESTEPPER_CONSOLELOGGER_HPP

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
			namespace timestepper {

				struct ConsoleLogger {

					// config
					std::string equationName;
					std::string timestepperName;
					Real t0;
					Real tf;
					Index interval;

					// consoleEnabled=false with an empty textFilePath prints nothing; use NullLogger directly for that case instead
					explicit ConsoleLogger(std::string equationNameIn, std::string timestepperNameIn, Real t0In, Real tfIn, Index reportInterval = 1, bool consoleEnabled = true, const std::string& textFilePath = "", const std::string& csvFilePath = "") :
						equationName(std::move(equationNameIn)), timestepperName(std::move(timestepperNameIn)), t0(t0In), tf(tfIn), interval(reportInterval) {

						// see solver::ConsoleLogger's identical note: teeBuf_/textFile_/out_ are
						// heap-allocated so a move of this struct (temporary -> variant) doesn't
						// dangle the pointer out_ holds into teeBuf_.
						teeBuf_ = std::make_unique<TeeStreamBuf>();

						if (consoleEnabled) teeBuf_->addTarget(std::cout.rdbuf());

						if (!textFilePath.empty()) {
							textFile_ = std::make_unique<std::ofstream>(textFilePath, std::ios::app);
							if (!textFile_->is_open()) {
								throw std::runtime_error("timestepper::ConsoleLogger: could not open '" + textFilePath + "' for writing");
							}
							teeBuf_->addTarget(textFile_->rdbuf());
						}

						out_ = std::make_unique<std::ostream>(teeBuf_.get());

						csv_ = std::make_unique<CsvWriter>(csvFilePath, std::vector<std::string>{"step", "time", "dt", "attempts", "residual_norm", "elapsed_s"});

						startTime_ = std::chrono::steady_clock::now();

					}

					// called before a step's nested nonlinear/linear solves begin, so their output
					// appears under this step's context rather than ahead of it
					void startStep(Index step, Real time) const {

						if (printHeader_) {
							printBanner();
							printHeader_ = false;
						}

						if (interval == 0) return;
						if ((step > 1) && ((step % interval) != 0)) return;

						*out_ << tag() << "step " << step << " starting: t=" << std::scientific << std::setprecision(4) << time << "\n";

					}

					void log(Index step, Real time, Real dt, Index attempts, Real residualNorm) const {

						const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime_).count();

						// CSV row: every accepted step, independent of any future console throttle
						if (csv_->enabled()) {
							csv_->writeRow({static_cast<Real>(step), time, dt, static_cast<Real>(attempts), residualNorm, static_cast<Real>(elapsed)});
						}

						lastStep_ = step;
						lastTime_ = time;
						totalAttempts_ += attempts;

						// console/text output throttle; the CSV row above is always written regardless
						if (interval == 0) return;
						if ((step > 1) && ((step % interval) != 0)) return;

						*out_ << tag() << "step " << step << " accepted: t=" << std::scientific << std::setprecision(4) << time;
						*out_ << ", dt=" << dt << ", " << attempts << (attempts == 1 ? " attempt" : " attempts") << ", ";
						*out_ << std::fixed << std::setprecision(3) << elapsed << "s\n";

					}

					void summary(bool converged) const {

						const double elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime_).count();

						*out_ << tag() << (converged ? "finished" : "did not reach tf") << ": " << lastStep_ << " steps, t=" << std::scientific << std::setprecision(4) << lastTime_ << " (" << totalAttempts_ << " total attempts), " << std::fixed << std::setprecision(3) << elapsed << "s\n";

					}

				private:

					std::string tag() const { return "  [TS:" + timestepperName + "] "; }

					void printBanner() const {

						*out_ << "\n";
						*out_ << tag() << "t0=" << t0 << "  tf=" << tf << "\n";

					}

					std::unique_ptr<TeeStreamBuf> teeBuf_;
					std::unique_ptr<std::ofstream> textFile_;
					std::unique_ptr<std::ostream> out_;
					std::unique_ptr<CsvWriter> csv_;

					mutable bool printHeader_ = true;
					mutable Index lastStep_ = 0;
					mutable Real lastTime_ = Real(0);
					mutable Index totalAttempts_ = 0;
					std::chrono::steady_clock::time_point startTime_;

				}; // struct ConsoleLogger

			} // namespace timestepper
		} // namespace logging
	} // namespace utils
} // namespace residuum

#endif
