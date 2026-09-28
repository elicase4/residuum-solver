#ifndef RESIDUUM_IO_FIELDIO_TIMESERIESVALUESOURCE_HPP
#define RESIDUUM_IO_FIELDIO_TIMESERIESVALUESOURCE_HPP

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <optional>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>

#include "core/Types.hpp"
#include "io/fieldio/FieldIO.hpp"
#include "mesh/Mesh.hpp"

namespace residuum {
	namespace io {
		namespace fieldio {

			template<Index numDOFs>
			class TimeSeriesValueSource {
			public:

				enum class Interpolation { Exact, Interpolate };

				static constexpr Index NumComponents = numDOFs;

				TimeSeriesValueSource(const mesh::Mesh& mesh, const std::string& pattern, Interpolation interpolation = Interpolation::Exact, Real timeTolerance = Real(1e-9)) : mesh_(mesh), pattern_(pattern), interpolation_(interpolation), timeTolerance_(timeTolerance) {}

				void value(Index nodeID, Index step, Real time, Real* out) {

					if (interpolation_ == Interpolation::Interpolate) {
						throw std::runtime_error("TimeSeriesValueSource: Interpolation::Interpolate is not yet implemented");
					}

					loadStepIfNeeded(step);

					const Real scale = std::max(Real(1), std::max(std::abs(currentFrame_.time), std::abs(time)));

					if (std::abs(currentFrame_.time - time) > timeTolerance_ * scale) {
						throw std::runtime_error("TimeSeriesValueSource: frame for step " + std::to_string(step) + " recorded time " + std::to_string(currentFrame_.time) + ", caller expected " + std::to_string(time) + " -- file series does not match the live solve's stepping");
					}

					for (Index c = 0; c < numDOFs; ++c) {
						out[c] = currentFrame_.values[nodeID * numDOFs + c];
					}

				}

				// time recorded in the currently-loaded frame's header
				Real currentFrameTime() const { return currentFrame_.time; }

			private:

				void loadStepIfNeeded(Index step) {

					if (lastStep_.has_value() && *lastStep_ == step) return;

					currentFrame_ = FieldIO::readBinary<numDOFs>(mesh_, formatStepFilename(pattern_, step));
					lastStep_ = step;

				}

				static std::string formatStepFilename(const std::string& pattern, Index step) {

					static const std::regex placeholder(R"(\{step:0(\d+)d\})");
					std::smatch match;

					if (!std::regex_search(pattern, match, placeholder)) {
						throw std::runtime_error("TimeSeriesValueSource: pattern '" + pattern + "' must contain a {step:0Nd} placeholder");
					}

					const int width = std::stoi(match[1].str());

					std::ostringstream oss;
					oss << std::setw(width) << std::setfill('0') << step;

					return std::regex_replace(pattern, placeholder, oss.str());

				}

				const mesh::Mesh& mesh_;
				std::string pattern_;
				Interpolation interpolation_;
				Real timeTolerance_;
				std::optional<Index> lastStep_;
				NodalFieldSnapshot currentFrame_;

			}; // class TimeSeriesValueSource

		} // namespace fieldio
	} // namespace io
} // namespace residuum

#endif
