#ifndef RESIDUUM_IO_FIELDIO_NODALFILEVALUESOURCE_HPP
#define RESIDUUM_IO_FIELDIO_NODALFILEVALUESOURCE_HPP

#include <string>
#include <vector>

#include "core/Types.hpp"
#include "io/fieldio/FieldIO.hpp"
#include "mesh/Mesh.hpp"

namespace residuum {
	namespace io {
		namespace fieldio {

			template<Index numDOFs>
			class NodalFileValueSource {
			public:

				static constexpr Index NumComponents = numDOFs;

				NodalFileValueSource(const mesh::Mesh& mesh, const std::string& file) {
					const auto snapshot = FieldIO::readBinary<numDOFs>(mesh, file);
					values_ = snapshot.values;
					fileTime_ = snapshot.time;
				}

				void eval(Index nodeID, Real* out) const {
					for (Index c = 0; c < numDOFs; ++c) {
						out[c] = values_[nodeID * numDOFs + c];
					}
				}

				Real fileTime() const { return fileTime_; }

			private:

				std::vector<Real> values_;
				Real fileTime_;

			}; // class NodalFileValueSource

		} // namespace fieldio
	} // namespace io
} // namespace residuum

#endif
