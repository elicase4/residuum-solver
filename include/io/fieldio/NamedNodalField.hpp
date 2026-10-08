#ifndef RESIDUUM_IO_FIELDIO_NAMEDNODALFIELD_HPP
#define RESIDUUM_IO_FIELDIO_NAMEDNODALFIELD_HPP

#include <string>

#include "core/Types.hpp"

namespace residuum {
	namespace io {
		namespace fieldio {

			// an already-computed, node-major point-data array for VTU/VTK output, additional to the
			// solved DOF field; numComponents == 1 writes as a scalar, 2 or 3 as a vector (see
			// VTUWriter::writeVector / VTKWriter::writeVector for the current component-count ceiling).
			// name should embed its unit (e.g. "Q[W/m^2]") -- neither VTK nor a ParaView-rendered VTU
			// surfaces a separate units attribute anywhere, so the array name is the only thing a
			// viewer actually displays
			struct NamedNodalField {

				std::string name;
				Index numComponents;
				const Real* data;

			}; // struct NamedNodalField

		} // namespace fieldio
	} // namespace io
} // namespace residuum

#endif
