#ifndef RESIDUUM_IO_FIELDIO_NAMEDNODALFIELD_HPP
#define RESIDUUM_IO_FIELDIO_NAMEDNODALFIELD_HPP

#include <string>

#include "core/Types.hpp"

namespace residuum {
	namespace io {
		namespace fieldio {

			// name must embed its own unit (e.g. "Q[W/m^2]"); neither VTK nor a ParaView-rendered VTU surfaces a separate units attribute, so the array name is the only thing a viewer displays
			struct NamedNodalField {

				std::string name;
				Index numComponents;
				const Real* data;

			}; // struct NamedNodalField

		} // namespace fieldio
	} // namespace io
} // namespace residuum

#endif
