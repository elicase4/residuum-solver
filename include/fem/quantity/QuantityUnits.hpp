#ifndef RESIDUUM_FEM_QUANTITY_QUANTITYUNITS_HPP
#define RESIDUUM_FEM_QUANTITY_QUANTITYUNITS_HPP

#include <string>

#include "core/Types.hpp"

#include "fem/quantity/Reduction.hpp"

namespace residuum::fem::quantity {

	template<typename FormT>
	std::string unitFor(Reduction mode, Index boundaryDim) {

		std::string unit = FormT::BaseUnit;
		if (mode == Reduction::Average) {
			unit += "/m^" + std::to_string(boundaryDim);
		}

		return unit;

	}

} // namespace residuum::fem::quantity

#endif
