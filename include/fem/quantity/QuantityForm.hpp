#ifndef RESIDUUM_FEM_QUANTITY_QUANTITYFORM_HPP
#define RESIDUUM_FEM_QUANTITY_QUANTITYFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace residuum {
	namespace fem {
		namespace quantity {

			template<typename FormT, typename QuadraturePointT>
			concept QuantityForm = requires (const FormT f, const QuadraturePointT& qp, const Real* Ue, Real* out) {
				{ FormT::NumComponents } -> std::convertible_to<Index>;
				{ f.computeElementLevelValue(qp, Ue, out) } -> std::same_as<void>;
			}; // concept QuantityForm

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#endif
