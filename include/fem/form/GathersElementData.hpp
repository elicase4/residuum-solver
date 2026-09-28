#ifndef RESIDUUM_FEM_FORM_GATHERSELEMENTDATA_HPP
#define RESIDUUM_FEM_FORM_GATHERSELEMENTDATA_HPP

#include <concepts>

#include "core/Types.hpp"

namespace residuum {
	namespace fem {
		namespace form {

			template<typename FormT>
			concept GathersElementData = requires (const FormT f, const Index* nodeIDs, Index nodesPerElement) {
				{ f.gatherElementData(nodeIDs, nodesPerElement) } -> std::same_as<void>;
			}; // concept GathersElementData

		} // namespace form
	} // namespace fem
} // namespace residuum

#endif
