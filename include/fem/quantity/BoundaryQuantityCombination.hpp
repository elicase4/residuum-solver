#ifndef RESIDUUM_FEM_QUANTITY_BOUNDARYQUANTITYCOMBINATION_HPP
#define RESIDUUM_FEM_QUANTITY_BOUNDARYQUANTITYCOMBINATION_HPP

#include <utility>
#include <vector>

#include "core/Types.hpp"

#include "fem/quantity/BoundaryQuantityRegistry.hpp"

namespace residuum {
	namespace fem {
		namespace quantity {

			template<Index NumComponents>
			class BoundaryQuantityCombination {
			public:

				void addTerm(Int tag, Real coefficient) {
					terms_.emplace_back(tag, coefficient);
				}

				// out has size NumComponents
				void evaluate(const BoundaryQuantityRegistry<NumComponents>& registry, Real* out) const {

					for (Index c = 0; c < NumComponents; ++c) out[c] = 0.0;

					for (const auto& [tag, coefficient] : terms_) {

						const Real* tagResult = registry.result(tag);

						for (Index c = 0; c < NumComponents; ++c) {
							out[c] += coefficient * tagResult[c];
						}

					}

				}

			private:

				std::vector<std::pair<Int, Real>> terms_;

			}; // class BoundaryQuantityCombination

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#endif
