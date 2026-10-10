#ifndef RESIDUUM_FEM_QUANTITY_BOUNDARYVECTORPROJECTION_HPP
#define RESIDUUM_FEM_QUANTITY_BOUNDARYVECTORPROJECTION_HPP

#include <array>
#include <utility>
#include <vector>

#include "core/Types.hpp"

#include "fem/quantity/BoundaryQuantityRegistry.hpp"

namespace residuum {
	namespace fem {
		namespace quantity {

			// collapses a multi-component per-tag result (e.g. a force vector) to one scalar via per-component weights summed across tags, unlike BoundaryQuantityCombination's single uniform coefficient
			template<Index NumComponents>
			class BoundaryVectorProjection {
			public:

				void addTerm(Int tag, const std::array<Real, NumComponents>& weights) {
					terms_.emplace_back(tag, weights);
				}

				Real evaluate(const BoundaryQuantityRegistry<NumComponents>& registry) const {

					Real out = 0.0;

					for (const auto& [tag, weights] : terms_) {

						const Real* tagResult = registry.result(tag);

						for (Index c = 0; c < NumComponents; ++c) {
							out += weights[c] * tagResult[c];
						}

					}

					return out;

				}

			private:

				std::vector<std::pair<Int, std::array<Real, NumComponents>>> terms_;

			}; // class BoundaryVectorProjection

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#endif
