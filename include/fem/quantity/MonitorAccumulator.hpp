#ifndef RESIDUUM_FEM_QUANTITY_MONITORACCUMULATOR_HPP
#define RESIDUUM_FEM_QUANTITY_MONITORACCUMULATOR_HPP

#include <functional>
#include <vector>

#include "core/Types.hpp"

#include "fem/quantity/BoundaryQuantityCombination.hpp"
#include "fem/quantity/BoundaryQuantityRegistry.hpp"

namespace residuum {
	namespace fem {
		namespace quantity {

			// registry + combination bookkeeping shared by every MonitorGroup implementation, regardless of how the registry actually gets filled
			template<Index NumComponents>
			class MonitorAccumulator {
			public:

				void registerTag(Int tag) { registry_.registerTag(tag); }

				Index addCombination() { combinations_.emplace_back(); return static_cast<Index>(combinations_.size() - 1); }

				void addTerm(Index combinationIndex, Int tag, Real coefficient) { combinations_[combinationIndex].addTerm(tag, coefficient); }

				BoundaryQuantityRegistry<NumComponents>& registry() { return registry_; }

				void forEachCombination(const std::function<void(Index, const Real*)>& callback) const {
					Real value[NumComponents];
					for (Index i = 0; i < static_cast<Index>(combinations_.size()); ++i) {
						combinations_[i].evaluate(registry_, value);
						callback(i, value);
					}
				}

			private:

				BoundaryQuantityRegistry<NumComponents> registry_;
				std::vector<BoundaryQuantityCombination<NumComponents>> combinations_;

			}; // class MonitorAccumulator

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#endif
