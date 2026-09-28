#ifndef RESIDUUM_FEM_BOUNDARY_ESSENTIALBOUNDARYREGISTRY_HPP
#define RESIDUUM_FEM_BOUNDARY_ESSENTIALBOUNDARYREGISTRY_HPP

#include <memory>
#include <unordered_map>
#include <vector>

#include "BoundaryCondition.hpp"
#include "EssentialBoundaryCondition.hpp"

namespace residuum {
	namespace fem {
		namespace boundary {

			class EssentialBoundaryRegistry {
			public:

				template<BoundaryFunction FunctionT>
				void registerBC(std::shared_ptr<BoundaryCondition<FunctionT>> bc) {
					entries_[bc->tag].push_back(std::make_unique<EssentialBoundaryCondition<FunctionT>>(std::move(bc)));
				}

				const std::vector<std::unique_ptr<EssentialBoundaryConditionBase>>* getEntries(Int tag) const {
					
					auto it = entries_.find(tag);

					if (it == entries_.end())
						return nullptr;

					return &it->second;
				}

				bool hasAny(Int tag) const {
					return entries_.contains(tag);
				}

				bool isEssential(Int tag, Index component) const {
					
					const auto* entries = getEntries(tag);
					
					if (!entries) {
						return false;
					}

					for (auto& entry : *entries) {
						if (entry->componentType(component) == BCCategory::Essential){
							return true;
						}
					}
				
					return false;

				}

			private:

				std::unordered_map<Int, std::vector<std::unique_ptr<EssentialBoundaryConditionBase>>> entries_;

			}; // class EssentialBoundaryRegistry

		} // namespace boundary
	} // namespace fem
} // namepspace residuum

#endif
