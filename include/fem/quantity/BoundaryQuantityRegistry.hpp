#ifndef RESIDUUM_FEM_QUANTITY_BOUNDARYQUANTITYREGISTRY_HPP
#define RESIDUUM_FEM_QUANTITY_BOUNDARYQUANTITYREGISTRY_HPP

#include <array>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "core/Types.hpp"

namespace residuum {
	namespace fem {
		namespace quantity {

			// per-tag accumulator; NumComponents is a plain count so any evaluation strategy can share
			// this (a QuantityFormsT's TotalComponents, a problem's numDOFs, etc.)
			template<Index NumComponents>
			class BoundaryQuantityRegistry {
			public:

				struct Entry {
					std::array<Real, NumComponents> result{};
					Real measureSum = 0;
				}; // struct Entry

				// setup-time only
				void registerTag(Int tag) {
					entries_.try_emplace(tag);
				}

				bool hasTag(Int tag) const {
					return entries_.contains(tag);
				}

				// finalized result for one registered tag
				const Real* result(Int tag) const {
					auto it = entries_.find(tag);
					if (it == entries_.end()) {
						throw std::runtime_error("BoundaryQuantityRegistry::result: tag " + std::to_string(tag) + " was never registered");
					}
					return it->second.result.data();
				}

				std::unordered_map<Int, Entry>& entries() { return entries_; }
				const std::unordered_map<Int, Entry>& entries() const { return entries_; }

			private:

				std::unordered_map<Int, Entry> entries_;

			}; // class BoundaryQuantityRegistry

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#endif
