#ifndef RESIDUUM_FEM_QUANTITY_BOUNDARYQUANTITYREGISTRY_HPP
#define RESIDUUM_FEM_QUANTITY_BOUNDARYQUANTITYREGISTRY_HPP

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "config/Platform.hpp"
#include "core/Types.hpp"

namespace residuum {
	namespace fem {
		namespace quantity {

			// per-tag accumulator; NumComponents is a plain count so any evaluation strategy can share this (a QuantityFormsT's TotalComponents, a problem's numDOFs, etc.)
			template<Index NumComponents>
			class BoundaryQuantityRegistry {
			public:

				struct Entry {
					std::array<Real, NumComponents> result{};
					Real measureSum = 0;
				}; // struct Entry

				void registerTag(Int tag){
					
					auto emplace_result = entries_.try_emplace(tag);
					
					if (static_cast<Index>(tag) >= tagToEntry_.size()){
						tagToEntry_.resize(tag + 1, nullptr);
					}
					
					tagToEntry_[tag] = &emplace_result.first->second;
				
				}

				bool hasTag(Int tag) const{
					
					if (tag < 0) {
						return false;
					}

					if (static_cast<Index>(tag) >= tagToEntry_.size()){
						return false;
					}

					return (tagToEntry_[tag] != nullptr);

				}

				// does this element have at least one face on a tag this registry is tracking
				bool touchesAny(const Int* rngTags, Index facesPerElement) const {
					for (Index f = 0; f < facesPerElement; ++f) {
						if (hasTag(rngTags[f])) return true;
					}
					return false;
				}

				// O(1) hot-path accessor for the per-face/per-element loop
				Entry& entryForTag(Int tag){
					
					return (*tagToEntry_[tag]);

				}

				// zeroes every registered tag's slot & measureSum, for a fresh evaluate() pass
				void reset() {
					for (auto& [tag, entry] : entries_) {
						entry.result.fill(0.0);
						entry.measureSum = 0.0;
					}
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

				std::unordered_map<Int, Entry> entries_;   // owning store; cold-path iteration (reset(), result())
				std::vector<Entry*> tagToEntry_;            // dense, direct-indexed by tag; nullptr if unregistered

			}; // class BoundaryQuantityRegistry

			// accumulates one boundary quadrature point's contribution to a tag's measureSum (length/area/volume)
			template<Index SpatialDim, typename EvalQPBdyT>
			PDE_HOST PDE_DEVICE void accumulateMeasure(const EvalQPBdyT& qp, Real& measureSum) {
				
				Real normalMag = 0;
				
				for (Index sD = 0; sD < SpatialDim; ++sD){
					normalMag += qp.normal[sD] * qp.normal[sD];
				}
				
				measureSum += std::sqrt(normalMag) * qp.w;
			}

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#endif
