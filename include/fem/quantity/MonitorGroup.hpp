#ifndef RESIDUUM_FEM_QUANTITY_MONITORGROUP_HPP
#define RESIDUUM_FEM_QUANTITY_MONITORGROUP_HPP

#include <functional>
#include <stdexcept>
#include <vector>

#include "core/Types.hpp"

namespace residuum {
	namespace fem {
		namespace quantity {

			class MonitorGroup {
			public:

				virtual ~MonitorGroup() = default;

				// setup-time only
				virtual void registerTag(Int tag) = 0;
				virtual Index addCombination() = 0;
				virtual void addTerm(Index combinationIndex, Int tag, Real coefficient) = 0;

				// setup-time only; separate index space from addCombination/addTerm, for monitors that
				// project a multi-component result (e.g. a force vector) to one scalar via per-component
				// weights rather than a single coefficient applied uniformly (e.g. lift). Default throws:
				// most monitors never need this, so it isn't pure virtual.
				virtual Index addWeightedCombination() {
					throw std::runtime_error("MonitorGroup: this monitor does not support weighted combinations");
				}

				virtual void addWeightedTerm(Index, Int, const std::vector<Real>&) {
					throw std::runtime_error("MonitorGroup: this monitor does not support per-component weighted terms");
				}

				virtual void evaluate(Real time, const std::function<void(Index, const Real*)>& callback) = 0;

			}; // class MonitorGroup

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#endif
