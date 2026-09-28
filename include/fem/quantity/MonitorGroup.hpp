#ifndef RESIDUUM_FEM_QUANTITY_MONITORGROUP_HPP
#define RESIDUUM_FEM_QUANTITY_MONITORGROUP_HPP

#include <functional>

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

				virtual void evaluate(Real time, const std::function<void(Index, const Real*)>& sink) = 0;

			}; // class MonitorGroup

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#endif
