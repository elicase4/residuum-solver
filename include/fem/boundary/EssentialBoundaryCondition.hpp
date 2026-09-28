#ifndef RESIDUUM_FEM_BOUNDARY_ESSENTIALBOUNDARYCONDITION_HPP
#define RESIDUUM_FEM_BOUNDARY_ESSENTIALBOUNDARYCONDITION_HPP

#include <memory>
#include <unordered_map>
#include <vector>

#include "BoundaryCondition.hpp"

namespace residuum {
	namespace fem {
		namespace boundary {

			class EssentialBoundaryConditionBase {
			public:

				virtual ~EssentialBoundaryConditionBase() = default;

				virtual void eval(Real time, const Real* x, Real* value) const = 0;

				virtual BCCategory componentType(Index component) const = 0;

			}; // class EssentialBoundaryConditionBase

			template<BoundaryFunction FunctionT>
			class EssentialBoundaryCondition final : public EssentialBoundaryConditionBase {
			public:

				explicit EssentialBoundaryCondition(std::shared_ptr<BoundaryCondition<FunctionT>> bc) : bc_(std::move(bc)) {}

				void eval(Real time, const Real* x, Real* value) const override {
					bc_->function.eval(time, x, value);
				}

				BCCategory componentType(Index component) const override {
					return bc_->componentType[component];
				}

			private:

				std::shared_ptr<BoundaryCondition<FunctionT>> bc_;

			}; // class EssentialBoundaryCondiiton

		} // namespace boundary
	} // namespace fem
} // namespace residuum

#endif
