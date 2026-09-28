#ifndef RESIDUUM_FEM_BOUNDARY_BOUNDARYCONDITION_HPP
#define RESIDUUM_FEM_BOUNDARY_BOUNDARYCONDITION_HPP

#include <concepts>

#include "core/Types.hpp"
#include "BoundaryCategory.hpp"

namespace residuum {
	namespace fem {
		namespace boundary {

			template<typename FunctionT>
			concept BoundaryFunction = requires(const FunctionT f, Real time, const Real* x, Real* value) {
				{ FunctionT::NumComponents } -> std::convertible_to<Index>;
				{ f.eval(time, x, value) } -> std::same_as<void>;
			}; // concept BoundaryFunction

			template<typename FunctionT>
			struct BoundaryCondition {

				static constexpr Index NumComponents = FunctionT::NumComponents;
				Int tag;
				BCCategory componentType[NumComponents];
				FunctionT function;
			
			}; // struct BoundaryCondition

		} // namespace boundary
	} // namespace fem
} // namespace residuum

#endif
