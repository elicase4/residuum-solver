#ifndef RESIDUUM_FEM_FORM_LINEARFORM_HPP
#define RESIDUUM_FEM_FORM_LINEARFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/evaluator/EvalQuadraturePointVolume.hpp"

namespace residuum {
	namespace fem {
		namespace form {
			
			template<typename FormT, typename QuadraturePointT>
			concept LinearForm = requires (const FormT f, const QuadraturePointT& qp, Real* Ue, Real* Fe) {
				{ f.computeElementLevelVector(qp, Ue, Fe) } -> std::same_as<void>;
			}; // concept LinearForm

		} // namespace form
	} // namespace fem
} // namespace residuum

#endif
