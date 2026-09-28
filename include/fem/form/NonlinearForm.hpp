#ifndef RESIDUUM_FEM_FORM_NONLINEARFORM_HPP
#define RESIDUUM_FEM_FORM_NONLINEARFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/evaluator/EvalQuadraturePointVolume.hpp"

namespace residuum {
	namespace fem {
		namespace form {
			
			template<typename FormT, typename QuadraturePointT>
			concept NonlinearForm = requires (const FormT f, const QuadraturePointT& qp, Real* Ue, Real* Re) {
				{ f.computeElementLevelVector(qp, Ue, Re) } -> std::same_as<void>;
			}; // concept NonlinearForm
		
		} // namespace form
	} // namespace fem
} // namespace residuum

#endif
