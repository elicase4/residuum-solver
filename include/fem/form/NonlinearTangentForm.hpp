#ifndef RESIDUUM_FEM_FORM_NONLINEARTANGENTFORM_HPP
#define RESIDUUM_FEM_FORM_NONLINEARTANGENTFORM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/evaluator/EvalQuadraturePointVolume.hpp"

namespace residuum {
	namespace fem {
		namespace form {
			
			template<typename FormT, typename QuadraturePointT>
			concept NonlinearTangentForm = requires (const FormT f, const QuadraturePointT& qp, Real* Ue, Real* Ke, Real* Oe) {
				{ f.computeElementLevelMatrix(qp, Ke) } -> std::same_as<void>;
				{ f.computeElementLevelVector(qp, Ue, Oe) } -> std::same_as<void>;
			}; // concept NonlinearTangentForm

		} // namespace form
	} // namespace fem
} // namespace residuum

#endif
