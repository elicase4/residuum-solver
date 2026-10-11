#ifndef RESIDUUM_FEM_FORM_FORMSWITHOPERAND_HPP
#define RESIDUUM_FEM_FORM_FORMSWITHOPERAND_HPP

namespace residuum {
	namespace fem {
		namespace form {

			// pairs a forms registry with a selector that reads one named buffer off a shared
			// StateVector, letting Assembler/QuantityEvaluator fuse multiple form registries into
			// one element/quadrature-point loop; storage is backend-specific (see backend/cpu)
			template<typename BackendT, typename FormsT, typename SelectorT>
			struct FormsWithOperand;

		} // namespace form
	} // namespace fem
} // namespace residuum

#include "backend/cpu/FormsWithOperand.tpp"

#endif
