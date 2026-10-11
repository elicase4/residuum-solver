#include "linalg/types/backend/CPU.hpp"

namespace residuum::fem::form {

	template<typename FormsT, typename SelectorT>
	struct FormsWithOperand<linalg::types::backend::CPU, FormsT, SelectorT> {
		const FormsT& forms; // host memory already lives where the assembler runs, so a zero-copy reference is correct; a CUDA specialization would cache lowered/device-resident values here instead, since expression-backed forms can never run on device regardless of storage
		SelectorT selector;
	}; // struct FormsWithOperand<linalg::types::backend::CPU, FormsT, SelectorT>

} // namespace residuum::fem::form
