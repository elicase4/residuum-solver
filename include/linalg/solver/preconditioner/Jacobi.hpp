#ifndef RESIDUUM_LINALG_SOLVER_PRECONDITIONER_JACOBI_HPP
#define RESIDUUM_LINALG_SOLVER_PRECONDITIONER_JACOBI_HPP

#include <memory>

#include "linalg/operations/VectorOps.hpp"
#include "linalg/operator/Operator.hpp"

namespace residuum {
	namespace linalg {
		namespace solver {
			namespace preconditioner {

				template<typename VectorT>
				class Jacobi {
				public:

					using DataType = typename VectorT::value_type;

					// default-constructible (like Identity) so CGRunner/GMRESRunner don't need to
					// thread a size through their constructors; sized lazily on the first update(),
					// using the operator's own size() -- the preconditioner pulls what it needs from
					// the operator itself, same as the update() contract everywhere else
					Jacobi() = default;

					template<typename OperatorT>
					requires linalg::op::DiagonalExtractable<OperatorT, VectorT>
					void update(const OperatorT& A) {

						if (!invDiag_) invDiag_ = std::make_unique<VectorT>(A.size());

						A.diagonal(*invDiag_);

						DataType* d = invDiag_->data();
						for (Index i = 0; i < invDiag_->size(); ++i) {
							d[i] = (std::abs(d[i]) > DataType(1e-300)) ? DataType(1) / d[i] : DataType(1);
						}

					}

					void apply(const VectorT& r, VectorT& z) {
						operations::multiply(r, *invDiag_, z); // z = r .* invDiag (elementwise)
					}

					Index flopsPerApply() const {
						return invDiag_ ? invDiag_->size() : Index(0);
					}

				private:
					std::unique_ptr<VectorT> invDiag_;

				}; // class Jacobi

			} // namespace preconditioner
		} // namespace solver
	} // namespace linalg
} // namespace residuum

#endif
