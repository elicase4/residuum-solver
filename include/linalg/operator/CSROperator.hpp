#ifndef RESIDUUM_LINALG_OPERATOR_CSROPERATOR_HPP
#define RESIDUUM_LINALG_OPERATOR_CSROPERATOR_HPP

#include "linalg/types/CSRMatrix.hpp"
#include "linalg/types/Vector.hpp"
#include "linalg/operations/MatrixOps.hpp"
#include "linalg/operator/Operator.hpp"

namespace residuum {
	namespace linalg {
		namespace op {

			template<typename MatrixT>
			class CSROperator {
			public:

				using value_type = typename MatrixT::value_type;

				const MatrixT& A;

				explicit CSROperator(const MatrixT& mat) : A(mat) {
					static_assert(LinearOperator<CSROperator, linalg::types::Vector<value_type, typename MatrixT::backend_type>>);
					static_assert(DiagonalExtractable<CSROperator, linalg::types::Vector<value_type, typename MatrixT::backend_type>>);
				}

				template<typename VectorT>
				void apply(const VectorT& x, VectorT& y) const {
					operations::matvec(A, x, y);
				}

				Index size() const {
					return A.nRows();
				}

				// SpMV: one multiply + one add per nonzero
				Index flopsPerApply() const {
					return 2 * A.rowPtr()[A.nRows()];
				}

				template<typename VectorT>
				void diagonal(VectorT& d) const {
					for (Index i = 0; i < A.nRows(); ++i) {
						d.data()[i] = A.data()[A.getDataIndex(i, i)];
					}
				}

			}; // class CSROperator

		} // namespace op
	} // namespace linalg
} // namespace residuum

#endif
