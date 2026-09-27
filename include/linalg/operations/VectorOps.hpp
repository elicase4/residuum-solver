#ifndef RESIDUUM_LINALG_OPERATIONS_VECTOROPS_HPP
#define RESIDUUM_LINALG_OPERATIONS_VECTOROPS_HPP

#include "linalg/types/Vector.hpp"

namespace residuum {
	namespace linalg {
		namespace operations {

			// c = y^T x
			template<typename VectorT>
			typename VectorT::value_type dot(const VectorT& a, const VectorT& b);

			// y = y + alpha*x
			template<typename VectorT>
			void axpy(const typename VectorT::value_type alpha, const VectorT& x, VectorT& y);

			// y = beta*y + alpha*x
			template<typename VectorT>
			void axpby(const typename VectorT::value_type alpha, const typename VectorT::value_type beta, const VectorT& x, VectorT& y);

			// x = alpha*x
			template<typename VectorT>
			void scal(const typename VectorT::value_type alpha, VectorT& x);

			// y = x
			template<typename VectorT>
			void copy(const VectorT& x, VectorT& y);

			// z = x .* y (elementwise)
			template<typename VectorT>
			void multiply(const VectorT& x, const VectorT& y, VectorT& z);

			// c = || x ||_2
			template<typename VectorT>
			typename VectorT::value_type norm(const VectorT& x);

		} // namespace operations
	} // namespace linalg
} // namespace residuum

#include "linalg/operations/backend/cpu/VectorOps.tpp"
//#include "linalg/operations/backend/cuda/VectorOps.tpp"

#endif
