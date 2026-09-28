#include <cmath>
#include <cassert>

namespace residuum::linalg::operations {

	// c = y^T x
	template<typename VectorT>
	typename VectorT::value_type dot(const VectorT& x, const VectorT& y){
		
		assert(x.size() == y.size());

		typename VectorT::value_type result = 0;
		for (Index i = 0; i < x.size(); ++i){
			result += x.data()[i] * y.data()[i];
		}

		return result;

	}

	// y = y + alpha*x
	template<typename VectorT>
	void axpy(const typename VectorT::value_type alpha, const VectorT& x, VectorT& y){

		assert(x.size() == y.size());

		for (Index i = 0; i < y.size(); ++i){
			y.data()[i] += alpha * x.data()[i];
		}

	}

	// y = y + alpha*x
	template<typename VectorT>
	void axpby(const typename VectorT::value_type alpha, const typename VectorT::value_type beta, const VectorT& x, VectorT& y){

		assert(x.size() == y.size());

		for (Index i = 0; i < y.size(); ++i){
			y.data()[i] = beta * y.data()[i] + alpha * x.data()[i];
		}

	}

	// x = alpha*x
	template<typename VectorT>
	void scal(const typename VectorT::value_type alpha, VectorT& x){

		for (Index i = 0; i < x.size(); ++i){
			x.data()[i] *= alpha;
		}

	}

	// y = x
	template<typename VectorT>
	void copy(const VectorT& x, VectorT& y){

		assert(x.size() == y.size());

		for (Index i = 0; i < y.size(); ++i){
			y.data()[i] = x.data()[i];
		}
	
	}

	// z = x .* y (elementwise)
	template<typename VectorT>
	void multiply(const VectorT& x, const VectorT& y, VectorT& z){

		assert(x.size() == y.size());
		assert(x.size() == z.size());

		for (Index i = 0; i < z.size(); ++i){
			z.data()[i] = x.data()[i] * y.data()[i];
		}

	}

	// c = || x ||_2
	template<typename VectorT>
	typename VectorT::value_type norm(const VectorT& x){

		return std::sqrt(dot(x,x));

	}

} // namespace residuum::linalg::operations
