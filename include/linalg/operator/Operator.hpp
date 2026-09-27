#ifndef RESIDUUM_LINALG_OPERATOR_OPERATOR_HPP
#define RESIDUUM_LINALG_OPERATOR_OPERATOR_HPP

#include <concepts>

namespace residuum {
	namespace linalg {
		namespace op {

			template<typename OperatorT, typename VectorT>
			concept LinearOperator = requires(const OperatorT A, const VectorT& x, VectorT& y){
				{ A.apply(x, y) } -> std::same_as<void>;
			};

			template<typename OperatorT, typename VectorT>
			concept DiagonalExtractable = requires(const OperatorT A, VectorT& d){
				{ A.diagonal(d) } -> std::same_as<void>;
			};

		} // namespace op
	} // namespace linalg
} // namespace residuum

#endif
