#ifndef RESIDUUM_UTILS_EXPRESSION_TENSOREXPRESSION_HPP
#define RESIDUUM_UTILS_EXPRESSION_TENSOREXPRESSION_HPP

#include <exprtk.hpp>
#include <string>
#include <stdexcept>
#include <vector>

#include "core/Types.hpp"

namespace residuum {
	namespace utils {
		namespace expression {

			class TensorExpression {
			public:

				explicit TensorExpression(const std::vector<std::vector<std::string>>& expressions);

				void operator()(Real t, const Real* x, Real* out) const;

				Index rows() const { return rows_; };
				
				Index cols() const { return cols_; };

			private:

				mutable Real t_ = 0.0;
				mutable Real x_ = 0.0;
				mutable Real y_ = 0.0;
				mutable Real z_ = 0.0;

				exprtk::symbol_table<Real> symbolTable_;
				std::vector<exprtk::expression<Real>> expressions_;
				exprtk::parser<Real> parser_;

				Index rows_;
				Index cols_;

			}; // class TensorExpression

		} // namespace expression
	} // namespace utils
} // namespace residuum

#endif
