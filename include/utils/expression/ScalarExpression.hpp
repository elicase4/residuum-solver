#ifndef RESIDUUM_UTILS_EXPRESSION_SCALAREXPRESSION_HPP
#define RESIDUUM_UTILS_EXPRESSION_SCALAREXPRESSION_HPP

#include <exprtk.hpp>
#include <string>
#include <stdexcept>

#include "core/Types.hpp"

namespace residuum {
	namespace utils {
		namespace expression {

			class ScalarExpression {
			public:

				explicit ScalarExpression(const std::string& expression);

				void operator()(Real t, const Real* x, Real* out) const;

			private:

				mutable Real t_ = 0.0;
				mutable Real x_ = 0.0;
				mutable Real y_ = 0.0;
				mutable Real z_ = 0.0;

				exprtk::symbol_table<Real> symbolTable_;
				exprtk::expression<Real> expression_;
				exprtk::parser<Real> parser_;

			}; // class ScalarExpression

		} // namespace expression
	} // namespace utils
} // namespace residuum

#endif
