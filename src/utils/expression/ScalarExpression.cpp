#include "utils/expression/ScalarExpression.hpp"

residuum::utils::expression::ScalarExpression::ScalarExpression(const std::string& expression) {

	symbolTable_.add_variable("t", t_);
	symbolTable_.add_variable("x", x_);
	symbolTable_.add_variable("y", y_);
	symbolTable_.add_variable("z", z_);

	symbolTable_.add_constants();

	expression_.register_symbol_table(symbolTable_);

	if (!parser_.compile(expression, expression_)) {
		throw std::runtime_error("Failed to compile scalar expression:\n" + expression + "\n" + parser_.error());
	}

}

void residuum::utils::expression::ScalarExpression::operator()(Real t, const Real* x, Real* out) const {
	
	t_ = t;

	x_ = x[0];
	y_ = x[1];
	z_ = x[2];

	out[0] = expression_.value();

}
