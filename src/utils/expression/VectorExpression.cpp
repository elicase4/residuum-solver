#include "utils/expression/VectorExpression.hpp"

residuum::utils::expression::VectorExpression::VectorExpression(const std::vector<std::string>& expressions) {

	symbolTable_.add_variable("t", t_);
	symbolTable_.add_variable("x", x_);
	symbolTable_.add_variable("y", y_);
	symbolTable_.add_variable("z", z_);

	symbolTable_.add_constants();

	expressions_.resize(expressions.size());

	for (Index i = 0; i < expressions.size(); ++i) {
		
		expressions_[i].register_symbol_table(symbolTable_);

		if (!parser_.compile(expressions[i], expressions_[i])) {
			throw std::runtime_error("Failed to compile vector expression:\n" + expressions[i] + "\n" + parser_.error());
		}

	}

}

void residuum::utils::expression::VectorExpression::operator()(Real t, const Real* x, Real* out) const {
	
	t_ = t;

	x_ = x[0];
	y_ = x[1];
	z_ = x[2];
	
	for (Index i = 0; i < expressions_.size(); ++i) {
		out[i] = expressions_[i].value();
	}

}
