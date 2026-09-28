#include "utils/expression/TensorExpression.hpp"

residuum::utils::expression::TensorExpression::TensorExpression(const std::vector<std::vector<std::string>>& expressions) {

	symbolTable_.add_variable("t", t_);
	symbolTable_.add_variable("x", x_);
	symbolTable_.add_variable("y", y_);
	symbolTable_.add_variable("z", z_);

	symbolTable_.add_constants();

	rows_ = expressions.size();
	cols_ = rows_ ? expressions.front().size() : 0;

	expressions_.resize(rows_ * cols_);

	Index k = 0;

	for (Index i = 0; i < rows_; ++i) {

		if (expressions[i].size() != cols_) {
			throw std::runtime_error("TensorExpression: inconsistent row sizes");
		}
	
		for (Index j = 0; j < cols_; ++j) {

			expressions_[k].register_symbol_table(symbolTable_);

			if (!parser_.compile(expressions[i][j], expressions_[k])) {
				throw std::runtime_error("Failed to compile vector expression:\n" + expressions[i][j] + "\n" + parser_.error());
			}

			++k;

		}

	}

}

void residuum::utils::expression::TensorExpression::operator()(Real t, const Real* x, Real* out) const {
	
	t_ = t;

	x_ = x[0];
	y_ = x[1];
	z_ = x[2];
	
	for (Index i = 0; i < expressions_.size(); ++i) {
		out[i] = expressions_[i].value();
	}

}
