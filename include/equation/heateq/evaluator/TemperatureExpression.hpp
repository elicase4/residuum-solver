#ifndef RESIDUUM_EQUATION_HEATEQ_EVAL_TEMPERATUREEXPRESSION_HPP
#define RESIDUUM_EQUATION_HEATEQ_EVAL_TEMPERATUREEXPRESSION_HPP

#include <exprtk.hpp>
#include <stdexcept>
#include <string>

#include "core/Types.hpp"

namespace residuum::equation::heateq::evaluator {

	class TemperatureExpression {
	public:

		explicit TemperatureExpression(const std::string& expression) : source_(expression) {
			compile();
		}

		TemperatureExpression(const TemperatureExpression& other) : source_(other.source_) {
			compile();
		}

		TemperatureExpression& operator=(const TemperatureExpression& other) {
			if (this != &other) {
				source_ = other.source_;
				compile();
			}
			return *this;
		}

		Real operator()(Real T) const {
			T_ = T;
			return expression_.value();
		}

	private:

		void compile() {

			symbolTable_ = exprtk::symbol_table<Real>();
			symbolTable_.add_variable("T", T_);
			symbolTable_.add_constants();

			expression_ = exprtk::expression<Real>();
			expression_.register_symbol_table(symbolTable_);

			if (!parser_.compile(source_, expression_)) {
				throw std::runtime_error("Failed to compile temperature expression:\n" + source_ + "\n" + parser_.error());
			}

		}

		std::string source_;
		mutable Real T_ = 0.0;

		exprtk::symbol_table<Real> symbolTable_;
		exprtk::expression<Real> expression_;
		exprtk::parser<Real> parser_;

	}; // class TemperatureExpression

} // namespace residuum::equation::heateq::evaluator

#endif
