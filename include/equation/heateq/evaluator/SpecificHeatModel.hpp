#ifndef RESIDUUM_EQUATION_HEATEQ_EVAL_SPECIFICHEATMODEL_HPP
#define RESIDUUM_EQUATION_HEATEQ_EVAL_SPECIFICHEATMODEL_HPP

#include <optional>
#include <string>

#include "equation/heateq/evaluator/TemperatureExpression.hpp"
#include "fem/evaluator/EvalModel.hpp"

namespace residuum::equation::heateq::evaluator {

	template<typename QuadraturePointT>
	struct SpecificHeatModel {

		enum class Dependence { Constant, TemperatureDependent };

		Dependence dependence = Dependence::Constant;

		Real value = Real(0);

		std::optional<TemperatureExpression> valueExpr;
		std::optional<TemperatureExpression> gradientExpr;

		void setConstant(Real cp) {
			dependence = Dependence::Constant;
			value = cp;
		}

		void setTemperatureDependent(const std::string& valueExpression, const std::string& gradientExpression) {
			dependence = Dependence::TemperatureDependent;
			valueExpr.emplace(valueExpression);
			gradientExpr.emplace(gradientExpression);
		}

		void eval(QuadraturePointT& qp) const {

			Real cpVal = value;

			if (dependence == Dependence::TemperatureDependent) {
				if constexpr (requires { qp.T.value; }) {
					cpVal = (*valueExpr)(qp.T.value);
				}
			}

			qp.cp = cpVal;

		}

		void evalGradient(QuadraturePointT& qp) const {

			Real dcp = Real(0);

			if (dependence == Dependence::TemperatureDependent) {
				if constexpr (requires { qp.T.value; }) {
					dcp = (*gradientExpr)(qp.T.value);
				}
			}

			qp.dcpdT = dcp;

		}

	}; // struct SpecificHeatModel

} // namespace residuum::equation::heateq::evaluator

#endif
