#ifndef RESIDUUM_EQUATION_HEATEQ_EVAL_DENSITYMODEL_HPP
#define RESIDUUM_EQUATION_HEATEQ_EVAL_DENSITYMODEL_HPP

#include "fem/evaluator/EvalModel.hpp"

namespace residuum::equation::heateq::evaluator {

	template<typename QuadraturePointT>
	struct DensityModel {

		Real value = Real(0);

		void eval(QuadraturePointT& qp) const { qp.rho = value; }
		void evalGradient(QuadraturePointT&) const {}

	}; // struct DensityModel

} // namespace residuum::equation::heateq::evaluator

#endif
