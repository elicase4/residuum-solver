#ifndef RESIDUUM_EQUATION_HEATEQ_EVAL_DEFAULTMODEL_HPP
#define RESIDUUM_EQUATION_HEATEQ_EVAL_DEFAULTMODEL_HPP

#include "fem/evaluator/EvalModel.hpp"

namespace residuum::equation::heateq::evaluator {

	template<typename QuadraturePointT>
	struct DefaultModel {
		
		void eval(QuadraturePointT&) const {}

		void evalGradient(QuadraturePointT&) const {}
		
	}; // struct DefaultModel

} // residuum::equation::heateq

#endif
