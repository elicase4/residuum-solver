#ifndef RESIDUUM_EQUATION_HEATEQ_EVAL_EVALELEMENT_HPP
#define RESIDUUM_EQUATION_HEATEQ_EVAL_EVALELEMENT_HPP

#include <utility>

#include "fem/evaluator/EvalElement.hpp"

namespace residuum::equation::heateq::evaluator {

	template<typename BasisT, Index SD>
	class EvalElement {
	public:

		static constexpr Index SpatialDim = SD;
		static constexpr Index ParametricDim = BasisT::ParametricDim;

		explicit EvalElement(BasisT basis) : basis_(std::move(basis)) {}

		Index nodesPerElement() const { return basis_.nodesPerElement(); }
		const BasisT& basis() const { return basis_; }

		// node coordinates
		const Real* nodeCoords;

		// time coordinate
		Real t;

		PDE_HOST PDE_DEVICE void bindElement(const Real* coords, const Real time){

			nodeCoords = coords;
			t = time;

		}

	private:
		BasisT basis_;

	}; // class EvalElement

} // namespace residuum::equation::heateq::evaluator

#endif
