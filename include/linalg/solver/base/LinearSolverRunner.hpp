#ifndef RESIDUUM_LINALG_SOLVER_BASE_LINEARSOLVERRUNNER_HPP
#define RESIDUUM_LINALG_SOLVER_BASE_LINEARSOLVERRUNNER_HPP

#include "linalg/solver/base/SolverReport.hpp"

namespace residuum {
	namespace linalg {
		namespace solver {

			template<typename VectorT>
			class LinearSolverRunner {
			public:

				virtual ~LinearSolverRunner() = default;

				virtual bool solve(const VectorT& b, VectorT& x, SolverReport<VectorT>& report) = 0;

			}; // class LinearSolverRunner

		} // namespace solver
	} // namespace linalg
} // namespace residuum

#endif
