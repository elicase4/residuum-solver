#ifndef RESIDUUM_SOLVER_NONLINEAR_NONLINEARSOLVERRUNNER_HPP
#define RESIDUUM_SOLVER_NONLINEAR_NONLINEARSOLVERRUNNER_HPP

#include "linalg/solver/base/SolverReport.hpp"

namespace residuum {
	namespace solver {
		namespace nonlinear {

			template<typename VectorT>
			class NonlinearSolverRunner {
			public:

				virtual ~NonlinearSolverRunner() = default;

				virtual bool solve(linalg::solver::SolverReport<VectorT>& report) = 0;

			}; // class NonlinearSolverRunner

		} // namespace nonlinear
	} // namespace solver
} // namespace residuum

#endif
