#ifndef RESIDUUM_LINALG_SOLVER_ITERATIVE_CG_SOLVER_HPP
#define RESIDUUM_LINALG_SOLVER_ITERATIVE_CG_SOLVER_HPP

#include <cmath>
#include <cassert>
#include <vector>

#include "linalg/solver/iterative/cg/Config.hpp"
#include "linalg/solver/iterative/cg/Workspace.hpp"

#include "linalg/operations/VectorOps.hpp"
#include "linalg/operator/Operator.hpp"
#include "linalg/solver/base/SolverReport.hpp"

namespace residuum {
	namespace linalg {
		namespace solver {
			namespace iterative {
				namespace cg {

					template<typename OperatorT, typename VectorT, typename PreconditionerT, typename LoggerT>
					requires linalg::op::LinearOperator<OperatorT, VectorT>
					class Solver {
					public:

						using Config = residuum::linalg::solver::iterative::cg::Config<VectorT>;
						using Workspace = residuum::linalg::solver::iterative::cg::Workspace<VectorT>;

						Solver(const Config& cfg) : config(cfg) {}

						bool solve(solver::SolverReport<VectorT>& report, LoggerT& logger, Workspace& W, PreconditionerT& M, const OperatorT& A, const VectorT& b, VectorT& x);

					private:
						Config config;

					}; // class Solver

				} // namespace cg
			} // namespace iterative
		} // namespace solver
	} // namespace linalg
} // namespace residuum

#include "Solver.tpp"

#endif
