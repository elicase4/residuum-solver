#ifndef RESIDUUM_LINALG_SOLVER_ITERATIVE_GMRES_SOLVER_HPP
#define RESIDUUM_LINALG_SOLVER_ITERATIVE_GMRES_SOLVER_HPP

#include <cmath>
#include <cassert>
#include <vector>

#include "linalg/solver/iterative/gmres/Config.hpp"
#include "linalg/solver/iterative/gmres/Workspace.hpp"

#include "linalg/operations/VectorOps.hpp"
#include "linalg/operator/Operator.hpp"
#include "linalg/solver/base/SolverReport.hpp"

namespace residuum {
	namespace linalg {
		namespace solver {
			namespace iterative {
				namespace gmres {

					// Right-preconditioned, restarted GMRES(m): builds an orthonormal Arnoldi basis for
					// K_m(A*M^-1, r0) and minimizes the true residual ||b - A*x|| over it every restart cycle.
					template<typename OperatorT, typename VectorT, typename PreconditionerT, typename LoggerT>
					requires linalg::op::LinearOperator<OperatorT, VectorT>
					class Solver {
					public:

						using Config = residuum::linalg::solver::iterative::gmres::Config<VectorT>;
						using Workspace = residuum::linalg::solver::iterative::gmres::Workspace<VectorT>;

						Solver(const Config& cfg) : config(cfg) {}

						bool solve(solver::SolverReport<VectorT>& report, LoggerT& logger, Workspace& W, PreconditionerT& M, const OperatorT& A, const VectorT& b, VectorT& x);

					private:
						Config config;

					}; // class Solver

				} // namespace gmres
			} // namespace iterative
		} // namespace solver
	} // namespace linalg
} // namespace residuum

#include "Solver.tpp"

#endif
