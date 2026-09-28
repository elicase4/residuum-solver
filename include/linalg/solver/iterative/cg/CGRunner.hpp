#ifndef RESIDUUM_LINALG_SOLVER_ITERATIVE_CG_CGRUNNER_HPP
#define RESIDUUM_LINALG_SOLVER_ITERATIVE_CG_CGRUNNER_HPP

#include "core/Types.hpp"

#include "linalg/solver/base/SolverReport.hpp"
#include "linalg/solver/base/LinearSolverRunner.hpp"
#include "linalg/solver/base/LinearSolverAlgorithm.hpp"

#include "linalg/solver/iterative/cg/Config.hpp"
#include "linalg/solver/iterative/cg/Workspace.hpp"
#include "linalg/solver/iterative/cg/Solver.hpp"

namespace residuum {
	namespace linalg {
		namespace solver {
			namespace iterative {
				namespace cg {

					template<typename OperatorT, typename VectorT, typename PreconditionerT, typename LoggerT>
					class CGRunner : public linalg::solver::LinearSolverRunner<VectorT> {
					public:

						using AlgorithmT = Solver<OperatorT, VectorT, PreconditionerT, LoggerT>;

						static_assert( linalg::solver::LinearSolverAlgorithm<AlgorithmT, OperatorT, VectorT, PreconditionerT, LoggerT>, "CGRunner: cg::Solver no longer satisfies linalg::solver::LinearSolverAlgorithm. Check its Config/Workspace aliases and solve() signature");

						CGRunner(const OperatorT& op, Index n, const typename AlgorithmT::Config& cfg, LoggerT logger);

						bool solve(const VectorT& b, VectorT& x, linalg::solver::SolverReport<VectorT>& report) override;

					private:

						OperatorT op_;
						typename AlgorithmT::Workspace workspace_;
						PreconditionerT preconditioner_;
						LoggerT logger_;
						AlgorithmT solver_;

					}; // class CGRunner

				} // namespace cg
			} // namespace iterative
		} // namespace solver
	} // namespace linalg
} // namespace residuum

#include "linalg/solver/iterative/cg/CGRunner.tpp"

#endif
