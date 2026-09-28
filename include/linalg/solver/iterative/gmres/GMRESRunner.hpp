#ifndef RESIDUUM_LINALG_SOLVER_ITERATIVE_GMRES_GMRESRUNNER_HPP
#define RESIDUUM_LINALG_SOLVER_ITERATIVE_GMRES_GMRESRUNNER_HPP

#include "core/Types.hpp"

#include "linalg/solver/base/SolverReport.hpp"
#include "linalg/solver/base/LinearSolverRunner.hpp"
#include "linalg/solver/base/LinearSolverAlgorithm.hpp"

#include "linalg/solver/iterative/gmres/Config.hpp"
#include "linalg/solver/iterative/gmres/Workspace.hpp"
#include "linalg/solver/iterative/gmres/Solver.hpp"

namespace residuum {
	namespace linalg {
		namespace solver {
			namespace iterative {
				namespace gmres {

					template<typename OperatorT, typename VectorT, typename PreconditionerT, typename LoggerT>
					class GMRESRunner : public linalg::solver::LinearSolverRunner<VectorT> {
					public:

						using AlgorithmT = Solver<OperatorT, VectorT, PreconditionerT, LoggerT>;

						static_assert( linalg::solver::LinearSolverAlgorithm<AlgorithmT, OperatorT, VectorT, PreconditionerT, LoggerT>, "GMRESRunner: gmres::Solver no longer satisfies linalg::solver::LinearSolverAlgorithm. Check its Config/Workspace aliases and solve() signature");

						GMRESRunner(const OperatorT& op, Index n, const typename AlgorithmT::Config& cfg, LoggerT logger);

						bool solve(const VectorT& b, VectorT& x, linalg::solver::SolverReport<VectorT>& report) override;

					private:

						OperatorT op_;
						typename AlgorithmT::Workspace workspace_;
						PreconditionerT preconditioner_;
						LoggerT logger_;
						AlgorithmT solver_;

					}; // class GMRESRunner

				} // namespace gmres
			} // namespace iterative
		} // namespace solver
	} // namespace linalg
} // namespace residuum

#include "linalg/solver/iterative/gmres/GMRESRunner.tpp"

#endif
