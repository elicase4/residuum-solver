#ifndef RESIDUUM_LINALG_SOLVER_BASE_LINEARSOLVERALGORITHM_HPP
#define RESIDUUM_LINALG_SOLVER_BASE_LINEARSOLVERALGORITHM_HPP

#include <concepts>

#include "linalg/solver/base/SolverReport.hpp"

namespace residuum {
	namespace linalg {
		namespace solver {

			template<typename ST, typename OperatorT, typename VectorT, typename PreconditionerT, typename LoggerT>
			concept LinearSolverAlgorithm = requires(ST& s, const typename ST::Config& cfg, SolverReport<VectorT>& report, LoggerT& logger, typename ST::Workspace& workspace, PreconditionerT& preconditioner, const OperatorT& op, const VectorT& b, VectorT& x) {
				{ ST(cfg) };
				{ s.solve(report, logger, workspace, preconditioner, op, b, x) } -> std::same_as<bool>;
			}; // concept LinearSolverAlgorithm

		} // namespace solver
	} // namespace linalg
} // namespace residuum

#endif
