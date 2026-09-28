namespace residuum {
	namespace linalg {
		namespace solver {
			namespace iterative {
				namespace cg {

					template<typename OperatorT, typename VectorT, typename PreconditionerT, typename LoggerT>
					CGRunner<OperatorT, VectorT, PreconditionerT, LoggerT>::CGRunner(const OperatorT& op, Index n, const typename AlgorithmT::Config& cfg, LoggerT logger) : op_(op), workspace_(n), preconditioner_() , logger_(std::move(logger)) , solver_(cfg) {}

					template<typename OperatorT, typename VectorT, typename PreconditionerT, typename LoggerT>
					bool CGRunner<OperatorT, VectorT, PreconditionerT, LoggerT>::solve(const VectorT& b, VectorT& x, linalg::solver::SolverReport<VectorT>& report) {
						return solver_.solve(report, logger_, workspace_, preconditioner_, op_, b, x);
					}

				} // namespace cg
			} // namespace iterative
		} // namespace solver
	} // namespace linalg
} // namespace residuum
