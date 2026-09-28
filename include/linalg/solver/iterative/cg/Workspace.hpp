#ifndef RESIDUUM_LINALG_SOLVER_ITERATIVE_CG_WORKSPACE_HPP
#define RESIDUUM_LINALG_SOLVER_ITERATIVE_CG_WORKSPACE_HPP

namespace residuum {
	namespace linalg {
		namespace solver {
			namespace iterative {
				namespace cg {

					template<typename VectorT>
					struct Workspace {

						VectorT r;
						VectorT p;
						VectorT Ap;
						VectorT z;

						explicit Workspace(Index n) : r(n), p(n), Ap(n), z(n) {}

					}; // struct Workspace

				} // namespace cg
			} // namespace iterative
		} // namespace solver
	} // namespace linalg
} // namespace residuum

#endif
