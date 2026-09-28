#ifndef RESIDUUM_LINALG_SOLVER_ITERATIVE_GMRES_WORKSPACE_HPP
#define RESIDUUM_LINALG_SOLVER_ITERATIVE_GMRES_WORKSPACE_HPP

#include <vector>

namespace residuum {
	namespace linalg {
		namespace solver {
			namespace iterative {
				namespace gmres {

					template<typename VectorT>
					struct Workspace {

						using DataType = typename VectorT::value_type;

						std::vector<VectorT> Q; // m+1 orthonormal Krylov basis vectors
						VectorT w; // scratch: A * (M^-1 * q_j) each Arnoldi step
						VectorT z; // scratch: preconditioner apply, and the Krylov-space correction at restart

						// the reduced (m+1)x(m) least-squares problem is dense and tiny (m ~ 20-50); it
						// stays host-side regardless of backend, there is nothing to gain from device residence
						std::vector<DataType> H; // Hessenberg, column-major: H[j*(m+1)+i] = row i, col j
						std::vector<DataType> cs; // Givens rotation cosines, size m
						std::vector<DataType> sn; // Givens rotation sines, size m
						std::vector<DataType> g; // rotated RHS of the reduced problem, size m+1

						Workspace(Index n, Index m) : w(n), z(n), H(static_cast<std::size_t>(m + 1) * m), cs(m), sn(m), g(m + 1) {
							Q.reserve(m + 1);
							for (Index i = 0; i <= m; ++i) Q.emplace_back(n);
						}

					}; // struct Workspace

				} // namespace gmres
			} // namespace iterative
		} // namespace solver
	} // namespace linalg
} // namespace residuum

#endif
