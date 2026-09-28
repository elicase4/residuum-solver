#ifndef RESIDUUM_LINALG_SOLVER_ITERATIVE_GMRES_CONFIG_HPP
#define RESIDUUM_LINALG_SOLVER_ITERATIVE_GMRES_CONFIG_HPP

namespace residuum {
	namespace linalg {
		namespace solver {
			namespace iterative {
				namespace gmres {

					enum class ToleranceType {Relative, Absolute}; // enum class Tolerance Type

					template<typename VectorT>
					struct Config {

						using DataType = typename VectorT::value_type;

						DataType tol = 1e-8;
						ToleranceType tolType = ToleranceType::Relative;
						Index maxIters = 1000; // total Arnoldi steps across all restart cycles

						Index krylovDim = 50; // restart length m; a fresh Arnoldi basis is built every m steps

						Index dofsPerNode = 1;

					}; // struct Config

				} // namespace gmres
			} // namespace iterative
		} // namespace solver
	} // namespace linalg
} // namespace residuum

#endif
