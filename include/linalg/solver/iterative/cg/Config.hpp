#ifndef RESIDUUM_LINALG_SOLVER_ITERATIVE_CG_CONFIG_HPP
#define RESIDUUM_LINALG_SOLVER_ITERATIVE_CG_CONFIG_HPP

namespace residuum {
	namespace linalg {
		namespace solver {
			namespace iterative {
				namespace cg {

					enum class ToleranceType {Relative, Absolute}; // enum class Tolerance Type

					template<typename VectorT>
					struct Config {
						
						using DataType = typename VectorT::value_type;

						DataType tol = 1e-8;
						ToleranceType tolType = ToleranceType::Relative;
						Index maxIters = 1000;

						Index dofsPerNode = 1;

					}; // struct Config

				} // namespace cg
			} // namespace iterative
		} // namespace solver
	} // namespace linalg
} // namespace residuum

#endif
