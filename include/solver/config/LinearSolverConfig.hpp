#ifndef RESIDUUM_SOLVER_CONFIG_LINEARSOLVERCONFIG_HPP
#define RESIDUUM_SOLVER_CONFIG_LINEARSOLVERCONFIG_HPP

#include <variant>

#include "core/Types.hpp"

namespace residuum {
	namespace solver {
		namespace config {

			struct IdentityPreconditionerParams {};

			struct JacobiPreconditionerParams {};

			struct PreconditionerConfig {

				enum class Type {
					Identity,
					Jacobi
				};

				Type type = Type::Identity;

				std::variant<IdentityPreconditionerParams, JacobiPreconditionerParams> params = IdentityPreconditionerParams{};

			}; // struct PreconditionerConfig

			// each linear solver type's own parameters
			struct CGParams {};

			struct GMRESParams {
				Index krylovDim = 50; // restart length m
			};

			struct BiCGSTABParams {};

			struct LUParams {};

			struct LinearSolverConfig {

				enum class Type {
					CG,
					GMRES,
					BiCGSTAB,
					LU
				};

				enum class OperatorType {
					CSR,
					FEM
				};

				Type type = Type::CG;

				OperatorType operatorType = OperatorType::CSR;

				PreconditionerConfig preconditioner;

				Real tolerance = 1e-10;

				Index maxIterations = 1000;

				std::variant<CGParams, GMRESParams, BiCGSTABParams, LUParams> params = CGParams{};

			};

		} // namespace config
	} // namespace solver
} // namespace residuum

#endif
