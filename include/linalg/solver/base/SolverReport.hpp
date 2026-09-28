#ifndef RESIDUUM_LINALG_SOLVER_BASE_SOLVERREPORT_HPP
#define RESIDUUM_LINALG_SOLVER_BASE_SOLVERREPORT_HPP

#include <vector>

#include "core/Types.hpp"

namespace residuum {
	namespace linalg {
		namespace solver {
			
			template<typename VectorT>
			struct SolverReport {
				
				using DataType = typename VectorT::value_type;

				bool converged = false;
				Index iterations = 0;
				DataType initialResidual = 0.0; // || r_0 ||
				DataType finalResidual = 0.0; // || r_k ||
				DataType finalResidualRel = 0.0; // finalResidual / initialResidual
				std::vector<DataType> perFieldResidual;

				DataType reduction() const {
					return (initialResidual > DataType(0)) ? (finalResidual / initialResidual) : DataType(0);
				}

			}; // struct SolverReport

		} // namespace solver
	} // namespace linalg
} // namespace residuum

#endif
