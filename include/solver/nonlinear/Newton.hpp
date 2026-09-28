#ifndef RESIDUUM_SOLVER_NONLINEAR_NEWTON_HPP
#define RESIDUUM_SOLVER_NONLINEAR_NEWTON_HPP

#include <string>
#include <utility>
#include <vector>

#include "core/Types.hpp"

#include "fem/dof/DOFOrdering.hpp"
#include "fem/dof/DOFResidualNorms.hpp"

#include "linalg/solver/base/SolverReport.hpp"

#include "solver/config/NonlinearSolverConfig.hpp"
#include "solver/nonlinear/NonlinearSolverRunner.hpp"
#include "solver/stage/NonlinearCapableStage.hpp"

#include "utils/logging/nonlinear/Logger.hpp"

namespace residuum {
	namespace solver {
		namespace nonlinear {

			template<stage::NonlinearCapableStage StageT, typename VectorT>
			class NewtonRunner : public NonlinearSolverRunner<VectorT> {
			public:

				NewtonRunner(StageT& stage, const config::NonlinearSolverConfig& cfg, utils::logging::nonlinear::Logger logger, std::vector<std::string> dofNames, Index freeDOFsPerField, fem::dof::DOFOrdering dofOrdering)
					: stage_(stage), cfg_(cfg), logger_(std::move(logger)), dofNames_(std::move(dofNames)), freeDOFsPerField_(freeDOFsPerField), dofOrdering_(dofOrdering) {}

				bool solve(linalg::solver::SolverReport<VectorT>& report) override {

					// a fresh nonlinear solve (e.g. a new timestep), so tag CSV rows with a new outer tick
					logger_.reset();

					Real res0 = Real(0);
					std::vector<Real> res0PerDOF;

					for (Index iter = 0; iter < cfg_.maxIterations; ++iter) {

						stage_.assemble();
						const Real resNorm = stage_.residualNorm();
						const auto& R = stage_.residual();
						const auto resPerDOF = fem::dof::computePerDOFNorms(R.data(), R.size(), static_cast<Index>(dofNames_.size()), freeDOFsPerField_, dofOrdering_);

						if (iter == 0) { res0 = resNorm; res0PerDOF = resPerDOF; }
						const Real resRel = (res0 > Real(0)) ? (resNorm / res0) : Real(0);

						std::vector<Real> resRelPerDOF(resPerDOF.size());
						for (Index i = 0; i < static_cast<Index>(resPerDOF.size()); ++i) {
							resRelPerDOF[i] = (res0PerDOF[i] > Real(0)) ? (resPerDOF[i] / res0PerDOF[i]) : Real(0);
						}

						logger_.log(iter, resRelPerDOF);

						if (resNorm < cfg_.absoluteTolerance || resRel < cfg_.relativeTolerance) {

							report.converged = true;
							report.iterations = iter;
							report.initialResidual = res0;
							report.finalResidual = resNorm;
							report.finalResidualRel = resRel;

							logger_.summary(true);
							return true;

						}

						if (!stage_.solveLinearStep()) {
							report.converged = false;
							logger_.summary(false);
							return false;
						}

					}

					report.converged = false;
					report.iterations = cfg_.maxIterations;
					logger_.summary(false);
					return false;

				}

			private:

				StageT& stage_;
				config::NonlinearSolverConfig cfg_;
				utils::logging::nonlinear::Logger logger_;
				std::vector<std::string> dofNames_;
				Index freeDOFsPerField_;
				fem::dof::DOFOrdering dofOrdering_;

			}; // class NewtonRunner

		} // namespace nonlinear
	} // namespace solver
} // namespace residuum

#endif
