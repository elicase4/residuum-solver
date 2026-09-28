#ifndef RESIDUUM_SOLVER_STAGE_STEADYSTAGE_HPP
#define RESIDUUM_SOLVER_STAGE_STEADYSTAGE_HPP

#include <memory>
#include <type_traits>

#include "core/Types.hpp"

#include "fem/assembly/ElementMap.hpp"
#include "fem/form/FormRegistry.hpp"
#include "linalg/operations/VectorOps.hpp"
#include "linalg/solver/base/LinearSolverRunner.hpp"
#include "linalg/solver/base/SolverReport.hpp"

#include "solver/SolverInstance.hpp"
#include "solver/config/LinearSolverConfig.hpp"
#include "solver/nonlinear/NonlinearSolverFactory.hpp"
#include "solver/nonlinear/NonlinearSolverRunner.hpp"
#include "solver/problem/Problem.hpp"

namespace residuum {
	namespace solver {
		namespace stage {

			template<problem::Problem ProblemT>
			class SteadyStage {
			public:

				using VectorT = typename ProblemT::VectorT;
				using MatrixT = typename ProblemT::MatrixT;

				// derived from the accessors
				using StiffnessFormsT = std::decay_t<decltype(std::declval<ProblemT&>().stiffnessForms())>;
				using TangentStiffnessFormsT = std::decay_t<decltype(std::declval<ProblemT&>().tangentStiffnessForms())>;

				// J = K + K_T
				using JacobianFormsT = fem::form::FormRegistry<StiffnessFormsT, TangentStiffnessFormsT>;

				explicit SteadyStage(ProblemT& problem) : problem_(problem), K_(needsK() ? problem_.createMatrix() : MatrixT(0, 0)), R_(problem_.createVector()) {

					linearSolverRunner_ = problem_.template makeLinearRunner<fem::assembly::GatherMode::Free>(&currentTime_, problem_.stiffnessForms(), problem_.stiffnessModel(), nullptr, {}, K_);

					if (solver::isNonlinear(problem_.solverInstance().mode)) {

						J_ = std::make_unique<MatrixT>(isSparse() ? problem_.createMatrix() : MatrixT(0, 0));
						dU_ = std::make_unique<VectorT>(problem_.createVector());
						jacobianSolverRunner_ = problem_.template makeLinearRunner<fem::assembly::GatherMode::Free>(&currentTime_, jacobianForms_, problem_.stiffnessModel(), &problem_.U(), {}, *J_);

						nonlinearSolverRunner_ = nonlinear::makeNonlinearSolverRunner<SteadyStage, VectorT>(*this, *problem_.solverInstance().nonlinear, problem_.loggingConfig().nonlinear, problem_.equationLabel(), problem_.dofNames(), problem_.freeDOFsPerField(), problem_.dofOrdering());

					}

				}

				void initialize() {}

				void assemble() {

					if (needsK()) {
						problem_.template assembleMatrix<fem::assembly::GatherMode::Free>(currentTime_, problem_.stiffnessForms(), problem_.stiffnessModel(), {}, K_);
					}
					problem_.assembleLoad(currentTime_);
					problem_.applyNatural(currentTime_);

					if (!nonlinearSolverRunner_) {
						problem_.applyEssential(currentTime_, problem_.stiffnessForms(), problem_.stiffnessModel(), problem_.F());
					}

				}

				bool solve() {

					linalg::solver::SolverReport<VectorT> report;

					if (nonlinearSolverRunner_) {
						return nonlinearSolverRunner_->solve(report);
					}

					return linearSolverRunner_->solve(problem_.F(), problem_.U(), report);

				}

				void finalize() {}

				// NonlinearCapableStage
				Real residualNorm() {

					problem_.template assembleResidual<fem::assembly::GatherMode::Full>(currentTime_, problem_.stiffnessForms(), problem_.stiffnessModel(), {}, R_);
					linalg::operations::axpby(Real(1), Real(-1), problem_.F(), R_);

					return linalg::operations::norm(R_);

				}

				// NonlinearCapableStage
				const VectorT& residual() const { return R_; }

				bool solveLinearStep() {

					if (needsJ()) {
						problem_.template assembleMatrix<fem::assembly::GatherMode::Full>(currentTime_, jacobianForms_, problem_.stiffnessModel(), {}, *J_);
					}

					dU_->zero();
					linalg::solver::SolverReport<VectorT> report;

					if (!jacobianSolverRunner_->solve(R_, *dU_, report)) return false;

					linalg::operations::axpy(Real(1), *dU_, problem_.U());

					return true;

				}

				decltype(auto) solution() const { return problem_.U(); }

			private:

				bool isSparse() const { return problem_.solverInstance().linear->operatorType == config::LinearSolverConfig::OperatorType::CSR; }
				bool needsK() const { return !solver::isNonlinear(problem_.solverInstance().mode) && isSparse(); }
				bool needsJ() const { return solver::isNonlinear(problem_.solverInstance().mode) && isSparse(); }

				ProblemT& problem_;

				Real currentTime_ = Real(0);

				MatrixT K_;
				VectorT R_;
				std::unique_ptr<linalg::solver::LinearSolverRunner<VectorT>> linearSolverRunner_;
				std::unique_ptr<nonlinear::NonlinearSolverRunner<VectorT>> nonlinearSolverRunner_;

				// nonlinear-only
				JacobianFormsT jacobianForms_;
				std::unique_ptr<MatrixT> J_;
				std::unique_ptr<VectorT> dU_;
				std::unique_ptr<linalg::solver::LinearSolverRunner<VectorT>> jacobianSolverRunner_;

			}; // class SteadyStage

		} // namespace stage
	} // namespace solver
} // namespace residuum

#endif
