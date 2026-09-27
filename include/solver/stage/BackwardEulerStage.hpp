#ifndef RESIDUUM_SOLVER_STAGE_BACKWARDEULERSTAGE_HPP
#define RESIDUUM_SOLVER_STAGE_BACKWARDEULERSTAGE_HPP

#include <memory>
#include <type_traits>

#include "core/Types.hpp"

#include "fem/assembly/ElementMap.hpp"
#include "fem/evaluator/ModelRegistry.hpp"
#include "fem/form/FormRegistry.hpp"
#include "fem/form/ScaledForm.hpp"

#include "linalg/operations/VectorOps.hpp"
#include "linalg/solver/base/LinearSolverRunner.hpp"
#include "linalg/solver/base/SolverReport.hpp"

#include "solver/SolverInstance.hpp"
#include "solver/config/LinearSolverConfig.hpp"
#include "solver/nonlinear/NonlinearSolverFactory.hpp"
#include "solver/nonlinear/NonlinearSolverRunner.hpp"
#include "solver/problem/TransientCapableProblem.hpp"

namespace residuum {
	namespace solver {
		namespace stage {

			template<problem::TransientCapableProblem ProblemT>
			class BackwardEulerStage {
			public:

				using VectorT = typename ProblemT::VectorT;
				using MatrixT = typename ProblemT::MatrixT;

				// derived from the accessors
				using StiffnessFormsT = std::decay_t<decltype(std::declval<ProblemT&>().stiffnessForms())>;
				using StiffnessModelT = std::decay_t<decltype(std::declval<ProblemT&>().stiffnessModel())>;
				using MassFormsT = std::decay_t<decltype(std::declval<ProblemT&>().massForms())>;
				using MassModelT = std::decay_t<decltype(std::declval<ProblemT&>().massModel())>;
				using TangentStiffnessFormsT = std::decay_t<decltype(std::declval<ProblemT&>().tangentStiffnessForms())>;
				using TangentMassFormsT = std::decay_t<decltype(std::declval<ProblemT&>().tangentMassForms())>;

				using OperatorFormsT = fem::form::FormRegistry<StiffnessFormsT, fem::form::ScaledForm<ProblemT::NumDOFs, MassFormsT>>;
				using OperatorModelT = fem::evaluator::ModelRegistry<StiffnessModelT, MassModelT>;

				using JacobianFormsT = fem::form::FormRegistry<StiffnessFormsT, TangentStiffnessFormsT, fem::form::ScaledForm<ProblemT::NumDOFs, MassFormsT>, TangentMassFormsT>;

				explicit BackwardEulerStage(ProblemT& problem) : problem_(problem), operatorModel_(problem_.stiffnessModel(), problem_.massModel()), K_(needsK() ? problem_.createMatrix() : MatrixT(0, 0)), massScratch_(problem_.createVector()), Udot_(problem_.createVector()), R_(problem_.createVector()) {

					setDt(problem_.solverInstance().timestepper->stepSize.dt);

					linearSolverRunner_ = problem_.template makeLinearRunner<fem::assembly::GatherMode::Free>(&currentTime_, operatorForms_, operatorModel_, nullptr, {}, K_);

					if (solver::isNonlinear(problem_.solverInstance().mode)) {

						J_ = std::make_unique<MatrixT>(isSparse() ? problem_.createMatrix() : MatrixT(0, 0));
						dU_ = std::make_unique<VectorT>(problem_.createVector());
						jacobianSolverRunner_ = problem_.template makeLinearRunner<fem::assembly::GatherMode::Free>(&currentTime_, jacobianForms_, operatorModel_, &problem_.U(), {&Udot_}, *J_);

						nonlinearSolverRunner_ = nonlinear::makeNonlinearSolverRunner<BackwardEulerStage, VectorT>(*this, *problem_.solverInstance().nonlinear, problem_.loggingConfig().nonlinear, problem_.equationLabel(), problem_.dofNames(), problem_.freeDOFsPerField(), problem_.dofOrdering());

					}

				}

				void initialize() {}

				void assemble() {

					linalg::operations::copy(problem_.U(), Udot_);
					linalg::operations::axpby(-Real(1) / dt_, Real(1) / dt_, problem_.U_prev(), Udot_);

					if (needsK()) {
						problem_.template assembleMatrix<fem::assembly::GatherMode::Free>(currentTime_, operatorForms_, operatorModel_, {&Udot_}, K_);
					}

					problem_.assembleLoad(currentTime_);

					if (!nonlinearSolverRunner_) {
						problem_.template assembleVector<fem::assembly::GatherMode::Full>(currentTime_ - dt_, problem_.massForms(), problem_.massModel(), problem_.U_prev(), {}, massScratch_);
						linalg::operations::axpy(Real(1) / dt_, massScratch_, problem_.F());
					}

					problem_.applyNatural(currentTime_);

					if (!nonlinearSolverRunner_) {
						problem_.applyEssential(currentTime_, operatorForms_, operatorModel_, problem_.F());
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

					problem_.template assembleVector<fem::assembly::GatherMode::Full>(currentTime_, problem_.massForms(), problem_.massModel(), Udot_, {}, massScratch_, &problem_.U());
					linalg::operations::axpy(Real(1), massScratch_, R_);

					linalg::operations::axpby(Real(1), Real(-1), problem_.F(), R_);

					return linalg::operations::norm(R_);

				}

				// NonlinearCapableStage
				const VectorT& residual() const { return R_; }

				bool solveLinearStep() {

					if (needsJ()) {
						problem_.template assembleMatrix<fem::assembly::GatherMode::Full>(currentTime_, jacobianForms_, operatorModel_, {&Udot_}, *J_);
					}

					dU_->zero();
					linalg::solver::SolverReport<VectorT> report;

					if (!jacobianSolverRunner_->solve(R_, *dU_, report)) return false;

					linalg::operations::axpy(Real(1), *dU_, problem_.U());

					return true;

				}

				// TransientCapableStage
				void setDt(Real dt) {

					if (dt == dt_) return; // common case: the active step-size policy kept dt unchanged

					dt_ = dt;
					operatorForms_ = OperatorFormsT(problem_.stiffnessForms(), fem::form::ScaledForm<ProblemT::NumDOFs, MassFormsT>(Real(1) / dt_));
					jacobianForms_ = JacobianFormsT(problem_.stiffnessForms(), problem_.tangentStiffnessForms(), fem::form::ScaledForm<ProblemT::NumDOFs, MassFormsT>(Real(1) / dt_), problem_.tangentMassForms());

					if (needsK()) {
						problem_.template assembleMatrix<fem::assembly::GatherMode::Free>(currentTime_, operatorForms_, operatorModel_, {}, K_);
					}

				}

				void setTime(Real t) { currentTime_ = t; }

				void advance() { linalg::operations::copy(problem_.U(), problem_.U_prev()); }

				void onStepComplete(Index step, Real time) {
					problem_.writeOutput(step, time);
					problem_.evaluateMonitors(step, time);
				}

				decltype(auto) solution() const { return problem_.U(); }

			private:

				bool isSparse() const { return problem_.solverInstance().linear->operatorType == config::LinearSolverConfig::OperatorType::CSR; }
				bool needsK() const { return !solver::isNonlinear(problem_.solverInstance().mode) && isSparse(); }
				bool needsJ() const { return solver::isNonlinear(problem_.solverInstance().mode) && isSparse(); }

				ProblemT& problem_;

				Real currentTime_ = Real(0);
				Real dt_ = Real(0);

				OperatorFormsT operatorForms_;
				OperatorModelT operatorModel_;

				MatrixT K_;
				VectorT massScratch_;
				VectorT Udot_;
				VectorT R_;

				std::unique_ptr<linalg::solver::LinearSolverRunner<VectorT>> linearSolverRunner_;
				std::unique_ptr<nonlinear::NonlinearSolverRunner<VectorT>> nonlinearSolverRunner_;

				// nonlinear-only
				JacobianFormsT jacobianForms_;
				std::unique_ptr<MatrixT> J_;
				std::unique_ptr<VectorT> dU_;
				std::unique_ptr<linalg::solver::LinearSolverRunner<VectorT>> jacobianSolverRunner_;

			}; // class BackwardEulerStage

		} // namespace stage
	} // namespace solver
} // namespace residuum

#endif
