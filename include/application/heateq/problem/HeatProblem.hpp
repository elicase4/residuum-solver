#ifndef RESIDUUM_APPLICATION_HEATEQ_PROBLEM_HEATPROBLEM_HPP
#define RESIDUUM_APPLICATION_HEATEQ_PROBLEM_HEATPROBLEM_HPP

#include <array>
#include <cstring>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "application/heateq/config/HeatConfig.hpp"

#include "core/Types.hpp"

#include "io/visualization/VisualizationWriter.hpp"
#include "io/fieldio/FieldIO.hpp"
#include "io/fieldio/NamedNodalField.hpp"
#include "io/fieldio/NodalFileValueSource.hpp"
#include "io/fieldio/NodalValueSourceAdapter.hpp"

#include "equation/heateq/evaluator/HeatBoundaryStateVector.hpp"
#include "equation/heateq/evaluator/HeatStateVector.hpp"

#include "fem/assembly/Assembler.hpp"
#include "fem/boundary/BoundaryApplicator.hpp"
#include "fem/dispatch/DiscretizationLimits.hpp"
#include "fem/dof/DOFOrdering.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/boundary/NaturalBoundaryRegistry.hpp"
#include "fem/evaluator/ModelRegistry.hpp"
#include "fem/form/FormRegistry.hpp"
#include "fem/form/ScaledForm.hpp"
#include "fem/quantity/BoundaryQuadratureMonitorGroup.hpp"
#include "fem/quantity/BoundaryReactionMonitorGroup.hpp"
#include "fem/quantity/MonitorGroup.hpp"
#include "fem/quantity/QuantityEvaluator.hpp"
#include "fem/quantity/QuantityUnits.hpp"
#include "fem/quantity/recovery/SimpleAverageRecovery.hpp"

#include "utils/logging/core/CsvWriter.hpp"

#include "linalg/operations/VectorOps.hpp"
#include "linalg/operator/CSROperator.hpp"
#include "linalg/operator/FEMOperator.hpp"
#include "linalg/solver/base/LinearSolverRunner.hpp"
#include "linalg/types/CSRMatrix.hpp"
#include "linalg/types/Vector.hpp"

#include "mesh/Mesh.hpp"

#include "solver/SolverInstance.hpp"
#include "solver/config/EvaluationMethod.hpp"
#include "solver/config/NodalFieldReadConfig.hpp"
#include "solver/linear/LinearSolverFactory.hpp"
#include "solver/logging/LoggerFactory.hpp"

#include "topology/TopologicalDOF.hpp"

#include "utils/expression/ScalarExpression.hpp"
#include "utils/expression/VectorExpression.hpp"
#include "utils/logging/driver/Logger.hpp"

namespace residuum {
	namespace application {
		namespace heateq {
			namespace problem {

				template<typename BackendT, typename HeatEqBundleT>
				class HeatProblem {
				public:

					using VectorT = linalg::types::Vector<Real, BackendT>;
					using MatrixT = linalg::types::CSRMatrix<Real, BackendT>;

					static constexpr Index NumDOFs = HeatEqBundleT::NumDOFs;

					using StateVectorT = typename equation::heateq::evaluator::HeatStateVector<HeatEqBundleT::NumDOFs, HeatEqBundleT::SpatialDim, HeatEqBundleT::EvalEle::ParametricDim, VectorT>;
					using BoundaryStateVectorT = typename equation::heateq::evaluator::HeatBoundaryStateVector<HeatEqBundleT::NumDOFs, HeatEqBundleT::SpatialDim, HeatEqBundleT::EvalEle::ParametricDim, VectorT>;

					HeatProblem(const config::HeatConfig& config, mesh::Mesh mesh, typename HeatEqBundleT::Basis basis, typename HeatEqBundleT::QuadratureVolumeType quadratureVolume, typename HeatEqBundleT::QuadratureBoundaryType quadratureBoundary);

					HeatProblem(const HeatProblem&) = delete;
					HeatProblem& operator=(const HeatProblem&) = delete;
					HeatProblem(HeatProblem&&) = delete;
					HeatProblem& operator=(HeatProblem&&) = delete;

					// stateless
					typename HeatEqBundleT::StiffnessForms stiffnessForms() const { return {}; }
					typename HeatEqBundleT::MassForms massForms() const { return {}; }
					typename HeatEqBundleT::TangentDiffusionForms tangentStiffnessForms() const { return {}; }
					typename HeatEqBundleT::TangentMassForms tangentMassForms() const { return {}; }

					const typename HeatEqBundleT::ConductivityModel& stiffnessModel() const { return conductivityModel_; }
					const typename HeatEqBundleT::MassModel& massModel() const { return massModel_; }

					MatrixT createMatrix() const { return fem::assembly::Assembler<BackendT>::template createMatrix<HeatEqBundleT::NumDOFs>(mesh_, topoDOF_); }
					VectorT createVector() const { return fem::assembly::Assembler<BackendT>::template createVector<HeatEqBundleT::NumDOFs>(mesh_, topoDOF_); }

					void assembleLoad(Real time);
					void applyNatural(Real time);

					template<typename FormsT, typename ModelT>
					void assembleMatrix(Real time, const FormsT& forms, const ModelT& model, StateVectorT& state, MatrixT& K);

					// selector reads whichever of state's named buffers (valueBuffer/rateBuffer) forms is meant to be evaluated against
					template<typename FormsT, typename ModelT, typename SelectorT>
					void assembleVector(Real time, const FormsT& forms, const ModelT& model, StateVectorT& state, SelectorT selector, VectorT& V);

					// a residual is always evaluated against the current solution (state.valueBuffer())
					template<typename FormsT, typename ModelT>
					void assembleResidual(Real time, const FormsT& forms, const ModelT& model, StateVectorT& state, VectorT& R);

					template<typename FormsT, typename ModelT>
					void applyEssential(Real time, const FormsT& forms, const ModelT& model, VectorT& F);

					template<typename FormsT, typename ModelT>
					std::unique_ptr<linalg::solver::LinearSolverRunner<VectorT>> makeLinearRunner(const Real* time, const FormsT& forms, const ModelT& model, StateVectorT& state, MatrixT& K);

					Index numFreeDOFs() const { return topoDOF_.numFreeDOFs(); }
					std::vector<std::string> dofNames() const { return {"T"}; }
					Index freeDOFsPerField() const { return topoDOF_.numFreeDOFs() / static_cast<Index>(dofNames().size()); }
					fem::dof::DOFOrdering dofOrdering() const { return topoDOF_.ordering(); }
					const solver::SolverInstance& solverInstance() const { return solverInstance_; }
					std::string equationLabel() const { return "Heat Equation"; }
					const solver::config::LoggingConfig& loggingConfig() const { return config_.logging; }

					VectorT& U() { return *U_; }
					VectorT& F() { return *F_; }
					VectorT& U_prev() { return *U_prev_; }

					VectorT& Udot() { return *Udot_; }

					const VectorT& solution() const { return *U_; }

					void advance() { linalg::operations::copy(*U_, *U_prev_); }

					void writeOutput(Index step, Real time) const;

					void writeLog() const;

					void evaluateMonitors(Index tick, Real time) const;

				private:

					using OpType = solver::config::LinearSolverConfig::OperatorType;

					void loadNodalField(const solver::config::NodalFieldReadConfig& cfg, VectorT& target) const;

					using ExpressionSourceFormsT = typename HeatEqBundleT::template ExpressionSourceForms<utils::expression::ScalarExpression>;
					using NodalSourceFormsT = typename HeatEqBundleT::template NodalSourceForms<typename HeatEqBundleT::NodalScalarSource>;
					using ExpressionFluxFormsT = typename HeatEqBundleT::template ExpressionFluxForms<utils::expression::VectorExpression>;
					using NodalFluxFormsT = typename HeatEqBundleT::template NodalFluxForms<typename HeatEqBundleT::NodalFluxSource>;

					using NegatedStiffnessFormsT = fem::form::ScaledForm<HeatEqBundleT::NumDOFs, typename HeatEqBundleT::StiffnessForms>;
					using NegatedMassFormsT = fem::form::ScaledForm<HeatEqBundleT::NumDOFs, typename HeatEqBundleT::MassForms>;
					using WeakFluxModelT = fem::evaluator::ModelRegistry<typename HeatEqBundleT::ConductivityModel, typename HeatEqBundleT::MassModel>;

					using DirichletExpressionT = typename HeatEqBundleT::template DirichletExpression<utils::expression::ScalarExpression>;
					using DirichletNodalT = typename HeatEqBundleT::template DirichletNodal<typename HeatEqBundleT::NodalScalarSource>;
					using FluxFunctionExpressionT = typename HeatEqBundleT::template FluxBCExpression<utils::expression::VectorExpression>;
					using FluxFunctionNodalT = typename HeatEqBundleT::template FluxBCNodal<typename HeatEqBundleT::NodalFluxSource>;

					struct MonitorOutput {
						std::string name;
						bool toConsole;
						std::string unit;
						utils::logging::CsvWriter csv;
					}; // struct MonitorOutput

					std::pair<std::unique_ptr<fem::quantity::MonitorGroup>, std::string> makeMonitorGroup(config::MonitorConfig::Quantity quantity, fem::quantity::Reduction mode, solver::config::EvaluationMethod fluxMethod) const;

					std::pair<std::unique_ptr<fem::quantity::MonitorGroup>, std::string> makeWeakFluxMonitorGroup(fem::quantity::Reduction mode) const;

					template<typename InnerSourceFormT>
					std::pair<std::unique_ptr<fem::quantity::MonitorGroup>, std::string> makeWeakFluxMonitorGroupFor(const InnerSourceFormT& sourceForm, fem::quantity::Reduction mode) const;

					template<typename FieldT, typename ModelT>
					std::vector<Real> computeNodalField(const ModelT& model, FieldT field, Real time) const;

					template<typename FormT>
					std::pair<std::unique_ptr<fem::quantity::MonitorGroup>, std::string> makeMonitorGroupForForm(fem::quantity::Reduction mode) const;

					template<typename FormT, fem::quantity::Reduction Mode>
					std::pair<std::unique_ptr<fem::quantity::MonitorGroup>, std::string> makeMonitorGroupFor() const;

					config::HeatConfig config_;

					solver::SolverInstance solverInstance_;

					mesh::Mesh mesh_;
					topology::TopologicalDOF<HeatEqBundleT::NumDOFs> topoDOF_;

					utils::logging::driver::Logger driverLogger_;

					typename HeatEqBundleT::EvalEle evalEleTemplate_;
					typename HeatEqBundleT::QuadratureVolumeType quadratureVolume_;
					typename HeatEqBundleT::QuadratureBoundaryType quadratureBoundary_;

					fem::assembly::Assembler<BackendT> assembler_;

					fem::boundary::BoundaryApplicator<BackendT> bcApplicator_;

					fem::boundary::EssentialBoundaryRegistry essentialBCs_;
					fem::boundary::NaturalBoundaryRegistry<typename HeatEqBundleT::EvalQPBdy> naturalBCs_;

					typename HeatEqBundleT::ConductivityModel conductivityModel_;
					typename HeatEqBundleT::ConductivityModelBdy conductivityModelBdy_;

					typename HeatEqBundleT::DefaultModel defaultModel_;
					typename HeatEqBundleT::DefaultModelBdy defaultModelBdy_;

					mutable std::vector<std::unique_ptr<fem::quantity::MonitorGroup>> monitorGroups_;
					std::vector<std::vector<MonitorOutput>> monitorOutputsByGroup_;
					std::vector<std::string> groupUnits_;

					std::optional<ExpressionSourceFormsT> sourceForms_;
					std::optional<NodalSourceFormsT> nodalSourceForms_;

					std::vector<std::unique_ptr<ExpressionFluxFormsT>> expressionFluxForms_;
					std::vector<std::unique_ptr<NodalFluxFormsT>> nodalFluxForms_;

					typename HeatEqBundleT::DensityModel densityModel_;
					typename HeatEqBundleT::SpecificHeatModel specificHeatModel_;
					typename HeatEqBundleT::MassModel massModel_;

					NegatedStiffnessFormsT negatedStiffnessForms_{Real(-1)};
					NegatedMassFormsT negatedMassForms_{Real(-1)};
					WeakFluxModelT weakFluxModel_;

					std::unique_ptr<VectorT> F_;
					std::unique_ptr<VectorT> U_;
					std::unique_ptr<VectorT> U_prev_;
					std::unique_ptr<VectorT> Udot_;

					mutable StateVectorT weakFluxState_;
					mutable BoundaryStateVectorT strongFluxState_;

					mutable std::optional<io::visualization::VisualizationWriter> vizWriter_;

				}; // class HeatProblem

			} // namespace problem
		} // namespace heateq
	} // namespace application
} // namespace residuum

#include "application/heateq/problem/HeatProblem.tpp"

#endif
