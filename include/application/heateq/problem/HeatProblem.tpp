namespace residuum::application::heateq::problem {

	template<typename BackendT, typename HeatEqBundleT>
	HeatProblem<BackendT, HeatEqBundleT>::HeatProblem(const application::heateq::config::HeatConfig& config, mesh::Mesh mesh, typename HeatEqBundleT::Basis basis, typename HeatEqBundleT::QuadratureVolumeType quadratureVolume, typename HeatEqBundleT::QuadratureBoundaryType quadratureBoundary) : config_(config), solverInstance_(solver::resolveSolverInstance(config_.solver)), mesh_(std::move(mesh)), topoDOF_(mesh_, config_.discretization.dofOrdering), driverLogger_(solver::logging::makeDriverLogger(config_.logging.driver, "heateq")), evalEleTemplate_(std::move(basis)), quadratureVolume_(std::move(quadratureVolume)), quadratureBoundary_(std::move(quadratureBoundary)) {

		if (config_.conductivity.type == config::ConductivityConfig::Type::Constant) {
			conductivityModel_.setConstant(config_.conductivity.value);
			conductivityModelBdy_.setConstant(config_.conductivity.value);
		} else if (config_.conductivity.type == config::ConductivityConfig::Type::Anisotropic) {
			conductivityModel_.setAnisotropic(config_.conductivity.tensor);
			conductivityModelBdy_.setAnisotropic(config_.conductivity.tensor);
		} else if (config_.conductivity.type == config::ConductivityConfig::Type::TemperatureDependentIsotropic) {
			conductivityModel_.setTemperatureDependentIsotropic(config_.conductivity.valueExpression, config_.conductivity.gradientExpression);
			conductivityModelBdy_.setTemperatureDependentIsotropic(config_.conductivity.valueExpression, config_.conductivity.gradientExpression);
		} else if (config_.conductivity.type == config::ConductivityConfig::Type::TemperatureDependentAnisotropic) {
			conductivityModel_.setTemperatureDependentAnisotropic(config_.conductivity.tensor, config_.conductivity.valueExpression, config_.conductivity.gradientExpression);
			conductivityModelBdy_.setTemperatureDependentAnisotropic(config_.conductivity.tensor, config_.conductivity.valueExpression, config_.conductivity.gradientExpression);
		} else {
			throw std::runtime_error("HeatProblem: unsupported conductivity type");
		}

		const bool transient = solver::isTransient(solverInstance_.mode);

		// absent 'output:' section means no field-snapshot output at all
		if (config_.output.has_value()) {
			const auto vizFormat = config_.output->format == solver::config::OutputConfig::Format::VTU ? io::visualization::VisualizationWriter::Format::VTU : io::visualization::VisualizationWriter::Format::VTK;
			vizWriter_.emplace(config_.output->directory, config_.output->prefix, vizFormat, transient);
		}

		if (transient) {

			if (!config_.density.has_value() || !config_.specificHeat.has_value()) {
				throw std::runtime_error("HeatProblem: 'materials.density' and 'materials.specific_heat' are required for a transient driver");
			}

			if (!solverInstance_.timestepper) {
				throw std::runtime_error("HeatProblem: solver.timestepper config is required for a transient driver");
			}

			densityModel_.value = config_.density->value;

			if (config_.specificHeat->type == config::SpecificHeatConfig::Type::Constant) {
				specificHeatModel_.setConstant(config_.specificHeat->value);
			} else if (config_.specificHeat->type == config::SpecificHeatConfig::Type::TemperatureDependent) {
				specificHeatModel_.setTemperatureDependent(config_.specificHeat->valueExpression, config_.specificHeat->gradientExpression);
			} else {
				throw std::runtime_error("HeatProblem: unsupported specific heat type");
			}

		}

		massModel_ = typename HeatEqBundleT::MassModel(densityModel_, specificHeatModel_);

		// source
		if (config_.source.read.mode == solver::config::NodalFieldReadConfig::Mode::Expression) {
			sourceForms_.emplace(config_.source.read.expression);
		} else {
			nodalSourceForms_.emplace(typename HeatEqBundleT::NodalScalarSource(mesh_, config_.source.read.file));
		}

		// boundary conditions
		Index numFluxBCs = config_.boundaryConditions.size();
		expressionFluxForms_.reserve(numFluxBCs);
		nodalFluxForms_.reserve(numFluxBCs);

		for (const auto& bcCfg : config_.boundaryConditions) {

			for (auto form : bcCfg.forms) {

				if (form == application::heateq::config::BoundaryConditionConfig::Form::ValueBC) {

					if (bcCfg.mode == solver::config::NodalFieldReadConfig::Mode::File) {

						auto bc = std::shared_ptr<fem::boundary::BoundaryCondition<DirichletNodalT>>(new fem::boundary::BoundaryCondition<DirichletNodalT>{bcCfg.boundaryID, {fem::boundary::BCCategory::Essential}, DirichletNodalT{mesh_, bcCfg.file}});
						essentialBCs_.registerBC<DirichletNodalT>(bc);

					} else {

						auto bc = std::shared_ptr<fem::boundary::BoundaryCondition<DirichletExpressionT>>(new fem::boundary::BoundaryCondition<DirichletExpressionT>{bcCfg.boundaryID, {fem::boundary::BCCategory::Essential}, DirichletExpressionT{bcCfg.expression}});
						essentialBCs_.registerBC<DirichletExpressionT>(bc);

					}

				} else if (form == config::BoundaryConditionConfig::Form::FluxBC) {

					if (bcCfg.mode == solver::config::NodalFieldReadConfig::Mode::File) {

						typename HeatEqBundleT::NodalFluxSource nodalFluxSource(mesh_, bcCfg.file);
						nodalFluxForms_.push_back(std::make_unique<NodalFluxFormsT>(nodalFluxSource));

						auto bc = std::shared_ptr<fem::boundary::BoundaryCondition<FluxFunctionNodalT>>(new fem::boundary::BoundaryCondition<FluxFunctionNodalT>{bcCfg.boundaryID, {fem::boundary::BCCategory::Natural}, FluxFunctionNodalT{mesh_, bcCfg.file}});
						naturalBCs_.template registerBC<FluxFunctionNodalT>(bc, *nodalFluxForms_.back(), defaultModelBdy_);

					} else {

						if (bcCfg.fluxExpression.size() != HeatEqBundleT::SpatialDim) {
							throw std::runtime_error("HeatProblem: flux boundary 'expression' must have SpatialDim components");
						}

						expressionFluxForms_.push_back(std::make_unique<ExpressionFluxFormsT>(bcCfg.fluxExpression));

						auto bc = std::shared_ptr<fem::boundary::BoundaryCondition<FluxFunctionExpressionT>>(new fem::boundary::BoundaryCondition<FluxFunctionExpressionT>{bcCfg.boundaryID, {fem::boundary::BCCategory::Natural}, FluxFunctionExpressionT{bcCfg.fluxExpression}});
						naturalBCs_.template registerBC<FluxFunctionExpressionT>(bc, *expressionFluxForms_.back(), defaultModelBdy_);

					}

				}

			}

		}

		// constraints & linear system containers
		topoDOF_.buildConstraints(evalEleTemplate_.basis(), essentialBCs_);

		F_ = std::make_unique<VectorT>(createVector());
		U_ = std::make_unique<VectorT>(createVector());
		U_->zero();

		// U_prev for transient solver
		U_prev_ = std::make_unique<VectorT>(createVector());

		// monitors
		{
			std::map<std::pair<config::MonitorConfig::Quantity, fem::quantity::Reduction>, Index> groupIndexForKey;

			for (const auto& monitorCfg : config_.monitors) {

				const auto key = std::make_pair(monitorCfg.quantity, monitorCfg.reduction);

				Index groupIdx;
				auto it = groupIndexForKey.find(key);
				if (it != groupIndexForKey.end()) {

					groupIdx = it->second;

				} else {

					groupIdx = static_cast<Index>(monitorGroups_.size());
					auto [group, unit] = makeMonitorGroup(monitorCfg.quantity, monitorCfg.reduction);
					monitorGroups_.push_back(std::move(group));
					monitorOutputsByGroup_.emplace_back();
					groupUnits_.push_back(unit);
					groupIndexForKey.emplace(key, groupIdx);

				}

				fem::quantity::MonitorGroup& group = *monitorGroups_[groupIdx];
				const Index combinationIdx = group.addCombination();

				for (const auto& term : monitorCfg.terms) {
					group.registerTag(term.boundary);
					group.addTerm(combinationIdx, term.boundary, term.coefficient);
				}

				const std::string& unit = groupUnits_[groupIdx];
				const std::vector<std::string> columns{"tick", "time", monitorCfg.name + " [" + unit + "]"};
				monitorOutputsByGroup_[groupIdx].push_back(MonitorOutput{monitorCfg.name, monitorCfg.output.console, unit, io::visualization::VisualizationWriter::makeCsvWriter(monitorCfg.output.file, columns)});

			}
		}

		// initial condition
		loadNodalField(config_.initialCondition.read, *U_);
		if (transient) {
			linalg::operations::copy(*U_, *U_prev_);
		}

	}

	template<typename BackendT, typename HeatEqBundleT>
	void HeatProblem<BackendT, HeatEqBundleT>::loadNodalField(const solver::config::NodalFieldReadConfig& cfg, VectorT& target) const {

		auto scatter = [&](Index nodeID, const Real* val) {
			for (Index c = 0; c < HeatEqBundleT::NumDOFs; ++c) {
				Index tdof = topoDOF_.getNodeDOF(nodeID, c);
				if (topoDOF_.isConstrained(tdof)) continue;
				target.data()[topoDOF_.toAlgebraic(tdof)] = val[c];
			}
		};

		if (cfg.mode == solver::config::NodalFieldReadConfig::Mode::Expression) {

			utils::expression::ScalarExpression expr(cfg.expression);

			for (Index nodeID = 0; nodeID < mesh_.data.numNodes; ++nodeID) {
				Real coords[3] = {Real(0), Real(0), Real(0)};
				const Real* p = mesh_.getNodeCoord(nodeID);
				for (Index d = 0; d < HeatEqBundleT::SpatialDim; ++d) coords[d] = p[d];
				Real val[HeatEqBundleT::NumDOFs];
				expr(Real(0), coords, val);
				scatter(nodeID, val);
			}

		} else {

			typename HeatEqBundleT::NodalScalarSource src(mesh_, cfg.file);

			for (Index nodeID = 0; nodeID < mesh_.data.numNodes; ++nodeID) {
				Real val[HeatEqBundleT::NumDOFs];
				src.eval(nodeID, val);
				scatter(nodeID, val);
			}

		}

	}

	template<typename BackendT, typename HeatEqBundleT>
	void HeatProblem<BackendT, HeatEqBundleT>::assembleLoad(Real time) {

		F_->zero();

		if (config_.source.read.mode == solver::config::NodalFieldReadConfig::Mode::Expression) {
			fem::assembly::Assembler<BackendT>::template assembleVector<HeatEqBundleT::NumDOFs, typename HeatEqBundleT::EvalEle, typename HeatEqBundleT::EvalQPVol, typename HeatEqBundleT::DefaultModel, ExpressionSourceFormsT, typename HeatEqBundleT::QuadratureVolumeType, fem::assembly::GatherMode::Free>(mesh_, topoDOF_, time, defaultModel_, *sourceForms_, evalEleTemplate_, quadratureVolume_, *U_, nullptr, {nullptr}, *F_, nullptr);
		} else {
			fem::assembly::Assembler<BackendT>::template assembleVector<HeatEqBundleT::NumDOFs, typename HeatEqBundleT::EvalEle, typename HeatEqBundleT::EvalQPVol, typename HeatEqBundleT::DefaultModel, NodalSourceFormsT, typename HeatEqBundleT::QuadratureVolumeType, fem::assembly::GatherMode::Free>(mesh_, topoDOF_, time, defaultModel_, *nodalSourceForms_, evalEleTemplate_, quadratureVolume_, *U_, nullptr, {nullptr}, *F_, nullptr);
		}

	}

	template<typename BackendT, typename HeatEqBundleT>
	void HeatProblem<BackendT, HeatEqBundleT>::applyNatural(Real time) {

		bcApplicator_.template applyNaturalBCs<HeatEqBundleT::NumDOFs, typename HeatEqBundleT::EvalEle, typename HeatEqBundleT::EvalQPBdy, typename HeatEqBundleT::QuadratureBoundaryType>(mesh_, topoDOF_, naturalBCs_, time, evalEleTemplate_, quadratureBoundary_, *F_);

	}

	template<typename BackendT, typename HeatEqBundleT>
	template<fem::assembly::GatherMode Mode, typename FormsT, typename ModelT>
	void HeatProblem<BackendT, HeatEqBundleT>::assembleMatrix(Real time, const FormsT& forms, const ModelT& model, const std::array<const VectorT*, NumAuxStates>& auxStates, MatrixT& K) {

		fem::assembly::Assembler<BackendT>::template assembleMatrix<HeatEqBundleT::NumDOFs, typename HeatEqBundleT::EvalEle, typename HeatEqBundleT::EvalQPVol, ModelT, FormsT, typename HeatEqBundleT::QuadratureVolumeType, Mode>(mesh_, topoDOF_, time, model, forms, evalEleTemplate_, quadratureVolume_, *U_, auxStates, K, &essentialBCs_);

	}

	template<typename BackendT, typename HeatEqBundleT>
	template<fem::assembly::GatherMode Mode, typename FormsT, typename ModelT>
	void HeatProblem<BackendT, HeatEqBundleT>::assembleVector(Real time, const FormsT& forms, const ModelT& model, const VectorT& gatherSource, const std::array<const VectorT*, NumAuxStates>& auxStates, VectorT& V, const VectorT* fieldSource) {

		fem::assembly::Assembler<BackendT>::template assembleVector<HeatEqBundleT::NumDOFs, typename HeatEqBundleT::EvalEle, typename HeatEqBundleT::EvalQPVol, ModelT, FormsT, typename HeatEqBundleT::QuadratureVolumeType, Mode>(mesh_, topoDOF_, time, model, forms, evalEleTemplate_, quadratureVolume_, gatherSource, fieldSource, auxStates, V, &essentialBCs_);

	}

	template<typename BackendT, typename HeatEqBundleT>
	template<fem::assembly::GatherMode Mode, typename FormsT, typename ModelT>
	void HeatProblem<BackendT, HeatEqBundleT>::assembleResidual(Real time, const FormsT& forms, const ModelT& model, const std::array<const VectorT*, NumAuxStates>& auxStates, VectorT& R) {

		assembleVector<Mode>(time, forms, model, *U_, auxStates, R);

	}

	template<typename BackendT, typename HeatEqBundleT>
	template<typename FormsT, typename ModelT>
	void HeatProblem<BackendT, HeatEqBundleT>::applyEssential(Real time, const FormsT& forms, const ModelT& model, VectorT& F) {

		bcApplicator_.template applyEssentialBCs<HeatEqBundleT::NumDOFs, typename HeatEqBundleT::EvalEle, typename HeatEqBundleT::EvalQPVol, ModelT, FormsT, typename HeatEqBundleT::QuadratureVolumeType>(mesh_, topoDOF_, essentialBCs_, time, model, forms, evalEleTemplate_, quadratureVolume_, F);

	}

	template<typename BackendT, typename HeatEqBundleT>
	template<fem::assembly::GatherMode Mode, typename FormsT, typename ModelT>
	std::unique_ptr<linalg::solver::LinearSolverRunner<typename HeatProblem<BackendT, HeatEqBundleT>::VectorT>> HeatProblem<BackendT, HeatEqBundleT>::makeLinearRunner(const Real* time, const FormsT& forms, const ModelT& model, const VectorT* fieldSource, const std::array<const VectorT*, NumAuxStates>& auxStates, MatrixT& K) {

		using CSROperatorT = linalg::op::CSROperator<MatrixT>;
		using MatrixFreeOperatorT = linalg::op::FEMOperator<fem::assembly::Assembler<BackendT>, topology::TopologicalDOF<HeatEqBundleT::NumDOFs>, typename HeatEqBundleT::EvalEle, typename HeatEqBundleT::EvalQPVol, ModelT, FormsT, typename HeatEqBundleT::QuadratureVolumeType, Mode, VectorT>;

		const OpType opType = solverInstance_.linear->operatorType;

		if (opType == OpType::CSR) {
			return solver::linear::makeLinearSolverRunner<CSROperatorT, VectorT>(CSROperatorT(K), topoDOF_.numFreeDOFs(), *solverInstance_.linear, config_.logging.linear, equationLabel(), dofNames(), freeDOFsPerField(), dofOrdering());
		} else if (opType == OpType::FEM) {
			return solver::linear::makeLinearSolverRunner<MatrixFreeOperatorT, VectorT>(MatrixFreeOperatorT(assembler_, mesh_, topoDOF_, time, model, forms, evalEleTemplate_, quadratureVolume_, &essentialBCs_, fieldSource, auxStates), topoDOF_.numFreeDOFs(), *solverInstance_.linear, config_.logging.linear, equationLabel(), dofNames(), freeDOFsPerField(), dofOrdering());
		}

		throw std::runtime_error("HeatProblem: unsupported operator");

	}

	template<typename BackendT, typename HeatEqBundleT>
	void HeatProblem<BackendT, HeatEqBundleT>::writeOutput(Index step, Real time) const {

		if (!vizWriter_.has_value() || (step % config_.output->writeFrequency != 0)) return;

		const std::string filename = vizWriter_->template writeField<HeatEqBundleT::NumDOFs>(mesh_, topoDOF_, essentialBCs_, step, time, U_->data(), {"T"}, {"K"});

		const char* ordStr = (topoDOF_.ordering() == fem::dof::DOFOrdering::Interleaved) ? "Interleaved" : "Block";
		driverLogger_.event("wrote '" + filename + "' - " + std::to_string(HeatEqBundleT::NumDOFs) + " field(s), " + std::to_string(mesh_.data.numNodes) + " nodes, " + ordStr + " ordering");

	}

	template<typename BackendT, typename HeatEqBundleT>
	void HeatProblem<BackendT, HeatEqBundleT>::writeLog() const {

	}

	template<typename BackendT, typename HeatEqBundleT>
	template<typename FormT, fem::quantity::Reduction Mode>
	std::pair<std::unique_ptr<fem::quantity::MonitorGroup>, std::string> HeatProblem<BackendT, HeatEqBundleT>::makeMonitorGroupFor() const {

		using QuantityFormsT = fem::quantity::QuantityForms<fem::quantity::ReducedQuantity<FormT, Mode>>;
		using GroupT = fem::quantity::MonitorGroupImpl<BackendT, HeatEqBundleT::NumDOFs, typename HeatEqBundleT::EvalEle, typename HeatEqBundleT::EvalQPBdy, typename HeatEqBundleT::ConductivityModelBdy, QuantityFormsT, typename HeatEqBundleT::QuadratureBoundaryType>;

		auto group = std::make_unique<GroupT>(mesh_, topoDOF_, essentialBCs_, conductivityModelBdy_, QuantityFormsT{}, evalEleTemplate_, quadratureBoundary_, *U_);
		const std::string unit = fem::quantity::unitFor<FormT>(Mode, HeatEqBundleT::SpatialDim - 1);

		return {std::move(group), unit};

	}

	template<typename BackendT, typename HeatEqBundleT>
	template<typename FormT>
	std::pair<std::unique_ptr<fem::quantity::MonitorGroup>, std::string> HeatProblem<BackendT, HeatEqBundleT>::makeMonitorGroupForForm(fem::quantity::Reduction mode) const {

		switch (mode) {
			case fem::quantity::Reduction::Integral: return makeMonitorGroupFor<FormT, fem::quantity::Reduction::Integral>();
			case fem::quantity::Reduction::Average: return makeMonitorGroupFor<FormT, fem::quantity::Reduction::Average>();
		}

		throw std::runtime_error("HeatProblem: unknown monitor reduction mode");

	}

	template<typename BackendT, typename HeatEqBundleT>
	std::pair<std::unique_ptr<fem::quantity::MonitorGroup>, std::string> HeatProblem<BackendT, HeatEqBundleT>::makeMonitorGroup(config::MonitorConfig::Quantity quantity, fem::quantity::Reduction mode) const {

		switch (quantity) {
			case config::MonitorConfig::Quantity::HeatFlux:
				return makeMonitorGroupForForm<typename HeatEqBundleT::HeatFluxIntegrand>(mode);
		}

		throw std::runtime_error("HeatProblem: unknown monitor quantity");

	}

	template<typename BackendT, typename HeatEqBundleT>
	void HeatProblem<BackendT, HeatEqBundleT>::evaluateMonitors(Index tick, Real time) const {

		for (Index g = 0; g < static_cast<Index>(monitorGroups_.size()); ++g) {

			monitorGroups_[g]->evaluate(time, [&](Index combinationIdx, const Real* value) {

				const auto& out = monitorOutputsByGroup_[g][combinationIdx];
				if (out.toConsole) driverLogger_.event("monitor '" + out.name + "' = " + std::to_string(value[0]) + " " + out.unit);
				out.csv.writeRow({static_cast<Real>(tick), time, value[0]});

			});

		}

	}

} // namespace residuum::application::heateq::problem
