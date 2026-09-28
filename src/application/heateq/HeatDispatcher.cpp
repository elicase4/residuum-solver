#include <fstream>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include "application/heateq/HeatDispatcher.hpp"
#include "application/heateq/problem/HeatProblem.hpp"

#include "equation/heateq/HeatEquation.hpp"

#include "fem/dispatch/DiscretizationDispatch.hpp"

#include "io/MeshIO.hpp"

#include "linalg/types/backend/CPU.hpp"

#include "mesh/Mesh.hpp"

#include "solver/driver/Steady.hpp"
#include "solver/driver/Transient.hpp"
#include "solver/logging/LoggerFactory.hpp"
#include "solver/stage/SteadyStage.hpp"
#include "solver/timestepper/TimeStepperFactory.hpp"

bool residuum::application::heateq::HeatDispatcher::run(const residuum::application::heateq::config::HeatConfig& config) {

	namespace sconfig = residuum::solver::config;

	if (config.backend.type != sconfig::BackendConfig::Type::CPU) {
		throw std::runtime_error("HeatDispatcher: only the CPU backend is supported so far");
	}

	if (!config.solver.linear.has_value() && !config.solver.nonlinear.has_value()) {
		throw std::runtime_error("HeatDispatcher: solver.linear or solver.nonlinear config is required");
	}

	if (config.solver.driver.type == sconfig::DriverConfig::Type::PseudoTransient) {
		throw std::runtime_error("HeatDispatcher: the pseudo_transient driver is not yet implemented");
	}

	const bool steady = (config.solver.driver.type == sconfig::DriverConfig::Type::Steady);

	for (const std::string& path : {config.logging.driver.textFile, config.logging.linear.textFile}) {
		if (!path.empty()) std::ofstream(path, std::ios::trunc);
	}

	residuum::mesh::Mesh mesh;
	residuum::io::MeshIO::readBinary(mesh, config.mesh.file);

	if (mesh.data.basisType != residuum::mesh::BasisType::Lagrange) {
		throw std::runtime_error("HeatDispatcher: only the Lagrange basis is supported so far");
	}

	const Index nsd = mesh.data.spatialDim;
	const Index npd = mesh.data.parametricDim;
	const auto family = mesh.data.elementFamily;

	const Index basisOrderX = mesh.data.basisOrder.size() > 0 ? mesh.data.basisOrder[0] : 1;
	const Index basisOrderY = mesh.data.basisOrder.size() > 1 ? mesh.data.basisOrder[1] : 1;
	const Index basisOrderZ = mesh.data.basisOrder.size() > 2 ? mesh.data.basisOrder[2] : 1;

	const auto logger = residuum::solver::logging::makeDriverLogger(config.logging.driver, "heateq");

	logger.event("mesh: " + config.mesh.file + " -> " + std::to_string(mesh.data.numElements) + " elements, nsd=" + std::to_string(nsd) + " npd=" + std::to_string(npd));

	bool converged = false;

	auto visitor = [&]<Index NSD, Index NPD, residuum::mesh::ElementFamily Family>(auto basis, auto quadVol, auto quadBdy) -> bool {

		using HeatEqBundle = residuum::equation::HeatEquation<NSD, NPD, Family>;
		using ProblemT = residuum::application::heateq::problem::HeatProblem<residuum::linalg::types::backend::CPU, HeatEqBundle>;

		ProblemT heatProblem(config, std::move(mesh), std::move(basis), std::move(quadVol), std::move(quadBdy));

		if (steady) {

			using StageT = residuum::solver::stage::SteadyStage<ProblemT>;
			StageT stage(heatProblem);

			residuum::solver::driver::Steady<StageT> driver;
			converged = driver.solve(stage);

			heatProblem.writeOutput(0, 0.0);
			heatProblem.evaluateMonitors(0, 0.0);

		} else {

			const auto& tsCfg = *heatProblem.solverInstance().timestepper;

			// step 0 == the initial condition
			heatProblem.writeOutput(0, tsCfg.t0);
			heatProblem.evaluateMonitors(0, tsCfg.t0);

			converged = residuum::solver::timestepper::dispatchTransientStage(heatProblem, tsCfg, config.logging.timestepper, heatProblem.equationLabel(), [](auto& stage, auto& stepper) {

				using StageT = std::decay_t<decltype(stage)>;
				residuum::solver::driver::Transient<StageT> driver;
				return driver.solve(stage, stepper);

			});

		}

		heatProblem.writeLog();

		if (converged) {
			logger.event("converged");
		} else {
			logger.warn("did not converge");
		}

		return converged;

	};

	const auto& d = config.discretization;
	const bool matched = residuum::fem::dispatch::dispatch(nsd, npd, family, basisOrderX, basisOrderY, basisOrderZ, d.quadrature.xi, d.quadrature.eta, d.quadrature.zeta, visitor);

	if (!matched) {
		throw std::runtime_error("HeatDispatcher: unsupported (nsd, npd, basis order, quadrature order) combination for this mesh/config");
	}

	return converged;

}
