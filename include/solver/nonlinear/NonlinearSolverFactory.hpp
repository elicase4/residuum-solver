#ifndef RESIDUUM_SOLVER_NONLINEAR_NONLINEARSOLVERFACTORY_HPP
#define RESIDUUM_SOLVER_NONLINEAR_NONLINEARSOLVERFACTORY_HPP

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "fem/dof/DOFOrdering.hpp"

#include "solver/config/LoggingConfig.hpp"
#include "solver/config/NonlinearSolverConfig.hpp"
#include "solver/logging/LoggerFactory.hpp"
#include "solver/nonlinear/Newton.hpp"
#include "solver/nonlinear/NonlinearSolverRunner.hpp"
#include "solver/stage/NonlinearCapableStage.hpp"

#include "utils/logging/nonlinear/Logger.hpp"

namespace residuum {
	namespace solver {
		namespace nonlinear {

			template<stage::NonlinearCapableStage StageT, typename VectorT>
			std::unique_ptr<NonlinearSolverRunner<VectorT>> makeNonlinearSolverRunner(StageT& stage, const config::NonlinearSolverConfig& cfg, const config::NonlinearLoggerConfig& loggerCfg, const std::string& equationLabel, std::vector<std::string> dofNames, Index freeDOFsPerField, fem::dof::DOFOrdering dofOrdering) {

				switch (cfg.type) {

					case config::NonlinearSolverConfig::Type::Newton:
						return std::make_unique<NewtonRunner<StageT, VectorT>>(stage, cfg, logging::makeNonlinearLogger(loggerCfg, equationLabel, "Newton", dofNames), dofNames, freeDOFsPerField, dofOrdering);

					case config::NonlinearSolverConfig::Type::Picard:
						throw std::runtime_error("NonlinearSolverFactory[" + equationLabel + "]: Picard solver not yet implemented");

				}

				throw std::runtime_error("NonlinearSolverFactory[" + equationLabel + "]: unknown nonlinear solver type");

			}

		} // namespace nonlinear
	} // namespace solver
} // namespace residuum

#endif
