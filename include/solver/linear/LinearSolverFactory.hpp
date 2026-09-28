#ifndef RESIDUUM_SOLVER_LINEAR_LINEARSOLVERFACTORY_HPP
#define RESIDUUM_SOLVER_LINEAR_LINEARSOLVERFACTORY_HPP

#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "core/Types.hpp"

#include "fem/dof/DOFOrdering.hpp"

#include "solver/config/LinearSolverConfig.hpp"
#include "solver/config/LoggingConfig.hpp"
#include "solver/logging/LoggerFactory.hpp"

#include "linalg/solver/base/LinearSolverRunner.hpp"
#include "linalg/solver/iterative/cg/CGRunner.hpp"
#include "linalg/solver/iterative/cg/Config.hpp"
#include "linalg/solver/iterative/gmres/GMRESRunner.hpp"
#include "linalg/solver/iterative/gmres/Config.hpp"

#include "linalg/solver/preconditioner/Identity.hpp"
#include "linalg/solver/preconditioner/Jacobi.hpp"

#include "utils/logging/core/NullLogger.hpp"
#include "utils/logging/linear/ConsoleLogger.hpp"
#include "utils/logging/linear/Logger.hpp"

namespace residuum {
	namespace solver {
		namespace linear {

			template<typename OperatorT, typename VectorT, typename PreconditionerT>
			std::unique_ptr<linalg::solver::LinearSolverRunner<VectorT>> makeLinearSolverRunnerWithPreconditioner(const OperatorT& op, Index n, const config::LinearSolverConfig& cfg, const config::LinearLoggerConfig& loggerCfg, const std::string& equationName, const std::string& preconditionerName, const std::vector<std::string>& dofNames, Index freeDOFsPerField, fem::dof::DOFOrdering ordering) {

				using LoggerT = utils::logging::linear::Logger;

				switch (cfg.type) {

					case config::LinearSolverConfig::Type::CG: {

						using CGConfigT = linalg::solver::iterative::cg::Config<VectorT>;
						CGConfigT cgCfg;
						cgCfg.tol = static_cast<typename CGConfigT::DataType>(cfg.tolerance);
						cgCfg.maxIters = cfg.maxIterations;

						std::ostringstream tolStream;
						tolStream << std::scientific << std::setprecision(1) << cfg.tolerance;

						std::vector<std::pair<std::string, std::string>> extraParams = {
							{"Maximum Iterations", std::to_string(cfg.maxIterations)},
							{"Tolerance", tolStream.str() + " (" + (cgCfg.tolType == linalg::solver::iterative::cg::ToleranceType::Relative ? "Relative" : "Absolute") + ")"}
						};

						LoggerT logger = logging::makeLinearLogger(loggerCfg, equationName, "CG", preconditionerName, dofNames, extraParams, freeDOFsPerField, ordering);

						return std::make_unique<linalg::solver::iterative::cg::CGRunner<OperatorT, VectorT, PreconditionerT, LoggerT>>(op, n, cgCfg, std::move(logger));

					}

					case config::LinearSolverConfig::Type::GMRES: {

						const auto& p = std::get<config::GMRESParams>(cfg.params);

						using GMRESConfigT = linalg::solver::iterative::gmres::Config<VectorT>;
						GMRESConfigT gmresCfg;
						gmresCfg.tol = static_cast<typename GMRESConfigT::DataType>(cfg.tolerance);
						gmresCfg.maxIters = cfg.maxIterations;
						gmresCfg.krylovDim = p.krylovDim;

						std::ostringstream tolStream;
						tolStream << std::scientific << std::setprecision(1) << cfg.tolerance;

						std::vector<std::pair<std::string, std::string>> extraParams = {
							{"Maximum Iterations", std::to_string(cfg.maxIterations)},
							{"Krylov Dimension", std::to_string(p.krylovDim)},
							{"Tolerance", tolStream.str() + " (" + (gmresCfg.tolType == linalg::solver::iterative::gmres::ToleranceType::Relative ? "Relative" : "Absolute") + ")"}
						};

						LoggerT logger = logging::makeLinearLogger(loggerCfg, equationName, "GMRES", preconditionerName, dofNames, extraParams, freeDOFsPerField, ordering);

						return std::make_unique<linalg::solver::iterative::gmres::GMRESRunner<OperatorT, VectorT, PreconditionerT, LoggerT>>(op, n, gmresCfg, std::move(logger));

					}

					case config::LinearSolverConfig::Type::BiCGSTAB:
						throw std::runtime_error("LinearSolverFactory: BiCGSTAB not yet implemented");

					case config::LinearSolverConfig::Type::LU:
						throw std::runtime_error("LinearSolverFactory: LU not yet implemented");

				}

				throw std::runtime_error("LinearSolverFactory: unknown linear solver type");

			}

			template<typename OperatorT, typename VectorT>
			std::unique_ptr<linalg::solver::LinearSolverRunner<VectorT>> makeLinearSolverRunner(const OperatorT& op, Index n, const config::LinearSolverConfig& cfg, const config::LinearLoggerConfig& loggerCfg, const std::string& equationName, const std::vector<std::string>& dofNames, Index freeDOFsPerField, fem::dof::DOFOrdering ordering) {

				switch (cfg.preconditioner.type) {

					case config::PreconditionerConfig::Type::Identity:
						return makeLinearSolverRunnerWithPreconditioner<OperatorT, VectorT, linalg::solver::preconditioner::Identity<VectorT>>(op, n, cfg, loggerCfg, equationName, "Identity", dofNames, freeDOFsPerField, ordering);

					case config::PreconditionerConfig::Type::Jacobi:
						return makeLinearSolverRunnerWithPreconditioner<OperatorT, VectorT, linalg::solver::preconditioner::Jacobi<VectorT>>(op, n, cfg, loggerCfg, equationName, "Jacobi", dofNames, freeDOFsPerField, ordering);

				}

				throw std::runtime_error("LinearSolverFactory: unknown preconditioner type");

			}

		} // namespace linear
	} // namespace solver
} // namespace residuum

#endif
