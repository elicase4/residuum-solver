#ifndef RESIDUUM_APPLICATION_HEATEQ_CONFIG_BOUNDARYCONDITIONCONFIG_HPP
#define RESIDUUM_APPLICATION_HEATEQ_CONFIG_BOUNDARYCONDITIONCONFIG_HPP

#include <string>
#include <vector>

#include "core/Types.hpp"
#include "application/heateq/config/ConductivityConfig.hpp"
#include "solver/config/NodalFieldReadConfig.hpp"

namespace residuum {
	namespace application {
		namespace heateq {
			namespace config {

				struct BoundaryConditionConfig {

					// Essential (Dirichlet) vs natural (flux)
					enum class Type {
						Value,
						Flux
					}; // enum class Type

					enum class Form {
						FluxBC,
						ValueBC
					}; // enum class Form

					Int boundaryID;

					Type type;

					// reuses NodalFieldReadConfig::Mode only; Flux needs a per-component vector expression, not its scalar shape
					solver::config::NodalFieldReadConfig::Mode mode = solver::config::NodalFieldReadConfig::Mode::Expression;

					// mode == Expression: Value takes a single scalar expression, Flux takes
					// one expression per spatial component. mode == File: file, regardless of type.
					std::string expression;
					std::vector<std::string> fluxExpression;
					std::string file;
					std::string unit; // required, SI unit label, e.g. "K", "W/m^2"

					std::vector<Form> forms;

					ConductivityConfig::Type model;

				}; // struct BoundaryConditionConfig

			} // namespace config
		} // namespace heateq
	} // namespace application
} // namespace residuum

#endif
