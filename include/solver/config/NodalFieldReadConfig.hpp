#ifndef RESIDUUM_SOLVER_CONFIG_NODALFIELDREADCONFIG_HPP
#define RESIDUUM_SOLVER_CONFIG_NODALFIELDREADCONFIG_HPP

#include <string>

namespace residuum {
	namespace solver {
		namespace config {

			struct NodalFieldReadConfig {

				enum class Mode { 
					Expression, 
					File 
				}; // enum class Mode

				Mode mode = Mode::Expression;
				std::string expression;
				std::string file;
				std::string unit;

			}; // struct NodalFieldReadConfig

		} // namespace config
	} // namespace solver
} // namespace residuum

#endif
