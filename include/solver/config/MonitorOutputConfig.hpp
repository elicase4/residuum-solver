#ifndef RESIDUUM_SOLVER_CONFIG_MONITOROUTPUTCONFIG_HPP
#define RESIDUUM_SOLVER_CONFIG_MONITOROUTPUTCONFIG_HPP

#include <string>

namespace residuum {
	namespace solver {
		namespace config {

			struct MonitorOutputConfig {

				bool console = true;
				std::string file; // optional CSV path

			}; // struct MonitorOutputConfig

		} // namespace config
	} // namespace solver
} // namespace residuum

#endif
