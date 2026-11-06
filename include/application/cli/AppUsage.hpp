#ifndef RESIDUUM_APPLICATION_CLI_APPUSAGE_HPP
#define RESIDUUM_APPLICATION_CLI_APPUSAGE_HPP

#include "core/Version.hpp"

#include <iostream>
#include <string>

namespace residuum {
	namespace application {
		namespace cli {

			// shared by each subcommand's own -h/--help, and by the top-level residuum dispatcher
			inline void printAppUsage(const std::string& subcommand, const std::string& description) {

				std::cout << description << "\n\n";
				std::cout << "Usage: residuum " << subcommand << " <config.yaml>\n";

			}

			inline void printVersion() {

				std::cout << "residuum " << RESIDUUM_VERSION << "\n";

			}

		} // namespace cli
	} // namespace application
} // namespace residuum

#endif
