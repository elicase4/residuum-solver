#ifndef RESIDUUM_APPLICATION_CLI_CLIAPP_HPP
#define RESIDUUM_APPLICATION_CLI_CLIAPP_HPP

#include "application/cli/AppUsage.hpp"
#include "application/cli/Subcommand.hpp"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace residuum {
	namespace application {
		namespace cli {

			inline void printTopLevelUsage(const std::vector<Subcommand>& subcommands) {

				std::cout << "Usage: residuum <subcommand> <config.yaml>\n\n";
				std::cout << "Subcommands:\n";

				for (const auto& cmd : subcommands) {
					std::cout << "  " << cmd.name << "\t" << cmd.description << "\n";
				}

				std::cout << "\nRun 'residuum <subcommand> --help' for subcommand-specific usage.\n";

			}

			// argv[1] selects the subcommand; the remaining argv (shifted by one) is forwarded to
			// its run(), so a subcommand's own argv[1] is unchanged -- still its config path
			inline int runSubcommandDispatch(int argc, char** argv, const std::vector<Subcommand>& subcommands) {

				if (argc < 2) {
					printTopLevelUsage(subcommands);
					return EXIT_FAILURE;
				}

				const std::string first = argv[1];

				if (first == "-h" || first == "--help") {
					printTopLevelUsage(subcommands);
					return EXIT_SUCCESS;
				}

				if (first == "--version") {
					printVersion();
					return EXIT_SUCCESS;
				}

				for (const auto& cmd : subcommands) {
					if (cmd.name == first) {
						return cmd.run(argc - 1, argv + 1);
					}
				}

				std::cerr << "residuum: unknown subcommand '" << first << "'\n\n";
				printTopLevelUsage(subcommands);
				return EXIT_FAILURE;

			}

		} // namespace cli
	} // namespace application
} // namespace residuum

#endif
