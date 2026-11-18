#ifndef RESIDUUM_APPLICATION_CLI_SUBCOMMAND_HPP
#define RESIDUUM_APPLICATION_CLI_SUBCOMMAND_HPP

#include <string>

namespace residuum {
	namespace application {
		namespace cli {

			struct Subcommand {

				std::string name;
				std::string description;
				int (*run)(int argc, char** argv);

			}; // struct Subcommand

		} // namespace cli
	} // namespace application
} // namespace residuum

#endif
