#ifndef RESIDUUM_UTILS_LOGGING_BANNER_HPP
#define RESIDUUM_UTILS_LOGGING_BANNER_HPP

#include <iostream>
#include <string>

namespace residuum {
	namespace utils {
		namespace logging {

			// prints unconditionally since no LoggingConfig exists yet at this point in main()
			inline void printStartupBanner(const std::string& appName, const std::string& description) {

				const int width = 60;

				std::cout << "\n";
				std::cout << std::string(width, '=') << "\n";
				std::cout << "  Residuum: " << appName << "\n";
				std::cout << "  " << description << "\n";
				std::cout << std::string(width, '=') << "\n";
				std::cout << "\n";

			}

		} // namespace logging
	} // namespace utils
} // namespace residuum

#endif
