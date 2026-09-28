#ifndef RESIDUUM_UTILS_LOGGING_CORE_ANSICOLOR_HPP
#define RESIDUUM_UTILS_LOGGING_CORE_ANSICOLOR_HPP

#include <cstdio>
#include <string>
#include <unistd.h>

namespace residuum {
	namespace utils {
		namespace logging {

			enum class Color {
				Default,
				Red,
				Yellow,
				Green,
				Cyan,
				Bold
			};

			inline bool colorEnabled() {
				static const bool enabled = (isatty(fileno(stdout)) != 0);
				return enabled;
			}

			inline const char* colorCode(Color c) {

				switch (c) {
					case Color::Red:    return "\033[31m";
					case Color::Yellow: return "\033[33m";
					case Color::Green:  return "\033[32m";
					case Color::Cyan:   return "\033[36m";
					case Color::Bold:   return "\033[1m";
					default:            return "";
				}

			}

			inline std::string colorize(const std::string& text, Color c) {

				if (!colorEnabled() || c == Color::Default) return text;
				return std::string(colorCode(c)) + text + "\033[0m";

			}

		} // namespace logging
	} // namespace utils
} // namespace residuum

#endif
