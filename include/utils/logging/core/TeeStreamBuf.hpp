#ifndef RESIDUUM_UTILS_LOGGING_CORE_TEESTREAMBUF_HPP
#define RESIDUUM_UTILS_LOGGING_CORE_TEESTREAMBUF_HPP

#include <ostream>
#include <streambuf>
#include <vector>

namespace residuum {
	namespace utils {
		namespace logging {

			class TeeStreamBuf : public std::streambuf {
			public:

				void addTarget(std::streambuf* buf) { targets_.push_back(buf); }

			protected:

				int overflow(int ch) override {
					if (ch == traits_type::eof()) return traits_type::not_eof(ch);
					for (auto* t : targets_) t->sputc(static_cast<char>(ch));
					return ch;
				}

				std::streamsize xsputn(const char* s, std::streamsize n) override {
					for (auto* t : targets_) t->sputn(s, n);
					return n;
				}

			private:

				std::vector<std::streambuf*> targets_;

			}; // class TeeStreamBuf

		} // namespace logging
	} // namespace utils
} // namespace residuum

#endif
