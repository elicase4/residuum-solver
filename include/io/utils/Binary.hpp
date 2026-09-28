#ifndef RESIDUUM_IO_UTILS_BINARY_HPP
#define RESIDUUM_IO_UTILS_BINARY_HPP

#include <cassert>
#include <cstdint>
#include <cstring>
#include <fstream>

#include "core/Types.hpp"

namespace residuum {
	namespace io {
		namespace binary {
			
			inline bool hostIsLittleEndian() {
				const uint32_t probe = 1u;
				uint8_t byte0;
				std::memcpy(&byte0, &probe, 1);
				return (byte0 == 1u);
			}

			template<typename T>
			T swapBytes(T val){
				
				static_assert((sizeof(T) == 4 || sizeof(T) == 8), "swapBytes: only 4-byte or 8-byte types are supported");

				T out;
				const char* src = reinterpret_cast<const char*>(&val);
				char* dst = reinterpret_cast<char*>(&out);
				for (std::size_t i = 0; i < sizeof(T); ++i){
					dst[i] = src[sizeof(T) - 1 - i];
				}

				return out;

			}
			
			template<typename T>
			void writeBE(std::ofstream& ofs, T val){
				if (hostIsLittleEndian()) val = swapBytes<T>(val);
				ofs.write(reinterpret_cast<const char*>(&val), sizeof(T));
			}

			template<typename T>
			T readBE(std::istream& is){
				T val;
				is.read(reinterpret_cast<char*>(&val), sizeof(T));
				if (!is) throw std::runtime_error("MeshIO: unexpected end of file");
				if (hostIsLittleEndian()) val = swapBytes<T>(val);
				return val;
			}

			template<typename T>
			void writeLE(std::ofstream& ofs, T val){
				if (!hostIsLittleEndian()) val = swapBytes<T>(val);
				ofs.write(reinterpret_cast<const char*>(&val), sizeof(T));
			}

			template<typename T>
			T readLE(std::istream& is){
				T val;
				is.read(reinterpret_cast<char*>(&val), sizeof(T));
				if (!is) throw std::runtime_error("MeshIO: unexpected end of file");
				if (!hostIsLittleEndian()) val = swapBytes<T>(val);
				return val;
			}

		} // namespace binary
	} // namespace io
} // namespace residuum

#endif
