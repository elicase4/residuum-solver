#ifndef RESIDUUM_IO_VISUALIZATION_VTUSERIESWRITER_HPP
#define RESIDUUM_IO_VISUALIZATION_VTUSERIESWRITER_HPP

#include <string>
#include <utility>
#include <vector>

#include "core/Types.hpp"

namespace residuum {
	namespace io {
		namespace visualization {

			class VTUSeriesWriter {
			public:

				VTUSeriesWriter(std::string directory, std::string prefix, Index digitWidth = 6);

				std::string addStep(Index step, Real time);

			private:

				void writePVD() const;

				std::string directory_;
				std::string prefix_;
				Index digitWidth_;

				std::vector<std::pair<std::string, Real>> entries_;

			}; // class VTUSeriesWriter

		} // namespace visualization
	} // namespace io
} // namespace residuum

#endif
