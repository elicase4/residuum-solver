#ifndef RESIDUUM_APPLICATION_HEATEQ_HEATAPPLICATION_HPP
#define RESIDUUM_APPLICATION_HEATEQ_HEATAPPLICATION_HPP

#include "application/heateq/config/HeatConfig.hpp"

namespace residuum {
	namespace application {
		namespace heateq {

			class HeatApplication {
			public:

				explicit HeatApplication(const config::HeatConfig& config);

				int run();

			private:

				config::HeatConfig config_;

			}; // class HeatApplication

			// entry point called by the residuum subcommand dispatcher; argv[1] is the config path
			int runApp(int argc, char** argv);

		} // namespace heateq
	} // namespace application
} // namespace residuum

#endif
