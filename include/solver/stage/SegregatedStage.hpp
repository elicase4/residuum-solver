#ifndef RESIDUUM_SOLVER_STAGE_SEGREGATEDSTAGE_HPP
#define RESIDUUM_SOLVER_STAGE_SEGREGATEDSTAGE_HPP

#include <tuple>
#include <cmath>
#include <utility>

#include "core/Types.hpp"
#include "solver/stage/Stage.hpp"

namespace residuum {
	namespace solver {
		namespace stage {

			template<typename ResidualFunctorT, Stage... StageTypes>
			class SegregatedStage {
			public:

				struct Config {
					Index maxOuterIterations = 50;
					Real absoluteTolerance = 1e-8;
					Real relativeTolerance = 1e-6;
				}; // struct Config

				explicit SegregatedStage(Config cfg, ResidualFunctorT residual, StageTypes&... stages) : cfg_(cfg), residual_(residual), stages_(stages...) {}

				void initialize() {
					initializeAll(std::index_sequence_for<StageTypes...>{});
				}

				void assemble() {}

				bool solve() {

					Real residual0 = -1.0;

					for (Index iter = 0; iter < cfg_.maxOuterIterations; ++iter) {

						// segregated solve for all stages in specified order
						bool allConverged = sweepAll(std::index_sequence_for<StageTypes...>{});
						if (!allConverged) {
							return false;
						}

						// evaluate inter-field coupling residual
						Real res = residual_(iter);
						if (iter == 0) {
							residual0 = res;
						}

						// check absolute and relative convergence
						if (res < cfg_.absoluteTolerance) {
							return true;
						}

						if (residual0 > 0.0 && (res / residual0) < cfg_.relativeTolerance) {
							return true;
						}

					}

					return false;

				}

				void finalize() {
					finalizeAll(std::index_sequence_for<StageTypes...>{});
				}

			private:

				template<Index... Is>
				void initializeAll(std::index_sequence<Is...>) {
					(std::get<Is>(stages_).initialize(), ...);
				}

				template<Index... Is>
				void finalizeAll(std::index_sequence<Is...>) {
					(std::get<Is>(stages_).finalize(), ...);
				}

				template<Index... Is>
				bool sweepAll(std::index_sequence<Is...>) {
					return (... && solveStage(std::get<Is>(stages_)));
				}

				template<typename ST>
				bool solveStage(ST& stage) {
					stage.assemble();
					return stage.solve();
				}

				Config cfg_;
				ResidualFunctorT residual_;
				std::tuple<StageTypes&...> stages_;

			}; // class SegregatedStage

		} // namespace stage
	} // namespace solver
} // namespace residuum

#endif
