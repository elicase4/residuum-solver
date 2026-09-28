#ifndef RESIDUUM_SOLVER_TIMESTEPPER_TIMESTEPPER_HPP
#define RESIDUUM_SOLVER_TIMESTEPPER_TIMESTEPPER_HPP

#include <concepts>

#include "core/Types.hpp"

namespace residuum {
	namespace solver {

		enum class TemporalOrder { First, Second };

		template<typename T>
		concept TimeStepper = requires(T& ts) {

			{ ts.step() } -> std::same_as<bool>;

			{ ts.finished() } -> std::same_as<bool>;

			{ ts.currentStep() } -> std::same_as<Index>;

			{ ts.currentTime() } -> std::same_as<Real>;

			{ T::Order } -> std::convertible_to<TemporalOrder>;

		}; // concept TimeStepper

	} // namespace solver
} // namespace residuum

#endif
