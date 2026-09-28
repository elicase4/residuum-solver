#ifndef RESIDUUM_FEM_EVAL_MODELREGISTRY_HPP
#define RESIDUUM_FEM_EVAL_MODELREGISTRY_HPP

#include <tuple>
#include <utility>

#include "config/Platform.hpp"

namespace residuum {
	namespace fem {
		namespace evaluator {

			template<typename... Models>
			class ModelRegistry {
			public:

				constexpr ModelRegistry() = default;

				constexpr explicit ModelRegistry(Models... models) : models_(std::forward<Models>(models)...) {}

				template<typename QuadraturePointT>
				PDE_HOST PDE_DEVICE void eval(QuadraturePointT& qp) const {
					std::apply([&](const auto&... model) {
						(model.eval(qp), ...);
					}, models_);
				}

				template<typename QuadraturePointT>
				PDE_HOST PDE_DEVICE void evalGradient(QuadraturePointT& qp) const {
					std::apply([&](const auto&... model) {
						(model.evalGradient(qp), ...);
					}, models_);
				}

			private:

				std::tuple<Models...> models_;

			}; // class ModelRegistry

		} // namespace evaluator
	} // namespace fem
} // namespace residuum

#endif
