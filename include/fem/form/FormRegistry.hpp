#ifndef RESIDUUM_FEM_FORM_FORMREGISTRY_HPP
#define RESIDUUM_FEM_FORM_FORMREGISTRY_HPP

#include <tuple>
#include <type_traits>
#include <utility>

#include "core/Types.hpp"
#include "config/Platform.hpp"
#include "fem/form/BilinearForm.hpp"
#include "fem/form/GathersElementData.hpp"
#include "fem/form/LinearForm.hpp"
#include "fem/form/NonlinearForm.hpp"
#include "fem/form/NonlinearTangentForm.hpp"

namespace residuum {
	namespace fem {
		namespace form {

			template<typename... Forms>
			class FormRegistry {
			public:

				constexpr FormRegistry() = default;

				constexpr explicit FormRegistry(Forms... forms) : forms_(std::forward<Forms>(forms)...) {}

				template<typename... Args>
				requires (sizeof...(Forms) == 1 && !(sizeof...(Args) == 1 && (std::is_same_v<std::remove_cvref_t<Args>, FormRegistry> && ...)))
				constexpr explicit FormRegistry(Args&&... args) : forms_(std::forward<Args>(args)...) {}

				PDE_HOST PDE_DEVICE void gatherElementData(const Index* nodeIDs, Index nodesPerElement) const {
					std::apply([&](const auto&... form) {
						([&]{
							if constexpr (GathersElementData<std::decay_t<decltype(form)>>) {
								form.gatherElementData(nodeIDs, nodesPerElement);
							}
						}(), ...);
					}, forms_);
				}

				template<typename QuadraturePointT>
				requires ((LinearForm<Forms, QuadraturePointT> || NonlinearForm<Forms, QuadraturePointT>) && ...)
				PDE_HOST PDE_DEVICE void computeElementLevelVector(const QuadraturePointT& qp, Real* Ue, Real* Fe) const {
					std::apply([&](const auto&... form) {
						(form.computeElementLevelVector(qp, Ue, Fe), ...);
					}, forms_);
				}

				template<typename QuadraturePointT>
				requires ((BilinearForm<Forms, QuadraturePointT> || NonlinearTangentForm<Forms, QuadraturePointT>) && ...)
				PDE_HOST PDE_DEVICE void computeElementLevelMatrix(const QuadraturePointT& qp, Real* Ke) const {
					std::apply([&](const auto&... form) {
						(form.computeElementLevelMatrix(qp, Ke), ...);
					}, forms_);
				}

				static constexpr Index numForms() { return sizeof...(Forms); }

				template<std::size_t I>
				constexpr const auto& get() const {
					return std::get<I>(forms_);
				}

			private:

				std::tuple<Forms...> forms_;

			}; // class FormRegistry

		} // namespace form
	} // namespace fem
} // namespace residuum

#endif
