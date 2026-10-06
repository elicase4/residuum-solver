#ifndef RESIDUUM_FEM_QUANTITY_QUANTITYFORMS_HPP
#define RESIDUUM_FEM_QUANTITY_QUANTITYFORMS_HPP

#include <array>
#include <tuple>
#include <type_traits>
#include <utility>

#include "config/Platform.hpp"
#include "core/Types.hpp"

#include "fem/quantity/QuantityForm.hpp"
#include "fem/quantity/Reduction.hpp"

namespace residuum {
	namespace fem {
		namespace quantity {

			// pairs one QuantityForm with its Reduction mode at compile time
			template<typename FormT, Reduction ReductionMode>
			struct ReducedQuantity {
				using FormType = FormT;
				static constexpr Reduction Mode = ReductionMode;
			}; // struct ReducedQuantity

			template<typename... ReducedQuantities>
			class QuantityForms {
			public:

				static constexpr Index NumQuantities = sizeof...(ReducedQuantities);
				static constexpr Index TotalComponents = (ReducedQuantities::FormType::NumComponents + ... + 0);

				template<typename... Args>
				requires (!(sizeof...(Args) == 1 && (std::is_same_v<std::remove_cvref_t<Args>, QuantityForms> && ...)))
				constexpr explicit QuantityForms(Args&&... args) : forms_(std::forward<Args>(args)...) {}

				// each form writes its own contribution, then gets added into its own slice of out
				template<typename QuadraturePointT>
				requires (QuantityForm<typename ReducedQuantities::FormType, QuadraturePointT> && ...)
				PDE_HOST PDE_DEVICE void computeElementLevelValue(const QuadraturePointT& qp, const Real* Ue, Real* out) const {
					Index offset = 0;
					std::apply([&](const auto&... form) {
						([&] {
							using FormT = std::decay_t<decltype(form)>;
							Real contribution[FormT::NumComponents] = {0};
							form.computeElementLevelValue(qp, Ue, contribution);
							for (Index c = 0; c < FormT::NumComponents; ++c) out[offset + c] += contribution[c];
							offset += FormT::NumComponents;
						}(), ...);
					}, forms_);
				}

				// divides each Average-mode slice by measureSum; Integral-mode slices are untouched
				void finalize(Real* out, Real measureSum) const {
					Index offset = 0;
					([&] {
						if constexpr (ReducedQuantities::Mode == Reduction::Average) {
							for (Index c = 0; c < ReducedQuantities::FormType::NumComponents; ++c) out[offset + c] /= measureSum;
						}
						offset += ReducedQuantities::FormType::NumComponents;
					}(), ...);
				}

				// compile-time offset of the I-th quantity's slice
				template<std::size_t I>
				static constexpr Index offset() {
					return Offsets[I];
				}

				template<std::size_t I>
				const auto& get() const { return std::get<I>(forms_); }

			private:

				static constexpr std::array<Index, NumQuantities> Offsets = [] {
					std::array<Index, NumQuantities> offs{};
					Index running = 0;
					std::size_t i = 0;
					((offs[i++] = running, running += ReducedQuantities::FormType::NumComponents), ...);
					return offs;
				}();

				std::tuple<typename ReducedQuantities::FormType...> forms_;

			}; // class QuantityForms

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#endif
