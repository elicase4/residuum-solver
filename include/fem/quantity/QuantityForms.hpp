#ifndef RESIDUUM_FEM_QUANTITY_QUANTITYFORMS_HPP
#define RESIDUUM_FEM_QUANTITY_QUANTITYFORMS_HPP

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

				template<typename QuadraturePointT>
				requires (QuantityForm<typename ReducedQuantities::FormType, QuadraturePointT> && ...)
				PDE_HOST PDE_DEVICE void computeElementLevelValue(const QuadraturePointT& qp, const Real* Ue, Real* out) const {
					computeElementLevelValueImpl(qp, Ue, out, std::index_sequence_for<ReducedQuantities...>{});
				}

				// finalizes Average-mode slices
				void finalize(Real* out, Real measureSum) const {
					finalizeImpl(out, measureSum, std::index_sequence_for<ReducedQuantities...>{});
				}

				// compile-time offset of the I-th quantity's slice
				template<std::size_t I>
				static constexpr Index offset() {
					return offsetImpl(std::make_index_sequence<I>{});
				}

				template<std::size_t I>
				const auto& get() const { return std::get<I>(forms_); }

			private:

				template<typename QuadraturePointT, std::size_t... Is>
				PDE_HOST PDE_DEVICE void computeElementLevelValueImpl(const QuadraturePointT& qp, const Real* Ue, Real* out, std::index_sequence<Is...>) const {
					(computeOne<Is>(qp, Ue, out), ...);
				}

				template<std::size_t I, typename QuadraturePointT>
				PDE_HOST PDE_DEVICE void computeOne(const QuadraturePointT& qp, const Real* Ue, Real* out) const {
					using Q = std::tuple_element_t<I, std::tuple<ReducedQuantities...>>;
					Real contribution[Q::FormType::NumComponents] = {0};
					std::get<I>(forms_).computeElementLevelValue(qp, Ue, contribution);
					for (Index c = 0; c < Q::FormType::NumComponents; ++c) {
						out[offset<I>() + c] += contribution[c];
					}
				}

				template<std::size_t... Is>
				void finalizeImpl(Real* out, Real measureSum, std::index_sequence<Is...>) const {
					(finalizeOne<Is>(out, measureSum), ...);
				}

				template<std::size_t I>
				void finalizeOne(Real* out, Real measureSum) const {
					using Q = std::tuple_element_t<I, std::tuple<ReducedQuantities...>>;
					if constexpr (Q::Mode == Reduction::Average) {
						for (Index c = 0; c < Q::FormType::NumComponents; ++c) {
							out[offset<I>() + c] /= measureSum;
						}
					}
				}

				template<std::size_t... Is>
				static constexpr Index offsetImpl(std::index_sequence<Is...>) {
					return (std::tuple_element_t<Is, std::tuple<ReducedQuantities...>>::FormType::NumComponents + ... + 0);
				}

				std::tuple<typename ReducedQuantities::FormType...> forms_;

			}; // class QuantityForms

		} // namespace quantity
	} // namespace fem
} // namespace residuum

#endif
