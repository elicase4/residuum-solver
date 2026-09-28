#include <gtest/gtest.h>

#include <string>
#include <type_traits>
#include <utility>

#include "core/Types.hpp"
#include "fem/form/FormRegistry.hpp"

#include "equation/heateq/HeatEquation.hpp"
#include "mesh/ElementFamily.hpp"
#include "utils/expression/ScalarExpression.hpp"

using namespace residuum;

namespace {

	// A minimal, purpose-built Form wrapping a resource that is neither movable
	// nor copyable, mirroring what SourceForm<...ScalarExpression> looks like
	// structurally but without dragging in exprtk/the FEM stack. Isolates
	// FormRegistry's own construction mechanics from everything downstream.
	struct Immovable {

		int value;
		explicit Immovable(int v) : value(v) {}
		Immovable(const Immovable&) = delete;
		Immovable& operator=(const Immovable&) = delete;
		Immovable(Immovable&&) = delete;
		Immovable& operator=(Immovable&&) = delete;

	}; // struct Immovable

	struct NonMovableForm {

		Immovable data;

		template<typename... Args>
			requires (!(sizeof...(Args) == 1 && (std::is_same_v<std::remove_cvref_t<Args>, NonMovableForm> && ...)))
		constexpr explicit NonMovableForm(Args&&... args) : data(std::forward<Args>(args)...) {}

		template<typename QuadraturePoint>
		void computeElementLevelVector(const QuadraturePoint&, const Real*, Real* Fe) const {
			Fe[0] += static_cast<Real>(data.value);
		}

		template<typename QuadraturePoint>
		void computeElementLevelMatrix(const QuadraturePoint&, Real*) const {}

	}; // struct NonMovableForm

	struct StatelessForm {

		template<typename QuadraturePoint>
		void computeElementLevelVector(const QuadraturePoint&, const Real*, Real*) const {}

		template<typename QuadraturePoint>
		void computeElementLevelMatrix(const QuadraturePoint&, Real*) const {}

	}; // struct StatelessForm

	struct AnotherStatelessForm {

		template<typename QuadraturePoint>
		void computeElementLevelVector(const QuadraturePoint&, const Real*, Real*) const {}

		template<typename QuadraturePoint>
		void computeElementLevelMatrix(const QuadraturePoint&, Real*) const {}

	}; // struct AnotherStatelessForm

	struct DummyQP {};

} // namespace

static_assert(!std::is_move_constructible_v<NonMovableForm>);
static_assert(!std::is_copy_constructible_v<NonMovableForm>);

TEST(FormRegistry, SingleFormForwardingConstructorBuildsInPlace) {

	fem::form::FormRegistry<NonMovableForm> registry(42);

	Real Fe[1] = {0.0};
	DummyQP qp;
	registry.computeElementLevelVector(qp, nullptr, Fe);

	EXPECT_DOUBLE_EQ(Fe[0], 42.0);

}

TEST(FormRegistry, DefaultConstructionStillWorksForStatelessForms) {

	fem::form::FormRegistry<StatelessForm> registry;
	EXPECT_EQ(registry.numForms(), 1u);

}

TEST(FormRegistry, CopyAndMoveStillWorkWhenFormIsCopyable) {

	// Regression check for the forwarding ctor's self-type exclusion guard: without
	// it, these would resolve to the forwarding ctor instead of the implicit
	// copy/move ctor and fail to compile.
	fem::form::FormRegistry<StatelessForm> a;
	fem::form::FormRegistry<StatelessForm> b(a);
	fem::form::FormRegistry<StatelessForm> c(std::move(a));

	EXPECT_EQ(b.numForms(), 1u);
	EXPECT_EQ(c.numForms(), 1u);

}

TEST(FormRegistry, MultiFormByValueConstructorStillWorks) {

	fem::form::FormRegistry<StatelessForm, AnotherStatelessForm> registry{StatelessForm{}, AnotherStatelessForm{}};
	EXPECT_EQ(registry.numForms(), 2u);

}

TEST(FormRegistry, SingleFormForwardingConstructorBuildsRealScalarExpressionSourceFormInPlace) {

	using Callable = utils::expression::ScalarExpression;
	using Bundle = equation::HeatEquation<2, 2, mesh::ElementFamily::Quad>;
	using SourceFormT = Bundle::SourceForm<Callable>;
	using SourceFormsT = fem::form::FormRegistry<SourceFormT>;

	static_assert(!std::is_move_constructible_v<SourceFormT>);
	static_assert(!std::is_copy_constructible_v<SourceFormT>);

	// Constructed straight from the raw expression string -- exercising the exact
	// path HeatProblem's constructor now uses for sourceForms_.
	SourceFormsT forms{std::string("x*y")};

	Real x[3] = {2.0, 5.0, 0.0};
	Real out[1] = {0.0};
	forms.get<0>().sourceFunction.eval(0.0, x, out);

	EXPECT_DOUBLE_EQ(out[0], 10.0);

}
