#include <gtest/gtest.h>

#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "equation/heateq/evaluator/SourceFunction.hpp"
#include "equation/heateq/boundary/BoundaryValueFunction.hpp"
#include "equation/heateq/boundary/BoundaryFluxFunction.hpp"
#include "utils/expression/ScalarExpression.hpp"
#include "utils/expression/VectorExpression.hpp"

using namespace residuum;

using Callable = utils::expression::ScalarExpression;
using VectorCallable = utils::expression::VectorExpression;

using SourceFunctionT = equation::heateq::evaluator::SourceFunction<2, 1, Callable>;
using BoundaryValueFunctionT = equation::heateq::BoundaryValueFunction<2, 1, Callable>;
using BoundaryFluxFunctionT = equation::heateq::BoundaryFluxFunction<2, 1, VectorCallable>;

// These wrappers only matter to test with a Callable that is genuinely neither
// movable nor copyable (ScalarExpression/VectorExpression, via exprtk::parser) --
// that's the exact case the forwarding constructors exist for. Confirm we're
// actually exercising that case, not a movable stand-in like a lambda.
static_assert(!std::is_move_constructible_v<SourceFunctionT>);
static_assert(!std::is_copy_constructible_v<SourceFunctionT>);
static_assert(!std::is_move_constructible_v<BoundaryValueFunctionT>);
static_assert(!std::is_copy_constructible_v<BoundaryValueFunctionT>);
static_assert(!std::is_move_constructible_v<BoundaryFluxFunctionT>);
static_assert(!std::is_copy_constructible_v<BoundaryFluxFunctionT>);

TEST(SourceFunction, ConstructsInPlaceFromRawExpressionAndEvaluates) {

	SourceFunctionT f{std::string("x*y + t")};

	Real x[3] = {2.0, 3.0, 0.0};
	Real out[1] = {0.0};
	f.eval(1.0, x, out);

	EXPECT_DOUBLE_EQ(out[0], 2.0 * 3.0 + 1.0);

}

TEST(BoundaryValueFunction, ConstructsInPlaceFromRawExpressionAndEvaluates) {

	BoundaryValueFunctionT f{std::string("(1-x)*(1-y)")};

	Real x[3] = {0.25, 0.5, 0.0};
	Real out[1] = {0.0};
	f.eval(0.0, x, out);

	EXPECT_DOUBLE_EQ(out[0], (1.0 - 0.25) * (1.0 - 0.5));

}

TEST(BoundaryFluxFunction, ConstructsInPlaceFromRawExpressionAndEvaluates) {

	BoundaryFluxFunctionT f{std::vector<std::string>{"2*x", "3*y"}};

	Real x[3] = {4.0, 5.0, 0.0};
	Real out[2] = {0.0, 0.0};
	f.eval(0.0, x, out);

	EXPECT_DOUBLE_EQ(out[0], 8.0);
	EXPECT_DOUBLE_EQ(out[1], 15.0);

}

TEST(SourceFunction, ConstructorPropagatesExpressionCompileErrors) {

	EXPECT_THROW((SourceFunctionT{std::string("not an expression (")}), std::runtime_error);

}
