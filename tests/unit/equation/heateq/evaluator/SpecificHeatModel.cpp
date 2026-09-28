#include <gtest/gtest.h>

#include "equation/heateq/evaluator/SpecificHeatModel.hpp"

using namespace residuum;

namespace {

	struct MockQP {
		Real cp = Real(0);
		Real dcpdT = Real(0);
		struct { Real value = Real(0); } T;
	};

	using ModelT = equation::heateq::evaluator::SpecificHeatModel<MockQP>;

} // namespace

TEST(SpecificHeatModel, SetConstantEvaluatesToConstantValue) {

	ModelT model;
	model.setConstant(2.0);

	MockQP qp;
	model.eval(qp);

	EXPECT_DOUBLE_EQ(qp.cp, 2.0);

}

TEST(SpecificHeatModel, SetConstantGradientIsZero) {

	ModelT model;
	model.setConstant(2.0);

	MockQP qp;
	model.evalGradient(qp);

	EXPECT_DOUBLE_EQ(qp.dcpdT, 0.0);

}

TEST(SpecificHeatModel, SetTemperatureDependentEvaluatesExpression) {

	ModelT model;
	model.setTemperatureDependent("1.0 + 0.05*T", "0.05");

	MockQP qp;
	qp.T.value = 100.0;

	model.eval(qp);
	EXPECT_DOUBLE_EQ(qp.cp, 6.0); // cp(100) = 1 + 0.05*100 = 6

	model.evalGradient(qp);
	EXPECT_DOUBLE_EQ(qp.dcpdT, 0.05);

}

TEST(SpecificHeatModel, SetTemperatureDependentTracksTemperatureChanges) {

	ModelT model;
	model.setTemperatureDependent("1.0 + 0.05*T", "0.05");

	MockQP qp;

	qp.T.value = 0.0;
	model.eval(qp);
	EXPECT_DOUBLE_EQ(qp.cp, 1.0);

	qp.T.value = 200.0;
	model.eval(qp);
	EXPECT_DOUBLE_EQ(qp.cp, 11.0);

}
