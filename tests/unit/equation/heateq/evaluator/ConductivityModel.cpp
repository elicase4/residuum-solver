#include <gtest/gtest.h>

#include <array>

#include "equation/heateq/evaluator/ConductivityModel.hpp"

using namespace residuum;

namespace {

	template<Index NSD>
	struct MockQP {
		static constexpr Index SpatialDim = NSD;
		std::array<Real, NSD * NSD> K{};
		std::array<Real, NSD * NSD> dKdT{};
		struct { Real value = Real(0); } T;
	};

	using ModelT2D = equation::heateq::evaluator::ConductivityModel<MockQP<2>>;
	using ModelT3D = equation::heateq::evaluator::ConductivityModel<MockQP<3>>;

} // namespace

TEST(ConductivityModel, SetConstantEvaluatesToScaledIdentity) {

	ModelT2D model;
	model.setConstant(3.0);

	MockQP<2> qp;
	model.eval(qp);

	EXPECT_DOUBLE_EQ(qp.K[0 * 2 + 0], 3.0);
	EXPECT_DOUBLE_EQ(qp.K[0 * 2 + 1], 0.0);
	EXPECT_DOUBLE_EQ(qp.K[1 * 2 + 0], 0.0);
	EXPECT_DOUBLE_EQ(qp.K[1 * 2 + 1], 3.0);

}

TEST(ConductivityModel, SetConstantGradientIsZero) {

	ModelT2D model;
	model.setConstant(3.0);

	MockQP<2> qp;
	model.evalGradient(qp);

	for (Real v : qp.dKdT) EXPECT_DOUBLE_EQ(v, 0.0);

}

// Regression test: setAnisotropic() used to leave `scalar` at its default-constructed
// 0, so eval()'s k*tensor[i,j] silently zeroed out the entire tensor regardless of the
// values passed in -- caught via examples/heateq/steady/anisotropic_conductivity trivially
// "converging" with zero residual because the assembled operator was the zero matrix.
TEST(ConductivityModel, SetAnisotropicEvaluatesFullTensorUnscaled) {

	ModelT2D model;
	model.setAnisotropic({{2.0, 0.5}, {0.5, 1.5}});

	MockQP<2> qp;
	model.eval(qp);

	EXPECT_DOUBLE_EQ(qp.K[0 * 2 + 0], 2.0);
	EXPECT_DOUBLE_EQ(qp.K[0 * 2 + 1], 0.5);
	EXPECT_DOUBLE_EQ(qp.K[1 * 2 + 0], 0.5);
	EXPECT_DOUBLE_EQ(qp.K[1 * 2 + 1], 1.5);

}

TEST(ConductivityModel, SetAnisotropicGradientIsZero) {

	ModelT2D model;
	model.setAnisotropic({{2.0, 0.5}, {0.5, 1.5}});

	MockQP<2> qp;
	model.evalGradient(qp);

	for (Real v : qp.dKdT) EXPECT_DOUBLE_EQ(v, 0.0);

}

TEST(ConductivityModel, SetAnisotropic3DEvaluatesFullTensorUnscaled) {

	ModelT3D model;
	model.setAnisotropic({{2.0, 0.5, 0.3}, {0.5, 1.5, 0.2}, {0.3, 0.2, 1.0}});

	MockQP<3> qp;
	model.eval(qp);

	const std::array<Real, 9> expected = {2.0, 0.5, 0.3, 0.5, 1.5, 0.2, 0.3, 0.2, 1.0};
	for (Index i = 0; i < 9; ++i) EXPECT_DOUBLE_EQ(qp.K[i], expected[i]);

}

TEST(ConductivityModel, SetTemperatureDependentIsotropicEvaluatesExpressionScaledIdentity) {

	ModelT2D model;
	model.setTemperatureDependentIsotropic("1.0 + 0.05*T", "0.05");

	MockQP<2> qp;
	qp.T.value = 100.0;

	model.eval(qp);
	EXPECT_DOUBLE_EQ(qp.K[0 * 2 + 0], 6.0); // k(100) = 1 + 0.05*100 = 6
	EXPECT_DOUBLE_EQ(qp.K[0 * 2 + 1], 0.0);
	EXPECT_DOUBLE_EQ(qp.K[1 * 2 + 1], 6.0);

	model.evalGradient(qp);
	EXPECT_DOUBLE_EQ(qp.dKdT[0 * 2 + 0], 0.05);
	EXPECT_DOUBLE_EQ(qp.dKdT[0 * 2 + 1], 0.0);
	EXPECT_DOUBLE_EQ(qp.dKdT[1 * 2 + 1], 0.05);

}

TEST(ConductivityModel, SetTemperatureDependentAnisotropicEvaluatesExpressionScaledTensor) {

	ModelT2D model;
	model.setTemperatureDependentAnisotropic({{2.0, 0.5}, {0.5, 1.5}}, "1.0 + 0.05*T", "0.05");

	MockQP<2> qp;
	qp.T.value = 100.0; // k(100) = 6, dk(100) = 0.05

	model.eval(qp);
	EXPECT_DOUBLE_EQ(qp.K[0 * 2 + 0], 12.0);
	EXPECT_DOUBLE_EQ(qp.K[0 * 2 + 1], 3.0);
	EXPECT_DOUBLE_EQ(qp.K[1 * 2 + 1], 9.0);

	model.evalGradient(qp);
	EXPECT_DOUBLE_EQ(qp.dKdT[0 * 2 + 0], 0.1);
	EXPECT_DOUBLE_EQ(qp.dKdT[0 * 2 + 1], 0.025);
	EXPECT_DOUBLE_EQ(qp.dKdT[1 * 2 + 1], 0.075);

}

TEST(ConductivityModel, SetAnisotropicWrongSizeTensorThrows) {

	ModelT2D model;
	EXPECT_THROW(model.setAnisotropic({{2.0, 0.5, 0.1}, {0.5, 1.5, 0.1}}), std::runtime_error);
	EXPECT_THROW(model.setAnisotropic({{2.0, 0.5}}), std::runtime_error);

}
