#include <gtest/gtest.h>

#include "core/Types.hpp"
#include "equation/heateq/quantity/HeatFluxIntegrand.hpp"
#include "fem/quantity/QuantityUnits.hpp"
#include "fem/quantity/Reduction.hpp"

using namespace residuum;

namespace {

	// minimal Form stand-in -- unitFor<Form> only ever touches Form::BaseUnit
	struct MassQuantity {
		static constexpr const char* BaseUnit = "kg";
	};

} // namespace

TEST(QuantityUnitsTest, IntegralReductionReturnsBareBaseUnit){

	const std::string unit = fem::quantity::unitFor<MassQuantity>(fem::quantity::Reduction::Integral, 2);
	EXPECT_EQ(unit, "kg");

}

TEST(QuantityUnitsTest, IntegralReductionIgnoresBoundaryDim){

	const std::string unit2D = fem::quantity::unitFor<MassQuantity>(fem::quantity::Reduction::Integral, 1);
	const std::string unit3D = fem::quantity::unitFor<MassQuantity>(fem::quantity::Reduction::Integral, 2);

	EXPECT_EQ(unit2D, "kg");
	EXPECT_EQ(unit3D, unit2D);

}

TEST(QuantityUnitsTest, AverageReductionAppendsPerBoundaryDim){

	const std::string unit1D = fem::quantity::unitFor<MassQuantity>(fem::quantity::Reduction::Average, 1);
	const std::string unit2D = fem::quantity::unitFor<MassQuantity>(fem::quantity::Reduction::Average, 2);

	EXPECT_EQ(unit1D, "kg/m^1");
	EXPECT_EQ(unit2D, "kg/m^2");

}

// real Form, not just the synthetic stand-in above -- HeatFluxIntegrand::BaseUnit
// is what actually reaches the monitor CSV/console for a heat-flux monitor
TEST(QuantityUnitsTest, HeatFluxIntegrandUnitsMatchExpectedSIForm){

	using Form = equation::heateq::quantity::HeatFluxIntegrand<int>; // QuadraturePointBoundary unused by BaseUnit/unitFor

	EXPECT_EQ(fem::quantity::unitFor<Form>(fem::quantity::Reduction::Integral, 2), "W");
	EXPECT_EQ(fem::quantity::unitFor<Form>(fem::quantity::Reduction::Average, 2), "W/m^2");

}
