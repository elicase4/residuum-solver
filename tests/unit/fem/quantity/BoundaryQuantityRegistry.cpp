#include <gtest/gtest.h>

#include "core/Types.hpp"
#include "fem/quantity/BoundaryQuantityRegistry.hpp"

using namespace residuum;

// Verifies the dense tagToIndex_ fast path agrees with the old map-based lookup, using sparse non-contiguous tags (RNG tags are rarely 0..N-1).
class BoundaryQuantityRegistryTest : public ::testing::Test {
protected:

	static constexpr Index NumComponents = 2;
	fem::quantity::BoundaryQuantityRegistry<NumComponents> registry;

	void SetUp() override {
		registry.registerTag(0);
		registry.registerTag(2);
		registry.registerTag(5);
	}

}; // class BoundaryQuantityRegistryTest

TEST_F(BoundaryQuantityRegistryTest, HasTagMatchesRegisteredSet) {

	EXPECT_TRUE(registry.hasTag(0));
	EXPECT_TRUE(registry.hasTag(2));
	EXPECT_TRUE(registry.hasTag(5));

	EXPECT_FALSE(registry.hasTag(1));
	EXPECT_FALSE(registry.hasTag(3));
	EXPECT_FALSE(registry.hasTag(4));

}

TEST_F(BoundaryQuantityRegistryTest, TouchesAnyFindsARegisteredTagAmongSeveralFaceTags) {

	Int facesAllUnregistered[] = {1, 3, 4};
	EXPECT_FALSE(registry.touchesAny(facesAllUnregistered, 3));

	Int facesOneRegistered[] = {1, 3, 5};
	EXPECT_TRUE(registry.touchesAny(facesOneRegistered, 3));

	Int facesAllRegistered[] = {0, 2, 5};
	EXPECT_TRUE(registry.touchesAny(facesAllRegistered, 3));

}

TEST_F(BoundaryQuantityRegistryTest, EntryForTagReadsAndWritesThroughToTheSameStorageAsEntries) {

	auto& viaFastPath = registry.entryForTag(5);
	viaFastPath.result[0] = 3.0;
	viaFastPath.result[1] = 4.0;
	viaFastPath.measureSum = 7.0;

	// entries() is still the owning store -- a write through entryForTag must be visible there too
	const auto& viaMap = registry.entries().at(5);
	EXPECT_EQ(viaMap.result[0], 3.0);
	EXPECT_EQ(viaMap.result[1], 4.0);
	EXPECT_EQ(viaMap.measureSum, 7.0);

	EXPECT_EQ(registry.result(5)[0], 3.0);
	EXPECT_EQ(registry.result(5)[1], 4.0);

}

TEST_F(BoundaryQuantityRegistryTest, ResetZeroesEveryRegisteredTag) {

	registry.entryForTag(0).result[0] = 1.0;
	registry.entryForTag(2).measureSum = 2.0;
	registry.entryForTag(5).result[1] = 3.0;

	registry.reset();

	for (Int tag : {0, 2, 5}) {
		EXPECT_EQ(registry.result(tag)[0], 0.0);
		EXPECT_EQ(registry.result(tag)[1], 0.0);
		EXPECT_EQ(registry.entries().at(tag).measureSum, 0.0);
	}

}
