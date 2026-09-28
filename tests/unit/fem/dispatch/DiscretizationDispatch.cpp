#include <gtest/gtest.h>

#include "core/Types.hpp"
#include "fem/dispatch/DiscretizationDispatch.hpp"
#include "mesh/ElementFamily.hpp"

using namespace residuum;

namespace {

	// never actually reached when validation correctly throws first
	struct NoopVisitor {
		template<Index NSD, Index NPD, mesh::ElementFamily Family>
		bool operator()(auto&&, auto&&, auto&&) const { return true; }
	};

}

TEST(DiscretizationDispatch, InRangeQuadSucceeds) {

	NoopVisitor visitor;
	EXPECT_TRUE(fem::dispatch::dispatchQuad<2>(1, 1, 2, 2, visitor));

}

TEST(DiscretizationDispatch, BasisOrderTooHighThrows) {

	NoopVisitor visitor;
	EXPECT_THROW(fem::dispatch::dispatchQuad<2>(fem::dispatch::kMaxBasisOrder + 1, 1, 2, 2, visitor), std::runtime_error);

}

TEST(DiscretizationDispatch, QuadraturePointsTooHighThrows) {

	NoopVisitor visitor;
	EXPECT_THROW(fem::dispatch::dispatchQuad<2>(1, 1, fem::dispatch::kMaxQuadraturePoints1D + 1, 2, visitor), std::runtime_error);

}

TEST(DiscretizationDispatch, HexBasisOrderTooHighThrows) {

	NoopVisitor visitor;
	EXPECT_THROW(fem::dispatch::dispatchHex<3>(1, 1, fem::dispatch::kMaxBasisOrder + 1, 2, 2, 2, visitor), std::runtime_error);

}
