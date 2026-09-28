#include <gtest/gtest.h>

#include "linalg/types/Vector.hpp"
#include "linalg/operations/VectorOps.hpp"

using namespace residuum;

TEST(VectorOpsCpu, Dot){

	using Backend = linalg::types::backend::CPU;
	using Vec = linalg::types::Vector<Real, Backend>;

	Vec a(3), b(3);
	
	a.data()[0] = 1.0; a.data()[1] = 2.0; a.data()[2] = 3.0;
	b.data()[0] = 4.0; b.data()[1] = 5.0; b.data()[2] = 6.0;

	Real result = linalg::operations::dot(a, b);

	EXPECT_NEAR(result, 32.0, 1e-12);

}

TEST(VectorOpsCpu, Axpy){

	using Backend = linalg::types::backend::CPU;
	using Vec = linalg::types::Vector<Real, Backend>;

	Vec x(3), y(3);
	
	x.data()[0] = 1.0; x.data()[1] = 2.0; x.data()[2] = 3.0;
	y.data()[0] = 1.0; y.data()[1] = 1.0; y.data()[2] = 1.0;

	linalg::operations::axpy(2.0, x, y);

	EXPECT_NEAR(y.data()[0], 3.0, 1e-12);
	EXPECT_NEAR(y.data()[1], 5.0, 1e-12);
	EXPECT_NEAR(y.data()[2], 7.0, 1e-12);

}

TEST(VectorOpsCpu, Scal){

	using Backend = linalg::types::backend::CPU;
	using Vec = linalg::types::Vector<Real, Backend>;

	Vec x(3);
	
	x.data()[0] = 1.0; x.data()[1] = 2.0; x.data()[2] = 3.0;

	linalg::operations::scal(2.0, x);

	EXPECT_NEAR(x.data()[0], 2.0, 1e-12);
	EXPECT_NEAR(x.data()[1], 4.0, 1e-12);
	EXPECT_NEAR(x.data()[2], 6.0, 1e-12);

}

TEST(VectorOpsCpu, Copy){

	using Backend = linalg::types::backend::CPU;
	using Vec = linalg::types::Vector<Real, Backend>;

	Vec x(3), y(3);
	
	x.data()[0] = 7.0; x.data()[1] = 8.0; x.data()[2] = 9.0;

	linalg::operations::copy(x, y);

	EXPECT_NEAR(y.data()[0], 7.0, 1e-12);
	EXPECT_NEAR(y.data()[1], 8.0, 1e-12);
	EXPECT_NEAR(y.data()[2], 9.0, 1e-12);

}

TEST(VectorOpsCpu, Norm){

	using Backend = linalg::types::backend::CPU;
	using Vec = linalg::types::Vector<Real, Backend>;

	Vec x(3);
	
	x.data()[0] = 3.0; x.data()[1] = 4.0; x.data()[2] = 0.0;

	Real n = linalg::operations::norm(x);

	EXPECT_NEAR(n, 5.0, 1e-12);

}
