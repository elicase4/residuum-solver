#include <gtest/gtest.h>

#include "linalg/types/CSRMatrix.hpp"
#include "linalg/types/Vector.hpp"
#include "linalg/operations/MatrixOps.hpp"

using namespace residuum;

TEST(MatrixOps, MatVecSimple){

	using Backend = linalg::types::backend::CPU;
	using Mat = linalg::types::CSRMatrix<Real, Backend>;
	using Vec = linalg::types::Vector<Real, Backend>;

	Mat A(2, 2);
	A.resize(4);

    // rowPtr
    A.rowPtr()[0] = 0;
    A.rowPtr()[1] = 2;
    A.rowPtr()[2] = 4;

    // colIdx
    A.colIdx()[0] = 0; A.colIdx()[1] = 1;
    A.colIdx()[2] = 0; A.colIdx()[3] = 1;

    // data
    A.data()[0] = 2; A.data()[1] = 1;
    A.data()[2] = 1; A.data()[3] = 3;

	Vec x(2), y(2);

	x.data()[0] = 1;
	x.data()[1] = 2;

	linalg::operations::matvec(A, x, y);

	EXPECT_NEAR(y.data()[0], 4.0, 1e-12);
	EXPECT_NEAR(y.data()[1], 7.0, 1e-12);

}
