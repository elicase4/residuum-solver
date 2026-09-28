#include <gtest/gtest.h>
#include <cmath>

#include "fem/dispatch/DiscretizationLimits.hpp"
#include "fem/quadrature/GaussQuadrature1D.hpp"
#include "fem/quadrature/GaussQuadratureQuad.hpp"
#include "fem/quadrature/GaussQuadratureHex.hpp"

using namespace residuum::fem;
using namespace residuum::fem::dispatch;
using namespace residuum::fem::quadrature;

namespace {

	// analytic integral of x^p over [-1, 1]
	Real analyticMonomialIntegral(int p) {
		return (p % 2 == 0) ? (2.0 / (p + 1)) : 0.0;
	}

} // namespace

// ======================================================
// GaussQuadrature1D Tests
// ======================================================

TEST(GaussQuadrature1D, WeightsSumToMeasure) {
	for (Index n = 1; n <= 5; ++n){
		GaussQuadrature1D quad(n);
		Real w[kMaxQuadraturePoints1D];
		quad.getWeights(w);
		Real sum = 0.0;
		for (Index i = 0; i < n; ++i) sum += w[i];
		EXPECT_NEAR(sum, 2.0, 1e-13) << "numPoints = " << n;
	}
}

TEST(GaussQuadrature1D, ExactForPolynomialsUpToDegree) {
	for (Index n = 1; n <= 5; ++n){
		GaussQuadrature1D quad(n);
		Real xi[kMaxQuadraturePoints1D];
		Real w[kMaxQuadraturePoints1D];
		quad.getPoints(xi);
		quad.getWeights(w);

		// an n-point Gauss rule is exact for polynomials up to degree 2n-1
		const int exactDegree = 2 * static_cast<int>(n) - 1;
		for (int p = 0; p <= exactDegree; ++p){
			Real sum = 0.0;
			for (Index i = 0; i < n; ++i) sum += w[i] * std::pow(xi[i], p);
			EXPECT_NEAR(sum, analyticMonomialIntegral(p), 1e-12) << "numPoints = " << n << ", degree = " << p;
		}
	}
}

TEST(GaussQuadrature1D, Accessors) {
	for (Index n = 1; n <= 5; ++n){
		GaussQuadrature1D quad(n);
		EXPECT_EQ(quad.numPoints(), n);
		EXPECT_EQ(quad.numPointsTotal(), n);
	}
}

// ======================================================
// GaussQuadratureQuad Tests
// ======================================================

TEST(GaussQuadratureQuad, WeightsSumToMeasure) {
	for (Index nx = 1; nx <= 4; ++nx){
		for (Index ny = 1; ny <= 4; ++ny){
			GaussQuadratureQuad quad(nx, ny);
			Real w[kMaxQuadraturePointsTotal<2>];
			quad.getWeights(w);
			Real sum = 0.0;
			for (Index q = 0; q < quad.numPointsTotal(); ++q) sum += w[q];
			EXPECT_NEAR(sum, 4.0, 1e-13) << "nx = " << nx << ", ny = " << ny;
		}
	}
}

TEST(GaussQuadratureQuad, ExactForPolynomialsUpToDegree) {
	for (Index nx = 1; nx <= 4; ++nx){
		for (Index ny = 1; ny <= 4; ++ny){
			GaussQuadratureQuad quad(nx, ny);
			Real xi[2*kMaxQuadraturePointsTotal<2>];
			Real w[kMaxQuadraturePointsTotal<2>];
			quad.getPoints(xi);
			quad.getWeights(w);

			// a tensor-product (nx,ny)-point rule is exact for x^p*y^q whenever
			// p <= 2nx-1 and q <= 2ny-1 independently, since the integral factors.
			const int exactDegreeX = 2 * static_cast<int>(nx) - 1;
			const int exactDegreeY = 2 * static_cast<int>(ny) - 1;

			for (int p = 0; p <= exactDegreeX; ++p){
				for (int q = 0; q <= exactDegreeY; ++q){
					Real sum = 0.0;
					for (Index qp = 0; qp < quad.numPointsTotal(); ++qp){
						sum += w[qp] * std::pow(xi[2*qp], p) * std::pow(xi[2*qp+1], q);
					}
					Real expected = analyticMonomialIntegral(p) * analyticMonomialIntegral(q);
					EXPECT_NEAR(sum, expected, 1e-11) << "nx=" << nx << " ny=" << ny << " p=" << p << " q=" << q;
				}
			}
		}
	}
}

TEST(GaussQuadratureQuad, Accessors) {
	for (Index nx = 1; nx <= 4; ++nx){
		for (Index ny = 1; ny <= 4; ++ny){
			GaussQuadratureQuad quad(nx, ny);
			EXPECT_EQ(quad.numPointsXi(), nx);
			EXPECT_EQ(quad.numPointsEta(), ny);
			EXPECT_EQ(quad.numPointsTotal(), nx * ny);
		}
	}
}

// ======================================================
// GaussQuadratureHex Tests
// ======================================================

TEST(GaussQuadratureHex, WeightsSumToMeasure) {
	for (Index nx = 1; nx <= 3; ++nx){
		for (Index ny = 1; ny <= 3; ++ny){
			for (Index nz = 1; nz <= 3; ++nz){
				GaussQuadratureHex quad(nx, ny, nz);
				Real w[kMaxQuadraturePointsTotal<3>];
				quad.getWeights(w);
				Real sum = 0.0;
				for (Index q = 0; q < quad.numPointsTotal(); ++q) sum += w[q];
				EXPECT_NEAR(sum, 8.0, 1e-13) << "nx = " << nx << ", ny = " << ny << ", nz = " << nz;
			}
		}
	}
}

TEST(GaussQuadratureHex, ExactForPolynomialsUpToDegree) {
	for (Index nx = 1; nx <= 3; ++nx){
		for (Index ny = 1; ny <= 3; ++ny){
			for (Index nz = 1; nz <= 3; ++nz){
				GaussQuadratureHex quad(nx, ny, nz);
				Real xi[3*kMaxQuadraturePointsTotal<3>];
				Real w[kMaxQuadraturePointsTotal<3>];
				quad.getPoints(xi);
				quad.getWeights(w);

				const int exactDegreeX = 2 * static_cast<int>(nx) - 1;
				const int exactDegreeY = 2 * static_cast<int>(ny) - 1;
				const int exactDegreeZ = 2 * static_cast<int>(nz) - 1;

				for (int p = 0; p <= exactDegreeX; ++p){
					for (int q = 0; q <= exactDegreeY; ++q){
						for (int r = 0; r <= exactDegreeZ; ++r){
							Real sum = 0.0;
							for (Index qp = 0; qp < quad.numPointsTotal(); ++qp){
								sum += w[qp] * std::pow(xi[3*qp], p) * std::pow(xi[3*qp+1], q) * std::pow(xi[3*qp+2], r);
							}
							Real expected = analyticMonomialIntegral(p) * analyticMonomialIntegral(q) * analyticMonomialIntegral(r);
							EXPECT_NEAR(sum, expected, 1e-10) << "nx=" << nx << " ny=" << ny << " nz=" << nz << " p=" << p << " q=" << q << " r=" << r;
						}
					}
				}
			}
		}
	}
}

TEST(GaussQuadratureHex, Accessors) {
	for (Index nx = 1; nx <= 3; ++nx){
		for (Index ny = 1; ny <= 3; ++ny){
			for (Index nz = 1; nz <= 3; ++nz){
				GaussQuadratureHex quad(nx, ny, nz);
				EXPECT_EQ(quad.numPointsXi(), nx);
				EXPECT_EQ(quad.numPointsEta(), ny);
				EXPECT_EQ(quad.numPointsZeta(), nz);
				EXPECT_EQ(quad.numPointsTotal(), nx * ny * nz);
			}
		}
	}
}
