#ifndef RESIDUUM_FEM_GEOMETRY_JACOBIANTRANSFORM_HPP
#define RESIDUUM_FEM_GEOMETRY_JACOBIANTRANSFORM_HPP

#include <math.h>

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace residuum {
	namespace fem {
		namespace geometry {
			
			template<Int SpatialDimension, Int ParametricDimension>
			class JacobianTransform {
			public:
				PDE_HOST PDE_DEVICE static void mapToPhysical(const Real* nodeCoords, const Real* N, Real* x, const Index nodesPerElement);
				
				PDE_HOST PDE_DEVICE static void computeJacobian(const Real* nodeCoords, const Real* dNdxi, Real* J, const Index nodesPerElement);

				PDE_HOST PDE_DEVICE static Real computeMeasure(const Real* g);

				PDE_HOST PDE_DEVICE static void computeMetric(const Real* J, Real* g);

				PDE_HOST PDE_DEVICE static void transformGradient(const Real* J, const Real* g, const Real* dNdxi, Real* dNdx, const Index nodesPerElement);

				PDE_HOST PDE_DEVICE static void computeBoundaryNormal(const Real* J, const Real* nRef, Real* n);

				static constexpr Index SpatialDim = SpatialDimension;
				static constexpr Index ParametricDim = ParametricDimension;

			private:
				PDE_HOST PDE_DEVICE static void invertMatrix(const Real* A, const Real detA, Real* invA);
				
				PDE_HOST PDE_DEVICE static Real computeMatrixDeterminant(const Real* A);
				
				PDE_HOST PDE_DEVICE static void computeMatrixCofactor(const Real* A, Real* cofA);

			}; // class JacobianTransform
			
		} // namespace geometry
	} // namespace fem
} // namespace residuum

#include "JacobianTransform.tpp"

#endif
