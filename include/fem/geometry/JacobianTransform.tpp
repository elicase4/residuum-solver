namespace residuum::fem::geometry {

template<Int SpatialDimension, Int ParametricDimension>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDimension, ParametricDimension>::mapToPhysical(const Real* nodeCoords, const Real* N, Real* x, const Index nodesPerElement){

	// initialize
	for (Index i = 0; i < SpatialDimension; ++i){
		x[i] = 0.0;
	}

	// compute x entries
	for (Index a = 0; a < nodesPerElement; ++a){
		for (Index i = 0; i < SpatialDimension; ++i){
			x[i] += N[a]  * nodeCoords[a*SpatialDimension + i];
		}
	}

}

template<Int SpatialDimension, Int ParametricDimension>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDimension, ParametricDimension>::computeJacobian(const Real* nodeCoords, const Real* dNdxi, Real* J, const Index nodesPerElement){

	// initialize
	for (Index i = 0; i < SpatialDimension*ParametricDimension; ++i){
		J[i] = 0.0;
	}

	// compute J entries
	for (Index i = 0; i < SpatialDimension; ++i){
		for (Index alpha = 0; alpha < ParametricDimension; ++alpha){
			for (Index a = 0; a < nodesPerElement; ++a){
				J[i*ParametricDimension + alpha] += dNdxi[a*ParametricDimension + alpha]  * nodeCoords[a*SpatialDimension + i];
			}
		}
	}

}

template<Int SpatialDimension, Int ParametricDimension>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDimension, ParametricDimension>::computeMetric(const Real* J, Real* g){

	for (Index beta = 0; beta < ParametricDimension; ++beta){
		for (Index alpha = 0; alpha < ParametricDimension; ++alpha){
			
			Real sum = 0.0;
			for (Index i = 0; i < SpatialDimension; ++i){
				sum += J[i*ParametricDimension + alpha] * J[i*ParametricDimension + beta];
			}

			g[beta*ParametricDimension + alpha] = sum;
			g[alpha*ParametricDimension + beta] = sum;
		}
	}

}

template<Int SpatialDimension, Int ParametricDimension>
PDE_HOST PDE_DEVICE Real JacobianTransform<SpatialDimension, ParametricDimension>::computeMeasure(const Real* g){

	Real detg = computeMatrixDeterminant(g);
	Real measure = sqrt(detg);

	return measure;

}

template<Int SpatialDimension, Int ParametricDimension>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDimension, ParametricDimension>::transformGradient(const Real* J, const Real* g, const Real* dNdxi, Real* dNdx, const Index nodesPerElement){

	// initialize
	for (Index i = 0; i < nodesPerElement*SpatialDimension; ++i){
		dNdx[i] = 0.0;
	}

	if constexpr (ParametricDimension == SpatialDimension){

		Real invJ[ParametricDimension*ParametricDimension];
		Real detJ = computeMatrixDeterminant(J);
		invertMatrix(J, detJ, invJ);

		for (Index a = 0; a < nodesPerElement; ++a){
			for (Index i = 0; i < SpatialDimension; ++i){
				for (Index alpha = 0; alpha < ParametricDimension; ++alpha){
					dNdx[a*SpatialDimension + i] += (invJ[i*ParametricDimension + alpha] * dNdxi[a*ParametricDimension + alpha]);
				}
			}
		}

	} else if constexpr (ParametricDimension != SpatialDimension){

		Real invg[ParametricDimension*ParametricDimension];
		Real detg = computeMatrixDeterminant(g);
		invertMatrix(g, detg, invg);

		for (Index a = 0; a < nodesPerElement; ++a){
			for (Index i = 0; i < SpatialDimension; ++i){
				for (Index beta = 0; beta < ParametricDimension; ++beta){
					for (Index alpha = 0; alpha < ParametricDimension; ++alpha){
						dNdx[a*SpatialDimension + i] += (J[i*ParametricDimension + beta] * invg[beta*ParametricDimension + alpha] * dNdxi[a*ParametricDimension + alpha]);
					}
				}
			}
		}

	}
	
}

template<Int SpatialDimension, Int ParametricDimension>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDimension, ParametricDimension>::computeBoundaryNormal(const Real* J, const Real* nRef, Real* n){
	
	if constexpr (ParametricDimension == SpatialDimension) {

		Real cofJ[ParametricDimension*ParametricDimension];
		computeMatrixCofactor(J, cofJ);

		for (Index i = 0; i < ParametricDimension; ++i){
			n[i] = 0.0;
			for (Index j = 0; j < ParametricDimension; ++j){
				n[i] += cofJ[i*ParametricDimension + j] * nRef[j];
			}
		}

	} else if constexpr (ParametricDimension == 2 && SpatialDimension == 3) {

		// extract a1
		Real a10 = J[0];
		Real a11 = J[2];
		Real a12 = J[4];

		// extract a2
		Real a20 = J[1];
		Real a21 = J[3];
		Real a22 = J[5];

		// compute surface normal
		Real ns0 = a11*a22 - a12*a21;
		Real ns1 = a12*a20 - a10*a22;
		Real ns2 = a10*a21 - a11*a20;

		// compute nref_purp
		Real r0 = -nRef[1];
		Real r1 =  nRef[0];

		// compute t = J*r
		Real t0 = a10*r0 + a20*r1;
		Real t1 = a11*r0 + a21*r1;
		Real t2 = a12*r0 + a22*r1;

		// compute n = t x n_s
		n[0] = t1*ns2 - t2*ns1;
		n[1] = t2*ns0 - t0*ns2;
		n[2] = t0*ns1 - t1*ns0;

	}

}

template<Int SpatialDimension, Int ParametricDimension>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDimension, ParametricDimension>::invertMatrix(const Real* A, const Real detA, Real* invA){
	
	Real detInvA = 1.0 / detA;

	if constexpr (ParametricDimension == 2) {
		
		invA[0] = detInvA  * A[3];
		invA[1] = -detInvA * A[1];
		invA[2] = -detInvA * A[2];
		invA[3] = detInvA  * A[0];

	} else if constexpr (ParametricDimension == 3){

		invA[0] = detInvA  * ( (A[4] * A[8]) - (A[5] * A[7]) );
		invA[1] = -detInvA * ( (A[3] * A[8]) - (A[5] * A[6]) );
		invA[2] = detInvA  * ( (A[3] * A[7]) - (A[4] * A[6]) );
		invA[3] = -detInvA * ( (A[1] * A[8]) - (A[2] * A[7]) );
		invA[4] = detInvA  * ( (A[0] * A[8]) - (A[2] * A[6]) );
		invA[5] = -detInvA * ( (A[0] * A[7]) - (A[1] * A[6]) );
		invA[6] = detInvA  * ( (A[1] * A[5]) - (A[2] * A[4]) );
		invA[7] = -detInvA * ( (A[0] * A[5]) - (A[2] * A[3]) );
		invA[8] = detInvA  * ( (A[0] * A[4]) - (A[1] * A[3]) );

	}

}

template<Int SpatialDimension, Int ParametricDimension>
PDE_HOST PDE_DEVICE Real JacobianTransform<SpatialDimension, ParametricDimension>::computeMatrixDeterminant(const Real* A){

	Real detA;

	if constexpr (ParametricDimension == 2) {
		
		detA = (A[0] * A[3]) - (A[1] * A[2]);

	} else if constexpr (ParametricDimension == 3){
		
		detA = (A[0] * ( (A[4] * A[8]) - (A[5] * A[7]) ))
			  -(A[1] * ( (A[3] * A[8]) - (A[5] * A[6]) ))
			  +(A[2] * ( (A[3] * A[7]) - (A[4] * A[6]) ));

	}

	return detA;

}

template<Int SpatialDimension, Int ParametricDimension>
PDE_HOST PDE_DEVICE void JacobianTransform<SpatialDimension, ParametricDimension>::computeMatrixCofactor(const Real* A, Real* cofA){

	if constexpr (ParametricDimension == 2) {
		
		cofA[0] =  A[3];
		cofA[1] = -A[2];
		cofA[2] = -A[1];
		cofA[3] =  A[0];

	} else if constexpr (ParametricDimension == 3){

		cofA[0] =  (A[4]*A[8] - A[5]*A[7]);
		cofA[1] = -(A[3]*A[8] - A[5]*A[6]);
		cofA[2] =  (A[3]*A[7] - A[4]*A[6]);

		cofA[3] = -(A[1]*A[8] - A[2]*A[7]);
		cofA[4] =  (A[0]*A[8] - A[2]*A[6]);
		cofA[5] = -(A[0]*A[7] - A[1]*A[6]);

		cofA[6] =  (A[1]*A[5] - A[2]*A[4]);
		cofA[7] = -(A[0]*A[5] - A[2]*A[3]);
		cofA[8] =  (A[0]*A[4] - A[1]*A[3]);

	}

}

} // namespace residuum::fem::geometry
