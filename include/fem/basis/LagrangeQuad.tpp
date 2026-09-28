namespace residuum::fem::basis {

// Implementation: eval
PDE_HOST PDE_DEVICE PDE_INLINE void LagrangeQuad::eval(const Real* xi, Real* N) const {

	Real Nx[fem::dispatch::kMaxBasisOrder + 1];
	Real Ny[fem::dispatch::kMaxBasisOrder + 1];

	basisX_.eval(xi[0], Nx);
	basisY_.eval(xi[1], Ny);

	// compute tensor product
	Index a = 0;
	for (Index j = 0; j <= basisY_.order(); ++j){
		for (Index i = 0; i <= basisX_.order(); ++i){
			N[a] = Nx[i] * Ny[j];
			a++;
		}
	}
}

// Implementation: evalGradient
PDE_HOST PDE_DEVICE PDE_INLINE void LagrangeQuad::evalGradient(const Real* xi, Real* dNdxi) const {

	Index pD = 2;

	Real Nx[fem::dispatch::kMaxBasisOrder + 1], Ny[fem::dispatch::kMaxBasisOrder + 1];
	Real dNx[fem::dispatch::kMaxBasisOrder + 1], dNy[fem::dispatch::kMaxBasisOrder + 1];

	basisX_.eval(xi[0], Nx);
	basisX_.evalFirstDerivative(xi[0], dNx);
	basisY_.eval(xi[1], Ny);
	basisY_.evalFirstDerivative(xi[1], dNy);

	// evaluate tensor product & chain rule
	Index a = 0;
	for (Index j = 0; j <= basisY_.order(); ++j){
		for (Index i = 0; i <= basisX_.order(); ++i){
			dNdxi[a*pD    ] = dNx[i] * Ny[j];
			dNdxi[a*pD + 1] = Nx[i] * dNy[j];
			a++;
		}
	}
}

// Implementation: evalHessian
PDE_HOST PDE_DEVICE PDE_INLINE void LagrangeQuad::evalHessian(const Real* xi, Real* d2Nd2xi) const {

	Real Nx[fem::dispatch::kMaxBasisOrder + 1], Ny[fem::dispatch::kMaxBasisOrder + 1];
	Real dNx[fem::dispatch::kMaxBasisOrder + 1], dNy[fem::dispatch::kMaxBasisOrder + 1];
	Real d2Nx[fem::dispatch::kMaxBasisOrder + 1], d2Ny[fem::dispatch::kMaxBasisOrder + 1];

	basisX_.eval(xi[0], Nx);
	basisX_.evalFirstDerivative(xi[0], dNx);
	basisX_.evalSecondDerivative(xi[0], d2Nx);
	basisY_.eval(xi[1], Ny);
	basisY_.evalFirstDerivative(xi[1], dNy);
	basisY_.evalSecondDerivative(xi[1], d2Ny);

	const Index NumEntries = 3;

	// evaluate tensor product & chain rule
	Index a = 0;
	for (Index j = 0; j <= basisY_.order(); ++j){
		for (Index i = 0; i <= basisX_.order(); ++ i){
			d2Nd2xi[a*NumEntries    ] = d2Nx[i] * Ny[j];
			d2Nd2xi[a*NumEntries + 1] = Nx[i] * d2Ny[j];
			d2Nd2xi[a*NumEntries + 2] = dNx[i] * dNy[j];
			a++;
		}
	}
}

// Implementation: evalLaplacian
PDE_HOST PDE_DEVICE PDE_INLINE void LagrangeQuad::evalLaplacian(const Real* xi, Real* lapN) const {

	Real Nx[fem::dispatch::kMaxBasisOrder + 1], Ny[fem::dispatch::kMaxBasisOrder + 1];
	Real d2Nx[fem::dispatch::kMaxBasisOrder + 1], d2Ny[fem::dispatch::kMaxBasisOrder + 1];

	basisX_.eval(xi[0], Nx);
	basisX_.evalSecondDerivative(xi[0], d2Nx);
	basisY_.eval(xi[1], Ny);
	basisY_.evalSecondDerivative(xi[1], d2Ny);

	// evaluate tensor product & chain rule
	Index a = 0;
	for (Index j = 0; j <= basisY_.order(); ++j){
		for (Index i = 0; i <= basisX_.order(); ++ i){
			lapN[a] = d2Nx[i] * Ny[j] + Nx[i] * d2Ny[j];
			a++;
		}
	}
}

// Implementation: getFaceTopology
PDE_HOST PDE_DEVICE PDE_INLINE void LagrangeQuad::getFaceTopology(const Int rngID, Real* nRef) const {
	switch (rngID){
		case 0:
			nRef[0] = -1.0;
			nRef[1] = 0.0;
			break;
		case 1:
			nRef[0] = 1.0;
			nRef[1] = 0.0;
			break;
		case 2:
			nRef[0] = 0.0;
			nRef[1] = -1.0;
			break;
		case 3:
			nRef[0] = 0.0;
			nRef[1] = 1.0;
			break;
		default:
			nRef[0] = 0.0;
			nRef[1] = 0.0;
	}
}

// Implementation: nodesPerFace
PDE_HOST PDE_DEVICE PDE_INLINE Index LagrangeQuad::nodesPerFace(const Int rngID) const {
	switch (rngID) {
		case 0:
			return (basisY_.order() + 1);
		case 1:
			return (basisY_.order() + 1);
		case 2:
			return (basisX_.order() + 1);
		case 3:
			return (basisX_.order() + 1);
		default:
			return 0;
	}
}

// Implementation: getFaceNodes
PDE_HOST PDE_DEVICE PDE_INLINE void LagrangeQuad::getFaceNodes(const Int rngID, Index* nodeIDs) const {
	switch (rngID) {
		case 0:
			for (Index i = 0; i < (basisY_.order() + 1); ++i)
				nodeIDs[i] = (Index) i*(basisX_.order() + 1);
			break;
		case 1:
			for (Index i = 0; i < (basisY_.order() + 1); ++i)
				nodeIDs[i] = (Index) i*(basisX_.order() + 1) + basisX_.order();
			break;
		case 2:
			for (Index i = 0; i < (basisX_.order() + 1); ++i)
				nodeIDs[i] = (Index) i;
			break;
		case 3:
			for (Index i = 0; i < (basisX_.order() + 1); ++i)
				nodeIDs[i] = (Index) i + (basisX_.order() + 1)*basisY_.order();
			break;
	}
}

// Implementation: mapFaceToElement
PDE_HOST PDE_DEVICE PDE_INLINE void LagrangeQuad::mapFaceToElement(const Int rngID, const Real* xi_face, Real* xi_elem) const {

	const Real s = xi_face[0];

	switch (rngID) {
		case 0:
			xi_elem[0] = -1.0;
			xi_elem[1] = s;
			break;
		case 1:
			xi_elem[0] = 1.0;
			xi_elem[1] = s;
			break;
		case 2:
			xi_elem[0] = s;
			xi_elem[1] = -1.0;
			break;
		case 3:
			xi_elem[0] = s;
			xi_elem[1] = 1.0;
			break;
	}

}

} // namespace residuum::fem::basis
