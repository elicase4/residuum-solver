namespace residuum::fem::basis {

// Implementation: eval
PDE_HOST PDE_DEVICE PDE_INLINE void LagrangeHex::eval(const Real* xi, Real* N) const {

	Real Nx[fem::dispatch::kMaxBasisOrder + 1];
	Real Ny[fem::dispatch::kMaxBasisOrder + 1];
	Real Nz[fem::dispatch::kMaxBasisOrder + 1];

	basisX_.eval(xi[0], Nx);
	basisY_.eval(xi[1], Ny);
	basisZ_.eval(xi[2], Nz);

	// compute tensor product
	Index a = 0;
	for (Index k = 0; k <= basisZ_.order(); ++k){
		for (Index j = 0; j <= basisY_.order(); ++j){
			for (Index i = 0; i <= basisX_.order(); ++i){
				N[a] = Nx[i] * Ny[j] * Nz[k];
				a++;
			}
		}
	}

}

// Implementation: evalGradient
PDE_HOST PDE_DEVICE PDE_INLINE void LagrangeHex::evalGradient(const Real* xi, Real* dNdxi) const {

	Index pD = 3;

	Real Nx[fem::dispatch::kMaxBasisOrder + 1], Ny[fem::dispatch::kMaxBasisOrder + 1], Nz[fem::dispatch::kMaxBasisOrder + 1];
	Real dNx[fem::dispatch::kMaxBasisOrder + 1], dNy[fem::dispatch::kMaxBasisOrder + 1], dNz[fem::dispatch::kMaxBasisOrder + 1];

	basisX_.eval(xi[0], Nx);
	basisX_.evalFirstDerivative(xi[0], dNx);
	basisY_.eval(xi[1], Ny);
	basisY_.evalFirstDerivative(xi[1], dNy);
	basisZ_.eval(xi[2], Nz);
	basisZ_.evalFirstDerivative(xi[2], dNz);

	// compute tensor product & chain rule
	Index a = 0;
	for (Index k = 0; k <= basisZ_.order(); ++k){
		for (Index j = 0; j <= basisY_.order(); ++j){
			for (Index i = 0; i <= basisX_.order(); ++i){
				dNdxi[a*pD    ] = dNx[i] * Ny[j] * Nz[k];
				dNdxi[a*pD + 1] = Nx[i] * dNy[j] * Nz[k];
				dNdxi[a*pD + 2] = Nx[i] * Ny[j] * dNz[k];
				a++;
			}
		}
	}

}

// Implementation: evalHessian
PDE_HOST PDE_DEVICE PDE_INLINE void LagrangeHex::evalHessian(const Real* xi, Real* d2Nd2xi) const {

	Real Nx[fem::dispatch::kMaxBasisOrder + 1], Ny[fem::dispatch::kMaxBasisOrder + 1], Nz[fem::dispatch::kMaxBasisOrder + 1];
	Real dNx[fem::dispatch::kMaxBasisOrder + 1], dNy[fem::dispatch::kMaxBasisOrder + 1], dNz[fem::dispatch::kMaxBasisOrder + 1];
	Real d2Nx[fem::dispatch::kMaxBasisOrder + 1], d2Ny[fem::dispatch::kMaxBasisOrder + 1], d2Nz[fem::dispatch::kMaxBasisOrder + 1];

	basisX_.eval(xi[0], Nx);
	basisX_.evalFirstDerivative(xi[0], dNx);
	basisX_.evalSecondDerivative(xi[0], d2Nx);
	basisY_.eval(xi[1], Ny);
	basisY_.evalFirstDerivative(xi[1], dNy);
	basisY_.evalSecondDerivative(xi[1], d2Ny);
	basisZ_.eval(xi[2], Nz);
	basisZ_.evalFirstDerivative(xi[2], dNz);
	basisZ_.evalSecondDerivative(xi[2], d2Nz);

	const Index NumEntries = 6;

	// compute tensor product & chain rule
	Index a = 0;
	for (Index k = 0; k <= basisZ_.order(); ++k){
		for (Index j = 0; j <= basisY_.order(); ++j){
			for (Index i = 0; i <= basisX_.order(); ++i){
				d2Nd2xi[a*NumEntries    ] = d2Nx[i] * Ny[j] * Nz[k];
				d2Nd2xi[a*NumEntries + 1] = Nx[i] * d2Ny[j] * Nz[k];
				d2Nd2xi[a*NumEntries + 2] = Nx[i] * Ny[j] * d2Nz[k];
				d2Nd2xi[a*NumEntries + 3] = dNx[i] * dNy[j] * Nz[k];
				d2Nd2xi[a*NumEntries + 4] = dNx[i] * Ny[j] * dNz[k];
				d2Nd2xi[a*NumEntries + 5] = Nx[i] * dNy[j] * dNz[k];
				a++;
			}
		}
	}

}

// Implementation: evalLaplacian
PDE_HOST PDE_DEVICE PDE_INLINE void LagrangeHex::evalLaplacian(const Real* xi, Real* lapN) const {

	Real Nx[fem::dispatch::kMaxBasisOrder + 1], Ny[fem::dispatch::kMaxBasisOrder + 1], Nz[fem::dispatch::kMaxBasisOrder + 1];
	Real d2Nx[fem::dispatch::kMaxBasisOrder + 1], d2Ny[fem::dispatch::kMaxBasisOrder + 1], d2Nz[fem::dispatch::kMaxBasisOrder + 1];

	basisX_.eval(xi[0], Nx);
	basisX_.evalSecondDerivative(xi[0], d2Nx);
	basisY_.eval(xi[1], Ny);
	basisY_.evalSecondDerivative(xi[1], d2Ny);
	basisZ_.eval(xi[2], Nz);
	basisZ_.evalSecondDerivative(xi[2], d2Nz);

	// compute tensor product
	Index a = 0;
	for (Index k = 0; k <= basisZ_.order(); ++k){
		for (Index j = 0; j <= basisY_.order(); ++j){
			for (Index i = 0; i <= basisX_.order(); ++i){
				lapN[a] = d2Nx[i] * Ny[j] * Nz[k] + Nx[i] * d2Ny[j] * Nz[k] + Nx[i] * Ny[j] * d2Nz[k];
				a++;
			}
		}
	}

}

// Implementation: getFaceTopology
PDE_HOST PDE_DEVICE PDE_INLINE void LagrangeHex::getFaceTopology(const Int rngID, Real* nRef) const {

	switch (rngID){
		case 0:
			nRef[0] = -1.0;
			nRef[1] = 0.0;
			nRef[2] = 0.0;
			break;
		case 1:
			nRef[0] = 1.0;
			nRef[1] = 0.0;
			nRef[2] = 0.0;
			break;
		case 2:
			nRef[0] = 0.0;
			nRef[1] = -1.0;
			nRef[2] = 0.0;
			break;
		case 3:
			nRef[0] = 0.0;
			nRef[1] = 1.0;
			nRef[2] = 0.0;
			break;
		case 4:
			nRef[0] = 0.0;
			nRef[1] = 0.0;
			nRef[2] = -1.0;
			break;
		case 5:
			nRef[0] = 0.0;
			nRef[1] = 0.0;
			nRef[2] = 1.0;
			break;
		default:
			nRef[0] = 0.0;
			nRef[1] = 0.0;
			nRef[2] = 0.0;
	}

}

// Implementation: nodesPerFace
PDE_HOST PDE_DEVICE PDE_INLINE Index LagrangeHex::nodesPerFace(const Int rngID) const {

	switch (rngID){
		case 0:
			return (basisZ_.order() + 1)*(basisY_.order() + 1);
		case 1:
			return (basisZ_.order() + 1)*(basisY_.order() + 1);
		case 2:
			return (basisZ_.order() + 1)*(basisX_.order() + 1);
		case 3:
			return (basisZ_.order() + 1)*(basisX_.order() + 1);
		case 4:
			return (basisY_.order() + 1)*(basisX_.order() + 1);
		case 5:
			return (basisY_.order() + 1)*(basisX_.order() + 1);
		default:
			return 0;
	}

}

// Implementation: getFaceNodes
PDE_HOST PDE_DEVICE PDE_INLINE void LagrangeHex::getFaceNodes(const Int rngID, Index* nodeIDs) const {

	switch (rngID){
		case 0:
			for (Index j = 0; j < (basisZ_.order() + 1); ++j) {
				for (Index i = 0; i < (basisY_.order() + 1); ++i)
					nodeIDs[i + j*(basisY_.order() + 1)] = (Index) i*(basisX_.order() + 1) + j*(basisX_.order() + 1)*(basisY_.order() + 1);
			}
			break;
		case 1:
			for (Index j = 0; j < (basisZ_.order() + 1); ++j) {
				for (Index i = 0; i < (basisY_.order() + 1); ++i)
					nodeIDs[i + j*(basisY_.order() + 1)] = (Index) i*(basisX_.order() + 1) + j*(basisX_.order() + 1)*(basisY_.order() + 1) + basisX_.order();
			}
			break;
		case 2:
			for (Index j = 0; j < (basisZ_.order() + 1); ++j) {
				for (Index i = 0; i < (basisX_.order() + 1); ++i)
					nodeIDs[i + j*(basisX_.order() + 1)] = (Index) i + j*(basisX_.order() + 1)*(basisY_.order() + 1);
			}
			break;
		case 3:
			for (Index j = 0; j < (basisZ_.order() + 1); ++j) {
				for (Index i = 0; i < (basisX_.order() + 1); ++i)
					nodeIDs[i + j*(basisX_.order() + 1)] = (Index) i + j*(basisX_.order() + 1)*(basisY_.order() + 1) + (basisX_.order() + 1)*basisY_.order();
			}
			break;
		case 4:
			for (Index j = 0; j < (basisY_.order() + 1); ++j) {
				for (Index i = 0; i < (basisX_.order() + 1); ++i)
					nodeIDs[i + j*(basisX_.order() + 1)] = (Index) i + j*(basisX_.order() + 1);
			}
			break;
		case 5:
			for (Index j = 0; j < (basisY_.order() + 1); ++j) {
				for (Index i = 0; i < (basisX_.order() + 1); ++i)
					nodeIDs[i + j*(basisX_.order() + 1)] = (Index) i + j*(basisX_.order() + 1) + (basisX_.order() + 1)*(basisY_.order() + 1)*basisZ_.order();
			}
			break;
	}

}

// Implementation: mapFaceToElement
PDE_HOST PDE_DEVICE PDE_INLINE void LagrangeHex::mapFaceToElement(const Int rngID, const Real* xi_face, Real* xi_elem) const {

	const Real s = xi_face[0];
	const Real t = xi_face[1];

	switch (rngID) {
		case 0:
			xi_elem[0] = -1.0;
			xi_elem[1] = s;
			xi_elem[2] = t;
			break;
		case 1:
			xi_elem[0] = 1.0;
			xi_elem[1] = s;
			xi_elem[2] = t;
			break;
		case 2:
			xi_elem[0] = s;
			xi_elem[1] = -1.0;
			xi_elem[2] = t;
			break;
		case 3:
			xi_elem[0] = s;
			xi_elem[1] = 1.0;
			xi_elem[2] = t;
			break;
		case 4:
			xi_elem[0] = s;
			xi_elem[1] = t;
			xi_elem[2] = -1.0;
			break;
		case 5:
			xi_elem[0] = s;
			xi_elem[1] = t;
			xi_elem[2] = 1.0;
			break;
	}

}

} // namespace residuum::fem::basis
