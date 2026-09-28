namespace residuum::fem::quadrature {

// implementation Points
PDE_HOST PDE_DEVICE PDE_INLINE void GaussQuadratureHex::getPoints(Real* xi) const {

	Real xi_xi[fem::dispatch::kMaxQuadraturePoints1D];
	Real xi_eta[fem::dispatch::kMaxQuadraturePoints1D];
	Real xi_zeta[fem::dispatch::kMaxQuadraturePoints1D];

	quadX_.getPoints(xi_xi);
	quadY_.getPoints(xi_eta);
	quadZ_.getPoints(xi_zeta);

	// compute tensor product
	Index q = 0;
	for (Index k = 0; k < quadZ_.numPoints(); ++k) {
		for (Index j = 0; j < quadY_.numPoints(); ++j) {
			for (Index i = 0; i < quadX_.numPoints(); ++i) {
				xi[3*q]   = xi_xi[i];
				xi[3*q+1] = xi_eta[j];
				xi[3*q+2] = xi_zeta[k];
				q++;
			}
		}
	}
}

// implementation weights
PDE_HOST PDE_DEVICE PDE_INLINE void GaussQuadratureHex::getWeights(Real* w) const {

	Real w_xi[fem::dispatch::kMaxQuadraturePoints1D];
	Real w_eta[fem::dispatch::kMaxQuadraturePoints1D];
	Real w_zeta[fem::dispatch::kMaxQuadraturePoints1D];

	quadX_.getWeights(w_xi);
	quadY_.getWeights(w_eta);
	quadZ_.getWeights(w_zeta);

	// compute tensor product
	Index q = 0;
	for (Index k = 0; k < quadZ_.numPoints(); ++k) {
		for (Index j = 0; j < quadY_.numPoints(); ++j) {
			for (Index i = 0; i < quadX_.numPoints(); ++i) {
				w[q] = w_xi[i] * w_eta[j] * w_zeta[k];
				q++;
			}
		}
	}
}

} // namespace residuum::fem::quadrature
