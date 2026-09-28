namespace residuum::fem::quadrature {

// implementation Points
PDE_HOST PDE_DEVICE PDE_INLINE void GaussQuadratureQuad::getPoints(Real* xi) const {

	Real xi_xi[fem::dispatch::kMaxQuadraturePoints1D];
	Real xi_eta[fem::dispatch::kMaxQuadraturePoints1D];

	quadX_.getPoints(xi_xi);
	quadY_.getPoints(xi_eta);

	// compute tensor product
	Index q = 0;
	for (Index j = 0; j < quadY_.numPoints(); ++j) {
		for (Index i = 0; i < quadX_.numPoints(); ++i) {
			xi[2*q]   = xi_xi[i];
			xi[2*q+1] = xi_eta[j];
			q++;
		}
	}
}

// implementation weights
PDE_HOST PDE_DEVICE PDE_INLINE void GaussQuadratureQuad::getWeights(Real* w) const {

	Real w_xi[fem::dispatch::kMaxQuadraturePoints1D];
	Real w_eta[fem::dispatch::kMaxQuadraturePoints1D];

	quadX_.getWeights(w_xi);
	quadY_.getWeights(w_eta);

	// compute tensor product
	Index q = 0;
	for (Index j = 0; j < quadY_.numPoints(); ++j) {
		for (Index i = 0; i < quadX_.numPoints(); ++i) {
			w[q] = w_xi[i] * w_eta[j];
			q++;
		}
	}
}

} // namespace residuum::fem::quadrature
