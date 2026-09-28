namespace residuum::fem::quadrature {

// Implementation: getPoints
PDE_HOST PDE_DEVICE PDE_INLINE void GaussQuadrature1D::getPoints(Real* xi) const {

	if (numPoints_ == 1){
		xi[0] = 0.0;
	}
	else if (numPoints_ == 2){
		xi[0] = -0.5773502691896257;  // -1/sqrt(3)
		xi[1] =  0.5773502691896257;  //  1/sqrt(3)
	}
	else if (numPoints_ == 3){
		xi[0] = -0.7745966692414834;  // -sqrt(3/5)
		xi[1] =  0.0;
		xi[2] =  0.7745966692414834;  //  sqrt(3/5)
	}
	else if (numPoints_ == 4){
		xi[0] = -0.8611363115940526;
		xi[1] = -0.3399810435848563;
		xi[2] =  0.3399810435848563;
		xi[3] =  0.8611363115940526;
	}
	else if (numPoints_ == 5){
		xi[0] = -0.9061798459386640;
		xi[1] = -0.5384693101056831;
		xi[2] =  0.0;
		xi[3] =  0.5384693101056831;
		xi[4] =  0.9061798459386640;
	}
	else {
		// unreachable: numPoints_ is bounds-checked at fem::dispatch::validateDiscretizationLimits
		PDE_UNREACHABLE();
	}

}

// Implementation: getWeights
PDE_HOST PDE_DEVICE PDE_INLINE void GaussQuadrature1D::getWeights(Real* w) const {

	if (numPoints_ == 1){
		w[0] = 2.0;
	}
	else if (numPoints_ == 2){
		w[0] = 1.0;
		w[1] = 1.0;
	}
	else if (numPoints_ == 3){
		w[0] = 0.5555555555555556;  // 5/9
		w[1] = 0.8888888888888889;  // 8/9
		w[2] = 0.5555555555555556;  // 5/9
	}
	else if (numPoints_ == 4){
		w[0] = 0.3478548451374538;
		w[1] = 0.6521451548625461;
		w[2] = 0.6521451548625461;
		w[3] = 0.3478548451374538;
	}
	else if (numPoints_ == 5){
		w[0] = 0.2369268850561891;
		w[1] = 0.4786286704993665;
		w[2] = 0.5688888888888889;
		w[3] = 0.4786286704993665;
		w[4] = 0.2369268850561891;
	}
	else {
		// unreachable: numPoints_ is bounds-checked at fem::dispatch::validateDiscretizationLimits
		PDE_UNREACHABLE();
	}

}

} // namespace residuum::fem::quadrature
