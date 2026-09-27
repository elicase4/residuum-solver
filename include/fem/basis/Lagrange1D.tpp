namespace residuum::fem::basis {

// Implementation: eval
PDE_HOST PDE_DEVICE PDE_INLINE void Lagrange1D::eval(Real xi, Real* N) const {

	if (order_ == 1){
		N[0] = 0.5 * (1.0 - xi);
		N[1] = 0.5 * (1.0 + xi);
	}
	else if (order_ == 2){
		N[0] = 0.5 * xi * (xi - 1.0);
		N[1] = (1.0 - xi) * (1.0 + xi);
		N[2] = 0.5 * xi * (xi + 1.0);
	}
	else if (order_ == 3){
		N[0] = (-9.0/16.0) * (xi - (1.0/3.0)) * (xi - 1.0) * (xi + (1.0/3.0));
		N[1] = (27.0/16.0) * (xi - 1.0) * (xi - (1.0/3.0)) * (xi + 1.0);
		N[2] = (-27.0/16.0) * (xi + (1.0/3.0)) * (xi - 1.0) * (xi + 1.0);
		N[3] = (9.0/16.0) * (xi - (1.0/3.0)) * (xi + 1.0) * (xi + (1.0/3.0));
	}
	else {
		// unreachable: order_ is bounds-checked at fem::dispatch::validateDiscretizationLimits
		PDE_UNREACHABLE();
	}

}

// Implementation: evalFirstDerivative
PDE_HOST PDE_DEVICE PDE_INLINE void Lagrange1D::evalFirstDerivative(Real xi, Real* N) const {

	if (order_ == 1){
		N[0] = -0.5;
		N[1] = 0.5;
	}
	else if (order_ == 2){
		N[0] = xi - 0.5;
		N[1] = -2.0 * xi;
		N[2] = xi + 0.5;
	}
	else if (order_ == 3){
		N[0] = (-9.0/16.0) * (3.0*xi*xi - 2.0*xi - (1.0/9.0));
		N[1] = (27.0/16.0) * (3.0*xi*xi - (2.0/3.0)*xi - 1.0);
		N[2] = (-27.0/16.0) * (3.0*xi*xi + (2.0/3.0)*xi - 1.0);
		N[3] = (9.0/16.0) * (3.0*xi*xi + 2.0*xi - (1.0/9.0));
	}
	else {
		// unreachable: order_ is bounds-checked at fem::dispatch::validateDiscretizationLimits
		PDE_UNREACHABLE();
	}

}

// Implementation: evalSecondDerivative
PDE_HOST PDE_DEVICE PDE_INLINE void Lagrange1D::evalSecondDerivative(Real xi, Real* N) const {

	if (order_ == 1){
		N[0] = 0.0;
		N[1] = 0.0;
	}
	else if (order_ == 2){
		N[0] = 1.0;
		N[1] = -2.0;
		N[2] = 1.0;
	}
	else if (order_ == 3){
		N[0] = (-9.0/16.0) * (6.0*xi - 2.0);
		N[1] = (27.0/16.0) * (6.0*xi - (2.0/3.0));
		N[2] = (-27.0/16.0) * (6.0*xi + (2.0/3.0));
		N[3] = (9.0/16.0) * (6.0*xi + 2.0);
	}
	else {
		// unreachable: order_ is bounds-checked at fem::dispatch::validateDiscretizationLimits
		PDE_UNREACHABLE();
	}

}

} // namespace residuum::fem::basis
