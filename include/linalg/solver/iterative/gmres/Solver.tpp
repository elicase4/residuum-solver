namespace residuum::linalg::solver::iterative::gmres {

	template<typename OperatorT, typename VectorT, typename PreconditionerT, typename LoggerT>
	requires linalg::op::LinearOperator<OperatorT, VectorT>
	bool Solver<OperatorT, VectorT, PreconditionerT, LoggerT>::solve(solver::SolverReport<VectorT>& report, LoggerT& logger, Workspace& W, PreconditionerT& M, const OperatorT& A, const VectorT& b, VectorT& x){

		// reset for fresh solve
		logger.reset();

		// let the preconditioner refresh itself against the current operator
		M.update(A);

		using DataType = typename VectorT::value_type;
		const bool relMode = (config.tolType == ToleranceType::Relative);
		const Index m = config.krylovDim;
		const DataType breakdownTol = DataType(1e-14);
		const DataType eps = DataType(1e-50);

		// column-major Hessenberg access: H(i,j) = row i, col j
		auto Hat = [&](Index i, Index j) -> DataType& { return W.H[static_cast<std::size_t>(j) * (m + 1) + i]; };

		// initial true residual: r0 = b - A*x, built directly into Q[0]
		A.apply(x, W.w); // w = A*x
		operations::copy(b, W.Q[0]); // Q[0] = b
		operations::axpy(DataType(-1), W.w, W.Q[0]); // Q[0] = b - A*x

		const DataType res0 = operations::norm(W.Q[0]);
		report.initialResidual = res0;

		auto converged = [&](DataType res) -> bool {
			DataType crit = relMode ? (res / (res0 + eps)) : res;
			return crit < config.tol;
		};

		// log iteration 0 from the genuine residual vector, matching cg::Solver's convention
		{
			auto perDOF = logger.template computePerDOFNorms<DataType>(W.Q[0].data(), W.Q[0].size());
			logger.log(Index(0), perDOF, DataType(0));
		}

		if (converged(res0)) {
			report.converged = true;
			report.finalResidual = res0;
			report.finalResidualRel = DataType(1);
			report.iterations = 0;
			logger.summary(true);
			return true;
		}

		Index totalIters = 0;
		DataType beta = res0;

		while (totalIters < config.maxIters) {

			// Q[0] currently holds the unnormalized residual for this restart cycle
			operations::scal(DataType(1) / beta, W.Q[0]);

			W.g[0] = beta;
			for (Index i = 1; i <= m; ++i) W.g[i] = DataType(0);

			Index mLocal = 0;
			bool breakdown = false;
			DataType cycleFlops = DataType(0);

			for (Index j = 0; j < m; ++j) {

				// w = A * (M^-1 * Q[j])  -- right preconditioning
				M.apply(W.Q[j], W.z);
				A.apply(W.z, W.w);

				// modified Gram-Schmidt against Q[0..j], recording the Hessenberg column
				for (Index i = 0; i <= j; ++i) {
					DataType hij = operations::dot(W.Q[i], W.w);
					Hat(i, j) = hij;
					operations::axpy(-hij, W.Q[i], W.w);
				}

				DataType hNext = operations::norm(W.w);

				++totalIters;
				breakdown = (hNext <= breakdownTol);

				if (!breakdown) {
					operations::copy(W.w, W.Q[j + 1]);
					operations::scal(DataType(1) / hNext, W.Q[j + 1]);
				}

				// apply the previously-computed Givens rotations to the new column
				for (Index i = 0; i < j; ++i) {
					DataType h1 = Hat(i, j);
					DataType h2 = Hat(i + 1, j);
					Hat(i, j)     =  W.cs[i] * h1 + W.sn[i] * h2;
					Hat(i + 1, j) = -W.sn[i] * h1 + W.cs[i] * h2;
				}

				// compute and apply the new rotation, zeroing the subdiagonal entry
				DataType hjj = Hat(j, j);
				DataType denom = std::sqrt(hjj * hjj + hNext * hNext);
				if (denom <= breakdownTol) {
					W.cs[j] = DataType(1);
					W.sn[j] = DataType(0);
				} else {
					W.cs[j] = hjj / denom;
					W.sn[j] = hNext / denom;
				}
				Hat(j, j) = W.cs[j] * hjj + W.sn[j] * hNext;

				DataType gj = W.g[j];
				W.g[j]     = W.cs[j] * gj;
				W.g[j + 1] = -W.sn[j] * gj;

				DataType resEstimate = std::abs(W.g[j + 1]);

				// this cheap scalar estimate drives convergence/breakdown/budget checks only -- a
				// genuine per-DOF residual isn't available until the true vector is reformed at
				// the end of this cycle (or at convergence), so logging happens there instead of
				// every Arnoldi step; accumulate this step's flops for that eventual log() call
				cycleFlops += static_cast<DataType>(A.flopsPerApply()) + static_cast<DataType>(M.flopsPerApply()) + DataType(4 * (j + 1)) * static_cast<DataType>(x.size());

				mLocal = j + 1;

				if (converged(resEstimate) || breakdown || totalIters >= config.maxIters) break;

			}

			// back-substitute H(0..mLocal-1, 0..mLocal-1) y = g(0..mLocal-1)
			std::vector<DataType> y(mLocal);
			for (Index ii = 0; ii < mLocal; ++ii) {
				Index i = mLocal - 1 - ii;
				DataType sum = W.g[i];
				for (Index k = i + 1; k < mLocal; ++k) sum -= Hat(i, k) * y[k];
				y[i] = sum / Hat(i, i);
			}

			// x += M^-1 * (sum_i y[i] * Q[i])
			W.z.zero();
			for (Index i = 0; i < mLocal; ++i) operations::axpy(y[i], W.Q[i], W.z);
			M.apply(W.z, W.w);
			operations::axpy(DataType(1), W.w, x);

			// fresh true residual for the next restart cycle (or the final report) -- also the
			// only point within a cycle where a genuine residual vector exists, so the per-DOF
			// breakdown is logged here: once per cycle, not once per Arnoldi step
			A.apply(x, W.w);
			operations::copy(b, W.Q[0]);
			operations::axpy(DataType(-1), W.w, W.Q[0]);
			beta = operations::norm(W.Q[0]);

			auto perDOF = logger.template computePerDOFNorms<DataType>(W.Q[0].data(), W.Q[0].size());
			logger.log(totalIters, perDOF, cycleFlops);

			if (converged(beta) || breakdown || totalIters >= config.maxIters) {
				report.converged = converged(beta);
				report.finalResidual = beta;
				report.finalResidualRel = beta / (res0 + eps);
				report.iterations = totalIters;
				logger.summary(report.converged);
				return report.converged;
			}

		}

		report.converged = false;
		report.finalResidual = beta;
		report.finalResidualRel = beta / (res0 + eps);
		report.iterations = totalIters;
		logger.summary(false);
		return false;

	}

} // namespace residuum::linalg::solver::iterative::gmres
