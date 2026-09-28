#include <cmath>
#include <cassert>
#include <vector>

#include "linalg/solver/iterative/cg/Config.hpp"
#include "linalg/solver/iterative/cg/Workspace.hpp"

#include "linalg/operations/VectorOps.hpp"
#include "linalg/solver/base/SolverReport.hpp"

namespace residuum::linalg::solver::iterative::cg {

	template<typename OperatorT, typename VectorT, typename PreconditionerT, typename LoggerT>
	requires linalg::op::LinearOperator<OperatorT, VectorT>
	bool Solver<OperatorT, VectorT, PreconditionerT, LoggerT>::solve(solver::SolverReport<VectorT>& report, LoggerT& logger, Workspace& W, PreconditionerT& M, const OperatorT& A, const VectorT& b, VectorT& x){

		// a fresh solve, so the logger's relative-residual baseline and per-call totals must not carry over
		logger.reset();

		// let the preconditioner refresh itself against the current operator (e.g. a Jacobi
		// preconditioner re-extracting the diagonal after T changed between Newton iterations)
		M.update(A);

		// get config info
		using DataType = typename VectorT::value_type;
		const bool relMode = (config.tolType == ToleranceType::Relative);

		// per-iteration flop cost is constant for CG
		const DataType flopsPerIter = static_cast<DataType>(A.flopsPerApply()) + static_cast<DataType>(M.flopsPerApply()) + DataType(10) * static_cast<DataType>(x.size());

		// compute intial residual
		A.apply(x, W.Ap); // Ap = A*x
		operations::copy(b, W.r); // r = b
		operations::axpy(DataType(-1.0), W.Ap, W.r); // r = b - Ap

		// compute absolute & relative residual magnitude
		const DataType res0 = operations::norm(W.r); // ||r||
		report.initialResidual = res0;

		// log iteration 0
		auto perDOF = logger.template computePerDOFNorms<DataType>(W.r.data(), W.r.size());
		logger.log(Index(0), perDOF, DataType(0));

		// convergence lambda
		auto converged = [&](DataType res) -> bool {
			DataType crit = relMode ? (res / (res0 + DataType(1e-50))) : res;
			return crit < config.tol;
		};

		// check convergence
		if (converged(res0)){
			report.converged = true;
			report.finalResidual = res0;
			report.finalResidualRel = DataType(1);
			report.iterations = 0;
			logger.summary(true);
			return true;
		}

		// apply preconditioner
		M.apply(W.r, W.z); // z = M*r
		operations::copy(W.z, W.p); // p = z
		DataType rz_old = operations::dot(W.r, W.z); // rz_old = r * z

		// pcg solver loop
		for (Index k = 1; k < config.maxIters; ++k){
			
			// operator application
			A.apply(W.p, W.Ap); // Ap = A * p

			// compute alpha
			DataType pAp = operations::dot(W.p, W.Ap); // pAp = p*Ap
			DataType alpha = rz_old / pAp;// alpha = r*z / p*Ap
			operations::axpy(alpha, W.p, x); // x += alpha * p
			operations::axpy(-alpha, W.Ap, W.r); // r -= alpha * Ap

			// compute absolute & relative residual magnitude
			DataType res = operations::norm(W.r); // res = ||r||
			DataType rel = res / (res0 + DataType(1e-50));

			// log iteration
			auto perDOF = logger.template computePerDOFNorms<DataType>(W.r.data(), W.r.size());
			logger.log(k, perDOF, flopsPerIter);

			// check convergence
			if (converged(res)){
				report.converged = true;
				report.finalResidual = res;
				report.finalResidualRel = rel;
				report.iterations = k;
				report.perFieldResidual = std::vector<DataType>(perDOF.begin(), perDOF.end());
				logger.summary(true);
				return true;
			}

			// apply preconditioner
			M.apply(W.r, W.z);

			// update quantities for new iter
			DataType rz_new = operations::dot(W.r, W.z); // rz_new = r*z
			DataType beta = rz_new / rz_old;

			operations::scal(beta, W.p); // p *= beta
			operations::axpy(DataType(1), W.z, W.p); // p += z

			rz_old = rz_new;

		}

		// report if max iterations are reached
		const DataType res_final = operations::norm(W.r);
		report.converged = false;
		report.finalResidual = res_final;
		report.finalResidualRel = res_final / (res0 + DataType(1e-50));
		report.iterations = config.maxIters;
		logger.summary(false);
		return false;

	}

} // namespace residuum::linalg::solver::iterative::cg
