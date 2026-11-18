// Mesh-refinement (h) convergence study for the steady heat equation on bilinear (Q1) quads.
//
// Manufactured solution on [-1, 1]^2 with k = 1:
//
//   T(x, y) = sin(pi x) sin(pi y) + x y
//   f       = -div(k grad T) = 2 pi^2 sin(pi x) sin(pi y)
//   g       = T on every boundary (reduces to x y, since the sine term vanishes at +-1)
//
// The x y term is harmonic, so it adds nothing to the source but makes the Dirichlet data
// nonzero, which exercises the essential-BC elimination path as well as the interior solve.
//
// For Q1 elements and a smooth solution, a priori FEM theory gives
//
//   || T - T_h ||_L2 = O(h^2)        | T - T_h |_H1 = O(h)
//
// The error norms are integrated element by element with a 4x4 Gauss rule (independent of the
// 2x2 rule used for assembly), using the bilinear interpolant of the nodal solution. Each
// element's corners are identified from their coordinates, so the post-processing does not
// depend on the mesh's local node ordering.
//
// Output: CSV of n, h, free DOFs, CG iterations, L2 error, H1 seminorm error, solve seconds,
// plus observed rates printed between successive refinements.

#include <array>
#include <chrono>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <string>
#include <vector>

#include "core/Config.hpp"
#include "core/FEM.hpp"
#include "core/LinAlg.hpp"
#include "core/Mesh.hpp"
#include "core/Topology.hpp"
#include "core/Types.hpp"
#include "core/Utils.hpp"

#include "mesh/ElementFamily.hpp"
#include "mesh/generator/BlockMesh2D.hpp"

#include "equation/heateq/HeatEquation.hpp"

using namespace residuum;

namespace {

	constexpr Real pi = std::numbers::pi_v<Real>;

	// manufactured solution and its gradient
	Real exactT(Real x, Real y) { return std::sin(pi * x) * std::sin(pi * y) + x * y; }
	Real exactTx(Real x, Real y) { return pi * std::cos(pi * x) * std::sin(pi * y) + y; }
	Real exactTy(Real x, Real y) { return pi * std::sin(pi * x) * std::cos(pi * y) + x; }

	// source and Dirichlet data in the (time, x, out) shape the heat-equation forms expect
	constexpr auto sourceFn = [](Real, const Real* x, Real* out) { out[0] = 2.0 * pi * pi * std::sin(pi * x[0]) * std::sin(pi * x[1]); };
	constexpr auto dirichletFn = [](Real, const Real* x, Real* out) { out[0] = std::sin(pi * x[0]) * std::sin(pi * x[1]) + x[0] * x[1]; };

	// free-DOF gather state, same shape as the one the integration tests use for raw Assembler calls
	template<Index numDOFs, Index SpatialDim, Index ParametricDim, typename VectorT>
	struct FreeGatherState {

		const VectorT* valueSource = nullptr;
		const VectorT* operandSource = nullptr;
		Real Ue[fem::dispatch::kMaxNodesPerElement<ParametricDim> * numDOFs];
		Real Ue_operand[fem::dispatch::kMaxNodesPerElement<ParametricDim> * numDOFs];

		Real* valueBuffer() { return Ue; }
		Real* operandBuffer() { return Ue_operand; }

		template<typename TopoDOFT>
		void gatherElement(const Index* nodeIDs, Index nodesPerElement, const Real* nodeCoords, const TopoDOFT& topoDOF, const Real* constrainedValues, Real time) {
			(void) nodeCoords; (void) constrainedValues; (void) time;
			std::memset(Ue, 0, sizeof(Ue));
			fem::assembly::gatherElementVector<numDOFs, fem::assembly::GatherMode::Free>(nodeIDs, nodesPerElement, topoDOF, Ue, valueSource);
			std::memset(Ue_operand, 0, sizeof(Ue_operand));
			if (operandSource != nullptr) {
				fem::assembly::gatherElementVector<numDOFs, fem::assembly::GatherMode::Free>(nodeIDs, nodesPerElement, topoDOF, Ue_operand, operandSource);
			}
		}

		template<typename QuadraturePointT>
		void interpolate(QuadraturePointT&) const {}

	}; // struct FreeGatherState

	struct CaseResult {
		Index n = 0;
		Real h = 0.0;
		Index freeDOFs = 0;
		Index iterations = 0;
		bool converged = false;
		Real errL2 = 0.0;
		Real errH1 = 0.0;
		Real seconds = 0.0;
	};

	// L2 and H1-seminorm error of the bilinear interpolant of the nodal field, 4x4 Gauss per element
	template<typename TopoDOFT>
	void integrateError(const mesh::Mesh& mesh, const TopoDOFT& topoDOF, const Real* U, CaseResult& out) {

		const std::array<Real, 4> gp = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
		const std::array<Real, 4> gw = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};

		auto nodalValue = [&](Index node) -> Real {
			const Index dof = topoDOF.getNodeDOF(node, 0);
			if (topoDOF.isConstrained(dof)) {
				const Real* xy = mesh.getNodeCoord(node);
				return exactT(xy[0], xy[1]);
			}
			return U[topoDOF.toAlgebraic(dof)];
		};

		Real sumL2 = 0.0;
		Real sumH1 = 0.0;

		for (Index e = 0; e < mesh.data.numElements; ++e) {

			const Index* nodes = mesh.getElementNodes(e);

			// element bounding box from its four corners
			Real xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
			for (Index a = 0; a < 4; ++a) {
				const Real* xy = mesh.getNodeCoord(nodes[a]);
				xmin = std::min(xmin, xy[0]); xmax = std::max(xmax, xy[0]);
				ymin = std::min(ymin, xy[1]); ymax = std::max(ymax, xy[1]);
			}

			// corner values ordered (xmin,ymin), (xmax,ymin), (xmax,ymax), (xmin,ymax)
			std::array<Real, 4> Tc{};
			for (Index a = 0; a < 4; ++a) {
				const Real* xy = mesh.getNodeCoord(nodes[a]);
				const bool right = std::abs(xy[0] - xmax) < std::abs(xy[0] - xmin);
				const bool top = std::abs(xy[1] - ymax) < std::abs(xy[1] - ymin);
				const Index corner = top ? (right ? 2 : 3) : (right ? 1 : 0);
				Tc[corner] = nodalValue(nodes[a]);
			}

			const Real hx = xmax - xmin;
			const Real hy = ymax - ymin;
			const Real jac = 0.25 * hx * hy;

			for (Index i = 0; i < 4; ++i) {
				for (Index j = 0; j < 4; ++j) {

					const Real xi = gp[i];
					const Real eta = gp[j];
					const Real w = gw[i] * gw[j] * jac;

					const Real N0 = 0.25 * (1 - xi) * (1 - eta);
					const Real N1 = 0.25 * (1 + xi) * (1 - eta);
					const Real N2 = 0.25 * (1 + xi) * (1 + eta);
					const Real N3 = 0.25 * (1 - xi) * (1 + eta);

					const Real Th = N0 * Tc[0] + N1 * Tc[1] + N2 * Tc[2] + N3 * Tc[3];
					const Real dTdxi = 0.25 * (-(1 - eta) * Tc[0] + (1 - eta) * Tc[1] + (1 + eta) * Tc[2] - (1 + eta) * Tc[3]);
					const Real dTdeta = 0.25 * (-(1 - xi) * Tc[0] - (1 + xi) * Tc[1] + (1 + xi) * Tc[2] + (1 - xi) * Tc[3]);
					const Real Thx = dTdxi * 2.0 / hx;
					const Real Thy = dTdeta * 2.0 / hy;

					const Real x = xmin + 0.5 * (xi + 1.0) * hx;
					const Real y = ymin + 0.5 * (eta + 1.0) * hy;

					const Real eT = exactT(x, y) - Th;
					const Real ex = exactTx(x, y) - Thx;
					const Real ey = exactTy(x, y) - Thy;

					sumL2 += w * eT * eT;
					sumH1 += w * (ex * ex + ey * ey);
				}
			}
		}

		out.errL2 = std::sqrt(sumL2);
		out.errH1 = std::sqrt(sumH1);
	}

	CaseResult runCase(Index n) {

		using BackendType = linalg::types::backend::CPU;
		using HeatEqBundle = equation::HeatEquation<2, 2, mesh::ElementFamily::Quad>;

		constexpr Real x0 = -1.0, x1 = 1.0, y0 = -1.0, y1 = 1.0;
		constexpr Index Px = 1, Py = 1;      // bilinear (Q1)
		constexpr Index numQuadPoint = 2;    // 2x2 Gauss for assembly
		const Real t = 0.0;

		HeatEqBundle::Basis basis{Px, Py};
		HeatEqBundle::QuadratureVolumeType quadVol{numQuadPoint, numQuadPoint};
		HeatEqBundle::EvalEle evalEle{basis};

		mesh::generator::BlockMesh2D gen{n, n, x0, x1, y0, y1, Px, Py};
		mesh::Mesh mesh2D = gen.generate();
		auto topoDOF = std::make_unique<topology::TopologicalDOF<HeatEqBundle::NumDOFs>>(mesh2D, fem::dof::DOFOrdering::Interleaved);

		// Dirichlet data on all four sides of the block (boundary tags 0-3)
		using DirichletT = HeatEqBundle::DirichletBC<decltype(dirichletFn)>;
		fem::boundary::EssentialBoundaryRegistry bcRegistry;
		std::vector<std::shared_ptr<fem::boundary::BoundaryCondition<DirichletT>>> bcs;
		for (Int tag = 0; tag < 4; ++tag) {
			bcs.push_back(std::make_shared<fem::boundary::BoundaryCondition<DirichletT>>(fem::boundary::BoundaryCondition<DirichletT>{tag, {fem::boundary::BCCategory::Essential}, DirichletT{dirichletFn}}));
			bcRegistry.registerBC<DirichletT>(bcs.back());
		}
		topoDOF->buildConstraints(basis, bcRegistry);

		fem::assembly::Assembler<BackendType> assembler;
		fem::boundary::BoundaryApplicator<BackendType> bcApplicator;

		HeatEqBundle::DefaultModel defaultModel;
		HeatEqBundle::ConductivityModel conductivity;
		conductivity.setConstant(1.0);

		HeatEqBundle::DiffusionForm diffusionForm;
		fem::form::FormRegistry<HeatEqBundle::DiffusionForm> operatorForms{diffusionForm};

		HeatEqBundle::SourceFunction<decltype(sourceFn)> sourceFunction{sourceFn};
		HeatEqBundle::SourceForm<decltype(sourceFn)> sourceForm{sourceFunction};
		fem::form::FormRegistry<HeatEqBundle::SourceForm<decltype(sourceFn)>> rhsForms{sourceForm};

		auto K = assembler.createMatrix<HeatEqBundle::NumDOFs>(mesh2D, *topoDOF);
		auto U = assembler.createVector<HeatEqBundle::NumDOFs>(mesh2D, *topoDOF);
		auto F = assembler.createVector<HeatEqBundle::NumDOFs>(mesh2D, *topoDOF);
		U.zero();

		FreeGatherState<HeatEqBundle::NumDOFs, HeatEqBundle::SpatialDim, HeatEqBundle::EvalEle::ParametricDim, decltype(U)> state;
		state.valueSource = &U;

		const auto start = std::chrono::steady_clock::now();

		assembler.assembleMatrix<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::ConductivityModel, decltype(operatorForms), HeatEqBundle::QuadratureVolumeType, decltype(state)>(mesh2D, *topoDOF, t, conductivity, operatorForms, evalEle, quadVol, state, K, nullptr);

		assembler.assembleVector<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::DefaultModel, HeatEqBundle::QuadratureVolumeType, decltype(state), fem::form::FormsWithOperand<BackendType, decltype(rhsForms), fem::evaluator::SelectValue>>(mesh2D, *topoDOF, t, defaultModel, evalEle, quadVol, state, F, nullptr, fem::form::FormsWithOperand<BackendType, decltype(rhsForms), fem::evaluator::SelectValue>{rhsForms, {}});

		// lift the Dirichlet data into the right-hand side (constrained DOFs are eliminated)
		bcApplicator.applyEssentialBCs<HeatEqBundle::NumDOFs, HeatEqBundle::EvalEle, HeatEqBundle::EvalQPVol, HeatEqBundle::ConductivityModel, decltype(operatorForms), HeatEqBundle::QuadratureVolumeType>(mesh2D, *topoDOF, bcRegistry, t, conductivity, operatorForms, evalEle, quadVol, F);

		// CG with Jacobi; tolerance far below the discretization error so it can't pollute the rates
		using VectorT = linalg::types::Vector<Real, BackendType>;
		linalg::op::CSROperator<linalg::types::CSRMatrix<Real, BackendType>> op(K);
		linalg::solver::iterative::cg::Workspace<VectorT> W(topoDOF->numFreeDOFs());
		linalg::solver::SolverReport<VectorT> report;
		linalg::solver::preconditioner::Jacobi<VectorT> M;
		M.update(op);
		utils::logging::NullLogger logger;
		linalg::solver::iterative::cg::Config<VectorT> cfg{1e-13, linalg::solver::iterative::cg::ToleranceType::Relative, 100000};
		linalg::solver::iterative::cg::Solver<decltype(op), VectorT, decltype(M), decltype(logger)> solver(cfg);
		solver.solve(report, logger, W, M, op, F, U);

		const auto stop = std::chrono::steady_clock::now();

		CaseResult r;
		r.n = n;
		r.h = (x1 - x0) / static_cast<Real>(n);
		r.freeDOFs = topoDOF->numFreeDOFs();
		r.iterations = report.iterations;
		r.converged = report.converged;
		r.seconds = std::chrono::duration<Real>(stop - start).count();
		integrateError(mesh2D, *topoDOF, U.data(), r);
		return r;
	}

} // namespace

int main() {

	const std::vector<Index> resolutions = {4, 8, 16, 32, 64, 128, 256};

	const std::string csvPath = std::string(BENCHMARK_RESULTS_PATH) + "/HeatEq_ConvergenceStudy_Q1.csv";
	std::ofstream csv(csvPath);
	csv << "n,h,free_dofs,cg_iterations,converged,error_L2,error_H1_semi,solve_seconds\n";
	csv << std::setprecision(10);

	std::cout << "Q1 h-convergence, steady heat equation, T = sin(pi x) sin(pi y) + x y on [-1,1]^2\n";
	std::cout << std::left << std::setw(6) << "n" << std::setw(12) << "h" << std::setw(10) << "DOFs" << std::setw(8) << "CG it"
	          << std::setw(14) << "L2 error" << std::setw(8) << "rate" << std::setw(14) << "H1 error" << std::setw(8) << "rate" << "\n";

	CaseResult prev;
	bool havePrev = false;

	for (Index n : resolutions) {

		const CaseResult r = runCase(n);

		csv << r.n << "," << r.h << "," << r.freeDOFs << "," << r.iterations << "," << (r.converged ? 1 : 0) << ","
		    << r.errL2 << "," << r.errH1 << "," << r.seconds << "\n";

		std::cout << std::setw(6) << r.n << std::setw(12) << std::setprecision(4) << r.h << std::setw(10) << r.freeDOFs << std::setw(8) << r.iterations
		          << std::setw(14) << std::scientific << std::setprecision(3) << r.errL2 << std::defaultfloat;
		if (havePrev) {
			std::cout << std::setw(8) << std::fixed << std::setprecision(2) << std::log(prev.errL2 / r.errL2) / std::log(prev.h / r.h) << std::defaultfloat;
		} else {
			std::cout << std::setw(8) << "-";
		}
		std::cout << std::setw(14) << std::scientific << std::setprecision(3) << r.errH1 << std::defaultfloat;
		if (havePrev) {
			std::cout << std::setw(8) << std::fixed << std::setprecision(2) << std::log(prev.errH1 / r.errH1) / std::log(prev.h / r.h) << std::defaultfloat;
		} else {
			std::cout << std::setw(8) << "-";
		}
		std::cout << (r.converged ? "" : "  (CG did not converge)") << "\n";

		prev = r;
		havePrev = true;
	}

	std::cout << "Expected: L2 rate -> 2, H1 rate -> 1\nCSV: " << csvPath << std::endl;
	return 0;
}
