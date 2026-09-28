#ifndef RESIDUUM_SOLVER_PROBLEM_PROBLEM_HPP
#define RESIDUUM_SOLVER_PROBLEM_PROBLEM_HPP

#include <concepts>
#include <memory>

#include "core/Types.hpp"

#include "fem/assembly/ElementMap.hpp"
#include "linalg/solver/base/LinearSolverRunner.hpp"

namespace residuum {
	namespace solver {
		namespace problem {

			template<typename PT>
			concept Problem = requires(PT& p, Real time, typename PT::MatrixT& K, typename PT::VectorT& V) {

				typename PT::VectorT;
				typename PT::MatrixT;

				{ PT::NumDOFs } -> std::convertible_to<Index>;
				{ PT::NumAuxStates } -> std::convertible_to<Index>;

				p.stiffnessForms();
				p.stiffnessModel();
				p.tangentStiffnessForms();

				{ p.createMatrix() } -> std::same_as<typename PT::MatrixT>;
				{ p.createVector() } -> std::same_as<typename PT::VectorT>;

				{ p.assembleLoad(time) };
				{ p.applyNatural(time) };

				{ p.template assembleMatrix<fem::assembly::GatherMode::Free>(time, p.stiffnessForms(), p.stiffnessModel(), {}, K) };
				{ p.template assembleVector<fem::assembly::GatherMode::Free>(time, p.stiffnessForms(), p.stiffnessModel(), V, {}, V) };
				{ p.template assembleResidual<fem::assembly::GatherMode::Free>(time, p.stiffnessForms(), p.stiffnessModel(), {}, V) };
				{ p.applyEssential(time, p.stiffnessForms(), p.stiffnessModel(), V) };

				{ p.template makeLinearRunner<fem::assembly::GatherMode::Free>(&time, p.stiffnessForms(), p.stiffnessModel(), nullptr, {}, K) } -> std::same_as<std::unique_ptr<linalg::solver::LinearSolverRunner<typename PT::VectorT>>>;

				{ p.numFreeDOFs() } -> std::convertible_to<Index>;
				{ p.solverInstance() };
				{ p.equationLabel() };
				{ p.loggingConfig() };

				{ p.U() } -> std::convertible_to<typename PT::VectorT&>;
				{ p.F() } -> std::convertible_to<typename PT::VectorT&>;

			}; // concept Problem

		} // namespace problem
	} // namespace solver
} // namespace residuum

#endif
