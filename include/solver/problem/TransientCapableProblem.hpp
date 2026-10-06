#ifndef RESIDUUM_SOLVER_PROBLEM_TRANSIENTCAPABLEPROBLEM_HPP
#define RESIDUUM_SOLVER_PROBLEM_TRANSIENTCAPABLEPROBLEM_HPP

#include <concepts>

#include "core/Types.hpp"
#include "solver/problem/Problem.hpp"

namespace residuum {
	namespace solver {
		namespace problem {

			// U_prev/Udot as separate named accessors is a stopgap, not general -- to be replaced by a
			// compile-time state-vector package (see the state-vector-unification issue) once a second
			// equation needs more than one rate-like quantity
			template<typename PT>
			concept TransientCapableProblem = Problem<PT> && requires(PT& p) {

				p.massForms();
				p.massModel();
				p.tangentMassForms();

				{ p.U_prev() } -> std::convertible_to<typename PT::VectorT&>;
				{ p.Udot() } -> std::convertible_to<typename PT::VectorT&>;

			}; // concept TransientCapableProblem

		} // namespace problem
	} // namespace solver
} // namespace residuum

#endif
