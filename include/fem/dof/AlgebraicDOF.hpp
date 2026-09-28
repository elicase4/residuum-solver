#ifndef RESIDUUM_FEM_DOF_ALGEBRAICDOF_HPP
#define RESIDUUM_FEM_DOF_ALGEBRAICDOF_HPP

#include "core/Types.hpp"
#include "config/Platform.hpp"

namespace residuum {
	namespace fem {
		namespace dof {

			class AlgebraicDOF {
			public:
				PDE_HOST PDE_DEVICE static Index toAlgebraic(Index topoDOF, const Int* topoToAlg) { return topoToAlg[topoDOF]; }
				PDE_HOST PDE_DEVICE static Index toTopological(Index algDOF, const Int* algToTopo) { return algToTopo[algDOF]; }
				
				PDE_HOST PDE_DEVICE static Index getElementDOFs(const Index* elemTopoDOFs, const Index numTopoDOFs, const Int* topoToAlg, Index* elemAlgDOFs);
			}; // class AlgebraicDOF

		} // namespace dof
	} // namespace fem
} // namespace residuum

#include "AlgebraicDOF.tpp"

#endif
