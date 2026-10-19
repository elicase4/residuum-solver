#ifndef RESIDUUM_FEM_BOUNDARY_ESSENTIALBOUNDARYVALUELOWERING_HPP
#define RESIDUUM_FEM_BOUNDARY_ESSENTIALBOUNDARYVALUELOWERING_HPP

#include "core/Types.hpp"

#include "fem/boundary/EssentialBoundaryRegistry.hpp"

#include "mesh/Mesh.hpp"

#include "topology/TopologicalDOF.hpp"

namespace residuum {
	namespace fem {
		namespace boundary {

			// lowers every constrained DOF's Dirichlet value into a flat buffer, once per
			// assembleMatrix/Vector/Diagonal or applyEssentialBCs call rather than once per
			// element -- stays host-only by design (no PDE_HOST PDE_DEVICE), since it still goes
			// through the registry's virtual dispatch; the hot per-element loop in
			// gatherElementVector only ever reads the resulting flat buffer, never the registry.
			//
			// constrainedValues must point to a buffer of size topoDOF.numGlobalDOFs(); only the
			// slots for constrained DOFs are written -- free-DOF slots are left untouched, since
			// gatherElementVector never reads them.
			//
			// Implementation sketch (tutor mode -- fill this in):
			//   for each topoDOF in [0, topoDOF.numGlobalDOFs()):
			//     if (!topoDOF.isConstrained(topoDOF)) continue;
			//     resolve this DOF's node id and component (inverse of getNodeDOF -- see
			//     TopologicalDOF::getDOFNode, and component = topoDOF % dofsPerNode)
			//     look up its constraint tag via topoDOF.getConstraintTag(...)
			//     fetch that tag's entries via bcRegistry.getEntries(tag) (same call the old
			//     gatherElementVector made) and evaluate them at this node's physical
			//     coordinate + time, same as EssentialBoundaryCondition::eval does today
			//     write the resulting component value into constrainedValues[topoDOF]
			template<Index numDOFs>
			void lowerEssentialBoundaryValues(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const EssentialBoundaryRegistry& bcRegistry, Real time, Real* constrainedValues){

				for (Index TdofIDi = 0; TdofIDi < topoDOF.numGlobalDOFs(); ++TdofIDi){

					if (!topoDOF.isConstrained(TdofIDi)) continue;

					Index nodeIDi = topoDOF.getDOFNode(TdofIDi);
					Index componentIDi = (TdofIDi % numDOFs);

					Int rngTag = topoDOF.getConstraintTag(TdofIDi);

					Real bcVal[numDOFs];
					const Real* nodeCoords = mesh.getNodeCoord(nodeIDi);
					const auto* entries = bcRegistry.getEntries(rngTag);

					if (entries) {
						for (const auto& entry : *entries){
							entry->eval(time, nodeCoords, bcVal);
						}
					}

					constrainedValues[TdofIDi] = bcVal[componentIDi];

				}

			}

		} // namespace boundary
	} // namespace fem
} // namespace residuum

#endif
