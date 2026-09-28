#ifndef RESIDUUM_TOPOLOGY_TOPOLOGICALDOF_HPP
#define RESIDUUM_TOPOLOGY_TOPOLOGICALDOF_HPP

#include <numeric>
#include <vector>
#include <unordered_map>
#include <set>

#include "core/Types.hpp"
#include "fem/dispatch/DiscretizationLimits.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/dof/DOFOrdering.hpp"
#include "mesh/Mesh.hpp"

namespace residuum {
	namespace topology {

		template<Index numDOFs>
		class TopologicalDOF {
		public:
			
			inline TopologicalDOF(const mesh::Mesh& mesh, fem::dof::DOFOrdering ordering);
			
			// size
			static constexpr Index dofsPerNode = numDOFs;
			Index numGlobalDOFs() const { return numGlobalDOFs_; }
			Index numFreeDOFs() const { return numFreeDOFs_; }
			Index numFreeDOFsPerField() const { return numFreeDOFsPerField_; }
			fem::dof::DOFOrdering ordering() const { return ordering_; }
			Index fieldOffset(Index component) const { return component * numFreeDOFsPerField_; }
			
			//Index numFreeDOFsPerField(Index component) const { return numFreeDOFsPerField_[component]; }
			//Index fieldOffset(Index component) { return std::acculumate(numFreeDOFsPerField_.begin(), numFreeDOFsPerField_.begin() + component, (Index) 0); }			

			// node mappings
			Index getNodeDOF(Index nodeId, Index component) const { return (nodeId * dofsPerNode + component); }
			Index getDOFNode(Index topoDOF) const { return (topoDOF / dofsPerNode); }
			
			// element mappings
			inline void getElementDOFs(Index elemId, Index* dofs) const;
			
			// constraints
			// basis is a runtime instance; only its nodesPerFace()/getFaceNodes() methods are used here
			template<typename BasisT>
			void buildConstraints(const BasisT& basis, const fem::boundary::EssentialBoundaryRegistry& bcRegistry);
			bool isConstrained(Index topoDOF) const { return topoToAlg_[topoDOF] == -1; }
			inline Int getConstraintTag(Index topoDOF) const;

			// mappings
			Index toAlgebraic(Index topoDOF) const { return topoToAlg_[topoDOF]; }
			Int* getTopoToAlg() { return topoToAlg_.data(); }
			Index toTopological(Index algDOF) const { return algToTopo_[algDOF]; }
			Int* getAlgToTopo() { return algToTopo_.data(); }

		private:
			// mesh
			const mesh::Mesh& mesh_;
			
			// dof counts
			Index numGlobalDOFs_;
			Index numFreeDOFs_;
			fem::dof::DOFOrdering ordering_;
			Index numFreeDOFsPerField_;

			//std::vector<Index> numFreeDOFsPerField_;
			
			// mappings
			std::vector<Int> topoToAlg_; // Size: numGlobalDOFs
			std::vector<Int> algToTopo_; // Size: numFreeDOFs

			// constraints
			std::unordered_map<Index, Int> constraintTags_;

		}; // class TopologicalDOF

	} // namespace topology
} // namespace residuum

#include "TopologicalDOF.tpp"

#endif
