#ifndef RESIDUUM_FEM_ASSEMBLY_ASSEMBLER_HPP
#define RESIDUUM_FEM_ASSEMBLY_ASSEMBLER_HPP

#include <algorithm>
#include <array>
#include <cstring>

#include "config/Platform.hpp"

#include "core/Types.hpp"

#include "fem/assembly/ElementMap.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/dispatch/DiscretizationLimits.hpp"
#include "fem/evaluator/EvalElement.hpp"
#include "fem/evaluator/EvalQuadraturePointVolume.hpp"
#include "fem/evaluator/EvalModel.hpp"
#include "fem/form/BilinearForm.hpp"
#include "fem/form/LinearForm.hpp"
#include "fem/form/NonlinearTangentForm.hpp"
#include "fem/form/NonlinearForm.hpp"

#include "mesh/Mesh.hpp"

#include "linalg/types/Matrix.hpp"
#include "linalg/types/CSRMatrix.hpp"
#include "linalg/types/Vector.hpp"
#include "linalg/types/DistributedCSRMatrix.hpp"
#include "linalg/types/DistributedVector.hpp"

#include "topology/TopologicalDOF.hpp"

namespace residuum {
	namespace fem {
		namespace assembly {

			namespace evaluator = residuum::fem::evaluator;

			template<typename BackendT>
			class Assembler {
			public:
			
				// allocation function
				template<Index numDOFs>
				static linalg::types::CSRMatrix<Real, BackendT> createMatrix(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF);
				
				// allocation function
				template<Index numDOFs>
				static linalg::types::Vector<Real, BackendT> createVector(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF);

				template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FormsT, typename QuadratureT, GatherMode Mode>
				requires evaluator::EvalModel<ModelT, EvalQPT>
				static void assembleMatrix(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const ModelT& model, const FormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, BackendT>& U, const std::array<const linalg::types::Vector<Real, BackendT>*, EvalQPT::NumAuxStates>& auxStates, linalg::types::CSRMatrix<Real, BackendT>& K, const fem::boundary::EssentialBoundaryRegistry* bcRegistry);

				template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FormsT, typename QuadratureT, GatherMode Mode>
				requires evaluator::EvalModel<ModelT, EvalQPT>
				static void assembleVector(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const ModelT& model, const FormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, BackendT>& U, const linalg::types::Vector<Real, BackendT>* fieldSource, const std::array<const linalg::types::Vector<Real, BackendT>*, EvalQPT::NumAuxStates>& auxStates, linalg::types::Vector<Real, BackendT>& F, const fem::boundary::EssentialBoundaryRegistry* bcRegistry);

				// scatters only each element's local Ke[a,a] into a plain vector -- no sparsity
				// pattern involved, so this needs no CSRMatrix/colIdx at all
				template<Index numDOFs, evaluator::EvalElement EvalEleT, evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FormsT, typename QuadratureT, GatherMode Mode>
				requires evaluator::EvalModel<ModelT, EvalQPT>
				static void assembleDiagonal(const mesh::Mesh& mesh, const topology::TopologicalDOF<numDOFs>& topoDOF, const Real time, const ModelT& model, const FormsT& forms, const EvalEleT& evalEle, const QuadratureT& quadrature, const linalg::types::Vector<Real, BackendT>* fieldSource, const std::array<const linalg::types::Vector<Real, BackendT>*, EvalQPT::NumAuxStates>& auxStates, linalg::types::Vector<Real, BackendT>& diag, const fem::boundary::EssentialBoundaryRegistry* bcRegistry);

			}; // class Assembler

		} // namespace assembly
	} // namespace fem
} // namespace residuum

#include "backend/cpu/Assembler.tpp"

#endif
