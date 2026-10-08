#ifndef RESIDUUM_LINALG_OPERATOR_FEMOPERATOR_HPP
#define RESIDUUM_LINALG_OPERATOR_FEMOPERATOR_HPP

#include "fem/assembly/Assembler.hpp"
#include "fem/assembly/ElementMap.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"

#include "fem/evaluator/EvalElement.hpp"
#include "fem/evaluator/EvalQuadraturePointVolume.hpp"
#include "fem/evaluator/EvalModel.hpp"
#include "fem/evaluator/StateVector.hpp"

#include "fem/form/BilinearForm.hpp"
#include "fem/form/NonlinearTangentForm.hpp"

#include "mesh/Mesh.hpp"

#include "linalg/operator/Operator.hpp"
#include "linalg/types/Vector.hpp"

#include "topology/TopologicalDOF.hpp"

namespace residuum {
	namespace linalg {
		namespace op {

			template<typename AssemblerT, typename TopologicalDOFT, fem::evaluator::EvalElement EvalEleT, fem::evaluator::EvalQuadraturePointVolume EvalQPT, typename ModelT, typename FormsT, typename QuadratureT, typename StateVectorT, typename VectorT>
			requires fem::evaluator::EvalModel<ModelT, EvalQPT> && fem::evaluator::StateVector<StateVectorT, TopologicalDOFT::dofsPerNode, EvalQPT>
			class FEMOperator {
			public:

				// state is a reference to whichever StateVector instance the caller already gathers/interpolates for this solve
				FEMOperator(const AssemblerT& assembler_, const mesh::Mesh& mesh_, const TopologicalDOFT& topoDOF_, const Real* time_, const ModelT& model_, const FormsT& forms_, const EvalEleT& evalEle_, const QuadratureT& quadrature_, const fem::boundary::EssentialBoundaryRegistry* bcRegistry_, StateVectorT& state_) : assembler(assembler_), mesh(mesh_), topoDOF(topoDOF_), time(time_), model(model_), forms(forms_), quadrature(quadrature_), evalEle(evalEle_), bcRegistry(bcRegistry_), state(state_) {
					static_assert(LinearOperator<FEMOperator, VectorT>);
					static_assert(DiagonalExtractable<FEMOperator, VectorT>);
				}

				// a Krylov basis vector x is always free-dof-only, so it is bound as state's operandSource
				// and gathered via GatherMode::Free inside StateVector::gatherElement -- never Full/Constrained
				void apply(const VectorT& x, VectorT& y) const {
					state.operandSource = &x;
					auto selectOperand = [](StateVectorT& s) -> Real* { return s.operandBuffer(); };
					assembler.template assembleVector<TopologicalDOFT::dofsPerNode, EvalEleT, EvalQPT, ModelT, QuadratureT, StateVectorT, fem::assembly::FormsWithOperand<FormsT, decltype(selectOperand)>>(mesh, topoDOF, *time, model, evalEle, quadrature, state, y, bcRegistry, fem::assembly::FormsWithOperand<FormsT, decltype(selectOperand)>{forms, selectOperand});
				}

				Index size() const {
					return topoDOF.numFreeDOFs();
				}

				void diagonal(VectorT& d) const {
					assembler.template assembleDiagonal<TopologicalDOFT::dofsPerNode, EvalEleT, EvalQPT, ModelT, FormsT, QuadratureT, StateVectorT>(mesh, topoDOF, *time, model, forms, evalEle, quadrature, state, d, bcRegistry);
				}

				// Approximation by dominant local dense-matvec term only
				Index flopsPerApply() const {
					return 2 * mesh.data.numElements * mesh.data.nodesPerElement * mesh.data.nodesPerElement;
				}

			private:
				const AssemblerT& assembler;
				const mesh::Mesh& mesh;
				const TopologicalDOFT& topoDOF;
				const Real* time;
				const ModelT& model;
				const FormsT& forms;
				const QuadratureT& quadrature;
				const EvalEleT& evalEle;

				const fem::boundary::EssentialBoundaryRegistry* bcRegistry;
				StateVectorT& state;

			}; // class FEMOperator

		} // namespace op
	} // namespace linalg
} // namespace residuum

#endif
