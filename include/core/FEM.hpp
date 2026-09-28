#ifndef RESIDUUM_CORE_FEM_HPP
#define RESIDUUM_CORE_FEM_HPP

#include "fem/assembly/Assembler.hpp"

#include "fem/boundary/BoundaryApplicator.hpp"
#include "fem/boundary/BoundaryCategory.hpp"
#include "fem/boundary/BoundaryCondition.hpp"
#include "fem/boundary/EssentialBoundaryCondition.hpp"
#include "fem/boundary/NaturalBoundaryOperator.hpp"
#include "fem/boundary/EssentialBoundaryRegistry.hpp"
#include "fem/boundary/NaturalBoundaryRegistry.hpp"

#include "fem/dof/AlgebraicDOF.hpp"
#include "fem/dof/DOFOrdering.hpp"

#include "fem/evaluator/EvalElement.hpp"
#include "fem/evaluator/EvalField.hpp"
#include "fem/evaluator/EvalFunction.hpp"
#include "fem/evaluator/EvalModel.hpp"
#include "fem/evaluator/EvalNodalData.hpp"
#include "fem/evaluator/EvalQuadraturePointVolume.hpp"
#include "fem/evaluator/EvalQuadraturePointBoundary.hpp"

#include "fem/form/FormRegistry.hpp"
#include "fem/form/BilinearForm.hpp"
#include "fem/form/LinearForm.hpp"
#include "fem/form/NonlinearForm.hpp"
#include "fem/form/NonlinearTangentForm.hpp"

#include "fem/geometry/JacobianTransform.hpp"

#include "fem/basis/Lagrange1D.hpp"
#include "fem/basis/LagrangeQuad.hpp"
#include "fem/basis/LagrangeHex.hpp"

#include "fem/quadrature/GaussQuadrature1D.hpp"
#include "fem/quadrature/GaussQuadratureQuad.hpp"
#include "fem/quadrature/GaussQuadratureHex.hpp"

#endif
