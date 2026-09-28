#ifndef RESIDUUM_CORE_LINALG_HPP
#define RESIDUUM_CORE_LINALG_HPP

#include "linalg/operations/MatrixOps.hpp"
#include "linalg/operations/VectorOps.hpp"

#include "linalg/operator/Operator.hpp"
#include "linalg/operator/CSROperator.hpp"
#include "linalg/operator/FEMOperator.hpp"

#include "linalg/types/Matrix.hpp"
#include "linalg/types/CSRMatrix.hpp"
#include "linalg/types/Vector.hpp"

#include "linalg/solver/base/SolverReport.hpp"
#include "linalg/solver/base/LinearSolverAlgorithm.hpp"
#include "linalg/solver/base/LinearSolverRunner.hpp"

#include "linalg/solver/iterative/cg/Config.hpp"
#include "linalg/solver/iterative/cg/Solver.hpp"
#include "linalg/solver/iterative/cg/Workspace.hpp"
#include "linalg/solver/iterative/cg/CGRunner.hpp"

#include "linalg/solver/preconditioner/Identity.hpp"

#endif
