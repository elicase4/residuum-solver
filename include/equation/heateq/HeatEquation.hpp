#ifndef RESIDUUM_EQUATION_HEATEQ_HEATEQUATION_HPP
#define RESIDUUM_EQUATION_HEATEQ_HEATEQUATION_HPP

#include "core/Types.hpp"
#include "core/FEM.hpp"

#include "fem/dispatch/DiscretizationDispatch.hpp"
#include "mesh/ElementFamily.hpp"

#include "equation/heateq/boundary/BoundaryFluxFunction.hpp"
#include "equation/heateq/boundary/BoundaryValueFunction.hpp"

#include "equation/heateq/evaluator/SourceFunction.hpp"
#include "equation/heateq/evaluator/DefaultModel.hpp"
#include "equation/heateq/evaluator/ConductivityModel.hpp"
#include "equation/heateq/evaluator/DensityModel.hpp"
#include "equation/heateq/evaluator/SpecificHeatModel.hpp"
#include "equation/heateq/evaluator/EvalElement.hpp"
#include "equation/heateq/evaluator/EvalQuadraturePointVolume.hpp"
#include "equation/heateq/evaluator/EvalQuadraturePointBoundary.hpp"

#include "equation/heateq/form/DiffusionForm.hpp"
#include "equation/heateq/form/MassForm.hpp"
#include "equation/heateq/form/TangentDiffusionForm.hpp"
#include "equation/heateq/form/TangentMassForm.hpp"
#include "equation/heateq/form/SourceForm.hpp"
#include "equation/heateq/form/NodalSourceForm.hpp"
#include "equation/heateq/form/FluxBoundaryForm.hpp"
#include "equation/heateq/form/NodalFluxForm.hpp"

#include "equation/heateq/quantity/HeatFluxIntegrand.hpp"
#include "equation/heateq/quantity/QuantityLimits.hpp"

#include "fem/evaluator/ModelRegistry.hpp"
#include "fem/form/FormRegistry.hpp"
#include "io/fieldio/NodalFileValueSource.hpp"
#include "io/fieldio/NodalValueSourceAdapter.hpp"

namespace residuum {
	namespace equation {

		template<Index NSD, Index NPD, mesh::ElementFamily Family>
		struct HeatEquation {

			// Primary traits
			static constexpr Index NumDOFs = 1;
			static constexpr Index SpatialDim = NSD;

			// Discretization Info
			using ElementTraits = fem::dispatch::ElementTypeTraits<Family>;
			using Basis = typename ElementTraits::BasisType;
			using QuadratureVolumeType = typename ElementTraits::QuadratureVolumeType;
			using QuadratureBoundaryType = typename ElementTraits::QuadratureBoundaryType;

			// Geometry
			using Transform = fem::geometry::JacobianTransform<NSD, NPD>;

			// Element Evaluator
			using EvalEle = heateq::evaluator::EvalElement<Basis, NSD>;
			
			// Quadrature Point Evaluator
			using EvalQPVol = heateq::evaluator::EvalQuadraturePointVolume<EvalEle, Basis, Transform>;
			using EvalQPBdy = heateq::evaluator::EvalQuadraturePointBoundary<EvalEle, Basis, Transform>;

			// Constitutive models
			using DefaultModel = heateq::evaluator::DefaultModel<EvalQPVol>;
			using DefaultModelBdy = heateq::evaluator::DefaultModel<EvalQPBdy>;
			using ConductivityModel = heateq::evaluator::ConductivityModel<EvalQPVol>;
			using ConductivityModelBdy = heateq::evaluator::ConductivityModel<EvalQPBdy>;
			using DensityModel = heateq::evaluator::DensityModel<EvalQPVol>;
			using SpecificHeatModel = heateq::evaluator::SpecificHeatModel<EvalQPVol>;
			using MassModel = fem::evaluator::ModelRegistry<DensityModel, SpecificHeatModel>;

			// Operator forms
			using DiffusionForm = heateq::DiffusionForm<EvalQPVol>;
			using StiffnessForms = fem::form::FormRegistry<DiffusionForm>;
			using MassForm = heateq::MassForm<EvalQPVol>;
			using MassForms = fem::form::FormRegistry<MassForm>;
			using TangentDiffusionForm = heateq::TangentDiffusionForm<EvalQPVol>;
			using TangentDiffusionForms = fem::form::FormRegistry<TangentDiffusionForm>;
			using TangentMassForm = heateq::TangentMassForm<EvalQPVol>;
			using TangentMassForms = fem::form::FormRegistry<TangentMassForm>;

			// Nodal-data sources
			using NodalScalarSource = io::fieldio::NodalFileValueSource<NumDOFs>;
			using NodalFluxSource = io::fieldio::NodalFileValueSource<NumDOFs * SpatialDim>;

			// SourceT: Expression callable vs Nodal data
			template<typename CallableT>
			using SourceFunction = heateq::evaluator::SourceFunction<NSD, NumDOFs, CallableT>;
			template<typename CallableT>
			using SourceForm = heateq::SourceForm<EvalQPVol, SourceFunction<CallableT>>;
			template<typename CallableT>
			using ExpressionSourceForms = fem::form::FormRegistry<SourceForm<CallableT>>;
			template<typename SourceT>
			using NodalSourceForm = heateq::NodalSourceForm<EvalQPVol, SourceT, NumDOFs>;
			template<typename SourceT>
			using NodalSourceForms = fem::form::FormRegistry<NodalSourceForm<SourceT>>;

			// Boundary Flux: Expression callable vs Nodal
			template<typename CallableT>
			using FluxBC = heateq::BoundaryFluxFunction<NSD, NumDOFs, CallableT>;
			template<typename CallableT>
			using FluxBCExpression = FluxBC<CallableT>;
			template<typename SourceT>
			using FluxBCNodal = FluxBC<io::fieldio::NodalValueSourceAdapter<SourceT>>;
			template<typename CallableT>
			using FluxForm = heateq::FluxBoundaryForm<EvalQPBdy, FluxBC<CallableT>>;
			template<typename CallableT>
			using ExpressionFluxForms = fem::form::FormRegistry<FluxForm<CallableT>>;
			template<typename SourceT>
			using NodalFluxForm = heateq::NodalFluxForm<EvalQPBdy, SourceT, NumDOFs, SpatialDim>;
			template<typename SourceT>
			using NodalFluxForms = fem::form::FormRegistry<NodalFluxForm<SourceT>>;

			// Boundary Value: Expression callable vs Nodal data
			template<typename CallableT>
			using DirichletBC = heateq::BoundaryValueFunction<NSD, NumDOFs, CallableT>;
			template<typename CallableT>
			using DirichletExpression = DirichletBC<CallableT>;
			template<typename SourceT>
			using DirichletNodal = DirichletBC<io::fieldio::NodalValueSourceAdapter<SourceT>>;

			// Derived quantities
			using HeatFluxIntegrand = heateq::quantity::HeatFluxIntegrand<EvalQPBdy>;

		}; // struct HeatEquation

	} // namespace equation
} // namespace residuum

#endif
