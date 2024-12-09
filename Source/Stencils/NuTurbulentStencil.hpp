#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "StencilFunctions.hpp"
#include "WallDistanceStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "BoundaryStencil.hpp"


/** Stencil for computing turbulent viscosity (nuT) using Prandtl's mixing length model.
 *
 * This stencil computes the turbulent viscosity based on the mixing length model,
 * which considers local wall distance, boundary layer thickness, and shear rate.
 */
namespace Stencils{
class NuTurbulentStencil : public FieldStencil<TurbulentFlowField> {
private:
    /** Computes the shear rate tensor magnitude S_ij S_ij */
    RealType localVelocity_[27*3];
    RealType localMeshSize_[27*3];
    
    /** Computes the boundary layer thickness delta */
    RealType computeBoundaryLayerThickness(RealType x, RealType U, RealType nu) const;

    /** Computes the local Reynolds number based on x, velocity, and viscosity */
    RealType computeLocalReynoldsNumber(RealType x, RealType U, RealType nu) const;

public:
    NuTurbulentStencil(const Parameters& parameters);
    ~NuTurbulentStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;       // For 2D
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;  // For 3D
};

}// namespace Stencils