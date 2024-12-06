#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "StencilFunctions.hpp"
#include "WallDistanceStencil.hpp"


/** Stencil for computing turbulent viscosity (nuT) using Prandtl's mixing length model.
 *
 * This stencil computes the turbulent viscosity based on the mixing length model,
 * which considers local wall distance, boundary layer thickness, and shear rate.
 */
class NuTurbulentStencil : public FieldStencil<FlowField> {
private:
    const RealType kappa_ = 0.41; //! von Kármán constant

    /** Computes the shear rate tensor magnitude S_ij S_ij */
    RealType computeShearRate(const RealType* lv, const RealType* lm) const;

    RealType interpolateVelocityToCellCenter(const FlowField& flowField, int i, int j, int k) const;
    
    /** Computes the boundary layer thickness delta */
    RealType computeBoundaryLayerThickness(RealType x, RealType U, RealType nu) const;

    /** Computes the local Reynolds number based on x, velocity, and viscosity */
    RealType computeLocalReynoldsNumber(RealType x, RealType U, RealType nu) const;

public:
    NuTurbulentStencil(const Parameters& parameters);
    ~NuTurbulentStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;       // For 2D
    void apply(FlowField& flowField, int i, int j, int k) override; // For 3D
};
