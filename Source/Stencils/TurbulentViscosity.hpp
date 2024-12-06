/* #pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "StencilFunctions.hpp"
#include "WallDistanceStencil.hpp"

namespace Stencils {

class NuTurbulentStencil : public FieldStencil<FlowField> {
private:
    const RealType kappa_ = 0.41; // von Kármán constant

    RealType computeShearRate(const RealType* lv, const RealType* lm) const;
    RealType computeBoundaryLayerThickness(RealType x, RealType U, RealType nu) const;
    RealType computeLocalReynoldsNumber(RealType x, RealType U, RealType nu) const;

    // Interpolates velocity to the cell center
    RealType interpolateVelocityToCellCenter(const FlowField& flowField, int i, int j, int k) const;

public:
    NuTurbulentStencil(const Parameters& parameters);
    ~NuTurbulentStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;
};

} // namespace Stencils

 */