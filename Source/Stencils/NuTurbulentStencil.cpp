#include "StdAfx.hpp"
#include "NuTurbulentStencil.hpp"
#include "StencilFunctions.hpp"
#include <cmath>

namespace Stencils{

NuTurbulentStencil::NuTurbulentStencil(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters) {}

// Compute the shear rate tensor magnitude (Sij Sij)
RealType NuTurbulentStencil::computeShearRate(const RealType* lv, const RealType* lm) const {
    RealType S11 = dudx(lv, lm); // du/dx
    RealType S22 = dvdy(lv, lm); // dv/dy
    RealType S33 = dwdz(lv, lm); // dw/dz (zero in 2D)

    RealType S12 = 0.5 * (dudy(lv, lm) + dvdx(lv, lm)); // (du/dy + dv/dx)
    RealType S13 = 0.5 * (dudz(lv, lm) + dwdx(lv, lm)); // (du/dz + dw/dx)
    RealType S23 = 0.5 * (dvdz(lv, lm) + dwdy(lv, lm)); // (dv/dz + dw/dy)

    return std::sqrt(2.0 * (S11 * S11 + S22 * S22 + S33 * S33 + 2.0 * (S12 * S12 + S13 * S13 + S23 * S23)));
}

// Interpolate velocity to the cell center
RealType NuTurbulentStencil::interpolateVelocityToCellCenter(FlowField& flowField, int i, int j, int k) const {
    RealType u = 0.5 * (flowField.getVelocity().getVector(i, j, k)[0] +
                        flowField.getVelocity().getVector(i - 1, j, k)[0]);
    RealType v = 0.5 * (flowField.getVelocity().getVector(i, j, k)[1] +
                        flowField.getVelocity().getVector(i, j - 1, k)[1]);
    RealType w = 0.5 * (flowField.getVelocity().getVector(i, j, k)[2] +
                        flowField.getVelocity().getVector(i, j, k - 1)[2]);
    return std::sqrt(u * u + v * v + w * w); // Magnitude of velocity at cell center
}

// Compute the boundary layer thickness delta
RealType NuTurbulentStencil::computeBoundaryLayerThickness(RealType x, RealType U, RealType Re) const {
    if (x == 0 || U == 0) return 0.0;

    RealType Rex = computeLocalReynoldsNumber(x, U, Re);

    const std::string& deltaType = parameters_.turbulenceModel.boundaryLayerType;
    if (deltaType == "laminar") {
        return 4.91 * x / std::sqrt(Rex); // Blasius laminar boundary layer thickness
    } else if (deltaType == "turbulent") {
        return 0.382 * x / std::pow(Rex, 1.0 / 5.0); // Turbulent flat plate boundary layer thickness
    } else {
        return 0.0; // No boundary layer
    }
}

// Compute local Reynolds number
RealType NuTurbulentStencil::computeLocalReynoldsNumber(RealType x, RealType U, RealType Re) const {
    return U * x * Re; // Re_x = U * x * Re (since nu = 1/Re)
}

// Apply stencil for 2D
void NuTurbulentStencil::apply(FlowField& flowField, int i, int j) {
    RealType x = parameters_.meshsize->getPosX(i, j); // x-coordinate
    RealType h = flowField.getWallDistance().getScalar(i, j); // Distance to the nearest wall

    // Interpolated velocity at the cell center
    RealType U = interpolateVelocityToCellCenter(flowField, i, j, 0);

    RealType Re = parameters_.flow.Re; // Global Reynolds number

    // Compute boundary layer thickness
    RealType delta = computeBoundaryLayerThickness(x, U, Re);

    // Compute mixing length scale
    RealType lmScale = std::min(kappa_ * h, parameters_.turbulenceModel.c0 * delta);

    // Compute shear rate
    const RealType* lv = flowField.getVelocity().getLocalVector(i, j);
    const RealType* lm = flowField.getMeshsize().getLocalMesh(i, j);
    RealType shearRate = computeShearRate(lv, lm);

    // Compute turbulent viscosity
    RealType nuT = lmScale * lmScale * shearRate;

    // Store turbulent viscosity //! redo this
    flowField.getTurbulentViscosity().getScalar(i, j) = nuT;
}

// Apply stencil for 3D
void NuTurbulentStencil::apply(FlowField& flowField, int i, int j, int k) {
    RealType x = parameters_.meshsize->getPosX(i, j, k); // x-coordinate
    RealType h = flowField.getWallDistance().getScalar(i, j, k); // Distance to the nearest wall

    // Interpolated velocity at the cell center
    RealType U = interpolateVelocityToCellCenter(flowField, i, j, k);

    RealType Re = parameters_.flow.Re; // Global Reynolds number

    // Compute boundary layer thickness
    RealType delta = computeBoundaryLayerThickness(x, U, Re);

    // Compute mixing length scale
    RealType lmScale = std::min(kappa_ * h, parameters_.turbulenceModel.c0 * delta);

    // Compute shear rate
    const RealType* lv = flowField.getVelocity().getLocalVector(i, j, k);
    const RealType* lm = flowField.getMeshsize().getLocalMesh(i, j, k);
    RealType shearRate = computeShearRate(lv, lm);

    // Compute turbulent viscosity
    RealType nuT = lmScale * lmScale * shearRate;

    // Store turbulent viscosity //! redo this
    flowField.getTurbulentViscosity().getScalar(i, j, k) = nuT;
}

} // namespace Stencils
