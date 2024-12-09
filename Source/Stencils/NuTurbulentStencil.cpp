#include "StdAfx.hpp"
#include "NuTurbulentStencil.hpp"
#include "StencilFunctions.hpp"
#include <cmath>

namespace Stencils{

NuTurbulentStencil::NuTurbulentStencil(const Parameters& parameters)
    : FieldStencil<TurbulentFlowField>(parameters) {}

// Compute the shear rate tensor magnitude (Sij Sij)

// Interpolate velocity to the cell center
// RealType NuTurbulentStencil::interpolateVelocityToCellCenter(FlowField& flowField, int i, int j, int k) const {
//     RealType u = 0.5 * (flowField.getVelocity().getVector(i, j, k)[0] +
//                         flowField.getVelocity().getVector(i - 1, j, k)[0]);
//     RealType v = 0.5 * (flowField.getVelocity().getVector(i, j, k)[1] +
//                         flowField.getVelocity().getVector(i, j - 1, k)[1]);
//     RealType w = 0.5 * (flowField.getVelocity().getVector(i, j, k)[2] +
//                         flowField.getVelocity().getVector(i, j, k - 1)[2]);
//     return std::sqrt(u * u + v * v + w * w); // Magnitude of velocity at cell center
// }

// Compute the boundary layer thickness delta
//! check if div/0
RealType NuTurbulentStencil::computeBoundaryLayerThickness(RealType x, RealType L , RealType Re) const {
    if (x == 0) return 0.0;

    RealType Rex = computeLocalReynoldsNumber(x, L, Re);

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
RealType NuTurbulentStencil::computeLocalReynoldsNumber(RealType x, RealType L , RealType Re) const {
    return x * Re / L; // Re_x = U * x * Re (since nu = 1/Re)
}

// Apply stencil for 2D
void NuTurbulentStencil::apply(TurbulentFlowField& flowField, int i, int j) {
    RealType x = parameters_.meshsize->getPosX(i, j) + 0.5 * parameters_.meshsize->getDx(i,j) ; // x-coordinate
    RealType h = flowField.getWallDistance().getScalar(i, j); // Distance to the nearest wall
    RealType L = parameters_.geometry.lengthX;

    // Interpolated velocity at the cell center
    // RealType U = interpolateVelocityToCellCenter(flowField, i, j, 0);

    RealType Re = parameters_.flow.Re; // Global Reynolds number

    // Compute boundary layer thickness
    RealType delta = computeBoundaryLayerThickness(x, L, Re);

    // Compute mixing length scale
    RealType lmScale = std::min(parameters_.turbulenceModel.kappa * h, parameters_.turbulenceModel.c0 * delta);

    // Compute shear rate
    RealType localVelocity[27*3];
    RealType localMeshSize[27*3];
    loadLocalVelocity2D(flowField, localVelocity, i, j);
    loadLocalMeshsize2D(parameters_, localMeshSize, i, j);

    RealType shearRate = computeShearRate(localVelocity, localMeshSize);

    // Compute turbulent viscosity
    RealType nuT = lmScale * lmScale * shearRate;

    // Store turbulent viscosity //! redo this
    flowField.getTurbulentViscosity().getScalar(i, j) = nuT;
}

// Apply stencil for 3D
void NuTurbulentStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {
    RealType x = parameters_.meshsize->getPosX(i, j, k) + 0.5 * parameters_.meshsize->getDx(i,j,k); // x-coordinate
    RealType h = flowField.getWallDistance().getScalar(i, j, k); // Distance to the nearest wall
    RealType L = parameters_.geometry.lengthX;

    // Interpolated velocity at the cell center
    // RealType U = interpolateVelocityToCellCenter(flowField, i, j, k);

    RealType Re = parameters_.flow.Re; // Global Reynolds number

    // Compute boundary layer thickness
    RealType delta = computeBoundaryLayerThickness(x, L, Re);

    // Compute mixing length scale
    RealType lmScale = std::min(parameters_.turbulenceModel.kappa * h, parameters_.turbulenceModel.c0 * delta);

    // Compute shear rate
    RealType localVelocity[27*3];
    RealType localMeshSize[27*3];
    loadLocalVelocity3D(flowField, localVelocity, i, j, k);
    loadLocalMeshsize3D(parameters_, localMeshSize, i, j, k);

    RealType shearRate = computeShearRate(localVelocity, localMeshSize);

    // Compute turbulent viscosity
    RealType nuT = lmScale * lmScale * shearRate;

    // Store turbulent viscosity //! redo this
    flowField.getTurbulentViscosity().getScalar(i, j, k) = nuT;
}

} // namespace Stencils
