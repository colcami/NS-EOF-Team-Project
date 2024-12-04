#include "StdAfx.hpp"
#include "TurbulentViscosity.hpp"
#include "StencilFunctions.hpp"
#include <cmath>
#include "WallDistanceStencil.hpp"

NuTurbulentStencil::NuTurbulentStencil(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters) {}

// Compute the shear rate tensor magnitude (Sij Sij)
RealType NuTurbulentStencil::computeShearRate(const RealType* lv, const RealType* lm) const {
    // Shear rate tensor components
    RealType S11 = dudx(lv, lm); // du/dx
    RealType S22 = dvdy(lv, lm); // dv/dy
    RealType S33 = dwdz(lv, lm); // dw/dz (zero in 2D)

    RealType S12 = 0.5 * (dudy(lv, lm) + dvdx(lv, lm)); // (du/dy + dv/dx)
    RealType S13 = 0.5 * (dudz(lv, lm) + dwdx(lv, lm)); // (du/dz + dw/dx)
    RealType S23 = 0.5 * (dvdz(lv, lm) + dwdy(lv, lm)); // (dv/dz + dw/dy)

    // Compute the magnitude of the shear rate tensor
    return std::sqrt(2.0 * (S11 * S11 + S22 * S22 + S33 * S33 + 2.0 * (S12 * S12 + S13 * S13 + S23 * S23)));
}

// Compute the boundary layer thickness delta
RealType NuTurbulentStencil::computeBoundaryLayerThickness(RealType x, RealType U, RealType nu) const {
    if (x == 0 || U == 0) return 0.0;

    // Compute local Reynolds number
    RealType Rex = computeLocalReynoldsNumber(x, U, nu);

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
RealType NuTurbulentStencil::computeLocalReynoldsNumber(RealType x, RealType U, RealType nu) const {
    return U * x / nu; // Re_x = Ux / nu
}

// Apply stencil for 2D
void NuTurbulentStencil::apply(FlowField& flowField, int i, int j) {
    RealType x = parameters_.meshsize->getPosX(i, j); // x-coordinate
    RealType h = flowField.getWallDistance().getScalar(i, j);    // Distance to nearest wall

    const RealType* lv = flowField.getVelocity().getLocalVector(i, j); // Local velocity
    const RealType* lm = flowField.getMeshsize().getLocalMesh(i, j);   // Local mesh sizes

    RealType U = std::sqrt(lv[0] * lv[0] + lv[1] * lv[1]); // Velocity magnitude
    RealType nu = flowField.getNu().getScalar(i, j);       // Kinematic viscosity (local)

    // Compute boundary layer thickness
    RealType delta = computeBoundaryLayerThickness(x, U, nu);

    // Compute mixing length
    RealType lmScale = std::min(kappa_ * h, parameters_.turbulenceModel.c0 * delta);

    // Compute shear rate
    RealType shearRate = computeShearRate(lv, lm);

    // Compute turbulent viscosity
    RealType nuT = lmScale * lmScale * shearRate;

    // Store turbulent viscosity
    flowField.getTurbulentViscosity().getScalar(i, j) = nuT;
}

// Apply stencil for 3D
void NuTurbulentStencil::apply(FlowField& flowField, int i, int j, int k) {
    RealType x = parameters_.meshsize->getPosX(i, j, k); // x-coordinate
    RealType h = flowField.getWallDistance().getScalar(i, j, k);    // Distance to nearest wall

    const RealType* lv = flowField.getVelocity().getLocalVector(i, j, k); // Local velocity
    const RealType* lm = flowField.getMeshsize().getLocalMesh(i, j, k);   // Local mesh sizes

    RealType U = std::sqrt(lv[0] * lv[0] + lv[1] * lv[1] + lv[2] * lv[2]); // Velocity magnitude
    RealType nu = flowField.getNu().getScalar(i, j, k);                    // Kinematic viscosity (local)

    // Compute boundary layer thickness
    RealType delta = computeBoundaryLayerThickness(x, U, nu);

    // Compute mixing length
    RealType lmScale = std::min(kappa_ * h, parameters_.turbulenceModel.c0 * delta);

    // Compute shear rate
    RealType shearRate = computeShearRate(lv, lm);

    // Compute turbulent viscosity
    RealType nuT = lmScale * lmScale * shearRate;

    // Store turbulent viscosity
    flowField.getTurbulentViscosity().getScalar(i, j, k) = nuT;
}
