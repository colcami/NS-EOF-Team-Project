/* 
#include "StdAfx.hpp"      
#include "TurbulentViscosity.hpp"
#include "StencilFunctions.hpp"
#include <cmath>
#include "WallDistanceStencil.hpp"
#include "FieldStencil.hpp"

namespace Stencils {

NuTurbulentStencil::NuTurbulentStencil(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters) {}

RealType NuTurbulentStencil::computeShearRate(const RealType* lv, const RealType* lm) const {
    RealType S11 = dudx(lv, lm); // du/dx
    RealType S22 = dvdy(lv, lm); // dv/dy
    RealType S33 = dwdz(lv, lm); // dw/dz

    RealType S12 = 0.5 * (dudy(lv, lm) + dvdx(lv, lm)); // (du/dy + dv/dx)
    RealType S13 = 0.5 * (dudz(lv, lm) + dwdx(lv, lm)); // (du/dz + dw/dx)
    RealType S23 = 0.5 * (dvdz(lv, lm) + dwdy(lv, lm)); // (dv/dz + dw/dy)

    return std::sqrt(2.0 * (S11 * S11 + S22 * S22 + S33 * S33 + 2.0 * (S12 * S12 + S13 * S13 + S23 * S23)));
}

RealType NuTurbulentStencil::interpolateVelocityToCellCenter(const FlowField& flowField, int i, int j, int k) const {
    const RealType* velocityHere = flowField.getVelocity().getVector(i, j, k);
    const RealType* velocityLeft = flowField.getVelocity().getVector(i - 1, j, k);
    const RealType* velocityDown = flowField.getVelocity().getVector(i, j - 1, k);
    const RealType* velocityBack = flowField.getVelocity().getVector(i, j, k - 1);

    RealType u = 0.5 * (velocityHere[0] + velocityLeft[0]);
    RealType v = 0.5 * (velocityHere[1] + velocityDown[1]);
    RealType w = 0.5 * (velocityHere[2] + velocityBack[2]);

    return std::sqrt(u * u + v * v + w * w);
}

RealType NuTurbulentStencil::computeBoundaryLayerThickness(RealType x, RealType U, RealType nu) const {
    if (x == 0 || U == 0) return 0.0;

    RealType Rex = computeLocalReynoldsNumber(x, U, nu);
    const std::string& deltaType = parameters_.turbulenceModel.boundaryLayerType;

    if (deltaType == "laminar") {
        return 4.91 * x / std::sqrt(Rex);
    } else if (deltaType == "turbulent") {
        return 0.382 * x / std::pow(Rex, 1.0 / 5.0);
    } else {
        return 0.0;
    }
}

RealType NuTurbulentStencil::computeLocalReynoldsNumber(RealType x, RealType U, RealType nu) const {
    return U * x / nu;
}

void NuTurbulentStencil::apply(FlowField& flowField, int i, int j, int k) {
    RealType x = parameters_.meshsize->getPosX(i, j, k);
    RealType h = flowField.getWallDistance().getScalar(i, j, k);

    RealType U = interpolateVelocityToCellCenter(flowField, i, j, k);
    RealType nu = 1.0 / parameters_.flow.Re;

    RealType delta = computeBoundaryLayerThickness(x, U, nu);
    RealType lmScale = std::min(kappa_ * h, parameters_.turbulenceModel.c0 * delta);

    const RealType* lv = flowField.getVelocity().getVector(i, j, k);
    const RealType* lm = flowField.getMeshsize()->getLocalMesh(i, j, k);

    RealType shearRate = computeShearRate(lv, lm);
    RealType nuT = lmScale * lmScale * shearRate;

    flowField.getTurbulentViscosity().getScalar(i, j, k) = nuT;
}

void NuTurbulentStencil::apply(FlowField& flowField, int i, int j) {
    RealType x = parameters_.meshsize->getPosX(i, j);
    RealType h = flowField.getWallDistance().getScalar(i, j);

    RealType U = interpolateVelocityToCellCenter(flowField, i, j, 0);
    RealType nu = 1.0 / parameters_.flow.Re;

    RealType delta = computeBoundaryLayerThickness(x, U, nu);
    RealType lmScale = std::min(kappa_ * h, parameters_.turbulenceModel.c0 * delta);

    const RealType* lv = flowField.getVelocity().getVector(i, j);
    const RealType* lm = flowField.getMeshsize()->getLocalMesh(i, j);

    RealType shearRate = computeShearRate(lv, lm);
    RealType nuT = lmScale * lmScale * shearRate;

    flowField.getTurbulentViscosity().getScalar(i, j) = nuT;
}

} // namespace Stencils
 */