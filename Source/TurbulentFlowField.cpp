#include "StdAfx.hpp"

#include "TurbulentFlowField.hpp"

TurbulentFlowField::TurbulentFlowField(int Nx, int Ny):
  FlowField(Nx, Ny),
  h_(ScalarField(Nx + 3, Ny + 3)),
  nu_t_(ScalarField(Nx + 3, Ny + 3)),
  nu_star_(ScalarField(Nx + 3, Ny + 3)) {

  ASSERTION(Nx > 0);
  ASSERTION(Ny > 0);
}

TurbulentFlowField::TurbulentFlowField(int Nx, int Ny, int Nz):
  FlowField(Nx, Ny, Nz),
  h_(ScalarField(Nx + 3, Ny + 3, Nz + 3)),
  nu_t_(ScalarField(Nx + 3, Ny + 3, Nz + 3)),
  nu_star_(ScalarField(Nx + 3, Ny + 3, Nz + 3)) {
  
  ASSERTION(Nx > 0);
  ASSERTION(Ny > 0);
  ASSERTION(Nz > 0);
}

TurbulentFlowField::TurbulentFlowField(const Parameters& parameters):
  FlowField(parameters),
  h_(
    parameters.geometry.dim == 2 ? ScalarField(getNx() + 3, getNy() + 3) : ScalarField(getNx() + 3, getNy() + 3, getNz() + 3)
  ),
  nu_t_(
    parameters.geometry.dim == 2 ? ScalarField(getNx() + 3, getNy() + 3) : ScalarField(getNx() + 3, getNy() + 3, getNz() + 3)
  ),
  nu_star_(
    parameters.geometry.dim == 2 ? ScalarField(getNx() + 3, getNy() + 3) : ScalarField(getNx() + 3, getNy() + 3, getNz() + 3)
  ){}


ScalarField& TurbulentFlowField::getWallDistance() { return h_; }

ScalarField& TurbulentFlowField::getTurbulentViscosity() { return nu_t_; }

ScalarField& TurbulentFlowField::getTotalViscosity() { return nu_star_; }

