#include "StdAfx.hpp"

#include "VTKWallDistanceStencil.hpp"
#include <cmath>
#include <limits>

// Constructor
Stencils::VTKWallDistanceStencil::VTKWallDistanceStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters), 
   written_(false), 
   prefix_(parameters.vtk.prefix) {}

// Write VTK header
void VTKWallDistanceStencil::writeVTKHeader(std::ostream& file) const {
    file << "# vtk DataFile Version 2.0\nNS-EOF: Wall Distance & Turbulent Viscosity\nASCII\n\n";
}

void Stencils::VTKWallDistanceStencil::writePoints(std::ostream& file, RealType simulationTime) const {
  // Number of points in every direction
  int px = parameters_.parallel.localSize[0] + 1;
  int py = parameters_.parallel.localSize[1] + 1;
  int pz = parameters_.geometry.dim == 2 ? 1 : parameters_.parallel.localSize[2] + 1;

  std::string grid;
  char        buffer[256];

  grid.reserve((file.precision() + 6) * px * py * pz * 3);

  sprintf(
    buffer,
    "DATASET STRUCTURED_GRID\nFIELD FieldData 1\nTIME 1 1 double\n%f\nDIMENSIONS %d %d %d\nPOINTS %d float\n",
    simulationTime,
    px,
    py,
    pz,
    px * py * pz
  );
  grid.append(buffer);

  if (parameters_.geometry.dim == 3) {
    for (int k = 2; k < 2 + pz; k++) {
      for (int j = 2; j < 2 + py; j++) {
        for (int i = 2; i < 2 + px; i++) {
          // Determine positions of grid points at lower/left/front corner of the respective grid cell (i,j,k) -> use
          // meshsize-ptr
          sprintf(
            buffer,
            "%f %f %f\n",
            parameters_.meshsize->getPosX(i, j, k),
            parameters_.meshsize->getPosY(i, j, k),
            parameters_.meshsize->getPosZ(i, j, k)
          );
          grid.append(buffer);
        }
      }
    }
  } else {
    for (int j = 2; j < 2 + py; j++) {
      for (int i = 2; i < 2 + px; i++) {
        sprintf(buffer, "%f %f 0.0\n", parameters_.meshsize->getPosX(i, j), parameters_.meshsize->getPosY(i, j));
        grid.append(buffer);
      }
    }
  }
  grid.append("\n");
  file << grid;
}

// Apply stencil for 2D
void Stencils::VTKWallDistanceStencil::apply(FlowField& flowField, int i, int j) {
    RealType wallDistance = flowField.getWallDistance(i, j); 
    RealType viscosity = flowField.getTurbulentViscosity(i, j); 

    wallDistanceStream_ << wallDistance << "\n";
    viscosityStream_ << viscosity << "\n";
}

// Apply stencil for 3D
void Stencils::VTKWallDistanceStencil::apply(FlowField& flowField, int i, int j, int k) {
    RealType wallDistance = flowField.getWallDistance(i, j, k); 
    RealType viscosity = flowField.getTurbulentViscosity(i, j, k); 

    wallDistanceStream_ << wallDistance << "\n";
    viscosityStream_ << viscosity << "\n";
}

// Open VTK file
void Stencils::VTKWallDistanceStencil::openFile(int timeStep, RealType simulationTime) {
    written_ = false;

    std::stringstream fileNameStream;
    fileNameStream << "Output/" << prefix_ << "/" << prefix_ << "."
                   << parameters_.parallel.rank << "." << timeStep << ".vtk";

    ofile_.open(fileNameStream.str().c_str());
    writeVTKHeader(ofile_);
    writePoints(ofile_, simulationTime);
}

// Write the data to the VTK file
void Stencils::VTKWallDistanceStencil::write(FlowField& flowField, int timeStep, RealType simulationTime) {
    openFile(timeStep, simulationTime);

    // Write cell data for turbulent viscosity and wall distance
    ofile_ << "CELL_DATA " << flowField.getNumberOfCells() << "\n";

    // Wall Distance
    ofile_ << "SCALARS wall_distance float 1\nLOOKUP_TABLE default\n";
    ofile_ << wallDistanceStream_.str();
    wallDistanceStream_.str("");

    // Turbulent Viscosity
    ofile_ << "SCALARS turbulent_viscosity float 1\nLOOKUP_TABLE default\n";
    ofile_ << viscosityStream_.str();
    viscosityStream_.str("");

    closeFile();
}

// Close VTK file
void Stencils::VTKWallDistanceStencil::closeFile() {
    ofile_.close();
}
