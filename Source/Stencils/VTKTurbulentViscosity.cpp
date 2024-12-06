#include "VTKTurbulentViscosityStencil.hpp"

namespace Stencils {

VTKTurbulentViscosityStencil::VTKTurbulentViscosityStencil(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters), written_(false), prefix_("turbulentViscosity") {}

void VTKTurbulentViscosityStencil::apply(FlowField& flowField, int i, int j) {
    RealType nuT = flowField.getTurbulentViscosity().getScalar(i, j);
    viscosityStream_ << nuT << "\n";
}

void VTKTurbulentViscosityStencil::apply(FlowField& flowField, int i, int j, int k) {
    RealType nuT = flowField.getTurbulentViscosity().getScalar(i, j, k);
    viscosityStream_ << nuT << "\n";
}

void VTKTurbulentViscosityStencil::write(FlowField& flowField, int timeStep, RealType simulationTime) {
    openFile(timeStep, simulationTime);

    writeVTKHeader(ofile_);
    writePoints(ofile_, simulationTime);

    ofile_ << "SCALARS TurbulentViscosity float 1\n";
    ofile_ << "LOOKUP_TABLE default\n";
    ofile_ << viscosityStream_.str();

    closeFile();

    viscosityStream_.str("");
    viscosityStream_.clear();
}

void VTKTurbulentViscosityStencil::writeVTKHeader(std::ostream& file) const {
    file << "# vtk DataFile Version 3.0\n";
    file << "Turbulent Viscosity Data\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
}

void VTKTurbulentViscosityStencil::writePoints(std::ostream& file, RealType simulationTime) const {
    file << "POINTS " << parameters_.geometry.sizeX * parameters_.geometry.sizeY * parameters_.geometry.sizeZ
         << " float\n";
    // Add points logic here
}

void VTKTurbulentViscosityStencil::openFile(int timeStep, RealType simulationTime) {
    std::stringstream filename;
    filename << prefix_ << "_" << timeStep << ".vtk";
    ofile_.open(filename.str());
}

void VTKTurbulentViscosityStencil::closeFile() {
    if (ofile_.is_open()) {
        ofile_.close();
    }
}

}
