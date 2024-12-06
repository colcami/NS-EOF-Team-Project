#include "VTKWallDistanceStencil.hpp"

namespace Stencils {

VTKWallDistanceStencil::VTKWallDistanceStencil(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters), written_(false), prefix_("wallDistance") {}

void VTKWallDistanceStencil::apply(FlowField& flowField, int i, int j) {
    RealType h = flowField.getWallDistance().getScalar(i, j);
    wallDistanceStream_ << h << "\n";
}

void VTKWallDistanceStencil::apply(FlowField& flowField, int i, int j, int k) {
    RealType h = flowField.getWallDistance().getScalar(i, j, k);
    wallDistanceStream_ << h << "\n";
}

void VTKWallDistanceStencil::write(FlowField& flowField, int timeStep, RealType simulationTime) {
    openFile(timeStep, simulationTime);

    writeVTKHeader(ofile_);
    writePoints(ofile_, simulationTime);

    ofile_ << "SCALARS WallDistance float 1\n";
    ofile_ << "LOOKUP_TABLE default\n";
    ofile_ << wallDistanceStream_.str();

    closeFile();

    wallDistanceStream_.str("");
    wallDistanceStream_.clear();
}

void VTKWallDistanceStencil::writeVTKHeader(std::ostream& file) const {
    file << "# vtk DataFile Version 3.0\n";
    file << "Wall Distance Data\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
}

void VTKWallDistanceStencil::writePoints(std::ostream& file, RealType simulationTime) const {
    file << "POINTS " << parameters_.geometry.sizeX * parameters_.geometry.sizeY * parameters_.geometry.sizeZ
         << " float\n";
    // Add points logic here
}

void VTKWallDistanceStencil::openFile(int timeStep, RealType simulationTime) {
    std::stringstream filename;
    filename << prefix_ << "_" << timeStep << ".vtk";
    ofile_.open(filename.str());
}

void VTKWallDistanceStencil::closeFile() {
    if (ofile_.is_open()) {
        ofile_.close();
    }
}

}
