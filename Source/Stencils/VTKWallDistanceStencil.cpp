#include "StdAfx.hpp"
#include "VTKWallDistanceStencil.hpp"
#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "WallDistanceStencil.hpp"
#include <sstream>
#include <fstream>
#include <filesystem>
#include <stdexcept>

namespace Stencils {

VTKWallDistanceStencil::VTKWallDistanceStencil(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters), written_(false), prefix_("wallDistance") {}

void VTKWallDistanceStencil::apply(FlowField& flowField, int i, int j) {
    // Extract and append wall distance for 2D
    RealType h = flowField.getWallDistance().getScalar(i, j);
    wallDistanceStream_ << h << "\n";
}

void VTKWallDistanceStencil::apply(FlowField& flowField, int i, int j, int k) {
    // Extract and append wall distance for 3D
    RealType h = flowField.getWallDistance().getScalar(i, j, k);
    wallDistanceStream_ << h << "\n";
}

void VTKWallDistanceStencil::write(FlowField& flowField, int timeStep, RealType simulationTime) {
    // Open the output file
    openFile(timeStep, simulationTime);

    if (wallDistanceStream_.str().empty()) {
        std::cerr << "Warning: Wall distance stream is empty for time step " << timeStep << ". Skipping file writing.\n";
        return;
    }

    // Write VTK file header
    writeVTKHeader(ofile_);
    writePoints(ofile_, simulationTime);

    // Write the scalar data for wall distance
    int numCells = parameters_.parallel.localSize[0] * parameters_.parallel.localSize[1] *
                   (parameters_.geometry.dim == 3 ? parameters_.parallel.localSize[2] : 1);

    ofile_ << "CELL_DATA " << numCells << "\n";
    ofile_ << "SCALARS WallDistance float 1\n";
    ofile_ << "LOOKUP_TABLE default\n";
    ofile_ << wallDistanceStream_.str(); // Output wall distance values

    // Close the file
    closeFile();

    // Reset the stream
    wallDistanceStream_.str("");
    wallDistanceStream_.clear();
}

void VTKWallDistanceStencil::writeVTKHeader(std::ostream& file) const {
    file << "# vtk DataFile Version 3.0\n";
    file << "Wall Distance Data\n";
    file << "ASCII\n\n";
}

void VTKWallDistanceStencil::writePoints(std::ostream& file, RealType simulationTime) const {
    file << "DATASET STRUCTURED_GRID\n";

    int sizeX = parameters_.parallel.localSize[0];
    int sizeY = parameters_.parallel.localSize[1];
    int sizeZ = (parameters_.geometry.dim == 3) ? parameters_.parallel.localSize[2] : 1;

    file << "DIMENSIONS " << sizeX + 1 << " " << sizeY + 1 << " " << sizeZ + 1 << "\n";
    file << "POINTS " << (sizeX + 1) * (sizeY + 1) * (sizeZ + 1) << " float\n";

    for (int k = 0; k <= sizeZ; ++k) {
        for (int j = 0; j <= sizeY; ++j) {
            for (int i = 0; i <= sizeX; ++i) {
                file << parameters_.meshsize->getPosX(i, j, k) << " "
                     << parameters_.meshsize->getPosY(i, j, k) << " "
                     << (parameters_.geometry.dim == 3 ? parameters_.meshsize->getPosZ(i, j, k) : 0.0) << "\n";
            }
        }
    }
}

void VTKWallDistanceStencil::openFile(int timeStep, RealType simulationTime) {
    const std::string folder = "Outputh";

    // Ensure the directory exists
    createDirectory(folder);

    // Incorporate rank into the filename to ensure unique file names in parallel simulations
    std::string scenario = parameters_.simulation.scenario;
    std::replace(scenario.begin(), scenario.end(), ' ', '_');

    std::stringstream filename;
    filename << folder << "/" << prefix_ << "_" << scenario << "_ws2_24_" << parameters_.parallel.rank << "_"
             << timeStep << ".vtk";

    // Open file
    ofile_.open(filename.str());
    if (!ofile_.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename.str());
    }
}

void VTKWallDistanceStencil::closeFile() {
    if (ofile_.is_open()) {
        ofile_.close();
    }
}

void VTKWallDistanceStencil::createDirectory(const std::string& folder) {
#ifdef _MSC_VER
    const int success = _mkdir(folder.c_str());
    if (success != 0 && errno != EEXIST) {
        throw std::runtime_error("Failed to create directory: " + folder);
    }
#else
    const int success = system(("mkdir -p " + folder).c_str());
    if (success != 0) {
        throw std::runtime_error("Failed to create directory: " + folder);
    }
#endif
}

} // namespace Stencils
