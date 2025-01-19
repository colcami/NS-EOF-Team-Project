#pragma once

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "WallDistanceStencil.hpp"
#include "TurbulentFlowField.hpp"
#include <sstream>
#include <fstream>
#include <string>

namespace Stencils {

/** Class to export wall distance data to VTK format.
 *
 * This stencil handles the computation and export of wall distance data stored
 * in the flow field to VTK files for visualization.
 */
class VTKWallDistanceStencil : public FieldStencil<TurbulentFlowField> {
private:
    bool written_;               //! Indicates whether the file has been written
    std::string prefix_;         //! Prefix for the output VTK files
    std::ofstream ofile_;        //! Output file stream

    std::stringstream wallDistanceStream_; //! Stream to store wall distance data

    /** Writes the VTK file header */
    void writeVTKHeader(std::ostream& file) const;

    /** Writes the grid points to the VTK file */
    void writePoints(std::ostream& file, RealType simulationTime) const;

    /** Opens the output VTK file for writing */
    void openFile(int timeStep, RealType simulationTime);

    /** Closes the output VTK file */
    void closeFile();

    /** Ensures the output directory exists */
    void createDirectory(const std::string& folder);

public:
    /** Constructor */
    VTKWallDistanceStencil(const Parameters& parameters);

    /** Destructor */
    ~VTKWallDistanceStencil() override = default;

    /** Applies the stencil for a 2D grid */
    void apply(TurbulentFlowField& flowField, int i, int j) override;

    /** Applies the stencil for a 3D grid */
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override;

    /** Writes the wall distance data to a VTK file */
    void write(TurbulentFlowField& flowField, int timeStep, RealType simulationTime);
};

} // namespace Stencils
