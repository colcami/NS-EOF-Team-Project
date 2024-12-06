/* #pragma once

#include "Definitions.hpp"
#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include <sstream>
#include <fstream>

namespace Stencils {
class VTKTurbulentViscosityStencil : public FieldStencil<FlowField> {
private:
    bool          written_;
    std::string   prefix_;
    std::ofstream ofile_;

    std::stringstream viscosityStream_;

    void writeVTKHeader(std::ostream& file) const;
    void writePoints(std::ostream& file, RealType simulationTime) const;

    void openFile(int timeStep, RealType simulationTime);
    void closeFile();

public:
    VTKTurbulentViscosityStencil(const Parameters& parameters);
    ~VTKTurbulentViscosityStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;
    void apply(FlowField& flowField, int i, int j, int k) override;

    void write(FlowField& flowField, int timeStep, RealType simulationTime);
};
}
 */