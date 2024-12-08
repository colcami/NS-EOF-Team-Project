#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
#include "TurbulentFlowField.hpp"
#include <cmath>
#include <algorithm>
#include <string>
#include <limits>
#include "Definitions.hpp"

/** Stencil for calculating normal distances to the nearest wall
 *
 * This stencil computes the normal distance `h` from the center of each cell to the nearest wall
 * for scenarios like channel, cavity, and backward-facing step.
 */
namespace Stencils {
    
class WallDistanceStencil : public FieldStencil<TurbulentFlowField> {
private:
    std::string scenario_; //! Simulation scenario: channel, cavity, or backward-facing step
    RealType lengthX_, lengthY_, lengthZ_; //! Domain dimensions
    RealType xStep_, yStep_; //! Backward-facing step parameters
    bool stretchX_, stretchY_, stretchZ_; //! Flags for stretched grids
    int dimensions_; //! Number of dimensions (2 or 3)

    void detectScenario(); //! Detect the scenario and initialize relevant parameters

public:
    WallDistanceStencil(const Parameters& parameters);
    ~WallDistanceStencil() override = default;

    void apply(TurbulentFlowField& flowField, int i, int j) override;       // Dispatch to apply2D or apply3D
    void apply(TurbulentFlowField& flowField, int i, int j, int k) override; // Dispatch to apply3D
};
} // namespace Stencils
