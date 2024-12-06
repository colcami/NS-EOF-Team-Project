#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"
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
class WallDistanceStencil : public FieldStencil<FlowField> {
private:
    std::string scenario_; //! Simulation scenario: channel, cavity, or backward-facing step
    RealType lengthX_, lengthY_, lengthZ_; //! Domain dimensions
    RealType xStep_, yStep_; //! Backward-facing step parameters
    bool stretchX_, stretchY_, stretchZ_; //! Flags for stretched grids

    void detectScenario(); //! Detect the scenario and initialize relevant parameters

public:
    WallDistanceStencil(const Parameters& parameters);
    ~WallDistanceStencil() override = default;

    void apply(FlowField& flowField, int i, int j) override;       // For 2D
    void apply(FlowField& flowField, int i, int j, int k) override; // For 3D
};
