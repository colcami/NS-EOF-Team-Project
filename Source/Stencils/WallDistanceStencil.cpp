#include "StdAfx.hpp"
#include "WallDistanceStencil.hpp"
#include "Definitions.hpp"
#include "Parameters.hpp"
#include <cmath>
#include <algorithm>
#include <limits>

namespace Stencils {

// Constructor
WallDistanceStencil::WallDistanceStencil(const Parameters& parameters)
    : FieldStencil<FlowField>(parameters) {
    detectScenario(); // Initialize scenario-specific parameters
}

// Detect the simulation scenario
void WallDistanceStencil::detectScenario() {
    scenario_ = parameters_.simulation.scenario;

    lengthX_ = parameters_.geometry.lengthX;
    lengthY_ = parameters_.geometry.lengthY;
    lengthZ_ = parameters_.geometry.lengthZ;

    stretchX_ = parameters_.geometry.stretchX;
    stretchY_ = parameters_.geometry.stretchY;
    stretchZ_ = parameters_.geometry.stretchZ;

    if (scenario_ == "channel" && parameters_.bfStep.xRatio > 0.0 && parameters_.bfStep.yRatio > 0.0) {
        // Backward-facing step in a channel
        xStep_ = parameters_.bfStep.xRatio * lengthX_;
        yStep_ = parameters_.bfStep.yRatio * lengthY_;
        scenario_ = "channelWithStep";
    } else if (scenario_ != "channel" && scenario_ != "channelWithStep" && scenario_ != "cavity") {
        throw std::runtime_error("Unknown scenario: " + scenario_);
    }
}

// Apply stencil for 2D cases
void WallDistanceStencil::apply(FlowField& flowField, int i, int j) {
    // Get cell center position
    RealType x = parameters_.meshsize->getPosX(i, j);
    RealType y = parameters_.meshsize->getPosY(i, j);

    RealType distance = std::numeric_limits<RealType>::max();

    if (scenario_ == "channel") {
        // Distance to top and bottom walls
        RealType toBottomWall = y;
        RealType toTopWall = lengthY_ - y;
        distance = std::min(toBottomWall, toTopWall);

    } else if (scenario_ == "channelWithStep") {
        if (x <= xStep_ && y <= yStep_) {
            // Inside the step: distance is 0
            distance = 0.0;
        } else if (x <= xStep_ && y > yStep_) {
            // Before the vertical wall of the step
            RealType toBottomWall = y - yStep_;
            RealType toTopWall = lengthY_ - y;
            distance = std::min(std::abs(toBottomWall), std::abs(toTopWall));
        } else if (x > xStep_ && y <= yStep_) {
            // After the step: include the vertical wall of the step
            RealType toBottomWall = y;
            RealType toVerticalWall = x - xStep_; // Horizontal distance to the vertical wall
            RealType toTopWall = lengthY_ - y;
            distance = std::min({std::abs(toBottomWall), std::abs(toVerticalWall), std::abs(toTopWall)});
        } else if (x > xStep_ && y > yStep_ && y < -x + xStep_ + lengthY_) {
            // Diagonal distance to the step corner
            RealType toVerticalWall = x - xStep_;
            RealType toBottomWall = y - yStep_;
            RealType toDiagonalStep = std::sqrt(toBottomWall * toBottomWall + toVerticalWall * toVerticalWall);
            RealType toTopWall = lengthY_ - y;
            distance = std::min(std::abs(toDiagonalStep), std::abs(toTopWall));
        } else {
            // General case: distance to top and bottom walls
            RealType toBottomWall = y;
            RealType toTopWall = lengthY_ - y;
            distance = std::min(std::abs(toBottomWall), std::abs(toTopWall));
        }

    } else if (scenario_ == "cavity") {
        // Distance to cavity walls
        RealType toLeftWall = x;
        RealType toRightWall = lengthX_ - x;
        RealType toBottomWall = y;
        RealType toTopWall = lengthY_ - y;

        distance = std::min({toLeftWall, toRightWall, toBottomWall, toTopWall});
    }

    // Store the result in the scalar field `h`
    flowField.getWallDistance().getScalar(i, j) = distance;
}

// Apply stencil for 3D cases
void WallDistanceStencil::apply(FlowField& flowField, int i, int j, int k) {
    // Get cell center position
    RealType x = parameters_.meshsize->getPosX(i, j, k);
    RealType y = parameters_.meshsize->getPosY(i, j, k);
    RealType z = parameters_.meshsize->getPosZ(i, j, k);

    RealType distance = std::numeric_limits<RealType>::max();

    if (scenario_ == "channel") {
        // Distance to top and bottom walls
        RealType toBottomWall = y;
        RealType toTopWall = lengthY_ - y;

        // Distance to front and back walls
        RealType toFrontWall = z;
        RealType toBackWall = lengthZ_ - z;

        distance = std::min({toBottomWall, toTopWall, toFrontWall, toBackWall});

    } else if (scenario_ == "channelWithStep") {
        if (x <= xStep_ && y <= yStep_) {
            // Inside the step: distance is 0
            distance = 0.0;
        } else if (x <= xStep_ && y > yStep_) {
            // Before the vertical wall of the step
            RealType toBottomWall = y - yStep_;
            RealType toTopWall = lengthY_ - y;
            RealType toFrontWall = z;
            RealType toBackWall = lengthZ_ - z;
            distance = std::min({std::abs(toBottomWall), std::abs(toTopWall), std::abs(toFrontWall), std::abs(toBackWall)});
        } else if (x > xStep_ && y <= yStep_) {
            // After the step: distances to bottom wall, vertical step wall, and others
            RealType toBottomWall = y;
            RealType toVerticalWall = x - xStep_; // Horizontal distance to the vertical wall
            RealType toTopWall = lengthY_ - y;
            RealType toFrontWall = z;
            RealType toBackWall = lengthZ_ - z;

            distance = std::min({std::abs(toBottomWall), std::abs(toVerticalWall), std::abs(toTopWall),
                                 std::abs(toFrontWall), std::abs(toBackWall)});
        } else if (x > xStep_ && y > yStep_ && y < -x + xStep_ + lengthY_) {
            // After the step: diagonal distance to the corner of the step
            RealType toBottomWall = y - yStep_;
            RealType toVerticalWall = x - xStep_;
            RealType toDiagonalStep = std::sqrt(toBottomWall * toBottomWall + toVerticalWall * toVerticalWall);
            RealType toTopWall = lengthY_ - y;
            RealType toFrontWall = z;
            RealType toBackWall = lengthZ_ - z;

            distance = std::min({std::abs(toDiagonalStep), std::abs(toTopWall), std::abs(toFrontWall), std::abs(toBackWall)});
        } else {
            // General case: top and bottom walls
            RealType toBottomWall = y;
            RealType toTopWall = lengthY_ - y;
            RealType toFrontWall = z;
            RealType toBackWall = lengthZ_ - z;

            distance = std::min({std::abs(toBottomWall), std::abs(toTopWall), std::abs(toFrontWall), std::abs(toBackWall)});
        }

    } else if (scenario_ == "cavity") {
        // Distance to cavity walls
        RealType toLeftWall = x;
        RealType toRightWall = lengthX_ - x;
        RealType toBottomWall = y;
        RealType toTopWall = lengthY_ - y;
        RealType toFrontWall = z;
        RealType toBackWall = lengthZ_ - z;

        distance = std::min({toLeftWall, toRightWall, toBottomWall, toTopWall, toFrontWall, toBackWall});
    }

    // Store the result in the scalar field `h`
    flowField.getWallDistance().getScalar(i, j, k) = distance;
}

} // namespace Stencils
