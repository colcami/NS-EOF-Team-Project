#include "StdAfx.hpp"
#include "WallDistanceStencil.hpp"
#include "Definitions.hpp"
#include "Parameters.hpp"
#include <cmath>
#include <algorithm>


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

    if (scenario_ == "backward facing step") {
        xStep_ = parameters_.backwardFacingStep.xRatio * lengthX_;
        yStep_ = parameters_.backwardFacingStep.yRatio * lengthY_;
    } else if (scenario_ != "channel" && scenario_ != "cavity") {
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

    } else if (scenario_ == "cavity") {
        // Distance to cavity walls
        RealType toLeftWall = x;
        RealType toRightWall = lengthX_ - x;
        RealType toBottomWall = y;
        RealType toTopWall = lengthY_ - y;

        distance = std::min({toLeftWall, toRightWall, toBottomWall, toTopWall});

    } else if (scenario_ == "backward facing step") {
        if (x < xStep_) {
            // Before the step: distance to flat bottom wall
            distance = y;
        } else {
            // After the step: compute distances to the raised bottom wall and step edge
            RealType toBottomWall = std::abs(y - yStep_);
            RealType toStepEdge = std::sqrt(std::pow(x - xStep_, 2) + std::pow(y - yStep_, 2));
            distance = std::min(toBottomWall, toStepEdge);
        }

        // Also consider the top wall
        RealType toTopWall = lengthY_ - y;
        distance = std::min(distance, toTopWall);
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

    } else if (scenario_ == "cavity") {
        // Distance to cavity walls
        RealType toLeftWall = x;
        RealType toRightWall = lengthX_ - x;
        RealType toBottomWall = y;
        RealType toTopWall = lengthY_ - y;
        RealType toFrontWall = z;
        RealType toBackWall = lengthZ_ - z;

        distance = std::min({toLeftWall, toRightWall, toBottomWall, toTopWall, toFrontWall, toBackWall});

    } else if (scenario_ == "backward facing step") {
        if (x < xStep_) {
            // Before the step: distance to flat bottom wall
            distance = y;
        } else {
            // After the step: compute distances to the raised bottom wall and step edge
            RealType toBottomWall = std::abs(y - yStep_);
            RealType toStepEdge = std::sqrt(
                std::pow(x - xStep_, 2) + std::pow(y - yStep_, 2) + std::pow(z - lengthZ_ / 2, 2)
            );
            distance = std::min(toBottomWall, toStepEdge);
        }

        // Also consider distances to top, front, and back walls
        RealType toTopWall = lengthY_ - y;
        RealType toFrontWall = z;
        RealType toBackWall = lengthZ_ - z;

        distance = std::min({distance, toTopWall, toFrontWall, toBackWall});
    }

    // Store the result in the scalar field `h`
    flowField.getWallDistance().getScalar(i, j, k) = distance;
}
