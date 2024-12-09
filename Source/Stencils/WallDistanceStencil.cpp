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
    : FieldStencil<TurbulentFlowField>(parameters) {}

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

void WallDistanceStencil::Distance(TurbulentFlowField& flowField, int i, int j){
    RealType x = parameters_.meshsize->getPosX(i, j);
    RealType y = parameters_.meshsize->getPosY(i, j);

    RealType distance = std::numeric_limits<RealType>::max();

    if (scenario_ == "channel") {
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
        RealType toLeftWall = x;
        RealType toRightWall = lengthX_ - x;
        RealType toBottomWall = y;
        RealType toTopWall = lengthY_ - y;

        distance = std::min({toLeftWall, toRightWall, toBottomWall, toTopWall});
    }

    flowField.getWallDistance().getScalar(i, j) = distance;
}

void WallDistanceStencil::Distance(TurbulentFlowField& flowField, int i, int j, int k){
    RealType x = parameters_.meshsize->getPosX(i, j, k);
    RealType y = parameters_.meshsize->getPosY(i, j, k);
    RealType z = parameters_.meshsize->getPosZ(i, j, k);

    RealType distance = std::numeric_limits<RealType>::max();

    if (scenario_ == "channel") {
        RealType toBottomWall = y;
        RealType toTopWall = lengthY_ - y;
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
            // After the step: include the vertical wall of the step
            RealType toBottomWall = y;
            RealType toVerticalWall = x - xStep_; // Horizontal distance to the vertical wall
            RealType toTopWall = lengthY_ - y;
            RealType toFrontWall = z;
            RealType toBackWall = lengthZ_ - z;
            distance = std::min({std::abs(toBottomWall), std::abs(toVerticalWall), std::abs(toTopWall), std::abs(toFrontWall), std::abs(toBackWall)});
        } else if (x > xStep_ && y > yStep_ && y < -x + xStep_ + lengthY_) {
            // Diagonal distance to the step corner
            RealType toVerticalWall = x - xStep_;
            RealType toBottomWall = y - yStep_;
            RealType toDiagonalStep = std::sqrt(toBottomWall * toBottomWall + toVerticalWall * toVerticalWall);
            RealType toTopWall = lengthY_ - y;
            RealType toFrontWall = z;
            RealType toBackWall = lengthZ_ - z;
            distance = std::min({std::abs(toDiagonalStep), std::abs(toTopWall), std::abs(toFrontWall), std::abs(toBackWall)});
        } else {
            // General case: distance to top and bottom walls
            RealType toBottomWall = y;
            RealType toTopWall = lengthY_ - y;
            RealType toFrontWall = z;
            RealType toBackWall = lengthZ_ - z;
            distance = std::min({std::abs(toBottomWall), std::abs(toTopWall), std::abs(toFrontWall), std::abs(toBackWall)});
        }

    } else if (scenario_ == "cavity") {
        RealType toLeftWall = x;
        RealType toRightWall = lengthX_ - x;
        RealType toBottomWall = y;
        RealType toFrontWall = z;
        RealType toBackWall = lengthZ_ - z;

        distance = std::min({toLeftWall, toRightWall, toBottomWall, toFrontWall, toBackWall});
    }

    flowField.getWallDistance().getScalar(i, j, k) = distance;
}

// Apply stencil for 2D or 3D cases
  void WallDistanceStencil::apply(TurbulentFlowField& flowField, int i, int j) { Distance(flowField, i, j) ;}

  void WallDistanceStencil::apply(TurbulentFlowField& flowField, int i, int j, int k) {Distance(flowField, i, j) ;}

  void WallDistanceStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j)  { Distance(flowField, i, j) ;}
  
  void WallDistanceStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) { Distance(flowField, i, j) ;}
  
  void WallDistanceStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) { Distance(flowField, i, j) ;}
  
  void WallDistanceStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) { Distance(flowField, i, j) ;}
  
  void WallDistanceStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) 
  { Distance(flowField, i, j, k) ;}
  
  void WallDistanceStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) 
  { Distance(flowField, i, j, k) ;}
  
  void WallDistanceStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) 
  { Distance(flowField, i, j, k) ;}
  
  void WallDistanceStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) 
  { Distance(flowField, i, j, k) ;}
  
  void WallDistanceStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) 
  { Distance(flowField, i, j, k) ;}
  
  void WallDistanceStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) 
  { Distance(flowField, i, j, k) ;}  
    
}// namespace Stencils

