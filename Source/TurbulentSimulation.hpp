#pragma once
#include "StdAfx.hpp"
#include "Definitions.hpp"
#include "FlowField.hpp"
#include "TurbulentFlowField.hpp"
#include "GlobalBoundaryFactory.hpp"
#include "Iterators.hpp"
#include "Simulation.hpp"

#include "Solvers/LinearSolver.hpp"
#include "Stencils/BFInputStencils.hpp"
#include "Stencils/BFStepInitStencil.hpp"
#include "Stencils/FGHStencil.hpp"
#include "Stencils/RHSStencil.hpp"
#include "Stencils/InitTaylorGreenFlowFieldStencil.hpp"
#include "Stencils/MaxUStencil.hpp"
#include "Stencils/MovingWallStencils.hpp"
#include "Stencils/NeumannBoundaryStencils.hpp"
#include "Stencils/ObstacleStencil.hpp"
#include "Stencils/PeriodicBoundaryStencils.hpp"
#include "Stencils/VelocityStencil.hpp"
#include "Stencils/VTKStencil.hpp"
#include "Stencils/TurbulentFGHStencil.hpp"
#include "Stencils/NuTurbulentStencil.hpp"
#include "Stencils/VTKTurbulentViscosityStencil.hpp"
#include "Stencils/VTKWallDistanceStencil.hpp"
#include "Stencils/WallDistanceStencil.hpp"

class TurbulentSimulation : public Simulation {
protected:

  TurbulentFlowField& flowField_;

  //! CHANGE THE FGH STENCIL
  Stencils::TurbulentFGHStencil          TfghStencil_;
  FieldIterator<TurbulentFlowField>      fghIterator_;

  Stencils::WallDistanceStencil      wallDistanceStencil_;
  FieldIterator<TurbulentFlowField>  wallDistanceIterator_;

  Stencils::VTKWallDistanceStencil   vtkWallDistanceStencil_;
  FieldIterator<TurbulentFlowField>  vtkWallDistanceIterator_;

  Stencils::NuTurbulentStencil       NuTurbulentStencil_;
  FieldIterator<TurbulentFlowField>  NuTurbulentIterator_;

  Stencils::VTKTurbulentViscosityStencil vtkTurbulentViscosityStencil_;
  FieldIterator<TurbulentFlowField>       vtkTurbulentViscosityIterator_;
  
  void setTimeStep() override;

public:
  TurbulentSimulation(Parameters& parameters, TurbulentFlowField& TflowField);
  ~TurbulentSimulation() override = default;

  /** Initialises the flow field according to the scenario */
  void initializeFlowField() override;

  void solveTimestep() override;

  /** Plots the flow field */
  void plotVTK(int timeStep, RealType simulationTime) override;

};