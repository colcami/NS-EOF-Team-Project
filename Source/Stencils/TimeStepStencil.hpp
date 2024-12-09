#pragma once

#include "FieldStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"

namespace Stencils{
  class TimeStepStencil : public FieldStencil<TurbulentFlowField>, public BoundaryStencil<TurbulentFlowField> {
    private:
      RealType min_dt_;

      void cellMinValue(TurbulentFlowField& flowField, int i, int j);

      void cellMinValue(TurbulentFlowField& flowField, int i, int j, int k);

    public:
      TimeStepStencil(const Parameters& paramters);
      ~TimeStepStencil() override = default;

      void apply(TurbulentFlowField& flowField, int i, int j) override ;
      void apply(TurbulentFlowField& flowField, int i, int j, int k) override ;

      void applyLeftWall(TurbulentFlowField& flowField, int i, int j) ;
      void applyRightWall(TurbulentFlowField& flowField, int i, int j) ;
      void applyBottomWall(TurbulentFlowField& flowField, int i, int j) ;
      void applyTopWall(TurbulentFlowField& flowField, int i, int j) ;

      void applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) ;
      void applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) ;
      void applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) ;
      void applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) ;
      void applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) ;
      void applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) ;

      void reset();

      RealType getdt();
  };
} //namespace Stencils