#include "StdAfx.hpp"
#include "Definitions.hpp"
#include "StencilFunctions.hpp"

#include "TimeStepStencil.hpp"

namespace Stencils{
  TimeStepStencil::TimeStepStencil(const Parameters& parameters):
    FieldStencil<TurbulentFlowField>(parameters),
    min_dt_(MY_FLOAT_MAX){}

  void TimeStepStencil::cellMinValue(TurbulentFlowField& flowField, int i, int j) {
    RealType       nu_t =      flowField.getTurbulentViscosity().getScalar(i, j);
    RealType factor = 1.0 / (parameters_.meshsize->getDxMin() * parameters_.meshsize->getDxMin())
                    + 1.0/ (parameters_.meshsize->getDyMin() * parameters_.meshsize->getDyMin());

    RealType temporary = 1 / parameters_.flow.Re + nu_t ;
    min_dt_ = std::min(1 /(2 * factor * temporary), min_dt_);
    }
  
  void TimeStepStencil::cellMinValue(TurbulentFlowField& flowField, int i, int j, int k) {
    RealType       nu_t =      flowField.getTurbulentViscosity().getScalar(i, j, k);
    RealType factor = 1.0 / (parameters_.meshsize->getDxMin() * parameters_.meshsize->getDxMin())
                    + 1.0/ (parameters_.meshsize->getDyMin() * parameters_.meshsize->getDyMin())
                    + 1.0/ (parameters_.meshsize->getDzMin() * parameters_.meshsize->getDzMin());

    RealType temporary = 1 / parameters_.flow.Re + nu_t;
    min_dt_ = std::min(1 /(2 * factor * temporary), min_dt_);
  }

  void TimeStepStencil::apply(TurbulentFlowField& flowField, int i, int j){ cellMinValue(flowField, i, j); }

  void TimeStepStencil::apply(TurbulentFlowField& flowField, int i, int j, int k){ cellMinValue(flowField, i, j, k); }

  void TimeStepStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j) { cellMinValue(flowField, i, j) ;} 
  
  void TimeStepStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j) { cellMinValue(flowField, i, j) ;} 
  
  void TimeStepStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j) { cellMinValue(flowField, i, j) ;} 
  
  void TimeStepStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j) { cellMinValue(flowField, i, j) ;} 

  void TimeStepStencil::applyLeftWall(TurbulentFlowField& flowField, int i, int j, int k) 
  { cellMinValue(flowField, i, j, k) ;} 
  
  void TimeStepStencil::applyRightWall(TurbulentFlowField& flowField, int i, int j, int k) 
  { cellMinValue(flowField, i, j, k) ;} 
  
  void TimeStepStencil::applyBottomWall(TurbulentFlowField& flowField, int i, int j, int k) 
  { cellMinValue(flowField, i, j, k) ;} 
  
  void TimeStepStencil::applyTopWall(TurbulentFlowField& flowField, int i, int j, int k) 
  { cellMinValue(flowField, i, j, k) ;} 
  
  void TimeStepStencil::applyFrontWall(TurbulentFlowField& flowField, int i, int j, int k) 
  { cellMinValue(flowField, i, j, k) ;} 
  
  void TimeStepStencil::applyBackWall(TurbulentFlowField& flowField, int i, int j, int k) 
  { cellMinValue(flowField, i, j, k) ;} 

  void TimeStepStencil::reset(){
    const double eps = 10e-16;
    min_dt_ = eps;
  }

  RealType TimeStepStencil::getdt(){return min_dt_;}

}