#include "StdAfx.hpp"

#include "PressureBufferFillStencil.hpp"

#include "Definitions.hpp"

template <class FlowFieldType>
Stencils::PressureBufferFillStencil<FlowFieldType>::PressureBufferFillStencil(
  const Parameters&      parameters,
  std::vector<RealType>& leftBuffer,
  std::vector<RealType>& rightBuffer,
  std::vector<RealType>& topBuffer,
  std::vector<RealType>& bottomBuffer,
  std::vector<RealType>& frontBuffer,
  std::vector<RealType>& backBuffer
):
  BoundaryStencil<FlowFieldType>(parameters),
  leftBuffer_(leftBuffer),
  rightBuffer_(rightBuffer),
  topBuffer_(topBuffer),
  bottomBuffer_(bottomBuffer),
  frontBuffer_(frontBuffer),
  backBuffer_(backBuffer) {}

template <class FlowFieldType>
void Stencils::PressureBufferFillStencil<FlowFieldType>::applyLeftWall(FlowFieldType& flowField, int i, int j) {
  leftBuffer_.at(j) = flowField.getPressure().getScalar(i + 2, j);
}

template <class FlowFieldType>
void Stencils::PressureBufferFillStencil<FlowFieldType>::applyRightWall(FlowFieldType& flowField, int i, int j) {
  rightBuffer_.at(j) = flowField.getPressure().getScalar(i - 1, j);
}

template <class FlowFieldType>
void Stencils::PressureBufferFillStencil<FlowFieldType>::applyBottomWall(FlowFieldType& flowField, int i, int j) {
  bottomBuffer_.at(i) = flowField.getPressure().getScalar(i, j + 2);
}

template <class FlowFieldType>
void Stencils::PressureBufferFillStencil<FlowFieldType>::applyTopWall(FlowFieldType& flowField, int i, int j) {
  topBuffer_.at(i) = flowField.getPressure().getScalar(i, j - 1);
}

template <class FlowFieldType>
void Stencils::PressureBufferFillStencil<FlowFieldType>::applyLeftWall(FlowFieldType& flowField, int i, int j, int k) {
  leftBuffer_.at(j * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k) = flowField.getPressure().getScalar(i + 2, j, k);
}

template <class FlowFieldType>
void Stencils::PressureBufferFillStencil<FlowFieldType>::applyRightWall(FlowFieldType& flowField, int i, int j, int k) {
  rightBuffer_.at(j * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k) = flowField.getPressure().getScalar(i - 1, j, k);
}

template <class FlowFieldType>
void Stencils::PressureBufferFillStencil<FlowFieldType>::applyBottomWall(FlowFieldType& flowField, int i, int j, int k) {
  bottomBuffer_.at(i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k) = flowField.getPressure().getScalar(i, j + 2, k);
}

template <class FlowFieldType>
void Stencils::PressureBufferFillStencil<FlowFieldType>::applyTopWall(FlowFieldType& flowField, int i, int j, int k) {
  topBuffer_.at(i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k) = flowField.getPressure().getScalar(i, j - 1, k);
}

template <class FlowFieldType>
void Stencils::PressureBufferFillStencil<FlowFieldType>::applyFrontWall(FlowFieldType& flowField, int i, int j, int k) {
  frontBuffer_.at(i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) + j) = flowField.getPressure().getScalar(i, j, k + 2);
}

template <class FlowFieldType>
void Stencils::PressureBufferFillStencil<FlowFieldType>::applyBackWall(FlowFieldType& flowField, int i, int j, int k) {
  backBuffer_.at(i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) + j) = flowField.getPressure().getScalar(i, j, k - 1);
}
