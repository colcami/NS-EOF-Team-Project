#include "StdAfx.hpp"

#include "PressureBufferReadStencil.hpp"

#include "Definitions.hpp"

template <class FlowFieldType>
Stencils::PressureBufferReadStencil<FlowFieldType>::PressureBufferReadStencil(
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
void Stencils::PressureBufferReadStencil<FlowFieldType>::applyLeftWall(FlowFieldType& flowField, int i, int j) {
  flowField.getPressure().getScalar(i + 1, j) = leftBuffer_.at(j);
}

template <class FlowFieldType>
void Stencils::PressureBufferReadStencil<FlowFieldType>::applyRightWall(FlowFieldType& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = rightBuffer_.at(j);
}

template <class FlowFieldType>
void Stencils::PressureBufferReadStencil<FlowFieldType>::applyBottomWall(FlowFieldType& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j + 1) = bottomBuffer_.at(i);
}

template <class FlowFieldType>
void Stencils::PressureBufferReadStencil<FlowFieldType>::applyTopWall(FlowFieldType& flowField, int i, int j) {
  flowField.getPressure().getScalar(i, j) = topBuffer_.at(i);
}

template <class FlowFieldType>
void Stencils::PressureBufferReadStencil<FlowFieldType>::applyLeftWall(FlowFieldType& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i + 1, j, k) = leftBuffer_.at(j * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k);
}

template <class FlowFieldType>
void Stencils::PressureBufferReadStencil<FlowFieldType>::applyRightWall(FlowFieldType& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = rightBuffer_.at(j * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k);
}

template <class FlowFieldType>
void Stencils::PressureBufferReadStencil<FlowFieldType>::applyBottomWall(FlowFieldType& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j + 1, k) = bottomBuffer_.at(i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k);
}

template <class FlowFieldType>
void Stencils::PressureBufferReadStencil<FlowFieldType>::applyTopWall(FlowFieldType& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = topBuffer_.at(i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k);
}

template <class FlowFieldType>
void Stencils::PressureBufferReadStencil<FlowFieldType>::applyFrontWall(FlowFieldType& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k + 1) = frontBuffer_.at(i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) + j);
}

template <class FlowFieldType>
void Stencils::PressureBufferReadStencil<FlowFieldType>::applyBackWall(FlowFieldType& flowField, int i, int j, int k) {
  flowField.getPressure().getScalar(i, j, k) = backBuffer_.at(i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) + j);
}
