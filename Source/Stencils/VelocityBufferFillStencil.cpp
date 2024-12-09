#include "StdAfx.hpp"

#include "VelocityBufferFillStencil.hpp"

#include "Definitions.hpp"

template <class FlowFieldType>
Stencils::VelocityBufferFillStencil<FlowFieldType>::VelocityBufferFillStencil(
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
void Stencils::VelocityBufferFillStencil<FlowFieldType>::applyLeftWall(FlowFieldType& flowField, int i, int j) {
  leftBuffer_.at(j)                                                                         = flowField.getVelocity().getVector(i + 2, j)[0];
  leftBuffer_.at(BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3 + j) = flowField.getVelocity().getVector(i + 2, j)[1];
}

template <class FlowFieldType>
void Stencils::VelocityBufferFillStencil<FlowFieldType>::applyRightWall(FlowFieldType& flowField, int i, int j) {
  rightBuffer_.at(j)                                                                         = flowField.getVelocity().getVector(i - 2, j)[0];
  rightBuffer_.at(BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3 + j) = flowField.getVelocity().getVector(i - 1, j)[1];
}

template <class FlowFieldType>
void Stencils::VelocityBufferFillStencil<FlowFieldType>::applyBottomWall(FlowFieldType& flowField, int i, int j) {
  bottomBuffer_.at(i)                                                                         = flowField.getVelocity().getVector(i, j + 2)[0];
  bottomBuffer_.at(BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3 + i) = flowField.getVelocity().getVector(i, j + 2)[1];
}

template <class FlowFieldType>
void Stencils::VelocityBufferFillStencil<FlowFieldType>::applyTopWall(FlowFieldType& flowField, int i, int j) {
  topBuffer_.at(i)                                                                         = flowField.getVelocity().getVector(i, j - 1)[0];
  topBuffer_.at(BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3 + i) = flowField.getVelocity().getVector(i, j - 2)[1];
}

template <class FlowFieldType>
void Stencils::VelocityBufferFillStencil<FlowFieldType>::applyLeftWall(FlowFieldType& flowField, int i, int j, int k) {
  leftBuffer_.at(j * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k) = flowField.getVelocity().getVector(i + 2, j, k)[0];
  leftBuffer_.at(
    j * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k
    + (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3)
  ) = flowField.getVelocity().getVector(i + 2, j, k)[1];
  leftBuffer_.at(
    j * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k
    + 2 * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3)
  ) = flowField.getVelocity().getVector(i + 2, j, k)[2];
}

template <class FlowFieldType>
void Stencils::VelocityBufferFillStencil<FlowFieldType>::applyRightWall(FlowFieldType& flowField, int i, int j, int k) {
  rightBuffer_.at(j * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k) = flowField.getVelocity().getVector(i - 2, j, k)[0];
  rightBuffer_.at(
    j * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k
    + (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3)
  ) = flowField.getVelocity().getVector(i - 1, j, k)[1];
  rightBuffer_.at(
    j * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k
    + 2 * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3)
  ) = flowField.getVelocity().getVector(i - 1, j, k)[2];
}

template <class FlowFieldType>
void Stencils::VelocityBufferFillStencil<FlowFieldType>::applyBottomWall(FlowFieldType& flowField, int i, int j, int k) {
  bottomBuffer_.at(i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k) = flowField.getVelocity().getVector(i, j + 2, k)[0];
  bottomBuffer_.at(
    i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k
    + (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3)
  ) = flowField.getVelocity().getVector(i, j + 2, k)[1];
  bottomBuffer_.at(
    i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k
    + 2 * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3)
  ) = flowField.getVelocity().getVector(i, j + 2, k)[2];
}

template <class FlowFieldType>
void Stencils::VelocityBufferFillStencil<FlowFieldType>::applyTopWall(FlowFieldType& flowField, int i, int j, int k) {
  topBuffer_.at(i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k) = flowField.getVelocity().getVector(i, j - 1, k)[0];
  topBuffer_.at(
    i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k
    + (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3)
  ) = flowField.getVelocity().getVector(i, j - 2, k)[1];
  topBuffer_.at(
    i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k
    + 2 * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3)
  ) = flowField.getVelocity().getVector(i, j - 1, k)[2];
}

template <class FlowFieldType>
void Stencils::VelocityBufferFillStencil<FlowFieldType>::applyFrontWall(FlowFieldType& flowField, int i, int j, int k) {
  frontBuffer_.at(i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) + j) = flowField.getVelocity().getVector(i, j, k + 2)[0];
  frontBuffer_.at(
    i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) + j
    + (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3)
  ) = flowField.getVelocity().getVector(i, j, k + 2)[1];
  frontBuffer_.at(
    i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) + j
    + 2 * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3)
  ) = flowField.getVelocity().getVector(i, j, k + 2)[2];
}

template <class FlowFieldType>
void Stencils::VelocityBufferFillStencil<FlowFieldType>::applyBackWall(FlowFieldType& flowField, int i, int j, int k) {
  backBuffer_.at(i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) + j) = flowField.getVelocity().getVector(i, j, k - 1)[0];
  backBuffer_.at(
    i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) + j
    + (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3)
  ) = flowField.getVelocity().getVector(i, j, k - 1)[1];
  backBuffer_.at(
    i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) + j
    + 2 * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3)
  ) = flowField.getVelocity().getVector(i, j, k - 2)[2];
}
