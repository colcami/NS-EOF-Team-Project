#include "StdAfx.hpp"

#include "VelocityBufferReadStencil.hpp"

#include "Definitions.hpp"

template <class FlowFieldType>
Stencils::VelocityBufferReadStencil<FlowFieldType>::VelocityBufferReadStencil(
  const Parameters&     parameters,
  std::vector<RealType> leftBuffer,
  std::vector<RealType> rightBuffer,
  std::vector<RealType> topBuffer,
  std::vector<RealType> bottomBuffer,
  std::vector<RealType> frontBuffer,
  std::vector<RealType> backBuffer
):
  BoundaryStencil<FlowFieldType>(parameters),
  leftBuffer_(leftBuffer),
  rightBuffer_(rightBuffer),
  topBuffer_(topBuffer),
  bottomBuffer_(bottomBuffer),
  frontBuffer_(frontBuffer),
  backBuffer_(backBuffer) {}

template <class FlowFieldType>
void Stencils::VelocityBufferReadStencil<FlowFieldType>::applyLeftWall(FlowFieldType& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0]     = leftBuffer_.at(j);
  flowField.getVelocity().getVector(i + 1, j)[1] = leftBuffer_.at(BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3 + j);
}

template <class FlowFieldType>
void Stencils::VelocityBufferReadStencil<FlowFieldType>::applyRightWall(FlowFieldType& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = rightBuffer_.at(j);
  flowField.getVelocity().getVector(i, j)[1] = rightBuffer_.at(BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3 + j);
}

template <class FlowFieldType>
void Stencils::VelocityBufferReadStencil<FlowFieldType>::applyBottomWall(FlowFieldType& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j + 1)[0] = bottomBuffer_.at(i);
  flowField.getVelocity().getVector(i, j)[1]     = bottomBuffer_.at(BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3 + i);
}

template <class FlowFieldType>
void Stencils::VelocityBufferReadStencil<FlowFieldType>::applyTopWall(FlowFieldType& flowField, int i, int j) {
  flowField.getVelocity().getVector(i, j)[0] = topBuffer_.at(i);
  flowField.getVelocity().getVector(i, j)[1] = topBuffer_.at(BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3 + i);
}

template <class FlowFieldType>
void Stencils::VelocityBufferReadStencil<FlowFieldType>::applyLeftWall(FlowFieldType& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0]     = leftBuffer_.at(j * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k);
  flowField.getVelocity().getVector(i + 1, j, k)[1] = leftBuffer_.at(
    j * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k
    + (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3)
  );
  flowField.getVelocity().getVector(i + 1, j, k)[2] = leftBuffer_.at(
    j * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k
    + 2 * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3)
  );
}

template <class FlowFieldType>
void Stencils::VelocityBufferReadStencil<FlowFieldType>::applyRightWall(FlowFieldType& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = rightBuffer_.at(j * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k);
  flowField.getVelocity().getVector(i, j, k)[1] = rightBuffer_.at(
    j * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k
    + (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3)
  );
  flowField.getVelocity().getVector(i, j, k)[2] = rightBuffer_.at(
    j * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k
    + 2 * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3)
  );
}

template <class FlowFieldType>
void Stencils::VelocityBufferReadStencil<FlowFieldType>::applyBottomWall(FlowFieldType& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j + 1, k)[0] = bottomBuffer_.at(i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k);
  flowField.getVelocity().getVector(i, j, k)[1]     = bottomBuffer_.at(
    i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k
    + (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3)
  );
  flowField.getVelocity().getVector(i, j + 1, k)[2] = bottomBuffer_.at(
    i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k
    + 2 * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3)
  );
}

template <class FlowFieldType>
void Stencils::VelocityBufferReadStencil<FlowFieldType>::applyTopWall(FlowFieldType& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = topBuffer_.at(i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k);
  flowField.getVelocity().getVector(i, j, k)[1] = topBuffer_.at(
    i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k
    + (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3)
  );
  flowField.getVelocity().getVector(i, j, k)[2] = topBuffer_.at(
    i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3) + k
    + 2 * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[2] + 3)
  );
}

template <class FlowFieldType>
void Stencils::VelocityBufferReadStencil<FlowFieldType>::applyFrontWall(FlowFieldType& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k + 1)[0] = frontBuffer_.at(i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) + j);
  flowField.getVelocity().getVector(i, j, k + 1)[1] = frontBuffer_.at(
    i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) + j
    + (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3)
  );
  flowField.getVelocity().getVector(i, j, k)[2] = frontBuffer_.at(
    i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) + j
    + 2 * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3)
  );
}

template <class FlowFieldType>
void Stencils::VelocityBufferReadStencil<FlowFieldType>::applyBackWall(FlowFieldType& flowField, int i, int j, int k) {
  flowField.getVelocity().getVector(i, j, k)[0] = backBuffer_.at(i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) + j);
  flowField.getVelocity().getVector(i, j, k)[1] = backBuffer_.at(
    i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) + j
    + (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3)
  );
  flowField.getVelocity().getVector(i, j, k)[2] = backBuffer_.at(
    i * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3) + j
    + 2 * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[0] + 3) * (BoundaryStencil<FlowFieldType>::parameters_.parallel.localSize[1] + 3)
  );
}
