#include "PetscParallelManager.hpp"

#include "Iterators.hpp"
#include "Stencils/PressureBufferFillStencil.hpp"
#include "Stencils/PressureBufferReadStencil.hpp"
#include "Stencils/VelocityBufferFillStencil.hpp"
#include "Stencils/VelocityBufferReadStencil.hpp"

template <class FlowFieldType>
ParallelManagers::PetscParallelManager<FlowFieldType>::PetscParallelManager(Parameters& parameters, FlowFieldType& flowfield):
  parameters_(parameters),
  flowfield_(flowfield) {}

template <class FlowFieldType>
void ParallelManagers::PetscParallelManager<FlowFieldType>::initializeBuffers(
  std::vector<RealType>& leftBuffer,
  std::vector<RealType>& rightBuffer,
  std::vector<RealType>& topBuffer,
  std::vector<RealType>& bottomBuffer,
  std::vector<RealType>& frontBuffer,
  std::vector<RealType>& backBuffer,
  int                    componentMultiplier
) const {

  const int Nx = parameters_.parallel.localSize[0];
  const int Ny = parameters_.parallel.localSize[1];
  const int Nz = parameters_.parallel.localSize[2];

  if (parameters_.geometry.dim == 2) {
    leftBuffer.resize(componentMultiplier * (Ny + 3));
    rightBuffer.resize(componentMultiplier * (Ny + 3));
    topBuffer.resize(componentMultiplier * (Nx + 3));
    bottomBuffer.resize(componentMultiplier * (Nx + 3));
  } else if (parameters_.geometry.dim == 3) {
    leftBuffer.resize(componentMultiplier * ((Ny + 3) * (Nz + 3)));
    rightBuffer.resize(componentMultiplier * ((Ny + 3) * (Nz + 3)));
    topBuffer.resize(componentMultiplier * ((Nx + 3) * (Nz + 3)));
    bottomBuffer.resize(componentMultiplier * ((Nx + 3) * (Nz + 3)));
    frontBuffer.resize(componentMultiplier * ((Nx + 3) * (Ny + 3)));
    backBuffer.resize(componentMultiplier * ((Nx + 3) * (Ny + 3)));
  }
}

template <class FlowFieldType>
void ParallelManagers::PetscParallelManager<FlowFieldType>::communicatePressure() {

  std::vector<RealType> leftPressureBuffer, rightPressureBuffer, topPressureBuffer, bottomPressureBuffer, frontPressureBuffer, backPressureBuffer;
  std::vector<RealType> leftReceive, rightReceive, topReceive, bottomReceive, frontReceive, backReceive;

  initializeBuffers(leftPressureBuffer, rightPressureBuffer, topPressureBuffer, bottomPressureBuffer, frontPressureBuffer, backPressureBuffer, 1);
  initializeBuffers(leftReceive, rightReceive, topReceive, bottomReceive, frontReceive, backReceive, 1);

  Stencils::PressureBufferFillStencil<FlowFieldType> pressureBufferFillStencil_(
    parameters_, leftPressureBuffer, rightPressureBuffer, topPressureBuffer, bottomPressureBuffer, frontPressureBuffer, backPressureBuffer
  );
  ParallelBoundaryIterator<FlowFieldType> pressureFillIterator(flowfield_, parameters_, pressureBufferFillStencil_, 0, 0);
  pressureFillIterator.iterate();

  MPI_Status statusPressure;


  MPI_Sendrecv(
    leftPressureBuffer.data(),
    leftPressureBuffer.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.leftNb,
    66,
    rightReceive.data(),
    rightReceive.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.rightNb,
    66,
    PETSC_COMM_WORLD,
    &statusPressure
  );

  MPI_Sendrecv(
    rightPressureBuffer.data(),
    rightPressureBuffer.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.rightNb,
    47,
    leftReceive.data(),
    leftReceive.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.leftNb,
    47,
    PETSC_COMM_WORLD,
    &statusPressure
  );

  MPI_Sendrecv(
    topPressureBuffer.data(),
    topPressureBuffer.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.topNb,
    31,
    bottomReceive.data(),
    bottomReceive.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.bottomNb,
    31,
    PETSC_COMM_WORLD,
    &statusPressure
  );

  MPI_Sendrecv(
    bottomPressureBuffer.data(),
    bottomPressureBuffer.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.bottomNb,
    34,
    topReceive.data(),
    topReceive.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.topNb,
    34,
    PETSC_COMM_WORLD,
    &statusPressure
  );

  MPI_Sendrecv(
    frontPressureBuffer.data(),
    frontPressureBuffer.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.frontNb,
    73,
    backReceive.data(),
    backReceive.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.backNb,
    73,
    PETSC_COMM_WORLD,
    &statusPressure
  );

  MPI_Sendrecv(
    backPressureBuffer.data(),
    backPressureBuffer.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.backNb,
    41,
    frontReceive.data(),
    frontReceive.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.frontNb,
    41,
    PETSC_COMM_WORLD,
    &statusPressure
  );

  Stencils::PressureBufferReadStencil<FlowFieldType> pressureBufferReadStencil_(parameters_, leftReceive, rightReceive, topReceive, bottomReceive, frontReceive, backReceive);
  ParallelBoundaryIterator<FlowFieldType>            pressureReadIterator(flowfield_, parameters_, pressureBufferReadStencil_, 0, 0);
  pressureReadIterator.iterate();
}

template <class FlowFieldType>
void ParallelManagers::PetscParallelManager<FlowFieldType>::communicateVelocities() {

  std::vector<RealType> leftVelocityBuffer, rightVelocityBuffer, topVelocityBuffer, bottomVelocityBuffer, frontVelocityBuffer, backVelocityBuffer;
  std::vector<RealType> leftReceive, rightReceive, topReceive, bottomReceive, frontReceive, backReceive;

  const int velocityComponents = (parameters_.geometry.dim == 2) ? 2 : 3;
  initializeBuffers(leftVelocityBuffer, rightVelocityBuffer, topVelocityBuffer, bottomVelocityBuffer, frontVelocityBuffer, backVelocityBuffer, velocityComponents);
  initializeBuffers(leftReceive, rightReceive, topReceive, bottomReceive, frontReceive, backReceive, velocityComponents);

  Stencils::VelocityBufferFillStencil<FlowFieldType> velocityBufferFillStencil_(
    parameters_, leftVelocityBuffer, rightVelocityBuffer, topVelocityBuffer, bottomVelocityBuffer, frontVelocityBuffer, backVelocityBuffer
  );
  ParallelBoundaryIterator<FlowFieldType> velocityFillIterator(flowfield_, parameters_, velocityBufferFillStencil_, 0, 0);
  velocityFillIterator.iterate();

  MPI_Status statusVelocity;

  MPI_Sendrecv(
    leftVelocityBuffer.data(),
    leftVelocityBuffer.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.leftNb,
    1910,
    rightReceive.data(),
    rightReceive.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.rightNb,
    1910,
    PETSC_COMM_WORLD,
    &statusVelocity
  );

  MPI_Sendrecv(
    rightVelocityBuffer.data(),
    rightVelocityBuffer.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.rightNb,
    1915,
    leftReceive.data(),
    leftReceive.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.leftNb,
    1915,
    PETSC_COMM_WORLD,
    &statusVelocity
  );

  MPI_Sendrecv(
    topVelocityBuffer.data(),
    topVelocityBuffer.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.topNb,
    1920,
    bottomReceive.data(),
    bottomReceive.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.bottomNb,
    1920,
    PETSC_COMM_WORLD,
    &statusVelocity
  );

  MPI_Sendrecv(
    bottomVelocityBuffer.data(),
    bottomVelocityBuffer.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.bottomNb,
    1922,
    topReceive.data(),
    topReceive.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.topNb,
    1922,
    PETSC_COMM_WORLD,
    &statusVelocity
  );

  MPI_Sendrecv(
    frontVelocityBuffer.data(),
    frontVelocityBuffer.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.frontNb,
    1923,
    backReceive.data(),
    backReceive.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.backNb,
    1923,
    PETSC_COMM_WORLD,
    &statusVelocity
  );

  MPI_Sendrecv(
    backVelocityBuffer.data(),
    backVelocityBuffer.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.backNb,
    1938,
    frontReceive.data(),
    frontReceive.size(),
    MY_MPI_FLOAT,
    parameters_.parallel.frontNb,
    1938,
    PETSC_COMM_WORLD,
    &statusVelocity
  );

  Stencils::VelocityBufferReadStencil<FlowFieldType> velocityBufferReadStencil_(parameters_, leftReceive, rightReceive, topReceive, bottomReceive, frontReceive, backReceive);
  ParallelBoundaryIterator<FlowFieldType>            velocityReadIterator(flowfield_, parameters_, velocityBufferReadStencil_, 0, 0);
  velocityReadIterator.iterate();
}
