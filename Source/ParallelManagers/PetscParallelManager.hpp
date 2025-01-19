#pragma once

#include "Parameters.hpp"
#include "Stencils/PressureBufferFillStencil.cpp"
#include "Stencils/PressureBufferReadStencil.cpp"
#include "Stencils/VelocityBufferFillStencil.cpp"
#include "Stencils/VelocityBufferReadStencil.cpp"

namespace ParallelManagers {

  /**
   * PetscParallelManager handles parallel communication of pressure values between processes.
   * It uses buffer stencils to read and write pressure values at domain boundaries and
   * MPI communication to exchange these values between neighboring processes.
   */
  template <class FlowFieldType>
  class PetscParallelManager {
  private:
    Parameters&    parameters_;
    FlowFieldType& flowfield_;

    void initializeBuffers(
      std::vector<RealType>& leftBuffer,
      std::vector<RealType>& rightBuffer,
      std::vector<RealType>& topBuffer,
      std::vector<RealType>& bottomBuffer,
      std::vector<RealType>& frontBuffer,
      std::vector<RealType>& backBuffer,
      int                    componentMultiplier = 1
    ) const;

  public:
    PetscParallelManager(Parameters& parameters, FlowFieldType& flowfield);

    ~PetscParallelManager() = default;

    /** Communicates pressure values between processes:
     * 1. Fills buffer arrays with pressure values using PressureBufferFillStencil
     * 2. Communicates these values between neighboring processes using MPI_Sendrecv
     * 3. Reads received values from buffers and writes them to grid using PressureBufferReadStencil
     */
    void communicatePressure();

    /** Communicates velocity values between processes:
     * 1. Fills buffer arrays with velocity values using VelocityBufferFillStencil
     * 2. Communicates these values between neighboring processes using MPI_Sendrecv
     * 3. Reads received values from buffers and writes them to grid using VelocityBufferReadStencil
     */
    void communicateVelocities();
  };
} // namespace ParallelManagers
