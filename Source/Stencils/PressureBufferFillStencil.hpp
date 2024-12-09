#pragma once

#include <cassert>
#include <vector>

#include "BoundaryStencil.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /**
   * @brief Stencil for filling pressure buffer arrays from domain boundaries
   *
   * Boundary stencil PressureBufferFillStencil which reads the pressure values
   * in each of the six (3D) boundary faces of a sub-domain (i.e. the domain of one process) and
   * stores them consecutively in one-dimensional buffer arrays (one array for each of the six faces).
   */
  template <class FlowFieldType>
  class PressureBufferFillStencil: public BoundaryStencil<FlowFieldType> {
  private:
    std::vector<RealType>& leftBuffer_;
    std::vector<RealType>& rightBuffer_;
    std::vector<RealType>& topBuffer_;
    std::vector<RealType>& bottomBuffer_;
    std::vector<RealType>& frontBuffer_;
    std::vector<RealType>& backBuffer_;

  public:
    PressureBufferFillStencil(
      const Parameters&      parameters,
      std::vector<RealType>& leftBuffer,
      std::vector<RealType>& rightBuffer,
      std::vector<RealType>& topBuffer,
      std::vector<RealType>& bottomBuffer,
      std::vector<RealType>& frontBuffer,
      std::vector<RealType>& backBuffer
    );
    ~PressureBufferFillStencil() override = default;

    void applyLeftWall(FlowFieldType& flowField, int i, int j) override;
    void applyRightWall(FlowFieldType& flowField, int i, int j) override;
    void applyBottomWall(FlowFieldType& flowField, int i, int j) override;
    void applyTopWall(FlowFieldType& flowField, int i, int j) override;

    void applyLeftWall(FlowFieldType& flowField, int i, int j, int k) override;
    void applyRightWall(FlowFieldType& flowField, int i, int j, int k) override;
    void applyBottomWall(FlowFieldType& flowField, int i, int j, int k) override;
    void applyTopWall(FlowFieldType& flowField, int i, int j, int k) override;
    void applyFrontWall(FlowFieldType& flowField, int i, int j, int k) override;
    void applyBackWall(FlowFieldType& flowField, int i, int j, int k) override;
  };

} // namespace Stencils