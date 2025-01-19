#pragma once

#include <vector>

#include "BoundaryStencil.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /**
   * @brief Stencil for reading velocity values from buffer arrays into domain boundaries
   *
   * Boundary stencil VelocityBufferReadStencil which reads data from onedimensional arrays
   * (one array for each of the six faces) and writes them into the correct cells of
   * the boundary.
   */
  template <class FlowFieldType>
  class VelocityBufferReadStencil: public BoundaryStencil<FlowFieldType> {
  private:
    std::vector<RealType> leftBuffer_;
    std::vector<RealType> rightBuffer_;
    std::vector<RealType> topBuffer_;
    std::vector<RealType> bottomBuffer_;
    std::vector<RealType> frontBuffer_;
    std::vector<RealType> backBuffer_;

  public:
    VelocityBufferReadStencil(
      const Parameters&     parameters,
      std::vector<RealType> leftBuffer,
      std::vector<RealType> rightBuffer,
      std::vector<RealType> topBuffer,
      std::vector<RealType> bottomBuffer,
      std::vector<RealType> frontBuffer,
      std::vector<RealType> backBuffer
    );

    ~VelocityBufferReadStencil() override = default;

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