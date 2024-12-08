#pragma once

#include <vector>

#include "BoundaryStencil.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /**
   * @brief Stencil for filling velocity buffer arrays from domain boundaries
   *
   * Boundary stencil VelocityBufferFillStencil which reads the velocity values
   * in each of the six (3D) boundary faces of a sub-domain (i.e. the domain of one process) and
   * stores them consecutively in one-dimensional buffer arrays (one array for each of the six faces).
   */
  template <class FlowFieldType>
  class VelocityBufferFillStencil: public BoundaryStencil<FlowFieldType> {
  private:
    std::vector<RealType>& leftBuffer_;
    std::vector<RealType>& rightBuffer_;
    std::vector<RealType>& topBuffer_;
    std::vector<RealType>& bottomBuffer_;
    std::vector<RealType>& frontBuffer_;
    std::vector<RealType>& backBuffer_;

  public:
    VelocityBufferFillStencil(
      const Parameters&      parameters,
      std::vector<RealType>& leftBuffer,
      std::vector<RealType>& rightBuffer,
      std::vector<RealType>& topBuffer,
      std::vector<RealType>& bottomBuffer,
      std::vector<RealType>& frontBuffer,
      std::vector<RealType>& backBuffer
    );
    ~VelocityBufferFillStencil() override = default;

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