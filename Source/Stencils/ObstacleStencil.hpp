#pragma once

#include "FieldStencil.hpp"
#include "FlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  /** Compute all velocities on obstacle cells
   */
  class ObstacleStencil: public FieldStencil<FlowField> {
  public:
    ObstacleStencil(const Parameters& parameters);
    ~ObstacleStencil() override = default;

    // Apply existing no-slip conditions and compute minimal distances
    void apply(FlowField& flowField, int i, int j) override;       // For 2D
    void apply(FlowField& flowField, int i, int j, int k) override; // For 3D

    // Access minimal distance
    RealType getMinimalDistance(int i, int j) const;
    RealType getMinimalDistance(int i, int j, int k) const;
 
  private:
    // Minimal distances storage
    std::vector<std::vector<RealType>> minimalDistances2D;
    std::vector<std::vector<std::vector<RealType>>> minimalDistances3D;

    void initializeDistances2D(int sizeX, int sizeY);
    void initializeDistances3D(int sizeX, int sizeY, int sizeZ);
  };

} // namespace Stencils
