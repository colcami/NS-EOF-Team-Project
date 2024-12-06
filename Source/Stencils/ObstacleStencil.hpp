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

    
  };

} // namespace Stencils
