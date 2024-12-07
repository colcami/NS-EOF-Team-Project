#pragma once

#include "FieldStencil.hpp"
#include "TurbulentFlowField.hpp"
#include "Parameters.hpp"

namespace Stencils {

  class TurbulentFGHStencil : public FieldStencil<TurbulentFlowField> {
  private:
    RealType localVelocity_[27 * 3];
    RealType localMeshsize_[27 * 3];
  
  public:
    TurbulentFGHStencil(const Parameters& parameters);
    ~TurbulentFGHStencil() override = default;

    void apply(TurbulentFlowField& TflowField, int i, int j) override;
    void apply(TurbulentFlowField& TflowField, int i, int j, int k) override;
  };

} // namespace Stencils