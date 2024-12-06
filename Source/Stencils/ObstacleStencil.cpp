#include "StdAfx.hpp"

#include "ObstacleStencil.hpp"
#include <cmath>
#include <limits>

namespace Stencils {
ObstacleStencil::ObstacleStencil(const Parameters& parameters):
  FieldStencil<FlowField>(parameters) {}

// Initialize minimal distances for 2D
void ObstacleStencil::initializeDistances2D(int sizeX, int sizeY) {
    minimalDistances2D.resize(sizeX, std::vector<RealType>(sizeY, std::numeric_limits<RealType>::max()));
}

// Initialize minimal distances for 3D
void ObstacleStencil::initializeDistances3D(int sizeX, int sizeY, int sizeZ) {
    minimalDistances3D.resize(sizeX, std::vector<std::vector<RealType>>(sizeY, std::vector<RealType>(sizeZ, std::numeric_limits<RealType>::max())));
}

void Stencils::ObstacleStencil::apply(FlowField& flowField, int i, int j) {
  const int    obstacle = flowField.getFlags().getValue(i, j);
  VectorField& velocity = flowField.getVelocity();

  // Check if current cell is obstacle cell
  if ((obstacle & OBSTACLE_SELF) == 1) {
    // If top cell is fluid, then the no-slip boundary has to be enforced
    if ((obstacle & OBSTACLE_TOP) == 0) {
      const RealType dy_t         = parameters_.meshsize->getDy(i, j + 1);
      const RealType dy           = parameters_.meshsize->getDy(i, j);
      velocity.getVector(i, j)[0] = -dy / dy_t * velocity.getVector(i, j + 1)[0];
    }
    // Same for bottom
    if ((obstacle & OBSTACLE_BOTTOM) == 0) {
      const RealType dy_b         = parameters_.meshsize->getDy(i, j - 1);
      const RealType dy           = parameters_.meshsize->getDy(i, j);
      velocity.getVector(i, j)[0] = -dy / dy_b * velocity.getVector(i, j - 1)[0];
    }
    // If right cell is fluid, then the no-slip boundary has to be enforced
    if ((obstacle & OBSTACLE_RIGHT) == 0) {
      const RealType dx_r         = parameters_.meshsize->getDx(i + 1, j);
      const RealType dx           = parameters_.meshsize->getDx(i, j);
      velocity.getVector(i, j)[1] = -dx / dx_r * velocity.getVector(i + 1, j)[1];
    }
    // Same for left
    if ((obstacle & OBSTACLE_LEFT) == 0) {
      const RealType dx_l         = parameters_.meshsize->getDx(i - 1, j);
      const RealType dx           = parameters_.meshsize->getDx(i, j);
      velocity.getVector(i, j)[1] = -dx / dx_l * velocity.getVector(i - 1, j)[1];
    }

    // Set normal velocity to zero if right neighbour is not obstacle
    if ((obstacle & OBSTACLE_RIGHT) == 0) {
      velocity.getVector(i, j)[0] = 0.0;
    }

    // Set normal velocity to zero if top neighbour is not obstacle
    if ((obstacle & OBSTACLE_TOP) == 0) {
      velocity.getVector(i, j)[1] = 0.0;
    }
  }
// --- New Minimal Distance Calculation ---
    if (minimalDistances2D.empty()) {
        initializeDistances2D(flowField.getNx(), flowField.getNy());
    }

    if ((obstacle & OBSTACLE_SELF) == 1) {
        minimalDistances2D[i][j] = 0.0; // Obstacle cells have zero distance
    } else {
        for (int di = -1; di <= 1; ++di) {
            for (int dj = -1; dj <= 1; ++dj) {
                int ni = i + di;
                int nj = j + dj;

                // Check bounds
                if (ni >= 0 && ni < flowField.getNx() && nj >= 0 && nj < flowField.getNy()) {
                    if ((flowField.getFlags().getValue(ni, nj) & OBSTACLE_SELF) == 1) {
                        const RealType dx = parameters_.meshsize->getDx(i, j);
                        const RealType dy = parameters_.meshsize->getDy(i, j);
                        RealType dist = std::sqrt(di * di * dx * dx + dj * dj * dy * dy);
                        minimalDistances2D[i][j] = std::min(minimalDistances2D[i][j], dist);
                    }
                }
            }
        }
    }
}


void Stencils::ObstacleStencil::apply(FlowField& flowField, int i, int j, int k) {
  const int    obstacle = flowField.getFlags().getValue(i, j);
  VectorField& velocity = flowField.getVelocity();

  // Check if current cell is obstacle cell
  if ((obstacle & OBSTACLE_SELF) == 1) {
    // If top cell is fluid: two velocities have to be set: direction 0 and 2.
    if ((obstacle & OBSTACLE_TOP) == 0) {
      const RealType dy_t            = parameters_.meshsize->getDy(i, j + 1, k);
      const RealType dy              = parameters_.meshsize->getDy(i, j, k);
      velocity.getVector(i, j, k)[0] = -dy / dy_t * velocity.getVector(i, j + 1, k)[0];
      velocity.getVector(i, j, k)[2] = -dy / dy_t * velocity.getVector(i, j + 1, k)[2];
    }
    if ((obstacle & OBSTACLE_BOTTOM) == 0) {
      const RealType dy_b            = parameters_.meshsize->getDy(i, j - 1, k);
      const RealType dy              = parameters_.meshsize->getDy(i, j, k);
      velocity.getVector(i, j, k)[0] = -dy / dy_b * velocity.getVector(i, j - 1, k)[0];
      velocity.getVector(i, j, k)[2] = -dy / dy_b * velocity.getVector(i, j - 1, k)[2];
    }

    // If right cell is fluid: two velocities have to be set: direction 1 and 2.
    if ((obstacle & OBSTACLE_RIGHT) == 0) {
      const RealType dx_r            = parameters_.meshsize->getDx(i + 1, j, k);
      const RealType dx              = parameters_.meshsize->getDx(i, j, k);
      velocity.getVector(i, j, k)[1] = -dx / dx_r * velocity.getVector(i + 1, j, k)[1];
      velocity.getVector(i, j, k)[2] = -dx / dx_r * velocity.getVector(i + 1, j, k)[2];
    }
    if ((obstacle & OBSTACLE_LEFT) == 0) {
      const RealType dx_l            = parameters_.meshsize->getDx(i - 1, j, k);
      const RealType dx              = parameters_.meshsize->getDx(i, j, k);
      velocity.getVector(i, j, k)[1] = -dx / dx_l * velocity.getVector(i - 1, j, k)[1];
      velocity.getVector(i, j, k)[2] = -dx / dx_l * velocity.getVector(i - 1, j, k)[2];
    }

    // Same for fluid cell in front
    if ((obstacle & OBSTACLE_BACK) == 0) {
      const RealType dz_f            = parameters_.meshsize->getDx(i, j, k + 1);
      const RealType dz              = parameters_.meshsize->getDx(i, j, k);
      velocity.getVector(i, j, k)[1] = -dz / dz_f * velocity.getVector(i, j, k + 1)[1];
      velocity.getVector(i, j, k)[0] = -dz / dz_f * velocity.getVector(i, j, k + 1)[0];
    }
    if ((obstacle & OBSTACLE_FRONT) == 0) {
      const RealType dz_b            = parameters_.meshsize->getDx(i, j, k - 1);
      const RealType dz              = parameters_.meshsize->getDx(i, j, k);
      velocity.getVector(i, j, k)[1] = -dz / dz_b * velocity.getVector(i, j, k - 1)[1];
      velocity.getVector(i, j, k)[0] = -dz / dz_b * velocity.getVector(i, j, k - 1)[0];
    }

    // Now the normal velocities need to be set to zero to ensure no flow at interfaces between solid and fluid.
    if ((obstacle & OBSTACLE_RIGHT) == 0) {
      velocity.getVector(i, j, k)[0] = 0.0;
    }
    if ((obstacle & OBSTACLE_TOP) == 0) {
      velocity.getVector(i, j, k)[1] = 0.0;
    }
    if ((obstacle & OBSTACLE_BACK) == 0) {
      velocity.getVector(i, j, k)[2] = 0.0;
    }
  }
// Minimal distance computation
    if (minimalDistances3D.empty()) {
        initializeDistances3D(flowField.getNx(), flowField.getNy(), flowField.getNz());
    }

    if ((obstacle & OBSTACLE_SELF) == 1) {
        minimalDistances3D[i][j][k] = 0.0; // Obstacle cells have zero distance
    } else {
        for (int di = -1; di <= 1; ++di) {
            for (int dj = -1; dj <= 1; ++dj) {
                for (int dk = -1; dk <= 1; ++dk) {
                    int ni = i + di;
                    int nj = j + dj;
                    int nk = k + dk;

                    if (ni >= 0 && ni < flowField.getNx() &&
                        nj >= 0 && nj < flowField.getNy() &&
                        nk >= 0 && nk < flowField.getNz()) {
                        if ((flowField.getFlags().getValue(ni, nj, nk) & OBSTACLE_SELF) == 1) {
                            const RealType dx = parameters_.meshsize->getDx(i, j, k);
                            const RealType dy = parameters_.meshsize->getDy(i, j, k);
                            const RealType dz = parameters_.meshsize->getDz(i, j, k);
                            RealType dist = std::sqrt(di * di * dx * dx +
                                                      dj * dj * dy * dy +
                                                      dk * dk * dz * dz);
                            minimalDistances3D[i][j][k] = std::min(minimalDistances3D[i][j][k], dist);
                        }
                    }
                }
            }
        }
    }
}
// Get minimal distance
RealType ObstacleStencil::getMinimalDistance(int i, int j) const {
    return minimalDistances2D[i][j];
}

RealType ObstacleStencil::getMinimalDistance(int i, int j, int k) const {
    return minimalDistances3D[i][j][k];
}

}