#pragma once

#include "Definitions.hpp"
#include "Parameters.hpp"

namespace Stencils {

  // Load the local velocity cube with relevant velocities of the 2D plane
  inline void loadLocalVelocity2D(FlowField& flowField, RealType* const localVelocity, int i, int j) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        const RealType* const point                  = flowField.getVelocity().getVector(i + column, j + row);
        localVelocity[39 + 9 * row + 3 * column]     = point[0]; // x-component
        localVelocity[39 + 9 * row + 3 * column + 1] = point[1]; // y-component
      }
    }
  }

  // Load the local velocity cube with surrounding velocities
  inline void loadLocalVelocity3D(FlowField& flowField, RealType* const localVelocity, int i, int j, int k) {
    for (int layer = -1; layer <= 1; layer++) {
      for (int row = -1; row <= 1; row++) {
        for (int column = -1; column <= 1; column++) {
          const RealType* const point = flowField.getVelocity().getVector(i + column, j + row, k + layer);
          localVelocity[39 + 27 * layer + 9 * row + 3 * column]     = point[0]; // x-component
          localVelocity[39 + 27 * layer + 9 * row + 3 * column + 1] = point[1]; // y-component
          localVelocity[39 + 27 * layer + 9 * row + 3 * column + 2] = point[2]; // z-component
        }
      }
    }
  }

  // Load local meshsize for 2D -> same as loadLocalVelocity2D, but invoking call to meshsize-ptr
  inline void loadLocalMeshsize2D(const Parameters& parameters, RealType* const localMeshsize, int i, int j) {
    for (int row = -1; row <= 1; row++) {
      for (int column = -1; column <= 1; column++) {
        localMeshsize[39 + 9 * row + 3 * column]     = parameters.meshsize->getDx(i + column, j + row);
        localMeshsize[39 + 9 * row + 3 * column + 1] = parameters.meshsize->getDy(i + column, j + row);
      }
    }
  }

  // Load local meshsize for 3D
  inline void loadLocalMeshsize3D(const Parameters& parameters, RealType* const localMeshsize, int i, int j, int k) {
    for (int layer = -1; layer <= 1; layer++) {
      for (int row = -1; row <= 1; row++) {
        for (int column = -1; column <= 1; column++) {
          localMeshsize[39 + 27 * layer + 9 * row + 3 * column] = parameters.meshsize->getDx(
            i + column, j + row, k + layer
          );
          localMeshsize[39 + 27 * layer + 9 * row + 3 * column + 1] = parameters.meshsize->getDy(
            i + column, j + row, k + layer
          );
          localMeshsize[39 + 27 * layer + 9 * row + 3 * column + 2] = parameters.meshsize->getDz(
            i + column, j + row, k + layer
          );
        }
      }
    }
  }

  // Maps an index and a component to the corresponding value in the cube.
  inline int mapd(int i, int j, int k, int component) { return 39 + 27 * k + 9 * j + 3 * i + component; }

  // Derivative functions. They are applied to a cube of 3x3x3 cells. lv stands for the local velocity, lm represents
  // the local mesh sizes dudx <-> first derivative of u-component of velocity field w.r.t. x-direction.
  inline RealType dudx(const RealType* const lv, const RealType* const lm) {
    // Evaluate dudx in the cell center by a central difference
    const int index0 = mapd(0, 0, 0, 0);
    const int index1 = mapd(-1, 0, 0, 0);
    return (lv[index0] - lv[index1]) / lm[index0];
  }

  inline RealType dvdy(const RealType* const lv, const RealType* const lm) {
    const int index0 = mapd(0, 0, 0, 1);
    const int index1 = mapd(0, -1, 0, 1);
    return (lv[index0] - lv[index1]) / lm[index0];
  }

  inline RealType dwdz(const RealType* const lv, const RealType* const lm) {
    const int index0 = mapd(0, 0, 0, 2);
    const int index1 = mapd(0, 0, -1, 2);
    return (lv[index0] - lv[index1]) / lm[index0];
  }

  //Other first derivatives
 // Derivative of u with respect to y (du/dy) for non-uniform grids
  inline RealType dudy(const RealType* const lv, const RealType* const lm) {
    // Velocities at adjacent interfaces
    const RealType u1 = lv[mapd(0, -1, 0, 0)]; // u(i, j-1)
    const RealType u2 = lv[mapd(0,  0, 0, 0)]; // u(i, j)
    const RealType u3 = lv[mapd(0,  1, 0, 0)]; // u(i, j+1)
    const RealType u4 = lv[mapd(-1, 1, 0, 0)]; // u(i-1, j+1)
    const RealType u5 = lv[mapd(-1, 0, 0, 0)]; // u(i-1, j)
    const RealType u6 = lv[mapd(-1, -1, 0, 0)]; // u(i-1, j-1)

    // Delta y of each cell
    const RealType dy0 = lm[mapd(0, -1, 0, 1)]; // dy(j-1)
    const RealType dy1 = lm[mapd(0,  0, 0, 1)]; // dy(j)
    const RealType dy2 = lm[mapd(0,  1, 0, 1)]; // dy(j+1)

    // Compute interpolation weights
    const RealType a0 = dy1 / (dy1 + dy0); // Weight for j-1
    const RealType a1 = dy0 / (dy1 + dy0); // Weight for j
    const RealType a2 = dy2 / (dy1 + dy2); // Weight for j+1
    const RealType a3 = dy1 / (dy1 + dy2); // Weight for j

    // Interpolate corner velocities at half cells
    const RealType uint1 = u1 * a0 + u2 * a1; // u(i, j-1/2)
    const RealType uint2 = u6 * a0 + u5 * a1; // u(i-1, j-1/2)
    const RealType uint3 = u2 * a2 + u3 * a3; // u(i, j+1/2)
    const RealType uint4 = u5 * a2 + u4 * a3; // u(i-1, j+1/2)

    // Compute interface velocities
    const RealType uy1 = 0.5 * (uint1 + uint2); // u(i-1/2, j-1/2)
    const RealType uy2 = 0.5 * (uint3 + uint4); // u(i-1/2, j+1/2)

    // Derivative of u with respect to y
    return (uy2 - uy1) / dy1;
}

  // Derivative of u with respect to z (du/dz) for non-uniform grids
  inline RealType dudz(const RealType* const lv, const RealType* const lm) {
    // Velocities at adjacent interfaces
    const RealType u1 = lv[mapd(0, 0, -1, 0)]; // u(i, j, k-1)
    const RealType u2 = lv[mapd(0, 0,  0, 0)]; // u(i, j, k)
    const RealType u3 = lv[mapd(0, 0,  1, 0)]; // u(i, j, k+1)
    const RealType u4 = lv[mapd(-1, 0, 1, 0)]; // u(i-1, j, k+1)
    const RealType u5 = lv[mapd(-1, 0, 0, 0)]; // u(i-1, j, k)
    const RealType u6 = lv[mapd(-1, 0, -1, 0)]; // u(i-1, j, k-1)

    // Delta z of each cell
    const RealType dz0 = lm[mapd(0, 0, -1, 2)]; // dz(k-1)
    const RealType dz1 = lm[mapd(0, 0,  0, 2)]; // dz(k)
    const RealType dz2 = lm[mapd(0, 0,  1, 2)]; // dz(k+1)

    // Compute interpolation weights
    const RealType a0 = dz1 / (dz1 + dz0); // Weight for k-1
    const RealType a1 = dz0 / (dz1 + dz0); // Weight for k
    const RealType a2 = dz2 / (dz1 + dz2); // Weight for k+1
    const RealType a3 = dz1 / (dz1 + dz2); // Weight for k

    // Interpolate corner velocities at half cells
    const RealType uint1 = u1 * a0 + u2 * a1; // u(i, j, k-1/2)
    const RealType uint2 = u6 * a0 + u5 * a1; // u(i-1, j, k-1/2)
    const RealType uint3 = u2 * a2 + u3 * a3; // u(i, j, k+1/2)
    const RealType uint4 = u5 * a2 + u4 * a3; // u(i-1, j, k+1/2)

    // Compute interface velocities
    const RealType uz1 = 0.5 * (uint1 + uint2); // u(i-1/2, j, k-1/2)
    const RealType uz2 = 0.5 * (uint3 + uint4); // u(i-1/2, j, k+1/2)

    // Derivative of u with respect to z
    return (uz2 - uz1) / dz1;
}

  // Derivative of v with respect to x (dv/dx) for non-uniform grids
  inline RealType dvdx(const RealType* const lv, const RealType* const lm) {
    // Velocities at adjacent interfaces
    const RealType v1 = lv[mapd(-1, 0, 0, 1)]; // v(i-1, j, k)
    const RealType v2 = lv[mapd( 0, 0, 0, 1)]; // v(i, j, k)
    const RealType v3 = lv[mapd( 1, 0, 0, 1)]; // v(i+1, j, k)
    const RealType v4 = lv[mapd( 1, -1, 0, 1)]; // v(i+1, j-1, k)
    const RealType v5 = lv[mapd( 0, -1, 0, 1)]; // v(i, j-1, k)
    const RealType v6 = lv[mapd(-1, -1, 0, 1)]; // v(i-1, j-1, k)

    // Delta x of each cell
    const RealType dx0 = lm[mapd(-1, 0, 0, 0)]; // dx(i-1)
    const RealType dx1 = lm[mapd( 0, 0, 0, 0)]; // dx(i)
    const RealType dx2 = lm[mapd( 1, 0, 0, 0)]; // dx(i+1)

    // Compute interpolation weights
    const RealType a0 = dx1 / (dx1 + dx0); // Weight for i-1
    const RealType a1 = dx0 / (dx1 + dx0); // Weight for i
    const RealType a2 = dx2 / (dx1 + dx2); // Weight for i+1
    const RealType a3 = dx1 / (dx1 + dx2); // Weight for i

    // Interpolate corner velocities at half cells
    const RealType vint1 = v1 * a0 + v2 * a1; // v(i-1/2, j, k)
    const RealType vint2 = v6 * a0 + v5 * a1; // v(i-1/2, j-1, k)
    const RealType vint3 = v2 * a2 + v3 * a3; // v(i+1/2, j, k)
    const RealType vint4 = v5 * a2 + v4 * a3; // v(i+1/2, j-1, k)

    // Compute interface velocities
    const RealType vx1 = 0.5 * (vint1 + vint2); // v(i-1/2, j-1/2, k)
    const RealType vx2 = 0.5 * (vint3 + vint4); // v(i+1/2, j-1/2, k)

    // Derivative of v with respect to x
    return (vx2 - vx1) / dx1;
}

  // Derivative of v with respect to z (dv/dz) for non-uniform grids
  inline RealType dvdz(const RealType* const lv, const RealType* const lm) {
    // Velocities at adjacent interfaces
    const RealType v1 = lv[mapd(0, 0, -1, 1)]; // v(i, j, k-1)
    const RealType v2 = lv[mapd(0, 0,  0, 1)]; // v(i, j, k)
    const RealType v3 = lv[mapd(0, 0,  1, 1)]; // v(i, j, k+1)
    const RealType v4 = lv[mapd(0, -1, 1, 1)]; // v(i, j-1, k+1)
    const RealType v5 = lv[mapd(0, -1, 0, 1)]; // v(i, j-1, k)
    const RealType v6 = lv[mapd(0, -1, -1, 1)]; // v(i, j-1, k-1)

    // Delta z of each cell
    const RealType dz0 = lm[mapd(0, 0, -1, 2)]; // dz(k-1)
    const RealType dz1 = lm[mapd(0, 0,  0, 2)]; // dz(k)
    const RealType dz2 = lm[mapd(0, 0,  1, 2)]; // dz(k+1)

    // Compute interpolation weights
    const RealType a0 = dz1 / (dz1 + dz0); // Weight for k-1
    const RealType a1 = dz0 / (dz1 + dz0); // Weight for k
    const RealType a2 = dz2 / (dz1 + dz2); // Weight for k+1
    const RealType a3 = dz1 / (dz1 + dz2); // Weight for k

    // Interpolate corner velocities at half cells
    const RealType vint1 = v1 * a0 + v2 * a1; // v(i, j, k-1/2)
    const RealType vint2 = v6 * a0 + v5 * a1; // v(i, j-1, k-1/2)
    const RealType vint3 = v2 * a2 + v3 * a3; // v(i, j, k+1/2)
    const RealType vint4 = v5 * a2 + v4 * a3; // v(i, j-1, k+1/2)

    // Compute interface velocities
    const RealType vz1 = 0.5 * (vint1 + vint2); // v(i, j-1/2, k-1/2)
    const RealType vz2 = 0.5 * (vint3 + vint4); // v(i, j-1/2, k+1/2)

    // Derivative of v with respect to z
    return (vz2 - vz1) / dz1;
}

  // Derivative of w with respect to x (dw/dx) for non-uniform grids
  inline RealType dwdx(const RealType* const lv, const RealType* const lm) {
    // Velocities at adjacent interfaces
    const RealType w1 = lv[mapd(-1, 0, 0, 2)]; // w(i-1, j, k)
    const RealType w2 = lv[mapd( 0, 0, 0, 2)]; // w(i, j, k)
    const RealType w3 = lv[mapd( 1, 0, 0, 2)]; // w(i+1, j, k)
    const RealType w4 = lv[mapd( 1, -1, 0, 2)]; // w(i+1, j-1, k)
    const RealType w5 = lv[mapd( 0, -1, 0, 2)]; // w(i, j-1, k)
    const RealType w6 = lv[mapd(-1, -1, 0, 2)]; // w(i-1, j-1, k)

    // Delta x of each cell
    const RealType dx0 = lm[mapd(-1, 0, 0, 0)]; // dx(i-1)
    const RealType dx1 = lm[mapd( 0, 0, 0, 0)]; // dx(i)
    const RealType dx2 = lm[mapd( 1, 0, 0, 0)]; // dx(i+1)

    // Compute interpolation weights
    const RealType a0 = dx1 / (dx1 + dx0); // Weight for i-1
    const RealType a1 = dx0 / (dx1 + dx0); // Weight for i
    const RealType a2 = dx2 / (dx1 + dx2); // Weight for i+1
    const RealType a3 = dx1 / (dx1 + dx2); // Weight for i

    // Interpolate corner velocities at half cells
    const RealType wint1 = w1 * a0 + w2 * a1; // w(i-1/2, j, k)
    const RealType wint2 = w6 * a0 + w5 * a1; // w(i-1/2, j-1, k)
    const RealType wint3 = w2 * a2 + w3 * a3; // w(i+1/2, j, k)
    const RealType wint4 = w5 * a2 + w4 * a3; // w(i+1/2, j-1, k)

    // Compute interface velocities
    const RealType wx1 = 0.5 * (wint1 + wint2); // w(i-1/2, j-1/2, k)
    const RealType wx2 = 0.5 * (wint3 + wint4); // w(i+1/2, j-1/2, k)

    // Derivative of w with respect to x
    return (wx2 - wx1) / dx1;
}

  // Derivative of w with respect to y (dw/dy) for non-uniform grids
  inline RealType dwdy(const RealType* const lv, const RealType* const lm) {
    // Velocities at adjacent interfaces
    const RealType w1 = lv[mapd(0, -1, 0, 2)]; // w(i, j-1, k)
    const RealType w2 = lv[mapd(0,  0, 0, 2)]; // w(i, j, k)
    const RealType w3 = lv[mapd(0,  1, 0, 2)]; // w(i, j+1, k)
    const RealType w4 = lv[mapd(-1, 1, 0, 2)]; // w(i-1, j+1, k)
    const RealType w5 = lv[mapd(-1, 0, 0, 2)]; // w(i-1, j, k)
    const RealType w6 = lv[mapd(-1, -1, 0, 2)]; // w(i-1, j-1, k)

    // Delta y of each cell
    const RealType dy0 = lm[mapd(0, -1, 0, 1)]; // dy(j-1)
    const RealType dy1 = lm[mapd(0,  0, 0, 1)]; // dy(j)
    const RealType dy2 = lm[mapd(0,  1, 0, 1)]; // dy(j+1)

    // Compute interpolation weights
    const RealType a0 = dy1 / (dy1 + dy0); // Weight for j-1
    const RealType a1 = dy0 / (dy1 + dy0); // Weight for j
    const RealType a2 = dy2 / (dy1 + dy2); // Weight for j+1
    const RealType a3 = dy1 / (dy1 + dy2); // Weight for j

    // Interpolate corner velocities at half cells
    const RealType wint1 = w1 * a0 + w2 * a1; // w(i, j-1/2, k)
    const RealType wint2 = w6 * a0 + w5 * a1; // w(i-1, j-1/2, k)
    const RealType wint3 = w2 * a2 + w3 * a3; // w(i, j+1/2, k)
    const RealType wint4 = w5 * a2 + w4 * a3; // w(i-1, j+1/2, k)

    // Compute interface velocities
    const RealType wy1 = 0.5 * (wint1 + wint2); // w(i-1/2, j-1/2, k)
    const RealType wy2 = 0.5 * (wint3 + wint4); // w(i-1/2, j+1/2, k)

    // Derivative of w with respect to y
    return (wy2 - wy1) / dy1;
}




    // TODO WS1: Second derivatives

  // First derivative of product (u*v), evaluated at the location of the v-component.
  inline RealType duvdx(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 1, 0, 0)]) *
        (lv[mapd(0, 0, 0, 1)] + lv[mapd(1, 0, 0, 1)])) -
        ((lv[mapd(-1, 0, 0, 0)] + lv[mapd(-1, 1, 0, 0)]) *
            (lv[mapd(-1, 0, 0, 1)] + lv[mapd(0, 0, 0, 1)])))
        + parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 1, 0, 0)]) *
            (lv[mapd(0, 0, 0, 1)] - lv[mapd(1, 0, 0, 1)])) -
            (fabs(lv[mapd(-1, 0, 0, 0)] + lv[mapd(-1, 1, 0, 0)]) *
                (lv[mapd(-1, 0, 0, 1)] - lv[mapd(0, 0, 0, 1)])))
        ) / lm[mapd(0, 0, 0, 0)];
#endif

    const RealType hxShort = 0.5 * lm[mapd(0, 0, 0, 0)]; // Distance of corner points in x-direction from center v-value
    const RealType hxLong0 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(-1, 0, 0, 0)]); // Distance between center and west
                                                                                   // v-value
    const RealType hxLong1 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hyShort = 0.5 * lm[mapd(0, 0, 0, 1)]; // Distance of center u-value from upper edge of cell
    const RealType hyLong = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]); // Distance of north and center u-value

    const RealType u00 = lv[mapd(0, 0, 0, 0)];
    const RealType u01 = lv[mapd(0, 1, 0, 0)];
    const RealType v00 = lv[mapd(0, 0, 0, 1)];
    const RealType v10 = lv[mapd(1, 0, 0, 1)];

    const RealType uM10 = lv[mapd(-1, 0, 0, 0)];
    const RealType uM11 = lv[mapd(-1, 1, 0, 0)];
    const RealType vM10 = lv[mapd(-1, 0, 0, 1)];

    // This a central difference expression for the first-derivative. We therefore linearly interpolate u*v onto the
    // surface of the current cell (in 2D: upper left and upper right corner) and then take the central difference.
    const RealType secondOrder = (((hyLong - hyShort) / hyLong * u00 + hyShort / hyLong * u01
                                  ) * ((hxLong1 - hxShort) / hxLong1 * v00 + hxShort / hxLong1 * v10)
                                  - ((hyLong - hyShort) / hyLong * uM10 + hyShort / hyLong * uM11
                                    ) * ((hxLong0 - hxShort) / hxLong0 * v00 + hxShort / hxLong0 * vM10))
                                 / (2.0 * hxShort);

    // This is a forward-difference in donor-cell style. We apply donor cell and again interpolate the velocity values
    // (u-comp.) onto the surface of the cell. We then apply the standard donor cell scheme. This will, however, result
    // in non-equal mesh spacing evaluations (in case of stretched meshes).
    const RealType kr = (hyLong - hyShort) / hyLong * u00 + hyShort / hyLong * u01;
    const RealType kl = (hyLong - hyShort) / hyLong * uM10 + hyShort / hyLong * uM11;

    const RealType firstOrder
      = 1.0 / (4.0 * hxShort)
        * (kr * (v00 + v10) - kl * (vM10 + v00) + fabs(kr) * (v00 - v10) - fabs(kl) * (vM10 - v00));

    // Return linear combination of central and donor-cell difference
    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in duvdx");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. y for u*v at location of u-component. For details on implementation, see duvdx.
  inline RealType duvdy(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 1)] + lv[mapd(1, 0, 0, 1)]) *
        (lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 1, 0, 0)])) -
        ((lv[mapd(0, -1, 0, 1)] + lv[mapd(1, -1, 0, 1)]) *
            (lv[mapd(0, -1, 0, 0)] + lv[mapd(0, 0, 0, 0)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 1)] + lv[mapd(1, 0, 0, 1)]) *
            (lv[mapd(0, 0, 0, 0)] - lv[mapd(0, 1, 0, 0)])) -
            (fabs(lv[mapd(0, -1, 0, 1)] + lv[mapd(1, -1, 0, 1)]) *
                (lv[mapd(0, -1, 0, 0)] - lv[mapd(0, 0, 0, 0)])))) /
        lm[mapd(0, 0, 0, 1)];
#endif

    const RealType hyShort = 0.5 * lm[mapd(0, 0, 0, 1)]; // Distance of corner points in x-direction from center v-value
    const RealType hyLong0 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, -1, 0, 1)]); // Distance between center and west
                                                                                   // v-value
    const RealType hyLong1 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hxShort = 0.5 * lm[mapd(0, 0, 0, 0)]; // Distance of center u-value from upper edge of cell
    const RealType hxLong = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]); // Distance of north and center u-value

    const RealType v00 = lv[mapd(0, 0, 0, 1)];
    const RealType v10 = lv[mapd(1, 0, 0, 1)];
    const RealType u00 = lv[mapd(0, 0, 0, 0)];
    const RealType u01 = lv[mapd(0, 1, 0, 0)];

    const RealType v0M1 = lv[mapd(0, -1, 0, 1)];
    const RealType v1M1 = lv[mapd(1, -1, 0, 1)];
    const RealType u0M1 = lv[mapd(0, -1, 0, 0)];

    const RealType secondOrder = (((hxLong - hxShort) / hxLong * v00 + hxShort / hxLong * v10
                                  ) * ((hyLong1 - hyShort) / hyLong1 * u00 + hyShort / hyLong1 * u01)
                                  - ((hxLong - hxShort) / hxLong * v0M1 + hxShort / hxLong * v1M1
                                    ) * ((hyLong0 - hyShort) / hyLong0 * u00 + hyShort / hyLong0 * u0M1))
                                 / (2.0 * hyShort);

    const RealType kr = (hxLong - hxShort) / hxLong * v00 + hxShort / hxLong * v10;
    const RealType kl = (hxLong - hxShort) / hxLong * v0M1 + hxShort / hxLong * v1M1;

    const RealType firstOrder
      = 1.0 / (4.0 * hyShort)
        * (kr * (u00 + u01) - kl * (u0M1 + u00) + fabs(kr) * (u00 - u01) - fabs(kl) * (u0M1 - u00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in duvdy");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. x for u*w at location of w-component. For details on implementation, see duvdx.
  inline RealType duwdx(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 0, 1, 0)]) *
        (lv[mapd(0, 0, 0, 2)] + lv[mapd(1, 0, 0, 2)])) -
        ((lv[mapd(-1, 0, 0, 0)] + lv[mapd(-1, 0, 1, 0)]) *
            (lv[mapd(-1, 0, 0, 2)] + lv[mapd(0, 0, 0, 2)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 0, 1, 0)]) *
            (lv[mapd(0, 0, 0, 2)] - lv[mapd(1, 0, 0, 2)])) -
            (fabs(lv[mapd(-1, 0, 0, 0)] + lv[mapd(-1, 0, 1, 0)]) *
                (lv[mapd(-1, 0, 0, 2)] - lv[mapd(0, 0, 0, 2)])))) /
        lm[mapd(0, 0, 0, 0)];
#endif

    const RealType hxShort = 0.5 * lm[mapd(0, 0, 0, 0)]; // Distance of corner points in x-direction from center v-value
    const RealType hxLong0 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(-1, 0, 0, 0)]); // Distance between center and west
                                                                                   // v-value
    const RealType hxLong1 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hzShort = 0.5 * lm[mapd(0, 0, 0, 2)]; // Distance of center u-value from upper edge of cell
    const RealType hzLong = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]); // Distance of north and center u-value

    const RealType u00 = lv[mapd(0, 0, 0, 0)];
    const RealType u01 = lv[mapd(0, 0, 1, 0)];
    const RealType w00 = lv[mapd(0, 0, 0, 2)];
    const RealType w10 = lv[mapd(1, 0, 0, 2)];

    const RealType uM10 = lv[mapd(-1, 0, 0, 0)];
    const RealType uM11 = lv[mapd(-1, 0, 1, 0)];
    const RealType wM10 = lv[mapd(-1, 0, 0, 2)];

    const RealType secondOrder = (((hzLong - hzShort) / hzLong * u00 + hzShort / hzLong * u01
                                  ) * ((hxLong1 - hxShort) / hxLong1 * w00 + hxShort / hxLong1 * w10)
                                  - ((hzLong - hzShort) / hzLong * uM10 + hzShort / hzLong * uM11
                                    ) * ((hxLong0 - hxShort) / hxLong0 * w00 + hxShort / hxLong0 * wM10))
                                 / (2.0 * hxShort);

    const RealType kr = (hzLong - hzShort) / hzLong * u00 + hzShort / hzLong * u01;
    const RealType kl = (hzLong - hzShort) / hzLong * uM10 + hzShort / hzLong * uM11;

    const RealType firstOrder
      = 1.0 / (4.0 * hxShort)
        * (kr * (w00 + w10) - kl * (wM10 + w00) + fabs(kr) * (w00 - w10) - fabs(kl) * (wM10 - w00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in duwdx");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. z for u*w at location of u-component. For details on implementation, see duvdx.
  inline RealType duwdz(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 2)] + lv[mapd(1, 0, 0, 2)]) *
        (lv[mapd(0, 0, 0, 0)] + lv[mapd(0, 0, 1, 0)])) -
        ((lv[mapd(0, 0, -1, 2)] + lv[mapd(1, 0, -1, 2)]) *
            (lv[mapd(0, 0, -1, 0)] + lv[mapd(0, 0, 0, 0)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 2)] + lv[mapd(1, 0, 0, 2)]) *
            (lv[mapd(0, 0, 0, 0)] - lv[mapd(0, 0, 1, 0)])) -
            (fabs(lv[mapd(0, 0, -1, 2)] + lv[mapd(1, 0, -1, 2)]) *
                (lv[mapd(0, 0, -1, 0)] - lv[mapd(0, 0, 0, 0)])))) /
        lm[mapd(0, 0, 0, 2)];
#endif

    const RealType hzShort = 0.5 * lm[mapd(0, 0, 0, 2)]; // Distance of corner points in x-direction from center v-value
    const RealType hzLong0 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, -1, 2)]); // Distance between center and west
                                                                                   // v-value
    const RealType hzLong1 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hxShort = 0.5 * lm[mapd(0, 0, 0, 0)]; // Distance of center u-value from upper edge of cell
    const RealType hxLong = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]); // Distance of north and center u-value

    const RealType w00 = lv[mapd(0, 0, 0, 2)];
    const RealType w10 = lv[mapd(1, 0, 0, 2)];
    const RealType u00 = lv[mapd(0, 0, 0, 0)];
    const RealType u01 = lv[mapd(0, 0, 1, 0)];

    const RealType w0M1 = lv[mapd(0, 0, -1, 2)];
    const RealType w1M1 = lv[mapd(1, 0, -1, 2)];
    const RealType u0M1 = lv[mapd(0, 0, -1, 0)];

    const RealType secondOrder = (((hxLong - hxShort) / hxLong * w00 + hxShort / hxLong * w10
                                  ) * ((hzLong1 - hzShort) / hzLong1 * u00 + hzShort / hzLong1 * u01)
                                  - ((hxLong - hxShort) / hxLong * w0M1 + hxShort / hxLong * w1M1
                                    ) * ((hzLong0 - hzShort) / hzLong0 * u00 + hzShort / hzLong0 * u0M1))
                                 / (2.0 * hzShort);

    const RealType kr = (hxLong - hxShort) / hxLong * w00 + hxShort / hxLong * w10;
    const RealType kl = (hxLong - hxShort) / hxLong * w0M1 + hxShort / hxLong * w1M1;

    const RealType firstOrder
      = 1.0 / (4.0 * hzShort)
        * (kr * (u00 + u01) - kl * (u0M1 + u00) + fabs(kr) * (u00 - u01) - fabs(kl) * (u0M1 - u00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in duwdz");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. y for v*w at location of w-component. For details on implementation, see duvdx.
  inline RealType dvwdy(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 0, 1, 1)]) *
        (lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 1, 0, 2)])) -
        ((lv[mapd(0, -1, 0, 1)] + lv[mapd(0, -1, 1, 1)]) *
            (lv[mapd(0, -1, 0, 2)] + lv[mapd(0, 0, 0, 2)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 0, 1, 1)]) *
            (lv[mapd(0, 0, 0, 2)] - lv[mapd(0, 1, 0, 2)])) -
            (fabs(lv[mapd(0, -1, 0, 1)] + lv[mapd(0, -1, 1, 1)]) *
                (lv[mapd(0, -1, 0, 2)] - lv[mapd(0, 0, 0, 2)])))) /
        lm[mapd(0, 0, 0, 1)];
#endif

    const RealType hyShort = 0.5 * lm[mapd(0, 0, 0, 1)]; // Distance of corner points in x-direction from center v-value
    const RealType hyLong0 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, -1, 0, 1)]); // Distance between center and west
                                                                                   // v-value
    const RealType hyLong1 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hzShort = 0.5 * lm[mapd(0, 0, 0, 2)]; // Distance of center u-value from upper edge of cell
    const RealType hzLong = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]); // Distance of north and center u-value

    const RealType v00 = lv[mapd(0, 0, 0, 1)];
    const RealType v01 = lv[mapd(0, 0, 1, 1)];
    const RealType w00 = lv[mapd(0, 0, 0, 2)];
    const RealType w10 = lv[mapd(0, 1, 0, 2)];

    const RealType vM10 = lv[mapd(0, -1, 0, 1)];
    const RealType vM11 = lv[mapd(0, -1, 1, 1)];
    const RealType wM10 = lv[mapd(0, -1, 0, 2)];

    const RealType secondOrder = (((hzLong - hzShort) / hzLong * v00 + hzShort / hzLong * v01
                                  ) * ((hyLong1 - hyShort) / hyLong1 * w00 + hyShort / hyLong1 * w10)
                                  - ((hzLong - hzShort) / hzLong * vM10 + hzShort / hzLong * vM11
                                    ) * ((hyLong0 - hyShort) / hyLong0 * w00 + hyShort / hyLong0 * wM10))
                                 / (2.0 * hyShort);

    const RealType kr = (hzLong - hzShort) / hzLong * v00 + hzShort / hzLong * v01;
    const RealType kl = (hzLong - hzShort) / hzLong * vM10 + hzShort / hzLong * vM11;

    const RealType firstOrder
      = 1.0 / (4.0 * hyShort)
        * (kr * (w00 + w10) - kl * (wM10 + w00) + fabs(kr) * (w00 - w10) - fabs(kl) * (wM10 - w00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      throw std::runtime_error("Error in dvwdy");
    }
#endif

    return tmp2;
  }

  // Evaluates first derivative w.r.t. z for v*w at location of v-component. For details on implementation, see duvdx.
  inline RealType dvwdz(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 1, 0, 2)]) *
        (lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 0, 1, 1)])) -
        ((lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 1, -1, 2)]) *
            (lv[mapd(0, 0, -1, 1)] + lv[mapd(0, 0, 0, 1)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 1, 0, 2)]) *
            (lv[mapd(0, 0, 0, 1)] - lv[mapd(0, 0, 1, 1)])) -
            (fabs(lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 1, -1, 2)]) *
                (lv[mapd(0, 0, -1, 1)] - lv[mapd(0, 0, 0, 1)])))) /
        lm[mapd(0, 0, 0, 2)];
#endif

    const RealType hzShort = 0.5 * lm[mapd(0, 0, 0, 2)]; // Distance of corner points in x-direction from center v-value
    const RealType hzLong0 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, -1, 2)]); // Distance between center and west
                                                                                   // v-value
    const RealType hzLong1 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);  // Distance between center and east
                                                                                   // v-value
    const RealType hyShort = 0.5 * lm[mapd(0, 0, 0, 1)]; // Distance of center u-value from upper edge of cell
    const RealType hyLong = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]); // Distance of north and center u-value

    const RealType w00 = lv[mapd(0, 0, 0, 2)];
    const RealType w10 = lv[mapd(0, 1, 0, 2)];
    const RealType v00 = lv[mapd(0, 0, 0, 1)];
    const RealType v01 = lv[mapd(0, 0, 1, 1)];

    const RealType w0M1 = lv[mapd(0, 0, -1, 2)];
    const RealType w1M1 = lv[mapd(0, 1, -1, 2)];
    const RealType v0M1 = lv[mapd(0, 0, -1, 1)];

    const RealType secondOrder = (((hyLong - hyShort) / hyLong * w00 + hyShort / hyLong * w10
                                  ) * ((hzLong1 - hzShort) / hzLong1 * v00 + hzShort / hzLong1 * v01)
                                  - ((hyLong - hyShort) / hyLong * w0M1 + hyShort / hyLong * w1M1
                                    ) * ((hzLong0 - hzShort) / hzLong0 * v00 + hzShort / hzLong0 * v0M1))
                                 / (2.0 * hzShort);

    const RealType kr = (hyLong - hyShort) / hyLong * w00 + hyShort / hyLong * w10;
    const RealType kl = (hyLong - hyShort) / hyLong * w0M1 + hyShort / hyLong * w1M1;

    const RealType firstOrder
      = 1.0 / (4.0 * hzShort)
        * (kr * (v00 + v01) - kl * (v0M1 + v00) + fabs(kr) * (v00 - v01) - fabs(kl) * (v0M1 - v00));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      spdlog::error("{}, {}", tmp1, tmp2);
      throw std::runtime_error("Error in dvwdz");
    }
#endif

    return tmp2;
  }

  // First derivative of u*u w.r.t. x, evaluated at location of u-component.
  inline RealType du2dx(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 0)] + lv[mapd(1, 0, 0, 0)]) *
        (lv[mapd(0, 0, 0, 0)] + lv[mapd(1, 0, 0, 0)])) -
        ((lv[mapd(-1, 0, 0, 0)] + lv[mapd(0, 0, 0, 0)]) *
            (lv[mapd(-1, 0, 0, 0)] + lv[mapd(0, 0, 0, 0)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 0)] + lv[mapd(1, 0, 0, 0)]) *
            (lv[mapd(0, 0, 0, 0)] - lv[mapd(1, 0, 0, 0)])) -
            (fabs(lv[mapd(-1, 0, 0, 0)] + lv[mapd(0, 0, 0, 0)]) *
                (lv[mapd(-1, 0, 0, 0)] - lv[mapd(0, 0, 0, 0)])))) /
        lm[mapd(0, 0, 0, 0)];
#endif

    const RealType dxShort = 0.5 * lm[mapd(0, 0, 0, 0)];
    // const RealType dxLong0 = 0.5*(lm[mapd(-1,0,0,0)] + lm[mapd(0,0,0,0)]);
    const RealType dxLong1 = 0.5 * (lm[mapd(0, 0, 0, 0)] + lm[mapd(1, 0, 0, 0)]);

    const RealType u0  = lv[mapd(0, 0, 0, 0)];
    const RealType uM1 = lv[mapd(-1, 0, 0, 0)];
    const RealType u1  = lv[mapd(1, 0, 0, 0)];

    // const RealType kr = (dxLong1 - dxShort) / dxLong1 * u0 + dxShort / dxLong1 * u1;
    // const RealType kl = (dxLong0 - dxShort) / dxLong0 * u0 + dxShort / dxLong0 * uM1;
    const RealType kr = (u0 + u1) / 2;
    const RealType kl = (u0 + uM1) / 2;

    // Central difference expression which is second-order accurate for uniform meshes. We interpolate u half-way
    // between neighboured u-component values and afterwards build the central difference for u*u.

    /*const RealType secondOrder = (((dxLong1 - dxShort) / dxLong1 * u0 + dxShort / dxLong1 * u1) * ((dxLong1 - dxShort)
       / dxLong1 * u0 + dxShort / dxLong1 * u1)
        - ((dxLong0 - dxShort) / dxLong0 * u0 + dxShort / dxLong0 * uM1) * ((dxLong0 - dxShort) / dxLong0 * u0 + dxShort
       / dxLong0 * uM1) ) / (2.0 * dxShort);*/

    const RealType secondOrder = ((u0 + u1) * (u0 + u1) - (u0 + uM1) * (u0 + uM1)) / (4 * dxLong1);

    // Donor-cell like derivative expression. We evaluate u half-way between neighboured u-components and use this as a
    // prediction of the transport direction.
    const RealType firstOrder = 1.0 / (4.0 * dxShort)
                                * (kr * (u0 + u1) - kl * (uM1 + u0) + fabs(kr) * (u0 - u1) - fabs(kl) * (uM1 - u0));

    // Return linear combination of central- and upwind difference
    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      spdlog::error("{}, {}", tmp1, tmp2);
      throw std::runtime_error("Error in du2dx");
    }
#endif

    return tmp2;
  }

  // First derivative of v*v w.r.t. y, evaluated at location of v-component. For details, see du2dx.
  inline RealType dv2dy(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 1, 0, 1)]) *
        (lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 1, 0, 1)])) -
        ((lv[mapd(0, -1, 0, 1)] + lv[mapd(0, 0, 0, 1)]) *
            (lv[mapd(0, -1, 0, 1)] + lv[mapd(0, 0, 0, 1)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 1)] + lv[mapd(0, 1, 0, 1)]) *
            (lv[mapd(0, 0, 0, 1)] - lv[mapd(0, 1, 0, 1)])) -
            (fabs(lv[mapd(0, -1, 0, 1)] + lv[mapd(0, 0, 0, 1)]) *
                (lv[mapd(0, -1, 0, 1)] - lv[mapd(0, 0, 0, 1)])))) /
        lm[mapd(0, 0, 0, 1)];
#endif

    const RealType dyShort = 0.5 * lm[mapd(0, 0, 0, 1)];
    // const RealType dyLong0 = 0.5*(lm[mapd(0,-1,0,1)] + lm[mapd(0,0,0,1)]);
    const RealType dyLong1 = 0.5 * (lm[mapd(0, 0, 0, 1)] + lm[mapd(0, 1, 0, 1)]);

    const RealType v0  = lv[mapd(0, 0, 0, 1)];
    const RealType vM1 = lv[mapd(0, -1, 0, 1)];
    const RealType v1  = lv[mapd(0, 1, 0, 1)];

    // const RealType kr = (dyLong1 - dyShort) / dyLong1 * v0 + dyShort / dyLong1 * v1;
    // const RealType kl = (dyLong0 - dyShort) / dyLong0 * v0 + dyShort / dyLong0 * vM1;
    const RealType kr = (v0 + v1) / 2;
    const RealType kl = (v0 + vM1) / 2;

    /*const RealType secondOrder = (((dyLong1 - dyShort) / dyLong1 * v0 + dyShort / dyLong1 * v1) * ((dyLong1 - dyShort)
       / dyLong1 * v0 + dyShort / dyLong1 * v1)
        - ((dyLong0 - dyShort) / dyLong0 * v0 + dyShort / dyLong0 * vM1) * ((dyLong0 - dyShort) / dyLong0 * v0 + dyShort
       / dyLong0 * vM1) ) / (2.0 * dyShort);*/

    const RealType secondOrder = ((v0 + v1) * (v0 + v1) - (v0 + vM1) * (v0 + vM1)) / (4 * dyLong1);

    const RealType firstOrder = 1.0 / (4.0 * dyShort)
                                * (kr * (v0 + v1) - kl * (vM1 + v0) + fabs(kr) * (v0 - v1) - fabs(kl) * (vM1 - v0));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      spdlog::error("{}, {}", tmp1, tmp2);
      throw std::runtime_error("Error in dv2dy");
    }
#endif

    return tmp2;
  }

  // First derivative of w*w w.r.t. z, evaluated at location of w-component. For details, see du2dx.
  inline RealType dw2dz(const RealType* const lv, const Parameters& parameters, const RealType* const lm) {
#ifndef NDEBUG
    const RealType tmp1 = 1.0 / 4.0 * ((((lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 0, 1, 2)]) *
        (lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 0, 1, 2)])) -
        ((lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 0, 0, 2)]) *
            (lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 0, 0, 2)]))) +
        parameters.solver.gamma * ((fabs(lv[mapd(0, 0, 0, 2)] + lv[mapd(0, 0, 1, 2)]) *
            (lv[mapd(0, 0, 0, 2)] - lv[mapd(0, 0, 1, 2)])) -
            (fabs(lv[mapd(0, 0, -1, 2)] + lv[mapd(0, 0, 0, 2)]) *
                (lv[mapd(0, 0, -1, 2)] - lv[mapd(0, 0, 0, 2)])))) /
        lm[mapd(0, 0, 0, 2)];
#endif

    const RealType dzShort = 0.5 * lm[mapd(0, 0, 0, 2)];
    // const RealType dzLong0 = 0.5 * (lm[mapd(0, 0, -1, 2)] + lm[mapd(0, 0, 0, 2)]);
    const RealType dzLong1 = 0.5 * (lm[mapd(0, 0, 0, 2)] + lm[mapd(0, 0, 1, 2)]);

    const RealType w0  = lv[mapd(0, 0, 0, 2)];
    const RealType wM1 = lv[mapd(0, 0, -1, 2)];
    const RealType w1  = lv[mapd(0, 0, 1, 2)];

    // const RealType kr = (dzLong1 - dzShort) / dzLong1 * w0 + dzShort / dzLong1 * w1;
    // const RealType kl = (dzLong0 - dzShort) / dzLong0 * w0 + dzShort / dzLong0 * wM1;
    const RealType kr = (w0 + w1) / 2;
    const RealType kl = (w0 + wM1) / 2;

    /*const RealType secondOrder = (((dzLong1 - dzShort) / dzLong1 * w0 + dzShort / dzLong1 * w1) * ((dzLong1 - dzShort)
       / dzLong1 * w0 + dzShort / dzLong1 * w1)
        - ((dzLong0 - dzShort) / dzLong0 * w0 + dzShort / dzLong0 * wM1) * ((dzLong0 - dzShort) / dzLong0 * w0 + dzShort
       / dzLong0 * wM1) ) / (2.0 * dzShort);*/

    const RealType secondOrder = ((w0 + w1) * (w0 + w1) - (w0 + wM1) * (w0 + wM1)) / (4 * dzLong1);

    const RealType firstOrder = 1.0 / (4.0 * dzShort)
                                * (kr * (w0 + w1) - kl * (wM1 + w0) + fabs(kr) * (w0 - w1) - fabs(kl) * (wM1 - w0));

    const RealType tmp2 = (1.0 - parameters.solver.gamma) * secondOrder + parameters.solver.gamma * firstOrder;

#ifndef NDEBUG
    if (fabs(tmp1 - tmp2) > 1.0e-12) {
      spdlog::error("{}, {}", tmp1, tmp2);
      throw std::runtime_error("Error in dw2dz");
    }
#endif

    return tmp2;
  }

  // adding the second derivatives
  // derivative of u
  inline RealType d2udx2(const RealType* const lv, const RealType* const lm) {
    const RealType u0 = lv[mapd(0, 0, 0, 0)]; // uij
    const RealType u1 = lv[mapd(1, 0, 0, 0)]; // ui+1
    const RealType u2 = lv[mapd(-1, 0, 0, 0)]; // ui-1
    const RealType dx_n = lm[mapd(0, 0, 0, 0)]; // dx-
    const RealType dx_p = lm[mapd(1, 0, 0, 0)]; // dx+
    return 2*( u1/(dx_p*(dx_p + dx_n)) - u0/(dx_n*dx_p) + u2/(dx_n*(dx_p + dx_n)));
  }

  inline RealType d2udy2(const RealType* const lv, const RealType* const lm) {
    const RealType u0 = lv[mapd(0, 0, 0, 0)]; // uij
    const RealType u1 = lv[mapd(0, 1, 0, 0)]; // uj+1
    const RealType u2 = lv[mapd(0, -1, 0, 0)]; // uj-1
    const RealType dy_n = lm[mapd(0, 0, 0, 1)]; // dy-
    const RealType dy_p = lm[mapd(0, 1, 0, 1)]; // dy+
    return 2*( u1/(dy_p*(dy_p + dy_n)) - u0/(dy_n*dy_p) + u2/(dy_n*(dy_p + dy_n)));
  }

  inline RealType d2udz2(const RealType* const lv, const RealType* const lm) {
    const RealType u0 = lv[mapd(0, 0, 0, 0)]; // uijk
    const RealType u1 = lv[mapd(0, 0, 1, 0)]; // uk+1
    const RealType u2 = lv[mapd(0, 0, -1, 0)]; //uk-1
    const RealType dz_n = lm[mapd(0, 0, 0, 2)]; // dz-
    const RealType dz_p = lm[mapd(0, 0, 1, 2)]; //dz+
    return 2*( u1/(dz_p*(dz_p + dz_n)) - u0/(dz_n*dz_p) + u2/(dz_n*(dz_p + dz_n)));
  }
  // derivative of v
  inline RealType d2vdx2(const RealType* const lv, const RealType* const lm) {
    const RealType v0 = lv[mapd(0, 0, 0, 1)]; // vij
    const RealType v1 = lv[mapd(1, 0, 0, 1)]; // vi+1
    const RealType v2 = lv[mapd(-1, 0, 0, 1)]; // vi-1
    const RealType dx_n = lm[mapd(0, 0, 0, 0)]; //dx-
    const RealType dx_p = lm[mapd(1, 0, 0, 0)]; //dx+
    return 2*( v1/(dx_p*(dx_p + dx_n)) - v0/(dx_n*dx_p) + v2/(dx_n*(dx_p + dx_n)));
  }

  inline RealType d2vdy2(const RealType* const lv, const RealType* const lm) {
    const RealType v0 = lv[mapd(0, 0, 0, 1)]; // vij
    const RealType v1 = lv[mapd(0, 1, 0, 1)]; // vj+1
    const RealType v2 = lv[mapd(0, -1, 0, 1)]; // vj-1
    const RealType dy_n = lm[mapd(0, 0, 0, 1)]; //dy-
    const RealType dy_p = lm[mapd(0, 1, 0, 1)]; //dy+
    return 2*( v1/(dy_p*(dy_p + dy_n)) - v0/(dy_n*dy_p) + v2/(dy_n*(dy_p + dy_n)));
  }
  
  inline RealType d2vdz2(const RealType* const lv, const RealType* const lm) {
    const RealType v0 = lv[mapd(0, 0, 0, 1)]; // vijk
    const RealType v1 = lv[mapd(0, 0, 1, 1)]; // vk+1
    const RealType v2 = lv[mapd(0, 0, -1, 1)]; // vk-1
    const RealType dz_n = lm[mapd(0, 0, 0, 2)]; //dz-
    const RealType dz_p = lm[mapd(0, 0, 1, 2)]; //dz+
    return 2*( v1/(dz_p*(dz_p + dz_n)) - v0/(dz_n*dz_p) + v2/(dz_n*(dz_p + dz_n)));
  }
  // derivative of w
  inline RealType d2wdx2(const RealType* const lv, const RealType* const lm) {
    const RealType w0 = lv[mapd(0, 0, 0, 2)]; // wijk
    const RealType w1 = lv[mapd(1, 0, 0, 2)]; // wi+1
    const RealType w2 = lv[mapd(-1, 0, 0, 2)]; // wi-1
    const RealType dx_n = lm[mapd(0, 0, 0, 0)]; // dx-
    const RealType dx_p = lm[mapd(1, 0, 0, 0)]; // dx+
    return 2*( w1/(dx_p*(dx_p + dx_n)) - w0/(dx_n*dx_p) + w2/(dx_n*(dx_p + dx_n)));
  }

  inline RealType d2wdy2(const RealType* const lv, const RealType* const lm) {
    const RealType w0 = lv[mapd(0, 0, 0, 2)]; // wijk
    const RealType w1 = lv[mapd(0, 1, 0, 2)]; // wj+1
    const RealType w2 = lv[mapd(0, -1, 0, 2)]; // wj-1
    const RealType dy_n = lm[mapd(0, 0, 0, 1)]; // dy-
    const RealType dy_p = lm[mapd(0, 1, 0, 1)]; // dy+
    return 2*( w1/(dy_p*(dy_p + dy_n)) - w0/(dy_n*dy_p) + w2/(dy_n*(dy_p + dy_n)));
  }
  
  inline RealType d2wdz2(const RealType* const lv, const RealType* const lm) {
    const RealType w0 = lv[mapd(0, 0, 0, 2)]; // wijk
    const RealType w1 = lv[mapd(0, 0, 1, 2)]; // wk+1
    const RealType w2 = lv[mapd(0, 0, -1, 2)]; // wk-1
    const RealType dz_n = lm[mapd(0, 0, 0, 2)]; // dz-
    const RealType dz_p = lm[mapd(0, 0, 1, 2)]; // dz+
    return 2*( w1/(dz_p*(dz_p + dz_n)) - w0/(dz_n*dz_p) + w2/(dz_n*(dz_p + dz_n)));
  }

  //New turbulence model derivatives for non uniform grids
 

 

  inline RealType computeF2D(
    const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    // TODO WS1
    return localVelocity[mapd(0, 0, 0, 0)]
        + dt * (1.0/parameters.flow.Re*(d2udx2(localVelocity, localMeshsize) + d2udy2(localVelocity, localMeshsize)) -du2dx(localVelocity, parameters, localMeshsize) - duvdy(localVelocity, parameters, localMeshsize) + parameters.environment.gx);
  }

  inline RealType computeG2D(
    const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    // TODO WS1
    return localVelocity[mapd(0, 0, 0, 1)]
        + dt * (1.0/parameters.flow.Re*(d2vdx2(localVelocity, localMeshsize) + d2vdy2(localVelocity, localMeshsize))- dv2dy(localVelocity, parameters, localMeshsize) - duvdx(localVelocity, parameters, localMeshsize) + parameters.environment.gy);
  }

  inline RealType computeF3D(
    const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    // TODO WS1
    return localVelocity[mapd(0, 0, 0, 0)]
        + dt * (1.0/parameters.flow.Re*(d2udx2(localVelocity, localMeshsize) + d2udy2(localVelocity, localMeshsize) + d2udz2(localVelocity, localMeshsize)) -du2dx(localVelocity, parameters, localMeshsize) - duvdy(localVelocity, parameters, localMeshsize)
            - duwdz(localVelocity, parameters, localMeshsize) + parameters.environment.gx);
  }

  inline RealType computeG3D(
    const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    // TODO WS1
    return localVelocity[mapd(0, 0, 0, 1)]
        + dt * (1.0/parameters.flow.Re*(d2vdx2(localVelocity, localMeshsize) + d2vdy2(localVelocity, localMeshsize) + d2vdz2(localVelocity, localMeshsize)) -dv2dy(localVelocity, parameters, localMeshsize) - duvdx(localVelocity, parameters, localMeshsize)
            - dvwdz(localVelocity, parameters, localMeshsize) + parameters.environment.gy);
  }

  inline RealType computeH3D(
    const RealType* const localVelocity, const RealType* const localMeshsize, const Parameters& parameters, RealType dt
  ) {
    // TODO WS1
    return localVelocity[mapd(0, 0, 0, 2)]
        + dt * (1.0/parameters.flow.Re*(d2wdx2(localVelocity, localMeshsize) + d2wdy2(localVelocity, localMeshsize) + d2wdz2(localVelocity, localMeshsize)) -dw2dz(localVelocity, parameters, localMeshsize) - duwdx(localVelocity, parameters, localMeshsize)
            - dvwdy(localVelocity, parameters, localMeshsize) + parameters.environment.gz);
  }

} // namespace Stencils
