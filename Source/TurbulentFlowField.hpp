#pragma once

#include "DataStructures.hpp"
#include "Parameters.hpp"
#include "FlowField.hpp"

/** Flow field
 *
 * Class intended to contain the state of the domain.
 */
class TurbulentFlowField : public FlowField {
private:
  ScalarField h_;    //! Scalar field for the nearest wall distance
  ScalarField nu_t_; //! Scalar field for the turbulent viscosity
  ScalarField nu_star_; //! Scalar field for the total viscosity
  
public:
  /** Constructor for the 2D flow field
   *
   * Constructor for the flow field. Allocates all the fields and sets
   * the sizes. Currently, this contructor is only used for testing purposes.
   *
   * @param Nx Size of the fuild domain (non-ghost cells), in the X direction
   * @param Ny Size of the fuild domain (non-ghost cells), in the Y direction
   */
  TurbulentFlowField(int Nx, int Ny);

  /** Constructor for the 3D flow field
   *
   * Constructor for the flow field. Allocates all the fields and sets
   * the sizes. Currently, this contructor is only used for testing purposes.
   *
   * @param Nx Size of the fuild domain (non-ghost cells), in the X direction
   * @param Ny Size of the fuild domain (non-ghost cells), in the Y direction
   * @param Nz Size of the fuild domain (non-ghost cells), in the Z direction
   */
  TurbulentFlowField(int Nx, int Ny, int Nz);

  /** Constructs a field from parameters object
   *
   * Constructs a field from a parameters object, so that it dimensionality can be defined in
   * the configuration file.
   *
   * @param parameters Parameters object with geometric information
   */
  TurbulentFlowField(const Parameters& parameters);

  ~TurbulentFlowField() override = default;

  /** Obtain size in the X direction
   *
   * @return Number of cells in the X direction
   */
  int getNx() const;

  /** Obtain size in the Y direction
   *
   * @return Number of cells in the Y direction
   */
  int getNy() const;

  /** Obtain size in the Z direction
   *
   * @return Number of cells in the Z direction
   */
  int getNz() const;

  int getCellsX() const;
  int getCellsY() const;
  int getCellsZ() const;

  ScalarField& getPressure();
  VectorField& getVelocity();

  IntScalarField& getFlags();

  VectorField& getFGH();

  ScalarField& getRHS();

  // New getters for turbulence fields
  ScalarField& getWallDistance();
  ScalarField& getTurbulentViscosity();
  ScalarField& getTotalViscosity();
  
  void getPressureAndVelocity(RealType& pressure, RealType* const velocity, int i, int j);
  void getPressureAndVelocity(RealType& pressure, RealType* const velocity, int i, int j, int k);
};
