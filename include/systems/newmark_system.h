// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA



#ifndef LIBMESH_NEWMARK_SYSTEM_H
#define LIBMESH_NEWMARK_SYSTEM_H

// Local Includes
#include "libmesh/linear_implicit_system.h"

// C++ includes

namespace libMesh
{

/**
 * This class contains a specific system class.
 * It provides an implicit time integration scheme
 * known as the Newmark method.
 *
 * In the algorithm implemented here the system is solved for
 * displacements.
 * Curently the Newmark scheme is implemented for constant
 * time step sizes only. This time step is stored in the
 * \p EquationSystems parameter named \p "Newmark \p time \p step".
 * For the case of constant time steps the matrix only has to be
 * assembled once, whereas the rhs has to be updated in each timestep.
 * Default values of the Newmark parameters \p alpha and \p delta
 * used for time integration are provided.
 * For details refer to the examples section.
 */
class NewmarkSystem : public LinearImplicitSystem
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  NewmarkSystem (EquationSystems & es,
                 const std::string & name,
                 const unsigned int number);

  /**
   * Destructor.
   */
  ~NewmarkSystem ();

  /**
   * The type of system.
   */
  typedef NewmarkSystem sys_type;

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear () libmesh_override;

  /**
   * Reinitializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  virtual void reinit () libmesh_override;

  /**
   * Assemble the linear system.  Does not
   * actually call the solver.
   */
  virtual void assemble () libmesh_override;

  /**
   * @returns \p "Newmark".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const libmesh_override { return "Newmark"; }


  //---------------------------------------------------------
  // These members are specific to the Newmark system
  //

  /**
   * Apply initial conditions.
   */
  void initial_conditions ();

  /**
   * Compute the global matrix by adding up scaled
   * mass damping and stiffness matrix.
   */
  void compute_matrix ();

  /**
   * Update the rhs.
   */
  void update_rhs ();

  /**
   * Update displacement, velocity and acceleration.
   */
  void update_u_v_a ();

  /**
   * Set the time step size and the newmark parameter alpha and
   * delta and calculate the constant parameters used for
   * time integration.
   */
  void set_newmark_parameters (const Real delta_T = _default_timestep,
                               const Real alpha   = _default_alpha,
                               const Real delta   = _default_delta);

private:

  /**
   * Constants used for the time integration.
   */
  Real _a_0;
  Real _a_1;
  Real _a_2;
  Real _a_3;
  Real _a_4;
  Real _a_5;
  Real _a_6;
  Real _a_7;

  /**
   * Returns true if the matrix assembly is finished.
   */
  bool _finished_assemble;

  /**
   * Default Newmark \p alpha
   */
  static const Real _default_alpha;

  /**
   * Default Newmark \p delta
   */
  static const Real _default_delta;

  /**
   * Default Newmark time step
   */
  static const Real _default_timestep;
};

} // namespace libMesh

#endif // LIBMESH_NEWMARK_SYSTEM_H
