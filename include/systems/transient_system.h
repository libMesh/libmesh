// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_TRANSIENT_SYSTEM_H
#define LIBMESH_TRANSIENT_SYSTEM_H

// Local Includes
#include "libmesh/system.h"
#include "libmesh/libmesh_config.h"

namespace libMesh
{

// Forward declarations
class LinearImplicitSystem;
class NonlinearImplicitSystem;
class ExplicitSystem;
#ifdef LIBMESH_HAVE_SLEPC
class EigenSystem;
#endif

/**
 * \brief Manages storage and variables for transient systems.
 *
 * This class is a specialized system for solving transient systems,
 * e.g., systems with a time dependency. It provides appropriate storage,
 * manages variables and ensures consistency for these systems.
 *
 * This template class adds the functionality to manage, in addition to the
 * solution vector, an old solution and an older solution vector. These
 * vectors are useful to simulate transient systems.
 *
 * \note Additional vectors/matrices can be added via parent class
 * interfaces.
 *
 * \author Benjamin S. Kirk
 * \date 2004
 * \brief Used for solving transient systems of equations.
 */
template <class Base>
class TransientSystem : public Base
{
public:

  /**
   * Constructor.
   */
  TransientSystem (EquationSystems & es,
                   const std::string & name,
                   const unsigned int number);

  /**
   * Special functions.
   * - This class has the same restrictions as the union of its
   *   potential base classes (currently LinearImplicitSystem,
   *   NonlinearImplicitSystem, ExplicitSystem, and System).
   * - Destructor is defaulted out-of-line.
   */
  TransientSystem (const TransientSystem &) = delete;
  TransientSystem & operator= (const TransientSystem &) = delete;
  TransientSystem (TransientSystem &&) = default;
  TransientSystem & operator= (TransientSystem &&) = delete;
  virtual ~TransientSystem ();

  /**
   * The type of system.
   */
  typedef TransientSystem<Base> sys_type;

  /**
   * \returns A reference to *this.
   */
  sys_type & system () { return *this; }

  /**
   * Clear all the data structures associated with
   * the system.
   */
  virtual void clear () override;

  /**
   * \returns \p "Transient" prepended to T::system_type().
   * Helps in identifying the system type in an equation
   * system file.
   */
  virtual std::string system_type () const override;


  //-----------------------------------------------------------------
  // access to the solution data fields

  /**
   * \returns The old solution (at the previous timestep)
   * for the specified global DOF.
   */
  Number old_solution (const dof_id_type global_dof_number) const;

  /**
   * \returns The older solution (two timesteps ago)
   * for the specified global DOF.
   */
  Number older_solution (const dof_id_type global_dof_number) const;

  /**
   * All the values I need to compute my contribution
   * to the simulation at hand.  Think of this as the
   * current solution with any ghost values needed from
   * other processors.
   */
  NumericVector<Number> * old_local_solution;

  /**
   * All the values I need to compute my contribution
   * to the simulation at hand.  Think of this as the
   * current solution with any ghost values needed from
   * other processors.
   */
  NumericVector<Number> * older_local_solution;

protected:

  /**
   * Re-update the local values when the mesh has changed.
   * This method takes the data updated by \p update() and
   * makes it up-to-date on the current mesh.
   */
  virtual void re_update () override;

private:

  /**
   * Helper function for (re-)adding old and older solution vectors.
   */
  virtual void add_old_vectors ();
};



// -----------------------------------------------------------
// Useful typedefs
typedef TransientSystem<LinearImplicitSystem> TransientImplicitSystem;
typedef TransientSystem<LinearImplicitSystem> TransientLinearImplicitSystem;
typedef TransientSystem<NonlinearImplicitSystem> TransientNonlinearImplicitSystem;
typedef TransientSystem<ExplicitSystem> TransientExplicitSystem;
typedef TransientSystem<System> TransientBaseSystem;
#ifdef LIBMESH_HAVE_SLEPC
typedef TransientSystem<EigenSystem> TransientEigenSystem;
#endif



// ------------------------------------------------------------
// TransientSystem inline methods
template <class Base>
inline
std::string TransientSystem<Base>::system_type () const
{
  std::string type = "Transient";
  type += Base::system_type ();

  return type;
}



} // namespace libMesh




#endif // LIBMESH_TRANSIENT_SYSTEM_H
