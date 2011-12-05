// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef __adaptive_h__
#define __adaptive_h__

// C++ includes

// Local includes
#include "libmesh_common.h"
#include "linear.h"
#include "error_vector.h"
#include "kelly_error_estimator.h"
#include "mesh_refinement.h"
#include "mesh.h"

namespace libMesh
{

/**
 * This is a generic class that defines a adaptive to be used in a
 * simulation.  A user can define a adaptive by deriving from this
 * class and implementing certain functions.
 *
 * @author Benjamin S. Kirk, 2003-2004.
 */

// ------------------------------------------------------------
// Adaptive class definition

template <class T = Linear<> >
class Adaptive : public T
{
public:

  /**
   * Constructor. Requires a reference to a system to be solved.
   */
  Adaptive (EquationSystems& es);

  /**
   * Constructor.  Requires a referece to the \p EquationSystems object.
   */
  Adaptive (EquationSystems& es,
	    const std::string& name,
	    const unsigned int number);

  /**
   * Destructor.
   */
  ~Adaptive ();

  /**
   * Re-implement the solve member to do a fixed number of
   * linear solves
   */
  virtual void solve ();

  /**
   * @returns the current refinement step.
   */
  unsigned int refinement_step () const { return _refinement_step; }

  /**
   * @returns the number of refinement steps to take.
   */
  unsigned int n_refinement_steps () const { return _n_refinement_steps; }

  /**
   * Sets the number of refinement steps to take.
   */
  unsigned int & n_refinement_steps () { return _n_refinement_steps; }

  /**
   * @returns the maximum level for mesh refinement.
   */
  unsigned int max_refinement_level () const { return _max_refinement_level; }

  /**
   * Sets the maximum level for mesh refinement.
   */
  unsigned int & max_refinement_level () { return _max_refinement_level; }


protected:


  /**
   * Sets the current refinement step.
   */
  unsigned int & refinement_step () { return _refinement_step; }


private:


  /**
   * The current refinement step.
   */
  unsigned int _refinement_step;

  /**
   * The number of refinement steps to take.
   */
  unsigned int _n_refinement_steps;

  /**
   * The maximum allowable levels of refinement.
   */
  unsigned int _max_refinement_level;
};



// ------------------------------------------------------------
// Adaptive inline members
template <class T>
Adaptive<T>::Adaptive (EquationSystems& es) :
  T                    (es), // Call the base class constructor
  _refinement_step     (0),  // Solver parameters
  _n_refinement_steps  (1),
  _max_refinement_level(100)
{
libmesh_deprecated();
}



template <class T>
Adaptive<T>::Adaptive (EquationSystems& es,
		       const std::string& name,
		       const unsigned int number) :
  Adaptive (es),
  T        (es, name, number)
{
libmesh_deprecated();
}



template <class T>
Adaptive<T>::~Adaptive ()
{
}



template <class T>
void Adaptive<T>::solve ()
{
  // First solve the base system
  T::solve ();

  for (this->refinement_step()=0;
       this->refinement_step() < this->n_refinement_steps();
       this->refinement_step()++)
    {
      // Then estimate the error in the base system
      // and refine the mesh
      {
	ErrorVector error;

	KellyErrorEstimator error_estimator;

	error_estimator.estimate_error (this->system(), "incomp_ns", error);

	MeshRefinement mesh_refinement (this->mesh());

	mesh_refinement.flag_elements_by_error_fraction (error,
							 0.40,
							 0.40,
							 this->max_refinement_level());

	mesh_refinement.refine_and_coarsen_elements ();

	this->system().reinit ();
      }

      // Then re-solve the base system
      T::solve ();
    }
}

} // namespace libMesh


#endif // #define __adaptive_h__
