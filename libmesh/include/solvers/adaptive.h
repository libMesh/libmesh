// $Id: adaptive.h,v 1.1 2004-01-03 15:37:42 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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
#include "error_estimator.h"
#include "mesh_refinement.h"
#include "mesh.h"

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
};



// ------------------------------------------------------------
// Adaptive inline members
template <class T>
Adaptive<T>::Adaptive (EquationSystems& es) :
  T                   (es), // Call the base class constructor
  _refinement_step    (0),  // Solver parameters
  _n_refinement_steps (1)     
{
}



template <class T>
Adaptive<T>::Adaptive (EquationSystems& es,
		       const std::string& name,
		       const unsigned int number) :
  Adaptive (es),
  T        (es, name, number)
{
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
	
	ErrorEstimator error_estimator;
	
	error_estimator.flux_jump (this->system(), "Poisson", error);
	
	MeshRefinement mesh_refinement (this->mesh());
	
	mesh_refinement.flag_elements_by_error_fraction (error,
							 0.80,
							 0.07,
							 100);
	
	mesh_refinement.refine_and_coarsen_elements ();

	this->system().reinit ();
      }
  
      // Then re-solve the base system
      T::solve ();
    }
}


#endif // #define __adaptive_h__
