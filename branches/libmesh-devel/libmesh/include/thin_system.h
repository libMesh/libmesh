// $Id: thin_system.h,v 1.3 2003-05-04 23:58:53 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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



#ifndef __thin_system_h__
#define __thin_system_h__

// C++ includes

// Local Includes
#include "equation_systems.h"
#include "system_base.h"


// Forward Declarations
class ThinSystem;


/**
 * This class provides a specific system class.  It aims
 * at thin systems, offering nothing more than just
 * the essentials needed to solve a system.  Note
 * that still additional vectors/matrices may be added,
 * as offered in the parent class \p SystemBase.
 */

// ------------------------------------------------------------
// ThinSystem class definition

class ThinSystem : public SystemBase
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  ThinSystem (EquationSystems<ThinSystem>& es,
	      const std::string&           name,
	      const unsigned int           number);

  /**
   * Destructor.
   */
  ~ThinSystem ();
  
  /**
   * Clear all the data structures associated with
   * the system. 
   */
  void clear ();

  /**
   * Initializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  void init ();

  /**
   * Reinitializes the member data fields associated with
   * the system, so that, e.g., \p assemble() may be used.
   */
  void reinit ();
 
  /**
   * Assemble the linear system.  Does not
   * actually call the solver.
   */
  void assemble ();
  
  /**
   * Assemble & solve the linear system.
   */
  std::pair<unsigned int, Real> solve ();
  
  /**
   * @returns \p "Thin".  Helps in identifying
   * the system type in an equation system file.
   */
  static const std::string system_type () { return "Thin"; }

  /**
   * @returns \p true when the other system contains
   * identical data, up to the given threshold.
   */
  bool compare (const ThinSystem& other_system, 
		const Real threshold,
		const bool verbose) const;

  /**
   * Register a user function to use in initializing the system.
   */
  void attach_init_function(void fptr(EquationSystems<ThinSystem>& es,
				      const std::string& name));
  
  /**
   * Register a user function to use in assembling the system
   * matrix and RHS.
   */
  void attach_assemble_function(void fptr(EquationSystems<ThinSystem>& es,
					  const std::string& name));
  
  /**
   * Function that initializes the system.
   */
  void (* init_system_fptr) (EquationSystems<ThinSystem>& es,
			     const std::string& name);
  
  /**
   * Function that assembles the system.
   */
  void (* assemble_fptr) (EquationSystems<ThinSystem>& es,
			  const std::string& name);


protected:

  /**
   * Reference to the \p EquationSystems<ThinSystem> data structure 
   * that handles us.   
   */
  EquationSystems<ThinSystem>& _equation_systems;

};



// ------------------------------------------------------------
// ThinSystem inline methods



#endif
