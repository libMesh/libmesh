// $Id: explicit_system.h,v 1.2 2004-02-29 18:28:09 benkirk Exp $

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



#ifndef __explicit_system_h__
#define __explicit_system_h__

// C++ includes

// Local Includes
#include "system.h"
#include "numeric_vector.h"


// Forward Declarations


/**
 * This class provides a specific system class.  It aims
 * at explicit systems, offering nothing more than just
 * the essentials needed to solve a system.  Note
 * that still additional vectors/matrices may be added,
 * as offered in the parent class \p System.
 */

// ------------------------------------------------------------
// ExplicitSystem class definition

class ExplicitSystem : public System
{
public:

  /**
   * Constructor.  Optionally initializes required
   * data structures.
   */
  ExplicitSystem (EquationSystems& es,
		  const std::string& name,
		  const unsigned int number);

  /**
   * Destructor.
   */
  ~ExplicitSystem ();
 
  /**
   * The type of system.
   */
  typedef ExplicitSystem sys_type;

  /**
   * @returns a clever pointer to the system.
   */
  sys_type & system () { return *this; }
  
  /**
   * Assembles & solves the linear system Ax=b. 
   */
  virtual void solve ();
 
  /**
   * @returns \p "Explicit".  Helps in identifying
   * the system type in an equation system file.
   */
  virtual std::string system_type () const { return "Explicit"; }

  /**
   * The system matrix.  Implicit systems are characterized by
   * the need to solve the linear system Ax=b.  This is the
   * right-hand-side vector b.
   */
  NumericVector<Number> * rhs;
  

protected:
};



// ------------------------------------------------------------
// ExplicitSystem inline methods


#endif
