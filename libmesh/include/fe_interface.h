// $Id: fe_interface.h,v 1.5 2003-01-24 17:24:38 jwpeterson Exp $

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



#ifndef __fe_interface_h__
#define __fe_interface_h__

// C++ includes
#include <vector>

// Local includes
#include "point.h"
#include "enum_elem_type.h"


// forward declarations
class Elem;
class FEType;


/**
 * This class provides an encapsulated access to all static
 * public member functions of finite element classes.
 * Using this class, one need not worry about the correct
 * finite element class, or modify all the various switch 
 * statements to account for new elements.  Instead, these 
 * different elements are caught here, and passed to the appropriate 
 * static class member.  By default, the FEBase class methods are used.
 *
 * Currently, only static access is provided.  A possible
 * extension might be to enable the finite element classes
 * access to the EquationSystems object, without
 * much hassle of redefining a lot.
 *
 * @author Daniel Dreyer, 2002
 */

// ------------------------------------------------------------
// FEInterface class definition

class FEInterface
{
private:

  /**
   * Empty constructor. Do not create an object of this type.
   */
  FEInterface();

public:
  
  /**
   * Destructor.
   */
  virtual ~FEInterface() {return;};

  /**
   * @returns the number of shape functions associated with this
   * finite element of type \p fe_t. 
   * Automatically decides which finite element class to use.
   */
  static unsigned int n_shape_functions(const unsigned int dim,
					const FEType& fe_t,
					const ElemType t);

  /**
   * @returns the number of shape functions associated with this
   * finite element.
   * Automatically decides which finite element class to use.
   */
  static unsigned int n_dofs(const unsigned int dim,
			     const FEType& fe_t,
			     const ElemType t);

  /**
   * @returns the number of dofs at node n for a finite element
   * of type \p fe_t.
   * Automatically decides which finite element class to use.
   */
  static unsigned int n_dofs_at_node(const unsigned int dim,
				     const FEType& fe_t,
				     const ElemType t,
				     const unsigned int n);

  /**
   * @returns the number of dofs interior to the element,
   * not associated with any interior nodes.
   * Automatically decides which finite element class to use.
   */
  static unsigned int n_dofs_per_elem(const unsigned int dim,
				      const FEType& fe_t,
				      const ElemType t);
				     
  
  /**
   * Build the nodal soln from the element soln.
   * This is the solution that will be plotted.
   * Automatically passes the request to the appropriate
   * finite element class member.
   */
  static void nodal_soln(const unsigned int dim,
			 const FEType& fe_t,
			 const Elem* elem,
			 const std::vector<number>& elem_soln,
			 std::vector<number>& nodal_soln);

  /**
   * @returns the location (on the reference element) of the
   * point \p p located in physical space.  This function requires
   * inverting the (probably nonlinear) transformation map, so
   * it is not trivial.
   */
  static Point inverse_map (const unsigned int dim,
			    const FEType& fe_t,
			    const Elem* elem,
			    const Point& p);

  /**
   * @returns true if the point p is located on the reference element
   * for element type t, false otherwise.  Since we are doing floating
   * point comparisons here the parameter \p eps can be specified to
   * indicate a tolerance.  For example, \f$ x \le 1 \f$  becomes
   * \f$ x \le 1 + \epsilon \f$. 
   */
  static bool on_reference_element(const Point& p,
				   const ElemType t,
				   const real eps=1.e-6);
  /**
   * @returns the value of the \f$ i^{th} \f$ shape function at
   * point \p p. This method allows you to specify the dimension,
   * element type, and order directly. Automatically passes the
   * request to the appropriate finite element class member.
   */
  static real shape(const unsigned int dim,
		    const FEType& fe_t,
		    const ElemType t,
		    const unsigned int i,
		    const Point& p);

  /**
   * @returns the value of the \f$ i^{th} \f$ shape function at
   * point \p p. This method allows you to specify the dimension,
   * element type, and order directly. Automatically passes the
   * request to the appropriate finite element class member.
   */
  static real shape(const unsigned int dim,
		    const FEType& fe_t,
		    const Elem* elem,
		    const unsigned int i,
		    const Point& p);
  
};





// ------------------------------------------------------------
// FEInterface class inline members



#endif
