// $Id: fe_interface.h,v 1.18 2003-04-18 15:46:32 spetersen Exp $

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
#include <map>

// Local includes
#include "point.h"
#include "enum_elem_type.h"


// forward declarations
class Elem;
class FEType;
class FEComputeData;


/**
 * This class provides an encapsulated access to all @e static
 * public member functions of finite element classes.
 * Using this class, one need not worry about the correct
 * finite element class.
 *
 * @author Daniel Dreyer, 2002-2003
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
  virtual ~FEInterface() {return;}

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
   * finite element class member.  To indicate that
   * results from this specific implementation of
   * \p nodal_soln should not be used, the vector 
   * \p nodal_soln is returned empty.
   */
  static void nodal_soln(const unsigned int dim,
			 const FEType& fe_t,
			 const Elem* elem,
			 const std::vector<Number>& elem_soln,
			 std::vector<Number>& nodal_soln);

  /**
   * @returns the location (on the reference element) of the
   * point \p p located in physical space.  This function requires
   * inverting the (probably nonlinear) transformation map, so
   * it is not trivial. The optional parameter \p tolerance defines
   * how close is "good enough."  The map inversion iteration
   * computes the sequence \f$ \{ p_n \} \f$, and the iteration is
   * terminated when \f$ \|\|p - p_n\|\| < \mbox{\texttt{tolerance}} \f$
   */
  static Point inverse_map (const unsigned int dim,
			    const FEType& fe_t,
			    const Elem* elem,
			    const Point& p,
			    const Real tolerance = 1.e-4,
			    const bool secure = true);

  /**
   * @returns true if the point p is located on the reference element
   * for element type t, false otherwise.  Since we are doing floating
   * point comparisons here the parameter \p eps can be specified to
   * indicate a tolerance.  For example, \f$ x \le 1 \f$  becomes
   * \f$ x \le 1 + \epsilon \f$. 
   */
  static bool on_reference_element(const Point& p,
				   const ElemType t,
				   const Real eps=1.e-6);
  /**
   * @returns the value of the \f$ i^{th} \f$ shape function at
   * point \p p. This method allows you to specify the dimension,
   * element type, and order directly. Automatically passes the
   * request to the appropriate finite element class member.
   */
  static Real shape(const unsigned int dim,
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
  static Real shape(const unsigned int dim,
		    const FEType& fe_t,
		    const Elem* elem,
		    const unsigned int i,
		    const Point& p);

  /**
   * Lets the appropriate child of \p FEBase compute the requested 
   * data for the input specified in \p data, and returns the values
   * also through \p data.  See this as a generalization of \p shape().
   * Currently, with disabled infinite elements, returns a vector of
   * all shape functions of \p elem evaluated ap \p p.
   */
  static void compute_data(const unsigned int dim,
			   const FEType& fe_t,
			   const Elem* elem,
			   FEComputeData& data);


  /**
   * Computes the constraint matrix contributions (for
   * non-conforming adapted meshes) corresponding to 
   * variable number \p var_number.
   */
  static void compute_constraints (std::map<unsigned int,
				            std::map<unsigned int,
				                     float> > & constraints,
				   const unsigned int system_number,
				   const unsigned int variable_number,
				   const FEType& fe_t,
				   const Elem* elem);
private:


  /**
   * @returns true if \p et is an element to be processed by
   * class \p InfFE.  Otherwise, it returns false.
   * For compatibility with disabled infinite elements
   * it always returns false.
   */
  static bool is_InfFE_elem(const ElemType et);


#ifdef ENABLE_INFINITE_ELEMENTS

  // ------------------------------------------------------------
  /*
   * All these private members do the same as their public
   * counterparts, but for infinite elements. This dis-entangles
   * the calls to \p FE and \p InfFE.
   */

  static unsigned int ifem_n_shape_functions(const unsigned int dim,
					     const FEType& fe_t,
					     const ElemType t);

  static unsigned int ifem_n_dofs(const unsigned int dim,
				  const FEType& fe_t,
				  const ElemType t);

  static unsigned int ifem_n_dofs_at_node(const unsigned int dim,
					  const FEType& fe_t,
					  const ElemType t,
					  const unsigned int n);

  static unsigned int ifem_n_dofs_per_elem(const unsigned int dim,
					   const FEType& fe_t,
					   const ElemType t);
  
  static void ifem_nodal_soln(const unsigned int dim,
			      const FEType& fe_t,
			      const Elem* elem,
			      const std::vector<Number>& elem_soln,
			      std::vector<Number>& nodal_soln);

  static Point ifem_inverse_map (const unsigned int dim,
				 const FEType& fe_t,
				 const Elem* elem,
				 const Point& p,
				 const Real tolerance,
				 const bool secure = true);

  static bool ifem_on_reference_element(const Point& p,
					const ElemType t,
					const Real eps);

  static Real ifem_shape(const unsigned int dim,
			 const FEType& fe_t,
			 const ElemType t,
			 const unsigned int i,
			 const Point& p);

  static Real ifem_shape(const unsigned int dim,
			 const FEType& fe_t,
			 const Elem* elem,
			 const unsigned int i,
			 const Point& p);

  static void ifem_compute_data(const unsigned int dim,
				const FEType& fe_t,
				const Elem* elem,
				FEComputeData& data);

#endif


};





// ------------------------------------------------------------
// FEInterface class inline members
#ifndef ENABLE_INFINITE_ELEMENTS 

inline bool FEInterface::is_InfFE_elem(const ElemType)
{
  return false; 
}

#else

inline bool FEInterface::is_InfFE_elem(const ElemType et)
{

  switch (et)
    {
    case INFEDGE2:
    case INFQUAD4:
    case INFQUAD6:
    case INFHEX8:
    case INFHEX16:
    case INFHEX18:
    case INFPRISM6:
    case INFPRISM12:
      {
        return true;
      }
      
    default:
      { 
	return false;
      }
    }
}

#endif //ifndef ENABLE_INFINITE_ELEMENTS 






#endif // ifndef __fe_interface_h__
