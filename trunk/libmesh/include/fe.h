// $Id: fe.h,v 1.12 2003-02-20 04:59:58 benkirk Exp $

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



#ifndef __fe_h__
#define __fe_h__

// C++ includes

// Local includes
#include "fe_base.h"

// forward declarations




/**
 * A specific instatiation of the \p FEBase class. This
 * class is templated, and specific template instantiations
 * will result in different Finite Element families. Full specialization
 * of the template for specific dimensions(\p Dim) and families
 * (\p T) provide support for specific finite element types.
 * The use of templates allows for compile-time optimization,
 * however it requires that the specific finite element family
 * and dimension is also known at compile time.  If this is
 * too restricting for your application you can use the
 * \p FEBase::build() member to create abstract (but still optimized)
 * finite elements.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2003
 * \version $Revision: 1.12 $
 */

//-------------------------------------------------------------
// FE class definition
template <unsigned int Dim, FEFamily T>
class FE : public FEBase
{
public:
  
  /**
   * Constructor.
   */
  FE(const FEType& fet);
  
  /**
   * @returns the value of the \f$ i^{th} \f$ shape function at
   * point \p p.  This method allows you to specify the imension,
   * element type, and order directly.  This allows the method to
   * be static.
   */
  static Real shape(const ElemType t,
		    const Order o,
		    const unsigned int i,
		    const Point& p);

  /**
   * @returns the value of the \f$ i^{th} \f$ shape function at
   * point \p p.  This method allows you to specify the imension,
   * element type, and order directly.  This allows the method to
   * be static.
   */
  static Real shape(const Elem* elem,
		    const Order o,
		    const unsigned int i,
		    const Point& p);
  
  /**
   * @returns the \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape function at point \p p.  This method allows you to
   * specify the dimension, element type, and order directly.
   */
  static Real shape_deriv(const ElemType t,
			  const Order o,
			  const unsigned int i,
			  const unsigned int j,
			  const Point& p);

  /**
   * @returns the \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape functionelement type, and order directly.
   */
  static Real shape_deriv(const Elem* elem,
			  const Order o,
			  const unsigned int i,
			  const unsigned int j,
			  const Point& p);
  
  /**
   * Build the nodal soln from the element soln.
   * This is the solution that will be plotted.
   */
  static void nodal_soln(const Elem* elem, const Order o,
			 const std::vector<Number>& elem_soln,
			 std::vector<Number>& nodal_soln);

  /**
   * @returns the number of shape functions associated with
   * this finite element.
   */
  unsigned int n_shape_functions () const
  { return n_dofs (elem_type, fe_type.order); }

  /**
   * @returns the number of shape functions associated with
   * a finite element of type \p t and approximation order \p o.
   */
  static unsigned int n_shape_functions (const ElemType t,
					 const Order o)
  { return n_dofs (t,o); }

  /**
   * @returns the number of shape functions associated with this
   * finite element.
   */
  static unsigned int n_dofs(const ElemType t,
			     const Order o);

  /**
   * @returns the number of dofs at node \p n for a finite element
   * of type \p t and order \p o.
   */
  static unsigned int n_dofs_at_node(const ElemType t,
				     const Order o,
				     const unsigned int n);

  /**
   * @returns the number of dofs interior to the element,
   * not associated with any interior nodes.
   */
  static unsigned int n_dofs_per_elem(const ElemType t,
				      const Order o);
				       
  /**
   * @returns the location (on the reference element) of the
   * point \p p located in physical space.  This function requires
   * inverting the (possibly nonlinear) transformation map, so
   * it is not trivial. The optional parameter \p tolerance defines
   * how close is "good enough."  The map inversion iteration
   * computes the sequence \f$ \{ p_n \} \f$, and the iteration is
   * terminated when \f$ \|p - p_n\| < \mbox{\texttt{tolerance}} \f$
   */
  static Point inverse_map (const Elem* elem,
			    const Point& p,
			    const Real tolerance = 1.e-6);
  
  /**
   * This is at the core of this class. Use this for each
   * new element in the mesh.  Reinitializes all the physical 
   * element-dependent data based on the current element 
   * \p elem.
   */
  void reinit (const Elem* elem);
    
  /**
   * Reinitializes all the physical element-dependent data based on
   * the \p side of \p face.
   */
  void reinit (QBase* qside,
	       const Elem* elem,
	       const unsigned int side);

  /**
   * Provides the class with the quadrature rule, which provides the
   * locations (on a reference element) where the shape functions are
   * to be calculated.
   */
  void attach_quadrature_rule (QBase* q)
  { assert (q != NULL); qrule = q; return; }
  
  /**
   * @returns the total number of quadrature points.  Call this
   * to get an upper bound for the \p for loop in your simulation
   * for matrix assembly of the current element.
   */
  unsigned int n_quadrature_points () const;

    
#ifdef ENABLE_INFINITE_ELEMENTS


protected:

  /**
   * Initialize the data fields for the base of an
   * an infinite element.
   */
  void init_base_shape_functions(const QBase* q,
				 const Elem* e);


#endif


private:


  
  /** 
   * Update the various member data fields \p phi,
   * \p dphidxi, \p dphideta, \p dphidzeta, etc.
   * for the current element.
   */
  void init_shape_functions(const QBase* q,
			    const Elem* e);

  /** 
   * Same as before, but for a side.
   */  
  void init_shape_functions(const QBase* q,
			    const Elem* e,
			    const unsigned int s);

  /**
   * @returns the location (in physical space) of the point
   * \p p located on the reference element.
   */
  static Point map (const Elem* elem,
		    const Point& reference_point);
  
  /**
   * @returns d(xyz)/dxi (in physical space) of the point
   * \p p located on the reference element.
   */
  static Point map_xi (const Elem* elem,
		       const Point& reference_point);
  
  /**
   * @returns d(xyz)/deta (in physical space) of the point
   * \p p located on the reference element.
   */
  static Point map_eta (const Elem* elem,
			const Point& reference_point);

  /**
   * @returns d(xyz)/dzeta (in physical space) of the point
   * \p p located on the reference element.
   */
  static Point map_zeta (const Elem* elem,
			 const Point& reference_point);

};



/**
 * Hierarchic finite elements.  Still templated on the dimension,
 * \p Dim.  
 *
 * \author Benjamin S. Kirk
 * \date 2002-2003
 * \version $Revision: 1.12 $
 */

//-------------------------------------------------------------
// FEHierarchic class definition
template <unsigned int Dim>
class FEHierarchic : public FE<Dim,HIERARCHIC>
{
public:

  /**
   * Constructor. Creates a hierarchic finite element
   * to be used in dimension \p Dim.
   */
  FEHierarchic(const FEType& fet);
};



/**
 * Lagrange finite elements.  Still templated on the dimension,
 * \p Dim.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2003
 * \version $Revision: 1.12 $
 */

//-------------------------------------------------------------
// FELagrange class definition
template <unsigned int Dim>
class FELagrange : public FE<Dim,LAGRANGE>
{
public:

  /**
   * Constructor. Creates a Lagrange finite element
   * to be used in dimension \p Dim.
   */
  FELagrange(const FEType& fet);
};



/**
 * Monomial finite elements.  Still templated on the dimension,
 * \p Dim.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2003
 * \version $Revision: 1.12 $
 */

//-------------------------------------------------------------
// FEMonomial class definition
template <unsigned int Dim>
class FEMonomial : public FE<Dim,MONOMIAL>
{
public:

  /**
   * Constructor. Creates a monomial finite element
   * to be used in dimension \p Dim.
   */
  FEMonomial(const FEType& fet);
};



/**
 * Provide Typedefs for various element types.
 */
namespace FiniteElements
{
  /**
   * Convenient definition for a 1D
   * Hierarchic finite element.
   */
  typedef FE<1,HIERARCHIC> FEHierarchic1D;
  
  /**
   * Convenient definition for a 2D
   * Hierarchic finite element.
   */
  typedef FE<2,HIERARCHIC> FEHierarchic2D;
  
  /**
   * Convenient definition for a 3D
   * Hierarchic finite element.
   */
  typedef FE<3,HIERARCHIC> FEHierarchic3D;
  


  /**
   * Convenient definition for a 1D
   * Lagrange finite element.
   */
  typedef FE<1,LAGRANGE> FELagrange1D;

  /**
   * Convenient definition for a 2D
   * Lagrange finite element.
   */
  typedef FE<2,LAGRANGE> FELagrange2D;

  /**
   * Convenient definition for a 3D
   * Lagrange finite element.
   */
  typedef FE<3,LAGRANGE> FELagrange3D;


  
  /**
   * Convenient definition for a 1D
   * Monomial finite element.
   */
  typedef FE<1,MONOMIAL> FEMonomial1D;
  
  /**
   * Convenient definition for a 2D
   * Monomial finite element.
   */
  typedef FE<2,MONOMIAL> FEMonomial2D;

  /**
   * Convenient definition for a 3D
   * Monomial finite element.
   */
  typedef FE<3,MONOMIAL> FEMonomial3D;
};




// ------------------------------------------------------------
// FE class inline members
template <unsigned int Dim, FEFamily T>
inline
FE<Dim,T>::FE (const FEType& fet) :
  FEBase (Dim,fet)
{
  // Sanity check.  Make sure the
  // Family specified in the template instantiation
  // matches the one in the FEType object
  assert (T == fe_type.family);
};



// ------------------------------------------------------------
// FEHierarchic class inline members
template <unsigned int Dim>
inline
FEHierarchic<Dim>::FEHierarchic (const FEType& fet) :
  FE<Dim,HIERARCHIC> (fet)
{
};



// ------------------------------------------------------------
// FELagrange class inline members
template <unsigned int Dim>
inline
FELagrange<Dim>::FELagrange (const FEType& fet) :
  FE<Dim,LAGRANGE> (fet)
{
}



// ------------------------------------------------------------
// FEMonomial class inline members
template <unsigned int Dim>
inline
FEMonomial<Dim>::FEMonomial (const FEType& fet) :
  FE<Dim,MONOMIAL> (fet)
{
}

#endif
