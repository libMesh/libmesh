// $Id: fe.h,v 1.25 2006-09-12 07:14:40 roystgnr Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
class DofConstraints;
class DofMap;
#ifdef ENABLE_INFINITE_ELEMENTS

template <unsigned int friend_Dim, FEFamily friend_T_radial, InfMapType friend_T_map>
class InfFE;

#endif



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
 * \version $Revision: 1.25 $
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
   *
   * On a p-refined element, \p o should be the total order of the element.
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
   *
   * On a p-refined element, \p o should be the base order of the element.
   */
  static Real shape(const Elem* elem,
		    const Order o,
		    const unsigned int i,
		    const Point& p);
  
  /**
   * @returns the \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape function at point \p p.  This method allows you to
   * specify the dimension, element type, and order directly.
   *
   * On a p-refined element, \p o should be the total order of the element.
   */
  static Real shape_deriv(const ElemType t,
			  const Order o,
			  const unsigned int i,
			  const unsigned int j,
			  const Point& p);

  /**
   * @returns the \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape function.  You must specify element type, and order directly.
   *
   * On a p-refined element, \p o should be the base order of the element.
   */
  static Real shape_deriv(const Elem* elem,
			  const Order o,
			  const unsigned int i,
			  const unsigned int j,
			  const Point& p);

  /**
   * @returns the second \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape function at the point \p p.  Note that cross-derivatives are
   * also possible, i.e.
   * j = 0 ==> d^2 phi / dxi^2
   * j = 1 ==> d^2 phi / dxi deta
   * j = 2 ==> d^2 phi / deta^2
   * j = 3 ==> d^2 phi / dxi dzeta
   * j = 4 ==> d^2 phi / deta dzeta
   * j = 5 ==> d^2 phi / dzeta^2
   *
   * Note:  Computing second derivatives is not currently supported
   * for all element types: C1 (Clough and Hermite), Lagrange,
   * Hierarchic, and Monomial are supported.
   * All other element types return an error when asked for second derivatives.
   *
   * On a p-refined element, \p o should be the total order of the element.
   */
  static Real shape_second_deriv(const ElemType t,
				 const Order o,
				 const unsigned int i,
				 const unsigned int j,
				 const Point& p);

  /**
   * @returns the second \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape function at the point \p p.  Note that cross-derivatives are
   * also possible, i.e.
   * j = 0 ==> d^2 phi / dxi^2
   * j = 1 ==> d^2 phi / dxi deta
   * j = 2 ==> d^2 phi / deta^2
   * j = 3 ==> d^2 phi / dxi dzeta
   * j = 4 ==> d^2 phi / deta dzeta
   * j = 5 ==> d^2 phi / dzeta^2
   *
   * Note:  Computing second derivatives is not currently supported
   * for all element types: C1 (Clough and Hermite), Lagrange,
   * Hierarchic, and Monomial are supported.
   * All other element types return an error when asked for second derivatives.
   *
   * On a p-refined element, \p o should be the base order of the element.
   */
  static Real shape_second_deriv(const Elem* elem,
				 const Order o,
				 const unsigned int i,
				 const unsigned int j,
				 const Point& p);
  
  /**
   * Build the nodal soln from the element soln.
   * This is the solution that will be plotted.
   *
   * On a p-refined element, \p o should be the base order of the element.
   */
  static void nodal_soln(const Elem* elem, const Order o,
			 const std::vector<Number>& elem_soln,
			 std::vector<Number>& nodal_soln);

  /**
   * @returns the number of shape functions associated with
   * this finite element.
   */
  virtual unsigned int n_shape_functions () const;

  /**
   * @returns the number of shape functions associated with
   * a finite element of type \p t and approximation order \p o.
   *
   * On a p-refined element, \p o should be the total order of the element.
   */
  static unsigned int n_shape_functions (const ElemType t,
					 const Order o)
  { return FE<Dim,T>::n_dofs (t,o); }

  /**
   * @returns the number of shape functions associated with this
   * finite element.
   *
   * On a p-refined element, \p o should be the total order of the element.
   */
  static unsigned int n_dofs(const ElemType t,
			     const Order o);

  /**
   * @returns the number of dofs at node \p n for a finite element
   * of type \p t and order \p o.
   *
   * On a p-refined element, \p o should be the total order of the element.
   */
  static unsigned int n_dofs_at_node(const ElemType t,
				     const Order o,
				     const unsigned int n);

  /**
   * @returns the number of dofs interior to the element,
   * not associated with any interior nodes.
   *
   * On a p-refined element, \p o should be the total order of the element.
   */
  static unsigned int n_dofs_per_elem(const ElemType t,
				      const Order o);

  /**
   * @returns the continuity level of the finite element.
   */
  virtual FEContinuity get_continuity() const;
				       
  /**
   * @returns true if the finite element's higher order shape functions are
   * hierarchic
   */
  virtual bool is_hierarchic() const;
				       
  /**
   * Fills the vector di with the local degree of freedom indices
   * associated with side \p s of element \p elem
   *
   * On a p-refined element, \p o should be the base order of the element.
   */
  static void dofs_on_side(const Elem* const elem,
		           const Order o,
		           unsigned int s,
		           std::vector<unsigned int>& di);
  /**
   * Fills the vector di with the local degree of freedom indices
   * associated with edge \p e of element \p elem
   *
   * On a p-refined element, \p o should be the base order of the element.
   */
  static void dofs_on_edge(const Elem* const elem,
		           const Order o,
		           unsigned int e,
		           std::vector<unsigned int>& di);

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
			    const Real tolerance = TOLERANCE,
			    const bool secure = true);

  /**
   * Takes a number points in physical space (in the \p
   * physical_points vector) and finds their location on the reference
   * element for the input element \p elem.  The values on the
   * reference element are returned in the vector \p
   * reference_points. The optional parameter \p tolerance defines how
   * close is "good enough."  The map inversion iteration computes the
   * sequence \f$ \{ p_n \} \f$, and the iteration is terminated when
   * \f$ \|p - p_n\| < \mbox{\texttt{tolerance}} \f$
   */
  static void inverse_map (const Elem* elem,
			   const std::vector<Point>& physical_points,
			   std::vector<Point>&       reference_points,
			   const Real tolerance = TOLERANCE,
			   const bool secure = true);
  
  /**
   * This is at the core of this class. Use this for each
   * new element in the mesh.  Reinitializes all the physical 
   * element-dependent data based on the current element 
   * \p elem.  By default the shape functions and associated
   * data are computed at the quadrature points specified
   * by the quadrature rule \p qrule, but may be any points
   * specified on the reference element specified in the optional
   * argument \p pts.
   */
  virtual void reinit (const Elem* elem,
		       const std::vector<Point>* const pts = NULL);
    
  /**
   * Reinitializes all the physical element-dependent data based on
   * the \p side of \p face.  The \p tolerance paremeter is passed to
   * the involved call to \p inverse_map().
   */
  virtual void reinit (const Elem* elem,
		       const unsigned int side,
		       const Real tolerance = TOLERANCE);

  /**
   * Reinitializes all the physical element-dependent data based on
   * the \p edge.  The \p tolerance paremeter is passed to the
   * involved call to \p inverse_map().
   */
  virtual void edge_reinit (const Elem* elem,
		            const unsigned int edge,
			    const Real tolerance = TOLERANCE);

  /**
   * Provides the class with the quadrature rule, which provides the
   * locations (on a reference element) where the shape functions are
   * to be calculated.
   */
  virtual void attach_quadrature_rule (QBase* q);
  
  /**
   * @returns the total number of quadrature points.  Call this
   * to get an upper bound for the \p for loop in your simulation
   * for matrix assembly of the current element.
   */
  virtual unsigned int n_quadrature_points () const;

  /**
   * Computes the constraint matrix contributions (for
   * non-conforming adapted meshes) corresponding to 
   * variable number \p var_number, using element-specific
   * optimizations if possible.
   */
  static void compute_constraints (DofConstraints &constraints,
				   DofMap &dof_map,
				   const unsigned int variable_number,
				   const Elem* elem);
  
  /**
   * Computes the constraint matrix contributions (for
   * non-conforming adapted meshes) corresponding to 
   * variable number \p var_number, using generic
   * projections.
   */
  static void compute_proj_constraints (DofConstraints &constraints,
				        DofMap &dof_map,
				        const unsigned int variable_number,
				        const Elem* elem);
  
#ifdef ENABLE_INFINITE_ELEMENTS

  /**
   * Initialize the data fields for the base of an
   * an infinite element.
   */
  void init_base_shape_functions(const std::vector<Point>& qp,
				 const Elem* e);

#endif

  /**
   * @returns \p true when the shape functions (for
   * this \p FEFamily) depend on the particular 
   * element, and therefore needs to be re-initialized
   * for each new element.  \p false otherwise.
   */
  virtual bool shapes_need_reinit() const;
  
  /** 
   * Update the various member data fields \p phi,
   * \p dphidxi, \p dphideta, \p dphidzeta, etc.
   * for the current element.  These data will be computed
   * at the points \p qp, which are generally (but need not be)
   * the quadrature points.
   */
  virtual void init_shape_functions(const std::vector<Point>& qp,
				    const Elem* e);

  /** 
   * Same as before, but for a side. This is used for boundary
   * integration.
   */  
  void init_face_shape_functions(const std::vector<Point>& qp,
				 const Elem* side);

  /** 
   * Same as before, but for an edge. This is used for some projection
   * operators.
   */  
  void init_edge_shape_functions(const std::vector<Point>& qp,
				 const Elem* edge);

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

#ifdef ENABLE_INFINITE_ELEMENTS
  /**
   * make InfFE classes friends, so that these may access
   * the private \p map, map_xyz methods
   */
  template <unsigned int friend_Dim, FEFamily friend_T_radial, InfMapType friend_T_map>
  friend class InfFE;
#endif

};



/**
 * Clough-Tocher finite elements.  Still templated on the dimension,
 * \p Dim.  
 *
 * \author Roy Stogner
 * \date 2004
 * \version $Revision: 1.25 $
 */

//-------------------------------------------------------------
// FEHierarchic class definition
template <unsigned int Dim>
class FEClough : public FE<Dim,CLOUGH>
{
public:

  /**
   * Constructor. Creates a hierarchic finite element
   * to be used in dimension \p Dim.
   */
  FEClough(const FEType& fet);
};



/**
 * Hermite finite elements.  Still templated on the dimension,
 * \p Dim.  
 *
 * \author Roy Stogner
 * \date 2005
 * \version $Revision: 1.25 $
 */

//-------------------------------------------------------------
// FEHierarchic class definition
template <unsigned int Dim>
class FEHermite : public FE<Dim,HERMITE>
{
public:

  /**
   * Constructor. Creates a hierarchic finite element
   * to be used in dimension \p Dim.
   */
  FEHermite(const FEType& fet);

  /**
   * 1D hermite functions on unit interval
   */
  static Real hermite_raw_shape_second_deriv(const unsigned int basis_num,
					     const Real xi);
  static Real hermite_raw_shape_deriv(const unsigned int basis_num,
				      const Real xi);
  static Real hermite_raw_shape(const unsigned int basis_num,
				const Real xi);
};



/**
 * Hierarchic finite elements.  Still templated on the dimension,
 * \p Dim.  
 *
 * \author Benjamin S. Kirk
 * \date 2002-2003
 * \version $Revision: 1.25 $
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
 * \version $Revision: 1.25 $
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
 * \version $Revision: 1.25 $
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
 * XYZ finite elements.  These require specialization
 * because the shape functions are defined in terms of
 * physical XYZ coordinates rather than local coordinates.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2003
 * \version $Revision: 1.25 $
 */

//-------------------------------------------------------------
// FEXYZ class definition
template <unsigned int Dim>
class FEXYZ : public FE<Dim,XYZ>
{
public:
  
  /**
   * Constructor. Creates a monomial finite element
   * to be used in dimension \p Dim.
   */
  FEXYZ(const FEType& fet);
    
  /**
   * Explicitly call base class method.  This prevents some
   * compilers being confused by partially overriding this virtual function.
   */
  virtual void reinit (const Elem* elem,
		       const std::vector<Point>* const pts = NULL)
  { FE<Dim,XYZ>::reinit (elem, pts); }

  /**
   * Reinitializes all the physical element-dependent data based on
   * the \p side of \p face.
   */
  virtual void reinit (const Elem* elem,
		       const unsigned int side,
		       const Real tolerance = TOLERANCE);


protected:

  /** 
   * Update the various member data fields \p phi,
   * \p dphidxi, \p dphideta, \p dphidzeta, etc.
   * for the current element.  These data will be computed
   * at the points \p qp, which are generally (but need not be)
   * the quadrature points.
   */
  virtual void init_shape_functions(const std::vector<Point>& qp,
				    const Elem* e);

  /** 
   * After having updated the jacobian and the transformation
   * from local to global coordinates in \p FEBase::compute_map(),
   * the first derivatives of the shape functions are 
   * transformed to global coordinates, giving \p dphi,
   * \p dphidx, \p dphidy, and \p dphidz. This method
   * should rarely be re-defined in derived classes, but
   * still should be usable for children. Therefore, keep
   * it protected.
   */
  virtual void compute_shape_functions(const Elem*);

  /**
   * Compute the map & shape functions for this face.
   */
  void compute_face_values (const Elem* elem,
			    const Elem* side);
};



/**
 * Provide Typedefs for various element types.
 */
namespace FiniteElements
{
  /**
   * Convenient definition for a 2D
   * Clough-Tocher finite element.
   */
  typedef FEClough<2> FEClough2D;
  
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
}




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
}



// ------------------------------------------------------------
// FEClough class inline members
template <unsigned int Dim>
inline
FEClough<Dim>::FEClough (const FEType& fet) :
  FE<Dim,CLOUGH> (fet)
{
}



// ------------------------------------------------------------
// FEHermite class inline members
template <unsigned int Dim>
inline
FEHermite<Dim>::FEHermite (const FEType& fet) :
  FE<Dim,HERMITE> (fet)
{
}



// ------------------------------------------------------------
// FEHierarchic class inline members
template <unsigned int Dim>
inline
FEHierarchic<Dim>::FEHierarchic (const FEType& fet) :
  FE<Dim,HIERARCHIC> (fet)
{
}



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




// ------------------------------------------------------------
// FELagrange class inline members
template <unsigned int Dim>
inline
FEXYZ<Dim>::FEXYZ (const FEType& fet) :
  FE<Dim,XYZ> (fet)
{
}


#endif
