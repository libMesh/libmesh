// $Id: fe.h,v 1.3 2003-01-20 17:06:10 jwpeterson Exp $

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
#include <vector>

// Local includes
#include "reference_counted_object.h"
#include "point.h"
#include "enum_elem_type.h"
#include "fe_type.h"
#include "auto_ptr.h"


// forward declarations
class QBase;
class MeshBase;
class Elem;
class FEBase;



/**
 * This class forms the foundation from which generic finite
 * elements may be derived. The current implementation offers
 * a wide variety of commonly used finite element concepts.
 * When this class is not enough, everything may be overloaded,
 * but it generally suffices to define own public static shape functions,
 * to overload FEBase::reinit(), and FEBase::init_shape_functions(). 
 * To use different coordinate transformations, re-define
 * FEBase::compute_map(), FEBase::on_reference_element(), 
 * FEBase::inverse_map() etc. 
 * For weird dof counts in derived classes, re-define the
 * n_dofs_...() family and, if necessary, FEBase::nodal_soln().
 * Note that all of these member functions are not virtual,
 * so the decision which child class of FEBase to use, is
 * up to the user by identifying the appropriate element type,
 * see InfFE::is_inf_elem().
 *
 * All interaction of this and derived classes with other classes, 
 * like DofMap, are handled through the interface class FEInterface. 
 * When the static member functions like FEBase::n_dofs() etc are 
 * not sufficient to represent your derived class, add calls to your 
 * version in FEInterface, better do not modify the methods from FEBase.
 * Within the FEBase class, things should remain unchanged. 
 * Through static member functions this class implements 
 * Lagrange shape functions, monomials, and hierarchic shape functions.
 * These may be used, through the public methods FEBase::shape() and
 * FEBase::shape_deriv(), as building blocks for shape functions in
 * derived classes.
 *
 * @author Benjamin S. Kirk, 2002
 */

// ------------------------------------------------------------
// FEBase class definition
class FEBase : public ReferenceCountedObject<FEBase>
{
protected:

  /**
   * Constructor.  Optionally initializes required data
   * structures.  Protected so that this base class
   * cannot be explicitly instantiated.
   */
  FEBase (const unsigned int dim,
	  const FEType& fet);
  
public:
  
  /**
   * Destructor.
   */
  virtual ~FEBase();

  /**
   * Builds a specific finite element type.  A \p AutoPtr<Elem> is
   * returned to prevent a memory leak. This way the user need not
   * remember to delete the object.
   */
  static AutoPtr<FEBase> build (const unsigned int dim,
				const FEType& type); 

  /**
   * This is at the core of this class. Use this for each
   * new element in the mesh.  Reinitializes all the physical 
   * element-dependent data based on the current element 
   * \p elem.
   */
  virtual void reinit (const Elem* elem) = 0;
    
  /**
   * Reinitializes all the physical element-dependent data based on
   * the \p side of \p face.
   */
  virtual void reinit (QBase* qside,
		       const Elem* elem,
		       const unsigned int side) = 0;
  
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
   * @returns the \p xyz spatial locations of the quadrature
   * points on the element.
   */    
  const std::vector<Point>& get_xyz() const
  { return xyz; };
  
  /**
   * @returns the shape function values at the quadrature points
   * on the element.
   */    
  const std::vector<std::vector<real> >& get_phi() const
  { return phi; };
  
  /**
   * @returns the element Jacobian times the quadrature weight for
   * each quadrature point.
   */    
  const std::vector<real>& get_JxW() const
  { return JxW; };

  /**
   * @returns the shape function derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<Point> >& get_dphi() const
  { return dphi; };
  
  /**
   * @returns the shape function x-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<real> >& get_dphidx() const
  { return dphidx; };
  
  /**
   * @returns the shape function y-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<real> >& get_dphidy() const
  { return dphidy; };
  
  /**
   * @returns the shape function z-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<real> >& get_dphidz() const
  { return dphidz; };
  
  /**
   * @returns the tangent vectors for face integration.
   */
  const std::vector<std::vector<Point> >& get_tangents() const
  { return tangents; };
  
  /**
   * @returns the normal vectors for face integration.
   */
  const std::vector<Point>& get_normals() const
  { return normals; };
  
  /**
   * Provides the class with the quadrature rule, which provides the
   * locations (on a reference element) where the shape functions are
   * to be calculated.
   */
  void attach_quadrature_rule (QBase* q)
  { assert (q != NULL); qrule = q; return; };
  
  /**
   * @returns the element type that the current shape functions
   * have been calculated for.  Useful in determining when shape
   * functions must be recomputed.
   */
  ElemType get_type()  const { return elem_type; };

  /**
   * @returns the approximation order of the finite element.
   */
  Order get_order()  const { return fe_type.order; };

  /**
   * @returns the finite element family of this element.
   */
  FEFamily get_family()  const { return fe_type.family; };

  /**
   * Prints the Jacobian times the weight for each quadrature point.
   */ 
  void print_JxW() const; 
  
  /**
   * Prints the value of each shape function at each quadrature point.
   */ 
  void print_phi() const; 
  
  /**
   * Prints the value of each shape function's derivative
   * at each quadrature point.
   */ 
  void print_dphi() const;
  
  /**
   * Prints the spatial location of each quadrature point
   * (on the physical element).
   */ 
  void print_xyz() const;

  /**
   * Prints all the relevant information about the current element.
   */
  void print_info() const;
  
  
  
protected:

  /**
   * Compute the jacobian and some other additional
   * data fields.
   */
  void compute_map(const QBase* q,
		   const Elem* e);
  
  /** 
   * Same as before, but for a side.
   */  
  void compute_map(const QBase* q,
		   const Elem* e,
		   const unsigned int s);

  /** 
   * After having updated the jacobian and the transformation
   * from local to global coordinates in FEBase::compute_map(),
   * the first derivatives of the shape functions are 
   * transformed to global coordinates, giving \p dphi,
   * \p dphidx, \p dphidy, and \p dphidz. This method
   * should barely be re-defined in derived classes, but
   * still should be usable for children. Therefore, keep
   * it protected.
   */
  void compute_shape_functions(const QBase* q);
  

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the x value of the pth entry of the dxzydxi_map.
   */
  real dxdxi_map(const unsigned int p) const   { return dxyzdxi_map[p](0); };

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the y value of the pth entry of the dxzydxi_map.
   */
  real dydxi_map(const unsigned int p) const   { return dxyzdxi_map[p](1); };

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the z value of the pth entry of the dxzydxi_map.
   */
  real dzdxi_map(const unsigned int p) const   { return dxyzdxi_map[p](2); };

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the x value of the pth entry of the dxzydeta_map.
   */
  real dxdeta_map(const unsigned int p) const  { return dxyzdeta_map[p](0); };

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the y value of the pth entry of the dxzydeta_map.
   */
  real dydeta_map(const unsigned int p) const  { return dxyzdeta_map[p](1); }; 

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the z value of the pth entry of the dxzydeta_map.
   */
  real dzdeta_map(const unsigned int p) const  { return dxyzdeta_map[p](2); };

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the x value of the pth entry of the dxzydzeta_map.
   */
  real dxdzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](0); };

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the y value of the pth entry of the dxzydzeta_map.
   */
  real dydzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](1); };

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the z value of the pth entry of the dxzydzeta_map.
   */
  real dzdzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](2); };




  
  /**
   * The dimensionality of the object
   */
  const unsigned int dim;

  /**
   * The spatial locations of the quadrature points
   */
  std::vector<Point> xyz;


  
  /**
   * Vector of parital derivatives: 
   * d(x)/d(xi), d(y)/d(xi), d(z)/d(xi) 
   */
  std::vector<Point> dxyzdxi_map;

  /**
   * Vector of parital derivatives: 
   * d(x)/d(eta), d(y)/d(eta), d(z)/d(eta)
   */
  std::vector<Point> dxyzdeta_map;

  /**
   * Vector of parital derivatives: 
   * d(x)/d(zeta), d(y)/d(zeta), d(z)/d(zeta)
   */
  std::vector<Point> dxyzdzeta_map;
  

  /**
   * Map for partial derivatives:
   * d(xi)/d(x). Needed for the Jacobian.
   */
  std::vector<real>  dxidx_map;

  /**
   * Map for partial derivatives:
   * d(xi)/d(y). Needed for the Jacobian.
   */
  std::vector<real>  dxidy_map;

  /**
   * Map for partial derivatives:
   * d(xi)/d(z). Needed for the Jacobian.
   */
  std::vector<real>  dxidz_map;



  
  /**
   * Map for partial derivatives:
   * d(eta)/d(x). Needed for the Jacobian.
   */
  std::vector<real>  detadx_map;

  /**
   * Map for partial derivatives:
   * d(eta)/d(y). Needed for the Jacobian.
   */
  std::vector<real>  detady_map;

  /**
   * Map for partial derivatives:
   * d(eta)/d(z). Needed for the Jacobian.
   */
  std::vector<real>  detadz_map;




  
  /**
   * Map for partial derivatives:
   * d(zeta)/d(x). Needed for the Jacobian.
   */
  std::vector<real>  dzetadx_map;

  /**
   * Map for partial derivatives:
   * d(zeta)/d(y). Needed for the Jacobian.
   */
  std::vector<real>  dzetady_map;

  /**
   * Map for partial derivatives:
   * d(zeta)/d(z). Needed for the Jacobian.
   */
  std::vector<real>  dzetadz_map;


  
  /**
   * Shape function values.
   */
  std::vector<std::vector<real> >   phi;

  /**
   * Shape function derivative values.
   */
  std::vector<std::vector<Point> >  dphi;

  /**
   * Shape function derivatives in the xi direction.
   */
  std::vector<std::vector<real> >   dphidxi;

  /**
   * Shape function derivatives in the eta direction.
   */
  std::vector<std::vector<real> >   dphideta;
  
  /**
   * Shape function derivatives in the zeta direction.
   */
  std::vector<std::vector<real> >   dphidzeta;

  /**
   * Shape function derivatives in the x direction.
   */
  std::vector<std::vector<real> >   dphidx;

  /**
   * Shape function derivatives in the y direction.
   */
  std::vector<std::vector<real> >   dphidy;

  /**
   * Shape function derivatives in the z direction.
   */
  std::vector<std::vector<real> >   dphidz;




  
  /**
   * Map for the shape function phi.
   */
  std::vector<std::vector<real> >   phi_map;

  /**
   * Map for the derivative, d(phi)/d(xi).
   */
  std::vector<std::vector<real> >   dphidxi_map;

  /**
   * Map for the derivative, d(phi)/d(eta).
   */
  std::vector<std::vector<real> >   dphideta_map;

  /**
   * Map for the derivative, d(phi)/d(zeta).
   */
  std::vector<std::vector<real> >   dphidzeta_map;




  
  /**
   * Map for the side shape functions, psi. 
   */
  std::vector<std::vector<real> >   psi_map;

  /**
   * Map for the derivative of the side functions,
   * d(psi)/d(xi).
   */
  std::vector<std::vector<real> >   dpsidxi_map;

  /**
   * Map for the derivative of the side function,
   * d(psi)/d(eta).
   */
  std::vector<std::vector<real> >   dpsideta_map;



  
  /**
   * Tangent vectors on boundary at quadrature points.
   */
  std::vector<std::vector<Point> >  tangents;

  /**
   * Normal vectors on boundary at quadrature points
   */
  std::vector<Point>                normals;

  /**
   * Jacobian values at quadrature points
   */
  std::vector<real>                 jac;
  
  /**
   * Jacobian*Weight values at quadrature points
   */
  std::vector<real>                 JxW;

  /**
   * The finite element type for this object.  Note that this
   * should be constant for the object.
   */
  const FEType fe_type;
  
  /**
   * The element type the current data structures are
   * set up for.
   */
  ElemType elem_type;

  /**
   * A pointer to the quadrature rule employed
   */
  QBase* qrule;



private:


  /**
   * Make the mesh class a friend
   */
  friend class MeshBase;
};



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
 * \version $Revision: 1.3 $
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
  static real shape(const ElemType t,
		    const Order o,
		    const unsigned int i,
		    const Point& p);

  /**
   * @returns the value of the \f$ i^{th} \f$ shape function at
   * point \p p.  This method allows you to specify the imension,
   * element type, and order directly.  This allows the method to
   * be static.
   */
  static real shape(const Elem* elem,
		    const Order o,
		    const unsigned int i,
		    const Point& p);
  
  /**
   * @returns the \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape function at point \p p.  This method allows you to
   * specify the dimension, element type, and order directly.
   */
  static real shape_deriv(const ElemType t,
			  const Order o,
			  const unsigned int i,
			  const unsigned int j,
			  const Point& p);

  /**
   * @returns the \f$ j^{th} \f$ derivative of the \f$ i^{th} \f$
   * shape functionelement type, and order directly.
   */
  static real shape_deriv(const Elem* elem,
			  const Order o,
			  const unsigned int i,
			  const unsigned int j,
			  const Point& p);
  
  /**
   * Build the nodal soln from the element soln.
   * This is the solution that will be plotted.
   */
  static void nodal_soln(const Elem* elem, const Order o,
			 const std::vector<number>& elem_soln,
			 std::vector<number>& nodal_soln);

  /**
   * @returns the number of shape functions associated with
   * this finite element.
   */
  unsigned int n_shape_functions () const
  { return n_dofs (elem_type, fe_type.order); };

  /**
   * @returns the number of shape functions associated with
   * a finite element of type \p t and approximation order \p o.
   */
  static unsigned int n_shape_functions (const ElemType t,
					 const Order o)
  { return n_dofs (t,o); };

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
   * it is not trivial.
   */
  static Point inverse_map (const Elem* elem,
			    const Point& p);
  
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

  /**
   * Make the \p FEBase class a friend so that its
   * \p FEBase::build() member will work.
   */
  friend class FEBase;
};



/**
 * Hierarchic finite elements.  Still templated on the dimension,
 * \p Dim.  
 *
 * \author Benjamin S. Kirk
 * \date 2002-2003
 * \version $Revision: 1.3 $
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
 * \version $Revision: 1.3 $
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
 * \version $Revision: 1.3 $
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
// FEBase class inline members
inline
FEBase::FEBase(const unsigned int d,
	       const FEType& fet) :
  dim(d),
  fe_type(fet),
  elem_type(INVALID_ELEM),
  qrule(NULL)
{
};



inline
FEBase::~FEBase()
{
};



inline
void FEBase::print_JxW() const
{
  for (unsigned int i=0; i<JxW.size(); ++i) std::cout << JxW[i] << std::endl;
};



inline
void FEBase::print_phi() const
{
  for (unsigned int i=0; i<phi.size(); ++i)
    {
      for (unsigned int j=0; j<phi[i].size(); ++j)
	{
	  std::cout << " phi[" << i << "][" << j << "]=" << phi[i][j] << std::endl;
	}
    }
};



inline
void FEBase::print_dphi() const
{
  for (unsigned int i=0; i<dphi.size(); ++i)
    {
      for (unsigned int j=0; j<dphi[i].size(); ++j)
	{
	  std::cout << " dphi[" << i << "][" << j << "]=";
	  dphi[i][j].print();
	}
    }
};



inline
void FEBase::print_xyz() const
{
  for (unsigned int i=0; i<xyz.size(); ++i) xyz[i].print();
};



inline
void FEBase::print_info() const
{
  std::cout << "Shape functions at the Gauss pts." << std::endl;
  print_phi();
  std::cout << "Shape function gradients at the Gauss pts." << std::endl;
  print_dphi();
  std::cout << "XYZ locations of the Gauss pts." << std::endl;
  print_xyz();
  std::cout << "Values of JxW at the Gauss pts." << std::endl;
  print_JxW();
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
};



// ------------------------------------------------------------
// FEMonomial class inline members
template <unsigned int Dim>
inline
FEMonomial<Dim>::FEMonomial (const FEType& fet) :
  FE<Dim,MONOMIAL> (fet)
{
};

#endif
