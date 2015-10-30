// $Id: fe.h,v 1.1.1.1 2003-01-10 16:17:48 libmesh Exp $

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
#include <memory>

// Local includes
#include "mesh_config.h"
#include "reference_counter.h"
#include "mesh_base.h"
#include "point.h"
#include "order.h"
#include "fe_type.h"


// forward declarations
class QBase;



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

class FEBase : public ReferenceCounter
{
protected:

  /**
   * Constructor.  Optionally initializes required data
   * structures.  Protected so that this base class
   * cannot be explicitly instantiated.
   */
  FEBase (const MeshBase& m,
	  const unsigned int dim,
	  const FEType& fet);
  
public:
  
  /**
   * Destructor.
   */
  virtual ~FEBase();

  /**
   * @returns "FEBase"
   */
  std::string class_name () const { return "FEBase"; };
  
  /**
   * Builds a specific finite element type.  A \p std::auto_ptr<Elem> is
   * returned to prevent a memory leak. This way the user need not
   * remember to delete the object.
   */
  static std::auto_ptr<FEBase> build (const MeshBase& m,
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
  

  // All these nice short-hand inline functions
  // are used in \p FEBase::compute_map(), which should be
  // be usable in derived classes, and therefore protected.
  real dxdxi_map(const unsigned int p) const   { return dxyzdxi_map[p](0); };

  real dydxi_map(const unsigned int p) const   { return dxyzdxi_map[p](1); };

  real dzdxi_map(const unsigned int p) const   { return dxyzdxi_map[p](2); };

  real dxdeta_map(const unsigned int p) const  { return dxyzdeta_map[p](0); };

  real dydeta_map(const unsigned int p) const  { return dxyzdeta_map[p](1); }; 

  real dzdeta_map(const unsigned int p) const  { return dxyzdeta_map[p](2); };

  real dxdzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](0); };

  real dydzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](1); };

  real dzdzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](2); };


  /**
   * A reference to the mesh object
   */
  const MeshBase& mesh;

  /**
   * The dimensionality of the object
   */
  const unsigned int dim;

  /**
   * The spatial locations of the quadrature points
   */
  std::vector<Point> xyz;

  // Vectors used to build the map from/to
  // the reference element
  std::vector<Point> dxyzdxi_map;
  std::vector<Point> dxyzdeta_map;
  std::vector<Point> dxyzdzeta_map;
  		     
  std::vector<real>  dxidx_map;
  std::vector<real>  dxidy_map;
  std::vector<real>  dxidz_map;
		     
  std::vector<real>  detadx_map;
  std::vector<real>  detady_map;
  std::vector<real>  detadz_map;
		     
  std::vector<real>  dzetadx_map;
  std::vector<real>  dzetady_map;
  std::vector<real>  dzetadz_map;
				    
  // FE Approximation shape functions    
  std::vector<std::vector<real> >   phi;
  std::vector<std::vector<Point> >  dphi;
  std::vector<std::vector<real> >   dphidxi;
  std::vector<std::vector<real> >   dphideta;
  std::vector<std::vector<real> >   dphidzeta;
  std::vector<std::vector<real> >   dphidx;
  std::vector<std::vector<real> >   dphidy;
  std::vector<std::vector<real> >   dphidz;

  // Mapping shape functions
  std::vector<std::vector<real> >   phi_map;
  std::vector<std::vector<real> >   dphidxi_map;
  std::vector<std::vector<real> >   dphideta_map;
  std::vector<std::vector<real> >   dphidzeta_map;

  // Shape functions for mapping side values
  std::vector<std::vector<real> >   psi_map;
  std::vector<std::vector<real> >   dpsidxi_map;
  std::vector<std::vector<real> >   dpsideta_map;

  // Normal vectors on boundary at quadrature points
  std::vector<std::vector<Point> >  tangents;
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
 * \version $Revision: 1.1.1.1 $
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
  FE(const MeshBase& m,
     const FEType& fet);
  
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
  static void nodal_soln(const MeshBase& mesh,
			 const Elem* elem, const Order o,
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
  static Point inverse_map (const MeshBase& mesh,
			    const Elem* elem,
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
  static Point map (const MeshBase& mesh,
		    const Elem* elem,
		    const Point& reference_point);
  
  /**
   * @returns d(xyz)/dxi (in physical space) of the point
   * \p p located on the reference element.
   */
  static Point map_xi (const MeshBase& mesh,
		       const Elem* elem,
		       const Point& reference_point);
  
  /**
   * @returns d(xyz)/deta (in physical space) of the point
   * \p p located on the reference element.
   */
  static Point map_eta (const MeshBase& mesh,
			const Elem* elem,
			const Point& reference_point);

  /**
   * @returns d(xyz)/dzeta (in physical space) of the point
   * \p p located on the reference element.
   */
  static Point map_zeta (const MeshBase& mesh,
			 const Elem* elem,
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
 * \version $Revision: 1.1.1.1 $
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
  FEHierarchic(const MeshBase& m,
	       const FEType&   fet);
};



/**
 * Lagrange finite elements.  Still templated on the dimension,
 * \p Dim.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2003
 * \version $Revision: 1.1.1.1 $
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
  FELagrange(const MeshBase& m,
	     const FEType&   fet);
};



/**
 * Monomial finite elements.  Still templated on the dimension,
 * \p Dim.
 *
 * \author Benjamin S. Kirk
 * \date 2002-2003
 * \version $Revision: 1.1.1.1 $
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
  FEMonomial(const MeshBase& m,
	     const FEType&   fet);
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
FEBase::FEBase(const MeshBase& m,
	       const unsigned int d,
	       const FEType& fet) :
  mesh(m),
  dim(d),
  fe_type(fet),
  elem_type(INVALID_ELEM),
  qrule(NULL)
{
  increment_constructor_count ();
};



inline
FEBase::~FEBase()
{
  increment_destructor_count ();
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
FE<Dim,T>::FE (const MeshBase& m,
	       const FEType& fet) :
  FEBase (m,Dim,fet)
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
FEHierarchic<Dim>::FEHierarchic (const MeshBase& m,
				 const FEType&   fet) :
  FE<Dim,HIERARCHIC> (m,fet)
{
};



// ------------------------------------------------------------
// FELagrange class inline members
template <unsigned int Dim>
inline
FELagrange<Dim>::FELagrange (const MeshBase& m,
			     const FEType&   fet) :
  FE<Dim,LAGRANGE> (m,fet)
{
};



// ------------------------------------------------------------
// FEMonomial class inline members
template <unsigned int Dim>
inline
FEMonomial<Dim>::FEMonomial (const MeshBase& m,
			     const FEType&   fet) :
  FE<Dim,MONOMIAL> (m,fet)
{
};

#endif
