// $Id: fe_base.h,v 1.1 2003-01-24 17:24:38 jwpeterson Exp $

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



#ifndef __fe_base_h__
#define __fe_base_h__

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

#ifdef ENABLE_INFINITE_ELEMENTS

template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
class InfFE;

#endif



/**
 * This class forms the foundation from which generic finite
 * elements may be derived. The current implementation offers
 * a wide variety of commonly used finite element concepts.
 * For actual use in simulations, take the templated FE
 * class.  Through static member functions this templated class implements 
 * Lagrange shape functions, monomials, and hierarchic shape functions.
 * These may be used, through the public methods FEBase::shape() and
 * FEBase::shape_deriv(), as building blocks for shape functions in
 * derived classes.
 *
 * All interaction of this and derived classes with other classes, 
 * like DofMap, are handled through the interface class FEInterface. 
 * When the static member functions like FEBase::n_dofs() etc are 
 * not sufficient to represent your derived class, add calls to your 
 * version in FEInterface, better do not modify the methods from FEBase.
 * Within the FEBase class, things should remain unchanged. 
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

  
#ifdef ENABLE_INFINITE_ELEMENTS

  /**
   * @returns the global first derivative of the phase term in
   * infinite elements, evaluated at the quadrature points.
   * To be implemented in derived classes.
   */
  virtual const std::vector<Point>& get_dphase() const = 0;

  /**
   * @returns the multiplicative weight at each quadrature point.
   * To be implemented in derived classes.
   */
  virtual const std::vector<real>& get_Sobolev_weight() const = 0;

  /**
   * @returns the first global derivative of the multiplicative 
   * weight at each quadrature point.  To be implemented in 
   * derived classes.
   */
  virtual const std::vector<Point>& get_Sobolev_dweight() const = 0;

#endif

  
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
   * Provides the class with the quadrature rule.  Implement
   * this in derived classes.
   */
  virtual void attach_quadrature_rule (QBase* q) = 0;

  /**
   * @returns the total number of approximation shape functions
   * for the current element.  Useful during matrix assembly.  
   * Implement this in derived classes.
   */
  virtual unsigned int n_shape_functions () const = 0;

  /**
   * @returns the total number of quadrature points.  Useful
   * during matrix assembly.  Implement this in derived classes.
   */
  virtual unsigned int n_quadrature_points () const = 0;
  
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


    
#ifdef ENABLE_INFINITE_ELEMENTS

  /**
   * Initialize the data fields for the base of an
   * an infinite element.  Implement this in the derived 
   * class \p FE<Dim,T>.
   */
  virtual void init_base_shape_functions(const QBase* q,
					 const Elem* e) = 0;

#endif

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
   * from local to global coordinates in \p FEBase::compute_map(),
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


#ifdef ENABLE_INFINITE_ELEMENTS

  /**
   * Make all \p InfFE<Dim,T_radial,T_map> classes friends
   * so that they can safely used \p FE<Dim-1,T_base> through
   * a \p FEBase* as base approximation.
   */
  template <unsigned int friend_Dim, FEFamily friend_T_radial, InfMapType friend_T_map>
  friend class InfFE;

#endif


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



#endif
