// $Id: fe_base.h,v 1.4 2004-02-18 23:04:08 benkirk Exp $

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



#ifndef __fe_base_h__
#define __fe_base_h__

// C++ includes
#include <vector>

// Local includes
#include "reference_counted_object.h"
#include "point.h"
#include "vector_value.h"
#include "enum_elem_type.h"
#include "fe_type.h"
#include "auto_ptr.h"


// forward declarations
class QBase;
class MeshBase;
class Elem;

#ifdef ENABLE_INFINITE_ELEMENTS

template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
class InfFE;

#endif



/**
 * This class forms the foundation from which generic finite
 * elements may be derived.  In the current implementation the
 * templated derived class \p FE offers a wide variety of commonly 
 * used finite element concepts.  Check there for details.
 * Use the \p FEBase::build() method to create an object of any of 
 * the derived classes.
 * Note that the amount of virtual members is kept to a minimum,
 * and the sophisticated template scheme of \p FE is quite
 * likely to offer acceptably fast code.
 *
 * All calls to static members of the \p FE classes should be
 * requested through the \p FEInterface.  This interface class
 * offers sort-of runtime polymorphism for the templated finite
 * element classes.  Even internal library classes, like \p DofMap,
 * request the number of dof's through this interface class.
 * Note that this also enables the co-existence of various
 * element-based schemes.  
 * This class is well 'at the heart' of the library, so 
 * things in here should better remain unchanged. 
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
   * Builds a specific finite element type.  A \p AutoPtr<FEBase> is
   * returned to prevent a memory leak. This way the user need not
   * remember to delete the object.
   */
  static AutoPtr<FEBase> build (const unsigned int dim,
				const FEType& type); 
  
#ifdef ENABLE_INFINITE_ELEMENTS

  /**
   * Builds a specific infinite element type.  A \p AutoPtr<FEBase> is
   * returned to prevent a memory leak. This way the user need not
   * remember to delete the object.
   */
  static AutoPtr<FEBase> build_InfFE (const unsigned int dim,
				      const FEType& type); 

#endif

  /**
   * This is at the core of this class. Use this for each
   * new element in the mesh.  Reinitializes all the physical 
   * element-dependent data based on the current element 
   * \p elem. By default the shape functions and associated
   * data are computed at the quadrature points specified
   * by the quadrature rule \p qrule, but may be any points
   * specified on the reference element specified in the optional
   * argument \p pts.
   */
  virtual void reinit (const Elem* elem,
		       const std::vector<Point>* const pts = NULL) = 0;
    
  /**
   * Reinitializes all the physical element-dependent data based on
   * the \p side of the element \p elem.
   */
  virtual void reinit (const Elem* elem,
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
				   const Real eps = TOLERANCE);
  
  /**
   * @returns the \p xyz spatial locations of the quadrature
   * points on the element.
   */    
  const std::vector<Point>& get_xyz() const
  { return xyz; }
  
  /**
   * @returns the shape function values at the quadrature points
   * on the element.
   */    
  const std::vector<std::vector<Real> >& get_phi() const
  { return phi; }
  
  /**
   * @returns the element Jacobian times the quadrature weight for
   * each quadrature point.
   */    
  const std::vector<Real>& get_JxW() const
  { return JxW; }

  /**
   * @returns the shape function derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<RealGradient> >& get_dphi() const
  { return dphi; }
  
  /**
   * @returns the shape function x-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<Real> >& get_dphidx() const
  { return dphidx; }
  
  /**
   * @returns the shape function y-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<Real> >& get_dphidy() const
  { return dphidy; }
  
  /**
   * @returns the shape function z-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<Real> >& get_dphidz() const
  { return dphidz; }

  
#ifdef ENABLE_INFINITE_ELEMENTS

  /**
   * @returns the global first derivative of the phase term 
   * which is used in infinite elements, evaluated at the 
   * quadrature points.  
   *
   * In case of the general finite element class \p FE this 
   * field is initialized to all zero, so that the variational 
   * formulation for an @e infinite element returns correct element 
   * matrices for a mesh using both finite and infinite elements.
   */
  const std::vector<RealGradient>& get_dphase() const
      { return dphase; }


  /**
   * @returns the multiplicative weight at each quadrature point.
   * This weight is used for certain infinite element weak 
   * formulations, so that @e weighted Sobolev spaces are
   * used for the trial function space.  This renders the
   * variational form easily computable. 
   *
   * In case of the general finite element class \p FE this 
   * field is initialized to all ones, so that the variational 
   * formulation for an @e infinite element returns correct element 
   * matrices for a mesh using both finite and infinite elements.
   */
  const std::vector<Real>& get_Sobolev_weight() const
      { return weight; }

  /**
   * @returns the first global derivative of the multiplicative 
   * weight at each quadrature point. See \p get_Sobolev_weight()
   * for details.  In case of \p FE initialized to all zero.
   */
  const std::vector<RealGradient>& get_Sobolev_dweight() const
      { return dweight; }

#endif

  
  /**
   * @returns the tangent vectors for face integration.
   */
  const std::vector<std::vector<Point> >& get_tangents() const
  { return tangents; }
  
  /**
   * @returns the normal vectors for face integration.
   */
  const std::vector<Point>& get_normals() const
  { return normals; }

  /**
   * @returns the curvatures for use in face integration.
   */
  const std::vector<Real>& get_curvatures() const
  { return curvatures;}
  
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
  ElemType get_type()  const { return elem_type; }

  /**
   * @returns the approximation order of the finite element.
   */
  Order get_order()  const { return fe_type.order; }

  /**
   * @returns the finite element family of this element.
   */
  FEFamily get_family()  const { return fe_type.family; }

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
  virtual void init_base_shape_functions(const std::vector<Point>& qp,
					 const Elem* e) = 0;

#endif

  /**
   * Compute the jacobian and some other additional
   * data fields. Takes the integration weights
   * as input, along with a pointer to the element.
   */
  void compute_map(const std::vector<Real>& qw,
		   const Elem* e);
  
  /** 
   * Same as before, but for a side.  Useful for boundary integration.
   */  
  void compute_face_map(const std::vector<Real>& qw,
			const Elem* side);

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
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the x value of the pth entry of the dxzydxi_map.
   */
  Real dxdxi_map(const unsigned int p) const   { return dxyzdxi_map[p](0); }

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the y value of the pth entry of the dxzydxi_map.
   */
  Real dydxi_map(const unsigned int p) const   { return dxyzdxi_map[p](1); }

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the z value of the pth entry of the dxzydxi_map.
   */
  Real dzdxi_map(const unsigned int p) const   { return dxyzdxi_map[p](2); }

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the x value of the pth entry of the dxzydeta_map.
   */
  Real dxdeta_map(const unsigned int p) const  { return dxyzdeta_map[p](0); }

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the y value of the pth entry of the dxzydeta_map.
   */
  Real dydeta_map(const unsigned int p) const  { return dxyzdeta_map[p](1); } 

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the z value of the pth entry of the dxzydeta_map.
   */
  Real dzdeta_map(const unsigned int p) const  { return dxyzdeta_map[p](2); }

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the x value of the pth entry of the dxzydzeta_map.
   */
  Real dxdzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](0); }

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the y value of the pth entry of the dxzydzeta_map.
   */
  Real dydzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](1); }

  /**
   * Used in \p FEBase::compute_map(), which should be
   * be usable in derived classes, and therefore protected.
   * Returns the z value of the pth entry of the dxzydzeta_map.
   */
  Real dzdzeta_map(const unsigned int p) const { return dxyzdzeta_map[p](2); }




  
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
  std::vector<RealGradient> dxyzdxi_map;

  /**
   * Vector of parital derivatives: 
   * d(x)/d(eta), d(y)/d(eta), d(z)/d(eta)
   */
  std::vector<RealGradient> dxyzdeta_map;

  /**
   * Vector of parital derivatives: 
   * d(x)/d(zeta), d(y)/d(zeta), d(z)/d(zeta)
   */
  std::vector<RealGradient> dxyzdzeta_map;
  
  /**
   * Vector of second partial derivatives in xi:
   * d^2(x)/d(xi)^2, d^2(y)/d(xi)^2, d^2(z)/d(xi)^2
   */
  std::vector<RealGradient> d2xyzdxi2_map;

  /**
   * Vector of mixed second partial derivatives in xi-eta:
   * d^2(x)/d(xi)d(eta) d^2(y)/d(xi)d(eta) d^2(z)/d(xi)d(eta)
   */
  std::vector<RealGradient> d2xyzdxideta_map;

  /**
   * Vector of second partial derivatives in eta:
   * d^2(x)/d(eta)^2 
   */
  std::vector<RealGradient> d2xyzdeta2_map;
  
  /**
   * Map for partial derivatives:
   * d(xi)/d(x). Needed for the Jacobian.
   */
  std::vector<Real>  dxidx_map;

  /**
   * Map for partial derivatives:
   * d(xi)/d(y). Needed for the Jacobian.
   */
  std::vector<Real>  dxidy_map;

  /**
   * Map for partial derivatives:
   * d(xi)/d(z). Needed for the Jacobian.
   */
  std::vector<Real>  dxidz_map;



  
  /**
   * Map for partial derivatives:
   * d(eta)/d(x). Needed for the Jacobian.
   */
  std::vector<Real>  detadx_map;

  /**
   * Map for partial derivatives:
   * d(eta)/d(y). Needed for the Jacobian.
   */
  std::vector<Real>  detady_map;

  /**
   * Map for partial derivatives:
   * d(eta)/d(z). Needed for the Jacobian.
   */
  std::vector<Real>  detadz_map;




  
  /**
   * Map for partial derivatives:
   * d(zeta)/d(x). Needed for the Jacobian.
   */
  std::vector<Real>  dzetadx_map;

  /**
   * Map for partial derivatives:
   * d(zeta)/d(y). Needed for the Jacobian.
   */
  std::vector<Real>  dzetady_map;

  /**
   * Map for partial derivatives:
   * d(zeta)/d(z). Needed for the Jacobian.
   */
  std::vector<Real>  dzetadz_map;


  
  /**
   * Shape function values.
   */
  std::vector<std::vector<Real> >   phi;

  /**
   * Shape function derivative values.
   */
  std::vector<std::vector<RealGradient> >  dphi;

  /**
   * Shape function derivatives in the xi direction.
   */
  std::vector<std::vector<Real> >   dphidxi;

  /**
   * Shape function derivatives in the eta direction.
   */
  std::vector<std::vector<Real> >   dphideta;
  
  /**
   * Shape function derivatives in the zeta direction.
   */
  std::vector<std::vector<Real> >   dphidzeta;

  /**
   * Shape function derivatives in the x direction.
   */
  std::vector<std::vector<Real> >   dphidx;

  /**
   * Shape function derivatives in the y direction.
   */
  std::vector<std::vector<Real> >   dphidy;

  /**
   * Shape function derivatives in the z direction.
   */
  std::vector<std::vector<Real> >   dphidz;




  
  /**
   * Map for the shape function phi.
   */
  std::vector<std::vector<Real> >   phi_map;

  /**
   * Map for the derivative, d(phi)/d(xi).
   */
  std::vector<std::vector<Real> >   dphidxi_map;

  /**
   * Map for the derivative, d(phi)/d(eta).
   */
  std::vector<std::vector<Real> >   dphideta_map;

  /**
   * Map for the derivative, d(phi)/d(zeta).
   */
  std::vector<std::vector<Real> >   dphidzeta_map;




  
  /**
   * Map for the side shape functions, psi. 
   */
  std::vector<std::vector<Real> >   psi_map;

  /**
   * Map for the derivative of the side functions,
   * d(psi)/d(xi).
   */
  std::vector<std::vector<Real> >   dpsidxi_map;

  /**
   * Map for the derivative of the side function,
   * d(psi)/d(eta).
   */
  std::vector<std::vector<Real> >   dpsideta_map;

  /**
   * Map for the second derivatives (in xi) of the
   * side shape functions.  Useful for computing
   * the curvature at the quadrature points.
   */
  std::vector<std::vector<Real> > d2psidxi2_map;

  /**
   * Map for the second (cross) derivatives in xi, eta
   * of the side shape functions.  Useful for
   * computing the curvature at the quadrature points.
   */
  std::vector<std::vector<Real> > d2psidxideta_map;

  /**
   * Map for the second derivatives (in eta) of the
   * side shape functions.  Useful for computing the
   * curvature at the quadrature points.
   */
  std::vector<std::vector<Real> > d2psideta2_map;
  
#ifdef ENABLE_INFINITE_ELEMENTS

  //--------------------------------------------------------------
  /* protected members for infinite elements, which are accessed 
   * from the outside through some inline functions
   */


  /**
   * Used for certain @e infinite element families:
   * the first derivatives of the phase term in global coordinates,
   * over @e all quadrature points.
   */
  std::vector<RealGradient> dphase;

  /**
   * Used for certain @e infinite element families:
   * the global derivative of the additional radial weight \f$ 1/{r^2} \f$,
   * over @e all quadrature points.
   */
  std::vector<RealGradient> dweight;

  /**
   * Used for certain @e infinite element families:
   * the additional radial weight \f$ 1/{r^2} \f$ in local coordinates,
   * over @e all quadrature points.
   */
  std::vector<Real>  weight;

#endif



  
  /**
   * Tangent vectors on boundary at quadrature points.
   */
  std::vector<std::vector<Point> >  tangents;

  /**
   * Normal vectors on boundary at quadrature points
   */
  std::vector<Point>                normals;

  /**
   * The mean curvature (= one half the sum of the principal
   * curvatures) on the boundary at the quadrature points.
   * The mean curvature is a scalar value.
   */
  std::vector<Real>                 curvatures;
  
  /**
   * Jacobian*Weight values at quadrature points
   */
  std::vector<Real>                 JxW;

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
   * @returns \p true when the shape functions (for
   * this \p FEFamily) depend on the particular
   * element, and therefore needs to be re-initialized
   * for each new element.  \p false otherwise.
   */
  virtual bool shapes_need_reinit() const = 0;


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
}



inline
FEBase::~FEBase()
{
}



inline
void FEBase::print_JxW() const
{
  for (unsigned int i=0; i<JxW.size(); ++i) std::cout << JxW[i] << std::endl;
}



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
}



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
}



inline
void FEBase::print_xyz() const
{
  for (unsigned int i=0; i<xyz.size(); ++i) xyz[i].print();
}



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
}



#endif
