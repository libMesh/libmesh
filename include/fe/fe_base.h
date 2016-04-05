// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_FE_BASE_H
#define LIBMESH_FE_BASE_H

// Local includes
#include "libmesh/libmesh_common.h"
#include "libmesh/auto_ptr.h"
#include "libmesh/compare_types.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/fe_abstract.h"
#include "libmesh/fe_transformation_base.h"
#include "libmesh/point.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/tensor_tools.h"
#include "libmesh/type_n_tensor.h"
#include "libmesh/vector_value.h"

// C++ includes
#include <cstddef>
#include <vector>


namespace libMesh
{


// forward declarations
template <typename T> class DenseMatrix;
template <typename T> class DenseVector;
class BoundaryInfo;
class DofConstraints;
class DofMap;
class Elem;
class MeshBase;
template <typename T> class NumericVector;
class QBase;
template <typename T> class FETransformationBase;
class FEType;

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
class NodeConstraints;
#endif

#ifdef LIBMESH_ENABLE_PERIODIC
class PeriodicBoundaries;
class PointLocatorBase;
#endif

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
class InfFE;
#endif

/**
 * This class forms the foundation from which generic finite
 * elements may be derived.  In the current implementation the
 * templated derived class \p FE offers a wide variety of commonly
 * used finite element concepts.  Check there for details.
 *
 * Use the \p FEGenericBase<OutputType>::build() method to create an
 * object of any of the derived classes which is compatible with
 * OutputType.
 *
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
 * \author Benjamin S. Kirk
 * \date 2002
 */
template <typename OutputType>
class FEGenericBase : public FEAbstract
{
protected:

  /**
   * Constructor.  Optionally initializes required data
   * structures.  Protected so that this base class
   * cannot be explicitly instantiated.
   */
  FEGenericBase (const unsigned int dim,
                 const FEType & fet);

public:

  /**
   * Destructor.
   */
  virtual ~FEGenericBase();

  /**
   * Builds a specific finite element type.  A \p
   * UniquePtr<FEGenericBase> is returned to prevent a memory leak. This
   * way the user need not remember to delete the object.
   *
   * The build call will fail if the OutputType of this class is not
   * compatible with the output required for the requested \p type
   */
  static UniquePtr<FEGenericBase> build (const unsigned int dim,
                                         const FEType & type);

  /**
   * Convenient typedefs for gradients of output, hessians of output,
   * and potentially-complex-valued versions of same.
   */
  typedef OutputType                                                      OutputShape;
  typedef typename TensorTools::IncrementRank<OutputShape>::type          OutputGradient;
  typedef typename TensorTools::IncrementRank<OutputGradient>::type       OutputTensor;
  typedef typename TensorTools::DecrementRank<OutputShape>::type          OutputDivergence;
  typedef typename TensorTools::MakeNumber<OutputShape>::type             OutputNumber;
  typedef typename TensorTools::IncrementRank<OutputNumber>::type         OutputNumberGradient;
  typedef typename TensorTools::IncrementRank<OutputNumberGradient>::type OutputNumberTensor;
  typedef typename TensorTools::DecrementRank<OutputNumber>::type         OutputNumberDivergence;



#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * Builds a specific infinite element type.  A \p
   * UniquePtr<FEGenericBase> is returned to prevent a memory leak. This
   * way the user need not remember to delete the object.
   *
   * The build call will fail if the OutputShape of this class is not
   * compatible with the output required for the requested \p type
   */
  static UniquePtr<FEGenericBase> build_InfFE (const unsigned int dim,
                                               const FEType & type);

#endif

#ifdef LIBMESH_ENABLE_AMR

  /**
   * Computes the constraint matrix contributions (for
   * non-conforming adapted meshes) corresponding to
   * variable number \p var_number, using generic
   * projections.
   */
  static void compute_proj_constraints (DofConstraints & constraints,
                                        DofMap & dof_map,
                                        const unsigned int variable_number,
                                        const Elem * elem);

  /**
   * Creates a local projection on \p coarse_elem, based on the
   * DoF values in \p global_vector for it's children.  Computes a
   * vector of coefficients corresponding to dof_indices for only the
   * single given \p var
   */

  static void coarsened_dof_values(const NumericVector<Number> & global_vector,
                                   const DofMap & dof_map,
                                   const Elem * coarse_elem,
                                   DenseVector<Number> & coarse_dofs,
                                   const unsigned int var,
                                   const bool use_old_dof_indices = false);

  /**
   * Creates a local projection on \p coarse_elem, based on the
   * DoF values in \p global_vector for it's children.  Computes a
   * vector of coefficients corresponding to all dof_indices.
   */

  static void coarsened_dof_values(const NumericVector<Number> & global_vector,
                                   const DofMap & dof_map,
                                   const Elem * coarse_elem,
                                   DenseVector<Number> & coarse_dofs,
                                   const bool use_old_dof_indices = false);

#endif // #ifdef LIBMESH_ENABLE_AMR

#ifdef LIBMESH_ENABLE_PERIODIC

  /**
   * Computes the constraint matrix contributions (for
   * meshes with periodic boundary conditions) corresponding to
   * variable number \p var_number, using generic projections.
   */
  static void compute_periodic_constraints (DofConstraints & constraints,
                                            DofMap & dof_map,
                                            const PeriodicBoundaries & boundaries,
                                            const MeshBase & mesh,
                                            const PointLocatorBase * point_locator,
                                            const unsigned int variable_number,
                                            const Elem * elem);

#endif // LIBMESH_ENABLE_PERIODIC

  /**
   * @returns the shape function values at the quadrature points
   * on the element.
   */
  const std::vector<std::vector<OutputShape> > & get_phi() const
  { libmesh_assert(!calculations_started || calculate_phi);
    calculate_phi = true; return phi; }

  /**
   * @returns the shape function derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<OutputGradient> > & get_dphi() const
  { libmesh_assert(!calculations_started || calculate_dphi);
    calculate_dphi = calculate_dphiref = true; return dphi; }

  /**
   * @returns the curl of the shape function at the quadrature
   * points.
   */
  const std::vector<std::vector<OutputShape> > & get_curl_phi() const
  { libmesh_assert(!calculations_started || calculate_curl_phi);
    calculate_curl_phi = calculate_dphiref = true; return curl_phi; }

  /**
   * @returns the divergence of the shape function at the quadrature
   * points.
   */
  const std::vector<std::vector<OutputDivergence> > & get_div_phi() const
  { libmesh_assert(!calculations_started || calculate_div_phi);
    calculate_div_phi = calculate_dphiref = true; return div_phi; }

  /**
   * @returns the shape function x-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<OutputShape> > & get_dphidx() const
  { libmesh_assert(!calculations_started || calculate_dphi);
    calculate_dphi = calculate_dphiref = true; return dphidx; }

  /**
   * @returns the shape function y-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<OutputShape> > & get_dphidy() const
  { libmesh_assert(!calculations_started || calculate_dphi);
    calculate_dphi = calculate_dphiref = true; return dphidy; }

  /**
   * @returns the shape function z-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<OutputShape> > & get_dphidz() const
  { libmesh_assert(!calculations_started || calculate_dphi);
    calculate_dphi = calculate_dphiref = true; return dphidz; }

  /**
   * @returns the shape function xi-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<OutputShape> > & get_dphidxi() const
  { libmesh_assert(!calculations_started || calculate_dphiref);
    calculate_dphiref = true; return dphidxi; }

  /**
   * @returns the shape function eta-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<OutputShape> > & get_dphideta() const
  { libmesh_assert(!calculations_started || calculate_dphiref);
    calculate_dphiref = true; return dphideta; }

  /**
   * @returns the shape function zeta-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<OutputShape> > & get_dphidzeta() const
  { libmesh_assert(!calculations_started || calculate_dphiref);
    calculate_dphiref = true; return dphidzeta; }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * @returns the shape function second derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<OutputTensor> > & get_d2phi() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = calculate_dphiref = true; return d2phi; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<OutputShape> > & get_d2phidx2() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = calculate_dphiref = true; return d2phidx2; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<OutputShape> > & get_d2phidxdy() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = calculate_dphiref = true; return d2phidxdy; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<OutputShape> > & get_d2phidxdz() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = calculate_dphiref = true; return d2phidxdz; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<OutputShape> > & get_d2phidy2() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi =  calculate_dphiref = true; return d2phidy2; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<OutputShape> > & get_d2phidydz() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = calculate_dphiref = true; return d2phidydz; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<OutputShape> > & get_d2phidz2() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = calculate_dphiref = true; return d2phidz2; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points, in reference coordinates
   */
  const std::vector<std::vector<OutputShape> > & get_d2phidxi2() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = calculate_dphiref = true; return d2phidxi2; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points, in reference coordinates
   */
  const std::vector<std::vector<OutputShape> > & get_d2phidxideta() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = calculate_dphiref = true; return d2phidxideta; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points, in reference coordinates
   */
  const std::vector<std::vector<OutputShape> > & get_d2phidxidzeta() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = calculate_dphiref = true; return d2phidxidzeta; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points, in reference coordinates
   */
  const std::vector<std::vector<OutputShape> > & get_d2phideta2() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = calculate_dphiref = true; return d2phideta2; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points, in reference coordinates
   */
  const std::vector<std::vector<OutputShape> > & get_d2phidetadzeta() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = calculate_dphiref = true; return d2phidetadzeta; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points, in reference coordinates
   */
  const std::vector<std::vector<OutputShape> > & get_d2phidzeta2() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = calculate_dphiref = true; return d2phidzeta2; }

#endif //LIBMESH_ENABLE_SECOND_DERIVATIVES

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

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
  const std::vector<OutputGradient> & get_dphase() const
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
  const std::vector<Real> & get_Sobolev_weight() const
  { return weight; }

  /**
   * @returns the first global derivative of the multiplicative
   * weight at each quadrature point. See \p get_Sobolev_weight()
   * for details.  In case of \p FE initialized to all zero.
   */
  const std::vector<RealGradient> & get_Sobolev_dweight() const
  { return dweight; }

#endif


  /**
   * Prints the value of each shape function at each quadrature point.
   */
  void print_phi(std::ostream & os) const;

  /**
   * Prints the value of each shape function's derivative
   * at each quadrature point.
   */
  void print_dphi(std::ostream & os) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Prints the value of each shape function's second derivatives
   * at each quadrature point.
   */
  void print_d2phi(std::ostream & os) const;

#endif


protected:



#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * Initialize the data fields for the base of an
   * an infinite element.  Implement this in the derived
   * class \p FE<Dim,T>.
   */
  virtual void init_base_shape_functions(const std::vector<Point> & qp,
                                         const Elem * e) = 0;

#endif

  /**
   * Determine which values are to be calculated, for both the FE
   * itself and for the FEMap.
   */
  void determine_calculations();

  /**
   * After having updated the jacobian and the transformation
   * from local to global coordinates in \p FEAbstract::compute_map(),
   * the first derivatives of the shape functions are
   * transformed to global coordinates, giving \p dphi,
   * \p dphidx, \p dphidy, and \p dphidz. This method
   * should rarely be re-defined in derived classes, but
   * still should be usable for children. Therefore, keep
   * it protected.
   */
  virtual void compute_shape_functions(const Elem * elem, const std::vector<Point> & qp);

  /**
   * Object that handles computing shape function values, gradients, etc
   * in the physical domain.
   */
  UniquePtr<FETransformationBase<OutputType> > _fe_trans;

  /**
   * Shape function values.
   */
  std::vector<std::vector<OutputShape> >   phi;

  /**
   * Shape function derivative values.
   */
  std::vector<std::vector<OutputGradient> >  dphi;

  /**
   * Shape function curl values. Only defined for vector types.
   */
  std::vector<std::vector<OutputShape> > curl_phi;

  /**
   * Shape function divergence values. Only defined for vector types.
   */
  std::vector<std::vector<OutputDivergence> > div_phi;

  /**
   * Shape function derivatives in the xi direction.
   */
  std::vector<std::vector<OutputShape> >   dphidxi;

  /**
   * Shape function derivatives in the eta direction.
   */
  std::vector<std::vector<OutputShape> >   dphideta;

  /**
   * Shape function derivatives in the zeta direction.
   */
  std::vector<std::vector<OutputShape> >   dphidzeta;

  /**
   * Shape function derivatives in the x direction.
   */
  std::vector<std::vector<OutputShape> >   dphidx;

  /**
   * Shape function derivatives in the y direction.
   */
  std::vector<std::vector<OutputShape> >   dphidy;

  /**
   * Shape function derivatives in the z direction.
   */
  std::vector<std::vector<OutputShape> >   dphidz;


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Shape function second derivative values.
   */
  std::vector<std::vector<OutputTensor> >  d2phi;

  /**
   * Shape function second derivatives in the xi direction.
   */
  std::vector<std::vector<OutputShape> >   d2phidxi2;

  /**
   * Shape function second derivatives in the xi-eta direction.
   */
  std::vector<std::vector<OutputShape> >   d2phidxideta;

  /**
   * Shape function second derivatives in the xi-zeta direction.
   */
  std::vector<std::vector<OutputShape> >   d2phidxidzeta;

  /**
   * Shape function second derivatives in the eta direction.
   */
  std::vector<std::vector<OutputShape> >   d2phideta2;

  /**
   * Shape function second derivatives in the eta-zeta direction.
   */
  std::vector<std::vector<OutputShape> >   d2phidetadzeta;

  /**
   * Shape function second derivatives in the zeta direction.
   */
  std::vector<std::vector<OutputShape> >   d2phidzeta2;

  /**
   * Shape function second derivatives in the x direction.
   */
  std::vector<std::vector<OutputShape> >   d2phidx2;

  /**
   * Shape function second derivatives in the x-y direction.
   */
  std::vector<std::vector<OutputShape> >   d2phidxdy;

  /**
   * Shape function second derivatives in the x-z direction.
   */
  std::vector<std::vector<OutputShape> >   d2phidxdz;

  /**
   * Shape function second derivatives in the y direction.
   */
  std::vector<std::vector<OutputShape> >   d2phidy2;

  /**
   * Shape function second derivatives in the y-z direction.
   */
  std::vector<std::vector<OutputShape> >   d2phidydz;

  /**
   * Shape function second derivatives in the z direction.
   */
  std::vector<std::vector<OutputShape> >   d2phidz2;

#endif


#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  //--------------------------------------------------------------
  /* protected members for infinite elements, which are accessed
   * from the outside through some inline functions
   */


  /**
   * Used for certain @e infinite element families:
   * the first derivatives of the phase term in global coordinates,
   * over @e all quadrature points.
   */
  std::vector<OutputGradient> dphase;

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

private:

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * Make all \p InfFE<Dim,T_radial,T_map> classes friends
   * so that they can safely used \p FE<Dim-1,T_base> through
   * a \p FEGenericBase * as base approximation.
   */
  template <unsigned int friend_Dim, FEFamily friend_T_radial, InfMapType friend_T_map>
  friend class InfFE;

#endif


};


// Typedefs for convenience and backwards compatibility
typedef FEGenericBase<Real> FEBase;
typedef FEGenericBase<RealGradient> FEVectorBase;




// ------------------------------------------------------------
// FEGenericBase class inline members
template <typename OutputType>
inline
FEGenericBase<OutputType>::FEGenericBase(const unsigned int d,
                                         const FEType & fet) :
  FEAbstract(d,fet),
  _fe_trans( FETransformationBase<OutputType>::build(fet) ),
  phi(),
  dphi(),
  curl_phi(),
  div_phi(),
  dphidxi(),
  dphideta(),
  dphidzeta(),
  dphidx(),
  dphidy(),
  dphidz()
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  ,d2phi(),
  d2phidxi2(),
  d2phidxideta(),
  d2phidxidzeta(),
  d2phideta2(),
  d2phidetadzeta(),
  d2phidzeta2(),
  d2phidx2(),
  d2phidxdy(),
  d2phidxdz(),
  d2phidy2(),
  d2phidydz(),
  d2phidz2()
#endif
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  ,dphase(),
  dweight(),
  weight()
#endif
{
}



template <typename OutputType>
inline
FEGenericBase<OutputType>::~FEGenericBase()
{
}

} // namespace libMesh

#endif // LIBMESH_FE_BASE_H
