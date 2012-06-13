// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef __fe_vector_base_h__
#define __fe_vector_base_h__

// Local includes
#include "reference_counted_object.h"
#include "fe_abstract.h"
#include "point.h"
#include "vector_value.h"
#include "enum_elem_type.h"
#include "fe_type.h"
#include "auto_ptr.h"

// C++ includes
#include <cstddef>
#include <vector>

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
#include "tensor_value.h"
#endif

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

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
class NodeConstraints;
#endif

#ifdef LIBMESH_ENABLE_PERIODIC
class PeriodicBoundaries;
class PointLocatorBase;
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
class FEVectorBase : public FEAbstract
{
protected:

  /**
   * Constructor.  Optionally initializes required data
   * structures.  Protected so that this base class
   * cannot be explicitly instantiated.
   */
  FEVectorBase (const unsigned int dim,
		const FEType& fet);

public:

  /**
   * Destructor.
   */
  virtual ~FEVectorBase();

  /**
   * Builds a specific finite element type.  A \p AutoPtr<FEVectorBase> is
   * returned to prevent a memory leak. This way the user need not
   * remember to delete the object.
   */
  static AutoPtr<FEVectorBase> build (const unsigned int dim,
				      const FEType& type);

#ifdef LIBMESH_ENABLE_AMR

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

  /**
   * Creates a local projection on \p coarse_elem, based on the
   * DoF values in \p global_vector for it's children.
   */

  static void coarsened_dof_values(const NumericVector<Number> &global_vector,
			           const DofMap &dof_map,
                                   const Elem *coarse_elem,
			           DenseVector<Number> &coarse_dofs,
			           const unsigned int var,
			           const bool use_old_dof_indices = false);

#endif // #ifdef LIBMESH_ENABLE_AMR

#ifdef LIBMESH_ENABLE_PERIODIC

  /**
   * Computes the constraint matrix contributions (for
   * meshes with periodic boundary conditions) corresponding to
   * variable number \p var_number, using generic projections.
   */
  static void compute_periodic_constraints (DofConstraints &constraints,
                                            DofMap &dof_map,
                                            const PeriodicBoundaries &boundaries,
                                            const MeshBase& mesh,
                                            const PointLocatorBase* point_locator,
                                            const unsigned int variable_number,
                                            const Elem* elem);

#endif // LIBMESH_ENABLE_PERIODIC

  /**
   * @returns the shape function values at the quadrature points
   * on the element.
   */
  const std::vector<std::vector<RealGradient> >& get_phi() const
  { libmesh_assert(!calculations_started || calculate_phi);
    calculate_phi = true; return phi; }

  /**
   * @returns the shape function derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<RealTensor> >& get_dphi() const
  { libmesh_assert(!calculations_started || calculate_dphi);
    calculate_dphi = true; return dphi; }

  /**
   * @returns the shape function x-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<RealGradient> >& get_dphidx() const
  { libmesh_assert(!calculations_started || calculate_dphi);
    calculate_dphi = true; return dphidx; }

  /**
   * @returns the shape function y-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<RealGradient> >& get_dphidy() const
  { libmesh_assert(!calculations_started || calculate_dphi);
    calculate_dphi = true; return dphidy; }

  /**
   * @returns the shape function z-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<RealGradient> >& get_dphidz() const
  { libmesh_assert(!calculations_started || calculate_dphi);
    calculate_dphi = true; return dphidz; }

  /**
   * @returns the shape function xi-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<RealGradient> >& get_dphidxi() const
  { libmesh_assert(!calculations_started || calculate_dphi);
    calculate_dphi = true; return dphidxi; }

  /**
   * @returns the shape function eta-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<RealGradient> >& get_dphideta() const
  { libmesh_assert(!calculations_started || calculate_dphi);
    calculate_dphi = true; return dphideta; }

  /**
   * @returns the shape function zeta-derivative at the quadrature
   * points.
   */
  const std::vector<std::vector<RealGradient> >& get_dphidzeta() const
  { libmesh_assert(!calculations_started || calculate_dphi);
    calculate_dphi = true; return dphidzeta; }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * @returns the shape function second derivatives at the quadrature
   * points.
   */
  /*
  const std::vector<std::vector<RealType3Tensor> >& get_d2phi() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = true; return d2phi; }
  */

  /**
   * @returns the shape function second derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<RealTensor> >& get_d2phidx2() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = true; return d2phidx2; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<RealTensor> >& get_d2phidxdy() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = true; return d2phidxdy; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<RealTensor> >& get_d2phidxdz() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = true; return d2phidxdz; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<RealTensor> >& get_d2phidy2() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = true; return d2phidy2; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<RealTensor> >& get_d2phidydz() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = true; return d2phidydz; }

  /**
   * @returns the shape function second derivatives at the quadrature
   * points.
   */
  const std::vector<std::vector<RealTensor> >& get_d2phidz2() const
  { libmesh_assert(!calculations_started || calculate_d2phi);
    calculate_d2phi = true; return d2phidz2; }

#endif

  /**
   * Prints the value of each shape function at each quadrature point.
   */
  void print_phi(std::ostream& os) const;

  /**
   * Prints the value of each shape function's derivative
   * at each quadrature point.
   */
  void print_dphi(std::ostream& os) const;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Prints the value of each shape function's second derivatives
   * at each quadrature point.
   */
  //void print_d2phi(std::ostream& os) const;

#endif


protected:

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
  virtual void compute_shape_functions(const Elem*);


  /**
   * Shape function values.
   */
  std::vector<std::vector<RealGradient> >   phi;

  /**
   * Shape function derivative values.
   */
  std::vector<std::vector<RealTensor> >  dphi;

  /**
   * Shape function derivatives in the xi direction.
   */
  std::vector<std::vector<RealGradient> >   dphidxi;

  /**
   * Shape function derivatives in the eta direction.
   */
  std::vector<std::vector<RealGradient> >   dphideta;

  /**
   * Shape function derivatives in the zeta direction.
   */
  std::vector<std::vector<RealGradient> >   dphidzeta;

  /**
   * Shape function derivatives in the x direction.
   */
  std::vector<std::vector<RealGradient> >   dphidx;

  /**
   * Shape function derivatives in the y direction.
   */
  std::vector<std::vector<RealGradient> >   dphidy;

  /**
   * Shape function derivatives in the z direction.
   */
  std::vector<std::vector<RealGradient> >   dphidz;


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

  /**
   * Shape function second derivative values.
   * TODO: Need to implement 3rd order tensor object to handle
   *       second derivatives of vector-valued elements
   */
  //std::vector<std::vector<RealType3Tensor> >  d2phi;

  /**
   * Shape function second derivatives in the xi direction.
   */
  std::vector<std::vector<RealTensor> >   d2phidxi2;

  /**
   * Shape function second derivatives in the xi-eta direction.
   */
  std::vector<std::vector<RealTensor> >   d2phidxideta;

  /**
   * Shape function second derivatives in the xi-zeta direction.
   */
  std::vector<std::vector<RealTensor> >   d2phidxidzeta;

  /**
   * Shape function second derivatives in the eta direction.
   */
  std::vector<std::vector<RealTensor> >   d2phideta2;

  /**
   * Shape function second derivatives in the eta-zeta direction.
   */
  std::vector<std::vector<RealTensor> >   d2phidetadzeta;

  /**
   * Shape function second derivatives in the zeta direction.
   */
  std::vector<std::vector<RealTensor> >   d2phidzeta2;

  /**
   * Shape function second derivatives in the x direction.
   */
  std::vector<std::vector<RealTensor> >   d2phidx2;

  /**
   * Shape function second derivatives in the x-y direction.
   */
  std::vector<std::vector<RealTensor> >   d2phidxdy;

  /**
   * Shape function second derivatives in the x-z direction.
   */
  std::vector<std::vector<RealTensor> >   d2phidxdz;

  /**
   * Shape function second derivatives in the y direction.
   */
  std::vector<std::vector<RealTensor> >   d2phidy2;

  /**
   * Shape function second derivatives in the y-z direction.
   */
  std::vector<std::vector<RealTensor> >   d2phidydz;

  /**
   * Shape function second derivatives in the z direction.
   */
  std::vector<std::vector<RealTensor> >   d2phidz2;

#endif

};




// ------------------------------------------------------------
// FEVectorBase class inline members
inline
FEVectorBase::FEVectorBase(const unsigned int d,
	       const FEType& fet) :
  FEAbstract(d,fet),
  phi(),
  dphi(),
  dphidxi(),
  dphideta(),
  dphidzeta(),
  dphidx(),
  dphidy(),
  dphidz()
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  //d2phi(),
  , d2phidxi2(),
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
{
}



inline
FEVectorBase::~FEVectorBase()
{
}

} // namespace libMesh

#endif
