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



#ifndef LIBMESH_HP_COARSENTEST_H
#define LIBMESH_HP_COARSENTEST_H

// Local Includes
#include "libmesh/auto_ptr.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/hp_selector.h"
#include "libmesh/libmesh_common.h"

#include "libmesh/fe.h"         // MipsPro requires fe.h and quadrature.h in order to
#include "libmesh/quadrature.h" //  delete AutoPtrs<> upon destruction

// C++ includes
#include <vector>

#ifdef LIBMESH_ENABLE_AMR

namespace libMesh
{

// Forward Declarations
class Elem;
class Point;
class System;
template <typename T> class TensorValue;
template <typename T> class VectorValue;
typedef VectorValue<Real> RealVectorValue;
typedef TensorValue<Real> RealTensorValue;
typedef RealVectorValue RealGradient;
typedef RealTensorValue RealTensor;


/**
 * This class uses the error estimate given by different types of
 * derefinement (h coarsening or p reduction) to choose between h
 * refining and p elevation.
 * Currently we assume that a set of elements has already been flagged
 * for h refinement, and we may want to change some of those elements
 * to be flagged for p refinement.
 *
 * This code is currently experimental and will not produce optimal
 * hp meshes without significant improvement.
 *
 * @author Roy H. Stogner, 2006.
 */
class HPCoarsenTest : public HPSelector
{
public:

  /**
   * Constructor.
   */
  HPCoarsenTest() : p_weight(1.0)
  {
    libmesh_experimental();
  }

  /**
   * Destructor.
   */
  virtual ~HPCoarsenTest() {}


  /**
   * This pure virtual function must be redefined
   * in derived classes to take a mesh flagged for h
   * refinement and potentially change the desired
   * refinement type.
   */
  virtual void select_refinement (System& system);

  /**
   * Because the coarsening test seems to always choose p refinement, we're
   * providing an option to make h refinement more likely
   */
  Real p_weight;

protected:
  /**
   * The helper function which adds individual fine element data to
   * the coarse element projection
   */
  void add_projection(const System&, const Elem*, unsigned int var);

  /**
   * The coarse element on which a solution projection is cached
   */
  const Elem *coarse;

  /**
   * Global DOF indices for fine elements
   */
  std::vector<unsigned int> dof_indices;

  /**
   * The finite element objects for fine and coarse elements
   */
  AutoPtr<FEBase> fe, fe_coarse;

  /**
   * The shape functions and their derivatives
   */
  const std::vector<std::vector<Real> > *phi, *phi_coarse;
  const std::vector<std::vector<RealGradient> > *dphi, *dphi_coarse;
  const std::vector<std::vector<RealTensor> > *d2phi, *d2phi_coarse;

  /**
   * Mapping jacobians
   */
  const std::vector<Real> *JxW;

  /**
   * Quadrature locations
   */
  const std::vector<Point> *xyz_values;
  std::vector<Point> coarse_qpoints;

  /**
   * The quadrature rule for the fine element
   */
  AutoPtr<QBase> qrule;

  /**
   * Linear system for projections
   */
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;
  /**
   * Coefficients for projected coarse and projected
   * p-derefined solutions
   */
  DenseVector<Number> Uc;
  DenseVector<Number> Up;
};

} // namespace libMesh

#endif // #ifdef LIBMESH_ENABLE_AMR

#endif // LIBMESH_HP_COARSENTEST_H

