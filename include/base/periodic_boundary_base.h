// The libMesh Finite Element Library.
// Copyright (C) 2002-2023 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_PERIODIC_BOUNDARY_BASE_H
#define LIBMESH_PERIODIC_BOUNDARY_BASE_H

// Local Includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_PERIODIC

// Local Includes
#include "libmesh/point.h"
#include "libmesh/dense_matrix.h"

// C++ Includes
#include <set>
#include <memory>

namespace libMesh
{

// Forward Declarations
class Elem;
class MeshBase;

/**
 * The type of enforcement of condition; strong(default) = constraint matrix
 */
enum class EnforcementType
  { STRONG_ENFORCEMENT=0,
    WEAK_ENFORCEMENT=1 };

/**
 * The base class for defining periodic boundaries.
 *
 * \author Roy Stogner
 * \date 2010
 * \brief Base class for all PeriodicBoundary implementations.
 */
class PeriodicBoundaryBase
{
public:
  enum TransformationType
    { FORWARD=0,
      INVERSE=1 };

  /**
   * The boundary ID of this boundary and its counterpart
   */
  boundary_id_type myboundary, pairedboundary;

  /**
   * Constructor
   */
  PeriodicBoundaryBase();

  /**
   * Copy constructor
   */
  PeriodicBoundaryBase(const PeriodicBoundaryBase & other);

  /**
   * Destructor
   */
  virtual ~PeriodicBoundaryBase() = default;

  /**
   * This function should be overridden by derived classes to
   * define how one finds corresponding nodes on the periodic
   * boundary pair.
   */
  virtual Point get_corresponding_pos(const Point & pt) const = 0;

  /**
   * If we want the DofMap to be able to make copies of references and
   * store them in the underlying map, this class must be clone'able,
   * i.e. have a kind of virtual construction mechanism.  The user can
   * also pass a flag to enable an 'inverse transformation' to be
   * cloned from a forward transformation. The simplest way to
   * implement a clone function like this is in terms of a copy
   * constructor, see periodic_boundary.h.
   *
   * \note Not every transformation needs to provide an automatic way
   * to clone an inverse: you can simply add a pair of
   * PeriodicBoundaryBase objects using the appropriate DofMap
   * interface instead.
   */
  virtual std::unique_ptr<PeriodicBoundaryBase> clone(TransformationType t = FORWARD) const = 0;

  void set_variable(unsigned int var);

  void merge(const PeriodicBoundaryBase & pb);

  bool is_my_variable(unsigned int var_num) const;

  /**
   * @return true if _transformation_matrix is not null.
   */
  bool has_transformation_matrix() const;

  /**
   * Get the transformation matrix, if it is defined.
   * Throw an error if it is not defined.
   */
  const DenseMatrix<Real> & get_transformation_matrix() const;

  /**
   * Set the transformation matrix. When calling this method we
   * require the following conditions:
   *  1) \p matrix is square with size that matches this->variables.size()
   *  2) the list of variables in this->variables set must all have the same FE type
   * Both of these conditions are asserted in DBG mode.
   */
  void set_transformation_matrix(const DenseMatrix<Real> & matrix);

  /**
   * Get the set of variables for this periodic boundary condition.
   */
  const std::set<unsigned int> & get_variables() const;

  /**
   * Get the enforcement type.
   */
  const EnforcementType & get_enforcement_type() const;

  /**
   * Set the enforcement type.
   */
  void set_enforcement_type(const EnforcementType & e_type);

protected:

  /**
   * Set of variables for this periodic boundary, empty means all variables possible
   */
  std::set<unsigned int> variables;

  /**
   * A DenseMatrix that defines the mapping of variables on this
   * boundary and the counterpart boundary. This is necessary for
   * periodic-boundaries with vector-valued quantities (e.g.
   * velocity or displacement) on a sector of a circular domain,
   * for example, since in that case we must map each variable
   * to a corresponding linear combination of all the variables.
   * We store the DenseMatrix via a unique_ptr, and an uninitialized
   * pointer is treated as equivalent to the identity matrix.
   */
  std::unique_ptr<DenseMatrix<Real>> _transformation_matrix;

  EnforcementType _enforcement_type;

};

} // namespace libmesh

#endif // LIBMESH_ENABLE_PERIODIC

#endif // LIBMESH_PERIODIC_BOUNDARY_BASE_H
