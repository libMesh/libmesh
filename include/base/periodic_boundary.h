// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_PERIODIC_BOUNDARY_H
#define LIBMESH_PERIODIC_BOUNDARY_H

// Local Includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_PERIODIC

// Local Includes
#include "libmesh/periodic_boundary_base.h"
#include "libmesh/vector_value.h" // RealVectorValue
#include "libmesh/auto_ptr.h" // libmesh_make_unique
#include "libmesh/libmesh.h" // libMesh::invalid_uinit

namespace libMesh
{

// Forward Declarations
template <typename> class ElemTempl;
template <typename> class MeshBaseTempl;

/**
 * The definition of a periodic boundary.
 *
 * \author Roy Stogner
 * \date 2010
 * \brief Used for implementing periodic BCs via constraints.
 */
template <typename RealType = Real>
class PeriodicBoundaryTempl : public PeriodicBoundaryBaseTempl<RealType>
{
public:
  typedef ElemTempl<RealType> Elem;
  typedef PointTempl<RealType> Point;
  typedef MeshBaseTempl<RealType> MeshBase;
  typedef PeriodicBoundaryTempl<RealType> PeriodicBoundary;
  typedef PeriodicBoundaryBaseTempl<RealType> PeriodicBoundaryBase;

  /**
   * Constructor
   */
  PeriodicBoundaryTempl();

  /**
   * Destructor
   */
  virtual ~PeriodicBoundaryTempl() {}

  /**
   * Copy constructor, with option for the copy to represent an inverse transformation.
   */
  PeriodicBoundaryTempl(const PeriodicBoundary & o, TransformationType t = FORWARD);

  /**
   * Constructor taking a reference to the translation vector.
   */
  PeriodicBoundaryTempl(const RealVectorValue & vector);

  /**
   * This function should be overridden by derived classes to
   * define how one finds corresponding nodes on the periodic
   * boundary pair.
   */
  virtual Point get_corresponding_pos(const Point & pt) const override;

  /**
   * If we want the DofMap to be able to make copies of references and
   * store them in the underlying map, this class must be clone'able,
   * i.e. have a kind of virtual construction mechanism.
   */
  virtual std::unique_ptr<PeriodicBoundaryBase> clone(TransformationType t = FORWARD) const override;

protected:

  // The vector which is added to points in myboundary
  // to produce corresponding points in pairedboundary
  RealVectorValue translation_vector;
};

template <typename RealType>
PeriodicBoundaryTempl<RealType>::PeriodicBoundaryTempl() :
  PeriodicBoundaryBase(),
  translation_vector()
{
}



template <typename RealType>
PeriodicBoundaryTempl<RealType>::PeriodicBoundaryTempl(const PeriodicBoundary & o, TransformationType t) :
  PeriodicBoundaryBase(o),
  translation_vector(o.translation_vector)
{
  if (t == INVERSE)
    {
      std::swap(myboundary, pairedboundary);
      translation_vector *= -1.0;
    }
}



template <typename RealType>
PeriodicBoundaryTempl<RealType>::PeriodicBoundaryTempl(const RealVectorValue & vector) :
  PeriodicBoundaryBase(),
  translation_vector(vector)
{
}



template <typename RealType>
PointTempl<RealType> PeriodicBoundaryTempl<RealType>::get_corresponding_pos(const Point & pt) const
{
  return pt + translation_vector;
}



template <typename RealType>
std::unique_ptr<PeriodicBoundaryBaseTempl<RealType>> PeriodicBoundaryTempl<RealType>::clone(TransformationType t) const
{
  return libmesh_make_unique<PeriodicBoundary>(*this, t);
}

typedef PeriodicBoundaryTempl<Real> PeriodicBoundary;

} // namespace libmesh

#endif // LIBMESH_ENABLE_PERIODIC

#endif // LIBMESH_PERIODIC_BOUNDARY_H
