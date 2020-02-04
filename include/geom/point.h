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



#ifndef LIBMESH_POINT_H
#define LIBMESH_POINT_H

// Local includes
#include "libmesh/vector_value.h"

namespace libMesh
{

template <typename> class NodeTempl;

/**
 * A \p Point defines a location in LIBMESH_DIM dimensional Real space.  Points
 * are always real-valued, even if the library is configured with
 * \p --enable-complex.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 * \brief A geometric point in (x,y,z) space.
 */
template <typename RealType>
class PointTempl : public VectorValue<RealType>
{
public:

  /**
   * Constructor.  By default sets all entries to 0.  Gives the point
   * 0 in \p LIBMESH_DIM dimensions.
   */
  PointTempl (const RealType x=0.,
         const RealType y=0.,
         const RealType z=0.) :
    VectorValue<RealType> (x,y,z)
  {}

  /**
   * Copy-constructor.
   */
  PointTempl (const PointTempl & p) :
    VectorValue<RealType> (p)
  {}

  /**
   * Copy-constructor.
   */
  PointTempl (const TypeVector<RealType> & p) :
    VectorValue<RealType> (p)
  {}

  /**
   * Copy from other point type constructor
   */
  template <typename RealType2>
  PointTempl (const PointTempl<RealType2> & p2) :
      VectorValue<RealType>(p2) {}

  /**
   * Copy-assignment operator.
   */
  PointTempl& operator=(const PointTempl & p) = default;

  /**
   * Copy-assignment from other point type
   */
  template <typename RealType2>
  PointTempl& operator=(const PointTempl<RealType2> & p2)
    {
      this->assign(p2);
    }

  /**
   * Disambiguate constructing from non-Real scalars
   */
  template <typename T,
            typename = typename
              boostcopy::enable_if_c<ScalarTraits<T>::value,void>::type>
  PointTempl (const T x) :
    VectorValue<RealType> (x,0,0)
  {}

  /**
   * Empty.
   */
  ~PointTempl() {}

protected:

  /**
   * Make the derived class a friend.
   */
  friend class NodeTempl<RealType>;
};

typedef PointTempl<Real> Point;

} // namespace libMesh

#endif // LIBMESH_POINT_H
