// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/raw_type.h"

namespace libMesh
{

/**
 * A \p Point defines a location in LIBMESH_DIM dimensional Real space.  Points
 * are always real-valued, even if the library is configured with
 * \p --enable-complex.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 * \brief A geometric point in (x,y,z) space.
 */
class Point : public VectorValue<GeomReal>
{
public:

  /**
   * Constructor.  By default sets all entries to 0.  Gives the point
   * 0 in \p LIBMESH_DIM dimensions.
   */
  Point (const GeomReal x=0.,
         const GeomReal y=0.,
         const GeomReal z=0.) :
    VectorValue<GeomReal> (x,y,z)
  {}

  /**
   * Trivial copy-constructor.
   */
  Point (const Point & p) = default;

  /**
   * Copy-constructor.
   */
  template <typename T>
  Point (const VectorValue<T> & p) :
    VectorValue<GeomReal> (p)
  {}

  /**
   * Copy-constructor from non-point Typevector.
   */
  template <typename T>
  Point (const TypeVector<T> & p) :
    VectorValue<GeomReal> (p)
  {}

  /**
   * Copy-assignment operator.
   */
  Point& operator=(const Point & p) = default;

  /**
   * Disambiguate constructing from non-Real scalars
   */
  template <typename T,
            typename = typename
              boostcopy::enable_if_c<ScalarTraits<T>::value,void>::type>
  Point (const T x) :
    VectorValue<GeomReal> (x,0,0)
  {}

  /**
   * Destructor.
   */
  ~Point() = default;

protected:

  /**
   * Make the derived class a friend.
   */
  friend class Node;
};

} // namespace libMesh

namespace MetaPhysicL
{
template <bool>
struct PointRawValueHelper;

template <>
struct PointRawValueHelper<true>
{
  typedef const libMesh::Point & value_type;
};

template <>
struct PointRawValueHelper<false>
{
  typedef libMesh::VectorValue<typename RawType<libMesh::GeomReal>::value_type> value_type;
};

template <>
struct RawType<libMesh::Point>
{
  typedef typename PointRawValueHelper<IsRawSame<libMesh::GeomReal>::value>::value_type value_type;

  template <typename T = libMesh::GeomReal,
            typename std::enable_if<!IsRawSame<T>::value, int>::type = 0>
  static value_type value (const libMesh::Point & in)
    {
      value_type ret;
      for (unsigned int i = 0; i < LIBMESH_DIM; ++i)
        ret(i) = raw_value(in(i));

      return ret;
    }

  template <typename T = libMesh::GeomReal,
            typename std::enable_if<IsRawSame<T>::value, int>::type = 0>
  static value_type value (const libMesh::Point & in)
    {
      return in;
    }
};
}


#endif // LIBMESH_POINT_H
