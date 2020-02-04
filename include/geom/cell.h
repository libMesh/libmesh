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



#ifndef LIBMESH_CELL_H
#define LIBMESH_CELL_H

// Local includes
#include "libmesh/elem.h"

namespace libMesh
{

/**
 * The \p Cell is an abstract element type that lives in
 * three dimensions.  A cell could be a tetrahedron, a hexahedron,
 * a pyramid, a prism, etc...
 *
 * \author Benjamin S. Kirk
 * \date 2002
 * \brief The base class for all 3D geometric element types.
 */
template <typename RealType = Real>
class CellTempl : public ElemTempl<RealType>
{
public:
  typedef CellTempl<RealType> Cell;
  typedef BoundingBoxTempl<RealType> BoundingBox;
  typedef ElemTempl<RealType> Elem;

  /**
   * Constructor.
   */
  CellTempl (const unsigned int nn,
             const unsigned int ns,
             Elem * p,
             Elem ** elemlinkdata,
             Node ** nodelinkdata) :
    Elem (nn, ns, p, elemlinkdata, nodelinkdata) {}

  CellTempl (Cell &&) = delete;
  CellTempl (const Cell &) = delete;
  Cell & operator= (const Cell &) = delete;
  Cell & operator= (Cell &&) = delete;
  virtual ~CellTempl() = default;

  /**
   * \returns 3, the dimensionality of the object.
   */
  virtual unsigned short dim () const override { return 3; }

  /**
   * \returns A bounding box (not necessarily the minimal bounding box)
   * containing the geometric element.
   */
  virtual BoundingBox loose_bounding_box () const override;

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * \returns \p false.  All classes derived from \p Cell
   * are finite elements.
   */
  virtual bool infinite () const override { return false; }

#endif
};

template <typename RealType>
BoundingBoxTempl<RealType> CellTempl<RealType>::loose_bounding_box () const
{
  // This might have curved sides, but it's definitely *not* curving
  // through 4-D space, so the full bounding box is just the merger of
  // the sides' bounding boxes.

  std::unique_ptr<const Elem> side_ptr { this->build_side_ptr(0, false) };
  BoundingBox bbox = side_ptr->loose_bounding_box();
  unsigned int my_n_sides = this->n_sides();
  for (unsigned s=1; s < my_n_sides; ++s)
    {
      this->build_side_ptr(side_ptr, s);
      bbox.union_with(side_ptr->loose_bounding_box());
    }

  return bbox;
}

typedef CellTempl<Real> Cell;

} // namespace libMesh


#endif // LIBMESH_CELL_H
