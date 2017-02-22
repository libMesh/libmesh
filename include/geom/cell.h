// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
class Cell : public Elem
{
public:

  /**
   * Constructor.
   */
  Cell (const unsigned int nn,
        const unsigned int ns,
        Elem * p,
        Elem ** elemlinkdata,
        Node ** nodelinkdata) :
    Elem (nn, ns, p, elemlinkdata, nodelinkdata) {}

  /**
   * @returns 3, the dimensionality of the object.
   */
  virtual unsigned int dim () const libmesh_override { return 3; }

  /**
   * @return a bounding box (not necessarily the minimal bounding box)
   * containing the geometric element.
   */
  virtual BoundingBox loose_bounding_box () const libmesh_override;

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * @returns \p false.  All classes derived from \p Cell
   * are finite elements.
   */
  virtual bool infinite () const libmesh_override { return false; }

#endif
};


} // namespace libMesh


#endif // LIBMESH_CELL_H
