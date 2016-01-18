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



#ifndef LIBMESH_FACE_H
#define LIBMESH_FACE_H

// Local includes
#include "libmesh/elem.h"


namespace libMesh
{

// Forward declarations
class Mesh;



/**
 * The \p Face is an abstract element type that lives in
 * two dimensions.  A face could be a triangle, a quadrilateral,
 * a pentagon, etc...
 */
class Face : public Elem
{
public:

  /**
   * Constructor.  Explicitly specifies the number of
   * nodes and neighbors for which storage will be allocated.
   */
  Face (const unsigned int nn,
        const unsigned int ns,
        Elem * p,
        Elem ** elemlinkdata,
        Node ** nodelinkdata) :
    Elem(nn, ns, p, elemlinkdata, nodelinkdata) {}

  /**
   * @returns 2, the dimensionality of the object.
   */
  virtual unsigned int dim () const libmesh_override { return 2; }

  /**
   * @returns 0.  All 2D elements have no faces, just
   * edges.
   */
  virtual unsigned int n_faces() const libmesh_override { return 0; }

  /**
   * build_side and build_edge are identical for faces
   */
  virtual UniquePtr<Elem> build_edge (const unsigned int i) const libmesh_override
  { return build_side(i); }

  /*
   * is_edge_on_side is trivial in 2D
   */
  virtual bool is_edge_on_side(const unsigned int e,
                               const unsigned int s) const libmesh_override
  { return (e == s); }

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

  /**
   * @returns \p false.  All classes derived from \p Face
   * are finite elements.
   */
  virtual bool infinite () const libmesh_override { return false; }

#endif

};

} // namespace libMesh

#endif // LIBMESH_FACE_H
