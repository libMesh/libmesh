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

#ifndef LIBMESH_FACE_QUAD4_SHELL_H
#define LIBMESH_FACE_QUAD4_SHELL_H

#include "face_quad4.h"

namespace libMesh
{

/**
 * QuadShell4 is almost identical to Quad4. The only difference is
 * with the type of boundary data we store for this case. We need this
 * "stub" class in order to differentiate between this class and other
 * classes when reading/writing Mesh files.
 *
 * \author David Knezevic
 * \date 2016
 * \brief A 2D quadrilateral shell element with 4 nodes.
 */
class QuadShell4 : public Quad4
{
public:
  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  QuadShell4 (Elem * p=nullptr) :
    Quad4(p) {}

  QuadShell4 (QuadShell4 &&) = delete;
  QuadShell4 (const QuadShell4 &) = delete;
  QuadShell4 & operator= (const QuadShell4 &) = delete;
  QuadShell4 & operator= (QuadShell4 &&) = delete;
  virtual ~QuadShell4() = default;

  /**
   * \returns \p QUADSHELL4.
   */
  virtual ElemType type () const override { return QUADSHELL4; }
};

} // namespace libMesh

#endif
