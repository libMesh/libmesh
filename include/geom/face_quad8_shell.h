// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_FACE_QUAD8_SHELL_H
#define LIBMESH_FACE_QUAD8_SHELL_H

#include "face_quad8.h"

namespace libMesh
{

/**
 * QuadShell8 is almost identical to Quad8. The only difference is
 * with the type of boundary data we store for this case. We need this
 * "stub" class in order to differentiate between this class and other
 * classes when reading/writing Mesh files.
 *
 * \author Sylvain Vallaghe
 * \date 2017
 * \brief A 2D quadrilateral shell element with 8 nodes.
 */
class QuadShell8 : public Quad8
{
public:
  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  QuadShell8 (Elem * p=libmesh_nullptr) :
    Quad8(p) {}

  QuadShell8 (QuadShell8 &&) = delete;
  QuadShell8 (const QuadShell8 &) = delete;
  QuadShell8 & operator= (const QuadShell8 &) = delete;
  QuadShell8 & operator= (QuadShell8 &&) = delete;
  virtual ~QuadShell8() = default;

  /**
   * \returns \p QUADSHELL8.
   */
  virtual ElemType type () const override { return QUADSHELL8; }
};

} // namespace libMesh

#endif
