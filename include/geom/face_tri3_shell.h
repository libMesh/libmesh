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

#ifndef LIBMESH_FACE_TRI3_SHELL_H
#define LIBMESH_FACE_TRI3_SHELL_H

#include "face_tri3.h"

namespace libMesh
{

/**
 * TriShell3 is almost identical to Tri3. The only difference is with
 * the type of boundary data we store for this case.  We need this
 * "stub" class in order to differentiate between this class and other
 * classes when reading/writing Mesh files.
 *
 * \author David Knezevic
 * \date 2016
 */
class TriShell3 : public Tri3
{
public:
  /**
   * Constructor.  By default this element has no parent.
   */
  explicit
  TriShell3 (Elem * p=nullptr) :
    Tri3(p) {}

  TriShell3 (TriShell3 &&) = delete;
  TriShell3 (const TriShell3 &) = delete;
  TriShell3 & operator= (const TriShell3 &) = delete;
  TriShell3 & operator= (TriShell3 &&) = delete;
  virtual ~TriShell3() = default;

  /**
   * \returns \p TRISHELL3.
   */
  virtual ElemType type () const override { return TRISHELL3; }
};

} // namespace libMesh

#endif
