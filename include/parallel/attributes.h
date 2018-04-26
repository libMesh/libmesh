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


#ifndef LIBMESH_ATTRIBUTES_H
#define LIBMESH_ATTRIBUTES_H

namespace libMesh
{

namespace Parallel
{
//-------------------------------------------------------------------

/*
 * The unspecialized class gives default, lowest-common-denominator
 * attributes, for values which can't be used with Parallel min/max.
 * Specialized classes can set this to true, and should define
 * the lowest and highest values possible for the type.
 */
template<typename T>
struct Attributes
{
  static const bool has_min_max = false;
  static void set_lowest(T &) {}
  static void set_highest(T &) {}
};

} // namespace Parallel

} // namespace libMesh

#endif // LIBMESH_ATTRIBUTES_H
