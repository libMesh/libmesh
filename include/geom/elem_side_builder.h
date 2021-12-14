// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_ELEM_SIDE_BUILDER_H
#define LIBMESH_ELEM_SIDE_BUILDER_H

#include "libmesh/elem.h"

namespace libMesh
{

/**
 * Helper for building element sides that minimizes the construction
 * of new elements.
 *
 * The building of element side pointers has the option to pass an
 * existing constructed element to avoid extraneous allocation.
 * If the element that is passed is of the correct type, the
 * build is done in place, that is, the nodes and subdomain are
 * done in place.
 *
 * This tool contains a cache with an entry per element type.
 * When a side is requested, the cached element of the desired
 * type is utilized. On first request of a side of a specific
 * type, the desired element is constructed. On all subsequent
 * calls for a side of the same type, the cached element is
 * used and the necessary members (nodes and subdomain)
 * are changed in place.
 *
 * NOTE: This tool is meant for on-the-fly use. Because
 * the cache is changed on each call, the references obtained
 * from previous calls should be considered invalid.
 */
class ElemSideBuilder
{
public:
  /**
   * \returns an element side for side \p s of element \p elem.
   */
  ///@{
  Elem & operator()(Elem & elem, const unsigned int s);
  const Elem & operator()(const Elem & elem, const unsigned int s);
  ///@}

private:
  /// Element cache for building sides; indexed by ElemType
  std::vector<std::unique_ptr<Elem>> _cached_elems;
};

} // namespace libMesh

#endif // LIBMESH_ELEM_SIDE_BUILDER_H
