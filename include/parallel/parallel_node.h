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


#ifndef LIBMESH_PARALLEL_NODE_H
#define LIBMESH_PARALLEL_NODE_H

// libMesh includes
#include "libmesh/id_types.h"
#include "libmesh/libmesh_config.h"

// TIMPI includes
#include "timpi/packing.h"


namespace libMesh
{

// Forward declarations
template <typename> class NodeTempl;
typedef NodeTempl<Real> Node;

namespace Parallel
{

template <>
class Packing<const Node *>
{
public:
  typedef largest_id_type buffer_type;

  template <typename OutputIter, typename Context>
  static void pack(const Node * const & object,
                   OutputIter data_out,
                   const Context * context);

  template <typename Context>
  static unsigned int packable_size(const Node * const & object,
                                    const Context * context);

  template <typename BufferIter>
  static unsigned int packed_size(BufferIter iter);

  template <typename BufferIter, typename Context>
  static const Node * unpack(BufferIter in, Context * ctx);
};


template <>
class Packing<Node *>
{
public:
  typedef largest_id_type buffer_type;

  template <typename OutputIter, typename Context>
  static void pack(Node * const & object,
                   OutputIter data_out,
                   const Context * context)
  { return Packing<const Node *>::pack(object, data_out, context); }

  template <typename Context>
  static unsigned int packable_size(Node * const & object,
                                    const Context * context)
  { return Packing<const Node*>::packable_size(object, context); }

  template <typename BufferIter>
  static unsigned int packed_size(BufferIter iter)
  { return Packing<const Node *>::packed_size(iter); }

  template <typename BufferIter, typename Context>
  static Node * unpack(BufferIter in, Context * ctx);
};


template <typename BufferIter, typename Context>
inline const Node *
Packing<const Node *>::unpack(BufferIter in, Context * ctx)
{ return Packing<Node *>::unpack(in, ctx); }



} // namespace Parallel
} // namespace libMesh

#endif // LIBMESH_PARALLEL_NODE_H
