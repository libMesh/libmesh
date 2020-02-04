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


#ifndef LIBMESH_PARALLEL_ELEM_H
#define LIBMESH_PARALLEL_ELEM_H

// libMesh includes
#include "libmesh/id_types.h"
#include "libmesh/libmesh_config.h"

// TIMPI includes
#include "timpi/packing.h"


namespace libMesh
{

// Forward declarations
template <typename> class ElemTempl;
template <typename> class MeshBaseTempl;
template <typename> class RemoteElemTempl;

namespace Parallel
{

template <typename RealType>
class Packing<const ElemTempl<RealType> *>
{
public:
  typedef largest_id_type buffer_type;
  typedef ElemTempl<RealType> Elem;
  typedef MeshBaseTempl<RealType> MeshBase;
  typedef RemoteElemTempl<RealType> RemoteElem;

  template <typename OutputIter, typename Context>
  static void pack(const Elem * const & object,
                   OutputIter data_out,
                   const Context * context);

  template <typename Context>
  static unsigned int packable_size(const Elem * const & object,
                                    const Context * context);

  template <typename BufferIter>
  static unsigned int packed_size(BufferIter iter);

  template <typename BufferIter, typename Context>
  static const Elem * unpack(BufferIter in, Context * ctx);
};


template <typename RealType>
class Packing<ElemTempl<RealType> *>
{
public:
  typedef largest_id_type buffer_type;
  typedef ElemTempl<RealType> Elem;
  typedef MeshBaseTempl<RealType> MeshBase;
  typedef RemoteElemTempl<RealType> RemoteElem;

  template <typename OutputIter, typename Context>
  static void pack(Elem * const & object,
                   OutputIter data_out,
                   const Context * context)
  { return Packing<const Elem *>::pack(object, data_out, context); }

  template <typename Context>
  static unsigned int packable_size(Elem * const & object,
                                    const Context * context)
  { return Packing<const Elem *>::packable_size(object, context); }

  template <typename BufferIter>
  static unsigned int packed_size(BufferIter iter)
  { return Packing<const Elem *>::packed_size(iter); }

  template <typename BufferIter, typename Context>
  static Elem * unpack(BufferIter in, Context * ctx);
};


template <typename RealType>
template <typename BufferIter, typename Context>
inline const ElemTempl<RealType>*
Packing<const ElemTempl<RealType>*>::unpack(BufferIter in, Context * ctx)
{ return Packing<ElemTempl<RealType>*>::unpack(in, ctx); }

} // namespace Parallel
} // namespace libMesh

#endif // LIBMESH_PARALLEL_ELEM_H
