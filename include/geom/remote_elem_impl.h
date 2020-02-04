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

#ifndef LIBMESH_REMOTE_ELEM_IMPL_H
#define LIBMESH_REMOTE_ELEM_IMPL_H

// Local includes
#include "libmesh/remote_elem.h"
#include "libmesh/libmesh_singleton.h"
#include "libmesh/threads.h"



namespace
{
using namespace libMesh;

typedef Threads::spin_mutex RemoteElemMutex;
RemoteElemMutex remote_elem_mtx;


// Class to be dispatched by Singleton::setup()
// to create the \p RemoteElem singleton.
// While this actual object has file-level static
// scope and will be initialized before main(),
// importantly the setup() method will not be invoked
// until after main().
template <typename RealType>
class RemoteElemSetup : public Singleton::Setup
{
  void setup () override
  {
    RemoteElemTempl<RealType>::create();
  }
};
}



namespace libMesh
{
template <typename RealType>
RemoteElemTempl<RealType>::~RemoteElemTempl()
{
  RemoteElemMutex::scoped_lock lock(remote_elem_mtx);

  remote_elem = nullptr;
}



template <typename RealType>
const ElemTempl<RealType> & RemoteElemTempl<RealType>::create ()
{
  if (remote_elem != nullptr)
    return *remote_elem;

  RemoteElemMutex::scoped_lock lock(remote_elem_mtx);

  // check again - object could have been created while waiting
  // for the lock to acquire!
  if (remote_elem == nullptr)
    remote_elem = new RemoteElemTempl<RealType>;

  return *remote_elem;
}


} // namespace libMesh

#endif // LIBMESH_REMOTE_ELEM_IMPL_H
