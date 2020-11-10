// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local includes
#include "libmesh/remote_elem.h"
#include "libmesh/libmesh_singleton.h"
#include "libmesh/threads.h"



namespace
{
using namespace libMesh;

// Class to be dispatched by Singleton::setup()
// to create the \p RemoteElem singleton.
// While this actual object has file-level static
// scope and will be initialized before main(),
// importantly the setup() method will not be invoked
// until after main().
class RemoteElemSetup : public Singleton::Setup
{
  void setup ()
  {
    RemoteElem::create();
  }
} remote_elem_setup;
}



namespace libMesh
{

// Pointer to singleton Remote Element (to be created in
// libMesh::init()
const RemoteElem * remote_elem;


RemoteElem::~RemoteElem()
{
  remote_elem = nullptr;
}



const Elem & RemoteElem::create ()
{
  if (remote_elem == nullptr)
    remote_elem = new RemoteElem;

  return *remote_elem;
}


} // namespace libMesh
