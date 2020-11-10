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
#include "libmesh/libmesh_singleton.h"
#include "libmesh/threads.h"
#include "libmesh/simple_range.h"

// C/C++ includes
#include <vector>


// --------------------------------------------------------
// Local anonymous namespace to hold miscellaneous bits
namespace
{
using namespace libMesh;

// Mutex object for required locking
typedef Threads::spin_mutex SingletonMutex;
SingletonMutex singleton_mtx, setup_mtx;

// global list of runtime Singleton objects - created dynamically,
// cleaned up in reverse order.
typedef std::vector<Singleton *> SingletonList;

SingletonList & get_singleton_cache()
{
  static SingletonList singleton_cache;
  return singleton_cache;
}

typedef std::vector<Singleton::Setup *> SetupList;
SetupList & get_setup_cache()
{
  static SetupList setup_cache;
  return setup_cache;
}

} // end anonymous namespace



// --------------------------------------------------------
// Local anonymous namespace to hold miscellaneous bits
namespace libMesh
{

Singleton::Singleton ()
{
  SingletonMutex::scoped_lock lock(singleton_mtx);

  get_singleton_cache().push_back (this);
}



Singleton::Setup::Setup ()
{
  SingletonMutex::scoped_lock lock(setup_mtx);

  get_setup_cache().push_back (this);
}



void Singleton::setup ()
{
  SingletonMutex::scoped_lock lock(setup_mtx);

  SetupList & setup_cache = get_setup_cache();

  for (auto & item : setup_cache)
    {
      libmesh_assert (item);
      item->setup();
    }
}



void Singleton::cleanup ()
{
  SingletonMutex::scoped_lock lock(singleton_mtx);

  SingletonList & singleton_cache = get_singleton_cache();

  for (auto & item : as_range(singleton_cache.rbegin(),
                              singleton_cache.rend()))
    {
      libmesh_assert (item);
      delete item;
      item = nullptr;
    }

  singleton_cache.clear();
}



} // namespace libMesh
