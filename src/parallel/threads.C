// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// System Includes

// Local Includes
#include "libmesh/threads.h"

#if LIBMESH_HAVE_OPENMP
#include <omp.h>
#endif

namespace libMesh
{

#if !defined(LIBMESH_HAVE_TBB_API) && defined(LIBMESH_HAVE_PTHREAD)
  std::map<pthread_t, unsigned int> Threads::_pthread_unique_ids;
  Threads::spin_mutex Threads::_pthread_unique_id_mutex;

  unsigned int Threads::pthread_unique_id()
  {
#if LIBMESH_HAVE_OPENMP
    return omp_get_thread_num();
#else
    spin_mutex::scoped_lock lock(_pthread_unique_id_mutex);
    return _pthread_unique_ids[pthread_self()];
#endif
  }
#endif

//-------------------------------------------------------------------------
// Threads:: object instantiation
Threads::spin_mutex Threads::spin_mtx;
Threads::recursive_mutex Threads::recursive_mtx;
bool Threads::in_threads = false;


} // namespace libMesh
