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


#ifndef LIBMESH_PARALLEL_ONLY_H
#define LIBMESH_PARALLEL_ONLY_H

// Macro to identify and debug functions which should only be called in
// parallel on every processor at once

#ifdef LIBMESH_ENABLE_DEPRECATED
#undef parallel_only
#ifndef NDEBUG
#define parallel_only() do {                                            \
    libmesh_deprecated();                                               \
    libmesh_assert(CommWorld.verify(std::string(__FILE__).size()));     \
    libmesh_assert(CommWorld.verify(std::string(__FILE__)));            \
    libmesh_assert(CommWorld.verify(__LINE__)); } while (0)
#else
#define parallel_only()  ((void) 0)
#endif
#endif

#undef libmesh_parallel_only
#ifndef NDEBUG
#define libmesh_parallel_only(comm_obj) do {                            \
    libmesh_assert((comm_obj).verify(std::string(__FILE__).size()));    \
    libmesh_assert((comm_obj).verify(std::string(__FILE__)));           \
    libmesh_assert((comm_obj).verify(__LINE__)); } while (0)
#else
#define libmesh_parallel_only(comm_obj)  ((void) 0)
#endif

// Macro to identify and debug functions which should only be called in
// parallel on every processor at once

#ifdef LIBMESH_ENABLE_DEPRECATED
#undef parallel_only_on
#ifndef NDEBUG
#define parallel_only_on(comm_arg) do {                                 \
    libmesh_deprecated();                                               \
    libmesh_assert(CommWorld.verify(std::string(__FILE__).size(), comm_arg)); \
    libmesh_assert(CommWorld.verify(std::string(__FILE__), comm_arg));  \
    libmesh_assert(CommWorld.verify(__LINE__), comm_arg); } while (0)
#else
#define parallel_only_on(comm_arg)  ((void) 0)
#endif
#endif

#undef libmesh_parallel_only_on
#ifndef NDEBUG
#define libmesh_parallel_only_on(comm_obj,comm_arg) do {                \
    libmesh_assert(comm_obj.verify(std::string(__FILE__).size(), comm_arg)); \
    libmesh_assert(comm_obj.verify(std::string(__FILE__), comm_arg));   \
    libmesh_assert(comm_obj.verify(__LINE__), comm_arg); } while (0)
#else
#define libmesh_parallel_only_on(comm_obj,comm_arg)  ((void) 0)
#endif

#endif // LIBMESH_PARALLEL_ONLY_H
