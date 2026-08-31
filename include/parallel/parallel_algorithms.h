
// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_PARALLEL_ALGORITHMS_H
#define LIBMESH_PARALLEL_ALGORITHMS_H

// The library configuration options
#include "libmesh/libmesh_config.h"

// C++ headers
#include <numeric>

// Workaround incomplete C++17 support in some compilers/libs

#ifdef LIBMESH_HAVE_CXX17_TRANSFORM_REDUCE
template <typename InputIt, class T, class BinaryOp, class UnaryOp>
T libmesh_transform_reduce(InputIt begin, InputIt end, T init, BinaryOp reduce, UnaryOp transform)
{
  return std::transform_reduce(begin, end, init, reduce, transform);
}
#else
template <typename InputIt, class T, class BinaryOp, class UnaryOp>
T libmesh_transform_reduce(InputIt begin, InputIt end, T init, BinaryOp reduce, UnaryOp transform)
{
  // Reuse init as returnval.
  //
  // Don't try to do any fancy reduce() reordering/vectorization;
  // people who want that can get a real compiler.
  for (auto it = begin; it != end; ++it)
    init = reduce(init, transform(*it));
  return init;
}
#endif // LIBMESH_HAVE_CXX17_TRANSFORM_REDUCE

#endif // LIBMESH_PARALLEL_ALGORITHMS_H
