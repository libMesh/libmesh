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

#ifndef LIBMESH_KOKKOS_STORAGE_H
#define LIBMESH_KOKKOS_STORAGE_H

// libMesh Kokkos storage helpers for dimension-aware vector/tensor views.
//
// The functions here convert between the two representations used by the
// Kokkos numerics layer:
//
// - "storage": a Kokkos View holding many vectors or tensors, with the
//   leading (runtime) extent indexing the entry and the trailing
//   LIBMESH_DIM extents indexing components, and
// - "values": individual owning objects (TypeVector, TypeTensor, Point,
//   ...) that kernels compute with in registers.
//
// load_*() reads one entry out of a view into an owning value;
// store_*() writes an owning value back into one entry of a view.  Both
// are usable inside device kernels.

#include "libmesh/kokkos_linalg_base.h"

#include "libmesh/libmesh_common.h"
#include "libmesh/libmesh_device.h"
#include "libmesh/type_tensor.h"
#include "libmesh/type_vector.h"

namespace libMesh::Kokkos
{

/**
 * \returns Entry \p i of the vector view \p view, materialized as an
 * owning vector of type \p VectorType (e.g. a TypeVector, VectorValue,
 * or Point).
 *
 * The view is expected to have one runtime extent indexing entries
 * followed by a LIBMESH_DIM extent indexing components, as created by
 * make_vector_storage().
 */
template <typename VectorType, typename ViewType>
LIBMESH_DEVICE_INLINE
VectorType load_vector(const ViewType & view, const unsigned int i)
{
  return materialize_vector<VectorType>(make_vector_ref(view, i));
}

/**
 * Copies the components of the vector \p v into entry \p i of the
 * vector view \p view, overwriting whatever that entry held.
 */
template <typename ViewType, typename VectorType>
LIBMESH_DEVICE_INLINE
void store_vector(const ViewType & view, const unsigned int i, const VectorType & v)
{
  auto out = make_vector_ref(view, i);

  for (unsigned int d = 0; d < LIBMESH_DIM; ++d)
    out(d) = v(d);
}

/**
 * \returns Entry \p i of the tensor view \p view, materialized as an
 * owning tensor of type \p TensorType (e.g. a TypeTensor or
 * TensorValue).
 *
 * The view is expected to have one runtime extent indexing entries
 * followed by LIBMESH_DIM x LIBMESH_DIM extents indexing components,
 * as created by make_tensor_storage().
 */
template <typename TensorType, typename ViewType>
LIBMESH_DEVICE_INLINE
TensorType load_tensor(const ViewType & view, const unsigned int i)
{
  return materialize_tensor<TensorType>(make_tensor_ref(view, i));
}

/**
 * Copies the components of the tensor \p T into entry \p i of the
 * tensor view \p view, overwriting whatever that entry held.
 */
template <typename ViewType, typename TensorType>
LIBMESH_DEVICE_INLINE
void store_tensor(const ViewType & view, const unsigned int i, const TensorType & T)
{
  auto out = make_tensor_ref(view, i);

  for (unsigned int row = 0; row < LIBMESH_DIM; ++row)
    for (unsigned int col = 0; col < LIBMESH_DIM; ++col)
      out(row, col) = T(row, col);
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_STORAGE_H
