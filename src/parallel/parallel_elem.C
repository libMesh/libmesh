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

// Local includes
#include "libmesh/parallel_elem_impl.h"

namespace libMesh
{
template class Parallel::Packing<const ElemTempl<Real> *>;
template class Parallel::Packing<ElemTempl<Real> *>;
template unsigned int Parallel::Packing<const ElemTempl<Real> *>::packed_size(
    typename std::vector<largest_id_type>::const_iterator);
template unsigned int Parallel::Packing<const ElemTempl<Real> *>::packed_size(
    typename std::vector<largest_id_type>::iterator);
template unsigned int
Parallel::Packing<const ElemTempl<Real> *>::packable_size(const ElemTempl<Real> * const &,
                                                          const MeshBaseTempl<Real> *);
template unsigned int
Parallel::Packing<const ElemTempl<Real> *>::packable_size(const ElemTempl<Real> * const &,
                                                          const DistributedMeshTempl<Real> *);
template unsigned int
Parallel::Packing<const ElemTempl<Real> *>::packable_size(const ElemTempl<Real> * const &,
                                                          const ParallelMeshTempl<Real> *);
template void Parallel::Packing<const ElemTempl<Real> *>::pack(
    const ElemTempl<Real> * const & elem,
    std::back_insert_iterator<std::vector<largest_id_type>> data_out,
    const MeshBaseTempl<Real> * mesh);
template void Parallel::Packing<const ElemTempl<Real> *>::pack(
    const ElemTempl<Real> * const & elem,
    std::back_insert_iterator<std::vector<largest_id_type>> data_out,
    const DistributedMeshTempl<Real> * mesh);
template void Parallel::Packing<const ElemTempl<Real> *>::pack(
    const ElemTempl<Real> * const & elem,
    std::back_insert_iterator<std::vector<largest_id_type>> data_out,
    const ParallelMeshTempl<Real> * mesh);
template ElemTempl<Real> *
Parallel::Packing<ElemTempl<Real> *>::unpack(std::vector<largest_id_type>::const_iterator in,
                                             MeshBaseTempl<Real> * mesh);
template ElemTempl<Real> *
Parallel::Packing<ElemTempl<Real> *>::unpack(std::vector<largest_id_type>::const_iterator in,
                                             DistributedMeshTempl<Real> * mesh);
template ElemTempl<Real> *
Parallel::Packing<ElemTempl<Real> *>::unpack(std::vector<largest_id_type>::const_iterator in,
                                             ParallelMeshTempl<Real> * mesh);
}
