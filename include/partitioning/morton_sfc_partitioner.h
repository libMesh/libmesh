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



#ifndef LIBMESH_MORTON_SFC_PARTITIONER_H
#define LIBMESH_MORTON_SFC_PARTITIONER_H

// Local Includes
#include "libmesh/sfc_partitioner.h"

namespace libMesh
{

/**
 * The \p MortonSFCPartitioner uses a Morton space filling curve to
 * partition the elements.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 * \brief Partitioner based on Morton's space filling curve algorithm.
 */
class MortonSFCPartitioner : public SFCPartitioner
{
public:

  /**
   * Constructor. Sets the underlying space filling curve type.
   */
  MortonSFCPartitioner ()
  {
    this->set_sfc_type ("Morton");
  }

  /**
   * Copy/move ctor, copy/move assignment operator, and destructor are
   * all explicitly defaulted for this class.
   */
  MortonSFCPartitioner (const MortonSFCPartitioner &) = default;
  MortonSFCPartitioner (MortonSFCPartitioner &&) = default;
  MortonSFCPartitioner & operator= (const MortonSFCPartitioner &) = default;
  MortonSFCPartitioner & operator= (MortonSFCPartitioner &&) = default;
  virtual ~MortonSFCPartitioner() = default;

  /**
   * \returns A copy of this partitioner wrapped in a smart pointer.
   */
  virtual std::unique_ptr<Partitioner> clone () const override
  {
    return libmesh_make_unique<MortonSFCPartitioner>(*this);
  }

protected:

  /**
   * Partition the \p MeshBase into \p n subdomains.
   */
  virtual void _do_partition (MeshBase & mesh,
                              const unsigned int n) override
  {
    SFCPartitioner::_do_partition (mesh, n);
  }
};

} // namespace libMesh

#endif // LIBMESH_MORTON_SFC_PARTITIONER_H
