// The libMesh Finite Element Library.
// Copyright (C) 2002-2022 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_METIS_PARTITIONER_H
#define LIBMESH_METIS_PARTITIONER_H

// Local Includes
#include "libmesh/partitioner.h"

// C++ Includes
#include <memory>

namespace libMesh
{

/**
 * The \p MetisPartitioner uses the Metis graph partitioner
 * to partition the elements.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 * \brief Partitioner which interfaces with the METIS library.
 */
class MetisPartitioner : public Partitioner
{
public:

  /**
   * Ctors, assignment operators, and destructor are all explicitly
   * defaulted for this class.
   */
  MetisPartitioner () = default;
  MetisPartitioner (const MetisPartitioner &) = default;
  MetisPartitioner (MetisPartitioner &&) = default;
  MetisPartitioner & operator= (const MetisPartitioner &) = default;
  MetisPartitioner & operator= (MetisPartitioner &&) = default;
  virtual ~MetisPartitioner() = default;

  /**
   * \returns A copy of this partitioner wrapped in a smart pointer.
   */
  virtual std::unique_ptr<Partitioner> clone () const override
  {
    return std::make_unique<MetisPartitioner>(*this);
  }

  virtual void attach_weights(ErrorVector * weights) override { _weights = weights; }

  /**
   * Called by the SubdomainPartitioner to partition elements in the range (it, end).
   */
  virtual void partition_range(MeshBase & mesh,
                               MeshBase::element_iterator it,
                               MeshBase::element_iterator end,
                               const unsigned int n) override;

protected:

  /**
   * Partition the \p MeshBase into \p n subdomains.
   */
  virtual void _do_partition (MeshBase & mesh,
                              const unsigned int n) override;
};

} // namespace libMesh

#endif // LIBMESH_METIS_PARTITIONER_H
