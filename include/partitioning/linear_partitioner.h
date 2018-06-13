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



#ifndef LIBMESH_LINEAR_PARTITIONER_H
#define LIBMESH_LINEAR_PARTITIONER_H

// Local Includes
#include "libmesh/partitioner.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

namespace libMesh
{

/**
 * The \p LinearPartitioner simply takes the element list and splits
 * it into equal-sized chunks assigned to each processor.  Warning:
 * the resulting domain decomposition can be arbitrarily bad in terms
 * of edge-cut and other communication-based metrics!
 *
 * \author Benjamin S. Kirk
 * \date 2003
 * \brief Partitions the elements based solely on their ids.
 */
class LinearPartitioner : public Partitioner
{
public:

  /**
   * Ctors, assignment operators, and destructor are
   * all explicitly defaulted for this class.
   */
  LinearPartitioner () = default;
  LinearPartitioner (const LinearPartitioner &) = default;
  LinearPartitioner (LinearPartitioner &&) = default;
  LinearPartitioner & operator= (const LinearPartitioner &) = default;
  LinearPartitioner & operator= (LinearPartitioner &&) = default;
  virtual ~LinearPartitioner() = default;

  /**
   * \returns A copy of this partitioner wrapped in a smart pointer.
   */
  virtual std::unique_ptr<Partitioner> clone () const override
  {
    return libmesh_make_unique<LinearPartitioner>(*this);
  }

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

#endif  // LIBMESH_LINEAR_PARTITIONER_H
