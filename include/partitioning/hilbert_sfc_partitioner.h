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



#ifndef LIBMESH_HILBERT_SFC_PARTITIONER_H
#define LIBMESH_HILBERT_SFC_PARTITIONER_H

// Local Includes
#include "libmesh/sfc_partitioner.h"

namespace libMesh
{

/**
 * The \p HilbertSFCPartitioner uses a Hilbert space
 * filling curve to partition the elements.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 * \brief Partitioner based on Hilbert's space filling curve algorithm.
 */
template <typename RealType = Real>
class HilbertSFCPartitionerTempl : public SFCPartitionerTempl<RealType>
{
public:
  typedef HilbertSFCPartitionerTempl<RealType> HilbertSFCPartitioner;
  typedef SFCPartitionerTempl<RealType> SFCPartitioner;
  typedef PartitionerTempl<RealType> Partitioner;
  typedef MeshBaseTempl<RealType> MeshBase;
  typedef PointTempl<RealType> Point;
  typedef NodeTempl<RealType> Node;

  /**
   * Constructor. Sets the underlying space filling curve type.
   */
  HilbertSFCPartitionerTempl ()
  {
    this->set_sfc_type ("Hilbert");
  }

  /**
   * Copy/move ctor, copy/move assignment operator, and destructor are
   * all explicitly defaulted for this class.
   */
  HilbertSFCPartitionerTempl (const HilbertSFCPartitioner &) = default;
  HilbertSFCPartitionerTempl (HilbertSFCPartitioner &&) = default;
  HilbertSFCPartitioner & operator= (const HilbertSFCPartitioner &) = default;
  HilbertSFCPartitioner & operator= (HilbertSFCPartitioner &&) = default;
  virtual ~HilbertSFCPartitionerTempl() = default;

  /**
   * \returns A copy of this partitioner wrapped in a smart pointer.
   */
  virtual std::unique_ptr<Partitioner> clone () const override
  {
    return libmesh_make_unique<HilbertSFCPartitioner>(*this);
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

typedef HilbertSFCPartitionerTempl<Real> HilbertSFCPartitioner;

} // namespace libMesh

#endif // LIBMESH_HILBERT_SFC_PARTITIONER_H
