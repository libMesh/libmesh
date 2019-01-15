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



#ifndef LIBMESH_PARMETIS_PARTITIONER_H
#define LIBMESH_PARMETIS_PARTITIONER_H

// Local Includes
#include "libmesh/id_types.h"
#include "libmesh/partitioner.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// C++ Includes
#include <cstddef>
#include <unordered_map>
#include <vector>

namespace libMesh
{

// Forward declarations
class ParmetisHelper;

/**
 * The \p ParmetisPartitioner uses the Parmetis graph partitioner
 * to partition the elements.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 * \brief Partitioner which provides an interface to ParMETIS.
 */
class ParmetisPartitioner : public Partitioner
{
public:

  /**
   * Default and copy ctors.
   */
  ParmetisPartitioner ();
  ParmetisPartitioner (const ParmetisPartitioner & other);

  /**
   * This class contains a unique_ptr member, so it can't be default
   * copy assigned.
   */
  ParmetisPartitioner & operator= (const ParmetisPartitioner &) = delete;

  /**
   * Move ctor, move assignment operator, and destructor are
   * all explicitly inline-defaulted for this class.
   */
  ParmetisPartitioner (ParmetisPartitioner &&) = default;
  ParmetisPartitioner & operator= (ParmetisPartitioner &&) = default;

  /**
   * The destructor is out-of-line-defaulted to play nice with forward
   * declarations.
   */
  virtual ~ParmetisPartitioner();

  /**
   * \returns A copy of this partitioner wrapped in a smart pointer.
   */
  virtual std::unique_ptr<Partitioner> clone () const override
  {
    return libmesh_make_unique<ParmetisPartitioner>(*this);
  }


protected:

  /**
   * Parmetis can handle dynamically repartitioning a mesh such
   * that the redistribution costs are minimized.  This method
   * takes a previously partitioned mesh (which may have
   * then been adaptively refined) and repartitions it.
   */
  virtual void _do_repartition (MeshBase & mesh,
                                const unsigned int n) override;

  /**
   * Partition the \p MeshBase into \p n subdomains.
   */
  virtual void _do_partition (MeshBase & mesh,
                              const unsigned int n) override;

#ifdef LIBMESH_HAVE_PARMETIS
  /**
  * Build the graph.
  */
  virtual void build_graph (const MeshBase & mesh) override;

private:

  // These methods and data only need to be available if the
  // ParMETIS library is available.

  /**
   * Initialize data structures.
   */
  void initialize (const MeshBase & mesh, const unsigned int n_sbdmns);

  /**
   * Pointer to the Parmetis-specific data structures.  Lets us avoid
   * including parmetis.h here.
   */
  std::unique_ptr<ParmetisHelper> _pmetis;

#endif
};

} // namespace libMesh

#endif // LIBMESH_PARMETIS_PARTITIONER_H
