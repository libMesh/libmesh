// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_SFC_PARTITIONER_H
#define LIBMESH_SFC_PARTITIONER_H

// Local Includes
#include "libmesh/partitioner.h"

// C++ Includes
#include <string>

namespace libMesh
{

/**
 * The \p SFCPartitioner uses a Hilbert or Morton-ordered space
 * filling curve to partition the elements.
 *
 * \author Benjamin S. Kirk
 * \date 2003
 * \brief Partitioner based on different types of space filling curves.
 */
class SFCPartitioner : public Partitioner
{
public:

  /**
   * Constructor.  Sets the default space filling
   * curve type to "Hilbert".
   */
  SFCPartitioner () :
    _sfc_type ("Hilbert")
  {}

  /**
   * Creates a new partitioner of this type and returns it in
   * a \p UniquePtr.
   */
  virtual UniquePtr<Partitioner> clone () const libmesh_override
  {
    return UniquePtr<Partitioner>(new SFCPartitioner());
  }

  /**
   * Sets the type of space-filling curve to use.  Valid types are
   * "Hilbert" (the default) and "Morton".
   */
  void set_sfc_type (const std::string & sfc_type)
  {
    libmesh_assert ((sfc_type == "Hilbert") ||
                    (sfc_type == "Morton"));

    _sfc_type = sfc_type;
  }

  /**
   * Called by the SubdomainPartitioner to partition elements in the range (it, end).
   */
  virtual void partition_range(MeshBase & mesh,
                               MeshBase::element_iterator it,
                               MeshBase::element_iterator end,
                               const unsigned int n) libmesh_override;

protected:

  /**
   * Partition the \p MeshBase into \p n subdomains.
   */
  virtual void _do_partition (MeshBase & mesh,
                              const unsigned int n) libmesh_override;


private:

  /**
   * The type of space-filling curve to use.  Hilbert by default.
   */
  std::string _sfc_type;
};

} // namespace libMesh

#endif // LIBMESH_SFC_PARTITIONER_H
