// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local Includes -----------------------------------
#include "libmesh/partitioner.h"

// C++ Includes   -----------------------------------
#include <string>

namespace libMesh
{



/**
 * The \p SFCPartitioner uses a Hilbert or Morton-ordered space
 * filling curve to partition the elements.
 */

// ------------------------------------------------------------
// SFCPartitioner class definition
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
   * an \p AutoPtr.
   */
  virtual AutoPtr<Partitioner> clone () const {
    AutoPtr<Partitioner> cloned_partitioner
      (new SFCPartitioner());
    return cloned_partitioner;
  }

  /**
   * Sets the type of space-filling curve to use.  Valid types are
   * "Hilbert" (the default) and "Morton"
   */
  void set_sfc_type (const std::string& sfc_type);


protected:

  /**
   * Partition the \p MeshBase into \p n subdomains.
   */
  virtual void _do_partition (MeshBase& mesh,
			      const unsigned int n);


private:


  /**
   * The type of space-filling curve to use.  Hilbert by default.
   */
  std::string _sfc_type;

};



// ------------------------------------------------------------
// LinearPartitioner inline members
inline
void SFCPartitioner::set_sfc_type (const std::string& sfc_type)
{
  libmesh_assert ((sfc_type == "Hilbert") ||
	  (sfc_type == "Morton"));

  _sfc_type = sfc_type;
}


} // namespace libMesh



#endif // LIBMESH_SFC_PARTITIONER_H
