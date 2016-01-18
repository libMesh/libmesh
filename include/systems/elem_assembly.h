// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef LIBMESH_ELEM_ASSEMBLY_H
#define LIBMESH_ELEM_ASSEMBLY_H

#include "libmesh/reference_counted_object.h"

namespace libMesh
{

class FEMContext;
class System;
class Node;

/**
 * ElemAssembly provides a per-element (interior and boundary) assembly
 * functionality.
 *
 * \author David J. Knezevic
 * \date 2011
 */
class ElemAssembly : public ReferenceCountedObject<ElemAssembly>
{
public:

  /**
   * Constructor.  Initializes required
   * data structures.
   */
  ElemAssembly () {}

  /**
   * Destructor.
   */
  virtual ~ElemAssembly () {}

  /**
   * Perform the element interior assembly.
   */
  virtual void interior_assembly(FEMContext &) { }

  /**
   * Perform the element boundary assembly.
   */
  virtual void boundary_assembly(FEMContext &) { }

  /**
   * Get values to add to the RHS vector based on \p node.
   * This allows one to impose point loads, for example.
   */
  virtual void
  get_nodal_rhs_values(std::map<numeric_index_type, Number> & values,
                       const System &,
                       const Node &)
  {
    // By default, just clear the values map
    values.clear();
  }
};

}

#endif // LIBMESH_ELEM_ASSEMBLY_H
