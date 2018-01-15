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

#ifndef LIBMESH_ELEM_ASSEMBLY_H
#define LIBMESH_ELEM_ASSEMBLY_H

#include "libmesh/reference_counted_object.h"
#include "libmesh/dense_matrix.h"

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
  ElemAssembly () :
    is_nodal_rhs_values_overriden(true)
  {}

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
   * Get values to add to the matrix or rhs vector based on \p node.
   * This allows one to impose point loads or springs, for example.
   */
  virtual void
  get_nodal_values(std::vector<dof_id_type>& ,
                   DenseMatrix<Number>& ,
                   DenseVector<Number>& ,
                   const System & ,
                   const Node &)
  {
    // Do nothing by default
  }

  /**
   * Get values to add to the RHS vector based on \p node.
   * This allows one to impose point loads, for example.
   * Get values to add to the matrix or rhs vector based on \p node.
   * This allows one to impose point loads or springs, for example.
   */
#ifdef LIBMESH_ENABLE_DEPRECATED
  virtual void
  get_nodal_rhs_values(std::map<numeric_index_type, Number> & values,
                       const System &,
                       const Node &)
  {
    // use get_nodal_values instead, since that enables matrix and vector assembly
    libmesh_deprecated();

    // By default, just clear the values map
    values.clear();

    // Set flag so that we know this is the default implementation
    is_nodal_rhs_values_overriden = false;
  }
#endif

  /**
   * Temporary flag to help us figure out if we should call the
   * deprecated get_nodal_rhs_values method or not.
   */
  bool is_nodal_rhs_values_overriden;

protected:

};

}

#endif // LIBMESH_ELEM_ASSEMBLY_H
