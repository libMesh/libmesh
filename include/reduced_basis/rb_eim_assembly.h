// rbOOmit: An implementation of the Certified Reduced Basis method.
// Copyright (C) 2009, 2010 David J. Knezevic

// This file is part of rbOOmit.

// rbOOmit is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// rbOOmit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#ifndef LIBMESH_RB_EIM_ASSEMBLY_H
#define LIBMESH_RB_EIM_ASSEMBLY_H

// rbOOmit includes
#include "libmesh/elem_assembly.h"

// libMesh includes
#include "libmesh/numeric_vector.h"
#include "libmesh/point.h"
#include "libmesh/fe.h"

// C++ includes
#include <memory>

namespace libMesh
{

class RBParameters;
class RBEIMConstruction;

/**
 * This class provides functionality required to define an assembly
 * object that arises from an "Empirical Interpolation Method" (EIM)
 * approximation.
 *
 * \author David J. Knezevic
 * \date 2012
 */
class RBEIMAssembly : public ElemAssembly
{
public:

  /**
   * Constructor.
   */
  RBEIMAssembly(RBEIMConstruction & rb_eim_eval_in,
                unsigned int basis_function_index_in);

  /**
   * Destructor.
   */
  virtual ~RBEIMAssembly();

  /**
   * Return the basis function values for all quadrature points for variable \p var
   * on element \p elem_id.
   */
  void evaluate_basis_function(dof_id_type elem_id,
                               unsigned int var,
                               std::vector<Number> & values);

  /**
   * Get a reference to the RBEIMEvaluation object.
   */
  RBEIMConstruction & get_rb_eim_construction();

private:

  /**
   * The RBEIMConstruction that the assembly data comes from.
   */
  RBEIMConstruction & _rb_eim_con;

  /**
   * The EIM basis function index (from _rb_eim_con's RBEIMEvaluation) for this assembly object.
   */
  unsigned int _basis_function_index;

};

}

#endif // LIBMESH_RB_EIM_ASSEMBLY_H
