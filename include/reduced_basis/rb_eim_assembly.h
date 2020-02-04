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

template <typename> class ElemTempl;
typedef ElemTempl<Real> Elem;
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
  RBEIMAssembly(RBEIMConstruction & rb_eim_con_in,
                unsigned int basis_function_index_in);

  /**
   * Destructor.
   */
  virtual ~RBEIMAssembly();

  /**
   * Evaluate variable \p var_number of this object's EIM basis function
   * at the points \p points, where the points are in reference coordinates.
   * Fill \p values with the basis function values.
   */
  virtual void evaluate_basis_function(unsigned int var,
                                       const Elem & element,
                                       const std::vector<Point> & points,
                                       std::vector<Number> & values);

  /**
   * Get a reference to the RBEIMConstruction object.
   */
  RBEIMConstruction & get_rb_eim_construction();

  /**
   * Get a reference to the ghosted_basis_function.
   */
  NumericVector<Number> & get_ghosted_basis_function();

  /**
   * Retrieve the FE object.
   */
  FEBase & get_fe();

private:

  /**
   * Initialize the FE object.
   */
  void initialize_fe();

  /**
   * The RBEIMConstruction object that this RBEIMAssembly is based on.
   */
  RBEIMConstruction & _rb_eim_con;

  /**
   * The EIM basis function index (from rb_eim_eval) for this assembly object.
   */
  unsigned int _basis_function_index;

  /**
   * The basis function that we sample to evaluate the
   * empirical interpolation approximation. This will be a GHOSTED
   * vector to facilitate interpolation in the case of multiple processors.
   */
  std::unique_ptr<NumericVector<Number>> _ghosted_basis_function;

  /**
   * We store an FE object so we can easily reinit in evaluate_basis_function.
   */
  std::unique_ptr<FEBase> _fe;

};

}

#endif // LIBMESH_RB_EIM_ASSEMBLY_H
