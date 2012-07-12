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

#ifndef __rb_eim_assembly_h__
#define __rb_eim_assembly_h__

// rbOOmit includes
#include "elem_assembly.h"

// libMesh includes
#include "auto_ptr.h"
#include "numeric_vector.h"
#include "point.h"

namespace libMesh
{

class Elem;
class RBParameters;
class RBEIMConstruction;

/**
 * This class provides functionality required to define an assembly
 * object that arises from an "Empirical Interpolation Method" (EIM)
 * approximation.
 *
 * @author David J. Knezevic, 2012
 */

class RBEIMAssembly : public ElemAssembly
{
public:

  /**
   * Constructor.
   */
  RBEIMAssembly(RBEIMConstruction& rb_eim_con_in,
                unsigned int basis_function_index_in);

  /**
   * Evaluate variable \p var_number of this object's EIM basis function
   * at the points \p qpoints. Fill \p values with the basis function values.
   */
  void evaluate_basis_function(unsigned int var,
                               const Elem& element,
                               const std::vector<Point>& qpoints,
                               std::vector<Number>& values);

private:

  /**
   * The RBEIMConstruction object that this RBEIMAssembly is based on.
   */
  RBEIMConstruction& _rb_eim_con;

  /**
   * The EIM basis function index (from rb_eim_eval) for this assembly object.
   */
  unsigned int _basis_function_index;

  /**
   * The basis function that we sample to evaluate the
   * empirical interpolation approximation. This will be a GHOSTED
   * vector to facilitate interpolation in the case of multiple processors.
   */
  AutoPtr< NumericVector<Number> > _ghosted_basis_function;

};

}

#endif
