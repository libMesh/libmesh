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

#ifndef LIBMESH_RB_PARAMETRIZED_FUNCTION_H
#define LIBMESH_RB_PARAMETRIZED_FUNCTION_H

#include "libmesh/libmesh_common.h"


namespace libMesh
{

class RBParameters;
class Point;

/**
 * A simple functor class that provides a RBParameter-dependent function.
 *
 * \author David Knezevic
 * \date 2012
 * \brief Provides a reduced basis parameterized function.
 */
class RBParametrizedFunction
{
public:

  /**
   * Virtual evaluate() gives us a vtable, so there's no cost in adding a
   * virtual destructor for safety's sake.
   */
  virtual ~RBParametrizedFunction();

  /**
   * Specify the number of components in this parametrized function.
   * A scalar-valued function has one component, a vector-valued
   * function has more than one component.
   */
  virtual unsigned int get_n_components() const = 0;

  /**
   * Evaluate the parametrized function at the specified point for
   * parameter \p mu.
   */
  virtual Number evaluate(const RBParameters & mu,
                          unsigned int comp,
                          Point xyz,
                          subdomain_id_type subdomain_id) = 0;

  /**
   * Vectorized version of evaluate.
   */
  virtual void vectorized_evaluate(const RBParameters & mu,
                                   const std::unordered_map<dof_id_type, std::vector<Point>> & all_xyz,
                                   std::unordered_map<dof_id_type, subdomain_id_type> sbd_ids,
                                   std::unordered_map<dof_id_type, std::vector<std::vector<Number>>> & output);

  /**
   * Store the result of vectorized_evaluate. This is helpful during EIM training,
   * since we can pre-evaluate and store the parameterized function for each training
   * sample.
   */
  virtual void preevaluate_parametrized_function(const RBParameters & mu,
                                                 const std::unordered_map<dof_id_type, std::vector<Point>> & all_xyz,
                                                 std::unordered_map<dof_id_type, subdomain_id_type> sbd_ids);

  /**
   * Look up the preevaluate values of the parametrized function for
   * component \p comp, element \p elem_id, and quadrature point \p qp.
   */
  virtual Number lookup_preevaluated_value(unsigned int comp,
                                           dof_id_type elem_id,
                                           unsigned int qp) const;

  /**
   * Storage for pre-evaluated values. The indexing here is:
   *   elem_id --> comp --> qp --> value
   */
  std::unordered_map<dof_id_type, std::vector<std::vector<Number>>> preevaluated_values;
};

}

#endif // LIBMESH_RB_PARAMETRIZED_FUNCTION_H
