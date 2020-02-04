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



#ifndef DTKEVALUATOR_H
#define DTKEVALUATOR_H


#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_TRILINOS_HAVE_DTK

#include "libmesh/mesh.h"

#include "libmesh/ignore_warnings.h"
#include <DTK_MeshContainer.hpp>
#include <DTK_FieldEvaluator.hpp>
#include <DTK_FieldContainer.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include "libmesh/restore_warnings.h"

#include <string>

namespace libMesh
{

// Forward declarations
class EquationSystems;
class System;
class DofMap;
template <typename> class MeshBaseTempl;
typedef MeshBaseTempl<Real> MeshBase;
class FEType;

template <typename T> class NumericVector;

/**
 * Implements the evaluate() function to compute FE solution values at
 * points requested by DTK.
 *
 * \author Derek Gaston
 * \date 2013
 */
class DTKEvaluator : public DataTransferKit::FieldEvaluator<int,DataTransferKit::FieldContainer<double>>
{
public:
  typedef DataTransferKit::MeshContainer<int> MeshContainerType;
  typedef DataTransferKit::FieldContainer<Number> FieldContainerType;

  DTKEvaluator(System & in_sys, std::string var_name);

  virtual FieldContainerType evaluate(const Teuchos::ArrayRCP<int> & elements,
                                      const Teuchos::ArrayRCP<double> & coords) override;

protected:
  System & sys;
  NumericVector<Number> & current_local_solution;
  EquationSystems & es;
  MeshBase & mesh;
  unsigned int dim;
  DofMap & dof_map;
  unsigned int var_num;
  const FEType & fe_type;
};

} // namespace libMesh

#endif // #ifdef LIBMESH_TRILINOS_HAVE_DTK

#endif // #define DTKEVALUATOR_H
