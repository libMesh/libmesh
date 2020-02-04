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

#ifndef LIBMESH_PERIODIC_BOUNDARY_BASE_IMPL_H
#define LIBMESH_PERIODIC_BOUNDARY_BASE_IMPL_H

// Local Includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_PERIODIC

#include "libmesh/periodic_boundary_base.h"

#include "libmesh/auto_ptr.h" // libmesh_make_unique
#include "libmesh/boundary_info.h" // BoundaryInfo::invalid_id

namespace libMesh
{

template <typename RealType>
PeriodicBoundaryBaseTempl<RealType>::PeriodicBoundaryBaseTempl() :
  myboundary(BoundaryInfo::invalid_id),
  pairedboundary(BoundaryInfo::invalid_id)
{
}



template <typename RealType>
PeriodicBoundaryBaseTempl<RealType>::PeriodicBoundaryBaseTempl(const PeriodicBoundaryBase & o) :
  myboundary(o.myboundary),
  pairedboundary(o.pairedboundary),
  variables(o.variables)
{
  // Make a deep copy of _transformation_matrix, if it's not null
  if(o._transformation_matrix)
  {
    this->_transformation_matrix = libmesh_make_unique<DenseMatrix<Real>>();
    *(this->_transformation_matrix) = *(o._transformation_matrix);
  }
}



template <typename RealType>
void PeriodicBoundaryBaseTempl<RealType>::set_variable(unsigned int var)
{
  variables.insert(var);
}



template <typename RealType>
void PeriodicBoundaryBaseTempl<RealType>::merge(const PeriodicBoundaryBase & pb)
{
  variables.insert(pb.variables.begin(), pb.variables.end());
}



template <typename RealType>
bool PeriodicBoundaryBaseTempl<RealType>::is_my_variable(unsigned int var_num) const
{
  bool a = variables.empty() || (!variables.empty() && variables.find(var_num) != variables.end());
  return a;
}



template <typename RealType>
bool PeriodicBoundaryBaseTempl<RealType>::has_transformation_matrix() const
{
  return bool(_transformation_matrix);
}



template <typename RealType>
const DenseMatrix<Real> & PeriodicBoundaryBaseTempl<RealType>::get_transformation_matrix() const
{
  if(!has_transformation_matrix())
  {
    libmesh_error_msg("Transformation matrix is not defined");
  }

  return *_transformation_matrix;
}



template <typename RealType>
void PeriodicBoundaryBaseTempl<RealType>::set_transformation_matrix(const DenseMatrix<Real> & matrix)
{
  // Make a deep copy of matrix
  this->_transformation_matrix = libmesh_make_unique<DenseMatrix<Real>>();
  *(this->_transformation_matrix) = matrix;

  // if _transformation_matrix is defined then it must be the same sie as variables.
  libmesh_assert_equal_to(_transformation_matrix->m(), variables.size());
  libmesh_assert_equal_to(_transformation_matrix->n(), variables.size());
}



template <typename RealType>
const std::set<unsigned int> & PeriodicBoundaryBaseTempl<RealType>::get_variables() const
{
  return variables;
}

} // namespace libMesh

#endif // LIBMESH_ENABLE_PERIODIC

#endif // LIBMESH_PERIODIC_BOUNDARY_BASE_IMPL_H
