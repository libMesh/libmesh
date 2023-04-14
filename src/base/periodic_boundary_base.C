// The libMesh Finite Element Library.
// Copyright (C) 2002-2023 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// Local Includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_PERIODIC

#include "libmesh/periodic_boundary_base.h"

#include "libmesh/boundary_info.h" // BoundaryInfo::invalid_id

// C++ Includes
#include <memory>


namespace libMesh
{

PeriodicBoundaryBase::PeriodicBoundaryBase() :
  myboundary(BoundaryInfo::invalid_id),
  pairedboundary(BoundaryInfo::invalid_id)
{
}



PeriodicBoundaryBase::PeriodicBoundaryBase(const PeriodicBoundaryBase & o) :
  myboundary(o.myboundary),
  pairedboundary(o.pairedboundary),
  variables(o.variables)
{
  // Make a deep copy of _transformation_matrix, if it's not null
  if(o._transformation_matrix)
  {
    this->_transformation_matrix = std::make_unique<DenseMatrix<Real>>();
    *(this->_transformation_matrix) = *(o._transformation_matrix);
  }
}



void PeriodicBoundaryBase::set_variable(unsigned int var)
{
  variables.insert(var);
}



void PeriodicBoundaryBase::merge(const PeriodicBoundaryBase & pb)
{
  variables.insert(pb.variables.begin(), pb.variables.end());
}



bool PeriodicBoundaryBase::is_my_variable(unsigned int var_num) const
{
  bool a = variables.empty() || (!variables.empty() && variables.find(var_num) != variables.end());
  return a;
}



bool PeriodicBoundaryBase::has_transformation_matrix() const
{
  return bool(_transformation_matrix);
}



const DenseMatrix<Real> & PeriodicBoundaryBase::get_transformation_matrix() const
{
  libmesh_error_msg_if(!has_transformation_matrix(),
                       "Transformation matrix is not defined");

  return *_transformation_matrix;
}



void PeriodicBoundaryBase::set_transformation_matrix(const DenseMatrix<Real> & matrix)
{
  // Make a deep copy of matrix
  this->_transformation_matrix = std::make_unique<DenseMatrix<Real>>();
  *(this->_transformation_matrix) = matrix;

  // if _transformation_matrix is defined then it must be the same sie as variables.
  libmesh_assert_equal_to(_transformation_matrix->m(), variables.size());
  libmesh_assert_equal_to(_transformation_matrix->n(), variables.size());
}



const std::set<unsigned int> & PeriodicBoundaryBase::get_variables() const
{
  return variables;
}


const EnforcementType & PeriodicBoundaryBase::get_enforcement_type() const
{
  return _enforcement_type;
}


void PeriodicBoundaryBase::set_enforcement_type(const EnforcementType & e_type)
{
  this->_enforcement_type = e_type;
}

} // namespace libMesh

#endif // LIBMESH_ENABLE_PERIODIC
