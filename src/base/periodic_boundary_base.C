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

// Local Includes -----------------------------------
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_PERIODIC

#include "libmesh/boundary_info.h" // BoundaryInfo::invalid_id
#include "libmesh/periodic_boundary_base.h"

namespace libMesh {

// ------------------------------------------------------------
// PeriodicBoundaryBase member functions


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


} // namespace libMesh

#endif // LIBMESH_ENABLE_PERIODIC
