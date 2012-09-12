// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh_config.h"

#ifdef LIBMESH_ENABLE_PERIODIC

#include "libmesh.h" // libMesh::invalid_uint
#include "periodic_boundary.h"

// ------------------------------------------------------------
// PeriodicBoundary member functions


PeriodicBoundary::PeriodicBoundary() :
  myboundary(libMesh::invalid_uint),
  pairedboundary(libMesh::invalid_uint),
  translation_vector()
{
}




PeriodicBoundary::PeriodicBoundary(const PeriodicBoundary & o, bool inverse) :
  myboundary(o.myboundary),
  pairedboundary(o.pairedboundary),
  translation_vector(o.translation_vector),
  variables(o.variables)
{
  if (inverse)
    {
      std::swap(myboundary, pairedboundary);
      translation_vector *= -1.0;
    }
}



PeriodicBoundary::PeriodicBoundary(const RealVectorValue & vector) :
  myboundary(libMesh::invalid_uint),
  pairedboundary(libMesh::invalid_uint),
  translation_vector (vector)
{
}



Point PeriodicBoundary::get_corresponding_pos(const Point & pt) const
{
  return pt + translation_vector;
}



void PeriodicBoundary::set_variable(unsigned int var)
{
  variables.insert(var);
}



void PeriodicBoundary::merge(const PeriodicBoundary & pb)
{
  variables.insert(pb.variables.begin(), pb.variables.end());
}



bool PeriodicBoundary::is_my_variable(unsigned int var_num) const
{
  bool a = variables.empty() || (!variables.empty() && variables.find(var_num) != variables.end());
  return a;
}


#endif // LIBMESH_ENABLE_PERIODIC
