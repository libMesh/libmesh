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

#include "fe_compute_data.h"
#include "equation_systems.h"

namespace libMesh
{



void FEComputeData::clear ()
{
  this->shape.clear();
#if defined(LIBMESH_ENABLE_INFINITE_ELEMENTS) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
  this->phase = 0.;
  this->speed = 0.;
#endif

#if defined (LIBMESH_ENABLE_INFINITE_ELEMENTS) && defined(LIBMESH_USE_COMPLEX_NUMBERS)
  this->speed = 0.;
  this->frequency = 0.;

#endif
}



void FEComputeData::init ()
{
  if (!(this->shape.empty()))
    std::fill (this->shape.begin(),   this->shape.end(),   0.);

#if defined(LIBMESH_ENABLE_INFINITE_ELEMENTS) && !defined(LIBMESH_USE_COMPLEX_NUMBERS)
  this->phase = 0.;

  if (equation_systems.parameters.have_parameter<Real>("speed"))
    this->speed = this->equation_systems.parameters.get<Real>("speed");
#endif

#if defined (LIBMESH_ENABLE_INFINITE_ELEMENTS) && defined(LIBMESH_USE_COMPLEX_NUMBERS)
  if (equation_systems.parameters.have_parameter<Real>("speed"))
    this->speed = this->equation_systems.parameters.get<Real>("speed");

  if (equation_systems.parameters.have_parameter<Real>("current frequency"))
    this->frequency = this->equation_systems.parameters.get<Real>("current frequency");

#endif
}


} // namespace libMesh
