// The libMesh Finite Element Library.
// Copyright (C) 2002-2021 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
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

#include "libmesh/fe_compute_data.h"
#include "libmesh/equation_systems.h"

namespace libMesh
{



void FEComputeData::clear ()
{
  this->shape.clear();
  this->dshape.clear();
  this->local_transform.clear();
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  this->phase = 0.;
  this->speed = 1.;
  this->frequency = 1.;

#endif
}



void FEComputeData::init ()
{
  if (!(this->shape.empty()))
    std::fill (this->shape.begin(),   this->shape.end(),   0.);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  if (equation_systems.parameters.have_parameter<Real>("speed"))
    this->speed = this->equation_systems.parameters.get<Real>("speed");

  // ensure that the wavenumber k=2. * libMesh::pi * frequency / speed is well-defined.
  libmesh_assert_not_equal_to(this->speed, 0);

  if (equation_systems.parameters.have_parameter<Number>("current frequency"))
    this->frequency = this->equation_systems.parameters.get<Number>("current frequency");

#if LIBMESH_USE_COMPLEX_NUMBERS
  else if (equation_systems.parameters.have_parameter<Real>("current frequency"))
    {
      // please use the type Number instead.
      libmesh_deprecated();
      this->frequency = static_cast<Number> (this->equation_systems.parameters.get<Real>("current frequency"));
    }
#endif

  this->phase = 0.;

#endif //LIBMESH_ENABLE_INFINITE_ELEMENTS

}


void FEComputeData::enable_derivative ()
{
  _need_dshape=true;
  if (!(this->dshape.empty()))
    std::fill (this->dshape.begin(),   this->dshape.end(),  0);
}


} // namespace libMesh
