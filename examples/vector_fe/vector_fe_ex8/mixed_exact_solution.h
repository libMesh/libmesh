// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#ifndef MIXED_EXACT_SOLUTION_H
#define MIXED_EXACT_SOLUTION_H

#include "libmesh/libmesh_common.h"
#include "libmesh/vector_value.h"

using namespace libMesh;

class MixedExactSolution
{
public:
  MixedExactSolution() = default;
  ~MixedExactSolution() = default;

  Real scalar(Real x, Real y)
  {
    return cos(.5*pi*x)*sin(.5*pi*y);
  }

  Real scalar(Real x, Real y, Real z)
  {
    return cos(.5*pi*x)*sin(.5*pi*y)*cos(.5*pi*z);
  }

  RealGradient operator() (Real x, Real y)
  {
    const Real ux =  .5*pi*sin(.5*pi*x)*sin(.5*pi*y);
    const Real uy = -.5*pi*cos(.5*pi*x)*cos(.5*pi*y);

    return RealGradient(ux, uy);
  }

  RealGradient operator() (Real x, Real y, Real z)
  {
    const Real ux =  .5*pi*sin(.5*pi*x)*sin(.5*pi*y)*cos(.5*pi*z);
    const Real uy = -.5*pi*cos(.5*pi*x)*cos(.5*pi*y)*cos(.5*pi*z);
    const Real uz =  .5*pi*cos(.5*pi*x)*sin(.5*pi*y)*sin(.5*pi*z);

    return RealGradient(ux, uy, uz);
  }

  RealTensor grad(Real x, Real y)
  {
    const Real dux_dx = .5*pi*.5*pi*cos(.5*pi*x)*sin(.5*pi*y);
    const Real dux_dy = .5*pi*.5*pi*sin(.5*pi*x)*cos(.5*pi*y);
    const Real duy_dx = .5*pi*.5*pi*sin(.5*pi*x)*cos(.5*pi*y);
    const Real duy_dy = .5*pi*.5*pi*cos(.5*pi*x)*sin(.5*pi*y);

    return RealTensor(dux_dx, dux_dy, Real(0), duy_dx, duy_dy);
  }

  RealTensor grad(Real x, Real y, Real z)
  {
    const Real dux_dx =  .5*pi*.5*pi*cos(.5*pi*x)*sin(.5*pi*y)*cos(.5*pi*z);
    const Real dux_dy =  .5*pi*.5*pi*sin(.5*pi*x)*cos(.5*pi*y)*cos(.5*pi*z);
    const Real dux_dz = -.5*pi*.5*pi*sin(.5*pi*x)*sin(.5*pi*y)*sin(.5*pi*z);
    const Real duy_dx =  .5*pi*.5*pi*sin(.5*pi*x)*cos(.5*pi*y)*cos(.5*pi*z);
    const Real duy_dy =  .5*pi*.5*pi*cos(.5*pi*x)*sin(.5*pi*y)*cos(.5*pi*z);
    const Real duy_dz =  .5*pi*.5*pi*cos(.5*pi*x)*cos(.5*pi*y)*sin(.5*pi*z);
    const Real duz_dx = -.5*pi*.5*pi*sin(.5*pi*x)*sin(.5*pi*y)*sin(.5*pi*z);
    const Real duz_dy =  .5*pi*.5*pi*cos(.5*pi*x)*cos(.5*pi*y)*sin(.5*pi*z);
    const Real duz_dz =  .5*pi*.5*pi*cos(.5*pi*x)*sin(.5*pi*y)*cos(.5*pi*z);

    return RealTensor(dux_dx, dux_dy, dux_dz, duy_dx, duy_dy, duy_dz, duz_dx, duz_dy, duz_dz);
  }

  Real div(Real x, Real y)
  {
    const Real dux_dx = grad(x, y)(0, 0);
    const Real duy_dy = grad(x, y)(1, 1);

    return dux_dx + duy_dy;
  }

  Real div(Real x, Real y, Real z)
  {
    const Real dux_dx = grad(x, y, z)(0, 0);
    const Real duy_dy = grad(x, y, z)(1, 1);
    const Real duz_dz = grad(x, y, z)(2, 2);

    return dux_dx + duy_dy + duz_dz;
  }

  Real forcing(Real x, Real y)
  {
    return div(x, y);
  }

  Real forcing(Real x, Real y, Real z)
  {
    return div(x, y, z);
  }
};

#endif // MIXED_EXACT_SOLUTION_H
