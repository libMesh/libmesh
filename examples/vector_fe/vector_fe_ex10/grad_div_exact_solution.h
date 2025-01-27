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

#ifndef GRAD_DIV_EXACT_SOLUTION_H
#define GRAD_DIV_EXACT_SOLUTION_H

#include "libmesh/libmesh_common.h"
#include "libmesh/vector_value.h"

using namespace libMesh;

class GradDivExactSolution
{
public:
  GradDivExactSolution() = default;
  ~GradDivExactSolution() = default;

  RealGradient operator() (Real x, Real y, Real z)
  {
    libmesh_ignore(z);

    const Real ux = cos(k*x)*sin(k*y);
    const Real uy = sin(k*x)*cos(k*y);

    return RealGradient(ux, uy);
  }

  RealTensor grad(Real x, Real y, Real z)
  {
    libmesh_ignore(z);

    const Real dux_dx = -k*sin(k*x)*sin(k*y);
    const Real dux_dy =  k*cos(k*x)*cos(k*y);
    const Real duy_dx = dux_dy;
    const Real duy_dy = dux_dx;

    return RealTensor(dux_dx, dux_dy, Real(0), duy_dx, duy_dy);
  }

  RealGradient forcing(Real x, Real y, Real z)
  {
    return (2*k*k + 1)*operator()(x, y, z);
  }

private:
  const Real k = pi;
};

#endif // GRAD_DIV_EXACT_SOLUTION_H
