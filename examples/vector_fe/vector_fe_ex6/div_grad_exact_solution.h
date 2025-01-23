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

#ifndef DIV_GRAD_EXACT_SOLUTION_H
#define DIV_GRAD_EXACT_SOLUTION_H

#include "libmesh/libmesh_common.h"
#include "libmesh/vector_value.h"

using namespace libMesh;

class DivGradExactSolution
{
public:
  DivGradExactSolution() = default;
  ~DivGradExactSolution() = default;

  Real scalar(Real x, Real y)
  {
    return cos(k*x)*sin(k*y);
  }

  Real scalar(Real x, Real y, Real z)
  {
    return cos(k*x)*sin(k*y)*cos(k*z);
  }

  RealGradient operator() (Real x, Real y)
  {
    const Real ux =  k*sin(k*x)*sin(k*y);
    const Real uy = -k*cos(k*x)*cos(k*y);

    return RealGradient(ux, uy);
  }

  RealGradient operator() (Real x, Real y, Real z)
  {
    const Real ux =  k*sin(k*x)*sin(k*y)*cos(k*z);
    const Real uy = -k*cos(k*x)*cos(k*y)*cos(k*z);
    const Real uz =  k*cos(k*x)*sin(k*y)*sin(k*z);

    return RealGradient(ux, uy, uz);
  }

  RealTensor grad(Real x, Real y)
  {
    const Real dux_dx = k*k*cos(k*x)*sin(k*y);
    const Real dux_dy = k*k*sin(k*x)*cos(k*y);
    const Real duy_dx = k*k*sin(k*x)*cos(k*y);
    const Real duy_dy = k*k*cos(k*x)*sin(k*y);

    return RealTensor(dux_dx, dux_dy, Real(0), duy_dx, duy_dy);
  }

  RealTensor grad(Real x, Real y, Real z)
  {
    const Real dux_dx =  k*k*cos(k*x)*sin(k*y)*cos(k*z);
    const Real dux_dy =  k*k*sin(k*x)*cos(k*y)*cos(k*z);
    const Real dux_dz = -k*k*sin(k*x)*sin(k*y)*sin(k*z);
    const Real duy_dx =  k*k*sin(k*x)*cos(k*y)*cos(k*z);
    const Real duy_dy =  k*k*cos(k*x)*sin(k*y)*cos(k*z);
    const Real duy_dz =  k*k*cos(k*x)*cos(k*y)*sin(k*z);
    const Real duz_dx = -k*k*sin(k*x)*sin(k*y)*sin(k*z);
    const Real duz_dy =  k*k*cos(k*x)*cos(k*y)*sin(k*z);
    const Real duz_dz =  k*k*cos(k*x)*sin(k*y)*cos(k*z);

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

private:
  const Real k = .5*pi;
};

#endif // DIV_GRAD_EXACT_SOLUTION_H
