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

#include "libmesh/libmesh_common.h"
#include "libmesh/vector_value.h"

using namespace libMesh;

#ifndef CURL_CURL_EXACT_SOLUTION_H
#define CURL_CURL_EXACT_SOLUTION_H

class CurlCurlExactSolution
{
public:
  CurlCurlExactSolution(){}

  ~CurlCurlExactSolution(){}

  RealGradient operator() (Real x, Real y)
  {
    const Real ux =  cos(pi*x)*sin(pi*y);
    const Real uy = -sin(pi*x)*cos(pi*y);

    return RealGradient(ux, uy);
  }

  RealTensor grad(Real x, Real y)
  {
    const Real dux_dx = -pi*sin(pi*x)*sin(pi*y);
    const Real dux_dy = pi*cos(pi*x)*cos(pi*y);
    const Real duy_dx = -pi*cos(pi*x)*cos(pi*y);
    const Real duy_dy = pi*sin(pi*x)*sin(pi*y);

    return RealTensor(dux_dx, dux_dy, Real(0), duy_dx, duy_dy);
  }

  RealGradient curl(Real x, Real y)
  {
    const Real dux_dy =  pi*cos(pi*x)*cos(pi*y);
    const Real duy_dx = -pi*cos(pi*x)*cos(pi*y);

    return RealGradient(Real(0), Real(0), duy_dx - dux_dy);
  }

  RealGradient forcing(Real x, Real y)
  {
    const Real fx =  (2*pi*pi + 1)*cos(pi*x)*sin(pi*y);
    const Real fy = -(2*pi*pi + 1)*sin(pi*x)*cos(pi*y);

    return RealGradient(fx, fy);
  }
};

#endif // CURL_CURL_EXACT_SOLUTION_H
