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

#ifndef CURL_CURL_EXACT_SOLUTION_H
#define CURL_CURL_EXACT_SOLUTION_H

#include "libmesh/libmesh_common.h"
#include "libmesh/vector_value.h"

using namespace libMesh;

class CurlCurlExactSolution
{
public:
  CurlCurlExactSolution() = default;
  ~CurlCurlExactSolution() = default;

  RealGradient operator() (Point p)
  {
    Point pp = R.transpose()*p;
    Real x = pp(0), y = pp(1);

    const Real ux =  cos(k*x)*sin(k*y);
    const Real uy = -sin(k*x)*cos(k*y);

    return R*RealGradient(ux, uy);
  }

  RealTensor grad(Point p)
  {
    Point pp = R.transpose()*p;
    Real x = pp(0), y = pp(1);

    const Real dux_dx = -k*sin(k*x)*sin(k*y);
    const Real dux_dy =  k*cos(k*x)*cos(k*y);
    const Real duy_dx = -k*cos(k*x)*cos(k*y);
    const Real duy_dy =  k*sin(k*x)*sin(k*y);

    return R*RealTensor(dux_dx, dux_dy, Real(0), duy_dx, duy_dy)*R.transpose();
  }

  RealGradient forcing(Point p)
  {
    return (2*k*k + 1)*operator()(p);
  }

  static void RM(RealTensor T)
  {
    R = T;
  }

private:
  static RealTensor R;
  const Real k = pi;
};

#endif // CURL_CURL_EXACT_SOLUTION_H
