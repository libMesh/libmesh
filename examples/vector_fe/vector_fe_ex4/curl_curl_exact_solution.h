/* The libMesh Finite Element Library. */
/* Copyright (C) 2003  Benjamin S. Kirk */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */

#include "libmesh/libmesh_common.h"
#include "libmesh/vector_value.h"

using namespace libMesh;

#ifndef __curl_curl_exact_solution_h__
#define __curl_curl_exact_solution_h__

class CurlCurlExactSolution
{
public:
  CurlCurlExactSolution(){}

  ~CurlCurlExactSolution(){}

  RealGradient operator()( Real x, Real y, Real z )
  {
    const Real ux = (1.0 - y*y)*(1.0 - z*z);
    const Real uy = (1.0 - x*x)*(1.0 - z*z);
    const Real uz = (1.0 - x*x)*(1.0 - y*y);

    return RealGradient( ux, uy, uz );
  }

  RealTensor grad( Real x, Real y, Real z )
  {
    const Real dux_dx = 0.0;
    const Real dux_dy = (1.0 - z*z)*(-2.0*y);
    const Real dux_dz = (1.0 - y*y)*(-2.0*z);

    const Real duy_dx = (1.0 - z*z)*(-2.0*x);
    const Real duy_dy = 0.0;
    const Real duy_dz = (1.0 - x*x)*(-2.0*z);

    const Real duz_dx = (1.0 - y*y)*(-2.0*x);
    const Real duz_dy = (1.0 - x*x)*(-2.0*y);
    const Real duz_dz = 0.0;

    return RealTensor( dux_dx, dux_dy, dux_dz, duy_dx, duy_dy, duy_dz, duz_dx, duz_dy, duz_dz );
  }

  RealGradient curl( Real x, Real y, Real z )
  {
    const Real duz_dy = (1.0 - x*x)*(-2.0*y);
    const Real duy_dz = (1.0 - x*x)*(-2.0*z);

    const Real dux_dz = (1.0 - y*y)*(-2.0*z);
    const Real duz_dx = (1.0 - y*y)*(-2.0*x);

    const Real dux_dy = (1.0 - z*z)*(-2.0*y);
    const Real duy_dx = (1.0 - z*z)*(-2.0*x);

    return RealGradient(duz_dy - duy_dz, dux_dz - duz_dx, duy_dx - dux_dy );
  }

  RealGradient forcing(  Real x, Real y, Real z )
  {
    const Real fx = 2.0*(1.0 - y*y) + 2.0*(1.0 - z*z) + (1.0 - y*y)*(1.0 - z*z);
    const Real fy = 2.0*(1.0 - x*x) + 2.0*(1.0 - z*z) + (1.0 - x*x)*(1.0 - z*z) ;
    const Real fz = 2.0*(1.0 - x*x) + 2.0*(1.0 - y*y) + (1.0 - x*x)*(1.0 - y*y);
    
    return RealGradient( fx, fy, fz );
  }

};

#endif // __curl_curl_exact_solution_h__
