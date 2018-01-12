// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

using namespace libMesh;

#ifndef LAPLACE_EXACT_SOLUTION_H
#define LAPLACE_EXACT_SOLUTION_H

class LaplaceExactSolution
{
public:
  LaplaceExactSolution(){}

  ~LaplaceExactSolution(){}

  Real operator() (unsigned int component,
                   Real x,
                   Real y,
                   Real z = 0.0)
  {
    const Real hp = 0.5*pi;

    switch(component)
      {
      case 0:
        return cos(hp*x)*sin(hp*y)*cos(hp*z);

      case 1:
        return sin(hp*x)*cos(hp*y)*cos(hp*z);

      case 2:
        return sin(hp*x)*cos(hp*y)*sin(hp*z);

      default:
        libmesh_error_msg("Invalid component = " << component);
      }
  }
};


class LaplaceExactGradient
{
public:
  LaplaceExactGradient(){}

  ~LaplaceExactGradient(){}

  RealGradient operator() (unsigned int component,
                           Real x,
                           Real y,
                           Real z = 0.0)
  {
    const Real hp = 0.5*pi;

    switch(component)
      {
      case 0:
        return RealGradient(-hp*sin(hp*x)*sin(hp*y)*cos(hp*z),
                            cos(hp*x)*(hp)*cos(hp*y)*cos(hp*z),
                            cos(hp*x)*sin(hp*y)*(-hp)*sin(hp*z));

      case 1:
        return RealGradient(hp*cos(hp*x)*cos(hp*y)*cos(hp*z),
                            sin(hp*x)*(-hp)*sin(hp*y)*cos(hp*z),
                            sin(hp*x)*cos(hp*y)*(-hp)*sin(hp*z));

      case 2:
        return RealGradient(hp*cos(hp*x)*cos(hp*y)*sin(hp*z),
                            sin(hp*x)*(-hp)*sin(hp*y)*sin(hp*z),
                            sin(hp*x)*cos(hp*y)*(hp)*cos(hp*z));

      default:
        libmesh_error_msg("Invalid component = " << component);
      }
  }
};

#endif // LAPLACE_EXACT_SOLUTION_H
