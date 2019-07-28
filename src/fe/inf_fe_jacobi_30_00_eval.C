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


// Local Includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

#include "libmesh/inf_fe.h"
#include "libmesh/jacobi_polynomials.h"

namespace libMesh
{

// Anonymous namespace for local helper functions
namespace {

Real jacobi_30_00_eval(unsigned n, Real x)
{
  if (n == 0)
    return 1.;

  Real val = JacobiPolynomials::value(n, /*alpha=*/3, /*beta=*/0, x);

  // For n>0, there is an even/odd shift of -1/+1 applied. I'm not
  // sure why this is done for the infinite elements, as it is not
  // part of the "standard" Jacobi polynomial definition, I'm just
  // copying what was done in the original implementation...
  return val + (n % 2 == 0 ? -1 : +1);
}

Real jacobi_30_00_eval_deriv(unsigned n, Real x)
{
  return JacobiPolynomials::deriv(n, /*alpha=*/3, /*beta=*/0, x);
}
} // anonymous namespace



// Specialize the eval() function for 1, 2, and 3 dimensions and the CARTESIAN mapping type
// to call the local helper function from the anonymous namespace.
template <> Real InfFE<1,JACOBI_30_00,CARTESIAN>::eval(Real x, Order, unsigned n) { return jacobi_30_00_eval(n, x); }
template <> Real InfFE<2,JACOBI_30_00,CARTESIAN>::eval(Real x, Order, unsigned n) { return jacobi_30_00_eval(n, x); }
template <> Real InfFE<3,JACOBI_30_00,CARTESIAN>::eval(Real x, Order, unsigned n) { return jacobi_30_00_eval(n, x); }

// Specialize the eval_deriv() function for 1, 2, and 3 dimensions and the CARTESIAN mapping type
// to call the local helper function from the anonymous namespace.
template <> Real InfFE<1,JACOBI_30_00,CARTESIAN>::eval_deriv(Real x, Order, unsigned n) { return jacobi_30_00_eval_deriv(n, x); }
template <> Real InfFE<2,JACOBI_30_00,CARTESIAN>::eval_deriv(Real x, Order, unsigned n) { return jacobi_30_00_eval_deriv(n, x); }
template <> Real InfFE<3,JACOBI_30_00,CARTESIAN>::eval_deriv(Real x, Order, unsigned n) { return jacobi_30_00_eval_deriv(n, x); }

} // namespace libMesh

#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
