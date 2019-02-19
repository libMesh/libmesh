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
#include "libmesh/inf_fe_macro.h"

using namespace libMesh;

// Anonymous namespace for local helper functions
namespace {

// Evaluate P_n(x) using the three term recurrence relation.  Note
// that using the three term recurrence is more accurate than the
// unrolled Horner scheme that was used here previously and requires a
// lot less code.
// TODO: C++17 will have std::legendre() in <cmath>, so eventually we
// can just use that.
Real legendre_eval(unsigned int n, Real x)
{
  Real p0 = 1;
  Real p1 = x;
  if (n == 0)
    return p0;

  unsigned int i = 1;
  while (i < n)
    {
      // Swapping saves a temporary, since this immediately updates p0
      // and then p1 is updated on the next line. Note that p0 and p1
      // appear in opposite positions than is usual in the update
      // formula because of this swap!
      std::swap(p0, p1);
      p1 = ((2*i + 1) * x * p0 - i * p1) / (i + 1);
      ++i;
    }

  // To match the original implementation, we add 1 to the result for
  // n odd, and subtract 1 for n even. I'm not sure why this is done
  // as it is not part of the "standard" Legendre polynomial definition.
  return p1 + (n % 2 == 0 ? -1 : +1);
}



// Evaluate d/dx P_n(x) using a recurrence relation.
// To evaluate the derivatives, we also need to evaluate the values,
// so this function duplicates some of legendre_eval(). That seems OK to me
// though since it is so simple. The recurrence relation we are using
// is from Wikipedia:
// (2n+1) * P_n = d/dx (P_{n+1} - P_{n-1})
// There are no plans for Legendre function derivatives to be in
// C++17, so we will need our own implementation for the forseeable
// future.  Finally, note that the shift which is applied to the
// values does not affect the derivatives, so it does not show up
// here.
Real legendre_eval_deriv(unsigned int n, Real x)
{
  if (n == 0)
    return 0;

  Real p0 = 1;
  Real p1 = x;

  // The recurrence relation we're using only requires one previous
  // derivative value, but it is from two iterations ago, so we still
  // need to keep track of two values.
  Real dp0 = 0;
  Real dp1 = 1;

  unsigned int i = 1;
  while (i < n)
    {
      // Swapping saves a temporary, since this immediately updates p0
      // and then p1 is updated on the next line. Note that p0 and p1
      // appear in opposite positions than is usual in the update
      // formula because of this swap!
      std::swap(p0, p1);
      p1 = ((2*i + 1) * x * p0 - i * p1) / (i + 1);

      // Note that dp1 appears in the formula below, but because we
      // swapped, it's really dp0. Also, we have already updated the
      // values, so:
      // p1 is now P_{i+1}(x)
      // p0 is now P_{i}(x)
      std::swap(dp0, dp1);
      dp1 = (2*i + 1) * p0 + dp1;
      ++i;
    }
  return dp1;
}

} // anonymous namespace


namespace libMesh
{

// Specialize the eval() function for 1, 2, and 3 dimensions and the CARTESIAN mapping type
// to call the local helper function from the anonymous namespace.
template <> Real InfFE<1,LEGENDRE,CARTESIAN>::eval(Real x, Order, unsigned n) { return legendre_eval(n, x); }
template <> Real InfFE<2,LEGENDRE,CARTESIAN>::eval(Real x, Order, unsigned n) { return legendre_eval(n, x); }
template <> Real InfFE<3,LEGENDRE,CARTESIAN>::eval(Real x, Order, unsigned n) { return legendre_eval(n, x); }

// Specialize the eval_deriv() function for 1, 2, and 3 dimensions and the CARTESIAN mapping type
// to call the local helper function from the anonymous namespace.
template <> Real InfFE<1,LEGENDRE,CARTESIAN>::eval_deriv(Real x, Order, unsigned n) { return legendre_eval_deriv(n, x); }
template <> Real InfFE<2,LEGENDRE,CARTESIAN>::eval_deriv(Real x, Order, unsigned n) { return legendre_eval_deriv(n, x); }
template <> Real InfFE<3,LEGENDRE,CARTESIAN>::eval_deriv(Real x, Order, unsigned n) { return legendre_eval_deriv(n, x); }

} // namespace libMesh


#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
