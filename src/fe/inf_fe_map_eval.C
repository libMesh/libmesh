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


// Local includes
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
#include "libmesh/inf_fe.h"

namespace libMesh
{

Real InfFEMap::eval(Real v, Order, unsigned i)
{
  libmesh_assert (-1.-1.e-5 <= v && v < 1.);

  switch (i)
    {
    case 0:
      return -2.*v/(1.-v);
    case 1:
      return (1.+v)/(1.-v);

    default:
      libmesh_error_msg("bad index i = " << i);
    }
}



Real InfFEMap::eval_deriv(Real v, Order, unsigned i)
{
  libmesh_assert (-1.-1.e-5 <= v && v < 1.);

  switch (i)
    {
    case 0:
      return -2./((1.-v)*(1.-v));
    case 1:
      return 2./((1.-v)*(1.-v));

    default:
      libmesh_error_msg("bad index i = " << i);
    }
}

// Specialize the eval() function for 1, 2, and 3 dimensions and the CARTESIAN mapping type
// to call the local helper function from the anonymous namespace.
template <> Real InfFE<1,INFINITE_MAP,CARTESIAN>::eval(Real v, Order o, unsigned i) { return InfFEMap::eval(v, o, i); }
template <> Real InfFE<2,INFINITE_MAP,CARTESIAN>::eval(Real v, Order o, unsigned i) { return InfFEMap::eval(v, o, i); }
template <> Real InfFE<3,INFINITE_MAP,CARTESIAN>::eval(Real v, Order o, unsigned i) { return InfFEMap::eval(v, o, i); }

// Specialize the eval_deriv() function for 1, 2, and 3 dimensions and the CARTESIAN mapping type
// to call the local helper function from the anonymous namespace.
template <> Real InfFE<1,INFINITE_MAP,CARTESIAN>::eval_deriv(Real v, Order o, unsigned i) { return InfFEMap::eval_deriv(v, o, i); }
template <> Real InfFE<2,INFINITE_MAP,CARTESIAN>::eval_deriv(Real v, Order o, unsigned i) { return InfFEMap::eval_deriv(v, o, i); }
template <> Real InfFE<3,INFINITE_MAP,CARTESIAN>::eval_deriv(Real v, Order o, unsigned i) { return InfFEMap::eval_deriv(v, o, i); }


} // namespace libMesh

#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
