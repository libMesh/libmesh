// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
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
#include "libmesh_config.h"
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
#include "inf_fe.h"

namespace libMesh
{



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
Real InfFE<Dim,T_radial,T_map>::eval(const Real v,
				     const Order, /* not used */
				     const unsigned int i)
{
  libmesh_assert (-1.-1.e-5 <= v && v < 1.);

  switch (i)
    {
    case 0:
      return -2.*v/(1.-v);
    case 1:
      return (1.+v)/(1.-v);

    default:
      libMesh::err << "bad index i = " << i << std::endl;
      libmesh_error();

    }

  // we never end up here.
  libmesh_error();
  return 0.;
}



template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
Real InfFE<Dim,T_radial,T_map>::eval_deriv(const Real v,
					   const Order, /* not used */
					   const unsigned int i)
{
  libmesh_assert (-1.-1.e-5 <= v && v < 1.);

  switch (i)
    {
    case 0:
      return -2./((1.-v)*(1.-v));
    case 1:
      return 2./((1.-v)*(1.-v));

    default:
      libMesh::err << "bad index i = " << i << std::endl;
      libmesh_error();

    }

  // we never end up here.
  libmesh_error();
  return 0.;
}



//--------------------------------------------------------------
// Explicit instantiations
// template class InfFE<1,INFINITE_MAP,CARTESIAN>;
// template class InfFE<2,INFINITE_MAP,CARTESIAN>;
// template class InfFE<3,INFINITE_MAP,CARTESIAN>;

// template class InfFE<1,INFINITE_MAP,SPHERICAL>;
// template class InfFE<2,INFINITE_MAP,SPHERICAL>;
// template class InfFE<3,INFINITE_MAP,SPHERICAL>;

// template class InfFE<1,INFINITE_MAP,ELLIPSOIDAL>;
// template class InfFE<2,INFINITE_MAP,ELLIPSOIDAL>;
// template class InfFE<3,INFINITE_MAP,ELLIPSOIDAL>;

template Real InfFE<1,INFINITE_MAP,CARTESIAN>::eval(const Real,const Order,const unsigned int);
template Real InfFE<1,INFINITE_MAP,CARTESIAN>::eval_deriv(const Real,const Order,const unsigned int);
template Real InfFE<2,INFINITE_MAP,CARTESIAN>::eval(const Real,const Order,const unsigned int);
template Real InfFE<2,INFINITE_MAP,CARTESIAN>::eval_deriv(const Real,const Order,const unsigned int);
template Real InfFE<3,INFINITE_MAP,CARTESIAN>::eval(const Real,const Order,const unsigned int);
template Real InfFE<3,INFINITE_MAP,CARTESIAN>::eval_deriv(const Real,const Order,const unsigned int);

} // namespace libMesh

#endif //ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

