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



// Local includes
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
#include "libmesh/inf_fe.h"
#include "libmesh/inf_fe_macro.h"
#include "libmesh/fe.h"
#include "libmesh/elem.h"

namespace libMesh
{



// ------------------------------------------------------------
// InfFE::Base class members
template <unsigned int Dim, FEFamily T_radial, InfMapType T_base>
Elem * InfFE<Dim,T_radial,T_base>::Base::build_elem (const Elem * inf_elem)
{
  std::unique_ptr<const Elem> ape(inf_elem->build_side_ptr(0));

  // The incoming inf_elem is const, but this function is required to
  // return a non-const Elem * so that it can be used by
  // update_base_elem().  Therefore a const_cast seems to be
  // unavoidable here.
  return const_cast<Elem *>(ape.release());
}




template <unsigned int Dim, FEFamily T_radial, InfMapType T_base>
ElemType InfFE<Dim,T_radial,T_base>::Base::get_elem_type (const ElemType type)
{
  switch (type)
    {
      // 3D infinite elements:
      // with Dim=3 -> infinite elements on their own
    case INFHEX8:
      return QUAD4;

    case INFHEX16:
      return QUAD8;

    case INFHEX18:
      return QUAD9;

    case INFPRISM6:
      return TRI3;

    case INFPRISM12:
      return TRI6;

      // 2D infinite elements:
      // with Dim=3 -> used as boundary condition,
      // with Dim=2 -> infinite elements on their own
    case INFQUAD4:
      return EDGE2;

    case INFQUAD6:
      return EDGE3;

      // 1D infinite elements:
      // with Dim=2 -> used as boundary condition,
      // with Dim=1 -> infinite elements on their own,
      //               but no base element!
    case INFEDGE2:
      return INVALID_ELEM;

    default:
      libmesh_error_msg("ERROR: Unsupported element type!: " << type);
    }
}





template <unsigned int Dim, FEFamily T_radial, InfMapType T_base>
unsigned int InfFE<Dim,T_radial,T_base>::Base::n_base_mapping_sf (const ElemType base_elem_type,
                                                                  const Order base_mapping_order)
{
  if (Dim == 1)
    return 1;

  else if (Dim == 2)
    return FE<1,LAGRANGE>::n_shape_functions (base_elem_type,
                                              base_mapping_order);
  else if (Dim == 3)
    return FE<2,LAGRANGE>::n_shape_functions (base_elem_type,
                                              base_mapping_order);
  else
    {
      libmesh_error_msg("Unsupported Dim = " << Dim);
      return 0;
    }
}





// ------------------------------------------------------------
// InfFE::Radial class members
template <unsigned int Dim, FEFamily T_radial, InfMapType T_map>
unsigned int InfFE<Dim,T_radial,T_map>::Radial::n_dofs_at_node (const Order o_radial,
                                                                const unsigned int n_onion)
{
  libmesh_assert_less (n_onion, 2);

  if (n_onion == 0)
    /*
     * in the base, no matter what, we have 1 node associated
     * with radial direction
     */
    return 1;
  else
    /*
     * this works, since for Order o_radial=CONST=0, we still
     * have the (1-v)/2 mode, associated to the base
     */
    return static_cast<unsigned int>(o_radial);
}






//--------------------------------------------------------------
// Explicit instantiations
INSTANTIATE_INF_FE_MBRF(1,CARTESIAN,Elem *,Base::build_elem(const Elem *));
INSTANTIATE_INF_FE_MBRF(2,CARTESIAN,Elem *,Base::build_elem(const Elem *));
INSTANTIATE_INF_FE_MBRF(3,CARTESIAN,Elem *,Base::build_elem(const Elem *));
INSTANTIATE_INF_FE_MBRF(1,CARTESIAN,ElemType,Base::get_elem_type(const ElemType type));
INSTANTIATE_INF_FE_MBRF(2,CARTESIAN,ElemType,Base::get_elem_type(const ElemType type));
INSTANTIATE_INF_FE_MBRF(3,CARTESIAN,ElemType,Base::get_elem_type(const ElemType type));
INSTANTIATE_INF_FE_MBRF(1,CARTESIAN,unsigned int,Base::n_base_mapping_sf(const ElemType,const Order));
INSTANTIATE_INF_FE_MBRF(2,CARTESIAN,unsigned int,Base::n_base_mapping_sf(const ElemType,const Order));
INSTANTIATE_INF_FE_MBRF(3,CARTESIAN,unsigned int,Base::n_base_mapping_sf(const ElemType,const Order));
INSTANTIATE_INF_FE_MBRF(1,CARTESIAN,unsigned int,Radial::n_dofs_at_node (const Order,const unsigned int));
INSTANTIATE_INF_FE_MBRF(2,CARTESIAN,unsigned int,Radial::n_dofs_at_node (const Order,const unsigned int));
INSTANTIATE_INF_FE_MBRF(3,CARTESIAN,unsigned int,Radial::n_dofs_at_node (const Order,const unsigned int));

} // namespace libMesh

#endif //ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
