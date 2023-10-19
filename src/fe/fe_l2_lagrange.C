// The libMesh Finite Element Library.
// Copyright (C) 2002-2023 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_macro.h"

namespace libMesh
{

// ------------------------------------------------------------
// Lagrange-specific implementations


// Instantiate (side_) nodal_soln() function for every dimension
LIBMESH_FE_NODAL_SOLN(L2_LAGRANGE, lagrange_nodal_soln)
LIBMESH_FE_SIDE_NODAL_SOLN(L2_LAGRANGE)

// L2 Lagrange elements have all dofs on elements (hence none on nodes)
template <> unsigned int FE<0,L2_LAGRANGE>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<1,L2_LAGRANGE>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<2,L2_LAGRANGE>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<3,L2_LAGRANGE>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }

// L2 Lagrange FEMs are DISCONTINUOUS
template <> FEContinuity FE<0,L2_LAGRANGE>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<1,L2_LAGRANGE>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<2,L2_LAGRANGE>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<3,L2_LAGRANGE>::get_continuity() const { return DISCONTINUOUS; }

// Lagrange FEMs are not hierarchic
template <> bool FE<0,L2_LAGRANGE>::is_hierarchic() const { return false; }
template <> bool FE<1,L2_LAGRANGE>::is_hierarchic() const { return false; }
template <> bool FE<2,L2_LAGRANGE>::is_hierarchic() const { return false; }
template <> bool FE<3,L2_LAGRANGE>::is_hierarchic() const { return false; }

// Lagrange FEM shapes do not need reinit (is this always true?)
template <> bool FE<0,L2_LAGRANGE>::shapes_need_reinit() const { return false; }
template <> bool FE<1,L2_LAGRANGE>::shapes_need_reinit() const { return false; }
template <> bool FE<2,L2_LAGRANGE>::shapes_need_reinit() const { return false; }
template <> bool FE<3,L2_LAGRANGE>::shapes_need_reinit() const { return false; }

// We don't need any constraints for this DISCONTINUOUS basis, hence these
// are NOOPs
#ifdef LIBMESH_ENABLE_AMR
template <>
void FE<2,L2_LAGRANGE>::compute_constraints (DofConstraints &,
                                             DofMap &,
                                             const unsigned int,
                                             const Elem *)
{ }

template <>
void FE<3,L2_LAGRANGE>::compute_constraints (DofConstraints &,
                                             DofMap &,
                                             const unsigned int,
                                             const Elem *)
{ }
#endif // LIBMESH_ENABLE_AMR

} // namespace libMesh
