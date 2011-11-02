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
#include "elem.h"
#include "fe.h"
#include "fe_interface.h"

namespace libMesh
{

  // ------------------------------------------------------------
  // Hierarchic-specific implementations

  // Anonymous namespace for local helper functions
  namespace {

    void l2_hierarchic_nodal_soln(const Elem* elem,
				  const Order order,
				  const std::vector<Number>& elem_soln,
				  std::vector<Number>&       nodal_soln,
				  unsigned Dim)
    {
      const unsigned int n_nodes = elem->n_nodes();

      const ElemType elem_type = elem->type();

      nodal_soln.resize(n_nodes);

      const Order totalorder = static_cast<Order>(order + elem->p_level());

      // FEType object to be passed to various FEInterface functions below.
      FEType fe_type(totalorder, L2_HIERARCHIC);

      switch (totalorder)
	{
	  // Constant shape functions
	case CONSTANT:
	  {
	    libmesh_assert (elem_soln.size() == 1);

	    const Number val = elem_soln[0];

	    for (unsigned int n=0; n<n_nodes; n++)
	      nodal_soln[n] = val;

	    return;
	  }


	  // For other orders do interpolation at the nodes
	  // explicitly.
	default:
	  {

	    const unsigned int n_sf =
	      // FE<Dim,T>::n_shape_functions(elem_type, totalorder);
	      FEInterface::n_shape_functions(Dim, fe_type, elem_type);

	    for (unsigned int n=0; n<n_nodes; n++)
	      {
		const Point mapped_point =
		  // FE<Dim,T>::inverse_map(elem, elem->point(n));
		  FEInterface::inverse_map(Dim, fe_type, elem, elem->point(n));

		libmesh_assert (elem_soln.size() == n_sf);

		// Zero before summation
		nodal_soln[n] = 0;

		// u_i = Sum (alpha_i phi_i)
		for (unsigned int i=0; i<n_sf; i++)
		  nodal_soln[n] += elem_soln[i] *
		    // FE<Dim,T>::shape(elem, order, i, mapped_point);
		    FEInterface::shape(Dim, fe_type, elem, i, mapped_point);
	      }

	    return;
	  }
	}
    } // l2_hierarchic_nodal_soln()




    unsigned int l2_hierarchic_n_dofs(const ElemType t, const Order o)
    {
      libmesh_assert (o > 0);
      switch (t)
	{
	case NODEELEM:
	  return 1;
	case EDGE2:
	case EDGE3:
	  return (o+1);
	case QUAD4:
	  libmesh_assert(o < 2);
	case QUAD8:
	case QUAD9:
	  return ((o+1)*(o+1));
	case HEX8:
	  libmesh_assert(o < 2);
	case HEX20:
	  libmesh_assert(o < 2);
	case HEX27:
	  return ((o+1)*(o+1)*(o+1));
	case TRI6:
	  return ((o+1)*(o+2)/2);
	default:
	  libmesh_error();
	}

      libmesh_error();
      return 0;
    } // l2_hierarchic_n_dofs()



    unsigned int l2_hierarchic_n_dofs_per_elem(const ElemType t,
					       const Order o)
    {
      switch (o)
	{
	  // constant shape functions always have 1 DOF per element
	case CONSTANT:
	  return 1;


	  // Discontinuous linear shape functions
	  // expressed in the monomials.
	case FIRST:
	  {
	    switch (t)
	      {
	      case NODEELEM:
		return 1;

		// 1D linears have 2 DOFs per element
	      case EDGE2:
	      case EDGE3:
	      case EDGE4:
		return 2;

		// 2D linears have 3 DOFs per element
	      case TRI3:
	      case TRI6:
	      case QUAD4:
	      case QUAD8:
	      case QUAD9:
		return 3;

		// 3D linears have 4 DOFs per element
	      case TET4:
	      case TET10:
	      case HEX8:
	      case HEX20:
	      case HEX27:
	      case PRISM6:
	      case PRISM15:
	      case PRISM18:
	      case PYRAMID5:
		return 4;

	      default:
		{
#ifdef DEBUG
		  libMesh::err << "ERROR: Bad ElemType = " << t
			       << " for " << o << "th order approximation!"
			       << std::endl;
#endif
		  libmesh_error();
		}
	      }
	  }

	  // Discontinuous quadratic shape functions
	  // expressed in the monomials.
	case SECOND:
	  {
	    switch (t)
	      {
	      case NODEELEM:
		return 1;

		// 1D quadratics have 3 DOFs per element
	      case EDGE2:
	      case EDGE3:
	      case EDGE4:
		return 3;

		// 2D quadratics have 6 DOFs per element
	      case TRI3:
	      case TRI6:
	      case QUAD4:
	      case QUAD8:
	      case QUAD9:
		return 6;

		// 3D quadratics have 10 DOFs per element
	      case TET4:
	      case TET10:
	      case HEX8:
	      case HEX20:
	      case HEX27:
	      case PRISM6:
	      case PRISM15:
	      case PRISM18:
	      case PYRAMID5:
		return 10;

	      default:
		{
#ifdef DEBUG
		  libMesh::err << "ERROR: Bad ElemType = " << t
			       << " for " << o << "th order approximation!"
			       << std::endl;
#endif
		  libmesh_error();
		}
	      }
	  }


	  // Discontinuous cubic shape functions
	  // expressed in the monomials.
	case THIRD:
	  {
	    switch (t)
	      {
	      case NODEELEM:
		return 1;

	      case EDGE2:
	      case EDGE3:
	      case EDGE4:
		return 4;

	      case TRI3:
	      case TRI6:
	      case QUAD4:
	      case QUAD8:
	      case QUAD9:
		return 10;

	      case TET4:
	      case TET10:
	      case HEX8:
	      case HEX20:
	      case HEX27:
	      case PRISM6:
	      case PRISM15:
	      case PRISM18:
	      case PYRAMID5:
		return 20;

	      default:
		{
#ifdef DEBUG
		  libMesh::err << "ERROR: Bad ElemType = " << t
			       << " for " << o << "th order approximation!"
			       << std::endl;
#endif
		  libmesh_error();
		}
	      }
	  }


	  // Discontinuous quartic shape functions
	  // expressed in the monomials.
	case FOURTH:
	  {
	    switch (t)
	      {
	      case NODEELEM:
		return 1;

	      case EDGE2:
	      case EDGE3:
	      case EDGE4:
		return 5;

	      case TRI3:
	      case TRI6:
	      case QUAD4:
	      case QUAD8:
	      case QUAD9:
		return 15;

	      case TET4:
	      case TET10:
	      case HEX8:
	      case HEX20:
	      case HEX27:
	      case PRISM6:
	      case PRISM15:
	      case PRISM18:
	      case PYRAMID5:
		return 35;

	      default:
		{
#ifdef DEBUG
		  libMesh::err << "ERROR: Bad ElemType = " << t
			       << " for " << o << "th order approximation!"
			       << std::endl;
#endif
		  libmesh_error();
		}
	      }
	  }

	default:
	  {
	    const unsigned int order = static_cast<unsigned int>(o);
	    switch (t)
	      {
	      case NODEELEM:
		return 1;

	      case EDGE2:
	      case EDGE3:
		return (order+1);

	      case TRI3:
	      case TRI6:
	      case QUAD4:
	      case QUAD8:
	      case QUAD9:
		return (order+1)*(order+2)/2;

	      case TET4:
	      case TET10:
	      case HEX8:
	      case HEX20:
	      case HEX27:
	      case PRISM6:
	      case PRISM15:
	      case PRISM18:
	      case PYRAMID5:
		return (order+1)*(order+2)*(order+3)/6;

	      default:
		{
#ifdef DEBUG
		  libMesh::err << "ERROR: Bad ElemType = " << t
			       << " for " << o << "th order approximation!"
			       << std::endl;
#endif
		  libmesh_error();
		}
	      }
	  }

	  return 0;
	}
    } // l2_hierarchic_n_dofs_per_elem()

  } // anonymous namespace




  // Do full-specialization of nodal_soln() function for every
  // dimension, instead of explicit instantiation at the end of this
  // file.
  // This could be macro-ified so that it fits on one line...
  template <>
  void FE<0,L2_HIERARCHIC>::nodal_soln(const Elem* elem,
				       const Order order,
				       const std::vector<Number>& elem_soln,
				       std::vector<Number>& nodal_soln)
  { l2_hierarchic_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/0); }

  template <>
  void FE<1,L2_HIERARCHIC>::nodal_soln(const Elem* elem,
				       const Order order,
				       const std::vector<Number>& elem_soln,
				       std::vector<Number>& nodal_soln)
  { l2_hierarchic_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/1); }

  template <>
  void FE<2,L2_HIERARCHIC>::nodal_soln(const Elem* elem,
				       const Order order,
				       const std::vector<Number>& elem_soln,
				       std::vector<Number>& nodal_soln)
  { l2_hierarchic_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/2); }

  template <>
  void FE<3,L2_HIERARCHIC>::nodal_soln(const Elem* elem,
				       const Order order,
				       const std::vector<Number>& elem_soln,
				       std::vector<Number>& nodal_soln)
  { l2_hierarchic_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/3); }

  // Full specialization of n_dofs() function for every dimension
  template <> unsigned int FE<0,L2_HIERARCHIC>::n_dofs(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }
  template <> unsigned int FE<1,L2_HIERARCHIC>::n_dofs(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }
  template <> unsigned int FE<2,L2_HIERARCHIC>::n_dofs(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }
  template <> unsigned int FE<3,L2_HIERARCHIC>::n_dofs(const ElemType t, const Order o) { return l2_hierarchic_n_dofs(t, o); }

  // Full specialization of n_dofs_at_node() function for every dimension.
  // Discontinuous L2 elements only have interior nodes
  template <> unsigned int FE<0,L2_HIERARCHIC>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
  template <> unsigned int FE<1,L2_HIERARCHIC>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
  template <> unsigned int FE<2,L2_HIERARCHIC>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
  template <> unsigned int FE<3,L2_HIERARCHIC>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }

  // Full specialization of n_dofs_per_elem() function for every dimension.
  template <> unsigned int FE<0,L2_HIERARCHIC>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_hierarchic_n_dofs_per_elem(t, o); }
  template <> unsigned int FE<1,L2_HIERARCHIC>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_hierarchic_n_dofs_per_elem(t, o); }
  template <> unsigned int FE<2,L2_HIERARCHIC>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_hierarchic_n_dofs_per_elem(t, o); }
  template <> unsigned int FE<3,L2_HIERARCHIC>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_hierarchic_n_dofs_per_elem(t, o); }

  // L2 Hierarchic FEMs are C^0 continuous
  template <> FEContinuity FE<0,L2_HIERARCHIC>::get_continuity() const { return C_ZERO; }
  template <> FEContinuity FE<1,L2_HIERARCHIC>::get_continuity() const { return C_ZERO; }
  template <> FEContinuity FE<2,L2_HIERARCHIC>::get_continuity() const { return C_ZERO; }
  template <> FEContinuity FE<3,L2_HIERARCHIC>::get_continuity() const { return C_ZERO; }

  // L2 Hierarchic FEMs are hierarchic (duh!)
  template <> bool FE<0,L2_HIERARCHIC>::is_hierarchic() const { return true; }
  template <> bool FE<1,L2_HIERARCHIC>::is_hierarchic() const { return true; }
  template <> bool FE<2,L2_HIERARCHIC>::is_hierarchic() const { return true; }
  template <> bool FE<3,L2_HIERARCHIC>::is_hierarchic() const { return true; }

#ifdef LIBMESH_ENABLE_AMR
  // compute_constraints() specializations are only needed for 2 and 3D
  template <>
  void FE<2,L2_HIERARCHIC>::compute_constraints (DofConstraints &constraints,
						 DofMap &dof_map,
						 const unsigned int variable_number,
						 const Elem* elem)
  { compute_proj_constraints(constraints, dof_map, variable_number, elem); }

  template <>
  void FE<3,L2_HIERARCHIC>::compute_constraints (DofConstraints &constraints,
						 DofMap &dof_map,
						 const unsigned int variable_number,
						 const Elem* elem)
  { compute_proj_constraints(constraints, dof_map, variable_number, elem); }
#endif // #ifdef LIBMESH_ENABLE_AMR

  // L2-Hierarchic FEM shapes need reinit
  template <> bool FE<0,L2_HIERARCHIC>::shapes_need_reinit() const { return true; }
  template <> bool FE<1,L2_HIERARCHIC>::shapes_need_reinit() const { return true; }
  template <> bool FE<2,L2_HIERARCHIC>::shapes_need_reinit() const { return true; }
  template <> bool FE<3,L2_HIERARCHIC>::shapes_need_reinit() const { return true; }

} // namespace libMesh
