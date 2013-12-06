// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"

namespace libMesh
{

  // ------------------------------------------------------------
  // Monomials-specific implementations


  // Anonymous namespace for local helper functions
  namespace {

    void monomial_nodal_soln(const Elem* elem,
			     const Order order,
			     const std::vector<Number>& elem_soln,
			     std::vector<Number>&       nodal_soln,
			     const unsigned Dim)
    {
      const unsigned int n_nodes = elem->n_nodes();

      const ElemType elem_type = elem->type();

      nodal_soln.resize(n_nodes);

      const Order totalorder = static_cast<Order>(order+elem->p_level());

      switch (totalorder)
	{
	  // Constant shape functions
	case CONSTANT:
	  {
	    libmesh_assert_equal_to (elem_soln.size(), 1);

	    const Number val = elem_soln[0];

	    for (unsigned int n=0; n<n_nodes; n++)
	      nodal_soln[n] = val;

	    return;
	  }


	  // For other orders, do interpolation at the nodes
	  // explicitly.
	default:
	  {
	    // FEType object to be passed to various FEInterface functions below.
	    FEType fe_type(totalorder, MONOMIAL);

	    const unsigned int n_sf =
	      // FE<Dim,T>::n_shape_functions(elem_type, totalorder);
	      FEInterface::n_shape_functions(Dim, fe_type, elem_type);

	    std::vector<Point> refspace_nodes;
	    FEBase::get_refspace_nodes(elem_type,refspace_nodes);
	    libmesh_assert_equal_to (refspace_nodes.size(), n_nodes);

	    for (unsigned int n=0; n<n_nodes; n++)
	      {
		libmesh_assert_equal_to (elem_soln.size(), n_sf);

		// Zero before summation
		nodal_soln[n] = 0;

		// u_i = Sum (alpha_i phi_i)
		for (unsigned int i=0; i<n_sf; i++)
		  nodal_soln[n] += elem_soln[i] *
		    // FE<Dim,T>::shape(elem, order, i, mapped_point);
		    FEInterface::shape(Dim, fe_type, elem, i, refspace_nodes[n]);
	      }

	    return;
	  } // default
	} // switch
    } // monomial_nodal_soln()




    unsigned int monomial_n_dofs(const ElemType t, const Order o)
    {
      switch (o)
	{

	  // constant shape functions
	  // no matter what shape there is only one DOF.
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

	      case EDGE2:
	      case EDGE3:
	      case EDGE4:
		return 2;

	      case TRI3:
	      case TRI6:
	      case QUAD4:
	      case QUAD8:
	      case QUAD9:
		return 3;

	      case TET4:
	      case TET10:
	      case HEX8:
	      case HEX20:
	      case HEX27:
	      case PRISM6:
	      case PRISM15:
	      case PRISM18:
	      case PYRAMID5:
              case PYRAMID14:
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

	      case EDGE2:
	      case EDGE3:
	      case EDGE4:
		return 3;

	      case TRI3:
	      case TRI6:
	      case QUAD4:
	      case QUAD8:
	      case QUAD9:
		return 6;

	      case TET4:
	      case TET10:
	      case HEX8:
	      case HEX20:
	      case HEX27:
	      case PRISM6:
	      case PRISM15:
	      case PRISM18:
	      case PYRAMID5:
	      case PYRAMID14:
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
	      case PYRAMID14:
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
	      case PYRAMID14:
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
	      case PYRAMID14:
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
	}

      libmesh_error();

      return 0;
    } // monomial_n_dofs()


  } // anonymous namespace





  // Do full-specialization for every dimension, instead
  // of explicit instantiation at the end of this file.
  // This could be macro-ified so that it fits on one line...
  template <>
  void FE<0,MONOMIAL>::nodal_soln(const Elem* elem,
				  const Order order,
				  const std::vector<Number>& elem_soln,
				  std::vector<Number>& nodal_soln)
  { monomial_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/0); }

  template <>
  void FE<1,MONOMIAL>::nodal_soln(const Elem* elem,
				  const Order order,
				  const std::vector<Number>& elem_soln,
				  std::vector<Number>& nodal_soln)
  { monomial_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/1); }

  template <>
  void FE<2,MONOMIAL>::nodal_soln(const Elem* elem,
				  const Order order,
				  const std::vector<Number>& elem_soln,
				  std::vector<Number>& nodal_soln)
  { monomial_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/2); }

  template <>
  void FE<3,MONOMIAL>::nodal_soln(const Elem* elem,
				  const Order order,
				  const std::vector<Number>& elem_soln,
				  std::vector<Number>& nodal_soln)
  { monomial_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/3); }


  // Full specialization of n_dofs() function for every dimension
  template <> unsigned int FE<0,MONOMIAL>::n_dofs(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
  template <> unsigned int FE<1,MONOMIAL>::n_dofs(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
  template <> unsigned int FE<2,MONOMIAL>::n_dofs(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
  template <> unsigned int FE<3,MONOMIAL>::n_dofs(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }

  // Full specialization of n_dofs_at_node() function for every dimension.
  // Monomials have no dofs at nodes, only element dofs.
  template <> unsigned int FE<0,MONOMIAL>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
  template <> unsigned int FE<1,MONOMIAL>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
  template <> unsigned int FE<2,MONOMIAL>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
  template <> unsigned int FE<3,MONOMIAL>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }

  // Full specialization of n_dofs_per_elem() function for every dimension.
  template <> unsigned int FE<0,MONOMIAL>::n_dofs_per_elem(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
  template <> unsigned int FE<1,MONOMIAL>::n_dofs_per_elem(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
  template <> unsigned int FE<2,MONOMIAL>::n_dofs_per_elem(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }
  template <> unsigned int FE<3,MONOMIAL>::n_dofs_per_elem(const ElemType t, const Order o) { return monomial_n_dofs(t, o); }


  // Full specialization of get_continuity() function for every dimension.
  template <> FEContinuity FE<0,MONOMIAL>::get_continuity() const { return DISCONTINUOUS; }
  template <> FEContinuity FE<1,MONOMIAL>::get_continuity() const { return DISCONTINUOUS; }
  template <> FEContinuity FE<2,MONOMIAL>::get_continuity() const { return DISCONTINUOUS; }
  template <> FEContinuity FE<3,MONOMIAL>::get_continuity() const { return DISCONTINUOUS; }

  // Full specialization of is_hierarchic() function for every dimension.
  // The monomials are hierarchic!
  template <> bool FE<0,MONOMIAL>::is_hierarchic() const { return true; }
  template <> bool FE<1,MONOMIAL>::is_hierarchic() const { return true; }
  template <> bool FE<2,MONOMIAL>::is_hierarchic() const { return true; }
  template <> bool FE<3,MONOMIAL>::is_hierarchic() const { return true; }

#ifdef LIBMESH_ENABLE_AMR

  // Full specialization of compute_constraints() function for 2D and
  // 3D only.  There are no constraints for discontinuous elements, so
  // we do nothing.
  template <> void FE<2,MONOMIAL>::compute_constraints (DofConstraints &, DofMap &, const unsigned int, const Elem*) {}
  template <> void FE<3,MONOMIAL>::compute_constraints (DofConstraints &, DofMap &, const unsigned int, const Elem*) {}

#endif // #ifdef LIBMESH_ENABLE_AMR

  // Full specialization of shapes_need_reinit() function for every dimension.
  template <> bool FE<0,MONOMIAL>::shapes_need_reinit() const { return false; }
  template <> bool FE<1,MONOMIAL>::shapes_need_reinit() const { return false; }
  template <> bool FE<2,MONOMIAL>::shapes_need_reinit() const { return false; }
  template <> bool FE<3,MONOMIAL>::shapes_need_reinit() const { return false; }

} // namespace libMesh
