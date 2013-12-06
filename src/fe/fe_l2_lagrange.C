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
#include "libmesh/dof_map.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/elem.h"
#include "libmesh/threads.h"

namespace libMesh
{

  // ------------------------------------------------------------
  // Lagrange-specific implementations


  // Anonymous namespace for local helper functions
  namespace {
    void l2_lagrange_nodal_soln(const Elem* elem,
                                const Order order,
			        const std::vector<Number>& elem_soln,
			        std::vector<Number>&       nodal_soln)
    {
      const unsigned int n_nodes = elem->n_nodes();
      const ElemType type        = elem->type();

      const Order totalorder = static_cast<Order>(order+elem->p_level());

      nodal_soln.resize(n_nodes);



      switch (totalorder)
	{
	  // linear Lagrange shape functions
	case FIRST:
	  {
	    switch (type)
	      {
	      case EDGE3:
		{
		  libmesh_assert_equal_to (elem_soln.size(), 2);
		  libmesh_assert_equal_to (nodal_soln.size(), 3);

		  nodal_soln[0] = elem_soln[0];
		  nodal_soln[1] = elem_soln[1];
		  nodal_soln[2] = .5*(elem_soln[0] + elem_soln[1]);

		  return;
		}

	      case EDGE4:
		{
		  libmesh_assert_equal_to (elem_soln.size(), 2);
		  libmesh_assert_equal_to (nodal_soln.size(), 4);

		  nodal_soln[0] = elem_soln[0];
		  nodal_soln[1] = elem_soln[1];
		  nodal_soln[2] = (2.*elem_soln[0] + elem_soln[1])/3.;
		  nodal_soln[3] = (elem_soln[0] + 2.*elem_soln[1])/3.;

		  return;
		}


	      case TRI6:
		{
		  libmesh_assert_equal_to (elem_soln.size(), 3);
		  libmesh_assert_equal_to (nodal_soln.size(), 6);

		  nodal_soln[0] = elem_soln[0];
		  nodal_soln[1] = elem_soln[1];
		  nodal_soln[2] = elem_soln[2];
		  nodal_soln[3] = .5*(elem_soln[0] + elem_soln[1]);
		  nodal_soln[4] = .5*(elem_soln[1] + elem_soln[2]);
		  nodal_soln[5] = .5*(elem_soln[2] + elem_soln[0]);

		  return;
		}


	      case QUAD8:
	      case QUAD9:
		{
		  libmesh_assert_equal_to (elem_soln.size(), 4);

		  if (type == QUAD8)
		    libmesh_assert_equal_to (nodal_soln.size(), 8);
		  else
		    libmesh_assert_equal_to (nodal_soln.size(), 9);


		  nodal_soln[0] = elem_soln[0];
		  nodal_soln[1] = elem_soln[1];
		  nodal_soln[2] = elem_soln[2];
		  nodal_soln[3] = elem_soln[3];
		  nodal_soln[4] = .5*(elem_soln[0] + elem_soln[1]);
		  nodal_soln[5] = .5*(elem_soln[1] + elem_soln[2]);
		  nodal_soln[6] = .5*(elem_soln[2] + elem_soln[3]);
		  nodal_soln[7] = .5*(elem_soln[3] + elem_soln[0]);

		  if (type == QUAD9)
		    nodal_soln[8] = .25*(elem_soln[0] + elem_soln[1] + elem_soln[2] + elem_soln[3]);

		  return;
		}


	      case TET10:
		{
		  libmesh_assert_equal_to (elem_soln.size(), 4);
		  libmesh_assert_equal_to (nodal_soln.size(), 10);

		  nodal_soln[0] = elem_soln[0];
		  nodal_soln[1] = elem_soln[1];
		  nodal_soln[2] = elem_soln[2];
		  nodal_soln[3] = elem_soln[3];
		  nodal_soln[4] = .5*(elem_soln[0] + elem_soln[1]);
		  nodal_soln[5] = .5*(elem_soln[1] + elem_soln[2]);
		  nodal_soln[6] = .5*(elem_soln[2] + elem_soln[0]);
		  nodal_soln[7] = .5*(elem_soln[3] + elem_soln[0]);
		  nodal_soln[8] = .5*(elem_soln[3] + elem_soln[1]);
		  nodal_soln[9] = .5*(elem_soln[3] + elem_soln[2]);

		  return;
		}


	      case HEX20:
	      case HEX27:
		{
		  libmesh_assert_equal_to (elem_soln.size(), 8);

		  if (type == HEX20)
		    libmesh_assert_equal_to (nodal_soln.size(), 20);
		  else
		    libmesh_assert_equal_to (nodal_soln.size(), 27);

		  nodal_soln[0]  = elem_soln[0];
		  nodal_soln[1]  = elem_soln[1];
		  nodal_soln[2]  = elem_soln[2];
		  nodal_soln[3]  = elem_soln[3];
		  nodal_soln[4]  = elem_soln[4];
		  nodal_soln[5]  = elem_soln[5];
		  nodal_soln[6]  = elem_soln[6];
		  nodal_soln[7]  = elem_soln[7];
		  nodal_soln[8]  = .5*(elem_soln[0] + elem_soln[1]);
		  nodal_soln[9]  = .5*(elem_soln[1] + elem_soln[2]);
		  nodal_soln[10] = .5*(elem_soln[2] + elem_soln[3]);
		  nodal_soln[11] = .5*(elem_soln[3] + elem_soln[0]);
		  nodal_soln[12] = .5*(elem_soln[0] + elem_soln[4]);
		  nodal_soln[13] = .5*(elem_soln[1] + elem_soln[5]);
		  nodal_soln[14] = .5*(elem_soln[2] + elem_soln[6]);
		  nodal_soln[15] = .5*(elem_soln[3] + elem_soln[7]);
		  nodal_soln[16] = .5*(elem_soln[4] + elem_soln[5]);
		  nodal_soln[17] = .5*(elem_soln[5] + elem_soln[6]);
		  nodal_soln[18] = .5*(elem_soln[6] + elem_soln[7]);
		  nodal_soln[19] = .5*(elem_soln[4] + elem_soln[7]);

		  if (type == HEX27)
		    {
		      nodal_soln[20] = .25*(elem_soln[0] + elem_soln[1] + elem_soln[2] + elem_soln[3]);
		      nodal_soln[21] = .25*(elem_soln[0] + elem_soln[1] + elem_soln[4] + elem_soln[5]);
		      nodal_soln[22] = .25*(elem_soln[1] + elem_soln[2] + elem_soln[5] + elem_soln[6]);
		      nodal_soln[23] = .25*(elem_soln[2] + elem_soln[3] + elem_soln[6] + elem_soln[7]);
		      nodal_soln[24] = .25*(elem_soln[3] + elem_soln[0] + elem_soln[7] + elem_soln[4]);
		      nodal_soln[25] = .25*(elem_soln[4] + elem_soln[5] + elem_soln[6] + elem_soln[7]);

		      nodal_soln[26] = .125*(elem_soln[0] + elem_soln[1] + elem_soln[2] + elem_soln[3] +
					     elem_soln[4] + elem_soln[5] + elem_soln[6] + elem_soln[7]);
		    }

		  return;
		}


	      case PRISM15:
	      case PRISM18:
		{
		  libmesh_assert_equal_to (elem_soln.size(), 6);

		  if (type == PRISM15)
		    libmesh_assert_equal_to (nodal_soln.size(), 15);
		  else
		    libmesh_assert_equal_to (nodal_soln.size(), 18);

		  nodal_soln[0]  = elem_soln[0];
		  nodal_soln[1]  = elem_soln[1];
		  nodal_soln[2]  = elem_soln[2];
		  nodal_soln[3]  = elem_soln[3];
		  nodal_soln[4]  = elem_soln[4];
		  nodal_soln[5]  = elem_soln[5];
		  nodal_soln[6]  = .5*(elem_soln[0] + elem_soln[1]);
		  nodal_soln[7]  = .5*(elem_soln[1] + elem_soln[2]);
		  nodal_soln[8]  = .5*(elem_soln[0] + elem_soln[2]);
		  nodal_soln[9]  = .5*(elem_soln[0] + elem_soln[3]);
		  nodal_soln[10] = .5*(elem_soln[1] + elem_soln[4]);
		  nodal_soln[11] = .5*(elem_soln[2] + elem_soln[5]);
		  nodal_soln[12] = .5*(elem_soln[3] + elem_soln[4]);
		  nodal_soln[13] = .5*(elem_soln[4] + elem_soln[5]);
		  nodal_soln[14] = .5*(elem_soln[3] + elem_soln[5]);

		  if (type == PRISM18)
		    {
		      nodal_soln[15] = .25*(elem_soln[0] + elem_soln[1] + elem_soln[4] + elem_soln[3]);
		      nodal_soln[16] = .25*(elem_soln[1] + elem_soln[2] + elem_soln[5] + elem_soln[4]);
		      nodal_soln[17] = .25*(elem_soln[2] + elem_soln[0] + elem_soln[3] + elem_soln[5]);
		    }

		  return;
		}



	      default:
		{
		  // By default the element solution _is_ nodal,
		  // so just copy it.
		  nodal_soln = elem_soln;

		  return;
		}
	      }
	  }

	case SECOND:
	  {
	    switch (type)
	      {
	      case EDGE4:
		{
		  libmesh_assert_equal_to (elem_soln.size(), 3);
		  libmesh_assert_equal_to (nodal_soln.size(), 4);

		  // Project quadratic solution onto cubic element nodes
		  nodal_soln[0] = elem_soln[0];
		  nodal_soln[1] = elem_soln[1];
		  nodal_soln[2] = (2.*elem_soln[0] - elem_soln[1] +
				   8.*elem_soln[2])/9.;
		  nodal_soln[3] = (-elem_soln[0] + 2.*elem_soln[1] +
				   8.*elem_soln[2])/9.;
		  return;
		}

	      default:
		{
		  // By default the element solution _is_ nodal,
		  // so just copy it.
		  nodal_soln = elem_soln;

		  return;
		}
	      }
	  }




	default:
	  {
	    // By default the element solution _is_ nodal,
	    // so just copy it.
	    nodal_soln = elem_soln;

	    return;
	  }
	}
    }


    // TODO: We should make this work, for example, for SECOND on a TRI3
    // (this is valid with L2_LAGRANGE, but not with LAGRANGE)
    unsigned int l2_lagrange_n_dofs(const ElemType t, const Order o)
    {
      switch (o)
	{

	  // linear Lagrange shape functions
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
		return 3;

	      case QUAD4:
	      case QUAD8:
	      case QUAD9:
		return 4;

	      case TET4:
	      case TET10:
		return 4;

	      case HEX8:
	      case HEX20:
	      case HEX27:
		return 8;

	      case PRISM6:
	      case PRISM15:
	      case PRISM18:
		return 6;

	      case PYRAMID5:
              case PYRAMID14:
		return 5;

	      default:
		{
#ifdef DEBUG
		  libMesh::err << "ERROR: Bad ElemType = " << t
			       << " for FIRST order approximation!"
			       << std::endl;
#endif
		  libmesh_error();
		}
	      }
	  }


	// quadratic Lagrange shape functions
	case SECOND:
	  {
	    switch (t)
	      {
	      case NODEELEM:
		return 1;

	      case EDGE3:
		return 3;

	      case TRI6:
		return 6;

	      case QUAD8:
		return 8;

	      case QUAD9:
		return 9;

	      case TET10:
		return 10;

	      case HEX20:
		return 20;

	      case HEX27:
		return 27;

	      case PRISM15:
		return 15;

	      case PRISM18:
		return 18;

	      default:
		{
#ifdef DEBUG
		  libMesh::err << "ERROR: Bad ElemType = " << t
			       << " for SECOND order approximation!"
			       << std::endl;
#endif
		  libmesh_error();
		}
	      }
	  }

	case THIRD:
	  {
	    switch (t)
	      {
	      case NODEELEM:
		return 1;

	      case EDGE4:
		return 4;

	      default:
		{
#ifdef DEBUG
		  libMesh::err << "ERROR: Bad ElemType = " << t
			       << " for THIRD order approximation!"
			       << std::endl;
#endif
		  libmesh_error();
		}
	      }
	  }

	default:
	  libmesh_error();
	}

      libmesh_error();
      return 0;
    }

  } // anonymous namespace


  // Do full-specialization for every dimension, instead
  // of explicit instantiation at the end of this file.
  // This could be macro-ified so that it fits on one line...
  template <>
  void FE<0,L2_LAGRANGE>::nodal_soln(const Elem* elem,
				  const Order order,
				  const std::vector<Number>& elem_soln,
				  std::vector<Number>& nodal_soln)
  { l2_lagrange_nodal_soln(elem, order, elem_soln, nodal_soln); }

  template <>
  void FE<1,L2_LAGRANGE>::nodal_soln(const Elem* elem,
				  const Order order,
				  const std::vector<Number>& elem_soln,
				  std::vector<Number>& nodal_soln)
  { l2_lagrange_nodal_soln(elem, order, elem_soln, nodal_soln); }

  template <>
  void FE<2,L2_LAGRANGE>::nodal_soln(const Elem* elem,
				  const Order order,
				  const std::vector<Number>& elem_soln,
				  std::vector<Number>& nodal_soln)
  { l2_lagrange_nodal_soln(elem, order, elem_soln, nodal_soln); }

  template <>
  void FE<3,L2_LAGRANGE>::nodal_soln(const Elem* elem,
				  const Order order,
				  const std::vector<Number>& elem_soln,
				  std::vector<Number>& nodal_soln)
  { l2_lagrange_nodal_soln(elem, order, elem_soln, nodal_soln); }


  // Do full-specialization for every dimension, instead
  // of explicit instantiation at the end of this function.
  // This could be macro-ified.
  template <> unsigned int FE<0,L2_LAGRANGE>::n_dofs(const ElemType t, const Order o) { return l2_lagrange_n_dofs(t, o); }
  template <> unsigned int FE<1,L2_LAGRANGE>::n_dofs(const ElemType t, const Order o) { return l2_lagrange_n_dofs(t, o); }
  template <> unsigned int FE<2,L2_LAGRANGE>::n_dofs(const ElemType t, const Order o) { return l2_lagrange_n_dofs(t, o); }
  template <> unsigned int FE<3,L2_LAGRANGE>::n_dofs(const ElemType t, const Order o) { return l2_lagrange_n_dofs(t, o); }


  // L2 Lagrange elements have all dofs on elements (hence none on nodes)
  template <> unsigned int FE<0,L2_LAGRANGE>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
  template <> unsigned int FE<1,L2_LAGRANGE>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
  template <> unsigned int FE<2,L2_LAGRANGE>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
  template <> unsigned int FE<3,L2_LAGRANGE>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }


  // L2 Lagrange elements have all dofs on elements
  template <> unsigned int FE<0,L2_LAGRANGE>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_lagrange_n_dofs(t, o); }
  template <> unsigned int FE<1,L2_LAGRANGE>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_lagrange_n_dofs(t, o); }
  template <> unsigned int FE<2,L2_LAGRANGE>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_lagrange_n_dofs(t, o); }
  template <> unsigned int FE<3,L2_LAGRANGE>::n_dofs_per_elem(const ElemType t, const Order o) { return l2_lagrange_n_dofs(t, o); }

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
					       const unsigned int ,
					       const Elem* )
  { }

  template <>
  void FE<3,L2_LAGRANGE>::compute_constraints (DofConstraints &,
					       DofMap &,
					       const unsigned int ,
					       const Elem* )
  { }
#endif // LIBMESH_ENABLE_AMR

} // namespace libMesh
