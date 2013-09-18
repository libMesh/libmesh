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
#include "libmesh/tensor_value.h"

namespace libMesh
{

  // ------------------------------------------------------------
  // Nedelec first kind specific implementations


  // Anonymous namespace for local helper functions
  namespace {
    void nedelec_one_nodal_soln(const Elem* elem,
				const Order order,
				const std::vector<Number>& elem_soln,
				const int dim,
				std::vector<Number>&       nodal_soln)
    {
      const unsigned int n_nodes = elem->n_nodes();
      const ElemType elem_type   = elem->type();

      const Order totalorder = static_cast<Order>(order+elem->p_level());

      nodal_soln.resize(n_nodes*dim);

      FEType fe_type(totalorder, NEDELEC_ONE);

      switch (totalorder)
	{
	case FIRST:
	  {
	    switch (elem_type)
	      {
	      case TRI6:
		{
		  libmesh_assert_equal_to (elem_soln.size(), 3);
		  libmesh_assert_equal_to (nodal_soln.size(), 6*2);
		  break;
		}
	      case QUAD8:
	      case QUAD9:
		{
		  libmesh_assert_equal_to (elem_soln.size(), 4);

		  if (elem_type == QUAD8)
		    libmesh_assert_equal_to (nodal_soln.size(), 8*2);
		  else
		    libmesh_assert_equal_to (nodal_soln.size(), 9*2);
		  break;
		}
	      case TET10:
		{
		  libmesh_assert_equal_to (elem_soln.size(), 6);
		  libmesh_assert_equal_to (nodal_soln.size(), 10*3);

		  libmesh_not_implemented();

		  break;
		}


	      case HEX20:
	      case HEX27:
		{
		  libmesh_assert_equal_to (elem_soln.size(), 12);

		  if (elem_type == HEX20)
		    libmesh_assert_equal_to (nodal_soln.size(), 20*3);
		  else
		    libmesh_assert_equal_to (nodal_soln.size(), 27*3);

		  break;
		}

	      default:
		{
		  libmesh_error();

		  break;
		}

	      } // switch(elem_type)

	    const unsigned int n_sf =
	      FEInterface::n_shape_functions(dim, fe_type, elem_type);

	    std::vector<Point> refspace_nodes;
	    FEVectorBase::get_refspace_nodes(elem_type,refspace_nodes);
	    libmesh_assert_equal_to (refspace_nodes.size(), n_nodes);


	    // Need to create new fe object so the shape function as the FETransformation
	    // applied to it.
	    AutoPtr<FEVectorBase> vis_fe = FEVectorBase::build(dim,fe_type);

	    const std::vector<std::vector<RealGradient> >& vis_phi = vis_fe->get_phi();

	    vis_fe->reinit(elem,&refspace_nodes);

	    for( unsigned int n = 0; n < n_nodes; n++ )
	      {
		libmesh_assert_equal_to (elem_soln.size(), n_sf);

		// Zero before summation
                for( int d = 0; d < dim; d++ )
                  {
                    nodal_soln[dim*n+d] = 0;
                  }

		// u = Sum (u_i phi_i)
		for (unsigned int i=0; i<n_sf; i++)
		  {
                    for( int d = 0; d < dim; d++ )
                      {
                        nodal_soln[dim*n+d]   += elem_soln[i]*(vis_phi[i][n](d));
                      }
		  }
	      }

	    return;
	  } // case FIRST

	default:
	  {
	    libmesh_error();
	  }

	}//switch (totalorder)

      return;
    } // nedelec_one_nodal_soln


    unsigned int nedelec_one_n_dofs(const ElemType t, const Order o)
    {
      switch (o)
	{
	case FIRST:
	  {
	    switch (t)
	      {
	      case TRI6:
		return 3;

	      case QUAD8:
	      case QUAD9:
		return 4;

	      case TET10:
		return 6;

	      case HEX20:
	      case HEX27:
		return 12;

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

	default:
	  libmesh_error();
	}

      libmesh_error();
      return 0;
    }




    unsigned int nedelec_one_n_dofs_at_node(const ElemType t,
					    const Order o,
					    const unsigned int n)
    {
      switch (o)
	{
	case FIRST:
	  {
	    switch (t)
	      {
	      case TRI6:
		{
		  switch (n)
		    {
		    case 0:
		    case 1:
		    case 2:
		      return 0;
		    case 3:
		    case 4:
		    case 5:
		      return 1;

		    default:
		      libmesh_error();
		    }
		}
	      case QUAD8:
		{
		  switch (n)
		    {
		    case 0:
		    case 1:
		    case 2:
		    case 3:
		      return 0;
		    case 4:
		    case 5:
		    case 6:
		    case 7:
		      return 1;

		    default:
		      libmesh_error();
		    }
		}
	      case QUAD9:
		{
		  switch (n)
		    {
		    case 0:
		    case 1:
		    case 2:
		    case 3:
		    case 8:
		      return 0;
		    case 4:
		    case 5:
		    case 6:
		    case 7:
		      return 1;

		    default:
		      libmesh_error();
		    }
		}
	      case TET10:
		{
		  switch (n)
		    {
		    case 0:
		    case 1:
		    case 2:
		    case 3:
		      return 0;
		    case 4:
		    case 5:
		    case 6:
		    case 7:
		    case 8:
		    case 9:
		      return 1;

		    default:
		      libmesh_error();
		    }
		}

	      case HEX20:
		{
		  switch (n)
		    {
		    case 0:
		    case 1:
		    case 2:
		    case 3:
		    case 4:
		    case 5:
		    case 6:
		    case 7:
		      return 0;
		    case 8:
		    case 9:
		    case 10:
		    case 11:
		    case 12:
		    case 13:
		    case 14:
		    case 15:
		    case 16:
		    case 17:
		    case 18:
		    case 19:
		      return 1;

		    default:
		      libmesh_error();
		    }
		}
		case HEX27:
		  {
		  switch (n)
		    {
		    case 0:
		    case 1:
		    case 2:
		    case 3:
		    case 4:
		    case 5:
		    case 6:
		    case 7:
		    case 20:
		    case 21:
		    case 22:
		    case 23:
		    case 24:
		    case 25:
		    case 26:
		      return 0;
		    case 8:
		    case 9:
		    case 10:
		    case 11:
		    case 12:
		    case 13:
		    case 14:
		    case 15:
		    case 16:
		    case 17:
		    case 18:
		    case 19:
		      return 1;

		    default:
		      libmesh_error();
		    }
		  }
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

	default:
	  libmesh_error();
	}

      libmesh_error();
      return 0;

    }



#ifdef LIBMESH_ENABLE_AMR
    void nedelec_one_compute_constraints (DofConstraints &/*constraints*/,
					  DofMap &/*dof_map*/,
					  const unsigned int /*variable_number*/,
					  const Elem* elem,
					  const unsigned Dim)
    {
      // Only constrain elements in 2,3D.
      if (Dim == 1)
	return;

      libmesh_assert(elem);

      libmesh_not_implemented();

      /*
      // Only constrain active and ancestor elements
      if (elem->subactive())
	return;

      FEType fe_type = dof_map.variable_type(variable_number);
      fe_type.order = static_cast<Order>(fe_type.order + elem->p_level());

      std::vector<unsigned int> my_dof_indices, parent_dof_indices;

      // Look at the element faces.  Check to see if we need to
      // build constraints.
      for (unsigned int s=0; s<elem->n_sides(); s++)
	if (elem->neighbor(s) != NULL)
	  if (elem->neighbor(s)->level() < elem->level()) // constrain dofs shared between
	    {                                                     // this element and ones coarser
	      // than this element.
	      // Get pointers to the elements of interest and its parent.
	      const Elem* parent = elem->parent();

	      // This can't happen...  Only level-0 elements have NULL
	      // parents, and no level-0 elements can be at a higher
	      // level than their neighbors!
	      libmesh_assert(parent);

	      const AutoPtr<Elem> my_side     (elem->build_side(s));
	      const AutoPtr<Elem> parent_side (parent->build_side(s));

	      // This function gets called element-by-element, so there
	      // will be a lot of memory allocation going on.  We can
	      // at least minimize this for the case of the dof indices
	      // by efficiently preallocating the requisite storage.
	      my_dof_indices.reserve (my_side->n_nodes());
	      parent_dof_indices.reserve (parent_side->n_nodes());

	      dof_map.dof_indices (my_side.get(), my_dof_indices,
				   variable_number);
	      dof_map.dof_indices (parent_side.get(), parent_dof_indices,
				   variable_number);

	      for (unsigned int my_dof=0;
		   my_dof<FEInterface::n_dofs(Dim-1, fe_type, my_side->type());
		   my_dof++)
		{
		  libmesh_assert_less (my_dof, my_side->n_nodes());

		  // My global dof index.
		  const unsigned int my_dof_g = my_dof_indices[my_dof];

		  // The support point of the DOF
		  const Point& support_point = my_side->point(my_dof);

		  // Figure out where my node lies on their reference element.
		  const Point mapped_point = FEInterface::inverse_map(Dim-1, fe_type,
								      parent_side.get(),
								      support_point);

		  // Compute the parent's side shape function values.
		  for (unsigned int their_dof=0;
		       their_dof<FEInterface::n_dofs(Dim-1, fe_type, parent_side->type());
		       their_dof++)
		    {
		      libmesh_assert_less (their_dof, parent_side->n_nodes());

		      // Their global dof index.
		      const unsigned int their_dof_g =
			parent_dof_indices[their_dof];

		      const Real their_dof_value = FEInterface::shape(Dim-1,
								      fe_type,
								      parent_side->type(),
								      their_dof,
								      mapped_point);

		      // Only add non-zero and non-identity values
		      // for Lagrange basis functions.
		      if ((std::abs(their_dof_value) > 1.e-5) &&
			  (std::abs(their_dof_value) < .999))
			{
			  // since we may be running this method concurretly
			  // on multiple threads we need to acquire a lock
			  // before modifying the shared constraint_row object.
			  Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

			  // A reference to the constraint row.
			  DofConstraintRow& constraint_row = constraints[my_dof_g].first;

			  constraint_row.insert(std::make_pair (their_dof_g,
								their_dof_value));
			}
#ifdef DEBUG
		      // Protect for the case u_i = 0.999 u_j,
		      // in which case i better equal j.
		      else if (their_dof_value >= .999)
			libmesh_assert_equal_to (my_dof_g, their_dof_g);
#endif
		    }
		}
	    }
      */
    } // nedelec_one_compute_constrants()
#endif // #ifdef LIBMESH_ENABLE_AMR

  } // anonymous namespace

#define NEDELEC_LOW_D_ERROR_MESSAGE \
  libMesh::err << "ERROR: This method makes no sense for low-D elements!" \
	        << std::endl;                      \
  libmesh_error();


  // Do full-specialization for every dimension, instead
  // of explicit instantiation at the end of this file.
  template <>
  void FE<0,NEDELEC_ONE>::nodal_soln(const Elem*,
				     const Order,
				     const std::vector<Number>&,
				     std::vector<Number>&)
  { NEDELEC_LOW_D_ERROR_MESSAGE }

  template <>
  void FE<1,NEDELEC_ONE>::nodal_soln(const Elem*,
				     const Order,
				     const std::vector<Number>&,
				     std::vector<Number>&)
  { NEDELEC_LOW_D_ERROR_MESSAGE }

  template <>
  void FE<2,NEDELEC_ONE>::nodal_soln(const Elem* elem,
				     const Order order,
				     const std::vector<Number>& elem_soln,
				     std::vector<Number>& nodal_soln)
  { nedelec_one_nodal_soln(elem, order, elem_soln, 2 /*dim*/, nodal_soln); }

  template <>
  void FE<3,NEDELEC_ONE>::nodal_soln(const Elem* elem,
				     const Order order,
				     const std::vector<Number>& elem_soln,
				     std::vector<Number>& nodal_soln)
  { nedelec_one_nodal_soln(elem, order, elem_soln, 3 /*dim*/, nodal_soln); }


  // Do full-specialization for every dimension, instead
  // of explicit instantiation at the end of this function.
  // This could be macro-ified.
  template <> unsigned int FE<0,NEDELEC_ONE>::n_dofs(const ElemType, const Order) { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <> unsigned int FE<1,NEDELEC_ONE>::n_dofs(const ElemType, const Order) { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <> unsigned int FE<2,NEDELEC_ONE>::n_dofs(const ElemType t, const Order o) { return nedelec_one_n_dofs(t, o); }
  template <> unsigned int FE<3,NEDELEC_ONE>::n_dofs(const ElemType t, const Order o) { return nedelec_one_n_dofs(t, o); }


  // Do full-specialization for every dimension, instead
  // of explicit instantiation at the end of this function.
  template <> unsigned int FE<0,NEDELEC_ONE>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <> unsigned int FE<1,NEDELEC_ONE>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <> unsigned int FE<2,NEDELEC_ONE>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return nedelec_one_n_dofs_at_node(t, o, n); }
  template <> unsigned int FE<3,NEDELEC_ONE>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return nedelec_one_n_dofs_at_node(t, o, n); }


  // Nedelec first type elements have no dofs per element
  // FIXME: Only for first order!
  template <> unsigned int FE<0,NEDELEC_ONE>::n_dofs_per_elem(const ElemType, const Order) { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <> unsigned int FE<1,NEDELEC_ONE>::n_dofs_per_elem(const ElemType, const Order) { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <> unsigned int FE<2,NEDELEC_ONE>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
  template <> unsigned int FE<3,NEDELEC_ONE>::n_dofs_per_elem(const ElemType, const Order) { return 0; }

  // Nedelec first type FEMs are always tangentially continuous
  template <> FEContinuity FE<0,NEDELEC_ONE>::get_continuity() const { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <> FEContinuity FE<1,NEDELEC_ONE>::get_continuity() const { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <> FEContinuity FE<2,NEDELEC_ONE>::get_continuity() const { return H_CURL; }
  template <> FEContinuity FE<3,NEDELEC_ONE>::get_continuity() const { return H_CURL; }

  // Nedelec first type FEMs are not hierarchic
  template <> bool FE<0,NEDELEC_ONE>::is_hierarchic() const { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <> bool FE<1,NEDELEC_ONE>::is_hierarchic() const { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <> bool FE<2,NEDELEC_ONE>::is_hierarchic() const { return false; }
  template <> bool FE<3,NEDELEC_ONE>::is_hierarchic() const { return false; }

  // Nedelec first type FEM shapes always need to be reinit'ed (because of orientation dependence)
  template <> bool FE<0,NEDELEC_ONE>::shapes_need_reinit() const { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <> bool FE<1,NEDELEC_ONE>::shapes_need_reinit() const { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <> bool FE<2,NEDELEC_ONE>::shapes_need_reinit() const { return true; }
  template <> bool FE<3,NEDELEC_ONE>::shapes_need_reinit() const { return true; }

#ifdef LIBMESH_ENABLE_AMR
  template <>
  void FE<0,NEDELEC_ONE>::compute_constraints (DofConstraints &,
					       DofMap &,
					       const unsigned int,
					       const Elem*)
  { NEDELEC_LOW_D_ERROR_MESSAGE }

  template <>
  void FE<1,NEDELEC_ONE>::compute_constraints (DofConstraints &,
					       DofMap &,
					       const unsigned int,
					       const Elem*)
  { NEDELEC_LOW_D_ERROR_MESSAGE }

  template <>
  void FE<2,NEDELEC_ONE>::compute_constraints (DofConstraints &constraints,
					       DofMap &dof_map,
					       const unsigned int variable_number,
					       const Elem* elem)
  { nedelec_one_compute_constraints(constraints, dof_map, variable_number, elem, /*Dim=*/2); }

  template <>
  void FE<3,NEDELEC_ONE>::compute_constraints (DofConstraints &constraints,
					       DofMap &dof_map,
					       const unsigned int variable_number,
					       const Elem* elem)
  { nedelec_one_compute_constraints(constraints, dof_map, variable_number, elem, /*Dim=*/3); }
#endif // LIBMESH_ENABLE_AMR

  // Specialize useless shape function methods
  template <>
  RealGradient FE<0,NEDELEC_ONE>::shape(const ElemType, const Order,const unsigned int,const Point&)
  { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <>
  RealGradient FE<0,NEDELEC_ONE>::shape(const Elem*,const Order,const unsigned int,const Point&)
  { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <>
  RealGradient FE<0,NEDELEC_ONE>::shape_deriv(const ElemType, const Order,const unsigned int,
					      const unsigned int,const Point&)
  { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <>
  RealGradient FE<0,NEDELEC_ONE>::shape_deriv(const Elem*,const Order,const unsigned int,
					      const unsigned int,const Point&)
  { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <>
  RealGradient FE<0,NEDELEC_ONE>::shape_second_deriv(const ElemType, const Order,const unsigned int,
						     const unsigned int,const Point&)
  { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <>
  RealGradient FE<0,NEDELEC_ONE>::shape_second_deriv(const Elem*,const Order,const unsigned int,
						     const unsigned int,const Point&)
  { NEDELEC_LOW_D_ERROR_MESSAGE }

  template <>
  RealGradient FE<1,NEDELEC_ONE>::shape(const ElemType, const Order,const unsigned int,const Point&)
  { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <>
  RealGradient FE<1,NEDELEC_ONE>::shape(const Elem*,const Order,const unsigned int,const Point&)
  { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <>
  RealGradient FE<1,NEDELEC_ONE>::shape_deriv(const ElemType, const Order,const unsigned int,
					      const unsigned int,const Point&)
  { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <>
  RealGradient FE<1,NEDELEC_ONE>::shape_deriv(const Elem*,const Order,const unsigned int,
					      const unsigned int,const Point&)
  { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <>
  RealGradient FE<1,NEDELEC_ONE>::shape_second_deriv(const ElemType, const Order,const unsigned int,
						     const unsigned int,const Point&)
  { NEDELEC_LOW_D_ERROR_MESSAGE }
  template <>
  RealGradient FE<1,NEDELEC_ONE>::shape_second_deriv(const Elem*,const Order,const unsigned int,
						     const unsigned int,const Point&)
  { NEDELEC_LOW_D_ERROR_MESSAGE }

} // namespace libMesh
