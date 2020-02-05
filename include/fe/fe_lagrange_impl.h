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

#ifndef LIBMESH_FE_LAGRANGE_IMPL_H
#define LIBMESH_FE_LAGRANGE_IMPL_H

// Local includes
#include "libmesh/dof_map.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/elem.h"
#include "libmesh/remote_elem.h"
#include "libmesh/threads.h"
#include "libmesh/enum_to_string.h"

namespace libMesh
{
// Anonymous namespace for local helper functions
namespace {
template <typename RealType>
void lagrange_nodal_soln(const ElemTempl<RealType> * elem,
                         const Order order,
                         const std::vector<Number> & elem_soln,
                         std::vector<Number> &       nodal_soln)
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

          case PYRAMID13:
            {
              libmesh_assert_equal_to (elem_soln.size(), 5);
              libmesh_assert_equal_to (nodal_soln.size(), 13);

              nodal_soln[0]  = elem_soln[0];
              nodal_soln[1]  = elem_soln[1];
              nodal_soln[2]  = elem_soln[2];
              nodal_soln[3]  = elem_soln[3];
              nodal_soln[4]  = elem_soln[4];
              nodal_soln[5]  = .5*(elem_soln[0] + elem_soln[1]);
              nodal_soln[6]  = .5*(elem_soln[1] + elem_soln[2]);
              nodal_soln[7]  = .5*(elem_soln[2] + elem_soln[3]);
              nodal_soln[8]  = .5*(elem_soln[3] + elem_soln[0]);
              nodal_soln[9]  = .5*(elem_soln[0] + elem_soln[4]);
              nodal_soln[10] = .5*(elem_soln[1] + elem_soln[4]);
              nodal_soln[11] = .5*(elem_soln[2] + elem_soln[4]);
              nodal_soln[12] = .5*(elem_soln[3] + elem_soln[4]);

              return;
            }

          case PYRAMID14:
            {
              libmesh_assert_equal_to (elem_soln.size(), 5);
              libmesh_assert_equal_to (nodal_soln.size(), 14);

              nodal_soln[0]  = elem_soln[0];
              nodal_soln[1]  = elem_soln[1];
              nodal_soln[2]  = elem_soln[2];
              nodal_soln[3]  = elem_soln[3];
              nodal_soln[4]  = elem_soln[4];
              nodal_soln[5]  = .5*(elem_soln[0] + elem_soln[1]);
              nodal_soln[6]  = .5*(elem_soln[1] + elem_soln[2]);
              nodal_soln[7]  = .5*(elem_soln[2] + elem_soln[3]);
              nodal_soln[8]  = .5*(elem_soln[3] + elem_soln[0]);
              nodal_soln[9]  = .5*(elem_soln[0] + elem_soln[4]);
              nodal_soln[10] = .5*(elem_soln[1] + elem_soln[4]);
              nodal_soln[11] = .5*(elem_soln[2] + elem_soln[4]);
              nodal_soln[12] = .5*(elem_soln[3] + elem_soln[4]);
              nodal_soln[13] = .25*(elem_soln[0] + elem_soln[1] + elem_soln[2] + elem_soln[3]);

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



unsigned int lagrange_n_dofs(const ElemType t, const Order o)
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
          case TRISHELL3:
          case TRI3SUBDIVISION:
          case TRI6:
            return 3;

          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
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
          case PYRAMID13:
          case PYRAMID14:
            return 5;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
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
          case QUADSHELL8:
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

          case PYRAMID13:
            return 13;

          case PYRAMID14:
            return 14;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
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

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }

    default:
      libmesh_error_msg("ERROR: Invalid Order " << Utility::enum_to_string(o) << " selected for LAGRANGE FE family!");
    }
}




unsigned int lagrange_n_dofs_at_node(const ElemType t,
                                     const Order o,
                                     const unsigned int n)
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
            {
              switch (n)
                {
                case 0:
                case 1:
                  return 1;

                default:
                  return 0;
                }
            }

          case TRI3:
          case TRISHELL3:
          case TRI3SUBDIVISION:
          case TRI6:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                  return 1;

                default:
                  return 0;
                }
            }

          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                case 3:
                  return 1;

                default:
                  return 0;
                }
            }


          case TET4:
          case TET10:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                case 3:
                  return 1;

                default:
                  return 0;
                }
            }

          case HEX8:
          case HEX20:
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
                  return 1;

                default:
                  return 0;
                }
            }

          case PRISM6:
          case PRISM15:
          case PRISM18:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                  return 1;

                default:
                  return 0;
                }
            }

          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                case 3:
                case 4:
                  return 1;

                default:
                  return 0;
                }
            }

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }

      // quadratic Lagrange shape functions
    case SECOND:
      {
        switch (t)
          {
            // quadratic lagrange has one dof at each node
          case NODEELEM:
          case EDGE3:
          case TRI6:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
          case TET10:
          case HEX20:
          case HEX27:
          case PRISM15:
          case PRISM18:
          case PYRAMID13:
          case PYRAMID14:
            return 1;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }

    case THIRD:
      {
        switch (t)
          {
          case NODEELEM:
          case EDGE4:
            return 1;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }

    default:
      libmesh_error_msg("Unsupported order: " << o );
    }
}



#ifdef LIBMESH_ENABLE_AMR
template <typename RealType>
void lagrange_compute_constraints (DofConstraints & constraints,
                                   DofMap & dof_map,
                                   const unsigned int variable_number,
                                   const ElemTempl<RealType> * elem,
                                   const unsigned Dim)
{
  typedef PointTempl<RealType> Point;

  // Only constrain elements in 2,3D.
  if (Dim == 1)
    return;

  libmesh_assert(elem);

  // Only constrain active and ancestor elements
  if (elem->subactive())
    return;

  FEType fe_type = dof_map.variable_type(variable_number);
  fe_type.order = static_cast<Order>(fe_type.order + elem->p_level());

  // Pull objects out of the loop to reduce heap operations
  std::vector<dof_id_type> my_dof_indices, parent_dof_indices;
  std::unique_ptr<const Elem> my_side, parent_side;

  // Look at the element faces.  Check to see if we need to
  // build constraints.
  for (auto s : elem->side_index_range())
    if (elem->neighbor_ptr(s) != nullptr &&
        elem->neighbor_ptr(s) != RemoteElemTempl<RealType>::get_instance())
      if (elem->neighbor_ptr(s)->level() < elem->level()) // constrain dofs shared between
        {                                                 // this element and ones coarser
          // than this element.
          // Get pointers to the elements of interest and its parent.
          const Elem * parent = elem->parent();

          // This can't happen...  Only level-0 elements have nullptr
          // parents, and no level-0 elements can be at a higher
          // level than their neighbors!
          libmesh_assert(parent);

          elem->build_side_ptr(my_side, s);
          parent->build_side_ptr(parent_side, s);

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

          const unsigned int n_side_dofs =
            FEInterface::n_dofs(Dim-1, fe_type, my_side->type());
          const unsigned int n_parent_side_dofs =
            FEInterface::n_dofs(Dim-1, fe_type, parent_side->type());
          for (unsigned int my_dof=0; my_dof != n_side_dofs; my_dof++)
            {
              libmesh_assert_less (my_dof, my_side->n_nodes());

              // My global dof index.
              const dof_id_type my_dof_g = my_dof_indices[my_dof];

              // Hunt for "constraining against myself" cases before
              // we bother creating a constraint row
              bool self_constraint = false;
              for (unsigned int their_dof=0;
                   their_dof != n_parent_side_dofs; their_dof++)
                {
                  libmesh_assert_less (their_dof, parent_side->n_nodes());

                  // Their global dof index.
                  const dof_id_type their_dof_g =
                    parent_dof_indices[their_dof];

                  if (their_dof_g == my_dof_g)
                    {
                      self_constraint = true;
                      break;
                    }
                }

              if (self_constraint)
                continue;

              DofConstraintRow * constraint_row;

              // we may be running constraint methods concurrently
              // on multiple threads, so we need a lock to
              // ensure that this constraint is "ours"
              {
                Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

                if (dof_map.is_constrained_dof(my_dof_g))
                  continue;

                constraint_row = &(constraints[my_dof_g]);
                libmesh_assert(constraint_row->empty());
              }

              // The support point of the DOF
              const Point & support_point = my_side->point(my_dof);

              // Figure out where my node lies on their reference element.
              const Point mapped_point = FEMap::inverse_map(Dim-1,
                                                            parent_side.get(),
                                                            support_point);

              // Compute the parent's side shape function values.
              for (unsigned int their_dof=0;
                   their_dof != n_parent_side_dofs; their_dof++)
                {
                  libmesh_assert_less (their_dof, parent_side->n_nodes());

                  // Their global dof index.
                  const dof_id_type their_dof_g =
                    parent_dof_indices[their_dof];

                  const auto their_dof_value = FEInterface::shape(Dim-1,
                                                                  fe_type,
                                                                  parent_side->type(),
                                                                  their_dof,
                                                                  mapped_point);

                  // Only add non-zero and non-identity values
                  // for Lagrange basis functions.
                  if ((std::abs(their_dof_value) > 1.e-5) &&
                      (std::abs(their_dof_value) < .999))
                    {
                      constraint_row->insert(std::make_pair (their_dof_g,
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
} // lagrange_compute_constraints()
#endif // #ifdef LIBMESH_ENABLE_AMR

} // anonymous namespace


  // Do full-specialization for every dimension, instead
  // of explicit instantiation at the end of this file.
  // This could be macro-ified so that it fits on one line...
template <typename RealType>
void FEShim<0,LAGRANGE,RealType>::nodal_soln(const Elem * elem,
                                const Order order,
                                const std::vector<Number> & elem_soln,
                                std::vector<Number> & nodal_soln)
{ lagrange_nodal_soln(elem, order, elem_soln, nodal_soln); }

template <typename RealType>
void FEShim<1,LAGRANGE,RealType>::nodal_soln(const Elem * elem,
                                const Order order,
                                const std::vector<Number> & elem_soln,
                                std::vector<Number> & nodal_soln)
{ lagrange_nodal_soln(elem, order, elem_soln, nodal_soln); }

template <typename RealType>
void FEShim<2,LAGRANGE,RealType>::nodal_soln(const Elem * elem,
                                const Order order,
                                const std::vector<Number> & elem_soln,
                                std::vector<Number> & nodal_soln)
{ lagrange_nodal_soln(elem, order, elem_soln, nodal_soln); }

template <typename RealType>
void FEShim<3,LAGRANGE,RealType>::nodal_soln(const Elem * elem,
                                const Order order,
                                const std::vector<Number> & elem_soln,
                                std::vector<Number> & nodal_soln)
{ lagrange_nodal_soln(elem, order, elem_soln, nodal_soln); }


// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this function.
// This could be macro-ified.
template <typename RealType> unsigned int FEShim<0,LAGRANGE,RealType>::n_dofs(const ElemType t, const Order o) { return lagrange_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<1,LAGRANGE,RealType>::n_dofs(const ElemType t, const Order o) { return lagrange_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<2,LAGRANGE,RealType>::n_dofs(const ElemType t, const Order o) { return lagrange_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<3,LAGRANGE,RealType>::n_dofs(const ElemType t, const Order o) { return lagrange_n_dofs(t, o); }


// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this function.
template <typename RealType> unsigned int FEShim<0,LAGRANGE,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return lagrange_n_dofs_at_node(t, o, n); }
template <typename RealType> unsigned int FEShim<1,LAGRANGE,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return lagrange_n_dofs_at_node(t, o, n); }
template <typename RealType> unsigned int FEShim<2,LAGRANGE,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return lagrange_n_dofs_at_node(t, o, n); }
template <typename RealType> unsigned int FEShim<3,LAGRANGE,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return lagrange_n_dofs_at_node(t, o, n); }


// Lagrange elements have no dofs per element
// (just at the nodes)
template <typename RealType> unsigned int FEShim<0,LAGRANGE,RealType>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <typename RealType> unsigned int FEShim<1,LAGRANGE,RealType>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <typename RealType> unsigned int FEShim<2,LAGRANGE,RealType>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <typename RealType> unsigned int FEShim<3,LAGRANGE,RealType>::n_dofs_per_elem(const ElemType, const Order) { return 0; }

// Lagrange FEMs are always C^0 continuous
template <typename RealType> FEContinuity FEShim<0,LAGRANGE,RealType>::get_continuity() { return C_ZERO; }
template <typename RealType> FEContinuity FEShim<1,LAGRANGE,RealType>::get_continuity() { return C_ZERO; }
template <typename RealType> FEContinuity FEShim<2,LAGRANGE,RealType>::get_continuity() { return C_ZERO; }
template <typename RealType> FEContinuity FEShim<3,LAGRANGE,RealType>::get_continuity() { return C_ZERO; }

// Lagrange FEMs are not hierarchic
template <typename RealType> bool FEShim<0,LAGRANGE,RealType>::is_hierarchic() { return false; }
template <typename RealType> bool FEShim<1,LAGRANGE,RealType>::is_hierarchic() { return false; }
template <typename RealType> bool FEShim<2,LAGRANGE,RealType>::is_hierarchic() { return false; }
template <typename RealType> bool FEShim<3,LAGRANGE,RealType>::is_hierarchic() { return false; }

// Lagrange FEM shapes do not need reinit (is this always true?)
template <typename RealType> bool FEShim<0,LAGRANGE,RealType>::shapes_need_reinit() { return false; }
template <typename RealType> bool FEShim<1,LAGRANGE,RealType>::shapes_need_reinit() { return false; }
template <typename RealType> bool FEShim<2,LAGRANGE,RealType>::shapes_need_reinit() { return false; }
template <typename RealType> bool FEShim<3,LAGRANGE,RealType>::shapes_need_reinit() { return false; }

// Methods for computing Lagrange constraints.  Note: we pass the
// dimension as the last argument to the anonymous helper function.
// Also note: we only need instantiations of this function for
// Dim==2 and 3.
#ifdef LIBMESH_ENABLE_AMR
template <typename RealType>
void FEShim<2,LAGRANGE,RealType>::compute_constraints (DofConstraints & constraints,
                                          DofMap & dof_map,
                                          const unsigned int variable_number,
                                                       const ElemTempl<Real> * elem)
{ lagrange_compute_constraints(constraints, dof_map, variable_number, elem, /*Dim=*/2); }

template <typename RealType>
void FEShim<3,LAGRANGE,RealType>::compute_constraints (DofConstraints & constraints,
                                          DofMap & dof_map,
                                          const unsigned int variable_number,
                                                       const ElemTempl<Real> * elem)
{ lagrange_compute_constraints(constraints, dof_map, variable_number, elem, /*Dim=*/3); }
#endif // LIBMESH_ENABLE_AMR

} // namespace libMesh

#endif // LIBMESH_FE_LAGRANGE_IMPL_H
