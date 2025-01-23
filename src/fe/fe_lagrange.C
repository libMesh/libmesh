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
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/fe.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_macro.h"
#include "libmesh/remote_elem.h"
#include "libmesh/threads.h"


#ifdef LIBMESH_ENABLE_AMR
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
// to have the macro 'inf_fe_switch' available.
#include "libmesh/fe_interface_macros.h"
#include "libmesh/inf_fe.h"
#endif
#endif

namespace libMesh
{

// global helper function
void lagrange_nodal_soln(const Elem * elem,
                         const Order order,
                         const std::vector<Number> & elem_soln,
                         std::vector<Number> &       nodal_soln,
                         const bool add_p_level)
{
  const unsigned int n_nodes = elem->n_nodes();
  const ElemType type        = elem->type();

  const Order totalorder = order + add_p_level*elem->p_level();

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


          case TRI7:
            libmesh_assert_equal_to (nodal_soln.size(), 7);
            nodal_soln[6] = (elem_soln[0] + elem_soln[1] + elem_soln[2])/3.;
            libmesh_fallthrough();
          case TRI6:
            {
              libmesh_assert (type == TRI7 || nodal_soln.size() == 6);
              libmesh_assert_equal_to (elem_soln.size(), 3);

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


          case TET14:
            libmesh_assert_equal_to (nodal_soln.size(), 14);
            nodal_soln[10] = (elem_soln[0] + elem_soln[1] + elem_soln[2])/3.;
            nodal_soln[11] = (elem_soln[0] + elem_soln[1] + elem_soln[3])/3.;
            nodal_soln[12] = (elem_soln[1] + elem_soln[2] + elem_soln[3])/3.;
            nodal_soln[13] = (elem_soln[0] + elem_soln[2] + elem_soln[3])/3.;
            libmesh_fallthrough();
          case TET10:
            {
              libmesh_assert_equal_to (elem_soln.size(), 4);
              libmesh_assert (type == TET14 || nodal_soln.size() == 10);

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


          case PRISM21:
            nodal_soln[20]  = (elem_soln[9] + elem_soln[10] + elem_soln[11])/Real(3);
            libmesh_fallthrough();
          case PRISM20:
            if (type == PRISM20)
              libmesh_assert_equal_to (nodal_soln.size(), 20);
            nodal_soln[18]  = (elem_soln[0] + elem_soln[1] + elem_soln[2])/Real(3);
            nodal_soln[19]  = (elem_soln[3] + elem_soln[4] + elem_soln[5])/Real(3);
            libmesh_fallthrough();
          case PRISM18:
            if (type == PRISM18)
              libmesh_assert_equal_to (nodal_soln.size(), 18);
            nodal_soln[15] = .25*(elem_soln[0] + elem_soln[1] + elem_soln[4] + elem_soln[3]);
            nodal_soln[16] = .25*(elem_soln[1] + elem_soln[2] + elem_soln[5] + elem_soln[4]);
            nodal_soln[17] = .25*(elem_soln[2] + elem_soln[0] + elem_soln[3] + elem_soln[5]);
            libmesh_fallthrough();
          case PRISM15:
            {
              libmesh_assert_equal_to (elem_soln.size(), 6);

              if (type == PRISM15)
                libmesh_assert_equal_to (nodal_soln.size(), 15);

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

              return;
            }

          case PYRAMID18:
            {
              libmesh_assert_equal_to (nodal_soln.size(), 18);

              nodal_soln[14] = (elem_soln[0] + elem_soln[1] + elem_soln[4])/Real(3);
              nodal_soln[15] = (elem_soln[1] + elem_soln[2] + elem_soln[4])/Real(3);
              nodal_soln[16] = (elem_soln[2] + elem_soln[3] + elem_soln[4])/Real(3);
              nodal_soln[17] = (elem_soln[0] + elem_soln[3] + elem_soln[4])/Real(3);

              libmesh_fallthrough();
            }

          case PYRAMID14:
            {
              if (type == PYRAMID14)
                libmesh_assert_equal_to (nodal_soln.size(), 14);

              nodal_soln[13] = .25*(elem_soln[0] + elem_soln[1] + elem_soln[2] + elem_soln[3]);

              libmesh_fallthrough();
            }

          case PYRAMID13:
            {
              libmesh_assert_equal_to (elem_soln.size(), 5);

              if (type == PYRAMID13)
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

          case TRI7:
            {
              libmesh_assert_equal_to (elem_soln.size(), 6);
              libmesh_assert_equal_to (nodal_soln.size(), 7);

              for (int i=0; i != 6; ++i)
                nodal_soln[i] = elem_soln[i];

              nodal_soln[6] = -1./9. * (elem_soln[0] + elem_soln[1] + elem_soln[2])
                              +4./9. * (elem_soln[3] + elem_soln[4] + elem_soln[5]);

              return;
            }

          case TET14:
            {
              libmesh_assert_equal_to (elem_soln.size(), 10);
              libmesh_assert_equal_to (nodal_soln.size(), 14);

              for (int i=0; i != 10; ++i)
                nodal_soln[i] = elem_soln[i];

              nodal_soln[10] = -1./9. * (elem_soln[0] + elem_soln[1] + elem_soln[2])
                               +4./9. * (elem_soln[4] + elem_soln[5] + elem_soln[6]);
              nodal_soln[11] = -1./9. * (elem_soln[0] + elem_soln[1] + elem_soln[3])
                               +4./9. * (elem_soln[4] + elem_soln[7] + elem_soln[8]);
              nodal_soln[12] = -1./9. * (elem_soln[1] + elem_soln[2] + elem_soln[3])
                               +4./9. * (elem_soln[5] + elem_soln[8] + elem_soln[9]);
              nodal_soln[13] = -1./9. * (elem_soln[0] + elem_soln[2] + elem_soln[3])
                               +4./9. * (elem_soln[6] + elem_soln[7] + elem_soln[9]);

              return;
            }

          case PRISM21:
            {
              nodal_soln[20]  = (elem_soln[9] + elem_soln[10] + elem_soln[11])/Real(3);
              libmesh_fallthrough();
            }
          case PRISM20:
            {
              if (type == PRISM20)
                libmesh_assert_equal_to (nodal_soln.size(), 20);

              for (int i=0; i != 18; ++i)
                nodal_soln[i] = elem_soln[i];

              nodal_soln[18]  = (elem_soln[0] + elem_soln[1] + elem_soln[2])/Real(3);
              nodal_soln[19]  = (elem_soln[3] + elem_soln[4] + elem_soln[5])/Real(3);
              return;
            }

          case PYRAMID18:
            {
              libmesh_assert_equal_to (nodal_soln.size(), 18);

              for (int i=0; i != 14; ++i)
                nodal_soln[i] = elem_soln[i];

              nodal_soln[14] = (elem_soln[0] + elem_soln[1] + elem_soln[4])/Real(3);
              nodal_soln[15] = (elem_soln[1] + elem_soln[2] + elem_soln[4])/Real(3);
              nodal_soln[16] = (elem_soln[2] + elem_soln[3] + elem_soln[4])/Real(3);
              nodal_soln[17] = (elem_soln[0] + elem_soln[3] + elem_soln[4])/Real(3);

              libmesh_fallthrough();
            }

          default:
            {
              // By default the element solution _is_ nodal, so just
              // copy the portion relevant to the nodal solution.
              // (this is the whole nodal solution for true Lagrange
              // elements, but a smaller part for "L2 Lagrange"
              libmesh_assert_less_equal(nodal_soln.size(), elem_soln.size());
              for (const auto i : index_range(nodal_soln))
                nodal_soln[i] = elem_soln[i];

              return;
            }
          }
      }




    default:
      {
        // By default the element solution _is_ nodal, so just copy
        // the portion relevant to the nodal solution.  (this is the
        // whole nodal solution for true Lagrange elements, but a
        // smaller part for "L2 Lagrange"
        libmesh_assert_less_equal(nodal_soln.size(), elem_soln.size());
        for (const auto i : index_range(nodal_soln))
          nodal_soln[i] = elem_soln[i];

        return;
      }
    }
}


// Anonymous namespace for local helper functions
namespace {
unsigned int lagrange_n_dofs(const ElemType t, const Order o)
{
  switch (o)
    {
      // lagrange can only be constant on a single node
    case CONSTANT:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }

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
          case TRI7:
            return 3;

          case QUAD4:
          case QUADSHELL4:
          case QUAD8:
          case QUADSHELL8:
          case QUAD9:
          case QUADSHELL9:
            return 4;

          case TET4:
          case TET10:
          case TET14:
            return 4;

          case HEX8:
          case HEX20:
          case HEX27:
            return 8;

          case PRISM6:
          case PRISM15:
          case PRISM18:
          case PRISM20:
          case PRISM21:
            return 6;

          case PYRAMID5:
          case PYRAMID13:
          case PYRAMID14:
          case PYRAMID18:
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
          case TRI7:
            return 6;

          case QUAD8:
          case QUADSHELL8:
            return 8;

          case QUAD9:
          case QUADSHELL9:
            return 9;

          case TET10:
          case TET14:
            return 10;

          case HEX20:
            return 20;

          case HEX27:
            return 27;

          case PRISM15:
            return 15;

          case PRISM18:
          case PRISM20:
          case PRISM21:
            return 18;

          case PYRAMID13:
            return 13;

          case PYRAMID14:
          case PYRAMID18:
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

          case PRISM20:
            return 20;

          case PRISM21:
            return 21;

          case PYRAMID18:
            return 18;

          case TRI7:
            return 7;

          case TET14:
            return 14;

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
      // lagrange can only be constant on a single node
    case CONSTANT:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          default:
            libmesh_error_msg("ERROR: Bad ElemType = " << Utility::enum_to_string(t) << " for " << Utility::enum_to_string(o) << " order approximation!");
          }
      }

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
          case TRI7:
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
          case QUADSHELL9:
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
          case TET14:
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
          case PRISM20:
          case PRISM21:
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
          case PYRAMID18:
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
          case QUADSHELL9:
          case TET10:
          case HEX20:
          case HEX27:
          case PRISM15:
          case PRISM18:
          case PYRAMID13:
          case PYRAMID14:
            return 1;

          case PRISM20:
          case PRISM21:
            return (n < 18);

          case PYRAMID18:
            return (n < 14);

          case TRI7:
            return (n < 6);

          case TET14:
            return (n < 10);

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
          case PRISM20:
          case PRISM21:
          case PYRAMID18:
          case TRI7:
          case TET14:
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
void lagrange_compute_constraints (DofConstraints & constraints,
                                   DofMap & dof_map,
                                   const unsigned int variable_number,
                                   const Elem * elem,
                                   const unsigned Dim)
{
  // Only constrain elements in 2,3D.
  if (Dim == 1)
    return;

  libmesh_assert(elem);

  // Only constrain active and ancestor elements
  if (elem->subactive())
    return;

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
   if (elem->infinite())
   {
      const FEType fe_t = dof_map.variable_type(variable_number);

      // expand the infinite_compute_constraint in its template-arguments.
      switch(Dim)
      {
         case 2:
            {
            inf_fe_family_mapping_switch(2, inf_compute_constraints (constraints, dof_map, variable_number, elem) , ,; break;);
            break;
          }
         case 3:
            {
            inf_fe_family_mapping_switch(3, inf_compute_constraints (constraints, dof_map, variable_number, elem) , ,; break;);
            break;
            }
         default:
           libmesh_error_msg("Invalid dim = " << Dim);
      }
      return;
   }
#endif

  const Variable & var = dof_map.variable(variable_number);
  const FEType fe_type = var.type();
  const bool add_p_level = dof_map.should_p_refine_var(variable_number);

  // Pull objects out of the loop to reduce heap operations
  std::vector<dof_id_type> my_dof_indices, parent_dof_indices;
  std::unique_ptr<const Elem> my_side, parent_side;

  // Look at the element faces.  Check to see if we need to
  // build constraints.
  for (auto s : elem->side_index_range())
    if (const Elem * const neigh = elem->neighbor_ptr(s);
        neigh && neigh != remote_elem)
      {
        // constrain dofs shared between this element and ones coarser
        if (neigh->level() < elem->level() &&
            var.active_on_subdomain(neigh->subdomain_id()))
          {
            // than this element.
            // Get pointers to the elements of interest and its parent.
            const Elem * parent = elem->parent();

            // This can't happen...  Only level-0 elements have nullptr
            // parents, and no level-0 elements can be at a higher
            // level than their neighbors!
            libmesh_assert(parent);

            elem->build_side_ptr(my_side, s);
            parent->build_side_ptr(parent_side, s);

            // We have some elements with interior bubble function bases
            // and lower-order sides.  We need to only ask for the lower
            // order in those cases.
            FEType side_fe_type = fe_type;
            int extra_order = 0;
            if (my_side->default_order() < fe_type.order)
              {
                side_fe_type.order = my_side->default_order();
                extra_order = (int)(side_fe_type.order) -
                              (int)(fe_type.order);
              }

            // This function gets called element-by-element, so there
            // will be a lot of memory allocation going on.  We can
            // at least minimize this for the case of the dof indices
            // by efficiently preallocating the requisite storage.
            my_dof_indices.reserve (my_side->n_nodes());
            parent_dof_indices.reserve (parent_side->n_nodes());

            dof_map.dof_indices (my_side.get(), my_dof_indices,
                                 variable_number, extra_order);
            dof_map.dof_indices (parent_side.get(), parent_dof_indices,
                                 variable_number, extra_order);

            const unsigned int n_side_dofs =
              FEInterface::n_dofs(side_fe_type, my_side.get());
            const unsigned int n_parent_side_dofs =
              FEInterface::n_dofs(side_fe_type, parent_side.get());
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

                    const Real their_dof_value = FEInterface::shape(Dim-1,
                                                                    side_fe_type,
                                                                    parent_side.get(),
                                                                    their_dof,
                                                                    mapped_point);

                    // Only add non-zero and non-identity values
                    // for Lagrange basis functions.
                    if ((std::abs(their_dof_value) > 1.e-5) &&
                        (std::abs(their_dof_value) < .999))
                      {
                        constraint_row->emplace(their_dof_g, their_dof_value);
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

        if (elem->active() && add_p_level)
          {
            // p refinement constraints:
            // constrain dofs shared between
            // active elements and neighbors with
            // lower polynomial degrees
            const unsigned int min_p_level =
              neigh->min_p_level_by_neighbor(elem, elem->p_level());
            if (min_p_level < elem->p_level())
              libmesh_error_msg(
                "Mismatched p-refinement levels for LAGRANGE is not currently supported");
          }
      }
} // lagrange_compute_constraints()
#endif // #ifdef LIBMESH_ENABLE_AMR

} // anonymous namespace


// Instantiate (side_) nodal_soln() function for every dimension
LIBMESH_FE_NODAL_SOLN(LAGRANGE, lagrange_nodal_soln)
LIBMESH_FE_SIDE_NODAL_SOLN(LAGRANGE)


// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this function.
// This could be macro-ified.
template <> unsigned int FE<0,LAGRANGE>::n_dofs(const ElemType t, const Order o) { return lagrange_n_dofs(t, o); }
template <> unsigned int FE<1,LAGRANGE>::n_dofs(const ElemType t, const Order o) { return lagrange_n_dofs(t, o); }
template <> unsigned int FE<2,LAGRANGE>::n_dofs(const ElemType t, const Order o) { return lagrange_n_dofs(t, o); }
template <> unsigned int FE<3,LAGRANGE>::n_dofs(const ElemType t, const Order o) { return lagrange_n_dofs(t, o); }


// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this function.
template <> unsigned int FE<0,LAGRANGE>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return lagrange_n_dofs_at_node(t, o, n); }
template <> unsigned int FE<1,LAGRANGE>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return lagrange_n_dofs_at_node(t, o, n); }
template <> unsigned int FE<2,LAGRANGE>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return lagrange_n_dofs_at_node(t, o, n); }
template <> unsigned int FE<3,LAGRANGE>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return lagrange_n_dofs_at_node(t, o, n); }


// Lagrange elements have no dofs per element
// (just at the nodes)
template <> unsigned int FE<0,LAGRANGE>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <> unsigned int FE<1,LAGRANGE>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <> unsigned int FE<2,LAGRANGE>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <> unsigned int FE<3,LAGRANGE>::n_dofs_per_elem(const ElemType, const Order) { return 0; }

// Lagrange FEMs are always C^0 continuous
template <> FEContinuity FE<0,LAGRANGE>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<1,LAGRANGE>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<2,LAGRANGE>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<3,LAGRANGE>::get_continuity() const { return C_ZERO; }

// Lagrange FEMs are not hierarchic
template <> bool FE<0,LAGRANGE>::is_hierarchic() const { return false; }
template <> bool FE<1,LAGRANGE>::is_hierarchic() const { return false; }
template <> bool FE<2,LAGRANGE>::is_hierarchic() const { return false; }
template <> bool FE<3,LAGRANGE>::is_hierarchic() const { return false; }

// Lagrange FEM shapes do not need reinit (is this always true?)
template <> bool FE<0,LAGRANGE>::shapes_need_reinit() const { return false; }
template <> bool FE<1,LAGRANGE>::shapes_need_reinit() const { return false; }
template <> bool FE<2,LAGRANGE>::shapes_need_reinit() const { return false; }
template <> bool FE<3,LAGRANGE>::shapes_need_reinit() const { return false; }

// Methods for computing Lagrange constraints.  Note: we pass the
// dimension as the last argument to the anonymous helper function.
// Also note: we only need instantiations of this function for
// Dim==2 and 3.
#ifdef LIBMESH_ENABLE_AMR
template <>
void FE<2,LAGRANGE>::compute_constraints (DofConstraints & constraints,
                                          DofMap & dof_map,
                                          const unsigned int variable_number,
                                          const Elem * elem)
{ lagrange_compute_constraints(constraints, dof_map, variable_number, elem, /*Dim=*/2); }

template <>
void FE<3,LAGRANGE>::compute_constraints (DofConstraints & constraints,
                                          DofMap & dof_map,
                                          const unsigned int variable_number,
                                          const Elem * elem)
{ lagrange_compute_constraints(constraints, dof_map, variable_number, elem, /*Dim=*/3); }
#endif // LIBMESH_ENABLE_AMR

} // namespace libMesh
