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

#ifndef LIBMESH_FE_SZABAB_IMPL_H
#define LIBMESH_FE_SZABAB_IMPL_H

// Local includes
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/enum_to_string.h"

namespace libMesh
{

// ------------------------------------------------------------
// Szabo-Babuska-specific implementations, Steffen Petersen 2004

// Anonymous namespace for local helper functions
namespace {

template <typename RealType>
void szabab_nodal_soln(const ElemTempl<RealType> * elem,
                       const Order order,
                       const std::vector<Number> & elem_soln,
                       std::vector<Number> &       nodal_soln,
                       unsigned Dim)
{
  const unsigned int n_nodes = elem->n_nodes();

  const ElemType elem_type = elem->type();

  nodal_soln.resize(n_nodes);

  const Order totalorder = static_cast<Order>(order+elem->p_level());

  // FEType object to be passed to various FEInterface functions below.
  FEType fe_type(totalorder, SZABAB);

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


      // For other bases do interpolation at the nodes
      // explicitly.
    case FIRST:
    case SECOND:
    case THIRD:
    case FOURTH:
    case FIFTH:
    case SIXTH:
    case SEVENTH:
      {

        const unsigned int n_sf =
          // FE<Dim,T>::n_shape_functions(elem_type, totalorder);
          FEInterface::n_shape_functions(Dim, fe_type, elem_type);

        std::vector<PointTempl<RealType>> refspace_nodes;
        FEGenericBase<Real,RealType>::get_refspace_nodes(elem_type,refspace_nodes);
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
      }

    default:
      libmesh_error_msg("ERROR: Invalid total order " << totalorder);
    }
} // szabab_nodal_soln()




unsigned int szabab_n_dofs(const ElemType t, const Order o)
{
  switch (o)
    {
      // Szabo-Babuska 1st-order polynomials.
    case FIRST:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
            return 2;

          case TRI3:
          case TRI6:
            return 3;

          case QUAD4:
          case QUAD8:
          case QUAD9:
            return 4;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }


      // Szabo-Babuska 2nd-order polynomials.
    case SECOND:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
            return 3;

          case TRI6:
            return 6;

          case QUAD8:
            return 8;
          case QUAD9:
            return 9;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }


      // Szabo-Babuska 3rd-order polynomials.
    case THIRD:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
            return 4;

          case TRI6:
            return 10;

          case QUAD8:
          case QUAD9:
            return 16;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }


      // Szabo-Babuska 4th-order polynomials.
    case FOURTH:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
            return 5;

          case TRI6:
            return 15;

          case QUAD8:
          case QUAD9:
            return 25;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }


      // Szabo-Babuska 5th-order polynomials.
    case FIFTH:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
            return 6;

          case TRI6:
            return 21;

          case QUAD8:
          case QUAD9:
            return 36;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }


      // Szabo-Babuska 6th-order polynomials.
    case SIXTH:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
            return 7;

          case TRI6:
            return 28;

          case QUAD8:
          case QUAD9:
            return 49;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }

      // Szabo-Babuska 7th-order polynomials.
    case SEVENTH:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

          case EDGE2:
          case EDGE3:
            return 8;

          case TRI6:
            return 36;

          case QUAD8:
          case QUAD9:
            return 64;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }


    default:
      libmesh_error_msg("ERROR: Invalid Order " << Utility::enum_to_string(o) << " selected for SZABAB FE family!");
    }
} // szabab_n_dofs()





unsigned int szabab_n_dofs_at_node(const ElemType t,
                                   const Order o,
                                   const unsigned int n)
{
  switch (o)
    {
      // The first-order Szabo-Babuska shape functions
    case FIRST:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

            // The 1D Szabo-Babuska defined on a three-noded edge
          case EDGE2:
          case EDGE3:
            {
              switch (n)
                {
                case 0:
                case 1:
                  return 1;

                case 2:
                  return 0;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for EDGE2/3!");
                }
            }


            // The 2D Szabo-Babuska defined on a 6-noded triangle
          case TRI3:
          case TRI6:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                  return 1;

                case 3:
                case 4:
                case 5:
                  return 0;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for TRI3/6!");
                }
            }


            // The 2D tensor-product Szabo-Babuska defined on a
            // nine-noded quadrilateral.
          case QUAD4:
          case QUAD8:
          case QUAD9:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                case 3:
                  return 1;

                case 4:
                case 5:
                case 6:
                case 7:
                  return 0;

                case 8:
                  return 0;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for QUAD4/8/9!");
                }
            }

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }



      // The second-order Szabo-Babuska shape functions
    case SECOND:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

            // The 1D Szabo-Babuska defined on a three-noded edge
          case EDGE2:
          case EDGE3:
            {
              switch (n)
                {
                case 0:
                case 1:
                  return 1;

                case 2:
                  return 1;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for EDGE2/3!");
                }
            }


            // The 2D Szabo-Babuska defined on a 6-noded triangle
          case TRI6:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                  return 1;

                case 3:
                case 4:
                case 5:
                  return 1;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for TRI6!");
                }
            }


            // The 2D tensor-product Szabo-Babuska defined on a
            // nine-noded quadrilateral.
          case QUAD8:
          case QUAD9:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                case 3:
                  return 1;

                case 4:
                case 5:
                case 6:
                case 7:
                  return 1;

                case 8:
                  return 1;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for QUAD8/9!");
                }
            }

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }



      // The third-order Szabo-Babuska shape functions
    case THIRD:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

            // The 1D Szabo-Babuska defined on a three-noded edge
          case EDGE2:
          case EDGE3:
            {
              switch (n)
                {
                case 0:
                case 1:
                  return 1;

                case 2:
                  return 2;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for EDGE2/3!");
                }
            }


            // The 2D Szabo-Babuska defined on a 6-noded triangle
          case TRI6:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                  return 1;

                case 3:
                case 4:
                case 5:
                  return 2;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for TRI6!");
                }
            }


            // The 2D tensor-product Szabo-Babuska defined on a
            // nine-noded quadrilateral.
          case QUAD8:
          case QUAD9:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                case 3:
                  return 1;

                case 4:
                case 5:
                case 6:
                case 7:
                  return 2;

                case 8:
                  return 4;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for QUAD8/9!");
                }
            }

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }



      // The fourth-order Szabo-Babuska shape functions
    case FOURTH:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

            // The 1D Szabo-Babuska defined on a three-noded edge
          case EDGE2:
          case EDGE3:
            {
              switch (n)
                {
                case 0:
                case 1:
                  return 1;

                case 2:
                  return 3;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for EDGE2/3!");
                }
            }


            // The 2D Szabo-Babuska defined on a 6-noded triangle
          case TRI6:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                  return 1;

                case 3:
                case 4:
                case 5:
                  return 3;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for TRI6!");
                }
            }


            // The 2D tensor-product Szabo-Babuska defined on a
            // nine-noded quadrilateral.
          case QUAD8:
          case QUAD9:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                case 3:
                  return 1;

                case 4:
                case 5:
                case 6:
                case 7:
                  return 3;

                case 8:
                  return 9;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for QUAD8/9!");
                }
            }

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }



      // The fifth-order Szabo-Babuska shape functions
    case FIFTH:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

            // The 1D Szabo-Babuska defined on a three-noded edge
          case EDGE2:
          case EDGE3:
            {
              switch (n)
                {
                case 0:
                case 1:
                  return 1;

                case 2:
                  return 4;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for EDGE2/3!");
                }
            }


            // The 2D Szabo-Babuska defined on a 6-noded triangle
          case TRI6:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                  return 1;

                case 3:
                case 4:
                case 5:
                  return 4;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for TRI6!");
                }
            }


            // The 2D tensor-product Szabo-Babuska defined on a
            // nine-noded quadrilateral.
          case QUAD8:
          case QUAD9:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                case 3:
                  return 1;

                case 4:
                case 5:
                case 6:
                case 7:
                  return 4;

                case 8:
                  return 16;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for QUAD8/9!");
                }
            }

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }



      // The sixth-order Szabo-Babuska shape functions
    case SIXTH:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

            // The 1D Szabo-Babuska defined on a three-noded edge
          case EDGE2:
          case EDGE3:
            {
              switch (n)
                {
                case 0:
                case 1:
                  return 1;

                case 2:
                  return 5;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for EDGE2/3!");
                }
            }


            // The 2D Szabo-Babuska defined on a 6-noded triangle
          case TRI6:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                  return 1;

                case 3:
                case 4:
                case 5:
                  return 5;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for TRI6!");
                }
            }


            // The 2D tensor-product Szabo-Babuska defined on a
            // nine-noded quadrilateral.
          case QUAD8:
          case QUAD9:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                case 3:
                  return 1;

                case 4:
                case 5:
                case 6:
                case 7:
                  return 5;

                case 8:
                  return 25;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for QUAD8/9!");
                }
            }

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }


      // The seventh-order Szabo-Babuska shape functions
    case SEVENTH:
      {
        switch (t)
          {
          case NODEELEM:
            return 1;

            // The 1D Szabo-Babuska defined on a three-noded edge
          case EDGE2:
          case EDGE3:
            {
              switch (n)
                {
                case 0:
                case 1:
                  return 1;

                case 2:
                  return 6;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for EDGE2/3!");
                }
            }


            // The 2D Szabo-Babuska defined on a 6-noded triangle
          case TRI6:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                  return 1;

                case 3:
                case 4:
                case 5:
                  return 6;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for TRI6!");
                }
            }


            // The 2D tensor-product Szabo-Babuska defined on a
            // nine-noded quadrilateral.
          case QUAD8:
          case QUAD9:
            {
              switch (n)
                {
                case 0:
                case 1:
                case 2:
                case 3:
                  return 1;

                case 4:
                case 5:
                case 6:
                case 7:
                  return 6;

                case 8:
                  return 36;

                default:
                  libmesh_error_msg("ERROR: Invalid node ID " << n << " selected for QUAD8/9!");
                }
            }

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }


    default:
      libmesh_error_msg("ERROR: Invalid Order " << Utility::enum_to_string(o) << " selected for SZABAB FE family!");
    }
} // szabab_n_dofs_at_node()



unsigned int szabab_n_dofs_per_elem(const ElemType t, const Order o)
{
  switch (o)
    {
      // The first-order Szabo-Babuska shape functions
    case FIRST:
      {
        switch (t)
          {
          case NODEELEM:
            return 0;

            // The 1D Szabo-Babuska defined on an edge
          case EDGE2:
          case EDGE3:
            return 0;

            // The 2D Szabo-Babuska defined on a triangle
          case TRI3:
          case TRI6:
            return 0;

            // The 2D tensor-product Szabo-Babuska defined on a
            // quadrilateral.
          case QUAD4:
          case QUAD8:
          case QUAD9:
            return 0;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }



      // The second-order Szabo-Babuska shape functions
    case SECOND:
      {
        switch (t)
          {
          case NODEELEM:
            return 0;

            // The 1D Szabo-Babuska defined on a two-noded edge
          case EDGE2:
            return 1;

            // The 1D Szabo-Babuska defined on a three-noded edge
          case EDGE3:
            return 0;

            // The 2D Szabo-Babuska defined on a 6-noded triangle
          case TRI6:
            return 0;

            // The 2D tensor-product Szabo-Babuska defined on a
            // eight-noded quadrilateral.
          case QUAD8:
            return 0;

            // The 2D tensor-product Szabo-Babuska defined on a
            // nine-noded quadrilateral.
          case QUAD9:
            return 0;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }



      // The third-order Szabo-Babuska shape functions
    case THIRD:
      {
        switch (t)
          {
          case NODEELEM:
            return 0;

            // The 1D Szabo-Babuska defined on a two-noded edge
          case EDGE2:
            return 2;

            // The 1D Szabo-Babuska defined on a three-noded edge
          case EDGE3:
            return 0;

            // The 2D Szabo-Babuska defined on a 6-noded triangle
          case TRI6:
            return 1;

            // The 2D tensor-product Szabo-Babuska defined on a
            // eight-noded quadrilateral.
          case QUAD8:
            return 4;

            // The 2D tensor-product Szabo-Babuska defined on a
            // nine-noded quadrilateral.
          case QUAD9:
            return 0;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }



      // The fourth-order Szabo-Babuska shape functions
    case FOURTH:
      {
        switch (t)
          {
          case NODEELEM:
            return 0;

            // The 1D Szabo-Babuska defined on a two-noded edge
          case EDGE2:
            return 3;

            // The 1D Szabo-Babuska defined on a three-noded edge
          case EDGE3:
            return 0;

            // The 2D Szabo-Babuska defined on a 6-noded triangle
          case TRI6:
            return 3;

            // The 2D tensor-product Szabo-Babuska defined on a
            // eight-noded quadrilateral.
          case QUAD8:
            return 9;

            // The 2D tensor-product Szabo-Babuska defined on a
            // nine-noded quadrilateral.
          case QUAD9:
            return 0;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }



      // The fifth-order Szabo-Babuska shape functions
    case FIFTH:
      {
        switch (t)
          {
          case NODEELEM:
            return 0;

            // The 1D Szabo-Babuska defined on a two-noded edge
          case EDGE2:
            return 4;

            // The 1D Szabo-Babuska defined on a three-noded edge
          case EDGE3:
            return 0;

            // The 2D Szabo-Babuska defined on a 6-noded triangle
          case TRI6:
            return 6;

            // The 2D tensor-product Szabo-Babuska defined on a
            // eight-noded quadrilateral.
          case QUAD8:
            return 16;

            // The 2D tensor-product Szabo-Babuska defined on a
            // nine-noded quadrilateral.
          case QUAD9:
            return 0;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }


      // The sixth-order Szabo-Babuska shape functions
    case SIXTH:
      {
        switch (t)
          {
          case NODEELEM:
            return 0;

            // The 1D Szabo-Babuska defined on a two-noded edge
          case EDGE2:
            return 5;

            // The 1D Szabo-Babuska defined on a three-noded edge
          case EDGE3:
            return 0;

            // The 2D Szabo-Babuska defined on a 6-noded triangle
          case TRI6:
            return 10;

            // The 2D tensor-product Szabo-Babuska defined on a
            // eight-noded quadrilateral.
          case QUAD8:
            return 25;

            // The 2D tensor-product Szabo-Babuska defined on a
            // nine-noded quadrilateral.
          case QUAD9:
            return 0;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }


      // The seventh-order Szabo-Babuska shape functions
    case SEVENTH:
      {
        switch (t)
          {
          case NODEELEM:
            return 0;

            // The 1D Szabo-Babuska defined on a two-noded edge
          case EDGE2:
            return 6;

            // The 1D Szabo-Babuska defined on a three-noded edge
          case EDGE3:
            return 0;

            // The 2D Szabo-Babuska defined on a 6-noded triangle
          case TRI6:
            return 15;

            // The 2D tensor-product Szabo-Babuska defined on a
            // eight-noded quadrilateral.
          case QUAD8:
            return 36;

            // The 2D tensor-product Szabo-Babuska defined on a
            // nine-noded quadrilateral.
          case QUAD9:
            return 0;

          case INVALID_ELEM:
            return 0;

          default:
            libmesh_error_msg("ERROR: Invalid ElemType " << Utility::enum_to_string(t) << " selected for SZABAB FE family!");
          }
      }


      // Otherwise no DOFS per element
    default:
      return 0;
    }
} // szabab_n_dofs_per_elem

} // anonymous namespace




  // Do full-specialization of nodal_soln() function for every
  // dimension, instead of explicit instantiation at the end of this
  // file.
  // This could be macro-ified so that it fits on one line...
template <typename RealType>
void FEShim<0,SZABAB,RealType>::nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln)
{ szabab_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/0); }

template <typename RealType>
void FEShim<1,SZABAB,RealType>::nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln)
{ szabab_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/1); }

template <typename RealType>
void FEShim<2,SZABAB,RealType>::nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln)
{ szabab_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/2); }

template <typename RealType>
void FEShim<3,SZABAB,RealType>::nodal_soln(const Elem * elem,
                              const Order order,
                              const std::vector<Number> & elem_soln,
                              std::vector<Number> & nodal_soln)
{ szabab_nodal_soln(elem, order, elem_soln, nodal_soln, /*Dim=*/3); }


// Full specialization of n_dofs() function for every dimension
template <typename RealType> unsigned int FEShim<0,SZABAB,RealType>::n_dofs(const ElemType t, const Order o) { return szabab_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<1,SZABAB,RealType>::n_dofs(const ElemType t, const Order o) { return szabab_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<2,SZABAB,RealType>::n_dofs(const ElemType t, const Order o) { return szabab_n_dofs(t, o); }
template <typename RealType> unsigned int FEShim<3,SZABAB,RealType>::n_dofs(const ElemType t, const Order o) { return szabab_n_dofs(t, o); }

// Full specialization of n_dofs_at_node() function for every dimension.
template <typename RealType> unsigned int FEShim<0,SZABAB,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return szabab_n_dofs_at_node(t, o, n); }
template <typename RealType> unsigned int FEShim<1,SZABAB,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return szabab_n_dofs_at_node(t, o, n); }
template <typename RealType> unsigned int FEShim<2,SZABAB,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return szabab_n_dofs_at_node(t, o, n); }
template <typename RealType> unsigned int FEShim<3,SZABAB,RealType>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return szabab_n_dofs_at_node(t, o, n); }

// Full specialization of n_dofs_per_elem() function for every dimension.
template <typename RealType> unsigned int FEShim<0,SZABAB,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return szabab_n_dofs_per_elem(t, o); }
template <typename RealType> unsigned int FEShim<1,SZABAB,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return szabab_n_dofs_per_elem(t, o); }
template <typename RealType> unsigned int FEShim<2,SZABAB,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return szabab_n_dofs_per_elem(t, o); }
template <typename RealType> unsigned int FEShim<3,SZABAB,RealType>::n_dofs_per_elem(const ElemType t, const Order o) { return szabab_n_dofs_per_elem(t, o); }

// Szabab FEMs are C^0 continuous
template <typename RealType> FEContinuity FEShim<0,SZABAB,RealType>::get_continuity() { return C_ZERO; }
template <typename RealType> FEContinuity FEShim<1,SZABAB,RealType>::get_continuity() { return C_ZERO; }
template <typename RealType> FEContinuity FEShim<2,SZABAB,RealType>::get_continuity() { return C_ZERO; }
template <typename RealType> FEContinuity FEShim<3,SZABAB,RealType>::get_continuity() { return C_ZERO; }

// Szabab FEMs are hierarchic
template <typename RealType> bool FEShim<0,SZABAB,RealType>::is_hierarchic() { return true; }
template <typename RealType> bool FEShim<1,SZABAB,RealType>::is_hierarchic() { return true; }
template <typename RealType> bool FEShim<2,SZABAB,RealType>::is_hierarchic() { return true; }
template <typename RealType> bool FEShim<3,SZABAB,RealType>::is_hierarchic() { return true; }

#ifdef LIBMESH_ENABLE_AMR
// compute_constraints() specializations are only needed for 2 and 3D
template <typename RealType>
void FEShim<2,SZABAB,RealType>::compute_constraints (DofConstraints & constraints,
                                        DofMap & dof_map,
                                        const unsigned int variable_number,
                                                     const ElemTempl<Real> * elem)
{ FEGenericBase<typename FEOutputType<SZABAB>::type, RealType>::compute_proj_constraints(constraints, dof_map, variable_number, elem); }

template <typename RealType>
void FEShim<3,SZABAB,RealType>::compute_constraints (DofConstraints & constraints,
                                        DofMap & dof_map,
                                        const unsigned int variable_number,
                                                     const ElemTempl<Real> * elem)
{ FEGenericBase<typename FEOutputType<SZABAB>::type, RealType>::compute_proj_constraints(constraints, dof_map, variable_number, elem); }
#endif // #ifdef LIBMESH_ENABLE_AMR

// Szabab shapes need reinit only for approximation orders >= 3,
// but we might reach that with p refinement
template <typename RealType> bool FEShim<0,SZABAB,RealType>::shapes_need_reinit() { return true; }
template <typename RealType> bool FEShim<1,SZABAB,RealType>::shapes_need_reinit() { return true; }
template <typename RealType> bool FEShim<2,SZABAB,RealType>::shapes_need_reinit() { return true; }
template <typename RealType> bool FEShim<3,SZABAB,RealType>::shapes_need_reinit() { return true; }

} // namespace libMesh

#endif //LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#endif // LIBMESH_FE_SZABAB_IMPL_H
