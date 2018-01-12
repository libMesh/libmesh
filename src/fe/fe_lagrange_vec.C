// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
// Lagrange-specific implementations


// Anonymous namespace for local helper functions
namespace {
void lagrange_vec_nodal_soln(const Elem * elem,
                             const Order order,
                             const std::vector<Number> & elem_soln,
                             const int dim,
                             std::vector<Number> &       nodal_soln)
{
  const unsigned int n_nodes = elem->n_nodes();
  const ElemType type        = elem->type();

  const Order totalorder = static_cast<Order>(order+elem->p_level());

  nodal_soln.resize(dim*n_nodes);

  switch (totalorder)
    {
      // linear Lagrange shape functions
    case FIRST:
      {
        switch (type)
          {
          case TRI6:
            {
              libmesh_assert_equal_to (elem_soln.size(), 2*3);
              libmesh_assert_equal_to (nodal_soln.size(), 2*6);

              // node 0 components
              nodal_soln[0] = elem_soln[0];
              nodal_soln[1] = elem_soln[1];

              // node 1 components
              nodal_soln[2] = elem_soln[2];
              nodal_soln[3] = elem_soln[3];

              // node 2 components
              nodal_soln[4] = elem_soln[4];
              nodal_soln[5] = elem_soln[5];

              // node 3 components
              nodal_soln[6] = .5*(elem_soln[0] + elem_soln[2]);
              nodal_soln[7] = .5*(elem_soln[1] + elem_soln[3]);

              // node 4 components
              nodal_soln[8] = .5*(elem_soln[2] + elem_soln[4]);
              nodal_soln[9] = .5*(elem_soln[3] + elem_soln[5]);

              // node 5 components
              nodal_soln[10] = .5*(elem_soln[0] + elem_soln[4]);
              nodal_soln[11] = .5*(elem_soln[1] + elem_soln[5]);

              return;
            }


          case QUAD8:
          case QUAD9:
            {
              libmesh_assert_equal_to (elem_soln.size(), 2*4);

              if (type == QUAD8)
                libmesh_assert_equal_to (nodal_soln.size(), 2*8);
              else
                libmesh_assert_equal_to (nodal_soln.size(), 2*9);

              // node 0 components
              nodal_soln[0] = elem_soln[0];
              nodal_soln[1] = elem_soln[1];

              // node 1 components
              nodal_soln[2] = elem_soln[2];
              nodal_soln[3] = elem_soln[3];

              // node 2 components
              nodal_soln[4] = elem_soln[4];
              nodal_soln[5] = elem_soln[5];

              // node 3 components
              nodal_soln[6] = elem_soln[6];
              nodal_soln[7] = elem_soln[7];

              // node 4 components
              nodal_soln[8] = .5*(elem_soln[0] + elem_soln[2]);
              nodal_soln[9] = .5*(elem_soln[1] + elem_soln[3]);

              // node 5 components
              nodal_soln[10] = .5*(elem_soln[2] + elem_soln[4]);
              nodal_soln[11] = .5*(elem_soln[3] + elem_soln[5]);

              // node 6 components
              nodal_soln[12] = .5*(elem_soln[4] + elem_soln[6]);
              nodal_soln[13] = .5*(elem_soln[5] + elem_soln[7]);

              // node 7 components
              nodal_soln[14] = .5*(elem_soln[6] + elem_soln[0]);
              nodal_soln[15] = .5*(elem_soln[7] + elem_soln[1]);

              if (type == QUAD9)
                {
                  // node 8 components
                  nodal_soln[16] = .25*(elem_soln[0] + elem_soln[2] + elem_soln[4] + elem_soln[6]);
                  nodal_soln[17] = .25*(elem_soln[1] + elem_soln[3] + elem_soln[5] + elem_soln[7]);
                }

              return;
            }


          case TET10:
            {
              libmesh_assert_equal_to (elem_soln.size(), 3*4);
              libmesh_assert_equal_to (nodal_soln.size(), 3*10);

              // node 0 components
              nodal_soln[0] = elem_soln[0];
              nodal_soln[1] = elem_soln[1];
              nodal_soln[2] = elem_soln[2];

              // node 1 components
              nodal_soln[3] = elem_soln[3];
              nodal_soln[4] = elem_soln[4];
              nodal_soln[5] = elem_soln[5];

              // node 2 components
              nodal_soln[6] = elem_soln[6];
              nodal_soln[7] = elem_soln[7];
              nodal_soln[8] = elem_soln[8];

              // node 3 components
              nodal_soln[9]  = elem_soln[9];
              nodal_soln[10] = elem_soln[10];
              nodal_soln[11] = elem_soln[11];

              // node 4 components
              nodal_soln[12] = .5*(elem_soln[0] + elem_soln[3]);
              nodal_soln[13] = .5*(elem_soln[1] + elem_soln[4]);
              nodal_soln[14] = .5*(elem_soln[2] + elem_soln[5]);

              // node 5 components
              nodal_soln[15] = .5*(elem_soln[3] + elem_soln[6]);
              nodal_soln[16] = .5*(elem_soln[4] + elem_soln[7]);
              nodal_soln[17] = .5*(elem_soln[5] + elem_soln[8]);

              // node 6 components
              nodal_soln[18] = .5*(elem_soln[6] + elem_soln[0]);
              nodal_soln[19] = .5*(elem_soln[7] + elem_soln[1]);
              nodal_soln[20] = .5*(elem_soln[8] + elem_soln[2]);

              // node 7 components
              nodal_soln[21] = .5*(elem_soln[9]  + elem_soln[0]);
              nodal_soln[22] = .5*(elem_soln[10] + elem_soln[1]);
              nodal_soln[23] = .5*(elem_soln[11] + elem_soln[2]);

              // node 8 components
              nodal_soln[24] = .5*(elem_soln[9]  + elem_soln[3]);
              nodal_soln[25] = .5*(elem_soln[10] + elem_soln[4]);
              nodal_soln[26] = .5*(elem_soln[11] + elem_soln[5]);

              // node 9 components
              nodal_soln[27] = .5*(elem_soln[9]  + elem_soln[6]);
              nodal_soln[28] = .5*(elem_soln[10] + elem_soln[7]);
              nodal_soln[29] = .5*(elem_soln[11] + elem_soln[8]);

              return;
            }


          case HEX20:
          case HEX27:
            {
              libmesh_assert_equal_to (elem_soln.size(), 3*8);

              if (type == HEX20)
                libmesh_assert_equal_to (nodal_soln.size(), 3*20);
              else
                libmesh_assert_equal_to (nodal_soln.size(), 3*27);

              // node 0 components
              nodal_soln[0]  = elem_soln[0];
              nodal_soln[1]  = elem_soln[1];
              nodal_soln[2]  = elem_soln[2];

              // node 1 components
              nodal_soln[3]  = elem_soln[3];
              nodal_soln[4]  = elem_soln[4];
              nodal_soln[5]  = elem_soln[5];

              // node 2 components
              nodal_soln[6]  = elem_soln[6];
              nodal_soln[7]  = elem_soln[7];
              nodal_soln[8]  = elem_soln[8];

              // node 3 components
              nodal_soln[9]   = elem_soln[9];
              nodal_soln[10]  = elem_soln[10];
              nodal_soln[11]  = elem_soln[11];

              // node 4 components
              nodal_soln[12]  = elem_soln[12];
              nodal_soln[13]  = elem_soln[13];
              nodal_soln[14]  = elem_soln[14];

              // node 5 components
              nodal_soln[15]  = elem_soln[15];
              nodal_soln[16]  = elem_soln[16];
              nodal_soln[17]  = elem_soln[17];

              // node 6 components
              nodal_soln[18]  = elem_soln[18];
              nodal_soln[19]  = elem_soln[19];
              nodal_soln[20]  = elem_soln[20];

              // node 7 components
              nodal_soln[21]  = elem_soln[21];
              nodal_soln[22]  = elem_soln[22];
              nodal_soln[23]  = elem_soln[23];

              // node 8 components
              nodal_soln[24]  = .5*(elem_soln[0] + elem_soln[3]);
              nodal_soln[25]  = .5*(elem_soln[1] + elem_soln[4]);
              nodal_soln[26]  = .5*(elem_soln[2] + elem_soln[5]);

              // node 9 components
              nodal_soln[27]  = .5*(elem_soln[3] + elem_soln[6]);
              nodal_soln[28]  = .5*(elem_soln[4] + elem_soln[7]);
              nodal_soln[29]  = .5*(elem_soln[5] + elem_soln[8]);

              // node 10 components
              nodal_soln[30]  = .5*(elem_soln[6] + elem_soln[9]);
              nodal_soln[31]  = .5*(elem_soln[7] + elem_soln[10]);
              nodal_soln[32]  = .5*(elem_soln[8] + elem_soln[11]);

              // node 11 components
              nodal_soln[33]  = .5*(elem_soln[9]  + elem_soln[0]);
              nodal_soln[34]  = .5*(elem_soln[10] + elem_soln[1]);
              nodal_soln[35]  = .5*(elem_soln[11] + elem_soln[2]);

              // node 12 components
              nodal_soln[36]  = .5*(elem_soln[0] + elem_soln[12]);
              nodal_soln[37]  = .5*(elem_soln[1] + elem_soln[13]);
              nodal_soln[38]  = .5*(elem_soln[2] + elem_soln[14]);

              // node 13 components
              nodal_soln[39]  = .5*(elem_soln[3] + elem_soln[15]);
              nodal_soln[40]  = .5*(elem_soln[4] + elem_soln[16]);
              nodal_soln[41]  = .5*(elem_soln[5] + elem_soln[17]);

              // node 14 components
              nodal_soln[42]  = .5*(elem_soln[6] + elem_soln[18]);
              nodal_soln[43]  = .5*(elem_soln[7] + elem_soln[19]);
              nodal_soln[44]  = .5*(elem_soln[8] + elem_soln[20]);

              // node 15 components
              nodal_soln[45]  = .5*(elem_soln[9]  + elem_soln[21]);
              nodal_soln[46]  = .5*(elem_soln[10] + elem_soln[22]);
              nodal_soln[47]  = .5*(elem_soln[11] + elem_soln[23]);

              // node 16 components
              nodal_soln[48]  = .5*(elem_soln[12] + elem_soln[15]);
              nodal_soln[49]  = .5*(elem_soln[13] + elem_soln[16]);
              nodal_soln[50]  = .5*(elem_soln[14] + elem_soln[17]);

              // node 17 components
              nodal_soln[51]  = .5*(elem_soln[15] + elem_soln[18]);
              nodal_soln[52]  = .5*(elem_soln[16] + elem_soln[19]);
              nodal_soln[53]  = .5*(elem_soln[17] + elem_soln[20]);

              // node 18 components
              nodal_soln[54]  = .5*(elem_soln[18] + elem_soln[21]);
              nodal_soln[55]  = .5*(elem_soln[19] + elem_soln[22]);
              nodal_soln[56]  = .5*(elem_soln[20] + elem_soln[23]);

              // node 19 components
              nodal_soln[57]  = .5*(elem_soln[12] + elem_soln[21]);
              nodal_soln[58]  = .5*(elem_soln[13] + elem_soln[22]);
              nodal_soln[59]  = .5*(elem_soln[14] + elem_soln[23]);

              if (type == HEX27)
                {
                  // node 20 components
                  nodal_soln[60]  = .25*(elem_soln[0] + elem_soln[3] + elem_soln[6] + elem_soln[9]);
                  nodal_soln[61]  = .25*(elem_soln[1] + elem_soln[4] + elem_soln[7] + elem_soln[10]);
                  nodal_soln[62]  = .25*(elem_soln[2] + elem_soln[5] + elem_soln[8] + elem_soln[11]);

                  // node 21 components
                  nodal_soln[63]  = .25*(elem_soln[0] + elem_soln[3] + elem_soln[12] + elem_soln[15]);
                  nodal_soln[64]  = .25*(elem_soln[1] + elem_soln[4] + elem_soln[13] + elem_soln[16]);
                  nodal_soln[65]  = .25*(elem_soln[2] + elem_soln[5] + elem_soln[14] + elem_soln[17]);

                  // node 22 components
                  nodal_soln[66]  = .25*(elem_soln[3] + elem_soln[6] + elem_soln[15] + elem_soln[18]);
                  nodal_soln[67]  = .25*(elem_soln[4] + elem_soln[7] + elem_soln[16] + elem_soln[19]);
                  nodal_soln[68]  = .25*(elem_soln[5] + elem_soln[8] + elem_soln[17] + elem_soln[20]);

                  // node 23 components
                  nodal_soln[69]  = .25*(elem_soln[6] + elem_soln[9]  + elem_soln[18] + elem_soln[21]);
                  nodal_soln[70]  = .25*(elem_soln[7] + elem_soln[10] + elem_soln[19] + elem_soln[22]);
                  nodal_soln[71]  = .25*(elem_soln[8] + elem_soln[11] + elem_soln[20] + elem_soln[23]);

                  // node 24 components
                  nodal_soln[72]  = .25*(elem_soln[9]  + elem_soln[0] + elem_soln[21] + elem_soln[12]);
                  nodal_soln[73]  = .25*(elem_soln[10] + elem_soln[1] + elem_soln[22] + elem_soln[13]);
                  nodal_soln[74]  = .25*(elem_soln[11] + elem_soln[2] + elem_soln[23] + elem_soln[14]);

                  // node 25 components
                  nodal_soln[75]  = .25*(elem_soln[12] + elem_soln[15] + elem_soln[18] + elem_soln[21]);
                  nodal_soln[76]  = .25*(elem_soln[13] + elem_soln[16] + elem_soln[19] + elem_soln[22]);
                  nodal_soln[77]  = .25*(elem_soln[14] + elem_soln[17] + elem_soln[20] + elem_soln[23]);

                  // node 26 components
                  nodal_soln[78]  = .125*(elem_soln[0]  + elem_soln[3]  + elem_soln[6]  + elem_soln[9] +
                                          elem_soln[12] + elem_soln[15] + elem_soln[18] + elem_soln[21]);

                  nodal_soln[79]  = .125*(elem_soln[1]  + elem_soln[4]  + elem_soln[7]  + elem_soln[10] +
                                          elem_soln[13] + elem_soln[16] + elem_soln[19] + elem_soln[22]);

                  nodal_soln[80]  = .125*(elem_soln[2]  + elem_soln[5]  + elem_soln[8]  + elem_soln[11] +
                                          elem_soln[14] + elem_soln[17] + elem_soln[20] + elem_soln[23]);
                }

              return;
            }


          case PRISM15:
          case PRISM18:
            {
              libmesh_assert_equal_to (elem_soln.size(), 3*6);

              if (type == PRISM15)
                libmesh_assert_equal_to (nodal_soln.size(), 3*15);
              else
                libmesh_assert_equal_to (nodal_soln.size(), 3*18);

              // node 0 components
              nodal_soln[0]  = elem_soln[0];
              nodal_soln[1]  = elem_soln[1];
              nodal_soln[2]  = elem_soln[2];

              // node 1 components
              nodal_soln[3]  = elem_soln[3];
              nodal_soln[4]  = elem_soln[4];
              nodal_soln[5]  = elem_soln[5];

              // node 2 components
              nodal_soln[6]  = elem_soln[6];
              nodal_soln[7]  = elem_soln[7];
              nodal_soln[8]  = elem_soln[8];

              // node 3 components
              nodal_soln[9]   = elem_soln[9];
              nodal_soln[10]  = elem_soln[10];
              nodal_soln[11]  = elem_soln[11];

              // node 4 components
              nodal_soln[12]  = elem_soln[12];
              nodal_soln[13]  = elem_soln[13];
              nodal_soln[14]  = elem_soln[14];

              // node 5 components
              nodal_soln[15]  = elem_soln[15];
              nodal_soln[16]  = elem_soln[16];
              nodal_soln[17]  = elem_soln[17];

              // node 6 components
              nodal_soln[18]  = .5*(elem_soln[0] + elem_soln[3]);
              nodal_soln[19]  = .5*(elem_soln[1] + elem_soln[4]);
              nodal_soln[20]  = .5*(elem_soln[2] + elem_soln[5]);

              // node 7 components
              nodal_soln[21]  = .5*(elem_soln[3] + elem_soln[6]);
              nodal_soln[22]  = .5*(elem_soln[4] + elem_soln[7]);
              nodal_soln[23]  = .5*(elem_soln[5] + elem_soln[8]);

              // node 8 components
              nodal_soln[24]  = .5*(elem_soln[0] + elem_soln[6]);
              nodal_soln[25]  = .5*(elem_soln[1] + elem_soln[7]);
              nodal_soln[26]  = .5*(elem_soln[2] + elem_soln[8]);

              // node 9 components
              nodal_soln[27]  = .5*(elem_soln[0] + elem_soln[9]);
              nodal_soln[28]  = .5*(elem_soln[1] + elem_soln[10]);
              nodal_soln[29]  = .5*(elem_soln[2] + elem_soln[11]);

              // node 10 components
              nodal_soln[30]  = .5*(elem_soln[3] + elem_soln[12]);
              nodal_soln[31]  = .5*(elem_soln[4] + elem_soln[13]);
              nodal_soln[32]  = .5*(elem_soln[5] + elem_soln[14]);

              // node 11 components
              nodal_soln[33]  = .5*(elem_soln[6] + elem_soln[15]);
              nodal_soln[34]  = .5*(elem_soln[7] + elem_soln[16]);
              nodal_soln[35]  = .5*(elem_soln[8] + elem_soln[17]);

              // node 12 components
              nodal_soln[36]  = .5*(elem_soln[9]  + elem_soln[12]);
              nodal_soln[37]  = .5*(elem_soln[10] + elem_soln[13]);
              nodal_soln[38]  = .5*(elem_soln[11] + elem_soln[14]);

              // node 13 components
              nodal_soln[39]  = .5*(elem_soln[12] + elem_soln[15]);
              nodal_soln[40]  = .5*(elem_soln[13] + elem_soln[16]);
              nodal_soln[41]  = .5*(elem_soln[14] + elem_soln[17]);

              // node 14 components
              nodal_soln[42]  = .5*(elem_soln[12] + elem_soln[15]);
              nodal_soln[43]  = .5*(elem_soln[13] + elem_soln[16]);
              nodal_soln[44]  = .5*(elem_soln[14] + elem_soln[17]);

              if (type == PRISM18)
                {
                  // node 15 components
                  nodal_soln[45]  = .25*(elem_soln[0] + elem_soln[3] + elem_soln[12] + elem_soln[9]);
                  nodal_soln[46]  = .25*(elem_soln[1] + elem_soln[4] + elem_soln[13] + elem_soln[10]);
                  nodal_soln[47]  = .25*(elem_soln[2] + elem_soln[5] + elem_soln[14] + elem_soln[11]);

                  // node 16 components
                  nodal_soln[48]  = .25*(elem_soln[3] + elem_soln[6] + elem_soln[15] + elem_soln[12]);
                  nodal_soln[49]  = .25*(elem_soln[4] + elem_soln[7] + elem_soln[16] + elem_soln[13]);
                  nodal_soln[50]  = .25*(elem_soln[5] + elem_soln[8] + elem_soln[17] + elem_soln[14]);

                  // node 17 components
                  nodal_soln[51]  = .25*(elem_soln[6] + elem_soln[0] + elem_soln[9]  + elem_soln[15]);
                  nodal_soln[52]  = .25*(elem_soln[7] + elem_soln[1] + elem_soln[10] + elem_soln[16]);
                  nodal_soln[53]  = .25*(elem_soln[8] + elem_soln[2] + elem_soln[11] + elem_soln[17]);
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

      }

    } // switch(totalorder)

}// void lagrange_vec_nodal_soln

} // anonymous namespace


  // Do full-specialization for every dimension, instead
  // of explicit instantiation at the end of this file.
  // This could be macro-ified so that it fits on one line...
template <>
void FE<0,LAGRANGE_VEC>::nodal_soln(const Elem * elem,
                                    const Order order,
                                    const std::vector<Number> & elem_soln,
                                    std::vector<Number> & nodal_soln)
{ FE<0,LAGRANGE>::nodal_soln(elem, order, elem_soln, nodal_soln); }

template <>
void FE<1,LAGRANGE_VEC>::nodal_soln(const Elem * elem,
                                    const Order order,
                                    const std::vector<Number> & elem_soln,
                                    std::vector<Number> & nodal_soln)
{ FE<1,LAGRANGE>::nodal_soln(elem, order, elem_soln, nodal_soln); }

template <>
void FE<2,LAGRANGE_VEC>::nodal_soln(const Elem * elem,
                                    const Order order,
                                    const std::vector<Number> & elem_soln,
                                    std::vector<Number> & nodal_soln)
{ lagrange_vec_nodal_soln(elem, order, elem_soln, 2 /*dimension*/, nodal_soln); }

template <>
void FE<3,LAGRANGE_VEC>::nodal_soln(const Elem * elem,
                                    const Order order,
                                    const std::vector<Number> & elem_soln,
                                    std::vector<Number> & nodal_soln)
{ lagrange_vec_nodal_soln(elem, order, elem_soln, 3 /*dimension*/, nodal_soln); }


// Specialize for shape function routines by leveraging scalar LAGRANGE elements

// 0-D
template <> RealGradient FE<0,LAGRANGE_VEC>::shape(const ElemType type, const Order order,
                                                   const unsigned int i, const Point & p)
{
  Real value = FE<0,LAGRANGE>::shape( type, order, i, p );
  return libMesh::RealGradient( value );
}
template <> RealGradient FE<0,LAGRANGE_VEC>::shape_deriv(const ElemType type, const Order order,
                                                         const unsigned int i, const unsigned int j,
                                                         const Point & p)
{
  Real value = FE<0,LAGRANGE>::shape_deriv( type, order, i, j, p );
  return libMesh::RealGradient( value );
}
template <> RealGradient FE<0,LAGRANGE_VEC>::shape_second_deriv(const ElemType type, const Order order,
                                                                const unsigned int i, const unsigned int j,
                                                                const Point & p)
{
  Real value = FE<0,LAGRANGE>::shape_second_deriv( type, order, i, j, p );
  return libMesh::RealGradient( value );
}

// 1-D
template <> RealGradient FE<1,LAGRANGE_VEC>::shape(const ElemType type, const Order order,
                                                   const unsigned int i, const Point & p)
{
  Real value = FE<1,LAGRANGE>::shape( type, order, i, p );
  return libMesh::RealGradient( value );
}
template <> RealGradient FE<1,LAGRANGE_VEC>::shape_deriv(const ElemType type, const Order order,
                                                         const unsigned int i, const unsigned int j,
                                                         const Point & p)
{
  Real value = FE<1,LAGRANGE>::shape_deriv( type, order, i, j, p );
  return libMesh::RealGradient( value );
}
template <> RealGradient FE<1,LAGRANGE_VEC>::shape_second_deriv(const ElemType type, const Order order,
                                                                const unsigned int i, const unsigned int j,
                                                                const Point & p)
{
  Real value = FE<1,LAGRANGE>::shape_second_deriv( type, order, i, j, p );
  return libMesh::RealGradient( value );
}

// 2-D
template <> RealGradient FE<2,LAGRANGE_VEC>::shape(const ElemType type, const Order order,
                                                   const unsigned int i, const Point & p)
{
  Real value = FE<2,LAGRANGE>::shape( type, order, i/2, p );

  switch( i%2 )
    {
    case 0:
      return libMesh::RealGradient( value );

    case 1:
      return libMesh::RealGradient( Real(0), value );

    default:
      libmesh_error_msg("i%2 must be either 0 or 1!");
    }

  //dummy
  return libMesh::RealGradient();
}
template <> RealGradient FE<2,LAGRANGE_VEC>::shape_deriv(const ElemType type, const Order order,
                                                         const unsigned int i, const unsigned int j,
                                                         const Point & p)
{
  Real value = FE<2,LAGRANGE>::shape_deriv( type, order, i/2, j, p );

  switch( i%2 )
    {
    case 0:
      return libMesh::RealGradient( value );

    case 1:
      return libMesh::RealGradient( Real(0), value );

    default:
      libmesh_error_msg("i%2 must be either 0 or 1!");
    }

  //dummy
  return libMesh::RealGradient();
}
template <> RealGradient FE<2,LAGRANGE_VEC>::shape_second_deriv(const ElemType type, const Order order,
                                                                const unsigned int i, const unsigned int j,
                                                                const Point & p)
{
  Real value = FE<2,LAGRANGE>::shape_second_deriv( type, order, i/2, j, p );

  switch( i%2 )
    {
    case 0:
      return libMesh::RealGradient( value );

    case 1:
      return libMesh::RealGradient( Real(0), value );

    default:
      libmesh_error_msg("i%2 must be either 0 or 1!");
    }

  //dummy
  return libMesh::RealGradient();
}


// 3-D
template <> RealGradient FE<3,LAGRANGE_VEC>::shape(const ElemType type, const Order order,
                                                   const unsigned int i, const Point & p)
{
  Real value = FE<3,LAGRANGE>::shape( type, order, i/3, p );

  switch( i%3 )
    {
    case 0:
      return libMesh::RealGradient( value );

    case 1:
      return libMesh::RealGradient( Real(0), value );

    case 2:
      return libMesh::RealGradient( Real(0), Real(0), value );

    default:
      libmesh_error_msg("i%3 must be 0, 1, or 2!");
    }

  //dummy
  return libMesh::RealGradient();
}
template <> RealGradient FE<3,LAGRANGE_VEC>::shape_deriv(const ElemType type, const Order order,
                                                         const unsigned int i, const unsigned int j,
                                                         const Point & p)
{
  Real value = FE<3,LAGRANGE>::shape_deriv( type, order, i/3, j, p );

  switch( i%3 )
    {
    case 0:
      return libMesh::RealGradient( value );

    case 1:
      return libMesh::RealGradient( Real(0), value );

    case 2:
      return libMesh::RealGradient( Real(0), Real(0), value );

    default:
      libmesh_error_msg("i%3 must be 0, 1, or 2!");
    }

  //dummy
  return libMesh::RealGradient();
}
template <> RealGradient FE<3,LAGRANGE_VEC>::shape_second_deriv(const ElemType type, const Order order,
                                                                const unsigned int i, const unsigned int j,
                                                                const Point & p)
{
  Real value = FE<3,LAGRANGE>::shape_second_deriv( type, order, i/3, j, p );

  switch( i%3 )
    {
    case 0:
      return libMesh::RealGradient( value );

    case 1:
      return libMesh::RealGradient( Real(0), value );

    case 2:
      return libMesh::RealGradient( Real(0), Real(0), value );

    default:
      libmesh_error_msg("i%3 must be 0, 1, or 2!");
    }

  //dummy
  return libMesh::RealGradient();
}



// 0-D
template <> RealGradient FE<0,LAGRANGE_VEC>::shape(const Elem * elem, const Order order,
                                                   const unsigned int i, const Point & p)
{
  Real value = FE<0,LAGRANGE>::shape( elem->type(), static_cast<Order>(order + elem->p_level()), i, p);
  return libMesh::RealGradient( value );
}
template <> RealGradient FE<0,LAGRANGE_VEC>::shape_deriv(const Elem * elem, const Order order,
                                                         const unsigned int i, const unsigned int j,
                                                         const Point & p)
{
  Real value = FE<0,LAGRANGE>::shape_deriv( elem->type(), static_cast<Order>(order + elem->p_level()), i, j, p);
  return libMesh::RealGradient( value );
}
template <> RealGradient FE<0,LAGRANGE_VEC>::shape_second_deriv(const Elem * elem, const Order order,
                                                                const unsigned int i, const unsigned int j,
                                                                const Point & p)
{
  Real value = FE<0,LAGRANGE>::shape_second_deriv( elem->type(), static_cast<Order>(order + elem->p_level()), i, j, p);
  return libMesh::RealGradient( value );
}

// 1-D
template <> RealGradient FE<1,LAGRANGE_VEC>::shape(const Elem * elem, const Order order,
                                                   const unsigned int i, const Point & p)
{
  Real value = FE<1,LAGRANGE>::shape( elem->type(), static_cast<Order>(order + elem->p_level()), i, p);
  return libMesh::RealGradient( value );
}
template <> RealGradient FE<1,LAGRANGE_VEC>::shape_deriv(const Elem * elem, const Order order,
                                                         const unsigned int i, const unsigned int j,
                                                         const Point & p)
{
  Real value = FE<1,LAGRANGE>::shape_deriv( elem->type(), static_cast<Order>(order + elem->p_level()), i, j, p);
  return libMesh::RealGradient( value );
}
template <> RealGradient FE<1,LAGRANGE_VEC>::shape_second_deriv(const Elem * elem, const Order order,
                                                                const unsigned int i, const unsigned int j,
                                                                const Point & p)
{
  Real value = FE<1,LAGRANGE>::shape_second_deriv( elem->type(), static_cast<Order>(order + elem->p_level()), i, j, p);
  return libMesh::RealGradient( value );
}

// 2-D
template <> RealGradient FE<2,LAGRANGE_VEC>::shape(const Elem * elem, const Order order,
                                                   const unsigned int i, const Point & p)
{
  Real value = FE<2,LAGRANGE>::shape( elem->type(), static_cast<Order>(order + elem->p_level()), i/2, p );

  switch( i%2 )
    {
    case 0:
      return libMesh::RealGradient( value );

    case 1:
      return libMesh::RealGradient( Real(0), value );

    default:
      libmesh_error_msg("i%2 must be either 0 or 1!");
    }

  //dummy
  return libMesh::RealGradient();
}
template <> RealGradient FE<2,LAGRANGE_VEC>::shape_deriv(const Elem * elem, const Order order,
                                                         const unsigned int i, const unsigned int j,
                                                         const Point & p)
{
  Real value = FE<2,LAGRANGE>::shape_deriv( elem->type(), static_cast<Order>(order + elem->p_level()), i/2, j, p );

  switch( i%2 )
    {
    case 0:
      return libMesh::RealGradient( value );

    case 1:
      return libMesh::RealGradient( Real(0), value );

    default:
      libmesh_error_msg("i%2 must be either 0 or 1!");
    }

  //dummy
  return libMesh::RealGradient();
}
template <> RealGradient FE<2,LAGRANGE_VEC>::shape_second_deriv(const Elem * elem, const Order order,
                                                                const unsigned int i, const unsigned int j,
                                                                const Point & p)
{
  Real value = FE<2,LAGRANGE>::shape_second_deriv( elem->type(), static_cast<Order>(order + elem->p_level()), i/2, j, p );

  switch( i%2 )
    {
    case 0:
      return libMesh::RealGradient( value );

    case 1:
      return libMesh::RealGradient( Real(0), value );

    default:
      libmesh_error_msg("i%2 must be either 0 or 1!");
    }

  //dummy
  return libMesh::RealGradient();
}

// 3-D
template <> RealGradient FE<3,LAGRANGE_VEC>::shape(const Elem * elem, const Order order,
                                                   const unsigned int i, const Point & p)
{
  Real value = FE<3,LAGRANGE>::shape( elem->type(), static_cast<Order>(order + elem->p_level()), i/3, p );

  switch( i%3 )
    {
    case 0:
      return libMesh::RealGradient( value );

    case 1:
      return libMesh::RealGradient( Real(0), value );

    case 2:
      return libMesh::RealGradient( Real(0), Real(0), value );

    default:
      libmesh_error_msg("i%3 must be 0, 1, or 2!");
    }

  //dummy
  return libMesh::RealGradient();
}
template <> RealGradient FE<3,LAGRANGE_VEC>::shape_deriv(const Elem * elem, const Order order,
                                                         const unsigned int i, const unsigned int j,
                                                         const Point & p)
{
  Real value = FE<3,LAGRANGE>::shape_deriv( elem->type(), static_cast<Order>(order + elem->p_level()), i/3, j, p );

  switch( i%3 )
    {
    case 0:
      return libMesh::RealGradient( value );

    case 1:
      return libMesh::RealGradient( Real(0), value );

    case 2:
      return libMesh::RealGradient( Real(0), Real(0), value );

    default:
      libmesh_error_msg("i%3 must be 0, 1, or 2!");
    }

  //dummy
  return libMesh::RealGradient();
}
template <> RealGradient FE<3,LAGRANGE_VEC>::shape_second_deriv(const Elem * elem, const Order order,
                                                                const unsigned int i, const unsigned int j,
                                                                const Point & p)
{
  Real value = FE<3,LAGRANGE>::shape_second_deriv( elem->type(), static_cast<Order>(order + elem->p_level()), i/3, j, p );

  switch( i%3 )
    {
    case 0:
      return libMesh::RealGradient( value );

    case 1:
      return libMesh::RealGradient( Real(0), value );

    case 2:
      return libMesh::RealGradient( Real(0), Real(0), value );

    default:
      libmesh_error_msg("i%3 must be 0, 1, or 2!");
    }

  //dummy
  return libMesh::RealGradient();
}

// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this function.
// This could be macro-ified.
template <> unsigned int FE<0,LAGRANGE_VEC>::n_dofs(const ElemType t, const Order o) { return FE<0,LAGRANGE>::n_dofs(t,o); }
template <> unsigned int FE<1,LAGRANGE_VEC>::n_dofs(const ElemType t, const Order o) { return FE<1,LAGRANGE>::n_dofs(t,o); }
template <> unsigned int FE<2,LAGRANGE_VEC>::n_dofs(const ElemType t, const Order o) { return 2*FE<2,LAGRANGE>::n_dofs(t,o); }
template <> unsigned int FE<3,LAGRANGE_VEC>::n_dofs(const ElemType t, const Order o) { return 3*FE<3,LAGRANGE>::n_dofs(t,o); }


// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this function.
template <> unsigned int FE<0,LAGRANGE_VEC>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return FE<0,LAGRANGE>::n_dofs_at_node(t,o,n); }
template <> unsigned int FE<1,LAGRANGE_VEC>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return FE<1,LAGRANGE>::n_dofs_at_node(t,o,n); }
template <> unsigned int FE<2,LAGRANGE_VEC>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return 2*FE<2,LAGRANGE>::n_dofs_at_node(t,o,n); }
template <> unsigned int FE<3,LAGRANGE_VEC>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return 3*FE<2,LAGRANGE>::n_dofs_at_node(t,o,n); }


// Lagrange elements have no dofs per element
// (just at the nodes)
template <> unsigned int FE<0,LAGRANGE_VEC>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <> unsigned int FE<1,LAGRANGE_VEC>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <> unsigned int FE<2,LAGRANGE_VEC>::n_dofs_per_elem(const ElemType, const Order) { return 0; }
template <> unsigned int FE<3,LAGRANGE_VEC>::n_dofs_per_elem(const ElemType, const Order) { return 0; }

// Lagrange FEMs are always C^0 continuous
template <> FEContinuity FE<0,LAGRANGE_VEC>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<1,LAGRANGE_VEC>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<2,LAGRANGE_VEC>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<3,LAGRANGE_VEC>::get_continuity() const { return C_ZERO; }

// Lagrange FEMs are not hierarchic
template <> bool FE<0,LAGRANGE_VEC>::is_hierarchic() const { return false; }
template <> bool FE<1,LAGRANGE_VEC>::is_hierarchic() const { return false; }
template <> bool FE<2,LAGRANGE_VEC>::is_hierarchic() const { return false; }
template <> bool FE<3,LAGRANGE_VEC>::is_hierarchic() const { return false; }

// Lagrange FEM shapes do not need reinit (is this always true?)
template <> bool FE<0,LAGRANGE_VEC>::shapes_need_reinit() const { return false; }
template <> bool FE<1,LAGRANGE_VEC>::shapes_need_reinit() const { return false; }
template <> bool FE<2,LAGRANGE_VEC>::shapes_need_reinit() const { return false; }
template <> bool FE<3,LAGRANGE_VEC>::shapes_need_reinit() const { return false; }

// Methods for computing Lagrange constraints.  Note: we pass the
// dimension as the last argument to the anonymous helper function.
// Also note: we only need instantiations of this function for
// Dim==2 and 3.
#ifdef LIBMESH_ENABLE_AMR
template <>
void FE<2,LAGRANGE_VEC>::compute_constraints (DofConstraints & constraints,
                                              DofMap & dof_map,
                                              const unsigned int variable_number,
                                              const Elem * elem)
{ //libmesh_not_implemented();
  FEVectorBase::compute_proj_constraints(constraints, dof_map, variable_number, elem);
}

template <>
void FE<3,LAGRANGE_VEC>::compute_constraints (DofConstraints & constraints,
                                              DofMap & dof_map,
                                              const unsigned int variable_number,
                                              const Elem * elem)
{ //libmesh_not_implemented();
  FEVectorBase::compute_proj_constraints(constraints, dof_map, variable_number, elem);
}
#endif // LIBMESH_ENABLE_AMR

} // namespace libMesh
