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

#ifndef LIBMESH_FE_CLOUGH_IMPL_SHAPE_2D_H
#define LIBMESH_FE_CLOUGH_IMPL_SHAPE_2D_H

// C++ includes

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"


// Anonymous namespace for persistent variables.
// This allows us to cache the global-to-local mapping transformation
// FIXME: This should also screw up multithreading royally
namespace
{
using namespace libMesh;

// Keep track of which element was most recently used to generate
// cached data
static dof_id_type old_elem_id = DofObject::invalid_id;
static const Elem * old_elem_ptr = nullptr;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <typename RealType>
RealType clough_raw_shape_second_deriv(const unsigned int basis_num,
                                       const unsigned int deriv_type,
                                       const PointTempl<RealType> & p);
#endif

template <typename RealType>
RealType clough_raw_shape_deriv(const unsigned int basis_num,
                            const unsigned int deriv_type,
                            const PointTempl<RealType> & p);
template <typename RealType>
RealType clough_raw_shape(const unsigned int basis_num,
                      const PointTempl<RealType> & p);
template <typename RealType>
unsigned char subtriangle_lookup(const PointTempl<RealType> & p);


// Compute the static coefficients for an element
template <typename RealType>
void clough_compute_coefs(const ElemTempl<RealType> * elem,
                          RealType & d1d2n,
                          RealType & d1d3n,
                          RealType & d2d3n,
                          RealType & d2d1n,
                          RealType & d3d1n,
                          RealType & d3d2n,
                          RealType & d1xd1x,
                          RealType & d1xd1y,
                          RealType & d1xd2n,
                          RealType & d1xd3n,
                          RealType & d1yd1x,
                          RealType & d1yd1y,
                          RealType & d1yd2n,
                          RealType & d1yd3n,
                          RealType & d2xd2x,
                          RealType & d2xd2y,
                          RealType & d2xd3n,
                          RealType & d2xd1n,
                          RealType & d2yd2x,
                          RealType & d2yd2y,
                          RealType & d2yd3n,
                          RealType & d2yd1n,
                          RealType & d3xd3x,
                          RealType & d3xd3y,
                          RealType & d3xd1n,
                          RealType & d3xd2n,
                          RealType & d3yd3x,
                          RealType & d3yd3y,
                          RealType & d3yd1n,
                          RealType & d3yd2n,
                          RealType & d1nd1n,
                          RealType & d2nd2n,
                          RealType & d3nd3n,
                          RealType & N01x,
                          RealType & N01y,
                          RealType & N10x,
                          RealType & N10y,
                          RealType & N02x,
                          RealType & N02y,
                          RealType & N20x,
                          RealType & N20y,
                          RealType & N21x,
                          RealType & N21y,
                          RealType & N12x,
                          RealType & N12y)
{
  typedef PointTempl<RealType> Point;

  // Using static globals for old_elem_id, etc. will fail
  // horribly with more than one thread.
  libmesh_assert_equal_to (libMesh::n_threads(), 1);

  // Coefficients are cached from old elements; we rely on that cache
  // except in dbg mode
#ifndef DEBUG
  if (elem->id() == old_elem_id &&
      elem == old_elem_ptr)
    return;
#endif

  const FEFamily mapping_family = FEMapTempl<RealType>::map_fe_type(*elem);
  const Order mapping_order        (elem->default_order());
  const ElemType mapping_elem_type (elem->type());

  const FEType map_fe_type(mapping_order, mapping_family);

  const int n_mapping_shape_functions =
    FEInterface::n_shape_functions(2, map_fe_type, mapping_elem_type);

  // Degrees of freedom are at vertices and edge midpoints
  std::vector<Point> dofpt;
  dofpt.push_back(Point(0,0));
  dofpt.push_back(Point(1,0));
  dofpt.push_back(Point(0,1));
  dofpt.push_back(Point(1/2.,1/2.));
  dofpt.push_back(Point(0,1/2.));
  dofpt.push_back(Point(1/2.,0));

  // Mapping functions - first derivatives at each dofpt
  std::vector<RealType> dxdxi(6), dxdeta(6), dydxi(6), dydeta(6);
  std::vector<RealType> dxidx(6), detadx(6), dxidy(6), detady(6);

  FEInterface::shape_deriv_ptr<RealType> shape_deriv_ptr =
    FEInterface::shape_deriv_function<RealType>(2, map_fe_type);

  for (int p = 0; p != 6; ++p)
    {
      //      libMesh::err << p << ' ' << dofpt[p];
      for (int i = 0; i != n_mapping_shape_functions; ++i)
        {
          const auto ddxi = shape_deriv_ptr
            (elem, mapping_order, i, 0, dofpt[p], false);
          const auto ddeta = shape_deriv_ptr
            (elem, mapping_order, i, 1, dofpt[p], false);

          //      libMesh::err << ddxi << ' ';
          //      libMesh::err << ddeta << std::endl;

          dxdxi[p] += elem->point(i)(0) * ddxi;
          dydxi[p] += elem->point(i)(1) * ddxi;
          dxdeta[p] += elem->point(i)(0) * ddeta;
          dydeta[p] += elem->point(i)(1) * ddeta;
        }

      //      for (int i = 0; i != 12; ++i)
      //          libMesh::err << i << ' ' << clough_raw_shape(i, dofpt[p]) << std::endl;

      //      libMesh::err << elem->point(p)(0) << ' ';
      //      libMesh::err << elem->point(p)(1) << ' ';
      //      libMesh::err << dxdxi[p] << ' ';
      //      libMesh::err << dydxi[p] << ' ';
      //      libMesh::err << dxdeta[p] << ' ';
      //      libMesh::err << dydeta[p] << std::endl << std::endl;

      const auto inv_jac = 1. / (dxdxi[p]*dydeta[p] -
                                 dxdeta[p]*dydxi[p]);
      dxidx[p] = dydeta[p] * inv_jac;
      dxidy[p] = - dxdeta[p] * inv_jac;
      detadx[p] = - dydxi[p] * inv_jac;
      detady[p] = dxdxi[p] * inv_jac;
    }

  // Calculate midpoint normal vectors
  auto N1x = dydeta[3] - dydxi[3];
  auto N1y = dxdxi[3] - dxdeta[3];
  auto Nlength = std::sqrt(static_cast<RealType>(N1x*N1x + N1y*N1y));
  N1x /= Nlength; N1y /= Nlength;

  auto N2x = - dydeta[4];
  auto N2y = dxdeta[4];
  Nlength = std::sqrt(static_cast<RealType>(N2x*N2x + N2y*N2y));
  N2x /= Nlength; N2y /= Nlength;

  auto N3x = dydxi[5];
  auto N3y = - dxdxi[5];
  Nlength = std::sqrt(static_cast<RealType>(N3x*N3x + N3y*N3y));
  N3x /= Nlength; N3y /= Nlength;

  // Calculate corner normal vectors (used for reduced element)
  N01x = dydxi[0];
  N01y = - dxdxi[0];
  Nlength = std::sqrt(static_cast<RealType>(N01x*N01x + N01y*N01y));
  N01x /= Nlength; N01y /= Nlength;

  N10x = dydxi[1];
  N10y = - dxdxi[1];
  Nlength = std::sqrt(static_cast<RealType>(N10x*N10x + N10y*N10y));
  N10x /= Nlength; N10y /= Nlength;

  N02x = - dydeta[0];
  N02y = dxdeta[0];
  Nlength = std::sqrt(static_cast<RealType>(N02x*N02x + N02y*N02y));
  N02x /= Nlength; N02y /= Nlength;

  N20x = - dydeta[2];
  N20y = dxdeta[2];
  Nlength = std::sqrt(static_cast<RealType>(N20x*N20x + N20y*N20y));
  N20x /= Nlength; N20y /= Nlength;

  N12x = dydeta[1] - dydxi[1];
  N12y = dxdxi[1] - dxdeta[1];
  Nlength = std::sqrt(static_cast<RealType>(N12x*N12x + N12y*N12y));
  N12x /= Nlength; N12y /= Nlength;

  N21x = dydeta[1] - dydxi[1];
  N21y = dxdxi[1] - dxdeta[1];
  Nlength = std::sqrt(static_cast<RealType>(N21x*N21x + N21y*N21y));
  N21x /= Nlength; N21y /= Nlength;

  //  for (int i=0; i != 6; ++i) {
  //    libMesh::err << elem->node_id(i) << ' ';
  //  }
  //  libMesh::err << std::endl;

  //  for (int i=0; i != 6; ++i) {
  //    libMesh::err << elem->point(i)(0) << ' ';
  //    libMesh::err << elem->point(i)(1) << ' ';
  //  }
  //  libMesh::err << std::endl;


  // give normal vectors a globally unique orientation

  if (elem->point(2) < elem->point(1))
    {
      //      libMesh::err << "Flipping nodes " << elem->node_id(2);
      //      libMesh::err << " and " << elem->node_id(1);
      //      libMesh::err << " around node " << elem->node_id(4);
      //      libMesh::err << std::endl;
      N1x = -N1x; N1y = -N1y;
      N12x = -N12x; N12y = -N12y;
      N21x = -N21x; N21y = -N21y;
    }
  else
    {
      //      libMesh::err << "Not flipping nodes " << elem->node_id(2);
      //      libMesh::err << " and " << elem->node_id(1);
      //      libMesh::err << " around node " << elem->node_id(4);
      //      libMesh::err << std::endl;
    }
  if (elem->point(0) < elem->point(2))
    {
      //      libMesh::err << "Flipping nodes " << elem->node_id(0);
      //      libMesh::err << " and " << elem->node_id(2);
      //      libMesh::err << " around node " << elem->node_id(5);
      //      libMesh::err << std::endl;
      //      libMesh::err << N2x << ' ' << N2y << std::endl;
      N2x = -N2x; N2y = -N2y;
      N02x = -N02x; N02y = -N02y;
      N20x = -N20x; N20y = -N20y;
      //      libMesh::err << N2x << ' ' << N2y << std::endl;
    }
  else
    {
      //      libMesh::err << "Not flipping nodes " << elem->node_id(0);
      //      libMesh::err << " and " << elem->node_id(2);
      //      libMesh::err << " around node " << elem->node_id(5);
      //      libMesh::err << std::endl;
    }
  if (elem->point(1) < elem->point(0))
    {
      //      libMesh::err << "Flipping nodes " << elem->node_id(1);
      //      libMesh::err << " and " << elem->node_id(0);
      //      libMesh::err << " around node " << elem->node_id(3);
      //      libMesh::err << std::endl;
      N3x = -N3x;
      N3y = -N3y;
      N01x = -N01x; N01y = -N01y;
      N10x = -N10x; N10y = -N10y;
    }
  else
    {
      //      libMesh::err << "Not flipping nodes " << elem->node_id(1);
      //      libMesh::err << " and " << elem->node_id(0);
      //      libMesh::err << " around node " << elem->node_id(3);
      //      libMesh::err << std::endl;
    }

  //  libMesh::err << N2x << ' ' << N2y << std::endl;

  // Cache basis function gradients
  // FIXME: the raw_shape calls shouldn't be done on every element!
  // FIXME: I should probably be looping, too...
  // Gradient naming: d(1)d(2n)d(xi) is the xi component of the
  // gradient of the
  // local basis function corresponding to value 1 at the node
  // corresponding to normal vector 2

  auto d1d2ndxi   = clough_raw_shape_deriv(0, 0, dofpt[4]);
  auto d1d2ndeta  = clough_raw_shape_deriv(0, 1, dofpt[4]);
  auto d1d2ndx = d1d2ndxi * dxidx[4] + d1d2ndeta * detadx[4];
  auto d1d2ndy = d1d2ndxi * dxidy[4] + d1d2ndeta * detady[4];
  auto d1d3ndxi   = clough_raw_shape_deriv(0, 0, dofpt[5]);
  auto d1d3ndeta  = clough_raw_shape_deriv(0, 1, dofpt[5]);
  auto d1d3ndx = d1d3ndxi * dxidx[5] + d1d3ndeta * detadx[5];
  auto d1d3ndy = d1d3ndxi * dxidy[5] + d1d3ndeta * detady[5];
  auto d2d3ndxi   = clough_raw_shape_deriv(1, 0, dofpt[5]);
  auto d2d3ndeta  = clough_raw_shape_deriv(1, 1, dofpt[5]);
  auto d2d3ndx = d2d3ndxi * dxidx[5] + d2d3ndeta * detadx[5];
  auto d2d3ndy = d2d3ndxi * dxidy[5] + d2d3ndeta * detady[5];
  auto d2d1ndxi   = clough_raw_shape_deriv(1, 0, dofpt[3]);
  auto d2d1ndeta  = clough_raw_shape_deriv(1, 1, dofpt[3]);
  auto d2d1ndx = d2d1ndxi * dxidx[3] + d2d1ndeta * detadx[3];
  auto d2d1ndy = d2d1ndxi * dxidy[3] + d2d1ndeta * detady[3];
  auto d3d1ndxi   = clough_raw_shape_deriv(2, 0, dofpt[3]);
  auto d3d1ndeta  = clough_raw_shape_deriv(2, 1, dofpt[3]);
  auto d3d1ndx = d3d1ndxi * dxidx[3] + d3d1ndeta * detadx[3];
  auto d3d1ndy = d3d1ndxi * dxidy[3] + d3d1ndeta * detady[3];
  auto d3d2ndxi   = clough_raw_shape_deriv(2, 0, dofpt[4]);
  auto d3d2ndeta  = clough_raw_shape_deriv(2, 1, dofpt[4]);
  auto d3d2ndx = d3d2ndxi * dxidx[4] + d3d2ndeta * detadx[4];
  auto d3d2ndy = d3d2ndxi * dxidy[4] + d3d2ndeta * detady[4];
  auto d1xd2ndxi  = clough_raw_shape_deriv(3, 0, dofpt[4]);
  auto d1xd2ndeta = clough_raw_shape_deriv(3, 1, dofpt[4]);
  auto d1xd2ndx = d1xd2ndxi * dxidx[4] + d1xd2ndeta * detadx[4];
  auto d1xd2ndy = d1xd2ndxi * dxidy[4] + d1xd2ndeta * detady[4];
  auto d1xd3ndxi  = clough_raw_shape_deriv(3, 0, dofpt[5]);
  auto d1xd3ndeta = clough_raw_shape_deriv(3, 1, dofpt[5]);
  auto d1xd3ndx = d1xd3ndxi * dxidx[5] + d1xd3ndeta * detadx[5];
  auto d1xd3ndy = d1xd3ndxi * dxidy[5] + d1xd3ndeta * detady[5];
  auto d1yd2ndxi  = clough_raw_shape_deriv(4, 0, dofpt[4]);
  auto d1yd2ndeta = clough_raw_shape_deriv(4, 1, dofpt[4]);
  auto d1yd2ndx = d1yd2ndxi * dxidx[4] + d1yd2ndeta * detadx[4];
  auto d1yd2ndy = d1yd2ndxi * dxidy[4] + d1yd2ndeta * detady[4];
  auto d1yd3ndxi  = clough_raw_shape_deriv(4, 0, dofpt[5]);
  auto d1yd3ndeta = clough_raw_shape_deriv(4, 1, dofpt[5]);
  auto d1yd3ndx = d1yd3ndxi * dxidx[5] + d1yd3ndeta * detadx[5];
  auto d1yd3ndy = d1yd3ndxi * dxidy[5] + d1yd3ndeta * detady[5];
  auto d2xd3ndxi  = clough_raw_shape_deriv(5, 0, dofpt[5]);
  auto d2xd3ndeta = clough_raw_shape_deriv(5, 1, dofpt[5]);
  auto d2xd3ndx = d2xd3ndxi * dxidx[5] + d2xd3ndeta * detadx[5];
  auto d2xd3ndy = d2xd3ndxi * dxidy[5] + d2xd3ndeta * detady[5];
  auto d2xd1ndxi  = clough_raw_shape_deriv(5, 0, dofpt[3]);
  auto d2xd1ndeta = clough_raw_shape_deriv(5, 1, dofpt[3]);
  auto d2xd1ndx = d2xd1ndxi * dxidx[3] + d2xd1ndeta * detadx[3];
  auto d2xd1ndy = d2xd1ndxi * dxidy[3] + d2xd1ndeta * detady[3];
  auto d2yd3ndxi  = clough_raw_shape_deriv(6, 0, dofpt[5]);
  auto d2yd3ndeta = clough_raw_shape_deriv(6, 1, dofpt[5]);
  auto d2yd3ndx = d2yd3ndxi * dxidx[5] + d2yd3ndeta * detadx[5];
  auto d2yd3ndy = d2yd3ndxi * dxidy[5] + d2yd3ndeta * detady[5];
  auto d2yd1ndxi  = clough_raw_shape_deriv(6, 0, dofpt[3]);
  auto d2yd1ndeta = clough_raw_shape_deriv(6, 1, dofpt[3]);
  auto d2yd1ndx = d2yd1ndxi * dxidx[3] + d2yd1ndeta * detadx[3];
  auto d2yd1ndy = d2yd1ndxi * dxidy[3] + d2yd1ndeta * detady[3];
  auto d3xd1ndxi  = clough_raw_shape_deriv(7, 0, dofpt[3]);
  auto d3xd1ndeta = clough_raw_shape_deriv(7, 1, dofpt[3]);
  auto d3xd1ndx = d3xd1ndxi * dxidx[3] + d3xd1ndeta * detadx[3];
  auto d3xd1ndy = d3xd1ndxi * dxidy[3] + d3xd1ndeta * detady[3];
  auto d3xd2ndxi  = clough_raw_shape_deriv(7, 0, dofpt[4]);
  auto d3xd2ndeta = clough_raw_shape_deriv(7, 1, dofpt[4]);
  auto d3xd2ndx = d3xd2ndxi * dxidx[4] + d3xd2ndeta * detadx[4];
  auto d3xd2ndy = d3xd2ndxi * dxidy[4] + d3xd2ndeta * detady[4];
  auto d3yd1ndxi  = clough_raw_shape_deriv(8, 0, dofpt[3]);
  auto d3yd1ndeta = clough_raw_shape_deriv(8, 1, dofpt[3]);
  auto d3yd1ndx = d3yd1ndxi * dxidx[3] + d3yd1ndeta * detadx[3];
  auto d3yd1ndy = d3yd1ndxi * dxidy[3] + d3yd1ndeta * detady[3];
  auto d3yd2ndxi  = clough_raw_shape_deriv(8, 0, dofpt[4]);
  auto d3yd2ndeta = clough_raw_shape_deriv(8, 1, dofpt[4]);
  auto d3yd2ndx = d3yd2ndxi * dxidx[4] + d3yd2ndeta * detadx[4];
  auto d3yd2ndy = d3yd2ndxi * dxidy[4] + d3yd2ndeta * detady[4];
  auto d1nd1ndxi  = clough_raw_shape_deriv(9, 0, dofpt[3]);
  auto d1nd1ndeta = clough_raw_shape_deriv(9, 1, dofpt[3]);
  auto d1nd1ndx = d1nd1ndxi * dxidx[3] + d1nd1ndeta * detadx[3];
  auto d1nd1ndy = d1nd1ndxi * dxidy[3] + d1nd1ndeta * detady[3];
  auto d2nd2ndxi  = clough_raw_shape_deriv(10, 0, dofpt[4]);
  auto d2nd2ndeta = clough_raw_shape_deriv(10, 1, dofpt[4]);
  auto d2nd2ndx = d2nd2ndxi * dxidx[4] + d2nd2ndeta * detadx[4];
  auto d2nd2ndy = d2nd2ndxi * dxidy[4] + d2nd2ndeta * detady[4];
  auto d3nd3ndxi  = clough_raw_shape_deriv(11, 0, dofpt[5]);
  auto d3nd3ndeta = clough_raw_shape_deriv(11, 1, dofpt[5]);
  auto d3nd3ndx = d3nd3ndxi * dxidx[3] + d3nd3ndeta * detadx[3];
  auto d3nd3ndy = d3nd3ndxi * dxidy[3] + d3nd3ndeta * detady[3];

  auto d1xd1dxi   = clough_raw_shape_deriv(3, 0, dofpt[0]);
  auto d1xd1deta  = clough_raw_shape_deriv(3, 1, dofpt[0]);
  auto d1xd1dx    = d1xd1dxi * dxidx[0] + d1xd1deta * detadx[0];
  auto d1xd1dy    = d1xd1dxi * dxidy[0] + d1xd1deta * detady[0];
  auto d1yd1dxi   = clough_raw_shape_deriv(4, 0, dofpt[0]);
  auto d1yd1deta  = clough_raw_shape_deriv(4, 1, dofpt[0]);
  auto d1yd1dx    = d1yd1dxi * dxidx[0] + d1yd1deta * detadx[0];
  auto d1yd1dy    = d1yd1dxi * dxidy[0] + d1yd1deta * detady[0];
  auto d2xd2dxi   = clough_raw_shape_deriv(5, 0, dofpt[1]);
  auto d2xd2deta  = clough_raw_shape_deriv(5, 1, dofpt[1]);
  auto d2xd2dx    = d2xd2dxi * dxidx[1] + d2xd2deta * detadx[1];
  auto d2xd2dy    = d2xd2dxi * dxidy[1] + d2xd2deta * detady[1];

  //  libMesh::err << dofpt[4](0) << ' ';
  //  libMesh::err << dofpt[4](1) << ' ';
  //  libMesh::err << (int)subtriangle_lookup(dofpt[5]) << ' ';
  //  libMesh::err << dxdxi[4] << ' ';
  //  libMesh::err << dxdeta[4] << ' ';
  //  libMesh::err << dydxi[4] << ' ';
  //  libMesh::err << dydeta[4] << ' ';
  //  libMesh::err << dxidx[4] << ' ';
  //  libMesh::err << dxidy[4] << ' ';
  //  libMesh::err << detadx[4] << ' ';
  //  libMesh::err << detady[4] << ' ';
  //  libMesh::err << N1x << ' ';
  //  libMesh::err << N1y << ' ';
  //  libMesh::err << d2yd1ndxi << ' ';
  //  libMesh::err << d2yd1ndeta << ' ';
  //  libMesh::err << d2yd1ndx << ' ';
  //  libMesh::err << d2yd1ndy << std::endl;

  auto d2yd2dxi   = clough_raw_shape_deriv(6, 0, dofpt[1]);
  auto d2yd2deta  = clough_raw_shape_deriv(6, 1, dofpt[1]);
  auto d2yd2dx    = d2yd2dxi * dxidx[1] + d2yd2deta * detadx[1];
  auto d2yd2dy    = d2yd2dxi * dxidy[1] + d2yd2deta * detady[1];
  auto d3xd3dxi   = clough_raw_shape_deriv(7, 0, dofpt[2]);
  auto d3xd3deta  = clough_raw_shape_deriv(7, 1, dofpt[2]);
  auto d3xd3dx    = d3xd3dxi * dxidx[1] + d3xd3deta * detadx[1];
  auto d3xd3dy    = d3xd3dxi * dxidy[1] + d3xd3deta * detady[1];
  auto d3yd3dxi   = clough_raw_shape_deriv(8, 0, dofpt[2]);
  auto d3yd3deta  = clough_raw_shape_deriv(8, 1, dofpt[2]);
  auto d3yd3dx    = d3yd3dxi * dxidx[1] + d3yd3deta * detadx[1];
  auto d3yd3dy    = d3yd3dxi * dxidy[1] + d3yd3deta * detady[1];

  // Calculate normal derivatives

  auto d1nd1ndn = d1nd1ndx * N1x + d1nd1ndy * N1y;
  auto d1xd2ndn = d1xd2ndx * N2x + d1xd2ndy * N2y;
  auto d1xd3ndn = d1xd3ndx * N3x + d1xd3ndy * N3y;
  auto d1yd2ndn = d1yd2ndx * N2x + d1yd2ndy * N2y;
  auto d1yd3ndn = d1yd3ndx * N3x + d1yd3ndy * N3y;

  auto d2nd2ndn = d2nd2ndx * N2x + d2nd2ndy * N2y;
  auto d2xd3ndn = d2xd3ndx * N3x + d2xd3ndy * N3y;
  auto d2xd1ndn = d2xd1ndx * N1x + d2xd1ndy * N1y;
  auto d2yd3ndn = d2yd3ndx * N3x + d2yd3ndy * N3y;
  auto d2yd1ndn = d2yd1ndx * N1x + d2yd1ndy * N1y;

  auto d3nd3ndn = d3nd3ndx * N3x + d3nd3ndy * N3y;
  auto d3xd1ndn = d3xd1ndx * N1x + d3xd1ndy * N1y;
  auto d3xd2ndn = d3xd2ndx * N2x + d3xd2ndy * N2y;
  auto d3yd1ndn = d3yd1ndx * N1x + d3yd1ndy * N1y;
  auto d3yd2ndn = d3yd2ndx * N2x + d3yd2ndy * N2y;

  // Calculate midpoint scaling factors

  d1nd1n = 1. / d1nd1ndn;
  d2nd2n = 1. / d2nd2ndn;
  d3nd3n = 1. / d3nd3ndn;

  // Calculate midpoint derivative adjustments to nodal value
  // interpolant functions

#ifdef DEBUG
  // The cached factors should equal our calculations if we're still
  // operating on the cached element
  if (elem->id() == old_elem_id &&
      elem == old_elem_ptr)
    {
      libmesh_assert_equal_to(d1d2n, -(d1d2ndx * N2x + d1d2ndy * N2y) / d2nd2ndn);
      libmesh_assert_equal_to(d1d3n, -(d1d3ndx * N3x + d1d3ndy * N3y) / d3nd3ndn);
      libmesh_assert_equal_to(d2d3n, -(d2d3ndx * N3x + d2d3ndy * N3y) / d3nd3ndn);
      libmesh_assert_equal_to(d2d1n, -(d2d1ndx * N1x + d2d1ndy * N1y) / d1nd1ndn);
      libmesh_assert_equal_to(d3d1n, -(d3d1ndx * N1x + d3d1ndy * N1y) / d1nd1ndn);
      libmesh_assert_equal_to(d3d2n, -(d3d2ndx * N2x + d3d2ndy * N2y) / d2nd2ndn);
      libmesh_assert_equal_to(d1xd1x, 1. / (d1xd1dx - d1xd1dy * d1yd1dx / d1yd1dy));
      libmesh_assert_equal_to(d1xd1y, 1. / (d1yd1dx - d1xd1dx * d1yd1dy / d1xd1dy));
      libmesh_assert_equal_to(d1yd1y, 1. / (d1yd1dy - d1yd1dx * d1xd1dy / d1xd1dx));
      libmesh_assert_equal_to(d1yd1x, 1. / (d1xd1dy - d1yd1dy * d1xd1dx / d1yd1dx));
      libmesh_assert_equal_to(d2xd2x, 1. / (d2xd2dx - d2xd2dy * d2yd2dx / d2yd2dy));
      libmesh_assert_equal_to(d2xd2y, 1. / (d2yd2dx - d2xd2dx * d2yd2dy / d2xd2dy));
      libmesh_assert_equal_to(d2yd2y, 1. / (d2yd2dy - d2yd2dx * d2xd2dy / d2xd2dx));
      libmesh_assert_equal_to(d2yd2x, 1. / (d2xd2dy - d2yd2dy * d2xd2dx / d2yd2dx));
      libmesh_assert_equal_to(d3xd3x, 1. / (d3xd3dx - d3xd3dy * d3yd3dx / d3yd3dy));
      libmesh_assert_equal_to(d3xd3y, 1. / (d3yd3dx - d3xd3dx * d3yd3dy / d3xd3dy));
      libmesh_assert_equal_to(d3yd3y, 1. / (d3yd3dy - d3yd3dx * d3xd3dy / d3xd3dx));
      libmesh_assert_equal_to(d3yd3x, 1. / (d3xd3dy - d3yd3dy * d3xd3dx / d3yd3dx));
      libmesh_assert_equal_to(d1xd2n, -(d1xd1x * d1xd2ndn + d1xd1y * d1yd2ndn) / d2nd2ndn);
      libmesh_assert_equal_to(d1yd2n, -(d1yd1y * d1yd2ndn + d1yd1x * d1xd2ndn) / d2nd2ndn);
      libmesh_assert_equal_to(d1xd3n, -(d1xd1x * d1xd3ndn + d1xd1y * d1yd3ndn) / d3nd3ndn);
      libmesh_assert_equal_to(d1yd3n, -(d1yd1y * d1yd3ndn + d1yd1x * d1xd3ndn) / d3nd3ndn);
      libmesh_assert_equal_to(d2xd3n, -(d2xd2x * d2xd3ndn + d2xd2y * d2yd3ndn) / d3nd3ndn);
      libmesh_assert_equal_to(d2yd3n, -(d2yd2y * d2yd3ndn + d2yd2x * d2xd3ndn) / d3nd3ndn);
      libmesh_assert_equal_to(d2xd1n, -(d2xd2x * d2xd1ndn + d2xd2y * d2yd1ndn) / d1nd1ndn);
      libmesh_assert_equal_to(d2yd1n, -(d2yd2y * d2yd1ndn + d2yd2x * d2xd1ndn) / d1nd1ndn);
      libmesh_assert_equal_to(d3xd1n, -(d3xd3x * d3xd1ndn + d3xd3y * d3yd1ndn) / d1nd1ndn);
      libmesh_assert_equal_to(d3yd1n, -(d3yd3y * d3yd1ndn + d3yd3x * d3xd1ndn) / d1nd1ndn);
      libmesh_assert_equal_to(d3xd2n, -(d3xd3x * d3xd2ndn + d3xd3y * d3yd2ndn) / d2nd2ndn);
      libmesh_assert_equal_to(d3yd2n, -(d3yd3y * d3yd2ndn + d3yd3x * d3xd2ndn) / d2nd2ndn);
    }
#endif

  d1d2n = -(d1d2ndx * N2x + d1d2ndy * N2y) / d2nd2ndn;
  d1d3n = -(d1d3ndx * N3x + d1d3ndy * N3y) / d3nd3ndn;
  d2d3n = -(d2d3ndx * N3x + d2d3ndy * N3y) / d3nd3ndn;
  d2d1n = -(d2d1ndx * N1x + d2d1ndy * N1y) / d1nd1ndn;
  d3d1n = -(d3d1ndx * N1x + d3d1ndy * N1y) / d1nd1ndn;
  d3d2n = -(d3d2ndx * N2x + d3d2ndy * N2y) / d2nd2ndn;

  // Calculate nodal derivative scaling factors

  d1xd1x = 1. / (d1xd1dx - d1xd1dy * d1yd1dx / d1yd1dy);
  d1xd1y = 1. / (d1yd1dx - d1xd1dx * d1yd1dy / d1xd1dy);
  //  d1xd1y = - d1xd1x * (d1xd1dy / d1yd1dy);
  d1yd1y = 1. / (d1yd1dy - d1yd1dx * d1xd1dy / d1xd1dx);
  d1yd1x = 1. / (d1xd1dy - d1yd1dy * d1xd1dx / d1yd1dx);
  //  d1yd1x = - d1yd1y * (d1yd1dx / d1xd1dx);
  d2xd2x = 1. / (d2xd2dx - d2xd2dy * d2yd2dx / d2yd2dy);
  d2xd2y = 1. / (d2yd2dx - d2xd2dx * d2yd2dy / d2xd2dy);
  //  d2xd2y = - d2xd2x * (d2xd2dy / d2yd2dy);
  d2yd2y = 1. / (d2yd2dy - d2yd2dx * d2xd2dy / d2xd2dx);
  d2yd2x = 1. / (d2xd2dy - d2yd2dy * d2xd2dx / d2yd2dx);
  //  d2yd2x = - d2yd2y * (d2yd2dx / d2xd2dx);
  d3xd3x = 1. / (d3xd3dx - d3xd3dy * d3yd3dx / d3yd3dy);
  d3xd3y = 1. / (d3yd3dx - d3xd3dx * d3yd3dy / d3xd3dy);
  //  d3xd3y = - d3xd3x * (d3xd3dy / d3yd3dy);
  d3yd3y = 1. / (d3yd3dy - d3yd3dx * d3xd3dy / d3xd3dx);
  d3yd3x = 1. / (d3xd3dy - d3yd3dy * d3xd3dx / d3yd3dx);
  //  d3yd3x = - d3yd3y * (d3yd3dx / d3xd3dx);

  //  libMesh::err << d1xd1dx << ' ';
  //  libMesh::err << d1xd1dy << ' ';
  //  libMesh::err << d1yd1dx << ' ';
  //  libMesh::err << d1yd1dy << ' ';
  //  libMesh::err << d2xd2dx << ' ';
  //  libMesh::err << d2xd2dy << ' ';
  //  libMesh::err << d2yd2dx << ' ';
  //  libMesh::err << d2yd2dy << ' ';
  //  libMesh::err << d3xd3dx << ' ';
  //  libMesh::err << d3xd3dy << ' ';
  //  libMesh::err << d3yd3dx << ' ';
  //  libMesh::err << d3yd3dy << std::endl;

  // Calculate midpoint derivative adjustments to nodal derivative
  // interpolant functions

  d1xd2n = -(d1xd1x * d1xd2ndn + d1xd1y * d1yd2ndn) / d2nd2ndn;
  d1yd2n = -(d1yd1y * d1yd2ndn + d1yd1x * d1xd2ndn) / d2nd2ndn;
  d1xd3n = -(d1xd1x * d1xd3ndn + d1xd1y * d1yd3ndn) / d3nd3ndn;
  d1yd3n = -(d1yd1y * d1yd3ndn + d1yd1x * d1xd3ndn) / d3nd3ndn;
  d2xd3n = -(d2xd2x * d2xd3ndn + d2xd2y * d2yd3ndn) / d3nd3ndn;
  d2yd3n = -(d2yd2y * d2yd3ndn + d2yd2x * d2xd3ndn) / d3nd3ndn;
  d2xd1n = -(d2xd2x * d2xd1ndn + d2xd2y * d2yd1ndn) / d1nd1ndn;
  d2yd1n = -(d2yd2y * d2yd1ndn + d2yd2x * d2xd1ndn) / d1nd1ndn;
  d3xd1n = -(d3xd3x * d3xd1ndn + d3xd3y * d3yd1ndn) / d1nd1ndn;
  d3yd1n = -(d3yd3y * d3yd1ndn + d3yd3x * d3xd1ndn) / d1nd1ndn;
  d3xd2n = -(d3xd3x * d3xd2ndn + d3xd3y * d3yd2ndn) / d2nd2ndn;
  d3yd2n = -(d3yd3y * d3yd2ndn + d3yd3x * d3xd2ndn) / d2nd2ndn;

  old_elem_id = elem->id();
  old_elem_ptr = elem;

  // Cross your fingers
  //  libMesh::err << d1nd1ndn << ' ';
  //  libMesh::err << d2xd1ndn << ' ';
  //  libMesh::err << d2yd1ndn << ' ';
  //  libMesh::err << std::endl;

  //  libMesh::err << "Transform variables: ";
  //  libMesh::err << d1nd1n << ' ';
  //  libMesh::err << d2nd2n << ' ';
  //  libMesh::err << d3nd3n << ' ';
  //  libMesh::err << d1d2n << ' ';
  //  libMesh::err << d1d3n << ' ';
  //  libMesh::err << d2d3n << ' ';
  //  libMesh::err << d2d1n << ' ';
  //  libMesh::err << d3d1n << ' ';
  //  libMesh::err << d3d2n << std::endl;
  //  libMesh::err << d1xd1x << ' ';
  //  libMesh::err << d1xd1y << ' ';
  //  libMesh::err << d1yd1x << ' ';
  //  libMesh::err << d1yd1y << ' ';
  //  libMesh::err << d2xd2x << ' ';
  //  libMesh::err << d2xd2y << ' ';
  //  libMesh::err << d2yd2x << ' ';
  //  libMesh::err << d2yd2y << ' ';
  //  libMesh::err << d3xd3x << ' ';
  //  libMesh::err << d3xd3y << ' ';
  //  libMesh::err << d3yd3x << ' ';
  //  libMesh::err << d3yd3y << std::endl;
  //  libMesh::err << d1xd2n << ' ';
  //  libMesh::err << d1yd2n << ' ';
  //  libMesh::err << d1xd3n << ' ';
  //  libMesh::err << d1yd3n << ' ';
  //  libMesh::err << d2xd3n << ' ';
  //  libMesh::err << d2yd3n << ' ';
  //  libMesh::err << d2xd1n << ' ';
  //  libMesh::err << d2yd1n << ' ';
  //  libMesh::err << d3xd1n << ' ';
  //  libMesh::err << d3yd1n << ' ';
  //  libMesh::err << d3xd2n << ' ';
  //  libMesh::err << d3yd2n << ' ';
  //  libMesh::err << std::endl;
  //  libMesh::err << std::endl;
}

template <typename RealType>
unsigned char subtriangle_lookup(const PointTempl<RealType> & p)
{
  if ((p(0) >= p(1)) && (p(0) + 2 * p(1) <= 1))
    return 0;
  if ((p(0) <= p(1)) && (p(1) + 2 * p(0) <= 1))
    return 2;
  return 1;
}


#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
// Return shape function second derivatives on the unit right
// triangle
template <typename RealType>
RealType clough_raw_shape_second_deriv(const unsigned int basis_num,
                                       const unsigned int deriv_type,
                                       const PointTempl<RealType> & p)
{
  RealType xi = p(0), eta = p(1);

  switch (deriv_type)
    {

      // second derivative in xi-xi direction
    case 0:
      switch (basis_num)
        {
        case 0:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -6 + 12*xi;
            case 1:
              return -30 + 42*xi + 42*eta;
            case 2:
              return -6 + 18*xi - 6*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 1:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 6 - 12*xi;
            case 1:
              return 18 - 27*xi - 21*eta;
            case 2:
              return 6 - 15*xi + 3*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 2:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 0;
            case 1:
              return 12 - 15*xi - 21*eta;
            case 2:
              return -3*xi + 3*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 3:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -4 + 6*xi;
            case 1:
              return -9 + 13*xi + 8*eta;
            case 2:
              return -1 - 7*xi + 4*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 4:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 4*eta;
            case 1:
              return 1 - 2*xi + 3*eta;
            case 2:
              return -3 + 14*xi - eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 5:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -2 + 6*xi;
            case 1:
              return -4 + 17./2.*xi + 7./2.*eta;
            case 2:
              return -2 + 13./2.*xi - 1./2.*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 6:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 4*eta;
            case 1:
              return 9 - 23./2.*xi - 23./2.*eta;
            case 2:
              return -1 + 5./2.*xi + 9./2.*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 7:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 0;
            case 1:
              return 7 - 17./2.*xi - 25./2.*eta;
            case 2:
              return 1 - 13./2.*xi + 7./2.*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 8:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 0;
            case 1:
              return -2 + 5./2.*xi + 7./2.*eta;
            case 2:
              return 1./2.*xi - 1./2.*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 9:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 0;
            case 1:
              return std::sqrt(2.) * (8 - 10*xi - 14*eta);
            case 2:
              return std::sqrt(2.) * (-2*xi + 2*eta);

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 10:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 0;
            case 1:
              return -4 + 4*xi + 8*eta;
            case 2:
              return -4 + 20*xi - 8*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 11:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -8*eta;
            case 1:
              return -12 + 16*xi + 12*eta;
            case 2:
              return 4 - 16*xi - 4*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }

        default:
          libmesh_error_msg("Invalid shape function index i = " <<
                            basis_num);
        }

      // second derivative in xi-eta direction
    case 1:
      switch (basis_num)
        {
        case 0:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -6*eta;
            case 1:
              return -30 +42*xi
                + 42*eta;
            case 2:
              return -6*xi;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 1:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return + 3*eta;
            case 1:
              return 15 - 21*xi - 21*eta;
            case 2:
              return 3*xi;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 2:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 3*eta;
            case 1:
              return 15 - 21*xi - 21*eta;
            case 2:
              return 3*xi;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 3:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -eta;
            case 1:
              return -4 + 8*xi + 3*eta;
            case 2:
              return -3 + 4*xi + 4*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 4:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -3 + 4*xi + 4*eta;
            case 1:
              return - 4 + 3*xi + 8*eta;
            case 2:
              return -xi;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 5:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return - 1./2.*eta;
            case 1:
              return -5./2. + 7./2.*xi + 7./2.*eta;
            case 2:
              return - 1./2.*xi;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 6:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -1 + 4*xi + 7./2.*eta;
            case 1:
              return 19./2. - 23./2.*xi - 25./2.*eta;
            case 2:
              return 9./2.*xi;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 7:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 9./2.*eta;
            case 1:
              return 19./2. - 25./2.*xi - 23./2.*eta;
            case 2:
              return -1 + 7./2.*xi + 4*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 8:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -1./2.*eta;
            case 1:
              return -5./2. + 7./2.*xi + 7./2.*eta;
            case 2:
              return -1./2.*xi;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 9:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return std::sqrt(2.) * (2*eta);
            case 1:
              return std::sqrt(2.) * (10 - 14*xi - 14*eta);
            case 2:
              return std::sqrt(2.) * (2*xi);

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 10:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -4*eta;
            case 1:
              return - 8 + 8*xi + 12*eta;
            case 2:
              return 4 - 8*xi - 8*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 11:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 4 - 8*xi - 8*eta;
            case 1:
              return -8 + 12*xi + 8*eta;
            case 2:
              return -4*xi;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }

        default:
          libmesh_error_msg("Invalid shape function index i = " <<
                            basis_num);
        }

      // second derivative in eta-eta direction
    case 2:
      switch (basis_num)
        {
        case 0:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -6 - 6*xi + 18*eta;
            case 1:
              return -30 + 42*xi + 42*eta;
            case 2:
              return -6 + 12*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 1:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 3*xi - 3*eta;
            case 1:
              return 12 - 21*xi - 15*eta;
            case 2:
              return 0;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 2:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 6 + 3*xi - 15*eta;
            case 1:
              return 18 - 21.*xi - 27*eta;
            case 2:
              return 6 - 12*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 3:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -3 - xi + 14*eta;
            case 1:
              return 1 + 3*xi - 2*eta;
            case 2:
              return 4*xi;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 4:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -1 + 4*xi - 7*eta;
            case 1:
              return -9 + 8*xi + 13*eta;
            case 2:
              return -4 + 6*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 5:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return - 1./2.*xi + 1./2.*eta;
            case 1:
              return -2 + 7./2.*xi + 5./2.*eta;
            case 2:
              return 0;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 6:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 1 + 7./2.*xi - 13./2.*eta;
            case 1:
              return 7 - 25./2.*xi - 17./2.*eta;
            case 2:
              return 0;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 7:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -1 + 9./2.*xi + 5./2.*eta;
            case 1:
              return 9 - 23./2.*xi - 23./2.*eta;
            case 2:
              return 4*xi;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 8:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -2 - 1./2.*xi + 13./2.*eta;
            case 1:
              return -4 + 7./2.*xi + 17./2.*eta;
            case 2:
              return -2 + 6*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 9:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return std::sqrt(2.) * (2*xi - 2*eta);
            case 1:
              return std::sqrt(2.) * (8 - 14*xi - 10*eta);
            case 2:
              return 0;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 10:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 4 - 4*xi - 16*eta;
            case 1:
              return -12 + 12*xi + 16*eta;
            case 2:
              return -8*xi;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 11:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -4 - 8*xi + 20*eta;
            case 1:
              return -4 + 8*xi + 4*eta;
            case 2:
              return 0;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }

        default:
          libmesh_error_msg("Invalid shape function index i = " <<
                            basis_num);
        }

    default:
      libmesh_error_msg("Invalid shape function derivative j = " <<
                        deriv_type);
    }
}

#endif



template <typename RealType>
RealType clough_raw_shape_deriv(const unsigned int basis_num,
                                const unsigned int deriv_type,
                                const PointTempl<RealType> & p)
{
  RealType xi = p(0), eta = p(1);

  switch (deriv_type)
    {
      // derivative in xi direction
    case 0:
      switch (basis_num)
        {
        case 0:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -6*xi + 6*xi*xi
                - 3*eta*eta;
            case 1:
              return 9 - 30*xi + 21*xi*xi
                - 30*eta + 42*xi*eta
                + 21*eta*eta;
            case 2:
              return -6*xi + 9*xi*xi
                - 6*xi*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 1:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 6*xi - 6*xi*xi
                + 3./2.*eta*eta;
            case 1:
              return - 9./2. + 18*xi - 27./2.*xi*xi
                + 15*eta - 21*xi*eta
                - 21./2.*eta*eta;
            case 2:
              return 6*xi - 15./2.*xi*xi
                + 3*xi*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 2:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 3./2.*eta*eta;
            case 1:
              return - 9./2. + 12*xi - 15./2.*xi*xi
                + 15*eta - 21*xi*eta
                - 21./2.*eta*eta;
            case 2:
              return -3./2.*xi*xi
                + 3*xi*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 3:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 1 - 4*xi + 3*xi*xi
                - 1./2.*eta*eta;
            case 1:
              return 5./2. - 9*xi + 13./2.*xi*xi
                - 4*eta + 8*xi*eta
                + 3./2.*eta*eta;
            case 2:
              return 1 - xi - 7./2.*xi*xi
                - 3*eta + 4*xi*eta
                + 2*eta*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 4:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return - 3*eta + 4*xi*eta
                + 2*eta*eta;
            case 1:
              return xi - xi*xi
                - 4*eta + 3*xi*eta
                + 4*eta*eta;
            case 2:
              return -3*xi + 7*xi*xi
                - xi*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 5:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -2*xi + 3*xi*xi
                - 1./4.*eta*eta;
            case 1:
              return 3./4. - 4*xi + 17./4.*xi*xi
                - 5./2.*eta + 7./2.*xi*eta
                + 7./4.*eta*eta;
            case 2:
              return -2*xi + 13./4.*xi*xi
                - 1./2.*xi*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 6:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -eta + 4*xi*eta
                + 7./4.*eta*eta;
            case 1:
              return -13./4. + 9*xi - 23./4.*xi*xi
                + 19./2.*eta - 23./2.*xi*eta
                - 25./4.*eta*eta;
            case 2:
              return -xi + 5./4.*xi*xi
                + 9./2.*xi*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 7:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 9./4.*eta*eta;
            case 1:
              return - 11./4. + 7*xi - 17./4.*xi*xi
                + 19./2.*eta - 25./2.*xi*eta
                - 23./4.*eta*eta;
            case 2:
              return xi - 13./4.*xi*xi
                - eta + 7./2.*xi*eta + 2*eta*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 8:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return - 1./4.*eta*eta;
            case 1:
              return 3./4. - 2*xi + 5./4.*xi*xi
                - 5./2.*eta + 7./2.*xi*eta
                + 7./4.*eta*eta;
            case 2:
              return 1./4.*xi*xi
                - 1./2.*xi*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 9:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return std::sqrt(2.) * eta*eta;
            case 1:
              return std::sqrt(2.) * (-3 + 8*xi - 5*xi*xi
                                      + 10*eta - 14*xi*eta
                                      - 7*eta*eta);
            case 2:
              return std::sqrt(2.) * (-xi*xi + 2*xi*eta);

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 10:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -2*eta*eta;
            case 1:
              return 2 - 4*xi + 2*xi*xi
                - 8*eta + 8*xi*eta
                + 6*eta*eta;
            case 2:
              return -4*xi + 10*xi*xi
                + 4*eta - 8*xi*eta
                - 4*eta*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 11:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 4*eta - 8*xi*eta
                - 4*eta*eta;
            case 1:
              return 4 - 12*xi + 8*xi*xi
                - 8*eta + 12*xi*eta
                + 4*eta*eta;
            case 2:
              return 4*xi - 8*xi*xi
                - 4*xi*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }

        default:
          libmesh_error_msg("Invalid shape function index i = " <<
                            basis_num);
        }

      // derivative in eta direction
    case 1:
      switch (basis_num)
        {
        case 0:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return - 6*eta - 6*xi*eta + 9*eta*eta;
            case 1:
              return 9 - 30*xi + 21*xi*xi
                - 30*eta + 42*xi*eta + 21*eta*eta;
            case 2:
              return - 3*xi*xi
                - 6*eta + 6*eta*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 1:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return + 3*xi*eta
                - 3./2.*eta*eta;
            case 1:
              return - 9./2. + 15*xi - 21./2.*xi*xi
                + 12*eta - 21*xi*eta - 15./2.*eta*eta;
            case 2:
              return + 3./2.*xi*xi;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 2:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 6*eta + 3*xi*eta - 15./2.*eta*eta;
            case 1:
              return - 9./2. + 15*xi - 21./2.*xi*xi
                + 18*eta - 21.*xi*eta - 27./2.*eta*eta;
            case 2:
              return 3./2.*xi*xi
                + 6*eta - 6*eta*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 3:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return - 3*eta - xi*eta + 7*eta*eta;
            case 1:
              return - 4*xi + 4*xi*xi
                + eta + 3*xi*eta - eta*eta;
            case 2:
              return - 3*xi + 2*xi*xi
                + 4*xi*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 4:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 1 - 3*xi + 2*xi*xi
                - eta + 4*xi*eta - 7./2.*eta*eta;
            case 1:
              return 5./2. - 4*xi + 3./2.*xi*xi
                - 9.*eta + 8*xi*eta + 13./2.*eta*eta;
            case 2:
              return 1 - 1./2.*xi*xi - 4*eta + 3*eta*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 5:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return - 1./2.*xi*eta + 1./4.*eta*eta;
            case 1:
              return 3./4. - 5./2.*xi + 7./4.*xi*xi
                - 2*eta + 7./2.*xi*eta + 5./4.*eta*eta;
            case 2:
              return - 1./4.*xi*xi;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 6:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -xi + 2*xi*xi
                + eta + 7./2.*xi*eta - 13./4.*eta*eta;
            case 1:
              return - 11./4. + 19./2.*xi - 23./4.*xi*xi
                + 7*eta - 25./2.*xi*eta - 17./4.*eta*eta;
            case 2:
              return 9./4.*xi*xi;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 7:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -eta + 9./2.*xi*eta + 5./4.*eta*eta;
            case 1:
              return - 13./4. + 19./2.*xi - 25./4.*xi*xi
                + 9*eta - 23./2.*xi*eta - 23./4.*eta*eta;
            case 2:
              return - xi + 7./4.*xi*xi + 4*xi*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 8:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return -2*eta - 1./2.*xi*eta + 13./4.*eta*eta;
            case 1:
              return 3./4. - 5./2.*xi + 7./4.*xi*xi
                - 4*eta + 7./2.*xi*eta + 17./4.*eta*eta;
            case 2:
              return - 1./4.*xi*xi
                - 2*eta + 3*eta*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 9:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return std::sqrt(2.) * (2*xi*eta - eta*eta);
            case 1:
              return std::sqrt(2.) * (- 3 + 10*xi - 7*xi*xi
                                      + 8*eta - 14*xi*eta - 5*eta*eta);
            case 2:
              return std::sqrt(2.) * (xi*xi);

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 10:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 4*eta - 4*xi*eta - 8*eta*eta;
            case 1:
              return 4 - 8*xi + 4*xi*xi
                - 12*eta + 12*xi*eta + 8*eta*eta;
            case 2:
              return 4*xi - 4*xi*xi
                - 8*xi*eta;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }
        case 11:
          switch (subtriangle_lookup(p))
            {
            case 0:
              return 4*xi - 4*xi*xi
                - 4*eta - 8*xi*eta + 10.*eta*eta;
            case 1:
              return 2 - 8*xi + 6*xi*xi
                - 4*eta + 8*xi*eta + 2*eta*eta;
            case 2:
              return - 2*xi*xi;

            default:
              libmesh_error_msg("Invalid subtriangle lookup = " <<
                                subtriangle_lookup(p));
            }

        default:
          libmesh_error_msg("Invalid shape function index i = " <<
                            basis_num);
        }

    default:
      libmesh_error_msg("Invalid shape function derivative j = " <<
                        deriv_type);
    }
}

template <typename RealType>
RealType clough_raw_shape(const unsigned int basis_num,
                          const PointTempl<RealType> & p)
{
  RealType xi = p(0), eta = p(1);

  switch (basis_num)
    {
    case 0:
      switch (subtriangle_lookup(p))
        {
        case 0:
          return 1 - 3*xi*xi + 2*xi*xi*xi
            - 3*eta*eta - 3*xi*eta*eta + 3*eta*eta*eta;
        case 1:
          return -1 + 9*xi - 15*xi*xi + 7*xi*xi*xi
            + 9*eta - 30*xi*eta +21*xi*xi*eta
            - 15*eta*eta + 21*xi*eta*eta + 7*eta*eta*eta;
        case 2:
          return 1 - 3*xi*xi + 3*xi*xi*xi
            - 3*xi*xi*eta
            - 3*eta*eta + 2*eta*eta*eta;

        default:
          libmesh_error_msg("Invalid subtriangle lookup = " <<
                            subtriangle_lookup(p));
        }
    case 1:
      switch (subtriangle_lookup(p))
        {
        case 0:
          return 3*xi*xi - 2*xi*xi*xi
            + 3./2.*xi*eta*eta
            - 1./2.*eta*eta*eta;
        case 1:
          return 1 - 9./2.*xi + 9*xi*xi - 9./2.*xi*xi*xi
            - 9./2.*eta + 15*xi*eta - 21./2.*xi*xi*eta
            + 6*eta*eta - 21./2.*xi*eta*eta - 5./2.*eta*eta*eta;
        case 2:
          return 3*xi*xi - 5./2.*xi*xi*xi
            + 3./2.*xi*xi*eta;

        default:
          libmesh_error_msg("Invalid subtriangle lookup = " <<
                            subtriangle_lookup(p));
        }
    case 2:
      switch (subtriangle_lookup(p))
        {
        case 0:
          return 3*eta*eta + 3./2.*xi*eta*eta - 5./2.*eta*eta*eta;
        case 1:
          return 1 - 9./2.*xi + 6*xi*xi - 5./2.*xi*xi*xi
            - 9./2.*eta + 15*xi*eta - 21./2.*xi*xi*eta
            + 9*eta*eta - 21./2.*xi*eta*eta - 9./2.*eta*eta*eta;
        case 2:
          return -1./2.*xi*xi*xi
            + 3./2.*xi*xi*eta
            + 3*eta*eta - 2*eta*eta*eta;

        default:
          libmesh_error_msg("Invalid subtriangle lookup = " <<
                            subtriangle_lookup(p));
        }
    case 3:
      switch (subtriangle_lookup(p))
        {
        case 0:
          return xi - 2*xi*xi + xi*xi*xi
            - 3./2.*eta*eta - 1./2.*xi*eta*eta + 7./3.*eta*eta*eta;
        case 1:
          return -1./6. + 5./2.*xi - 9./2.*xi*xi + 13./6.*xi*xi*xi
            - 4*xi*eta + 4*xi*xi*eta
            + 1./2.*eta*eta + 3./2.*xi*eta*eta - 1./3.*eta*eta*eta;
        case 2:
          return xi - 1./2.*xi*xi - 7./6.*xi*xi*xi
            - 3*xi*eta + 2*xi*xi*eta
            + 2*xi*eta*eta;

        default:
          libmesh_error_msg("Invalid subtriangle lookup = " <<
                            subtriangle_lookup(p));
        }
    case 4:
      switch (subtriangle_lookup(p))
        {
        case 0:
          return eta - 3*xi*eta + 2*xi*xi*eta
            - 1./2.*eta*eta + 2*xi*eta*eta - 7./6.*eta*eta*eta;
        case 1:
          return -1./6. + 1./2.*xi*xi - 1./3.*xi*xi*xi
            + 5./2.*eta - 4*xi*eta + 3./2.*xi*xi*eta
            - 9./2.*eta*eta + 4*xi*eta*eta + 13./6.*eta*eta*eta;
        case 2:
          return -3./2.*xi*xi + 7./3.*xi*xi*xi
            + eta - 1./2.*xi*xi*eta - 2*eta*eta + eta*eta*eta;

        default:
          libmesh_error_msg("Invalid subtriangle lookup = " <<
                            subtriangle_lookup(p));
        }
    case 5:
      switch (subtriangle_lookup(p))
        {
        case 0:
          return -xi*xi + xi*xi*xi
            - 1./4.*xi*eta*eta + 1./12.*eta*eta*eta;
        case 1:
          return -1./6. + 3./4.*xi - 2*xi*xi + 17./12.*xi*xi*xi
            + 3./4.*eta - 5./2.*xi*eta + 7./4.*xi*xi*eta
            - eta*eta + 7./4.*xi*eta*eta + 5./12.*eta*eta*eta;
        case 2:
          return -xi*xi + 13./12.*xi*xi*xi
            - 1./4.*xi*xi*eta;

        default:
          libmesh_error_msg("Invalid subtriangle lookup = " <<
                            subtriangle_lookup(p));
        }
    case 6:
      switch (subtriangle_lookup(p))
        {
        case 0:
          return -xi*eta + 2*xi*xi*eta
            + 1./2.*eta*eta + 7./4.*xi*eta*eta - 13./12.*eta*eta*eta;
        case 1:
          return 2./3. - 13./4.*xi + 9./2.*xi*xi - 23./12.*xi*xi*xi
            - 11./4.*eta + 19./2.*xi*eta - 23./4.*xi*xi*eta
            + 7./2.*eta*eta - 25./4.*xi*eta*eta - 17./12.*eta*eta*eta;
        case 2:
          return -1./2.*xi*xi + 5./12.*xi*xi*xi
            + 9./4.*xi*xi*eta;

        default:
          libmesh_error_msg("Invalid subtriangle lookup = " <<
                            subtriangle_lookup(p));
        }
    case 7:
      switch (subtriangle_lookup(p))
        {
        case 0:
          return -1./2.*eta*eta + 9./4.*xi*eta*eta + 5./12.*eta*eta*eta;
        case 1:
          return 2./3. - 11./4.*xi + 7./2.*xi*xi - 17./12.*xi*xi*xi
            - 13./4.*eta + 19./2.*xi*eta - 25./4.*xi*xi*eta
            + 9./2.*eta*eta - 23./4.*xi*eta*eta - 23./12.*eta*eta*eta;
        case 2:
          return 1./2.*xi*xi - 13./12.*xi*xi*xi
            - xi*eta + 7./4.*xi*xi*eta + 2*xi*eta*eta;

        default:
          libmesh_error_msg("Invalid subtriangle lookup = " <<
                            subtriangle_lookup(p));
        }
    case 8:
      switch (subtriangle_lookup(p))
        {
        case 0:
          return -eta*eta - 1./4.*xi*eta*eta + 13./12.*eta*eta*eta;
        case 1:
          return -1./6. + 3./4.*xi - xi*xi + 5./12.*xi*xi*xi
            + 3./4.*eta - 5./2.*xi*eta + 7./4.*xi*xi*eta
            - 2*eta*eta + 7./4.*xi*eta*eta + 17./12.*eta*eta*eta;
        case 2:
          return 1./12.*xi*xi*xi
            - 1./4.*xi*xi*eta
            - eta*eta + eta*eta*eta;

        default:
          libmesh_error_msg("Invalid subtriangle lookup = " <<
                            subtriangle_lookup(p));
        }
    case 9:
      switch (subtriangle_lookup(p))
        {
        case 0:
          return std::sqrt(2.) * (xi*eta*eta - 1./3.*eta*eta*eta);
        case 1:
          return std::sqrt(2.) * (2./3. - 3*xi + 4*xi*xi - 5./3.*xi*xi*xi
                                  - 3*eta + 10*xi*eta - 7*xi*xi*eta
                                  + 4*eta*eta - 7*xi*eta*eta - 5./3.*eta*eta*eta);
        case 2:
          return std::sqrt(2.) * (-1./3.*xi*xi*xi + xi*xi*eta);

        default:
          libmesh_error_msg("Invalid subtriangle lookup = " <<
                            subtriangle_lookup(p));
        }
    case 10:
      switch (subtriangle_lookup(p))
        {
        case 0:
          return 2*eta*eta - 2*xi*eta*eta - 8./3.*eta*eta*eta;
        case 1:
          return -2./3. + 2*xi - 2*xi*xi + 2./3.*xi*xi*xi
            + 4*eta - 8*xi*eta + 4*xi*xi*eta
            - 6*eta*eta + 6*xi*eta*eta + 8./3.*eta*eta*eta;
        case 2:
          return -2*xi*xi + 10./3.*xi*xi*xi
            + 4*xi*eta - 4*xi*xi*eta
            - 4*xi*eta*eta;

        default:
          libmesh_error_msg("Invalid subtriangle lookup = " <<
                            subtriangle_lookup(p));
        }
    case 11:
      switch (subtriangle_lookup(p))
        {
        case 0:
          return 4*xi*eta - 4*xi*xi*eta
            - 2*eta*eta - 4*xi*eta*eta + 10./3.*eta*eta*eta;
        case 1:
          return -2./3. + 4*xi - 6*xi*xi + 8./3.*xi*xi*xi
            + 2*eta - 8*xi*eta + 6*xi*xi*eta
            - 2*eta*eta + 4*xi*eta*eta + 2./3.*eta*eta*eta;
        case 2:
          return 2*xi*xi - 8./3.*xi*xi*xi
            - 2*xi*xi*eta;

        default:
          libmesh_error_msg("Invalid subtriangle lookup = " <<
                            subtriangle_lookup(p));
        }

    default:
      libmesh_error_msg("Invalid shape function index i = " <<
                        basis_num);
    }
}


} // end anonymous namespace


namespace libMesh
{


template <typename RealType>
typename FEShim<2,CLOUGH,RealType>::OutputShape FEShim<2,CLOUGH,RealType>::shape(const ElemType,
                         const Order,
                         const unsigned int,
                         const Point &)
{
  libmesh_error_msg("Clough-Tocher elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}



template <typename RealType>
typename FEShim<2,CLOUGH,RealType>::OutputShape FEShim<2,CLOUGH,RealType>::shape(const ElemTempl<RealType> * elem,
                                                                                 const Order order,
                                                                                 const unsigned int i,
                                                                                 const PointTempl<RealType> & p,
                                                                                 const bool add_p_level)
{
  libmesh_assert(elem);

  // Coefficient naming: d(1)d(2n) is the coefficient of the
  // global shape function corresponding to value 1 in terms of the
  // local shape function corresponding to normal derivative 2
  static RealType d1d2n, d1d3n, d2d3n, d2d1n, d3d1n, d3d2n;
  static RealType d1xd1x, d1xd1y, d1xd2n, d1xd3n;
  static RealType d1yd1x, d1yd1y, d1yd2n, d1yd3n;
  static RealType d2xd2x, d2xd2y, d2xd3n, d2xd1n;
  static RealType d2yd2x, d2yd2y, d2yd3n, d2yd1n;
  static RealType d3xd3x, d3xd3y, d3xd1n, d3xd2n;
  static RealType d3yd3x, d3yd3y, d3yd1n, d3yd2n;
  static RealType d1nd1n, d2nd2n, d3nd3n;
  // Normal vector naming: N01x is the x component of the
  // unit vector at point 0 normal to (possibly curved) side 01
  static RealType N01x, N01y, N10x, N10y;
  static RealType N02x, N02y, N20x, N20y;
  static RealType N21x, N21y, N12x, N12y;

  clough_compute_coefs(elem,
                       d1d2n,
                       d1d3n,
                       d2d3n,
                       d2d1n,
                       d3d1n,
                       d3d2n,
                       d1xd1x,
                       d1xd1y,
                       d1xd2n,
                       d1xd3n,
                       d1yd1x,
                       d1yd1y,
                       d1yd2n,
                       d1yd3n,
                       d2xd2x,
                       d2xd2y,
                       d2xd3n,
                       d2xd1n,
                       d2yd2x,
                       d2yd2y,
                       d2yd3n,
                       d2yd1n,
                       d3xd3x,
                       d3xd3y,
                       d3xd1n,
                       d3xd2n,
                       d3yd3x,
                       d3yd3y,
                       d3yd1n,
                       d3yd2n,
                       d1nd1n,
                       d2nd2n,
                       d3nd3n,
                       N01x,
                       N01y,
                       N10x,
                       N10y,
                       N02x,
                       N02y,
                       N20x,
                       N20y,
                       N21x,
                       N21y,
                       N12x,
                       N12y);

  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order + add_p_level * elem->p_level());

  switch (totalorder)
    {
      // 2nd-order restricted Clough-Tocher element
    case SECOND:
      {
        // There may be a bug in the 2nd order case; the 3rd order
        // Clough-Tocher elements are pretty uniformly better anyways
        // so use those instead.
        libmesh_experimental();
        switch (type)
          {
            // C1 functions on the Clough-Tocher triangle.
          case TRI6:
            {
              libmesh_assert_less (i, 9);
              // FIXME: it would be nice to calculate (and cache)
              // clough_raw_shape(j,p) only once per triangle, not 1-7
              // times
              switch (i)
                {
                  // Note: these DoF numbers are "scrambled" because my
                  // initial numbering conventions didn't match libMesh
                case 0:
                  return clough_raw_shape(0, p)
                    + d1d2n * clough_raw_shape(10, p)
                    + d1d3n * clough_raw_shape(11, p);
                case 3:
                  return clough_raw_shape(1, p)
                    + d2d3n * clough_raw_shape(11, p)
                    + d2d1n * clough_raw_shape(9, p);
                case 6:
                  return clough_raw_shape(2, p)
                    + d3d1n * clough_raw_shape(9, p)
                    + d3d2n * clough_raw_shape(10, p);
                case 1:
                  return d1xd1x * clough_raw_shape(3, p)
                    + d1xd1y * clough_raw_shape(4, p)
                    + d1xd2n * clough_raw_shape(10, p)
                    + d1xd3n * clough_raw_shape(11, p)
                    + 0.5 * N01x * d3nd3n * clough_raw_shape(11, p)
                    + 0.5 * N02x * d2nd2n * clough_raw_shape(10, p);
                case 2:
                  return d1yd1y * clough_raw_shape(4, p)
                    + d1yd1x * clough_raw_shape(3, p)
                    + d1yd2n * clough_raw_shape(10, p)
                    + d1yd3n * clough_raw_shape(11, p)
                    + 0.5 * N01y * d3nd3n * clough_raw_shape(11, p)
                    + 0.5 * N02y * d2nd2n * clough_raw_shape(10, p);
                case 4:
                  return d2xd2x * clough_raw_shape(5, p)
                    + d2xd2y * clough_raw_shape(6, p)
                    + d2xd3n * clough_raw_shape(11, p)
                    + d2xd1n * clough_raw_shape(9, p)
                    + 0.5 * N10x * d3nd3n * clough_raw_shape(11, p)
                    + 0.5 * N12x * d1nd1n * clough_raw_shape(9, p);
                case 5:
                  return d2yd2y * clough_raw_shape(6, p)
                    + d2yd2x * clough_raw_shape(5, p)
                    + d2yd3n * clough_raw_shape(11, p)
                    + d2yd1n * clough_raw_shape(9, p)
                    + 0.5 * N10y * d3nd3n * clough_raw_shape(11, p)
                    + 0.5 * N12y * d1nd1n * clough_raw_shape(9, p);
                case 7:
                  return d3xd3x * clough_raw_shape(7, p)
                    + d3xd3y * clough_raw_shape(8, p)
                    + d3xd1n * clough_raw_shape(9, p)
                    + d3xd2n * clough_raw_shape(10, p)
                    + 0.5 * N20x * d2nd2n * clough_raw_shape(10, p)
                    + 0.5 * N21x * d1nd1n * clough_raw_shape(9, p);
                case 8:
                  return d3yd3y * clough_raw_shape(8, p)
                    + d3yd3x * clough_raw_shape(7, p)
                    + d3yd1n * clough_raw_shape(9, p)
                    + d3yd2n * clough_raw_shape(10, p)
                    + 0.5 * N20y * d2nd2n * clough_raw_shape(10, p)
                    + 0.5 * N21y * d1nd1n * clough_raw_shape(9, p);
                default:
                  libmesh_error_msg("Invalid shape function index i = " << i);
                }
            }
          default:
            libmesh_error_msg("ERROR: Unsupported element type = " << type);
          }
      }
      // 3rd-order Clough-Tocher element
    case THIRD:
      {
        switch (type)
          {
            // C1 functions on the Clough-Tocher triangle.
          case TRI6:
            {
              libmesh_assert_less (i, 12);

              // FIXME: it would be nice to calculate (and cache)
              // clough_raw_shape(j,p) only once per triangle, not 1-7
              // times
              switch (i)
                {
                  // Note: these DoF numbers are "scrambled" because my
                  // initial numbering conventions didn't match libMesh
                case 0:
                  return clough_raw_shape(0, p)
                    + d1d2n * clough_raw_shape(10, p)
                    + d1d3n * clough_raw_shape(11, p);
                case 3:
                  return clough_raw_shape(1, p)
                    + d2d3n * clough_raw_shape(11, p)
                    + d2d1n * clough_raw_shape(9, p);
                case 6:
                  return clough_raw_shape(2, p)
                    + d3d1n * clough_raw_shape(9, p)
                    + d3d2n * clough_raw_shape(10, p);
                case 1:
                  return d1xd1x * clough_raw_shape(3, p)
                    + d1xd1y * clough_raw_shape(4, p)
                    + d1xd2n * clough_raw_shape(10, p)
                    + d1xd3n * clough_raw_shape(11, p);
                case 2:
                  return d1yd1y * clough_raw_shape(4, p)
                    + d1yd1x * clough_raw_shape(3, p)
                    + d1yd2n * clough_raw_shape(10, p)
                    + d1yd3n * clough_raw_shape(11, p);
                case 4:
                  return d2xd2x * clough_raw_shape(5, p)
                    + d2xd2y * clough_raw_shape(6, p)
                    + d2xd3n * clough_raw_shape(11, p)
                    + d2xd1n * clough_raw_shape(9, p);
                case 5:
                  return d2yd2y * clough_raw_shape(6, p)
                    + d2yd2x * clough_raw_shape(5, p)
                    + d2yd3n * clough_raw_shape(11, p)
                    + d2yd1n * clough_raw_shape(9, p);
                case 7:
                  return d3xd3x * clough_raw_shape(7, p)
                    + d3xd3y * clough_raw_shape(8, p)
                    + d3xd1n * clough_raw_shape(9, p)
                    + d3xd2n * clough_raw_shape(10, p);
                case 8:
                  return d3yd3y * clough_raw_shape(8, p)
                    + d3yd3x * clough_raw_shape(7, p)
                    + d3yd1n * clough_raw_shape(9, p)
                    + d3yd2n * clough_raw_shape(10, p);
                case 10:
                  return d1nd1n * clough_raw_shape(9, p);
                case 11:
                  return d2nd2n * clough_raw_shape(10, p);
                case 9:
                  return d3nd3n * clough_raw_shape(11, p);

                default:
                  libmesh_error_msg("Invalid shape function index i = " << i);
                }
            }
          default:
            libmesh_error_msg("ERROR: Unsupported element type = " << type);
          }
      }
      // by default throw an error
    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order = " << order);
    }
}



template <typename RealType>
typename FEShim<2,CLOUGH,RealType>::OutputShape FEShim<2,CLOUGH,RealType>::shape_deriv(const ElemType,
                               const Order,
                               const unsigned int,
                               const unsigned int,
                               const Point &)
{
  libmesh_error_msg("Clough-Tocher elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}



template <typename RealType>
typename FEShim<2,CLOUGH,RealType>::OutputShape FEShim<2,CLOUGH,RealType>::shape_deriv(const Elem * elem,
                               const Order order,
                               const unsigned int i,
                               const unsigned int j,
                               const Point & p,
                               const bool add_p_level)
{
  libmesh_assert(elem);

  // Coefficient naming: d(1)d(2n) is the coefficient of the
  // global shape function corresponding to value 1 in terms of the
  // local shape function corresponding to normal derivative 2
  static RealType d1d2n, d1d3n, d2d3n, d2d1n, d3d1n, d3d2n;
  static RealType d1xd1x, d1xd1y, d1xd2n, d1xd3n;
  static RealType d1yd1x, d1yd1y, d1yd2n, d1yd3n;
  static RealType d2xd2x, d2xd2y, d2xd3n, d2xd1n;
  static RealType d2yd2x, d2yd2y, d2yd3n, d2yd1n;
  static RealType d3xd3x, d3xd3y, d3xd1n, d3xd2n;
  static RealType d3yd3x, d3yd3y, d3yd1n, d3yd2n;
  static RealType d1nd1n, d2nd2n, d3nd3n;
  // Normal vector naming: N01x is the x component of the
  // unit vector at point 0 normal to (possibly curved) side 01
  static RealType N01x, N01y, N10x, N10y;
  static RealType N02x, N02y, N20x, N20y;
  static RealType N21x, N21y, N12x, N12y;

  clough_compute_coefs(elem,
                       d1d2n,
                       d1d3n,
                       d2d3n,
                       d2d1n,
                       d3d1n,
                       d3d2n,
                       d1xd1x,
                       d1xd1y,
                       d1xd2n,
                       d1xd3n,
                       d1yd1x,
                       d1yd1y,
                       d1yd2n,
                       d1yd3n,
                       d2xd2x,
                       d2xd2y,
                       d2xd3n,
                       d2xd1n,
                       d2yd2x,
                       d2yd2y,
                       d2yd3n,
                       d2yd1n,
                       d3xd3x,
                       d3xd3y,
                       d3xd1n,
                       d3xd2n,
                       d3yd3x,
                       d3yd3y,
                       d3yd1n,
                       d3yd2n,
                       d1nd1n,
                       d2nd2n,
                       d3nd3n,
                       N01x,
                       N01y,
                       N10x,
                       N10y,
                       N02x,
                       N02y,
                       N20x,
                       N20y,
                       N21x,
                       N21y,
                       N12x,
                       N12y);

  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order + add_p_level * elem->p_level());

  switch (totalorder)
    {
      // 2nd-order restricted Clough-Tocher element
    case SECOND:
      {
        // There may be a bug in the 2nd order case; the 3rd order
        // Clough-Tocher elements are pretty uniformly better anyways
        // so use those instead.
        libmesh_experimental();
        switch (type)
          {
            // C1 functions on the Clough-Tocher triangle.
          case TRI6:
            {
              libmesh_assert_less (i, 9);
              // FIXME: it would be nice to calculate (and cache)
              // clough_raw_shape(j,p) only once per triangle, not 1-7
              // times
              switch (i)
                {
                  // Note: these DoF numbers are "scrambled" because my
                  // initial numbering conventions didn't match libMesh
                case 0:
                  return clough_raw_shape_deriv(0, j, p)
                    + d1d2n * clough_raw_shape_deriv(10, j, p)
                    + d1d3n * clough_raw_shape_deriv(11, j, p);
                case 3:
                  return clough_raw_shape_deriv(1, j, p)
                    + d2d3n * clough_raw_shape_deriv(11, j, p)
                    + d2d1n * clough_raw_shape_deriv(9, j, p);
                case 6:
                  return clough_raw_shape_deriv(2, j, p)
                    + d3d1n * clough_raw_shape_deriv(9, j, p)
                    + d3d2n * clough_raw_shape_deriv(10, j, p);
                case 1:
                  return d1xd1x * clough_raw_shape_deriv(3, j, p)
                    + d1xd1y * clough_raw_shape_deriv(4, j, p)
                    + d1xd2n * clough_raw_shape_deriv(10, j, p)
                    + d1xd3n * clough_raw_shape_deriv(11, j, p)
                    + 0.5 * N01x * d3nd3n * clough_raw_shape_deriv(11, j, p)
                    + 0.5 * N02x * d2nd2n * clough_raw_shape_deriv(10, j, p);
                case 2:
                  return d1yd1y * clough_raw_shape_deriv(4, j, p)
                    + d1yd1x * clough_raw_shape_deriv(3, j, p)
                    + d1yd2n * clough_raw_shape_deriv(10, j, p)
                    + d1yd3n * clough_raw_shape_deriv(11, j, p)
                    + 0.5 * N01y * d3nd3n * clough_raw_shape_deriv(11, j, p)
                    + 0.5 * N02y * d2nd2n * clough_raw_shape_deriv(10, j, p);
                case 4:
                  return d2xd2x * clough_raw_shape_deriv(5, j, p)
                    + d2xd2y * clough_raw_shape_deriv(6, j, p)
                    + d2xd3n * clough_raw_shape_deriv(11, j, p)
                    + d2xd1n * clough_raw_shape_deriv(9, j, p)
                    + 0.5 * N10x * d3nd3n * clough_raw_shape_deriv(11, j, p)
                    + 0.5 * N12x * d1nd1n * clough_raw_shape_deriv(9, j, p);
                case 5:
                  return d2yd2y * clough_raw_shape_deriv(6, j, p)
                    + d2yd2x * clough_raw_shape_deriv(5, j, p)
                    + d2yd3n * clough_raw_shape_deriv(11, j, p)
                    + d2yd1n * clough_raw_shape_deriv(9, j, p)
                    + 0.5 * N10y * d3nd3n * clough_raw_shape_deriv(11, j, p)
                    + 0.5 * N12y * d1nd1n * clough_raw_shape_deriv(9, j, p);
                case 7:
                  return d3xd3x * clough_raw_shape_deriv(7, j, p)
                    + d3xd3y * clough_raw_shape_deriv(8, j, p)
                    + d3xd1n * clough_raw_shape_deriv(9, j, p)
                    + d3xd2n * clough_raw_shape_deriv(10, j, p)
                    + 0.5 * N20x * d2nd2n * clough_raw_shape_deriv(10, j, p)
                    + 0.5 * N21x * d1nd1n * clough_raw_shape_deriv(9, j, p);
                case 8:
                  return d3yd3y * clough_raw_shape_deriv(8, j, p)
                    + d3yd3x * clough_raw_shape_deriv(7, j, p)
                    + d3yd1n * clough_raw_shape_deriv(9, j, p)
                    + d3yd2n * clough_raw_shape_deriv(10, j, p)
                    + 0.5 * N20y * d2nd2n * clough_raw_shape_deriv(10, j, p)
                    + 0.5 * N21y * d1nd1n * clough_raw_shape_deriv(9, j, p);
                default:
                  libmesh_error_msg("Invalid shape function index i = " << i);
                }
            }
          default:
            libmesh_error_msg("ERROR: Unsupported element type = " << type);
          }
      }
      // 3rd-order Clough-Tocher element
    case THIRD:
      {
        switch (type)
          {
            // C1 functions on the Clough-Tocher triangle.
          case TRI6:
            {
              libmesh_assert_less (i, 12);

              // FIXME: it would be nice to calculate (and cache)
              // clough_raw_shape(j,p) only once per triangle, not 1-7
              // times
              switch (i)
                {
                  // Note: these DoF numbers are "scrambled" because my
                  // initial numbering conventions didn't match libMesh
                case 0:
                  return clough_raw_shape_deriv(0, j, p)
                    + d1d2n * clough_raw_shape_deriv(10, j, p)
                    + d1d3n * clough_raw_shape_deriv(11, j, p);
                case 3:
                  return clough_raw_shape_deriv(1, j, p)
                    + d2d3n * clough_raw_shape_deriv(11, j, p)
                    + d2d1n * clough_raw_shape_deriv(9, j, p);
                case 6:
                  return clough_raw_shape_deriv(2, j, p)
                    + d3d1n * clough_raw_shape_deriv(9, j, p)
                    + d3d2n * clough_raw_shape_deriv(10, j, p);
                case 1:
                  return d1xd1x * clough_raw_shape_deriv(3, j, p)
                    + d1xd1y * clough_raw_shape_deriv(4, j, p)
                    + d1xd2n * clough_raw_shape_deriv(10, j, p)
                    + d1xd3n * clough_raw_shape_deriv(11, j, p);
                case 2:
                  return d1yd1y * clough_raw_shape_deriv(4, j, p)
                    + d1yd1x * clough_raw_shape_deriv(3, j, p)
                    + d1yd2n * clough_raw_shape_deriv(10, j, p)
                    + d1yd3n * clough_raw_shape_deriv(11, j, p);
                case 4:
                  return d2xd2x * clough_raw_shape_deriv(5, j, p)
                    + d2xd2y * clough_raw_shape_deriv(6, j, p)
                    + d2xd3n * clough_raw_shape_deriv(11, j, p)
                    + d2xd1n * clough_raw_shape_deriv(9, j, p);
                case 5:
                  return d2yd2y * clough_raw_shape_deriv(6, j, p)
                    + d2yd2x * clough_raw_shape_deriv(5, j, p)
                    + d2yd3n * clough_raw_shape_deriv(11, j, p)
                    + d2yd1n * clough_raw_shape_deriv(9, j, p);
                case 7:
                  return d3xd3x * clough_raw_shape_deriv(7, j, p)
                    + d3xd3y * clough_raw_shape_deriv(8, j, p)
                    + d3xd1n * clough_raw_shape_deriv(9, j, p)
                    + d3xd2n * clough_raw_shape_deriv(10, j, p);
                case 8:
                  return d3yd3y * clough_raw_shape_deriv(8, j, p)
                    + d3yd3x * clough_raw_shape_deriv(7, j, p)
                    + d3yd1n * clough_raw_shape_deriv(9, j, p)
                    + d3yd2n * clough_raw_shape_deriv(10, j, p);
                case 10:
                  return d1nd1n * clough_raw_shape_deriv(9, j, p);
                case 11:
                  return d2nd2n * clough_raw_shape_deriv(10, j, p);
                case 9:
                  return d3nd3n * clough_raw_shape_deriv(11, j, p);

                default:
                  libmesh_error_msg("Invalid shape function index i = " << i);
                }
            }
          default:
            libmesh_error_msg("ERROR: Unsupported element type = " << type);
          }
      }
      // by default throw an error
    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order = " << order);
    }
}



#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <typename RealType>
typename FEShim<2,CLOUGH,RealType>::OutputShape FEShim<2,CLOUGH,RealType>::shape_second_deriv(const ElemType,
                                      const Order,
                                      const unsigned int,
                                      const unsigned int,
                                      const Point &)
{
  libmesh_error_msg("Clough-Tocher elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}


template <typename RealType>
typename FEShim<2,CLOUGH,RealType>::OutputShape FEShim<2,CLOUGH,RealType>::shape_second_deriv(const Elem * elem,
                                      const Order order,
                                      const unsigned int i,
                                      const unsigned int j,
                                      const Point & p,
                                      const bool add_p_level)
{
  libmesh_assert(elem);

  // Coefficient naming: d(1)d(2n) is the coefficient of the
  // global shape function corresponding to value 1 in terms of the
  // local shape function corresponding to normal derivative 2
  static RealType d1d2n, d1d3n, d2d3n, d2d1n, d3d1n, d3d2n;
  static RealType d1xd1x, d1xd1y, d1xd2n, d1xd3n;
  static RealType d1yd1x, d1yd1y, d1yd2n, d1yd3n;
  static RealType d2xd2x, d2xd2y, d2xd3n, d2xd1n;
  static RealType d2yd2x, d2yd2y, d2yd3n, d2yd1n;
  static RealType d3xd3x, d3xd3y, d3xd1n, d3xd2n;
  static RealType d3yd3x, d3yd3y, d3yd1n, d3yd2n;
  static RealType d1nd1n, d2nd2n, d3nd3n;
  // Normal vector naming: N01x is the x component of the
  // unit vector at point 0 normal to (possibly curved) side 01
  static RealType N01x, N01y, N10x, N10y;
  static RealType N02x, N02y, N20x, N20y;
  static RealType N21x, N21y, N12x, N12y;

  clough_compute_coefs(elem,
                       d1d2n,
                       d1d3n,
                       d2d3n,
                       d2d1n,
                       d3d1n,
                       d3d2n,
                       d1xd1x,
                       d1xd1y,
                       d1xd2n,
                       d1xd3n,
                       d1yd1x,
                       d1yd1y,
                       d1yd2n,
                       d1yd3n,
                       d2xd2x,
                       d2xd2y,
                       d2xd3n,
                       d2xd1n,
                       d2yd2x,
                       d2yd2y,
                       d2yd3n,
                       d2yd1n,
                       d3xd3x,
                       d3xd3y,
                       d3xd1n,
                       d3xd2n,
                       d3yd3x,
                       d3yd3y,
                       d3yd1n,
                       d3yd2n,
                       d1nd1n,
                       d2nd2n,
                       d3nd3n,
                       N01x,
                       N01y,
                       N10x,
                       N10y,
                       N02x,
                       N02y,
                       N20x,
                       N20y,
                       N21x,
                       N21y,
                       N12x,
                       N12y);

  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order + add_p_level * elem->p_level());

  switch (totalorder)
    {
      // 2nd-order restricted Clough-Tocher element
    case SECOND:
      {
        switch (type)
          {
            // C1 functions on the Clough-Tocher triangle.
          case TRI6:
            {
              libmesh_assert_less (i, 9);
              // FIXME: it would be nice to calculate (and cache)
              // clough_raw_shape(j,p) only once per triangle, not 1-7
              // times
              switch (i)
                {
                  // Note: these DoF numbers are "scrambled" because my
                  // initial numbering conventions didn't match libMesh
                case 0:
                  return clough_raw_shape_second_deriv(0, j, p)
                    + d1d2n * clough_raw_shape_second_deriv(10, j, p)
                    + d1d3n * clough_raw_shape_second_deriv(11, j, p);
                case 3:
                  return clough_raw_shape_second_deriv(1, j, p)
                    + d2d3n * clough_raw_shape_second_deriv(11, j, p)
                    + d2d1n * clough_raw_shape_second_deriv(9, j, p);
                case 6:
                  return clough_raw_shape_second_deriv(2, j, p)
                    + d3d1n * clough_raw_shape_second_deriv(9, j, p)
                    + d3d2n * clough_raw_shape_second_deriv(10, j, p);
                case 1:
                  return d1xd1x * clough_raw_shape_second_deriv(3, j, p)
                    + d1xd1y * clough_raw_shape_second_deriv(4, j, p)
                    + d1xd2n * clough_raw_shape_second_deriv(10, j, p)
                    + d1xd3n * clough_raw_shape_second_deriv(11, j, p)
                    + 0.5 * N01x * d3nd3n * clough_raw_shape_second_deriv(11, j, p)
                    + 0.5 * N02x * d2nd2n * clough_raw_shape_second_deriv(10, j, p);
                case 2:
                  return d1yd1y * clough_raw_shape_second_deriv(4, j, p)
                    + d1yd1x * clough_raw_shape_second_deriv(3, j, p)
                    + d1yd2n * clough_raw_shape_second_deriv(10, j, p)
                    + d1yd3n * clough_raw_shape_second_deriv(11, j, p)
                    + 0.5 * N01y * d3nd3n * clough_raw_shape_second_deriv(11, j, p)
                    + 0.5 * N02y * d2nd2n * clough_raw_shape_second_deriv(10, j, p);
                case 4:
                  return d2xd2x * clough_raw_shape_second_deriv(5, j, p)
                    + d2xd2y * clough_raw_shape_second_deriv(6, j, p)
                    + d2xd3n * clough_raw_shape_second_deriv(11, j, p)
                    + d2xd1n * clough_raw_shape_second_deriv(9, j, p)
                    + 0.5 * N10x * d3nd3n * clough_raw_shape_second_deriv(11, j, p)
                    + 0.5 * N12x * d1nd1n * clough_raw_shape_second_deriv(9, j, p);
                case 5:
                  return d2yd2y * clough_raw_shape_second_deriv(6, j, p)
                    + d2yd2x * clough_raw_shape_second_deriv(5, j, p)
                    + d2yd3n * clough_raw_shape_second_deriv(11, j, p)
                    + d2yd1n * clough_raw_shape_second_deriv(9, j, p)
                    + 0.5 * N10y * d3nd3n * clough_raw_shape_second_deriv(11, j, p)
                    + 0.5 * N12y * d1nd1n * clough_raw_shape_second_deriv(9, j, p);
                case 7:
                  return d3xd3x * clough_raw_shape_second_deriv(7, j, p)
                    + d3xd3y * clough_raw_shape_second_deriv(8, j, p)
                    + d3xd1n * clough_raw_shape_second_deriv(9, j, p)
                    + d3xd2n * clough_raw_shape_second_deriv(10, j, p)
                    + 0.5 * N20x * d2nd2n * clough_raw_shape_second_deriv(10, j, p)
                    + 0.5 * N21x * d1nd1n * clough_raw_shape_second_deriv(9, j, p);
                case 8:
                  return d3yd3y * clough_raw_shape_second_deriv(8, j, p)
                    + d3yd3x * clough_raw_shape_second_deriv(7, j, p)
                    + d3yd1n * clough_raw_shape_second_deriv(9, j, p)
                    + d3yd2n * clough_raw_shape_second_deriv(10, j, p)
                    + 0.5 * N20y * d2nd2n * clough_raw_shape_second_deriv(10, j, p)
                    + 0.5 * N21y * d1nd1n * clough_raw_shape_second_deriv(9, j, p);
                default:
                  libmesh_error_msg("Invalid shape function index i = " << i);
                }
            }
          default:
            libmesh_error_msg("ERROR: Unsupported element type = " << type);
          }
      }
      // 3rd-order Clough-Tocher element
    case THIRD:
      {
        switch (type)
          {
            // C1 functions on the Clough-Tocher triangle.
          case TRI6:
            {
              libmesh_assert_less (i, 12);

              // FIXME: it would be nice to calculate (and cache)
              // clough_raw_shape(j,p) only once per triangle, not 1-7
              // times
              switch (i)
                {
                  // Note: these DoF numbers are "scrambled" because my
                  // initial numbering conventions didn't match libMesh
                case 0:
                  return clough_raw_shape_second_deriv(0, j, p)
                    + d1d2n * clough_raw_shape_second_deriv(10, j, p)
                    + d1d3n * clough_raw_shape_second_deriv(11, j, p);
                case 3:
                  return clough_raw_shape_second_deriv(1, j, p)
                    + d2d3n * clough_raw_shape_second_deriv(11, j, p)
                    + d2d1n * clough_raw_shape_second_deriv(9, j, p);
                case 6:
                  return clough_raw_shape_second_deriv(2, j, p)
                    + d3d1n * clough_raw_shape_second_deriv(9, j, p)
                    + d3d2n * clough_raw_shape_second_deriv(10, j, p);
                case 1:
                  return d1xd1x * clough_raw_shape_second_deriv(3, j, p)
                    + d1xd1y * clough_raw_shape_second_deriv(4, j, p)
                    + d1xd2n * clough_raw_shape_second_deriv(10, j, p)
                    + d1xd3n * clough_raw_shape_second_deriv(11, j, p);
                case 2:
                  return d1yd1y * clough_raw_shape_second_deriv(4, j, p)
                    + d1yd1x * clough_raw_shape_second_deriv(3, j, p)
                    + d1yd2n * clough_raw_shape_second_deriv(10, j, p)
                    + d1yd3n * clough_raw_shape_second_deriv(11, j, p);
                case 4:
                  return d2xd2x * clough_raw_shape_second_deriv(5, j, p)
                    + d2xd2y * clough_raw_shape_second_deriv(6, j, p)
                    + d2xd3n * clough_raw_shape_second_deriv(11, j, p)
                    + d2xd1n * clough_raw_shape_second_deriv(9, j, p);
                case 5:
                  return d2yd2y * clough_raw_shape_second_deriv(6, j, p)
                    + d2yd2x * clough_raw_shape_second_deriv(5, j, p)
                    + d2yd3n * clough_raw_shape_second_deriv(11, j, p)
                    + d2yd1n * clough_raw_shape_second_deriv(9, j, p);
                case 7:
                  return d3xd3x * clough_raw_shape_second_deriv(7, j, p)
                    + d3xd3y * clough_raw_shape_second_deriv(8, j, p)
                    + d3xd1n * clough_raw_shape_second_deriv(9, j, p)
                    + d3xd2n * clough_raw_shape_second_deriv(10, j, p);
                case 8:
                  return d3yd3y * clough_raw_shape_second_deriv(8, j, p)
                    + d3yd3x * clough_raw_shape_second_deriv(7, j, p)
                    + d3yd1n * clough_raw_shape_second_deriv(9, j, p)
                    + d3yd2n * clough_raw_shape_second_deriv(10, j, p);
                case 10:
                  return d1nd1n * clough_raw_shape_second_deriv(9, j, p);
                case 11:
                  return d2nd2n * clough_raw_shape_second_deriv(10, j, p);
                case 9:
                  return d3nd3n * clough_raw_shape_second_deriv(11, j, p);

                default:
                  libmesh_error_msg("Invalid shape function index i = " << i);
                }
            }
          default:
            libmesh_error_msg("ERROR: Unsupported element type = " << type);
          }
      }
      // by default throw an error
    default:
      libmesh_error_msg("ERROR: Unsupported polynomial order = " << order);
    }
}

#endif

} // namespace libMesh

#endif // LIBMESH_FE_CLOUGH_IMPL_SHAPE_2D_H
