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


// C++ inlcludes

// Local includes
#include "fe.h"
#include "elem.h"
#include "number_lookups.h"


// Anonymous namespace for persistant variables.
// This allows us to cache the global-to-local mapping transformation
// This caching is turned off when TBB is enabled...
namespace
{
  using namespace libMesh;

#ifndef LIBMESH_HAVE_TBB_API
  static unsigned int old_elem_id = libMesh::invalid_uint;
  // Mapping functions - derivatives at each dofpt
  std::vector<std::vector<Real> > dxdxi(2, std::vector<Real>(2, 0));
#endif //LIBMESH_HAVE_TBB_API


#ifndef LIBMESH_HAVE_TBB_API
  // Compute the static coefficients for an element
  void hermite_compute_coefs(const Elem* elem)
  {
    // Coefficients are cached from old elements
    if (elem->id() == old_elem_id)
      return;

    old_elem_id = elem->id();
#else
  void hermite_compute_coefs(const Elem* elem, std::vector<std::vector<Real> > & dxdxi)
  {
#endif //LIBMESH_HAVE_TBB_API

#ifdef DEBUG
  std::vector<Real> dxdeta(2), dydxi(2);
#endif //DEBUG
  const Order mapping_order        (elem->default_order());
  const ElemType mapping_elem_type (elem->type());
  const int n_mapping_shape_functions =
    FE<2,LAGRANGE>::n_shape_functions(mapping_elem_type,
				      mapping_order);

  std::vector<Point> dofpt;
  dofpt.push_back(Point(-1,-1));
  dofpt.push_back(Point(1,1));

  for (int p = 0; p != 2; ++p)
    {
      dxdxi[0][p] = 0;
      dxdxi[1][p] = 0;
#ifdef DEBUG
      dxdeta[p] = 0;
      dydxi[p] = 0;
#endif
      for (int i = 0; i != n_mapping_shape_functions; ++i)
        {
          const Real ddxi = FE<2,LAGRANGE>::shape_deriv
            (mapping_elem_type, mapping_order, i, 0, dofpt[p]);
          const Real ddeta = FE<2,LAGRANGE>::shape_deriv
            (mapping_elem_type, mapping_order, i, 1, dofpt[p]);

          dxdxi[0][p] += elem->point(i)(0) * ddxi;
          dxdxi[1][p] += elem->point(i)(1) * ddeta;
// dxdeta and dydxi should be 0!
#ifdef DEBUG
          dxdeta[p] += elem->point(i)(0) * ddeta;
          dydxi[p] += elem->point(i)(1) * ddxi;
#endif
        }
      // No singular elements!
      libmesh_assert(dxdxi[0][p]);
      libmesh_assert(dxdxi[1][p]);
      // No non-rectilinear or non-axis-aligned elements!
#ifdef DEBUG
      libmesh_assert(std::abs(dxdeta[p]) < 1e-9);
      libmesh_assert(std::abs(dydxi[p]) < 1e-9);
#endif
    }
}



Real hermite_bases_2D
 (std::vector<unsigned int> &bases1D,
  const std::vector<std::vector<Real> > &dxdxi,
  const Order &o,
  unsigned int i)
{
  bases1D.clear();
  bases1D.resize(2,0);
  Real coef = 1.0;

  unsigned int e = o-3;

  // Nodes
  if (i < 16)
    {
      switch (i / 4)
        {
        case 0:
          break;
        case 1:
          bases1D[0] = 1;
          break;
        case 2:
          bases1D[0] = 1;
          bases1D[1] = 1;
          break;
        case 3:
          bases1D[1] = 1;
          break;
        }

      unsigned int basisnum = i%4;
      switch (basisnum)
        {
          case 0: // DoF = value at node
            coef = 1.0;
            break;
          case 1: // DoF = x derivative at node
            coef = dxdxi[0][bases1D[0]];
            bases1D[0] += 2; break;
          case 2: // DoF = y derivative at node
            coef = dxdxi[1][bases1D[1]];
            bases1D[1] += 2; break;
          case 3: // DoF = xy derivative at node
            coef = dxdxi[0][bases1D[0]] * dxdxi[1][bases1D[1]];
            bases1D[0] += 2; bases1D[1] += 2; break;
          default:
            libmesh_error();
        }
    }
  // Edges
  else if (i < 16 + 4*2*e)
    {
      unsigned int basisnum = (i - 16) % (2*e);
      switch ((i - 16) / (2*e))
        {
        case 0:
          bases1D[0] = basisnum/2 + 4;
          bases1D[1] = basisnum%2*2;
          if (basisnum%2)
            coef = dxdxi[1][0];
          break;
        case 1:
          bases1D[0] = basisnum%2*2 + 1;
          bases1D[1] = basisnum/2 + 4;
          if (basisnum%2)
            coef = dxdxi[0][1];
          break;
        case 2:
          bases1D[0] = basisnum/2 + 4;
          bases1D[1] = basisnum%2*2 + 1;
          if (basisnum%2)
            coef = dxdxi[1][1];
          break;
        case 3:
          bases1D[0] = basisnum%2*2;
          bases1D[1] = basisnum/2 + 4;
          if (basisnum%2)
            coef = dxdxi[0][0];
          break;
        default:
          libmesh_error();
        }
    }
  // Interior
  else
    {
      unsigned int basisnum = i - 16 - 8*e;
      bases1D[0] = square_number_row[basisnum]+4;
      bases1D[1] = square_number_column[basisnum]+4;
    }

  // No singular elements
  libmesh_assert(coef);
  return coef;
}

} // end anonymous namespace


namespace libMesh
{


template <>
Real FE<2,HERMITE>::shape(const ElemType,
			  const Order,
			  const unsigned int,
			  const Point&)
{
  libMesh::err << "Hermite elements require the real element\n"
	        << "to construct gradient-based degrees of freedom."
	        << std::endl;

  libmesh_error();
  return 0.;
}



template <>
Real FE<2,HERMITE>::shape(const Elem* elem,
			  const Order order,
			  const unsigned int i,
			  const Point& p)
{
  libmesh_assert (elem != NULL);

#ifndef LIBMESH_HAVE_TBB_API
  hermite_compute_coefs(elem);
#else
  std::vector<std::vector<Real> > dxdxi(2, std::vector<Real>(2, 0));

  hermite_compute_coefs(elem, dxdxi);
#endif //LIBMESH_HAVE_TBB_API

  const ElemType type = elem->type();

  const Order totalorder = static_cast<Order>(order + elem->p_level());

  switch (type)
    {
    case QUAD4:
      libmesh_assert (totalorder < 4);
    case QUAD8:
    case QUAD9:
      {
        libmesh_assert (i<(totalorder+1u)*(totalorder+1u));

        std::vector<unsigned int> bases1D;

        Real coef = hermite_bases_2D(bases1D, dxdxi, totalorder, i);

        return coef * FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) *
               FEHermite<1>::hermite_raw_shape(bases1D[1],p(1));
      }
    default:
      libMesh::err << "ERROR: Unsupported element type!" << std::endl;
      libmesh_error();
    }

  libmesh_error();
  return 0.;
}



template <>
Real FE<2,HERMITE>::shape_deriv(const ElemType,
				const Order,
				const unsigned int,
				const unsigned int,
				const Point&)
{
  libMesh::err << "Hermite elements require the real element\n"
	        << "to construct gradient-based degrees of freedom."
	        << std::endl;

  libmesh_error();
  return 0.;
}



template <>
Real FE<2,HERMITE>::shape_deriv(const Elem* elem,
				const Order order,
				const unsigned int i,
				const unsigned int j,
				const Point& p)
{
  libmesh_assert (elem != NULL);
  libmesh_assert (j == 0 || j == 1);

#ifndef LIBMESH_HAVE_TBB_API
  hermite_compute_coefs(elem);
#else
  std::vector<std::vector<Real> > dxdxi(2, std::vector<Real>(2, 0));

  hermite_compute_coefs(elem, dxdxi);
#endif //LIBMESH_HAVE_TBB_API

  const ElemType type = elem->type();

  const Order totalorder = static_cast<Order>(order + elem->p_level());

  switch (type)
    {
    case QUAD4:
      libmesh_assert (totalorder < 4);
    case QUAD8:
    case QUAD9:
      {
        libmesh_assert (i<(totalorder+1u)*(totalorder+1u));

        std::vector<unsigned int> bases1D;

        Real coef = hermite_bases_2D(bases1D, dxdxi, totalorder, i);

        switch (j)
          {
            case 0:
              return coef *
                FEHermite<1>::hermite_raw_shape_deriv(bases1D[0],p(0)) *
                FEHermite<1>::hermite_raw_shape(bases1D[1],p(1));
            case 1:
              return coef *
                FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) *
                FEHermite<1>::hermite_raw_shape_deriv(bases1D[1],p(1));
            default:
              libmesh_error();
          }
      }
    default:
      libMesh::err << "ERROR: Unsupported element type!" << std::endl;
      libmesh_error();
    }

  libmesh_error();
  return 0.;
}



template <>
Real FE<2,HERMITE>::shape_second_deriv(const Elem* elem,
                                       const Order order,
                                       const unsigned int i,
                                       const unsigned int j,
                                       const Point& p)
{
  libmesh_assert (elem != NULL);
  libmesh_assert (j == 0 || j == 1 || j == 2);

#ifndef LIBMESH_HAVE_TBB_API
  hermite_compute_coefs(elem);
#else
  std::vector<std::vector<Real> > dxdxi(2, std::vector<Real>(2, 0));

  hermite_compute_coefs(elem, dxdxi);
#endif //LIBMESH_HAVE_TBB_API

  const ElemType type = elem->type();

  const Order totalorder = static_cast<Order>(order + elem->p_level());

  switch (type)
    {
    case QUAD4:
      libmesh_assert (totalorder < 4);
    case QUAD8:
    case QUAD9:
      {
        libmesh_assert (i<(totalorder+1u)*(totalorder+1u));

        std::vector<unsigned int> bases1D;

        Real coef = hermite_bases_2D(bases1D, dxdxi, totalorder, i);

        switch (j)
          {
          case 0:
            return coef *
              FEHermite<1>::hermite_raw_shape_second_deriv(bases1D[0],p(0)) *
              FEHermite<1>::hermite_raw_shape(bases1D[1],p(1));
          case 1:
            return coef *
              FEHermite<1>::hermite_raw_shape_deriv(bases1D[0],p(0)) *
              FEHermite<1>::hermite_raw_shape_deriv(bases1D[1],p(1));
          case 2:
            return coef *
              FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) *
              FEHermite<1>::hermite_raw_shape_second_deriv(bases1D[1],p(1));
          default:
            libmesh_error();
          }
      }
    default:
      libMesh::err << "ERROR: Unsupported element type!" << std::endl;
      libmesh_error();
    }

  libmesh_error();
  return 0.;
}

} // namespace libMesh
