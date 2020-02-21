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


// C++ includes

// Local includes
#include "libmesh/fe.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/number_lookups.h"


namespace
{
using namespace libMesh;

// Compute the static coefficients for an element
void hermite_compute_coefs(const Elem * elem, std::vector<std::vector<Real>> & dxdxi
#ifdef DEBUG
                           , std::vector<Real> & dxdeta, std::vector<Real> & dydxi
#endif
                           )
{
  const FEFamily mapping_family = FEMap::map_fe_type(*elem);
  const ElemType mapping_elem_type (elem->type());

  const FEType map_fe_type(elem->default_order(), mapping_family);

  const int n_mapping_shape_functions =
    FEInterface::n_shape_functions (2, map_fe_type,
                                    mapping_elem_type);

  std::vector<Point> dofpt;
  dofpt.push_back(Point(-1,-1));
  dofpt.push_back(Point(1,1));

  FEInterface::shape_deriv_ptr shape_deriv_ptr =
    FEInterface::shape_deriv_function(2, map_fe_type, mapping_elem_type);

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
          const Real ddxi = shape_deriv_ptr
            (map_fe_type, elem, i, 0, dofpt[p], false);
          const Real ddeta = shape_deriv_ptr
            (map_fe_type, elem, i, 1, dofpt[p], false);

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
      libmesh_assert_less (std::abs(dxdeta[p]), 1e-9);
      libmesh_assert_less (std::abs(dydxi[p]), 1e-9);
#endif
    }
}



Real hermite_bases_2D (std::vector<unsigned int> & bases1D,
                       const std::vector<std::vector<Real>> & dxdxi,
                       const Order & o,
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
        default:
          libmesh_error_msg("Invalid basis node = " << i/4);
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
          libmesh_error_msg("Invalid basisnum = " << basisnum);
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
          libmesh_error_msg("Invalid basisnum = " << basisnum);
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


LIBMESH_DEFAULT_VECTORIZED_FE(2,HERMITE)


template <>
Real FE<2,HERMITE>::shape(const Elem * elem,
                          const Order order,
                          const unsigned int i,
                          const Point & p,
                          const bool add_p_level)
{
  libmesh_assert(elem);

  std::vector<std::vector<Real>> dxdxi(2, std::vector<Real>(2, 0));

#ifdef DEBUG
  std::vector<Real> dxdeta(2), dydxi(2);
#endif

  hermite_compute_coefs(elem,dxdxi
#ifdef DEBUG
                        ,dxdeta,dydxi
#endif
                        );

  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order + add_p_level * elem->p_level());

  switch (type)
    {
    case QUAD4:
    case QUADSHELL4:
      libmesh_assert_less (totalorder, 4);
      libmesh_fallthrough();
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
      {
        libmesh_assert_less (i, (totalorder+1u)*(totalorder+1u));

        std::vector<unsigned int> bases1D;

        Real coef = hermite_bases_2D(bases1D, dxdxi, totalorder, i);

        return coef * FEHermite<1>::hermite_raw_shape(bases1D[0],p(0)) *
          FEHermite<1>::hermite_raw_shape(bases1D[1],p(1));
      }
    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << type);
    }
}



template <>
Real FE<2,HERMITE>::shape(const ElemType,
                          const Order,
                          const unsigned int,
                          const Point &)
{
  libmesh_error_msg("Hermite elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}


template <>
Real FE<2,HERMITE>::shape(const FEType fet,
                          const Elem * elem,
                          const unsigned int i,
                          const Point & p,
                          const bool add_p_level)
{
  return FE<2,HERMITE>::shape(elem, fet.order, i, p, add_p_level);
}


template <>
Real FE<2,HERMITE>::shape_deriv(const Elem * elem,
                                const Order order,
                                const unsigned int i,
                                const unsigned int j,
                                const Point & p,
                                const bool add_p_level)
{
  libmesh_assert(elem);
  libmesh_assert (j == 0 || j == 1);

  std::vector<std::vector<Real>> dxdxi(2, std::vector<Real>(2, 0));

#ifdef DEBUG
  std::vector<Real> dxdeta(2), dydxi(2);
#endif

  hermite_compute_coefs(elem,dxdxi
#ifdef DEBUG
                        ,dxdeta,dydxi
#endif
                        );

  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order + add_p_level * elem->p_level());

  switch (type)
    {
    case QUAD4:
    case QUADSHELL4:
      libmesh_assert_less (totalorder, 4);
      libmesh_fallthrough();
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
      {
        libmesh_assert_less (i, (totalorder+1u)*(totalorder+1u));

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
            libmesh_error_msg("Invalid derivative index j = " << j);
          }
      }
    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << type);
    }
}


template <>
Real FE<2,HERMITE>::shape_deriv(const ElemType,
                                const Order,
                                const unsigned int,
                                const unsigned int,
                                const Point &)
{
  libmesh_error_msg("Hermite elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}


template <>
Real FE<2,HERMITE>::shape_deriv(const FEType fet,
                                const Elem * elem,
                                const unsigned int i,
                                const unsigned int j,
                                const Point & p,
                                const bool add_p_level)
{
  return FE<2,HERMITE>::shape_deriv(elem, fet.order, i, j, p, add_p_level);
}




#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES


template <>
Real FE<2,HERMITE>::shape_second_deriv(const Elem * elem,
                                       const Order order,
                                       const unsigned int i,
                                       const unsigned int j,
                                       const Point & p,
                                       const bool add_p_level)
{
  libmesh_assert(elem);
  libmesh_assert (j == 0 || j == 1 || j == 2);

  std::vector<std::vector<Real>> dxdxi(2, std::vector<Real>(2, 0));

#ifdef DEBUG
  std::vector<Real> dxdeta(2), dydxi(2);
#endif

  hermite_compute_coefs(elem,dxdxi
#ifdef DEBUG
                        ,dxdeta,dydxi
#endif
                        );

  const ElemType type = elem->type();

  const Order totalorder =
    static_cast<Order>(order + add_p_level * elem->p_level());

  switch (type)
    {
    case QUAD4:
    case QUADSHELL4:
      libmesh_assert_less (totalorder, 4);
      libmesh_fallthrough();
    case QUAD8:
    case QUADSHELL8:
    case QUAD9:
      {
        libmesh_assert_less (i, (totalorder+1u)*(totalorder+1u));

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
            libmesh_error_msg("Invalid derivative index j = " << j);
          }
      }
    default:
      libmesh_error_msg("ERROR: Unsupported element type = " << type);
    }
}


template <>
Real FE<2,HERMITE>::shape_second_deriv(const ElemType,
                                       const Order,
                                       const unsigned int,
                                       const unsigned int,
                                       const Point &)
{
  libmesh_error_msg("Hermite elements require the real element \nto construct gradient-based degrees of freedom.");
  return 0.;
}


template <>
Real FE<2,HERMITE>::shape_second_deriv(const FEType fet,
                                       const Elem * elem,
                                       const unsigned int i,
                                       const unsigned int j,
                                       const Point & p,
                                       const bool add_p_level)
{
  return FE<2,HERMITE>::shape_second_deriv(elem, fet.order, i, j, p, add_p_level);
}


#endif

} // namespace libMesh
