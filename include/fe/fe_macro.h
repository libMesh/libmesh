// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


#ifndef LIBMESH_FE_MACRO_H
#define LIBMESH_FE_MACRO_H




/**
 * These macros help in instantiating specific versions
 * of the \p FE class.  Simply include this file, and
 * instantiate at the end for the desired dimension(s).
 */

#define INSTANTIATE_MAPS(_dim,_type)                                    \
  template Point FE<_dim,_type>::map(const Elem*,const Point&);         \
  template Point FE<_dim,_type>::map_xi(const Elem*,const Point&);      \
  template Point FE<_dim,_type>::map_eta(const Elem*,const Point&);     \
  template Point FE<_dim,_type>::map_zeta(const Elem*,const Point&);    \
  template void  FE<_dim,_type>::inverse_map(const Elem*,const std::vector<Point>&,std::vector<Point>&,Real,bool); \
  template Point FE<_dim,_type>::inverse_map(const Elem*,const Point&,Real,bool)

#define INSTANTIATE_SUBDIVISION_MAPS                                    \
  template Point FE<2,SUBDIVISION>::map(const Elem*,const Point&);      \
  template Point FE<2,SUBDIVISION>::map_xi(const Elem*,const Point&);   \
  template Point FE<2,SUBDIVISION>::map_eta(const Elem*,const Point&);  \
  template Point FE<2,SUBDIVISION>::map_zeta(const Elem*,const Point&)

#define INSTANTIATE_SUBDIVISION_FE                                      \
  template unsigned int FE<2,SUBDIVISION>::n_shape_functions () const;  \
  template void         FE<2,SUBDIVISION>::attach_quadrature_rule (QBase*); \
  template unsigned int FE<2,SUBDIVISION>::n_quadrature_points () const; \
  template void         FE<2,SUBDIVISION>::reinit(const Elem*,const std::vector<Point>* const,const std::vector<Real>* const); \
  template void         FE<2,SUBDIVISION>::init_shape_functions(const std::vector<Point>&, const Elem*)


#ifndef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#define INSTANTIATE_FE(_dim)   template class FE< (_dim), CLOUGH>;      \
  template class FE< (_dim), HERMITE>;                                  \
  template class FE< (_dim), HIERARCHIC>;                               \
  template class FE< (_dim), L2_HIERARCHIC>;                            \
  template class FE< (_dim), LAGRANGE>;                                 \
  template class FE< (_dim), LAGRANGE_VEC>;                             \
  template class FE< (_dim), L2_LAGRANGE>;                              \
  template class FE< (_dim), MONOMIAL>;                                 \
  template class FE< (_dim), SCALAR>;                                   \
  template class FE< (_dim), XYZ>;                                      \
  template class FE< (_dim), NEDELEC_ONE>

#define INSTANTIATE_ALL_MAPS(_dim)              \
  INSTANTIATE_MAPS(_dim,CLOUGH);                \
  INSTANTIATE_MAPS(_dim,HERMITE);               \
  INSTANTIATE_MAPS(_dim,HIERARCHIC);            \
  INSTANTIATE_MAPS(_dim,L2_HIERARCHIC);         \
  INSTANTIATE_MAPS(_dim,LAGRANGE);              \
  INSTANTIATE_MAPS(_dim,LAGRANGE_VEC);          \
  INSTANTIATE_MAPS(_dim,L2_LAGRANGE);           \
  INSTANTIATE_MAPS(_dim,MONOMIAL);              \
  INSTANTIATE_MAPS(_dim,SCALAR);                \
  INSTANTIATE_MAPS(_dim,XYZ);                   \
  INSTANTIATE_MAPS(_dim,NEDELEC_ONE)

#else //LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#define INSTANTIATE_FE(_dim)   template class FE< (_dim), CLOUGH>;      \
  template class FE< (_dim), HERMITE>;                                  \
  template class FE< (_dim), HIERARCHIC>;                               \
  template class FE< (_dim), L2_HIERARCHIC>;                            \
  template class FE< (_dim), LAGRANGE>;                                 \
  template class FE< (_dim), LAGRANGE_VEC>;                             \
  template class FE< (_dim), L2_LAGRANGE>;                              \
  template class FE< (_dim), MONOMIAL>;                                 \
  template class FE< (_dim), SCALAR>;                                   \
  template class FE< (_dim), BERNSTEIN>;                                \
  template class FE< (_dim), SZABAB>;                                   \
  template class FE< (_dim), XYZ>;                                      \
  template class FE< (_dim), NEDELEC_ONE>

#define INSTANTIATE_ALL_MAPS(_dim)              \
  INSTANTIATE_MAPS(_dim,CLOUGH);                \
  INSTANTIATE_MAPS(_dim,HERMITE);               \
  INSTANTIATE_MAPS(_dim,HIERARCHIC);            \
  INSTANTIATE_MAPS(_dim,L2_HIERARCHIC);         \
  INSTANTIATE_MAPS(_dim,LAGRANGE);              \
  INSTANTIATE_MAPS(_dim,LAGRANGE_VEC);          \
  INSTANTIATE_MAPS(_dim,L2_LAGRANGE);           \
  INSTANTIATE_MAPS(_dim,MONOMIAL);              \
  INSTANTIATE_MAPS(_dim,SCALAR);                \
  INSTANTIATE_MAPS(_dim,BERNSTEIN);             \
  INSTANTIATE_MAPS(_dim,SZABAB);                \
  INSTANTIATE_MAPS(_dim,XYZ);                   \
  INSTANTIATE_MAPS(_dim,NEDELEC_ONE)

#endif //LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#endif // LIBMESH_FE_MACRO_H
