// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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




// These macros help in instantiating specific versions
// of the \p FE class.  Simply include this file, and
// instantiate at the end for the desired dimension(s).
#define INSTANTIATE_MAPS(_dim,_type)                                                    \
  template LIBMESH_EXPORT Point FE<_dim, _type>::map(const Elem *, const Point &);      \
  template LIBMESH_EXPORT Point FE<_dim, _type>::map_xi(const Elem *, const Point &);   \
  template LIBMESH_EXPORT Point FE<_dim, _type>::map_eta(const Elem *, const Point &);  \
  template LIBMESH_EXPORT Point FE<_dim, _type>::map_zeta(const Elem *, const Point &); \
  template LIBMESH_EXPORT void  FE<_dim, _type>::inverse_map(const Elem *, const std::vector<Point> &, std::vector<Point> &, Real, bool); \
  template LIBMESH_EXPORT Point FE<_dim, _type>::inverse_map(const Elem *, const Point &, Real, bool)

#define INSTANTIATE_SUBDIVISION_MAPS                                                      \
  template LIBMESH_EXPORT Point FE<2, SUBDIVISION>::map(const Elem *, const Point &);     \
  template LIBMESH_EXPORT Point FE<2, SUBDIVISION>::map_xi(const Elem *, const Point &);  \
  template LIBMESH_EXPORT Point FE<2, SUBDIVISION>::map_eta(const Elem *, const Point &); \
  template LIBMESH_EXPORT Point FE<2, SUBDIVISION>::map_zeta(const Elem *, const Point &)

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

#define INSTANTIATE_SUBDIVISION_FE                                      \
  template LIBMESH_EXPORT FE<2,SUBDIVISION>::FE(const FEType & fet); \
  template LIBMESH_EXPORT unsigned int FE<2,SUBDIVISION>::n_shape_functions () const;  \
  template LIBMESH_EXPORT void         FE<2,SUBDIVISION>::attach_quadrature_rule (QBase *); \
  template LIBMESH_EXPORT void         FE<2,SUBDIVISION>::reinit(const Elem *,const std::vector<Point> * const,const std::vector<Real> * const); \
  template LIBMESH_EXPORT void         FE<2,SUBDIVISION>::init_base_shape_functions(const std::vector<Point> &, const Elem *); \
  template LIBMESH_EXPORT void         FE<2,SUBDIVISION>::init_shape_functions(const std::vector<Point> &, const Elem *); \
  template LIBMESH_EXPORT void         FE<2,SUBDIVISION>::init_dual_shape_functions(unsigned int, unsigned int); \
  template LIBMESH_EXPORT void         FE<2,SUBDIVISION>::default_all_shape_derivs (const Elem * elem, const Order o, const std::vector<Point> & p, std::vector<std::vector<Real>> * comps[3], const bool add_p_level); \
  template LIBMESH_EXPORT void         FE<2,SUBDIVISION>::default_side_nodal_soln(const Elem * elem, const Order o, const unsigned int side, const std::vector<Number> & elem_soln, std::vector<Number> & nodal_soln_on_side, bool add_p_level)

#else // LIBMESH_ENABLE_INFINITE_ELEMENTS

#define INSTANTIATE_SUBDIVISION_FE                                      \
  template LIBMESH_EXPORT FE<2,SUBDIVISION>::FE(const FEType & fet);                     \
  template LIBMESH_EXPORT unsigned int FE<2,SUBDIVISION>::n_shape_functions () const;  \
  template LIBMESH_EXPORT void         FE<2,SUBDIVISION>::attach_quadrature_rule (QBase *); \
  template LIBMESH_EXPORT void         FE<2,SUBDIVISION>::reinit(const Elem *,const std::vector<Point> * const,const std::vector<Real> * const); \
  template LIBMESH_EXPORT void         FE<2,SUBDIVISION>::init_shape_functions(const std::vector<Point> &, const Elem *); \
  template LIBMESH_EXPORT void         FE<2,SUBDIVISION>::init_dual_shape_functions(unsigned int, unsigned int); \
  template LIBMESH_EXPORT void         FE<2,SUBDIVISION>::default_all_shape_derivs (const Elem * elem, const Order o, const std::vector<Point> & p, std::vector<std::vector<Real>> * comps[3], const bool add_p_level); \
  template LIBMESH_EXPORT void         FE<2,SUBDIVISION>::default_side_nodal_soln(const Elem * elem, const Order o, const unsigned int side, const std::vector<Number> & elem_soln, std::vector<Number> & nodal_soln_on_side, bool add_p_level)

#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS


#ifndef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#define INSTANTIATE_FE(_dim)                                            \
  template class LIBMESH_EXPORT FE< (_dim), CLOUGH>;                    \
  template class LIBMESH_EXPORT FE< (_dim), HERMITE>;                   \
  template class LIBMESH_EXPORT FE< (_dim), HIERARCHIC>;                \
  template class LIBMESH_EXPORT FE< (_dim), HIERARCHIC_VEC>;            \
  template class LIBMESH_EXPORT FE< (_dim), L2_HIERARCHIC>;             \
  template class LIBMESH_EXPORT FE< (_dim), L2_HIERARCHIC_VEC>;         \
  template class LIBMESH_EXPORT FE< (_dim), LAGRANGE>;                  \
  template class LIBMESH_EXPORT FE< (_dim), LAGRANGE_VEC>;              \
  template class LIBMESH_EXPORT FE< (_dim), L2_LAGRANGE>;               \
  template class LIBMESH_EXPORT FE< (_dim), L2_LAGRANGE_VEC>;           \
  template class LIBMESH_EXPORT FE< (_dim), MONOMIAL>;                  \
  template class LIBMESH_EXPORT FE< (_dim), SCALAR>;                    \
  template class LIBMESH_EXPORT FE< (_dim), SIDE_HIERARCHIC>;           \
  template class LIBMESH_EXPORT FE< (_dim), XYZ>;                       \
  template class LIBMESH_EXPORT FE< (_dim), NEDELEC_ONE>;               \
  template class LIBMESH_EXPORT FE< (_dim), MONOMIAL_VEC>;              \
  template class LIBMESH_EXPORT FE< (_dim), RAVIART_THOMAS>;            \
  template class LIBMESH_EXPORT FE< (_dim), L2_RAVIART_THOMAS>

#else //LIBMESH_ENABLE_HIGHER_ORDER_SHAPES

#define INSTANTIATE_FE(_dim)                                            \
  template class LIBMESH_EXPORT FE< (_dim), CLOUGH>;                    \
  template class LIBMESH_EXPORT FE< (_dim), HERMITE>;                   \
  template class LIBMESH_EXPORT FE< (_dim), HIERARCHIC>;                \
  template class LIBMESH_EXPORT FE< (_dim), HIERARCHIC_VEC>;            \
  template class LIBMESH_EXPORT FE< (_dim), L2_HIERARCHIC>;             \
  template class LIBMESH_EXPORT FE< (_dim), L2_HIERARCHIC_VEC>;         \
  template class LIBMESH_EXPORT FE< (_dim), LAGRANGE>;                  \
  template class LIBMESH_EXPORT FE< (_dim), LAGRANGE_VEC>;              \
  template class LIBMESH_EXPORT FE< (_dim), L2_LAGRANGE>;               \
  template class LIBMESH_EXPORT FE< (_dim), L2_LAGRANGE_VEC>;           \
  template class LIBMESH_EXPORT FE< (_dim), MONOMIAL>;                  \
  template class LIBMESH_EXPORT FE< (_dim), SCALAR>;                    \
  template class LIBMESH_EXPORT FE< (_dim), SIDE_HIERARCHIC>;           \
  template class LIBMESH_EXPORT FE< (_dim), BERNSTEIN>;                 \
  template class LIBMESH_EXPORT FE< (_dim), SZABAB>;                    \
  template class LIBMESH_EXPORT FE< (_dim), XYZ>;                       \
  template class LIBMESH_EXPORT FE< (_dim), RATIONAL_BERNSTEIN>;        \
  template class LIBMESH_EXPORT FE< (_dim), NEDELEC_ONE>;               \
  template class LIBMESH_EXPORT FE< (_dim), MONOMIAL_VEC>;              \
  template class LIBMESH_EXPORT FE< (_dim), RAVIART_THOMAS>;            \
  template class LIBMESH_EXPORT FE< (_dim), L2_RAVIART_THOMAS>

#endif //LIBMESH_ENABLE_HIGHER_ORDER_SHAPES


#define LIBMESH_FE_NODAL_SOLN_DIM(_fetype, _funcname, _dim)             \
template <>                                                             \
void FE<_dim,_fetype>::nodal_soln(const Elem * elem,                    \
                                  const Order order,                    \
                                  const std::vector<Number> & elem_soln,\
                                  std::vector<Number> & nodal_soln,     \
                                  const bool add_p_level)               \
{ _funcname(elem, order, elem_soln, nodal_soln, add_p_level); }

#define LIBMESH_FE_NODAL_SOLN(fetype, _funcname)                        \
LIBMESH_FE_NODAL_SOLN_DIM(fetype, _funcname, 0)                         \
LIBMESH_FE_NODAL_SOLN_DIM(fetype, _funcname, 1)                         \
LIBMESH_FE_NODAL_SOLN_DIM(fetype, _funcname, 2)                         \
LIBMESH_FE_NODAL_SOLN_DIM(fetype, _funcname, 3)


#define LIBMESH_FE_SIDE_NODAL_SOLN_DIM(_fetype, _dim)                   \
template <>                                                             \
void FE<_dim,_fetype>::side_nodal_soln(const Elem * elem,               \
                                       const Order order,               \
                                       const unsigned int side,         \
                                       const std::vector<Number> & elem_soln,\
                                       std::vector<Number> & nodal_soln, \
                                       const bool add_p_level)          \
{ default_side_nodal_soln(elem, order, side, elem_soln, nodal_soln, add_p_level); }

#define LIBMESH_FE_SIDE_NODAL_SOLN(fetype)                              \
LIBMESH_FE_SIDE_NODAL_SOLN_DIM(fetype, 0)                               \
LIBMESH_FE_SIDE_NODAL_SOLN_DIM(fetype, 1)                               \
LIBMESH_FE_SIDE_NODAL_SOLN_DIM(fetype, 2)                               \
LIBMESH_FE_SIDE_NODAL_SOLN_DIM(fetype, 3)


#endif // LIBMESH_FE_MACRO_H
