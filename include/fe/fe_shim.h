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

#ifndef LIBMESH_FE_SHIM_H
#define LIBMESH_FE_SHIM_H

#include "libmesh/libmesh_config.h"
#include "libmesh/enum_fe_family.h"
#include "libmesh/enum_order.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/fe_output_type.h"
#include "libmesh/fe_forward.h"

#include <vector>

namespace libMesh
{
template <typename>
class PointTempl;
template <typename>
class ElemTempl;
template <typename>
class FEMapTempl;
class DofMap;
class DofConstraints;

#define ReinitBody(Dim, Family)                                                                    \
  {                                                                                                \
    typedef PointTempl<RealType> Point;                                                            \
    typedef ElemTempl<RealType> Elem;                                                              \
    typedef FEMapTempl<RealType> FEMap;                                 \
                                                                                                   \
    static void reinit(FE<Dim, Family, RealType> & fe,                                             \
                       const Elem * elem,                                                          \
                       const unsigned int side,                                                    \
                       const Real tolerance = TOLERANCE,                                           \
                       const std::vector<Point> * const pts = nullptr,                             \
                       const std::vector<Real> * const weights = nullptr);                         \
  }

#define EdgeReinitBody(Dim, Family)                                                                \
  {                                                                                                \
    typedef PointTempl<RealType> Point;                                                            \
    typedef ElemTempl<RealType> Elem;                                                              \
    typedef FEMapTempl<RealType> FEMap;                                 \
                                                                                                   \
    static void edge_reinit(FE<Dim, Family, RealType> & fe,                                        \
                            const Elem * elem,                                                     \
                            const unsigned int edge,                                               \
                            const Real tolerance = TOLERANCE,                                      \
                            const std::vector<Point> * const pts = nullptr,                        \
                            const std::vector<Real> * const weights = nullptr);                    \
  }

#define SideMapBody(Dim, Family)                                                                   \
  {                                                                                                \
    typedef PointTempl<RealType> Point;                                                            \
    typedef ElemTempl<RealType> Elem;                                                              \
    typedef FEMapTempl<RealType> FEMap;                                 \
                                                                                                   \
    static void side_map(FE<Dim, Family, RealType> & fe,                                           \
                         const Elem * elem,                                                        \
                         const Elem * side,                                                        \
                         const unsigned int s,                                                     \
                         const std::vector<Point> & reference_side_points,                         \
                         std::vector<Point> & reference_points);                                   \
  }

#define InverseMapBody(Dim, Family)                                                                \
  {                                                                                                \
    typedef FEMapTempl<RealType> FEMap;                                                            \
    typedef PointTempl<RealType> Point;                                                            \
    typedef ElemTempl<RealType> Elem;                                                              \
                                                                                                   \
    static Point inverse_map(const Elem *, const Point &, const Real, const bool);                 \
                                                                                                   \
    static void inverse_map(                                                                       \
        const Elem *, const std::vector<Point> &, std::vector<Point> &, Real, bool);               \
  }

#define DofsOnSideBody(Dim, Family)                                                                \
  {                                                                                                \
    typedef ElemTempl<RealType> Elem;                                                              \
    typedef FEMapTempl<RealType> FEMap;                                 \
                                                                                                   \
    static void dofs_on_side(                                                                      \
        const Elem * const elem, const Order o, unsigned int s, std::vector<unsigned int> & di);   \
    static void dofs_on_edge(                                                                      \
        const Elem * const elem, const Order o, unsigned int e, std::vector<unsigned int> & di);   \
  }

template <unsigned int Dim, FEFamily T, typename RealType>
struct FEReinitShim ReinitBody(Dim, T);

#define PartialSpecReinit(Dim, Family)                                                             \
  template <typename RealType>                                                                     \
  struct FEReinitShim<Dim, Family, RealType> ReinitBody(Dim, Family)

template <unsigned int Dim, FEFamily T, typename RealType>
struct FEEdgeReinitShim EdgeReinitBody(Dim, T);

#define PartialSpecEdgeReinit(Dim, Family)                                                         \
  template <typename RealType>                                                                     \
  struct FEEdgeReinitShim<Dim, Family, RealType> EdgeReinitBody(Dim, Family)

template <unsigned int Dim, FEFamily T, typename RealType>
struct FESideMapShim SideMapBody(Dim, T);

#define PartialSpecSideMap(Dim, Family)                                                            \
  template <typename RealType>                                                                     \
  struct FESideMapShim<Dim, Family, RealType> SideMapBody(Dim, Family)

template <unsigned int Dim, FEFamily T, typename RealType>
struct FEInverseMapShim InverseMapBody(Dim, T);

#define PartialSpecInverseMap(Dim, Family)                                                         \
  template <typename RealType>                                                                     \
  struct FEInverseMapShim<Dim, Family, RealType> InverseMapBody(Dim, Family)

template <unsigned int Dim, FEFamily T, typename RealType>
struct FEDofsOnSideShim DofsOnSideBody(Dim, T);

#define PartialSpecDofsOnSide(Dim, Family)                                                         \
  template <typename RealType>                                                                     \
  struct FEDofsOnSideShim<Dim, Family, RealType> DofsOnSideBody(Dim, Family)

#define ALL_SIDES_0D(Family)                                                                       \
  PartialSpecReinit(0, Family);                                                                    \
  PartialSpecEdgeReinit(0, Family);                                                                \
  PartialSpecSideMap(0, Family)

ALL_SIDES_0D(CLOUGH);
ALL_SIDES_0D(HERMITE);
ALL_SIDES_0D(HIERARCHIC);
ALL_SIDES_0D(L2_HIERARCHIC);
ALL_SIDES_0D(LAGRANGE);
ALL_SIDES_0D(L2_LAGRANGE);
ALL_SIDES_0D(LAGRANGE_VEC);
ALL_SIDES_0D(MONOMIAL);
ALL_SIDES_0D(MONOMIAL_VEC);
ALL_SIDES_0D(NEDELEC_ONE);
ALL_SIDES_0D(SCALAR);
ALL_SIDES_0D(XYZ);
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
ALL_SIDES_0D(BERNSTEIN);
ALL_SIDES_0D(SZABAB);
ALL_SIDES_0D(RATIONAL_BERNSTEIN);
#endif

// 1D partial specs
PartialSpecEdgeReinit(1, CLOUGH);
PartialSpecEdgeReinit(1, HERMITE);
PartialSpecEdgeReinit(1, HIERARCHIC);
PartialSpecEdgeReinit(1, L2_HIERARCHIC);
PartialSpecEdgeReinit(1, LAGRANGE);
PartialSpecEdgeReinit(1, LAGRANGE_VEC);
PartialSpecEdgeReinit(1, L2_LAGRANGE);
PartialSpecReinit(1, XYZ);
PartialSpecEdgeReinit(1, XYZ);
PartialSpecEdgeReinit(1, MONOMIAL);
PartialSpecEdgeReinit(1, MONOMIAL_VEC);
PartialSpecEdgeReinit(1, SCALAR);
PartialSpecReinit(1, NEDELEC_ONE);
PartialSpecEdgeReinit(1, NEDELEC_ONE);
PartialSpecSideMap(1, NEDELEC_ONE);
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
PartialSpecEdgeReinit(1, BERNSTEIN);
PartialSpecEdgeReinit(1, SZABAB);
PartialSpecEdgeReinit(1, RATIONAL_BERNSTEIN);
#endif

// 2D partial specs
PartialSpecSideMap(2, SUBDIVISION);
PartialSpecEdgeReinit(2, SUBDIVISION);
PartialSpecInverseMap(2, SUBDIVISION);
PartialSpecDofsOnSide(2, SUBDIVISION);

template <unsigned int Dim, typename RealType>
struct FEReinitShim<Dim,XYZ,RealType> ReinitBody(Dim, XYZ);

/**
 * This struct is used for shimming when we want to partially specialize FE methods which is illegal
 * in C++.
 */
template <unsigned int Dim, FEFamily T, typename RealType>
struct FEShim;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

#define FEShimClassDefine(Dim, T)                                                                  \
  template <typename RealType>                                                                     \
  struct FEShim<Dim, T, RealType>                                                                  \
  {                                                                                                \
    typedef PointTempl<RealType> Point;                                                            \
    typedef ElemTempl<RealType> Elem;                                                              \
                                                                                                   \
    typedef typename FEOutputType<T, RealType>::type OutputShape;                                  \
                                                                                                   \
    static void nodal_soln(const Elem * elem,                                                      \
                           const Order order,                                                      \
                           const std::vector<Number> & elem_soln,                                  \
                           std::vector<Number> & nodal_soln);                                      \
                                                                                                   \
    static unsigned int n_dofs(const ElemType t, const Order o);                                   \
                                                                                                   \
    static unsigned int n_dofs_at_node(const ElemType t, const Order o, const unsigned int n);     \
                                                                                                   \
    static unsigned int n_dofs_per_elem(const ElemType t, const Order o);                          \
                                                                                                   \
    static FEContinuity get_continuity();                                                          \
                                                                                                   \
    static bool is_hierarchic();                                                                   \
                                                                                                   \
    static void compute_constraints(DofConstraints & constraints,                                  \
                                    DofMap & dof_map,                                              \
                                    const unsigned int variable_number,                            \
                                    const ElemTempl<Real> * elem);      \
                                                                                                   \
    static bool shapes_need_reinit();                                                              \
                                                                                                   \
    static OutputShape shape(const ElemType, const Order, const unsigned int, const Point &);      \
                                                                                                   \
    static OutputShape                                                                             \
    shape(const Elem *, const Order, const unsigned int, const Point &, const bool = true);        \
                                                                                                   \
    static OutputShape shape_deriv(                                                                \
        const ElemType, const Order, const unsigned int, const unsigned int, const Point &);       \
                                                                                                   \
    static OutputShape shape_deriv(const Elem *,                                                   \
                                   const Order,                                                    \
                                   const unsigned int,                                             \
                                   const unsigned int,                                             \
                                   const Point &,                                                  \
                                   const bool = true);                                             \
                                                                                                   \
    static OutputShape shape_second_deriv(                                                         \
        const ElemType, const Order, const unsigned int, const unsigned int, const Point &);       \
                                                                                                   \
    static OutputShape shape_second_deriv(const Elem *,                                            \
                                          const Order,                                             \
                                          const unsigned int,                                      \
                                          const unsigned int,                                      \
                                          const Point &,                                           \
                                          const bool = true);                                      \
  }

#else
#define FEShimClassDefine(Dim, T)                                                                  \
  template <typename RealType>                                                                     \
  struct FEShim<Dim, T, RealType>                                                                  \
  {                                                                                                \
    typedef PointTempl<RealType> Point;                                                            \
    typedef ElemTempl<RealType> Elem;                                                              \
                                                                                                   \
    typedef typename FEOutputType<T, RealType>::type OutputShape;                                  \
                                                                                                   \
    static void nodal_soln(const Elem * elem,                                                      \
                           const Order order,                                                      \
                           const std::vector<Number> & elem_soln,                                  \
                           std::vector<Number> & nodal_soln);                                      \
                                                                                                   \
    static unsigned int n_dofs(const ElemType t, const Order o);                                   \
                                                                                                   \
    static unsigned int n_dofs_at_node(const ElemType t, const Order o, const unsigned int n);     \
                                                                                                   \
    static unsigned int n_dofs_per_elem(const ElemType t, const Order o);                          \
                                                                                                   \
    static FEContinuity get_continuity();                                                          \
                                                                                                   \
    static bool is_hierarchic();                                                                   \
                                                                                                   \
    static void compute_constraints(DofConstraints & constraints,                                  \
                                    DofMap & dof_map,                                              \
                                    const unsigned int variable_number,                            \
                                    const ElemTempl<Real> * elem);      \
                                                                                                   \
    static bool shapes_need_reinit();                                                              \
                                                                                                   \
    static OutputShape shape(const ElemType, const Order, const unsigned int, const Point &);      \
                                                                                                   \
    static OutputShape                                                                             \
    shape(const Elem *, const Order, const unsigned int, const Point &, const bool = true);        \
                                                                                                   \
    static OutputShape shape_deriv(                                                                \
        const ElemType, const Order, const unsigned int, const unsigned int, const Point &);       \
                                                                                                   \
    static OutputShape shape_deriv(const Elem *,                                                   \
                                   const Order,                                                    \
                                   const unsigned int,                                             \
                                   const unsigned int,                                             \
                                   const Point &,                                                  \
                                   const bool = true);                                             \
  }
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

FEShimClassDefine(0, CLOUGH);
FEShimClassDefine(1, CLOUGH);
FEShimClassDefine(2, CLOUGH);
FEShimClassDefine(3, CLOUGH);
FEShimClassDefine(0, HERMITE);
FEShimClassDefine(1, HERMITE);
FEShimClassDefine(2, HERMITE);
FEShimClassDefine(3, HERMITE);
FEShimClassDefine(0, HIERARCHIC);
FEShimClassDefine(1, HIERARCHIC);
FEShimClassDefine(2, HIERARCHIC);
FEShimClassDefine(3, HIERARCHIC);
FEShimClassDefine(0, L2_HIERARCHIC);
FEShimClassDefine(1, L2_HIERARCHIC);
FEShimClassDefine(2, L2_HIERARCHIC);
FEShimClassDefine(3, L2_HIERARCHIC);
FEShimClassDefine(0, LAGRANGE);
FEShimClassDefine(1, LAGRANGE);
FEShimClassDefine(2, LAGRANGE);
FEShimClassDefine(3, LAGRANGE);
FEShimClassDefine(0, L2_LAGRANGE);
FEShimClassDefine(1, L2_LAGRANGE);
FEShimClassDefine(2, L2_LAGRANGE);
FEShimClassDefine(3, L2_LAGRANGE);
FEShimClassDefine(0, LAGRANGE_VEC);
FEShimClassDefine(1, LAGRANGE_VEC);
FEShimClassDefine(2, LAGRANGE_VEC);
FEShimClassDefine(3, LAGRANGE_VEC);
FEShimClassDefine(0, MONOMIAL);
FEShimClassDefine(1, MONOMIAL);
FEShimClassDefine(2, MONOMIAL);
FEShimClassDefine(3, MONOMIAL);
FEShimClassDefine(0, MONOMIAL_VEC);
FEShimClassDefine(1, MONOMIAL_VEC);
FEShimClassDefine(2, MONOMIAL_VEC);
FEShimClassDefine(3, MONOMIAL_VEC);
FEShimClassDefine(0, NEDELEC_ONE);
FEShimClassDefine(1, NEDELEC_ONE);
FEShimClassDefine(2, NEDELEC_ONE);
FEShimClassDefine(3, NEDELEC_ONE);
FEShimClassDefine(0, SCALAR);
FEShimClassDefine(1, SCALAR);
FEShimClassDefine(2, SCALAR);
FEShimClassDefine(3, SCALAR);
FEShimClassDefine(0, XYZ);
FEShimClassDefine(1, XYZ);
FEShimClassDefine(2, XYZ);
FEShimClassDefine(3, XYZ);
FEShimClassDefine(2, SUBDIVISION);
#ifdef LIBMESH_ENABLE_HIGHER_ORDER_SHAPES
FEShimClassDefine(0, BERNSTEIN);
FEShimClassDefine(1, BERNSTEIN);
FEShimClassDefine(2, BERNSTEIN);
FEShimClassDefine(3, BERNSTEIN);
FEShimClassDefine(0, RATIONAL_BERNSTEIN);
FEShimClassDefine(1, RATIONAL_BERNSTEIN);
FEShimClassDefine(2, RATIONAL_BERNSTEIN);
FEShimClassDefine(3, RATIONAL_BERNSTEIN);
FEShimClassDefine(0, SZABAB);
FEShimClassDefine(1, SZABAB);
FEShimClassDefine(2, SZABAB);
FEShimClassDefine(3, SZABAB);
#endif

template <unsigned int Dim, typename RealType>
struct FEHermiteShim;

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
#define FEHermiteShimClassDefine(Dim)                                                              \
  template <typename RealType>                                                                     \
  struct FEHermiteShim<Dim, RealType> : public FEShim<Dim, HERMITE, RealType>                      \
  {                                                                                                \
    /**                                                                                            \
     * 1D hermite functions on unit interval                                                       \
     */                                                                                            \
    static RealType hermite_raw_shape_second_deriv(const unsigned int basis_num,                   \
                                                   const RealType & xi);                           \
    static RealType hermite_raw_shape_deriv(const unsigned int basis_num, const RealType & xi);    \
    static RealType hermite_raw_shape(const unsigned int basis_num, const RealType & xi);          \
  }

#else

#define FEHermiteShimClassDefine(Dim)                                                              \
  template <typename RealType>                                                                     \
  struct FEHermiteShim<Dim, RealType> : public FEShim<Dim, HERMITE, RealType>                      \
  {                                                                                                \
    static RealType hermite_raw_shape_deriv(const unsigned int basis_num, const RealType & xi);    \
    static RealType hermite_raw_shape(const unsigned int basis_num, const RealType & xi);          \
  }

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

FEHermiteShimClassDefine(0);
FEHermiteShimClassDefine(1);
FEHermiteShimClassDefine(2);
FEHermiteShimClassDefine(3);

template <unsigned int Dim, typename RealType>
struct FEMapShapeFuncShim
{
  typedef PointTempl<RealType> Point;
  typedef ElemTempl<RealType> Elem;
  typedef FEMapTempl<RealType> FEMap;

  static void init_face_shape_functions(FEMap & map, const std::vector<Point> &, const Elem *);
  static void init_edge_shape_functions(FEMap & map, const std::vector<Point> &, const Elem *);
};

template <typename RealType>
struct FEMapShapeFuncShim<0, RealType>
{
  typedef PointTempl<RealType> Point;
  typedef ElemTempl<RealType> Elem;
  typedef FEMapTempl<RealType> FEMap;

  static void init_face_shape_functions(FEMap & map, const std::vector<Point> &, const Elem *);
  static void init_edge_shape_functions(FEMap & map, const std::vector<Point> &, const Elem *);
};

} // namespace libMesh

#endif // LIBMESH_FE_SHIM_H
