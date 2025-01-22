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
#include "libmesh/tensor_value.h"


namespace libMesh
{


LIBMESH_DEFAULT_VECTORIZED_FE(0,HIERARCHIC_VEC)
LIBMESH_DEFAULT_VECTORIZED_FE(1,HIERARCHIC_VEC)
LIBMESH_DEFAULT_VECTORIZED_FE(2,HIERARCHIC_VEC)
LIBMESH_DEFAULT_VECTORIZED_FE(3,HIERARCHIC_VEC)
LIBMESH_DEFAULT_VECTORIZED_FE(0,L2_HIERARCHIC_VEC)
LIBMESH_DEFAULT_VECTORIZED_FE(1,L2_HIERARCHIC_VEC)
LIBMESH_DEFAULT_VECTORIZED_FE(2,L2_HIERARCHIC_VEC)
LIBMESH_DEFAULT_VECTORIZED_FE(3,L2_HIERARCHIC_VEC)


// ------------------------------------------------------------
// Hierarchic-specific implementations


// Anonymous namespace for local helper functions
namespace {
void hierarchic_vec_nodal_soln(const Elem * elem,
                               const Order order,
                               const std::vector<Number> & elem_soln,
                               const int dim,
                               std::vector<Number> &       nodal_soln,
                               const bool add_p_level)
{
  const unsigned int n_nodes = elem->n_nodes();

  // Constant shape functions can't be supported, even for
  // L2_HIERARCHIC*, without breaking the "HIERARCHIC is
  // hierarchic" guarantee
  libmesh_assert(order != CONSTANT);

  nodal_soln.resize(dim*n_nodes);

  // Do interpolation at the nodes explicitly.
  FEType fe_type(order, HIERARCHIC);

  const unsigned int n_sf =
    FEInterface::n_shape_functions(fe_type, elem, add_p_level);

  std::vector<Point> refspace_nodes;
  FEBase::get_refspace_nodes(elem->type(),refspace_nodes);
  libmesh_assert_equal_to (refspace_nodes.size(), n_nodes);
  libmesh_assert_equal_to (elem_soln.size(), n_sf*dim);

  for (unsigned int n=0; n<n_nodes; n++)
    for (int d=0; d != dim; ++d)
      {
        const unsigned int ni=n*dim+d;
        nodal_soln[ni] = 0;

        // u_i = Sum (alpha_i phi_i); we're here only looking
        // at vector components in direction d
        for (unsigned int i=0; i<n_sf; i++)
          nodal_soln[ni] += elem_soln[i*dim+d] *
            FEInterface::shape(fe_type, elem, i, refspace_nodes[n]);
      }

}// void hierarchic_vec_nodal_soln

} // anonymous namespace


  // Do full-specialization for every dimension, instead
  // of explicit instantiation at the end of this file.
  // This could be macro-ified so that it fits on one line...
template <>
void FE<0,HIERARCHIC_VEC>::nodal_soln(const Elem * elem,
                                      const Order order,
                                      const std::vector<Number> & elem_soln,
                                      std::vector<Number> & nodal_soln,
                                      const bool add_p_level)
{ FE<0,HIERARCHIC>::nodal_soln(elem, order, elem_soln, nodal_soln, add_p_level); }

template <>
void FE<1,HIERARCHIC_VEC>::nodal_soln(const Elem * elem,
                                      const Order order,
                                      const std::vector<Number> & elem_soln,
                                      std::vector<Number> & nodal_soln,
                                      const bool add_p_level)
{ FE<1,HIERARCHIC>::nodal_soln(elem, order, elem_soln, nodal_soln, add_p_level); }

template <>
void FE<2,HIERARCHIC_VEC>::nodal_soln(const Elem * elem,
                                      const Order order,
                                      const std::vector<Number> & elem_soln,
                                      std::vector<Number> & nodal_soln,
                                      const bool add_p_level)
{ hierarchic_vec_nodal_soln(elem, order, elem_soln, 2 /*dimension*/, nodal_soln, add_p_level); }

template <>
void FE<3,HIERARCHIC_VEC>::nodal_soln(const Elem * elem,
                                      const Order order,
                                      const std::vector<Number> & elem_soln,
                                      std::vector<Number> & nodal_soln,
                                      const bool add_p_level)
{ hierarchic_vec_nodal_soln(elem, order, elem_soln, 3 /*dimension*/, nodal_soln, add_p_level); }

LIBMESH_FE_SIDE_NODAL_SOLN(HIERARCHIC_VEC)

template <>
void FE<0,L2_HIERARCHIC_VEC>::nodal_soln(const Elem * elem,
                                         const Order order,
                                         const std::vector<Number> & elem_soln,
                                         std::vector<Number> & nodal_soln,
                                         const bool add_p_level)
{ FE<0,HIERARCHIC_VEC>::nodal_soln(elem, order, elem_soln, nodal_soln, add_p_level); }

template <>
void FE<1,L2_HIERARCHIC_VEC>::nodal_soln(const Elem * elem,
                                         const Order order,
                                         const std::vector<Number> & elem_soln,
                                         std::vector<Number> & nodal_soln,
                                         const bool add_p_level)
{ FE<1,HIERARCHIC_VEC>::nodal_soln(elem, order, elem_soln, nodal_soln, add_p_level); }

template <>
void FE<2,L2_HIERARCHIC_VEC>::nodal_soln(const Elem * elem,
                                         const Order order,
                                         const std::vector<Number> & elem_soln,
                                         std::vector<Number> & nodal_soln,
                                         const bool add_p_level)
{ FE<2,HIERARCHIC_VEC>::nodal_soln(elem, order, elem_soln, nodal_soln, add_p_level); }

template <>
void FE<3,L2_HIERARCHIC_VEC>::nodal_soln(const Elem * elem,
                                         const Order order,
                                         const std::vector<Number> & elem_soln,
                                         std::vector<Number> & nodal_soln,
                                         const bool add_p_level)
{ FE<3,HIERARCHIC_VEC>::nodal_soln(elem, order, elem_soln, nodal_soln, add_p_level); }

LIBMESH_FE_SIDE_NODAL_SOLN(L2_HIERARCHIC_VEC)


// Specialize for shape function routines by leveraging scalar HIERARCHIC elements

// 0-D
template <> RealGradient FE<0,HIERARCHIC_VEC>::shape(const ElemType type, const Order order,
                                                     const unsigned int i, const Point & p)
{
  Real value = FE<0,HIERARCHIC>::shape( type, order, i, p );
  return libMesh::RealGradient( value );
}
template <> RealGradient FE<0,HIERARCHIC_VEC>::shape_deriv(const ElemType type, const Order order,
                                                           const unsigned int i, const unsigned int j,
                                                           const Point & p)
{
  Real value = FE<0,HIERARCHIC>::shape_deriv( type, order, i, j, p );
  return libMesh::RealGradient( value );
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <> RealGradient FE<0,HIERARCHIC_VEC>::shape_second_deriv(const ElemType type, const Order order,
                                                                  const unsigned int i, const unsigned int j,
                                                                  const Point & p)
{
  Real value = FE<0,HIERARCHIC>::shape_second_deriv( type, order, i, j, p );
  return libMesh::RealGradient( value );
}
#endif

template <> RealGradient FE<0,L2_HIERARCHIC_VEC>::shape(const ElemType type, const Order order,
                                                        const unsigned int i, const Point & p)
{
  return FE<0,HIERARCHIC_VEC>::shape(type, order, i, p);
}
template <> RealGradient FE<0,L2_HIERARCHIC_VEC>::shape_deriv(const ElemType type, const Order order,
                                                              const unsigned int i, const unsigned int j,
                                                              const Point & p)
{
  return FE<0,HIERARCHIC_VEC>::shape_deriv(type, order, i, j, p);
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <> RealGradient FE<0,L2_HIERARCHIC_VEC>::shape_second_deriv(const ElemType type, const Order order,
                                                                     const unsigned int i, const unsigned int j,
                                                                     const Point & p)
{
  return FE<0,HIERARCHIC_VEC>::shape_second_deriv(type, order, i, j, p);
}
#endif

// 1-D
template <> RealGradient FE<1,HIERARCHIC_VEC>::shape(const ElemType type, const Order order,
                                                     const unsigned int i, const Point & p)
{
  Real value = FE<1,HIERARCHIC>::shape( type, order, i, p );
  return libMesh::RealGradient( value );
}
template <> RealGradient FE<1,HIERARCHIC_VEC>::shape_deriv(const ElemType type, const Order order,
                                                           const unsigned int i, const unsigned int j,
                                                           const Point & p)
{
  Real value = FE<1,HIERARCHIC>::shape_deriv( type, order, i, j, p );
  return libMesh::RealGradient( value );
}
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <> RealGradient FE<1,HIERARCHIC_VEC>::shape_second_deriv(const ElemType type, const Order order,
                                                                  const unsigned int i, const unsigned int j,
                                                                  const Point & p)
{
  Real value = FE<1,HIERARCHIC>::shape_second_deriv( type, order, i, j, p );
  return libMesh::RealGradient( value );
}

#endif

template <> RealGradient FE<1,L2_HIERARCHIC_VEC>::shape(const ElemType type, const Order order,
                                                        const unsigned int i, const Point & p)
{
  return FE<1,HIERARCHIC_VEC>::shape(type, order, i, p);
}
template <> RealGradient FE<1,L2_HIERARCHIC_VEC>::shape_deriv(const ElemType type, const Order order,
                                                              const unsigned int i, const unsigned int j,
                                                              const Point & p)
{
  return FE<1,HIERARCHIC_VEC>::shape_deriv(type, order, i, j, p);
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <> RealGradient FE<1,L2_HIERARCHIC_VEC>::shape_second_deriv(const ElemType type, const Order order,
                                                                     const unsigned int i, const unsigned int j,
                                                                     const Point & p)
{
  return FE<1,HIERARCHIC_VEC>::shape_second_deriv(type, order, i, j, p);
}
#endif

// 2-D
template <> RealGradient FE<2,HIERARCHIC_VEC>::shape(const ElemType type, const Order order,
                                                     const unsigned int i, const Point & p)
{
  Real value = FE<2,HIERARCHIC>::shape( type, order, i/2, p );

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
template <> RealGradient FE<2,HIERARCHIC_VEC>::shape_deriv(const ElemType type, const Order order,
                                                           const unsigned int i, const unsigned int j,
                                                           const Point & p)
{
  Real value = FE<2,HIERARCHIC>::shape_deriv( type, order, i/2, j, p );

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

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <> RealGradient FE<2,HIERARCHIC_VEC>::shape_second_deriv(const ElemType type, const Order order,
                                                                  const unsigned int i, const unsigned int j,
                                                                  const Point & p)
{
  Real value = FE<2,HIERARCHIC>::shape_second_deriv( type, order, i/2, j, p );

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

#endif

template <> RealGradient FE<2,L2_HIERARCHIC_VEC>::shape(const ElemType type, const Order order,
                                                        const unsigned int i, const Point & p)
{
  return FE<2,HIERARCHIC_VEC>::shape(type, order, i, p);
}

template <> RealGradient FE<2,L2_HIERARCHIC_VEC>::shape_deriv(const ElemType type, const Order order,
                                                              const unsigned int i, const unsigned int j,
                                                              const Point & p)
{
  return FE<2,HIERARCHIC_VEC>::shape_deriv(type, order, i, j, p);
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <> RealGradient FE<2,L2_HIERARCHIC_VEC>::shape_second_deriv(const ElemType type, const Order order,
                                                                     const unsigned int i, const unsigned int j,
                                                                     const Point & p)
{
  return FE<2,HIERARCHIC_VEC>::shape_second_deriv(type, order, i, j, p);
}

#endif

// 3-D
template <> RealGradient FE<3,HIERARCHIC_VEC>::shape(const ElemType type, const Order order,
                                                     const unsigned int i, const Point & p)
{
  Real value = FE<3,HIERARCHIC>::shape( type, order, i/3, p );

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
template <> RealGradient FE<3,HIERARCHIC_VEC>::shape_deriv(const ElemType type, const Order order,
                                                           const unsigned int i, const unsigned int j,
                                                           const Point & p)
{
  Real value = FE<3,HIERARCHIC>::shape_deriv( type, order, i/3, j, p );

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

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <> RealGradient FE<3,HIERARCHIC_VEC>::shape_second_deriv(const ElemType type, const Order order,
                                                                  const unsigned int i, const unsigned int j,
                                                                  const Point & p)
{
  Real value = FE<3,HIERARCHIC>::shape_second_deriv( type, order, i/3, j, p );

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

#endif

template <> RealGradient FE<3,L2_HIERARCHIC_VEC>::shape(const ElemType type, const Order order,
                                                        const unsigned int i, const Point & p)
{
  return FE<3,HIERARCHIC_VEC>::shape(type, order, i, p);
}

template <> RealGradient FE<3,L2_HIERARCHIC_VEC>::shape_deriv(const ElemType type, const Order order,
                                                              const unsigned int i, const unsigned int j,
                                                              const Point & p)
{
  return FE<3,HIERARCHIC_VEC>::shape_deriv(type, order, i, j, p);
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <> RealGradient FE<3,L2_HIERARCHIC_VEC>::shape_second_deriv(const ElemType type, const Order order,
                                                                     const unsigned int i, const unsigned int j,
                                                                     const Point & p)
{
  return FE<3,HIERARCHIC_VEC>::shape_second_deriv(type, order, i, j, p);
}

#endif


// 0-D
template <> RealGradient FE<0,HIERARCHIC_VEC>::shape(const Elem * elem, const Order order,
                                                     const unsigned int i, const Point & p,
                                                     const bool add_p_level)
{
  const Real value = FE<0,HIERARCHIC>::shape(elem, order, i, p, add_p_level);
  return libMesh::RealGradient( value );
}
template <> RealGradient FE<0,HIERARCHIC_VEC>::shape_deriv(const Elem * elem, const Order order,
                                                           const unsigned int i, const unsigned int j,
                                                           const Point & p,
                                                           const bool add_p_level)
{
  const Real value = FE<0,HIERARCHIC>::shape_deriv(elem, order, i, j, p, add_p_level);
  return libMesh::RealGradient( value );
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <> RealGradient FE<0,HIERARCHIC_VEC>::shape_second_deriv(const Elem * elem, const Order order,
                                                                  const unsigned int i, const unsigned int j,
                                                                  const Point & p,
                                                                  const bool add_p_level)
{
  Real value = FE<0,HIERARCHIC>::shape_second_deriv(elem, order, i, j, p, add_p_level);
  return libMesh::RealGradient( value );
}

#endif

template <> RealGradient FE<0,L2_HIERARCHIC_VEC>::shape(const Elem * elem, const Order order,
                                                        const unsigned int i, const Point & p,
                                                        const bool add_p_level)
{
  return FE<0,HIERARCHIC_VEC>::shape(elem, order, i, p, add_p_level);
}

template <> RealGradient FE<0,L2_HIERARCHIC_VEC>::shape_deriv(const Elem * elem, const Order order,
                                                              const unsigned int i, const unsigned int j,
                                                              const Point & p,
                                                              const bool add_p_level)
{
  return FE<0,HIERARCHIC_VEC>::shape_deriv(elem, order, i, j, p, add_p_level);
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <> RealGradient FE<0,L2_HIERARCHIC_VEC>::shape_second_deriv(const Elem * elem, const Order order,
                                                                     const unsigned int i, const unsigned int j,
                                                                     const Point & p,
                                                                     const bool add_p_level)
{
  return FE<0,HIERARCHIC_VEC>::shape_second_deriv(elem, order, i, j, p, add_p_level);
}

#endif

// 1-D
template <> RealGradient FE<1,HIERARCHIC_VEC>::shape(const Elem * elem, const Order order,
                                                     const unsigned int i, const Point & p,
                                                     const bool add_p_level)
{
  Real value = FE<1,HIERARCHIC>::shape(elem, order, i, p, add_p_level);
  return libMesh::RealGradient( value );
}
template <> RealGradient FE<1,HIERARCHIC_VEC>::shape_deriv(const Elem * elem, const Order order,
                                                           const unsigned int i, const unsigned int j,
                                                           const Point & p,
                                                           const bool add_p_level)
{
  Real value = FE<1,HIERARCHIC>::shape_deriv(elem, order, i, j, p, add_p_level);
  return libMesh::RealGradient( value );
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <> RealGradient FE<1,HIERARCHIC_VEC>::shape_second_deriv(const Elem * elem, const Order order,
                                                                  const unsigned int i, const unsigned int j,
                                                                  const Point & p,
                                                                  const bool add_p_level)
{
  Real value = FE<1,HIERARCHIC>::shape_second_deriv(elem, order, i, j, p, add_p_level);
  return libMesh::RealGradient( value );
}

#endif

template <> RealGradient FE<1,L2_HIERARCHIC_VEC>::shape(const Elem * elem, const Order order,
                                                        const unsigned int i, const Point & p,
                                                        const bool add_p_level)
{
  return FE<1,HIERARCHIC_VEC>::shape(elem, order, i, p, add_p_level);
}

template <> RealGradient FE<1,L2_HIERARCHIC_VEC>::shape_deriv(const Elem * elem, const Order order,
                                                              const unsigned int i, const unsigned int j,
                                                              const Point & p,
                                                              const bool add_p_level)
{
  return FE<1,HIERARCHIC_VEC>::shape_deriv(elem, order, i, j, p, add_p_level);
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <> RealGradient FE<1,L2_HIERARCHIC_VEC>::shape_second_deriv(const Elem * elem, const Order order,
                                                                     const unsigned int i, const unsigned int j,
                                                                     const Point & p,
                                                                     const bool add_p_level)
{
  return FE<1,HIERARCHIC_VEC>::shape_second_deriv(elem, order, i, j, p, add_p_level);
}

#endif

// 2-D
template <> RealGradient FE<2,HIERARCHIC_VEC>::shape(const Elem * elem, const Order order,
                                                     const unsigned int i, const Point & p,
                                                     const bool add_p_level)
{
  const Real value = FE<2,HIERARCHIC>::shape(elem, order, i/2, p, add_p_level);

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
template <> RealGradient FE<2,HIERARCHIC_VEC>::shape_deriv(const Elem * elem, const Order order,
                                                           const unsigned int i, const unsigned int j,
                                                           const Point & p,
                                                           const bool add_p_level)
{
  const Real value = FE<2,HIERARCHIC>::shape_deriv(elem, order, i/2, j, p, add_p_level);

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

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
template <> RealGradient FE<2,HIERARCHIC_VEC>::shape_second_deriv(const Elem * elem, const Order order,
                                                                  const unsigned int i, const unsigned int j,
                                                                  const Point & p,
                                                                  const bool add_p_level)
{
  const Real value = FE<2,HIERARCHIC>::shape_second_deriv(elem, order, i/2, j, p, add_p_level);

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

#endif

template <> RealGradient FE<2,L2_HIERARCHIC_VEC>::shape(const Elem * elem, const Order order,
                                                        const unsigned int i, const Point & p,
                                                        const bool add_p_level)
{
  return FE<2,HIERARCHIC_VEC>::shape(elem, order, i, p, add_p_level);
}
template <> RealGradient FE<2,L2_HIERARCHIC_VEC>::shape_deriv(const Elem * elem, const Order order,
                                                              const unsigned int i, const unsigned int j,
                                                              const Point & p,
                                                              const bool add_p_level)
{
  return FE<2,HIERARCHIC_VEC>::shape_deriv(elem, order, i, j, p, add_p_level);
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <> RealGradient FE<2,L2_HIERARCHIC_VEC>::shape_second_deriv(const Elem * elem, const Order order,
                                                                     const unsigned int i, const unsigned int j,
                                                                     const Point & p,
                                                                     const bool add_p_level)
{
  return FE<2,HIERARCHIC_VEC>::shape_second_deriv(elem, order, i, j, p, add_p_level);
}

#endif

// 3-D
template <> RealGradient FE<3,HIERARCHIC_VEC>::shape(const Elem * elem, const Order order,
                                                     const unsigned int i, const Point & p,
                                                     const bool add_p_level)
{
  const Real value = FE<3,HIERARCHIC>::shape(elem, order, i/3, p, add_p_level);

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

template <> RealGradient FE<3,HIERARCHIC_VEC>::shape_deriv(const Elem * elem, const Order order,
                                                           const unsigned int i, const unsigned int j,
                                                           const Point & p,
                                                           const bool add_p_level)
{
  const Real value = FE<3,HIERARCHIC>::shape_deriv(elem, order, i/3, j, p, add_p_level);

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

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <> RealGradient FE<3,HIERARCHIC_VEC>::shape_second_deriv(const Elem * elem, const Order order,
                                                                  const unsigned int i, const unsigned int j,
                                                                  const Point & p,
                                                                  const bool add_p_level)
{
  const Real value = FE<3,HIERARCHIC>::shape_second_deriv(elem, order, i/3, j, p, add_p_level);

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

#endif

template <> RealGradient FE<3,L2_HIERARCHIC_VEC>::shape(const Elem * elem, const Order order,
                                                        const unsigned int i, const Point & p,
                                                        const bool add_p_level)
{
  return FE<3,HIERARCHIC_VEC>::shape(elem, order, i, p, add_p_level);
}

template <> RealGradient FE<3,L2_HIERARCHIC_VEC>::shape_deriv(const Elem * elem, const Order order,
                                                              const unsigned int i, const unsigned int j,
                                                              const Point & p,
                                                              const bool add_p_level)
{
  return FE<3,HIERARCHIC_VEC>::shape_deriv(elem, order, i, j, p, add_p_level);
}

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES

template <> RealGradient FE<3,L2_HIERARCHIC_VEC>::shape_second_deriv(const Elem * elem, const Order order,
                                                                   const unsigned int i, const unsigned int j,
                                                                   const Point & p,
                                                                   const bool add_p_level)
{
  return FE<3,HIERARCHIC_VEC>::shape_second_deriv(elem, order, i, j, p, add_p_level);
}

#endif

// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this function.
// This could be macro-ified.
template <> unsigned int FE<0,HIERARCHIC_VEC>::n_dofs(const ElemType t, const Order o) { return FE<0,HIERARCHIC>::n_dofs(t,o); }
template <> unsigned int FE<1,HIERARCHIC_VEC>::n_dofs(const ElemType t, const Order o) { return FE<1,HIERARCHIC>::n_dofs(t,o); }
template <> unsigned int FE<2,HIERARCHIC_VEC>::n_dofs(const ElemType t, const Order o) { return 2*FE<2,HIERARCHIC>::n_dofs(t,o); }
template <> unsigned int FE<3,HIERARCHIC_VEC>::n_dofs(const ElemType t, const Order o) { return 3*FE<3,HIERARCHIC>::n_dofs(t,o); }

template <> unsigned int FE<0,L2_HIERARCHIC_VEC>::n_dofs(const ElemType t, const Order o) { return FE<0,L2_HIERARCHIC>::n_dofs(t,o); }
template <> unsigned int FE<1,L2_HIERARCHIC_VEC>::n_dofs(const ElemType t, const Order o) { return FE<1,L2_HIERARCHIC>::n_dofs(t,o); }
template <> unsigned int FE<2,L2_HIERARCHIC_VEC>::n_dofs(const ElemType t, const Order o) { return 2*FE<2,L2_HIERARCHIC>::n_dofs(t,o); }
template <> unsigned int FE<3,L2_HIERARCHIC_VEC>::n_dofs(const ElemType t, const Order o) { return 3*FE<3,L2_HIERARCHIC>::n_dofs(t,o); }


// Do full-specialization for every dimension, instead
// of explicit instantiation at the end of this function.
template <> unsigned int FE<0,HIERARCHIC_VEC>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return FE<0,HIERARCHIC>::n_dofs_at_node(t,o,n); }
template <> unsigned int FE<1,HIERARCHIC_VEC>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return FE<1,HIERARCHIC>::n_dofs_at_node(t,o,n); }
template <> unsigned int FE<2,HIERARCHIC_VEC>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return 2*FE<2,HIERARCHIC>::n_dofs_at_node(t,o,n); }
template <> unsigned int FE<3,HIERARCHIC_VEC>::n_dofs_at_node(const ElemType t, const Order o, const unsigned int n) { return 3*FE<3,HIERARCHIC>::n_dofs_at_node(t,o,n); }

template <> unsigned int FE<0,L2_HIERARCHIC_VEC>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<1,L2_HIERARCHIC_VEC>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<2,L2_HIERARCHIC_VEC>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }
template <> unsigned int FE<3,L2_HIERARCHIC_VEC>::n_dofs_at_node(const ElemType, const Order, const unsigned int) { return 0; }

template <> unsigned int FE<0,HIERARCHIC_VEC>::n_dofs_per_elem(const ElemType t, const Order o) { return FE<0,HIERARCHIC>::n_dofs_per_elem(t,o); }
template <> unsigned int FE<1,HIERARCHIC_VEC>::n_dofs_per_elem(const ElemType t, const Order o) { return FE<1,HIERARCHIC>::n_dofs_per_elem(t,o); }
template <> unsigned int FE<2,HIERARCHIC_VEC>::n_dofs_per_elem(const ElemType t, const Order o) { return 2*FE<2,HIERARCHIC>::n_dofs_per_elem(t,o); }
template <> unsigned int FE<3,HIERARCHIC_VEC>::n_dofs_per_elem(const ElemType t, const Order o) { return 3*FE<3,HIERARCHIC>::n_dofs_per_elem(t,o); }

// L2 Hierarchic elements have all their dofs on the element
template <> unsigned int FE<0,L2_HIERARCHIC_VEC>::n_dofs_per_elem(const ElemType t, const Order o) { return FE<0,L2_HIERARCHIC_VEC>::n_dofs(t, o); }
template <> unsigned int FE<1,L2_HIERARCHIC_VEC>::n_dofs_per_elem(const ElemType t, const Order o) { return FE<1,L2_HIERARCHIC_VEC>::n_dofs(t, o); }
template <> unsigned int FE<2,L2_HIERARCHIC_VEC>::n_dofs_per_elem(const ElemType t, const Order o) { return FE<2,L2_HIERARCHIC_VEC>::n_dofs(t, o); }
template <> unsigned int FE<3,L2_HIERARCHIC_VEC>::n_dofs_per_elem(const ElemType t, const Order o) { return FE<3,L2_HIERARCHIC_VEC>::n_dofs(t, o); }

// Hierarchic FEMs are always C^0 continuous
template <> FEContinuity FE<0,HIERARCHIC_VEC>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<1,HIERARCHIC_VEC>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<2,HIERARCHIC_VEC>::get_continuity() const { return C_ZERO; }
template <> FEContinuity FE<3,HIERARCHIC_VEC>::get_continuity() const { return C_ZERO; }

// L2 Hierarchic FEMs are always discontinuous
template <> FEContinuity FE<0,L2_HIERARCHIC_VEC>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<1,L2_HIERARCHIC_VEC>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<2,L2_HIERARCHIC_VEC>::get_continuity() const { return DISCONTINUOUS; }
template <> FEContinuity FE<3,L2_HIERARCHIC_VEC>::get_continuity() const { return DISCONTINUOUS; }

// Hierarchic FEMs are hierarchic
template <> bool FE<0,HIERARCHIC_VEC>::is_hierarchic() const { return true; }
template <> bool FE<1,HIERARCHIC_VEC>::is_hierarchic() const { return true; }
template <> bool FE<2,HIERARCHIC_VEC>::is_hierarchic() const { return true; }
template <> bool FE<3,HIERARCHIC_VEC>::is_hierarchic() const { return true; }
template <> bool FE<0,L2_HIERARCHIC_VEC>::is_hierarchic() const { return true; }
template <> bool FE<1,L2_HIERARCHIC_VEC>::is_hierarchic() const { return true; }
template <> bool FE<2,L2_HIERARCHIC_VEC>::is_hierarchic() const { return true; }
template <> bool FE<3,L2_HIERARCHIC_VEC>::is_hierarchic() const { return true; }

// Hierarchic FEM shapes need reinit
template <> bool FE<0,HIERARCHIC_VEC>::shapes_need_reinit() const { return true; }
template <> bool FE<1,HIERARCHIC_VEC>::shapes_need_reinit() const { return true; }
template <> bool FE<2,HIERARCHIC_VEC>::shapes_need_reinit() const { return true; }
template <> bool FE<3,HIERARCHIC_VEC>::shapes_need_reinit() const { return true; }
template <> bool FE<0,L2_HIERARCHIC_VEC>::shapes_need_reinit() const { return true; }
template <> bool FE<1,L2_HIERARCHIC_VEC>::shapes_need_reinit() const { return true; }
template <> bool FE<2,L2_HIERARCHIC_VEC>::shapes_need_reinit() const { return true; }
template <> bool FE<3,L2_HIERARCHIC_VEC>::shapes_need_reinit() const { return true; }

// Methods for computing Hierarchic constraints.  Note: we pass the
// dimension as the last argument to the anonymous helper function.
// Also note: we only need instantiations of this function for
// Dim==2 and 3.
#ifdef LIBMESH_ENABLE_AMR
template <>
void FE<2,HIERARCHIC_VEC>::compute_constraints (DofConstraints & constraints,
                                                DofMap & dof_map,
                                                const unsigned int variable_number,
                                                const Elem * elem)
{ //libmesh_not_implemented();
  FEVectorBase::compute_proj_constraints(constraints, dof_map, variable_number, elem);
}

template <>
void FE<3,HIERARCHIC_VEC>::compute_constraints (DofConstraints & constraints,
                                                DofMap & dof_map,
                                                const unsigned int variable_number,
                                                const Elem * elem)
{ //libmesh_not_implemented();
  FEVectorBase::compute_proj_constraints(constraints, dof_map, variable_number, elem);
}

template <>
void FE<2,L2_HIERARCHIC_VEC>::compute_constraints (DofConstraints & constraints,
                                                   DofMap & dof_map,
                                                   const unsigned int variable_number,
                                                   const Elem * elem)
{ //libmesh_not_implemented();
  FEVectorBase::compute_proj_constraints(constraints, dof_map, variable_number, elem);
}

template <>
void FE<3,L2_HIERARCHIC_VEC>::compute_constraints (DofConstraints & constraints,
                                                   DofMap & dof_map,
                                                   const unsigned int variable_number,
                                                   const Elem * elem)
{ //libmesh_not_implemented();
  FEVectorBase::compute_proj_constraints(constraints, dof_map, variable_number, elem);
}
#endif // LIBMESH_ENABLE_AMR

} // namespace libMesh
