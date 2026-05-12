// Kokkos FEEvaluator specializations for MONOMIAL elements.
//
// MONOMIAL uses the complete total-degree polynomial space P_p. Following
// libMesh's FE<Dim, MONOMIAL>, the basis is parameterized by spatial
// dimension, not element class: TRI and QUAD share the 2-D implementation,
// while TET/HEX/PRISM/PYRAMID share the 3-D implementation.
//
// The implementation below mirrors the generic index-to-exponent decoding used
// by libMesh's host-side MONOMIAL FE code, so the Kokkos layer reuses the same
// basis ordering without hand-expanding every order into bespoke tables.

#ifndef LIBMESH_KOKKOS_FE_MONOMIAL_H
#define LIBMESH_KOKKOS_FE_MONOMIAL_H

#include "kokkos_fe_base.h"
#include "libmesh/enum_elem_type.h"

namespace libMesh::Kokkos
{

namespace detail
{

LIBMESH_DEVICE_INLINE Real
pow_unsigned(Real base, unsigned int exponent)
{
  Real value = 1;
  for (unsigned int i = 0; i < exponent; ++i)
    value *= base;
  return value;
}

template <unsigned int Dim>
struct monomial_exponents;

template <>
struct monomial_exponents<1>
{
  unsigned int nx;

  LIBMESH_DEVICE_INLINE static monomial_exponents decode(unsigned int i)
  {
    return {i};
  }
};

template <>
struct monomial_exponents<2>
{
  unsigned int nx;
  unsigned int ny;

  LIBMESH_DEVICE_INLINE static monomial_exponents decode(unsigned int i)
  {
    unsigned int degree = 0;
    for (; i >= (degree + 1) * (degree + 2) / 2; ++degree) {}

    const unsigned int ny = i - (degree * (degree + 1) / 2);
    const unsigned int nx = degree - ny;
    return {nx, ny};
  }
};

template <>
struct monomial_exponents<3>
{
  unsigned int nx;
  unsigned int ny;
  unsigned int nz;

  LIBMESH_DEVICE_INLINE static monomial_exponents decode(unsigned int i)
  {
    unsigned int degree = 0;
    for (; i >= (degree + 1) * (degree + 2) * (degree + 3) / 6; ++degree) {}

    const unsigned int degree_offset = degree * (degree + 1) * (degree + 2) / 6;
    const unsigned int local_index = i - degree_offset;

    unsigned int block = degree;
    unsigned int nz = 0;
    for (; block < local_index; block += (degree - nz + 1))
      ++nz;

    const unsigned int nx = block - local_index;
    const unsigned int ny = degree - nx - nz;
    return {nx, ny, nz};
  }
};

} // namespace detail

template <unsigned int N>
struct MonomialImpl1D
{
  static constexpr unsigned int n_dofs() { return N + 1; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real /*eta*/, Real /*zeta*/)
  {
    return detail::pow_unsigned(xi, i);
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real /*eta*/, Real /*zeta*/)
  {
    if (!i)
      return zero_vector();

    return make_vector(i * detail::pow_unsigned(xi, i - 1), 0.0, 0.0);
  }
};

template <unsigned int N>
struct MonomialImpl2D
{
  static constexpr unsigned int n_dofs() { return (N + 1) * (N + 2) / 2; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    const auto exponents = detail::monomial_exponents<2>::decode(i);
    return detail::pow_unsigned(xi, exponents.nx) *
           detail::pow_unsigned(eta, exponents.ny);
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real /*zeta*/)
  {
    const auto exponents = detail::monomial_exponents<2>::decode(i);
    const Real dx = exponents.nx
                      ? exponents.nx *
                          detail::pow_unsigned(xi, exponents.nx - 1) *
                          detail::pow_unsigned(eta, exponents.ny)
                      : 0.0;
    const Real dy = exponents.ny
                      ? exponents.ny *
                          detail::pow_unsigned(xi, exponents.nx) *
                          detail::pow_unsigned(eta, exponents.ny - 1)
                      : 0.0;
    return make_vector(dx, dy, 0.0);
  }
};

template <unsigned int N>
struct MonomialImpl3D
{
  static constexpr unsigned int n_dofs() { return (N + 1) * (N + 2) * (N + 3) / 6; }

  LIBMESH_DEVICE_INLINE static Real
  shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    const auto exponents = detail::monomial_exponents<3>::decode(i);
    return detail::pow_unsigned(xi, exponents.nx) *
           detail::pow_unsigned(eta, exponents.ny) *
           detail::pow_unsigned(zeta, exponents.nz);
  }

  LIBMESH_DEVICE_INLINE static RealVector
  grad_shape(unsigned int i, Real xi, Real eta, Real zeta)
  {
    const auto exponents = detail::monomial_exponents<3>::decode(i);
    const Real dx = exponents.nx
                      ? exponents.nx *
                          detail::pow_unsigned(xi, exponents.nx - 1) *
                          detail::pow_unsigned(eta, exponents.ny) *
                          detail::pow_unsigned(zeta, exponents.nz)
                      : 0.0;
    const Real dy = exponents.ny
                      ? exponents.ny *
                          detail::pow_unsigned(xi, exponents.nx) *
                          detail::pow_unsigned(eta, exponents.ny - 1) *
                          detail::pow_unsigned(zeta, exponents.nz)
                      : 0.0;
    const Real dz = exponents.nz
                      ? exponents.nz *
                          detail::pow_unsigned(xi, exponents.nx) *
                          detail::pow_unsigned(eta, exponents.ny) *
                          detail::pow_unsigned(zeta, exponents.nz - 1)
                      : 0.0;
    return make_vector(dx, dy, dz);
  }
};

// Per-topology FEEvaluator delegating specializations

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::EDGE2, N> : MonomialImpl1D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::EDGE3, N> : MonomialImpl1D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::EDGE4, N> : MonomialImpl1D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::TRI3, N> : MonomialImpl2D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::TRI6, N> : MonomialImpl2D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::TRI7, N> : MonomialImpl2D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::QUAD4, N> : MonomialImpl2D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::QUAD8, N> : MonomialImpl2D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::QUAD9, N> : MonomialImpl2D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::TET4, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::TET10, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::TET14, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::HEX8, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::HEX20, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::HEX27, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::PRISM6, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::PRISM15, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::PRISM18, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::PRISM20, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::PRISM21, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::PYRAMID5, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::PYRAMID13, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::PYRAMID14, N> : MonomialImpl3D<N> {};

template <unsigned int N>
struct FEEvaluator<libMesh::MONOMIAL, libMesh::PYRAMID18, N> : MonomialImpl3D<N> {};

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_FE_MONOMIAL_H
