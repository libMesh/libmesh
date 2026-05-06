// Primary FEEvaluator template for Kokkos device-compatible shape functions.
//
// Uses libMesh's own ElemType and FEFamily enums as non-type template
// parameters — no separate tag structs are needed.
//
// All uses must be explicit specializations defined in the kokkos_fe_lagrange_*.h
// and kokkos_fe_monomial.h headers.  Every specialization must provide:
//
//   static constexpr unsigned int n_dofs()
//
//   LIBMESH_DEVICE_INLINE
//   static Real shape(unsigned int i, Real xi, Real eta, Real zeta)
//
//   LIBMESH_DEVICE_INLINE
//   static RealVector grad_shape(unsigned int i, Real xi, Real eta, Real zeta)
//
// Reference-element coordinate conventions (matching libMesh):
//   Edge:  xi  in [-1, 1]
//   Quad:  (xi, eta)          in [-1,1]^2
//   Hex:   (xi, eta, zeta)    in [-1,1]^3
//   Tri:   (xi, eta)          in unit triangle, xi >= 0, eta >= 0, xi+eta <= 1
//   Tet:   (xi, eta, zeta)    in unit tetrahedron
//
// Unused coordinate arguments (e.g. zeta on a 2D element) are accepted but
// ignored, so call sites can always pass all three without special-casing.
//
#ifndef LIBMESH_KOKKOS_FE_BASE_H
#define LIBMESH_KOKKOS_FE_BASE_H

#include "kokkos_scalar_types.h"
#include "libmesh/libmesh_device.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/enum_fe_family.h"

namespace libMesh::Kokkos
{

template <libMesh::FEFamily family, libMesh::ElemType elem_type, unsigned int Order = 0>
struct FEEvaluator; // forward declaration only; instantiation requires a specialization

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_FE_BASE_H
