// Kokkos FE element and family tag types, and the primary FEEvaluator template.
//
// Tag types are zero-size, zero-cost structs used as compile-time template
// parameters to select FEEvaluator specializations at compile time.
//
// The primary FEEvaluator template is declared here; all uses must be explicit
// specializations defined in the fe_lagrange_*.h and fe_monomial.h headers.

#pragma once

#include "libmesh/kokkos/scalar_types.h"

namespace libMesh::Kokkos
{

// ── Element topology tags ─────────────────────────────────────────────────────

struct Edge2Tag  {}; // 1D, 2-node linear edge
struct Edge3Tag  {}; // 1D, 3-node quadratic edge
struct Tri3Tag   {}; // 2D, 3-node linear triangle
struct Tri6Tag   {}; // 2D, 6-node quadratic triangle
struct Quad4Tag  {}; // 2D, 4-node bilinear quadrilateral
struct Quad8Tag  {}; // 2D, 8-node serendipity quadrilateral
struct Quad9Tag  {}; // 2D, 9-node biquadratic quadrilateral
struct Tet4Tag   {}; // 3D, 4-node linear tetrahedron
struct Tet10Tag  {}; // 3D, 10-node quadratic tetrahedron
struct Hex8Tag   {}; // 3D, 8-node trilinear hexahedron
struct Hex20Tag  {}; // 3D, 20-node serendipity hexahedron
struct Hex27Tag  {}; // 3D, 27-node triquadratic hexahedron

// ── FE family tags ────────────────────────────────────────────────────────────

struct LagrangeTag    {};
struct LagrangeVecTag {};
struct HermiteTag     {};
struct MonomialTag    {};
struct MonomialVecTag {};

// ── Spatial dimension tags (for MONOMIAL) ─────────────────────────────────────
// MONOMIAL uses complete total-degree polynomials parameterised by spatial
// dimension, not element class: TRI and QUAD share Dim2Tag; TET/HEX/PRISM/
// PYRAMID share Dim3Tag.  This matches libMesh's FE<Dim, MONOMIAL> design.

struct Dim1Tag {};  // 1-D elements (EDGE)
struct Dim2Tag {};  // 2-D elements (TRI, QUAD)
struct Dim3Tag {};  // 3-D elements (TET, HEX, PRISM, PYRAMID)

// ── Primary FEEvaluator template ─────────────────────────────────────────────
//
// All uses must be explicit specializations. Every specialization must provide:
//
//   static constexpr unsigned int n_dofs()
//
//   KOKKOS_INLINE_FUNCTION
//   static Real shape(unsigned int i, Real xi, Real eta, Real zeta)
//
//   KOKKOS_INLINE_FUNCTION
//   static Real3 grad_shape(unsigned int i, Real xi, Real eta, Real zeta)
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

// Order is only meaningful for MONOMIAL specializations.
// Lagrange specializations always use Order = 0 (the default) because the
// element tag (e.g. Quad9Tag) already encodes the polynomial order.
template <typename FamilyTag, typename ElemTag, unsigned int Order = 0>
struct FEEvaluator; // forward declaration only; instantiation requires a specialization

} // namespace libMesh::Kokkos
