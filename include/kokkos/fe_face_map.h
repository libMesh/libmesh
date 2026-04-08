//* This file is part of the libMesh library.
//* https://www.libmesh.net
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#ifdef LIBMESH_HAVE_KOKKOS

#include "libmesh/kokkos/fe_evaluator.h"
#include "libmesh/elem.h"

namespace libMesh::Kokkos
{

/**
 * Map a face quadrature point from the side element's reference coordinate system
 * to the parent element's reference coordinate system.
 *
 * side_in_parent must be obtained via build_side_ptr() (not side_ptr()), so that
 * second-order sides carry their midpoint nodes.
 *
 * @param side_in_parent  The side element as embedded in the parent (from build_side_ptr())
 * @param side_topo       Topology of the side element
 * @param face_qpt        Quadrature point in the side element's reference coordinates
 * @returns               Corresponding point in the parent element's reference coordinates
 */
inline Real3
mapFaceQpToParent(const libMesh::Elem & side_in_parent,
                  FEElemTopology side_topo,
                  Real3 face_qpt)
{
  const unsigned int n = side_in_parent.n_nodes();
  Real3 parent_pt;

  // 1-D elements: the "side" is a single vertex node.  There are no reference
  // coordinates on a point element, so the parent reference coordinate is
  // simply the node's own position in the parent reference element.
  if (n == 1)
  {
    const auto & pt = side_in_parent.point(0);
    parent_pt.v[0] = pt(0);
    parent_pt.v[1] = pt(1);
    parent_pt.v[2] = pt(2);
    return parent_pt;
  }

  for (unsigned int k = 0; k < n; ++k)
  {
    Real psi;
    const Real s = face_qpt.v[0];
    const Real t = face_qpt.v[1];

    switch (side_topo)
    {
      case FEElemTopology::EDGE2:
        psi = FEEvaluator<LagrangeTag, Edge2Tag>::shape(k, s, 0.0, 0.0);
        break;
      case FEElemTopology::EDGE3:
        psi = FEEvaluator<LagrangeTag, Edge3Tag>::shape(k, s, 0.0, 0.0);
        break;
      case FEElemTopology::TRI3:
        psi = FEEvaluator<LagrangeTag, Tri3Tag>::shape(k, s, t, 0.0);
        break;
      case FEElemTopology::TRI6:
        psi = FEEvaluator<LagrangeTag, Tri6Tag>::shape(k, s, t, 0.0);
        break;
      case FEElemTopology::QUAD4:
        psi = FEEvaluator<LagrangeTag, Quad4Tag>::shape(k, s, t, 0.0);
        break;
      case FEElemTopology::QUAD8:
        psi = FEEvaluator<LagrangeTag, Quad8Tag>::shape(k, s, t, 0.0);
        break;
      case FEElemTopology::QUAD9:
        psi = FEEvaluator<LagrangeTag, Quad9Tag>::shape(k, s, t, 0.0);
        break;
      default:
        psi = 0.0;
        break;
    }

    const auto & pt = side_in_parent.point(k);
    parent_pt.v[0] += psi * pt(0);
    parent_pt.v[1] += psi * pt(1);
    parent_pt.v[2] += psi * pt(2);
  }

  return parent_pt;
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_HAVE_KOKKOS
