// Kokkos FE access to the shared libMesh Gauss quadrature rule tables.

#ifndef LIBMESH_KOKKOS_QUADRATURE_H
#define LIBMESH_KOKKOS_QUADRATURE_H

#include "kokkos_fe_base.h"
#include "libmesh/enum_elem_type.h"
#include "libmesh/quadrature_gauss_rules.h"

#include <vector>

namespace libMesh::Kokkos
{

struct GaussLegendre1D
{
  LIBMESH_DEVICE_INLINE static unsigned int n_points(unsigned int alg_order)
  {
    return Quadrature::Gauss::gauss_legendre_count(alg_order);
  }

  LIBMESH_DEVICE_INLINE static Real point(unsigned int alg_order, unsigned int i)
  {
    return Quadrature::Gauss::gauss_legendre_point(alg_order, i);
  }

  LIBMESH_DEVICE_INLINE static Real weight(unsigned int alg_order, unsigned int i)
  {
    return Quadrature::Gauss::gauss_legendre_weight(alg_order, i);
  }
};

struct GaussQuadrature
{
  LIBMESH_DEVICE_INLINE static unsigned int
  n_points(libMesh::ElemType topo, unsigned int order)
  {
    switch (topo)
    {
      case libMesh::EDGE2:
      case libMesh::EDGE3:
        return GaussLegendre1D::n_points(order);

      case libMesh::QUAD4:
      case libMesh::QUAD8:
      case libMesh::QUAD9:
      {
        const unsigned int n = GaussLegendre1D::n_points(order);
        return n * n;
      }

      case libMesh::HEX8:
      case libMesh::HEX20:
      case libMesh::HEX27:
      {
        const unsigned int n = GaussLegendre1D::n_points(order);
        return n * n * n;
      }

      case libMesh::TRI3:
      case libMesh::TRI6:
        return Quadrature::Gauss::triangle_count(order);

      case libMesh::TET4:
      case libMesh::TET10:
        return Quadrature::Gauss::tetrahedron_count(order);

      default:
        return 0;
    }
  }

  LIBMESH_DEVICE_INLINE static RealVector
  point(libMesh::ElemType topo, unsigned int order, unsigned int qp)
  {
    switch (topo)
    {
      case libMesh::EDGE2:
      case libMesh::EDGE3:
        return make_vector(GaussLegendre1D::point(order, qp), 0.0, 0.0);

      case libMesh::QUAD4:
      case libMesh::QUAD8:
      case libMesh::QUAD9:
      {
        const unsigned int n = Quadrature::Gauss::gauss_legendre_count(order);
        if (!n)
          return zero_vector();
        const unsigned int i = qp % n;
        const unsigned int j = qp / n;
        return make_vector(GaussLegendre1D::point(order, i),
                           GaussLegendre1D::point(order, j),
                           0.0);
      }

      case libMesh::HEX8:
      case libMesh::HEX20:
      case libMesh::HEX27:
      {
        const unsigned int n = Quadrature::Gauss::gauss_legendre_count(order);
        if (!n)
          return zero_vector();
        const unsigned int i = qp % n;
        const unsigned int j = (qp / n) % n;
        const unsigned int k = qp / (n * n);
        return make_vector(GaussLegendre1D::point(order, i),
                           GaussLegendre1D::point(order, j),
                           GaussLegendre1D::point(order, k));
      }

      case libMesh::TRI3:
      case libMesh::TRI6:
      {
        if (qp >= Quadrature::Gauss::triangle_count(order))
          return zero_vector();

        const auto point = Quadrature::Gauss::triangle_point(order, qp);
        return make_vector(point.x, point.y, 0.0);
      }

      case libMesh::TET4:
      case libMesh::TET10:
      {
        if (qp >= Quadrature::Gauss::tetrahedron_count(order))
          return zero_vector();

        const auto point = Quadrature::Gauss::tetrahedron_point(order, qp);
        return make_vector(point.x, point.y, point.z);
      }

      default:
        return zero_vector();
    }
  }

  LIBMESH_DEVICE_INLINE static Real
  weight(libMesh::ElemType topo, unsigned int order, unsigned int qp)
  {
    switch (topo)
    {
      case libMesh::EDGE2:
      case libMesh::EDGE3:
        return GaussLegendre1D::weight(order, qp);

      case libMesh::QUAD4:
      case libMesh::QUAD8:
      case libMesh::QUAD9:
      {
        const unsigned int n = Quadrature::Gauss::gauss_legendre_count(order);
        if (!n)
          return 0.0;
        return GaussLegendre1D::weight(order, qp % n) *
               GaussLegendre1D::weight(order, qp / n);
      }

      case libMesh::HEX8:
      case libMesh::HEX20:
      case libMesh::HEX27:
      {
        const unsigned int n = Quadrature::Gauss::gauss_legendre_count(order);
        if (!n)
          return 0.0;
        return GaussLegendre1D::weight(order, qp % n) *
               GaussLegendre1D::weight(order, (qp / n) % n) *
               GaussLegendre1D::weight(order, qp / (n * n));
      }

      case libMesh::TRI3:
      case libMesh::TRI6:
        return (qp < Quadrature::Gauss::triangle_count(order)) ? Quadrature::Gauss::triangle_weight(order, qp) : 0.0;

      case libMesh::TET4:
      case libMesh::TET10:
        return (qp < Quadrature::Gauss::tetrahedron_count(order)) ? Quadrature::Gauss::tetrahedron_weight(order, qp) : 0.0;

      default:
        return 0.0;
    }
  }
};

inline void
fill_quadrature(libMesh::ElemType topo,
                unsigned int order,
                std::vector<RealVector> & qpts,
                std::vector<Real> & weights)
{
  const unsigned int nqp = GaussQuadrature::n_points(topo, order);
  qpts.resize(nqp);
  weights.resize(nqp);
  for (unsigned int q = 0; q < nqp; ++q)
  {
    qpts[q] = GaussQuadrature::point(topo, order, q);
    weights[q] = GaussQuadrature::weight(topo, order, q);
  }
}

} // namespace libMesh::Kokkos

#endif // LIBMESH_KOKKOS_QUADRATURE_H
