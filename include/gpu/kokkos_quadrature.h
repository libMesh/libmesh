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
    return Quadrature::Gauss::gauss_legendre_rule(alg_order).count;
  }

  LIBMESH_DEVICE_INLINE static Real point(unsigned int alg_order, unsigned int i)
  {
    const auto rule = Quadrature::Gauss::gauss_legendre_rule(alg_order);
    return (i < rule.count) ? rule.points[i] : 0.0;
  }

  LIBMESH_DEVICE_INLINE static Real weight(unsigned int alg_order, unsigned int i)
  {
    const auto rule = Quadrature::Gauss::gauss_legendre_rule(alg_order);
    return (i < rule.count) ? rule.weights[i] : 0.0;
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
        return Quadrature::Gauss::triangle_rule(order).count;

      case libMesh::TET4:
      case libMesh::TET10:
        return Quadrature::Gauss::tetrahedron_rule(order).count;

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
        const auto rule = Quadrature::Gauss::gauss_legendre_rule(order);
        const unsigned int n = rule.count;
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
        const auto rule = Quadrature::Gauss::gauss_legendre_rule(order);
        const unsigned int n = rule.count;
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
        const auto rule = Quadrature::Gauss::triangle_rule(order);
        return (qp < rule.count) ? make_vector(rule.points[qp].x, rule.points[qp].y, 0.0) : zero_vector();
      }

      case libMesh::TET4:
      case libMesh::TET10:
      {
        const auto rule = Quadrature::Gauss::tetrahedron_rule(order);
        return (qp < rule.count)
                 ? make_vector(rule.points[qp].x, rule.points[qp].y, rule.points[qp].z)
                 : zero_vector();
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
        const auto rule = Quadrature::Gauss::gauss_legendre_rule(order);
        const unsigned int n = rule.count;
        if (!n)
          return 0.0;
        return GaussLegendre1D::weight(order, qp % n) *
               GaussLegendre1D::weight(order, qp / n);
      }

      case libMesh::HEX8:
      case libMesh::HEX20:
      case libMesh::HEX27:
      {
        const auto rule = Quadrature::Gauss::gauss_legendre_rule(order);
        const unsigned int n = rule.count;
        if (!n)
          return 0.0;
        return GaussLegendre1D::weight(order, qp % n) *
               GaussLegendre1D::weight(order, (qp / n) % n) *
               GaussLegendre1D::weight(order, qp / (n * n));
      }

      case libMesh::TRI3:
      case libMesh::TRI6:
      {
        const auto rule = Quadrature::Gauss::triangle_rule(order);
        return (qp < rule.count) ? rule.points[qp].w : 0.0;
      }

      case libMesh::TET4:
      case libMesh::TET10:
      {
        const auto rule = Quadrature::Gauss::tetrahedron_rule(order);
        return (qp < rule.count) ? rule.points[qp].w : 0.0;
      }

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
