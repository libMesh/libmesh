#ifndef LIBMESH_QUADRATURE_EXACTNESS_H
#define LIBMESH_QUADRATURE_EXACTNESS_H

#include <libmesh/enum_elem_type.h>
#include <libmesh/libmesh_common.h>
#include <libmesh/utility.h>

#include <algorithm>
#include <numeric>
#include <vector>

namespace quadrature_exactness
{

inline libMesh::Real
axis_integral(const unsigned int power)
{
  return (power % 2) ? libMesh::Real(0) : (libMesh::Real(2) / (power + 1));
}

inline libMesh::Real
edge_integral(const unsigned int x_power)
{
  return axis_integral(x_power);
}

inline libMesh::Real
quad_integral(const unsigned int x_power,
              const unsigned int y_power)
{
  return axis_integral(x_power) * axis_integral(y_power);
}

inline libMesh::Real
tri_integral(const unsigned int x_power,
             const unsigned int y_power)
{
  libMesh::Real analytical = 1.0;

  const unsigned int larger_power = std::max(x_power, y_power);
  const unsigned int smaller_power = std::min(x_power, y_power);

  std::vector<unsigned int> numerator(smaller_power > 1 ? smaller_power - 1 : 0);
  std::vector<unsigned int> denominator(2 + smaller_power);

  std::iota(numerator.begin(), numerator.end(), 2);
  std::iota(denominator.begin(), denominator.end(), larger_power + 1);

  for (std::size_t i = 0; i < denominator.size(); ++i)
  {
    if (i < numerator.size())
      analytical *= numerator[i];

    analytical /= denominator[i];
  }

  return analytical;
}

inline libMesh::Real
hex_integral(const unsigned int x_power,
             const unsigned int y_power,
             const unsigned int z_power)
{
  return axis_integral(x_power) * axis_integral(y_power) * axis_integral(z_power);
}

inline libMesh::Real
tet_integral(const unsigned int x_power,
             const unsigned int y_power,
             const unsigned int z_power)
{
  libMesh::Real analytical = 1.0;

  unsigned int sorted_powers[3] = {x_power, y_power, z_power};
  std::sort(sorted_powers, sorted_powers + 3);

  std::vector<unsigned int> numerator_1(sorted_powers[0] > 1 ? sorted_powers[0] - 1 : 0);
  std::vector<unsigned int> numerator_2(sorted_powers[1] > 1 ? sorted_powers[1] - 1 : 0);
  std::vector<unsigned int> denominator(3 + sorted_powers[0] + sorted_powers[1]);

  std::iota(numerator_1.begin(), numerator_1.end(), 2);
  std::iota(numerator_2.begin(), numerator_2.end(), 2);
  std::iota(denominator.begin(), denominator.end(), sorted_powers[2] + 1);

  for (std::size_t i = 0; i < denominator.size(); ++i)
  {
    if (i < numerator_1.size())
      analytical *= numerator_1[i];

    if (i < numerator_2.size())
      analytical *= numerator_2[i];

    analytical /= denominator[i];
  }

  return analytical;
}

inline libMesh::Real
prism_integral(const unsigned int x_power,
               const unsigned int y_power,
               const unsigned int z_power)
{
  return tri_integral(x_power, y_power) * axis_integral(z_power);
}

inline libMesh::Real
pyramid_integral(const unsigned int x_power,
                 const unsigned int y_power,
                 const unsigned int z_power)
{
  if (x_power % 2 || y_power % 2)
    return libMesh::Real(0);

  const unsigned int binom =
    libMesh::Utility::binomial(x_power + y_power + z_power + 3, z_power);

  return libMesh::Real(4) /
         ((x_power + 1) * (y_power + 1) * binom * (x_power + y_power + z_power + 3));
}

inline libMesh::Real
monomial_integral(const libMesh::ElemType elem_type,
                  const unsigned int x_power,
                  const unsigned int y_power = 0,
                  const unsigned int z_power = 0)
{
  switch (elem_type)
  {
    case libMesh::EDGE2:
    case libMesh::EDGE3:
    case libMesh::EDGE4:
      return edge_integral(x_power);

    case libMesh::TRI3:
    case libMesh::TRI6:
    case libMesh::TRI7:
      return tri_integral(x_power, y_power);

    case libMesh::QUAD4:
    case libMesh::QUAD8:
    case libMesh::QUAD9:
      return quad_integral(x_power, y_power);

    case libMesh::TET4:
    case libMesh::TET10:
    case libMesh::TET14:
      return tet_integral(x_power, y_power, z_power);

    case libMesh::HEX8:
    case libMesh::HEX20:
    case libMesh::HEX27:
      return hex_integral(x_power, y_power, z_power);

    case libMesh::PRISM6:
    case libMesh::PRISM15:
    case libMesh::PRISM18:
    case libMesh::PRISM20:
    case libMesh::PRISM21:
      return prism_integral(x_power, y_power, z_power);

    case libMesh::PYRAMID5:
    case libMesh::PYRAMID13:
    case libMesh::PYRAMID14:
    case libMesh::PYRAMID18:
      return pyramid_integral(x_power, y_power, z_power);

    default:
      return libMesh::Real(0);
  }
}

} // namespace quadrature_exactness

#endif // LIBMESH_QUADRATURE_EXACTNESS_H
