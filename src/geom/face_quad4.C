// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// C++ includes

// Local includes
#include "libmesh/side.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_quad4.h"

namespace libMesh
{




// ------------------------------------------------------------
// Quad class static member initialization
const unsigned int Quad4::side_nodes_map[4][2] =
  {
    {0, 1}, // Side 0
    {1, 2}, // Side 1
    {2, 3}, // Side 2
    {3, 0}  // Side 3
  };


#ifdef LIBMESH_ENABLE_AMR

const float Quad4::_embedding_matrix[4][4][4] =
  {
    // embedding matrix for child 0
    {
      // 0    1    2    3
      {1.0, 0.0, 0.0, 0.0}, // 0
      {0.5, 0.5, 0.0, 0.0}, // 1
      {.25, .25, .25, .25}, // 2
      {0.5, 0.0, 0.0, 0.5}  // 3
    },

    // embedding matrix for child 1
    {
      // 0    1    2    3
      {0.5, 0.5, 0.0, 0.0}, // 0
      {0.0, 1.0, 0.0, 0.0}, // 1
      {0.0, 0.5, 0.5, 0.0}, // 2
      {.25, .25, .25, .25}  // 3
    },

    // embedding matrix for child 2
    {
      // 0    1    2    3
      {0.5, 0.0, 0.0, 0.5}, // 0
      {.25, .25, .25, .25}, // 1
      {0.0, 0.0, 0.5, 0.5}, // 2
      {0.0, 0.0, 0.0, 1.0}  // 3
    },

    // embedding matrix for child 3
    {
      // 0    1    2    3
      {.25, .25, .25, .25}, // 0
      {0.0, 0.5, 0.5, 0.0}, // 1
      {0.0, 0.0, 1.0, 0.0}, // 2
      {0.0, 0.0, 0.5, 0.5}  // 3
    }
  };

#endif





// ------------------------------------------------------------
// Quad4 class member functions

bool Quad4::is_vertex(const unsigned int) const
{
  return true;
}

bool Quad4::is_edge(const unsigned int) const
{
  return false;
}

bool Quad4::is_face(const unsigned int) const
{
  return false;
}

bool Quad4::is_node_on_side(const unsigned int n,
                            const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  for (unsigned int i = 0; i != 2; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}



bool Quad4::has_affine_map() const
{
  Point v = this->point(3) - this->point(0);
  return (v.relative_fuzzy_equals(this->point(2) - this->point(1)));
}



std::unique_ptr<Elem> Quad4::build_side_ptr (const unsigned int i,
                                             bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    return libmesh_make_unique<Side<Edge2,Quad4>>(this,i);

  else
    {
      std::unique_ptr<Elem> edge = libmesh_make_unique<Edge2>();
      edge->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (unsigned n=0; n<edge->n_nodes(); ++n)
        edge->set_node(n) = this->node_ptr(Quad4::side_nodes_map[i][n]);

      return edge;
    }
}





void Quad4::connectivity(const unsigned int libmesh_dbg_var(sf),
                         const IOPackage iop,
                         std::vector<dof_id_type> & conn) const
{
  libmesh_assert_less (sf, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  // Create storage.
  conn.resize(4);

  switch (iop)
    {
    case TECPLOT:
      {
        conn[0] = this->node_id(0)+1;
        conn[1] = this->node_id(1)+1;
        conn[2] = this->node_id(2)+1;
        conn[3] = this->node_id(3)+1;
        return;
      }

    case VTK:
      {
        conn[0] = this->node_id(0);
        conn[1] = this->node_id(1);
        conn[2] = this->node_id(2);
        conn[3] = this->node_id(3);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



Real Quad4::volume () const
{
  // Make copies of our points.  It makes the subsequent calculations a bit
  // shorter and avoids dereferencing the same pointer multiple times.
  Point
    x0 = point(0), x1 = point(1),
    x2 = point(2), x3 = point(3);

  // Construct constant data vectors.
  // \vec{x}_{\xi}  = \vec{a1}*eta + \vec{b1}
  // \vec{x}_{\eta} = \vec{a2}*xi  + \vec{b2}
  // This is copy-pasted directly from the output of a Python script.
  Point
    a1 = x0/4 - x1/4 + x2/4 - x3/4,
    b1 = -x0/4 + x1/4 + x2/4 - x3/4,
    a2 = a1,
    b2 = -x0/4 - x1/4 + x2/4 + x3/4;

  // Check for quick return for parallelogram QUAD4.
  if (a1.relative_fuzzy_equals(Point(0,0,0)))
    return 4. * b1.cross(b2).norm();

  // Otherwise, use 2x2 quadrature to approximate the surface area.

  // 4-point rule, exact for bi-cubics.  The weights for this rule are
  // all equal to 1.
  const Real q[2] = {-std::sqrt(3.)/3, std::sqrt(3.)/3.};

  Real vol=0.;
  for (unsigned int i=0; i<2; ++i)
    for (unsigned int j=0; j<2; ++j)
      vol += cross_norm(q[j]*a1 + b1,
                        q[i]*a2 + b2);

  return vol;
}

} // namespace libMesh
