// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/cell_pyramid5.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_quad4.h"

namespace libMesh
{




// ------------------------------------------------------------
// Pyramid5 class static member initializations
const unsigned int Pyramid5::side_nodes_map[5][4] =
  {
    {0, 1, 4, 99}, // Side 0
    {1, 2, 4, 99}, // Side 1
    {2, 3, 4, 99}, // Side 2
    {3, 0, 4, 99}, // Side 3
    {0, 3, 2,  1}  // Side 4
  };

const unsigned int Pyramid5::edge_nodes_map[8][2] =
  {
    {0, 1}, // Side 0
    {1, 2}, // Side 1
    {2, 3}, // Side 2
    {0, 3}, // Side 3
    {0, 4}, // Side 4
    {1, 4}, // Side 5
    {2, 4}, // Side 6
    {3, 4}  // Side 7
  };



// ------------------------------------------------------------
// Pyramid5 class member functions

bool Pyramid5::is_vertex(const unsigned int) const
{
  return true;
}

bool Pyramid5::is_edge(const unsigned int) const
{
  return false;
}

bool Pyramid5::is_face(const unsigned int) const
{
  return false;
}

bool Pyramid5::is_node_on_side(const unsigned int n,
                               const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  for (unsigned int i = 0; i != 4; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool Pyramid5::is_node_on_edge(const unsigned int n,
                               const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  for (unsigned int i = 0; i != 2; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
}



bool Pyramid5::has_affine_map() const
{
  //  Point v = this->point(3) - this->point(0);
  //  return (v.relative_fuzzy_equals(this->point(2) - this->point(1)));
  return false;
}



UniquePtr<Elem> Pyramid5::build_side (const unsigned int i,
                                      bool proxy) const
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    {
      switch (i)
        {
        case 0:
        case 1:
        case 2:
        case 3:
          return UniquePtr<Elem>(new Side<Tri3,Pyramid5>(this,i));

        case 4:
          return UniquePtr<Elem>(new Side<Quad4,Pyramid5>(this,i));

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }
    }

  else
    {
      // Create NULL pointer to be initialized, returned later.
      Elem * face = libmesh_nullptr;

      switch (i)
        {
        case 0: // triangular face 1
        case 1: // triangular face 2
        case 2: // triangular face 3
        case 3: // triangular face 4
          {
            face = new Tri3;
            break;
          }
        case 4: // the quad face at z=0
          {
            face = new Quad4;
            break;
          }
        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      face->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (unsigned n=0; n<face->n_nodes(); ++n)
        face->set_node(n) = this->get_node(Pyramid5::side_nodes_map[i][n]);

      return UniquePtr<Elem>(face);
    }

  libmesh_error_msg("We'll never get here!");
  return UniquePtr<Elem>();
}



UniquePtr<Elem> Pyramid5::build_edge (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_edges());

  return UniquePtr<Elem>(new SideEdge<Edge2,Pyramid5>(this,i));
}



void Pyramid5::connectivity(const unsigned int libmesh_dbg_var(sc),
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const
{
  libmesh_assert(_nodes);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        conn.resize(8);
        conn[0] = this->node(0)+1;
        conn[1] = this->node(1)+1;
        conn[2] = this->node(2)+1;
        conn[3] = this->node(3)+1;
        conn[4] = this->node(4)+1;
        conn[5] = this->node(4)+1;
        conn[6] = this->node(4)+1;
        conn[7] = this->node(4)+1;
        return;
      }

    case VTK:
      {
        conn.resize(5);
        conn[0] = this->node(3);
        conn[1] = this->node(2);
        conn[2] = this->node(1);
        conn[3] = this->node(0);
        conn[4] = this->node(4);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}


Real Pyramid5::volume () const
{
  // Make copies of our points.  It makes the subsequent calculations a bit
  // shorter and avoids dereferencing the same pointer multiple times.
  Point
    x0 = point(0), x1 = point(1), x2 = point(2), x3 = point(3), x4 = point(4);

  // The number of components in the dx_dxi, dx_deta, and dx_dzeta arrays.
  const int n_components = 4;

  Point dx_dxi[n_components] =
    {
      -x0/4 + x1/4 + x2/4 - x3/4, // const
      x0/4 - x1/4 - x2/4 + x3/4,  // zeta
      x0/4 - x1/4 + x2/4 - x3/4   // eta
    };

  Point dx_deta[n_components] =
    {
      -x0/4 - x1/4 + x2/4 + x3/4, // const
      x0/4 + x1/4 - x2/4 - x3/4,  // zeta
      x0/4 - x1/4 + x2/4 - x3/4,  // xi
    };

  Point dx_dzeta[n_components] =
    {
      -x0/4 - x1/4 - x2/4 - x3/4 + x4,  // const
      x0/2 + x1/2 + x2/2 + x3/2 - 2*x4, // zeta
      -x0/4 - x1/4 - x2/4 - x3/4 + x4,  // zeta**2
      x0/4 - x1/4 + x2/4 - x3/4         // xi*eta
    };

  // Number of points in the 2D quadrature rule
  const int N = 8;

  // Parameters of the quadrature rule
  static const Real
    a1 = -5.0661630334978742377431469429037e-01L,
    a2 = -2.6318405556971359557121701304557e-01L,
    b1 = 1.2251482265544137786674043037115e-01L,
    b2 = 5.4415184401122528879992623629551e-01L,
    w1 = 2.3254745125350790274997694884235e-01L,
    w2 = 1.0078588207982543058335638449098e-01L;

  // The points and weights of the 2x2x2 quadrature rule
  static const Real xi[8]   = {a1,  a2,  a1,  a2, -a1, -a2, -a1, -a2};
  static const Real eta[8]  = {a1,  a2, -a1, -a2,  a1,  a2, -a1, -a2};
  static const Real zeta[8] = {b1,  b2,  b1,  b2,  b1,  b2,  b1,  b2};
  static const Real w[8] = {w1, w2, w1, w2, w1, w2, w1, w2};

  Real vol = 0.;
  for (int q=0; q<N; ++q)
    {
      // Note: we need to scale dx/dxi and dx/deta by (1-z), and dx/dzeta by (1-z)**2
      Point
        dx_dxi_q   = (dx_dxi[0] + zeta[q]*dx_dxi[1] + eta[q]*dx_dxi[2]) / (1. - zeta[q]),
        dx_deta_q  = (dx_deta[0] + zeta[q]*dx_deta[1] + xi[q]*dx_deta[2]) / (1. - zeta[q]),
        dx_dzeta_q = (dx_dzeta[0] + zeta[q]*dx_dzeta[1] + zeta[q]*zeta[q]*dx_dzeta[2] + xi[q]*eta[q]*dx_dzeta[3]) / (1. - zeta[q]) / (1. - zeta[q]);

      // Compute scalar triple product, multiply by weight, and accumulate volume.
      vol += w[q] * dx_dxi_q * dx_deta_q.cross(dx_dzeta_q);
    }

  return vol;
}

} // namespace libMesh
