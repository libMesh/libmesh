// The libMesh Finite Element Library.
// Copyright (C) 2002-2014 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



AutoPtr<Elem> Quad4::build_side (const unsigned int i,
                                 bool proxy) const
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    return AutoPtr<Elem>(new Side<Edge2,Quad4>(this,i));

  else
    {
      Elem * edge = new Edge2;
      edge->subdomain_id() = this->subdomain_id();

      switch (i)
        {
        case 0:
          {
            edge->set_node(0) = this->get_node(0);
            edge->set_node(1) = this->get_node(1);
            break;
          }
        case 1:
          {
            edge->set_node(0) = this->get_node(1);
            edge->set_node(1) = this->get_node(2);
            break;
          }
        case 2:
          {
            edge->set_node(0) = this->get_node(2);
            edge->set_node(1) = this->get_node(3);
            break;
          }
        case 3:
          {
            edge->set_node(0) = this->get_node(3);
            edge->set_node(1) = this->get_node(0);
            break;
          }
        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      return AutoPtr<Elem>(edge);
    }

  libmesh_error_msg("We'll never get here!");
  return AutoPtr<Elem>();
}





void Quad4::connectivity(const unsigned int libmesh_dbg_var(sf),
                         const IOPackage iop,
                         std::vector<dof_id_type>& conn) const
{
  libmesh_assert_less (sf, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  // Create storage.
  conn.resize(4);

  switch (iop)
    {
    case TECPLOT:
      {
        conn[0] = this->node(0)+1;
        conn[1] = this->node(1)+1;
        conn[2] = this->node(2)+1;
        conn[3] = this->node(3)+1;
        return;
      }

    case VTK:
      {
        conn[0] = this->node(0);
        conn[1] = this->node(1);
        conn[2] = this->node(2);
        conn[3] = this->node(3);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



Real Quad4::volume () const
{
  // The A,B,C,D naming scheme here corresponds exactly to the
  // libmesh counter-clockwise numbering scheme.

  //        3           2        D           C
  // QUAD4: o-----------o        o-----------o
  //        |           |        |           |
  //        |           |        |           |
  //        |           |        |           |
  //        |           |        |           |
  //        |           |        |           |
  //        o-----------o        o-----------o
  //        0           1        A           B

  // Vector pointing from A to C
  Point AC ( this->point(2) - this->point(0) );

  // Vector pointing from A to B
  Point AB ( this->point(1) - this->point(0) );

  // Vector pointing from A to D
  Point AD ( this->point(3) - this->point(0) );

  // The diagonal vector minus the side vectors
  Point AC_AB_AD (AC - AB - AD);

  // Check for quick return for planar QUAD4.  This will
  // be the most common case, occuring for all purely 2D meshes.
  if (AC_AB_AD == Point(0.,0.,0.))
    return AB.cross(AD).size();

  else
    {
      // Use 2x2 quadrature to approximate the surface area.  (The
      // true integral is too difficult to compute analytically.)  The
      // accuracy here is exactly the same as would be obtained via a
      // call to Elem::volume(), however it is a bit more optimized to
      // do it this way.  The technique used is to integrate the magnitude
      // of the normal vector over the whole area.  See for example,
      //
      // Y. Zhang, C. Bajaj, G. Xu. Surface Smoothing and Quality
      // Improvement of Quadrilateral/Hexahedral Meshes with Geometric
      // Flow. The special issue of the Journal Communications in
      // Numerical Methods in Engineering (CNME), submitted as an
      // invited paper, 2006.
      // http://www.ices.utexas.edu/~jessica/paper/quadhexgf/quadhex_geomflow_CNM.pdf

      // 4-point rule
      const Real q[2] = {0.5 - std::sqrt(3.) / 6.,
                         0.5 + std::sqrt(3.) / 6.};

      Real vol=0.;
      for (unsigned int i=0; i<2; ++i)
        for (unsigned int j=0; j<2; ++j)
          vol += (AB + q[i]*AC_AB_AD).cross(AD + q[j]*AC_AB_AD).size();

      return 0.25*vol;
    }
}

} // namespace libMesh
