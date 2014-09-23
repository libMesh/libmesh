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
#include "libmesh/cell_prism6.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_quad4.h"
#include "libmesh/face_tri3.h"

namespace libMesh
{



// ------------------------------------------------------------
// Prism6 class static member initializations
const unsigned int Prism6::side_nodes_map[5][4] =
  {
    {0, 2, 1, 99}, // Side 0
    {0, 1, 4,  3}, // Side 1
    {1, 2, 5,  4}, // Side 2
    {2, 0, 3,  5}, // Side 3
    {3, 4, 5, 99}  // Side 4
  };

const unsigned int Prism6::side_elems_map[5][4] =
  {
    {0, 1, 2, 3}, // Side 0
    {0, 1, 4, 5}, // Side 1
    {1, 2, 5, 6}, // Side 2
    {0, 2, 4, 6}, // Side 3
    {4, 5, 6, 7}  // Side 4
  };

const unsigned int Prism6::edge_nodes_map[9][2] =
  {
    {0, 1}, // Side 0
    {1, 2}, // Side 1
    {0, 2}, // Side 2
    {0, 3}, // Side 3
    {1, 4}, // Side 4
    {2, 5}, // Side 5
    {3, 4}, // Side 6
    {4, 5}, // Side 7
    {3, 5}  // Side 8
  };


// ------------------------------------------------------------
// Prism6 class member functions

bool Prism6::is_vertex(const unsigned int) const
{
  return true;
}

bool Prism6::is_edge(const unsigned int) const
{
  return false;
}

bool Prism6::is_face(const unsigned int) const
{
  return false;
}

bool Prism6::is_node_on_side(const unsigned int n,
                             const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  for (unsigned int i = 0; i != 4; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool Prism6::is_node_on_edge(const unsigned int n,
                             const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  for (unsigned int i = 0; i != 2; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
}



bool Prism6::has_affine_map() const
{
  // Make sure z edges are affine
  Point v = this->point(3) - this->point(0);
  if (!v.relative_fuzzy_equals(this->point(4) - this->point(1)) ||
      !v.relative_fuzzy_equals(this->point(5) - this->point(2)))
    return false;
  return true;
}



UniquePtr<Elem> Prism6::build_side (const unsigned int i,
                                    bool proxy) const
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    {
      switch(i)
        {
        case 0:
        case 4:
          return UniquePtr<Elem>(new Side<Tri3,Prism6>(this,i));

        case 1:
        case 2:
        case 3:
          return UniquePtr<Elem>(new Side<Quad4,Prism6>(this,i));

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }
    }

  else
    {
      // Create NULL pointer to be initialized, returned later.
      Elem* face = NULL;

      switch (i)
        {
        case 0:  // the triangular face at z=-1
          {
            face = new Tri3;

            face->set_node(0) = this->get_node(0);
            face->set_node(1) = this->get_node(2);
            face->set_node(2) = this->get_node(1);

            break;
          }
        case 1:  // the quad face at y=0
          {
            face = new Quad4;

            face->set_node(0) = this->get_node(0);
            face->set_node(1) = this->get_node(1);
            face->set_node(2) = this->get_node(4);
            face->set_node(3) = this->get_node(3);

            break;
          }
        case 2:  // the other quad face
          {
            face = new Quad4;

            face->set_node(0) = this->get_node(1);
            face->set_node(1) = this->get_node(2);
            face->set_node(2) = this->get_node(5);
            face->set_node(3) = this->get_node(4);

            break;
          }
        case 3: // the quad face at x=0
          {
            face = new Quad4;

            face->set_node(0) = this->get_node(2);
            face->set_node(1) = this->get_node(0);
            face->set_node(2) = this->get_node(3);
            face->set_node(3) = this->get_node(5);

            break;
          }
        case 4: // the triangular face at z=1
          {
            face = new Tri3;

            face->set_node(0) = this->get_node(3);
            face->set_node(1) = this->get_node(4);
            face->set_node(2) = this->get_node(5);

            break;
          }
        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      face->subdomain_id() = this->subdomain_id();
      return UniquePtr<Elem>(face);
    }

  libmesh_error_msg("We'll never get here!");
  return UniquePtr<Elem>();
}



UniquePtr<Elem> Prism6::build_edge (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_edges());

  return UniquePtr<Elem>(new SideEdge<Edge2,Prism6>(this,i));
}



void Prism6::connectivity(const unsigned int libmesh_dbg_var(sc),
                          const IOPackage iop,
                          std::vector<dof_id_type>& conn) const
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
        conn[3] = this->node(2)+1;
        conn[4] = this->node(3)+1;
        conn[5] = this->node(4)+1;
        conn[6] = this->node(5)+1;
        conn[7] = this->node(5)+1;
        return;
      }

    case VTK:
      {
        conn.resize(6);
        conn[0] = this->node(0);
        conn[1] = this->node(2);
        conn[2] = this->node(1);
        conn[3] = this->node(3);
        conn[4] = this->node(5);
        conn[5] = this->node(4);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}



#ifdef LIBMESH_ENABLE_AMR

const float Prism6::_embedding_matrix[8][6][6] =
  {
    // embedding matrix for child 0
    {
      //  0     1     2     3     4     5
      { 1.0,  0.0,  0.0,  0.0,  0.0,  0.0}, // 0
      { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 1
      { 0.5,  0.0,  0.5,  0.0,  0.0,  0.0}, // 2
      { 0.5,  0.0,  0.0,  0.5,  0.0,  0.0}, // 3
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 4
      { .25,  0.0,  .25,  .25,  0.0,  .25}  // 5
    },

    // embedding matrix for child 1
    {
      //  0     1     2     3     4     5
      { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 0
      { 0.0,  1.0,  0.0,  0.0,  0.0,  0.0}, // 1
      { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0}, // 2
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 3
      { 0.0,  0.5,  0.0,  0.0,  0.5,  0.0}, // 4
      { 0.0,  .25,  .25,  0.0,  .25,  .25}  // 5
    },

    // embedding matrix for child 2
    {
      //  0     1     2     3     4     5
      { 0.5,  0.0,  0.5,  0.0,  0.0,  0.0}, // 0
      { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0}, // 1
      { 0.0,  0.0,  1.0,  0.0,  0.0,  0.0}, // 2
      { .25,  0.0,  .25,  .25,  0.0,  .25}, // 3
      { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 4
      { 0.0,  0.0,  0.5,  0.0,  0.0,  0.5}  // 5
    },

    // embedding matrix for child 3
    {
      //  0     1     2     3     4     5
      { 0.5,  0.5,  0.0,  0.0,  0.0,  0.0}, // 0
      { 0.0,  0.5,  0.5,  0.0,  0.0,  0.0}, // 1
      { 0.5,  0.0,  0.5,  0.0,  0.0,  0.0}, // 2
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 3
      { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 4
      { .25,  0.0,  .25,  .25,  0.0,  .25}  // 5
    },

    // embedding matrix for child 4
    {
      //  0     1     2     3     4     5
      { 0.5,  0.0,  0.0,  0.5,  0.0,  0.0}, // 0
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 1
      { .25,  0.0,  .25,  .25,  0.0,  .25}, // 2
      { 0.0,  0.0,  0.0,  1.0,  0.0,  0.0}, // 3
      { 0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 4
      { 0.0,  0.0,  0.0,  0.5,  0.0,  0.5}  // 5
    },

    // embedding matrix for child 5
    {
      //  0     1     2     3     4     5
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 0
      { 0.0,  0.5,  0.0,  0.0,  0.5,  0.0}, // 1
      { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 2
      { 0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 3
      { 0.0,  0.0,  0.0,  0.0,  1.0,  0.0}, // 4
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5}  // 5
    },

    // embedding matrix for child 6
    {
      //  0     1     2     3     4     5
      { .25,  0.0,  .25,  .25,  0.0,  .25}, // 0
      { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 1
      { 0.0,  0.0,  0.5,  0.0,  0.0,  0.5}, // 2
      { 0.0,  0.0,  0.0,  0.5,  0.0,  0.5}, // 3
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5}, // 4
      { 0.0,  0.0,  0.0,  0.0,  0.0,  1.0}  // 5
    },

    // embedding matrix for child 7
    {
      //  0     1     2     3     4     5
      { .25,  .25,  0.0,  .25,  .25,  0.0}, // 0
      { 0.0,  .25,  .25,  0.0,  .25,  .25}, // 1
      { .25,  0.0,  .25,  .25,  0.0,  .25}, // 2
      { 0.0,  0.0,  0.0,  0.5,  0.5,  0.0}, // 3
      { 0.0,  0.0,  0.0,  0.0,  0.5,  0.5}, // 4
      { 0.0,  0.0,  0.0,  0.5,  0.0,  0.5}  // 5
    }
  };

#endif



Real Prism6::volume () const
{
  // The volume of the prism is computed by splitting
  // it into 2 tetrahedra and 3 pyramids with bilinear bases.
  // Then the volume formulae for the tetrahedron and pyramid
  // are applied and summed to obtain the prism's volume.

  static const unsigned char sub_pyr[3][4] =
    {
      {0, 1, 4, 3},
      {1, 2, 5, 4},
      {0, 3, 5, 2}
    };

  static const unsigned char sub_tet[2][3] =
    {
      {0, 1, 2},
      {5, 4, 3}
    };

  // The centroid is a convenient point to use
  // for the apex of all the pyramids.
  const Point R = this->centroid();

  // temporary storage for Nodes which form the base of the
  // subelements
  Node* base[4];

  // volume accumulation variable
  Real vol=0.;

  // Add up the sub-pyramid volumes
  for (unsigned int n=0; n<3; ++n)
    {
      // Set the nodes of the pyramid base
      for (unsigned int i=0; i<4; ++i)
        base[i] = this->_nodes[sub_pyr[n][i]];

      // Compute diff vectors
      Point a ( *base[0] - R );
      Point b ( *base[1] - *base[3] );
      Point c ( *base[2] - *base[0] );
      Point d ( *base[3] - *base[0] );
      Point e ( *base[1] - *base[0] );

      // Compute pyramid volume
      Real sub_vol = (1./6.)*(a*(b.cross(c))) + (1./12.)*(c*(d.cross(e)));

      libmesh_assert (sub_vol>0.);

      vol += sub_vol;
    }


  // Add up the sub-tet volumes
  for (unsigned int n=0; n<2; ++n)
    {
      // Set the nodes of the pyramid base
      for (unsigned int i=0; i<3; ++i)
        base[i] = this->_nodes[sub_tet[n][i]];

      // The volume of a tetrahedron is 1/6 the box product formed
      // by its base and apex vectors
      Point a ( R - *base[0] );

      // b is the vector pointing from 0 to 1
      Point b ( *base[1] - *base[0] );

      // c is the vector pointing from 0 to 2
      Point c ( *base[2] - *base[0] );

      Real sub_vol =  (1.0 / 6.0) * (a * (b.cross(c)));

      libmesh_assert (sub_vol>0.);

      vol += sub_vol;
    }


  // Done with all sub-volumes, so return
  return vol;
}

} // namespace libMesh
