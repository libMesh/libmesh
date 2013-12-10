// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/cell_pyramid14.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/face_quad9.h"

namespace libMesh
{




// ------------------------------------------------------------
// Pyramid14 class static member initializations
const unsigned int Pyramid14::side_nodes_map[5][9] =
{
  {0, 1, 4, 5, 10,  9, 99, 99, 99}, // Side 0 (front)
  {1, 2, 4, 6, 11, 10, 99, 99, 99}, // Side 1 (right)
  {2, 3, 4, 7, 12, 11, 99, 99, 99}, // Side 2 (back)
  {3, 0, 4, 8,  9, 12, 99, 99, 99}, // Side 3 (left)
  {0, 3, 2, 1,  8,  7,  6,  5, 13}  // Side 4 (base)
};

const unsigned int Pyramid14::edge_nodes_map[8][3] =
{
  {0, 1,  5}, // Edge 0
  {1, 2,  6}, // Edge 1
  {2, 3,  7}, // Edge 2
  {0, 3,  8}, // Edge 3
  {0, 4,  9}, // Edge 4
  {1, 4, 10}, // Edge 5
  {2, 4, 11}, // Edge 6
  {3, 4, 12}  // Edge 7
};



// ------------------------------------------------------------
// Pyramid14 class member functions

bool Pyramid14::is_vertex(const unsigned int i) const
{
  if (i < 5)
    return true;
  return false;
}



bool Pyramid14::is_edge(const unsigned int i) const
{
  if (i < 5)
    return false;
  if (i == 13)
    return false;
  return true;
}



bool Pyramid14::is_face(const unsigned int i) const
{
  if (i == 13)
    return true;
  return false;
}



bool Pyramid14::is_node_on_side(const unsigned int n,
                                const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  for (unsigned int i = 0; i != 9; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool Pyramid14::is_node_on_edge(const unsigned int n,
                                const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  for (unsigned int i = 0; i != 3; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
}



bool Pyramid14::has_affine_map() const
{
  // TODO: If the base is a parallelogram and all the triangular faces are planar,
  // the map should be linear, but I need to test this theory...
  return false;
}



AutoPtr<Elem> Pyramid14::build_side (const unsigned int i, bool proxy) const
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
	  {
	    AutoPtr<Elem> face(new Side<Tri6,Pyramid14>(this,i));
	    return face;
	  }

	case 4:
	  {
	    AutoPtr<Elem> face(new Side<Quad9,Pyramid14>(this,i));
	    return face;
	  }

	default:
	  {
	    libmesh_error();
	  }
	}
    }

  else
    {
      // Create NULL pointer to be initialized, returned later.
      AutoPtr<Elem> face(NULL);

      switch (i)
	{
	case 0:  // triangular face 1
	  {
            face.reset(new Tri6);

	    face->set_node(0) = this->get_node(0);
	    face->set_node(1) = this->get_node(1);
	    face->set_node(2) = this->get_node(4);
	    face->set_node(3) = this->get_node(5);
	    face->set_node(4) = this->get_node(10);
	    face->set_node(5) = this->get_node(9);

	    break;
	  }
	case 1:  // triangular face 2
	  {
            face.reset(new Tri6);

	    face->set_node(0) = this->get_node(1);
	    face->set_node(1) = this->get_node(2);
	    face->set_node(2) = this->get_node(4);
	    face->set_node(3) = this->get_node(6);
	    face->set_node(4) = this->get_node(11);
	    face->set_node(5) = this->get_node(10);

	    break;
	  }
	case 2:  // triangular face 3
	  {
            face.reset(new Tri6);

	    face->set_node(0) = this->get_node(2);
	    face->set_node(1) = this->get_node(3);
	    face->set_node(2) = this->get_node(4);
	    face->set_node(3) = this->get_node(7);
	    face->set_node(4) = this->get_node(12);
	    face->set_node(5) = this->get_node(11);

	    break;
	  }
	case 3:  // triangular face 4
	  {
            face.reset(new Tri6);

	    face->set_node(0) = this->get_node(3);
	    face->set_node(1) = this->get_node(0);
	    face->set_node(2) = this->get_node(4);
	    face->set_node(3) = this->get_node(8);
	    face->set_node(4) = this->get_node(9);
	    face->set_node(5) = this->get_node(12);

	    break;
	  }
	case 4:  // the quad face at z=0
	  {
            face.reset(new Quad9);

	    face->set_node(0) = this->get_node(0);
	    face->set_node(1) = this->get_node(3);
	    face->set_node(2) = this->get_node(2);
	    face->set_node(3) = this->get_node(1);
	    face->set_node(4) = this->get_node(8);
	    face->set_node(5) = this->get_node(7);
	    face->set_node(6) = this->get_node(6);
	    face->set_node(7) = this->get_node(5);
	    face->set_node(8) = this->get_node(13);

	    break;
	  }
	default:
	  {
	    libmesh_error();
	  }
	}

      face->subdomain_id() = this->subdomain_id();
      return face;
    }


  // We'll never get here.
  libmesh_error();
  AutoPtr<Elem> ap(NULL);  return ap;
}



AutoPtr<Elem> Pyramid14::build_edge (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_edges());

  return AutoPtr<Elem>(new SideEdge<Edge3,Pyramid14>(this,i));
}



void Pyramid14::connectivity(const unsigned int libmesh_dbg_var(sc),
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
	// FIXME
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
	// FIXME
	conn.resize(5);
	conn[0] = this->node(3);
	conn[1] = this->node(2);
	conn[2] = this->node(1);
	conn[3] = this->node(0);
	conn[4] = this->node(4);
	return;
      }

    default:
      libmesh_error();
    }

    libmesh_error();
}



  unsigned int Pyramid14::n_second_order_adjacent_vertices (const unsigned int n) const
  {
    switch (n)
      {
      case 5:
      case 6:
      case 7:
      case 8:
      case 9:
      case 10:
      case 11:
      case 12:
	return 2;

      case 13:
	return 4;

      default:
	libmesh_error();
      }

    // We'll never get here
    libmesh_error();
    return static_cast<unsigned int>(-1);
  }


  unsigned short int Pyramid14::second_order_adjacent_vertex (const unsigned int n,
                                                              const unsigned int v) const
  {
    libmesh_assert_greater_equal (n, this->n_vertices());
    libmesh_assert_less (n, this->n_nodes());

    switch (n)
      {
      case 5:
      case 6:
      case 7:
      case 8:
      case 9:
      case 10:
      case 11:
      case 12:
        {
          libmesh_assert_less (v, 2);

          // This is the analog of the static, const arrays
          // {Hex,Prism,Tet10}::_second_order_adjacent_vertices
          // defined in the respective source files... possibly treat
          // this similarly once the Pyramid13 has been added?
          unsigned short node_list[8][2] =
            {
              {0,1},
              {1,2},
              {2,3},
              {0,3},
              {0,4},
              {1,4},
              {2,4},
              {3,4}
            };

          return node_list[n-5][v];
        }

        // mid-face node on bottom
      case 13:
        {
          libmesh_assert_less (v, 4);

          // The vertex nodes surrounding node 13 are 0, 1, 2, and 3.
          // Thus, the v'th node is simply = v.
          return v;
        }

      default:
        {
          // We can't handle this n, throw an error!
          libmesh_error();
        }

      }

    // We'll never get here
    libmesh_error();
    return static_cast<unsigned short int>(-1);
  }

} // namespace libMesh
