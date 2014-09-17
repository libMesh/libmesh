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
#include "libmesh/face_tri3.h"

namespace libMesh
{



// ------------------------------------------------------------
// Tri3 class static member initializations
const unsigned int Tri3::side_nodes_map[3][2] =
  {
    {0, 1}, // Side 0
    {1, 2}, // Side 1
    {2, 0}  // Side 2
  };


#ifdef LIBMESH_ENABLE_AMR

const float Tri3::_embedding_matrix[4][3][3] =
  {
    // embedding matrix for child 0
    {
      // 0    1    2
      {1.0, 0.0, 0.0}, // 0
      {0.5, 0.5, 0.0}, // 1
      {0.5, 0.0, 0.5}  // 2
    },

    // embedding matrix for child 1
    {
      // 0    1    2
      {0.5, 0.5, 0.0}, // 0
      {0.0, 1.0, 0.0}, // 1
      {0.0, 0.5, 0.5}  // 2
    },

    // embedding matrix for child 2
    {
      // 0    1    2
      {0.5, 0.0, 0.5}, // 0
      {0.0, 0.5, 0.5}, // 1
      {0.0, 0.0, 1.0}  // 2
    },

    // embedding matrix for child 3
    {
      // 0    1    2
      {0.5, 0.5, 0.0}, // 0
      {0.0, 0.5, 0.5}, // 1
      {0.5, 0.0, 0.5}  // 2
    }
  };

#endif



// ------------------------------------------------------------
// Tri3 class member functions

bool Tri3::is_vertex(const unsigned int) const
{
  return true;
}

bool Tri3::is_edge(const unsigned int) const
{
  return false;
}

bool Tri3::is_face(const unsigned int) const
{
  return false;
}

bool Tri3::is_node_on_side(const unsigned int n,
                           const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  for (unsigned int i = 0; i != 2; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

AutoPtr<Elem> Tri3::build_side (const unsigned int i,
                                bool proxy) const
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    return AutoPtr<Elem>(new Side<Edge2,Tri3>(this,i));

  else
    {
      Elem* edge = new Edge2;
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


void Tri3::connectivity(const unsigned int libmesh_dbg_var(sf),
                        const IOPackage iop,
                        std::vector<dof_id_type>& conn) const
{
  libmesh_assert_less (sf, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        conn.resize(4);
        conn[0] = this->node(0)+1;
        conn[1] = this->node(1)+1;
        conn[2] = this->node(2)+1;
        conn[3] = this->node(2)+1;
        return;
      }

    case VTK:
      {
        conn.resize(3);
        conn[0] = this->node(0);
        conn[1] = this->node(1);
        conn[2] = this->node(2);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}






Real Tri3::volume () const
{
  // 3-node triangles have the following formula for computing the area
  Point v10 ( *(this->get_node(1)) - *(this->get_node(0)) );

  Point v20 ( *(this->get_node(2)) - *(this->get_node(0)) );

  return 0.5 * (v10.cross(v20)).size() ;
}



std::pair<Real, Real> Tri3::min_and_max_angle() const
{
  Point v10 ( this->point(1) - this->point(0) );
  Point v20 ( this->point(2) - this->point(0) );
  Point v21 ( this->point(2) - this->point(1) );

  const Real
    len_10=v10.size(),
    len_20=v20.size(),
    len_21=v21.size()
    ;

  const Real
    theta0=std::acos(( v10*v20)/len_10/len_20),
    theta1=std::acos((-v10*v21)/len_10/len_21),
    theta2=libMesh::pi - theta0 - theta1
    ;

  libmesh_assert_greater (theta0, 0.);
  libmesh_assert_greater (theta1, 0.);
  libmesh_assert_greater (theta2, 0.);

  return std::make_pair(std::min(theta0, std::min(theta1,theta2)),
                        std::max(theta0, std::max(theta1,theta2)));
}

} // namespace libMesh
