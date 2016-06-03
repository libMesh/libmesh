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

UniquePtr<Elem> Tri3::build_side_ptr (const unsigned int i,
                                      bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    return UniquePtr<Elem>(new Side<Edge2,Tri3>(this,i));

  else
    {
      Elem * edge = new Edge2;
      edge->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (unsigned n=0; n<edge->n_nodes(); ++n)
        edge->set_node(n) = this->node_ptr(Tri3::side_nodes_map[i][n]);

      return UniquePtr<Elem>(edge);
    }

  libmesh_error_msg("We'll never get here!");
  return UniquePtr<Elem>();
}


void Tri3::connectivity(const unsigned int libmesh_dbg_var(sf),
                        const IOPackage iop,
                        std::vector<dof_id_type> & conn) const
{
  libmesh_assert_less (sf, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
        conn.resize(4);
        conn[0] = this->node_id(0)+1;
        conn[1] = this->node_id(1)+1;
        conn[2] = this->node_id(2)+1;
        conn[3] = this->node_id(2)+1;
        return;
      }

    case VTK:
      {
        conn.resize(3);
        conn[0] = this->node_id(0);
        conn[1] = this->node_id(1);
        conn[2] = this->node_id(2);
        return;
      }

    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}






Real Tri3::volume () const
{
  // 3-node triangles have the following formula for computing the area
  return 0.5 * cross_norm(point(1) - point(0),
                          point(2) - point(0));
}



std::pair<Real, Real> Tri3::min_and_max_angle() const
{
  Point v10 ( this->point(1) - this->point(0) );
  Point v20 ( this->point(2) - this->point(0) );
  Point v21 ( this->point(2) - this->point(1) );

  const Real
    len_10=v10.norm(),
    len_20=v20.norm(),
    len_21=v21.norm();

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

bool Tri3::contains_point (const Point & p, Real tol) const
{
  // move test point relative to node 0
  const Point p0 ( p - this->point(0) );

  // define two in plane vectors
  Point v1 ( this->point(1) - this->point(0) );
  Point v2 ( this->point(2) - this->point(0) );

#if LIBMESH_DIM == 3
  // define out of plane vector
  Point oop ( v1.cross(v2) );

  // project moved test point onto out of plane component and bail
  // if it is farther out of plane than tol.
  if ( std::abs(p0 * oop) > tol || oop.norm_sq() < tol * tol)
    return false;
#endif

  const Real d01 = p0 * v1;
  const Real d02 = p0 * v2;
  const Real d11 = v1.norm_sq();
  const Real d12 = v1 * v2;
  const Real d22 = v2.norm_sq();

#if LIBMESH_DIM == 3
  // the denominator can only be 0 if v1 = v2, but then the cross product would be 0
  // and we'd have bailed out already
  const Real d = 1.0 / (d11 * d22 - d12 * d12);
#else
  Real d = d11 * d22 - d12 * d12;
  if (d < tol)
    return false;
  d = 1.0 / d;
#endif

  // the in plane check implements the barycentric technique from
  // http://www.blackpawn.com/texts/pointinpoly/
  const Real s1 = (d02 * d11 - d01 * d12) * d;
  if (s1 <= -tol)
    return false;

  const Real s2 = (d01 * d22 - d02 * d12) * d;
  return s2 > -tol && s1 + s2 < 1.0 + tol;
}

} // namespace libMesh
