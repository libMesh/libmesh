// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Local includes
#include "libmesh/side.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_tri3.h"
#include "libmesh/tensor_value.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"

namespace libMesh
{



// ------------------------------------------------------------
// Tet4 class static member initializations
const int Tet4::num_nodes;
const int Tet4::num_sides;
const int Tet4::num_edges;
const int Tet4::num_children;
const int Tet4::nodes_per_side;
const int Tet4::nodes_per_edge;

const unsigned int Tet4::side_nodes_map[Tet4::num_sides][Tet4::nodes_per_side] =
  {
    {0, 2, 1}, // Side 0
    {0, 1, 3}, // Side 1
    {1, 2, 3}, // Side 2
    {2, 0, 3}  // Side 3
  };

const unsigned int Tet4::edge_nodes_map[Tet4::num_edges][Tet4::nodes_per_edge] =
  {
    {0, 1}, // Edge 0
    {1, 2}, // Edge 1
    {0, 2}, // Edge 2
    {0, 3}, // Edge 3
    {1, 3}, // Edge 4
    {2, 3}  // Edge 5
  };


// ------------------------------------------------------------
// Tet4 class member functions

bool Tet4::is_vertex(const unsigned int) const
{
  return true;
}

bool Tet4::is_edge(const unsigned int) const
{
  return false;
}

bool Tet4::is_face(const unsigned int) const
{
  return false;
}

bool Tet4::is_node_on_edge(const unsigned int n,
                           const unsigned int e) const
{
  libmesh_assert_less (e, n_edges());
  return std::find(std::begin(edge_nodes_map[e]),
                   std::end(edge_nodes_map[e]),
                   n) != std::end(edge_nodes_map[e]);
}




#ifdef LIBMESH_ENABLE_AMR

// This function only works if LIBMESH_ENABLE_AMR...
bool Tet4::is_child_on_side(const unsigned int c,
                            const unsigned int s) const
{
  // OK, for the Tet4, this is pretty obvious... it is sets of nodes
  // not equal to the current node.  But if we want this algorithm to
  // be generic and work for Tet10 also it helps to do it this way.
  const unsigned int nodes_opposite[4][3] =
    {
      {1,2,3}, // nodes opposite node 0
      {0,2,3}, // nodes opposite node 1
      {0,1,3}, // nodes opposite node 2
      {0,1,2}  // nodes opposite node 3
    };

  // Call the base class helper function
  return Tet::is_child_on_side_helper(c, s, nodes_opposite);
}

#else

bool Tet4::is_child_on_side(const unsigned int /*c*/,
                            const unsigned int /*s*/) const
{
  libmesh_not_implemented();
  return false;
}

#endif //LIBMESH_ENABLE_AMR




bool Tet4::is_node_on_side(const unsigned int n,
                           const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
Tet4::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

Order Tet4::default_order() const
{
  return FIRST;
}

std::unique_ptr<Elem> Tet4::build_side_ptr (const unsigned int i,
                                            bool proxy)
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    return libmesh_make_unique<Side<Tri3,Tet4>>(this,i);

  else
    {
      std::unique_ptr<Elem> face = libmesh_make_unique<Tri3>();
      face->subdomain_id() = this->subdomain_id();

      for (unsigned n=0; n<face->n_nodes(); ++n)
        face->set_node(n) = this->node_ptr(Tet4::side_nodes_map[i][n]);

      return face;
    }
}



void Tet4::build_side_ptr (std::unique_ptr<Elem> & side,
                           const unsigned int i)
{
  this->simple_build_side_ptr<Tet4>(side, i, TRI3);
}



std::unique_ptr<Elem> Tet4::build_edge_ptr (const unsigned int i)
{
  libmesh_assert_less (i, this->n_edges());

  return libmesh_make_unique<SideEdge<Edge2,Tet4>>(this,i);
}


void Tet4::connectivity(const unsigned int libmesh_dbg_var(sc),
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
        conn[0] = this->node_id(0)+1;
        conn[1] = this->node_id(1)+1;
        conn[2] = this->node_id(2)+1;
        conn[3] = this->node_id(2)+1;
        conn[4] = this->node_id(3)+1;
        conn[5] = this->node_id(3)+1;
        conn[6] = this->node_id(3)+1;
        conn[7] = this->node_id(3)+1;
        return;
      }

    case VTK:
      {
        conn.resize(4);
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



#ifdef LIBMESH_ENABLE_AMR

const float Tet4::_embedding_matrix[Tet4::num_children][Tet4::num_nodes][Tet4::num_nodes] =
  {
    // embedding matrix for child 0
    {
      // 0    1    2    3
      {1.0, 0.0, 0.0, 0.0}, // 0
      {0.5, 0.5, 0.0, 0.0}, // 1
      {0.5, 0.0, 0.5, 0.0}, // 2
      {0.5, 0.0, 0.0, 0.5}  // 3
    },

    // embedding matrix for child 1
    {
      // 0    1    2    3
      {0.5, 0.5, 0.0, 0.0}, // 0
      {0.0, 1.0, 0.0, 0.0}, // 1
      {0.0, 0.5, 0.5, 0.0}, // 2
      {0.0, 0.5, 0.0, 0.5}  // 3
    },

    // embedding matrix for child 2
    {
      // 0    1    2    3
      {0.5, 0.0, 0.5, 0.0}, // 0
      {0.0, 0.5, 0.5, 0.0}, // 1
      {0.0, 0.0, 1.0, 0.0}, // 2
      {0.0, 0.0, 0.5, 0.5}  // 3
    },

    // embedding matrix for child 3
    {
      // 0    1    2    3
      {0.5, 0.0, 0.0, 0.5}, // 0
      {0.0, 0.5, 0.0, 0.5}, // 1
      {0.0, 0.0, 0.5, 0.5}, // 2
      {0.0, 0.0, 0.0, 1.0}  // 3
    },

    // embedding matrix for child 4
    {
      // 0    1    2    3
      {0.5, 0.5, 0.0, 0.0}, // 0
      {0.0, 0.5, 0.0, 0.5}, // 1
      {0.5, 0.0, 0.5, 0.0}, // 2
      {0.5, 0.0, 0.0, 0.5}  // 3
    },

    // embedding matrix for child 5
    {
      // 0    1    2    3
      {0.5, 0.5, 0.0, 0.0}, // 0
      {0.0, 0.5, 0.5, 0.0}, // 1
      {0.5, 0.0, 0.5, 0.0}, // 2
      {0.0, 0.5, 0.0, 0.5}  // 3
    },

    // embedding matrix for child 6
    {
      // 0    1    2    3
      {0.5, 0.0, 0.5, 0.0}, // 0
      {0.0, 0.5, 0.5, 0.0}, // 1
      {0.0, 0.0, 0.5, 0.5}, // 2
      {0.0, 0.5, 0.0, 0.5}  // 3
    },

    // embedding matrix for child 7
    {
      // 0    1    2    3
      {0.5, 0.0, 0.5, 0.0}, // 0
      {0.0, 0.5, 0.0, 0.5}, // 1
      {0.0, 0.0, 0.5, 0.5}, // 2
      {0.5, 0.0, 0.0, 0.5}  // 3
    }
  };

#endif // #ifdef LIBMESH_ENABLE_AMR





Real Tet4::volume () const
{
  // The volume of a tetrahedron is 1/6 the box product formed
  // by its base and apex vectors
  Point a = point(3) - point(0);

  // b is the vector pointing from 0 to 1
  Point b = point(1) - point(0);

  // c is the vector pointing from 0 to 2
  Point c = point(2) - point(0);

  return triple_product(a, b, c) / 6.;
}




std::pair<Real, Real> Tet4::min_and_max_angle() const
{
  Point n[4];

  // Compute the outward normal vectors on each face
  n[0] = (this->point(2) - this->point(0)).cross(this->point(1) - this->point(0));
  n[1] = (this->point(1) - this->point(0)).cross(this->point(3) - this->point(0));
  n[2] = (this->point(2) - this->point(1)).cross(this->point(3) - this->point(1));
  n[3] = (this->point(0) - this->point(2)).cross(this->point(3) - this->point(2));

  Real dihedral_angles[6]; // 01, 02, 03, 12, 13, 23

  // Compute dihedral angles
  for (unsigned int k=0,i=0; i<4; ++i)
    for (unsigned int j=i+1; j<4; ++j,k+=1)
      dihedral_angles[k] = std::acos(n[i]*n[j] / n[i].norm() / n[j].norm()); // return value is between 0 and PI

  // Return max/min dihedral angles
  return std::make_pair(*std::min_element(dihedral_angles, dihedral_angles+6),
                        *std::max_element(dihedral_angles, dihedral_angles+6));

}



dof_id_type Tet4::key () const
{
  return this->compute_key(this->node_id(0),
                           this->node_id(1),
                           this->node_id(2),
                           this->node_id(3));
}



bool Tet4::contains_point (const Point & p, Real tol) const
{
  // See the response by Tony Noe on this thread.
  // http://bbs.dartmouth.edu/~fangq/MATH/download/source/Point_in_Tetrahedron.htm
  Point
    col1 = point(1) - point(0),
    col2 = point(2) - point(0),
    col3 = point(3) - point(0);

  Point r;

  libmesh_try
    {
      RealTensorValue(col1(0), col2(0), col3(0),
                      col1(1), col2(1), col3(1),
                      col1(2), col2(2), col3(2)).solve(p - point(0), r);
    }
  libmesh_catch (ConvergenceFailure &)
    {
      // The Jacobian was singular, therefore the Tet had
      // zero volume.  In this degenerate case, it is not
      // possible for the Tet to contain any points.
      return false;
    }

  // The point p is inside the tetrahedron if r1 + r2 + r3 < 1 and
  // 0 <= ri <= 1 for i = 1,2,3.
  return
    r(0) > -tol &&
    r(1) > -tol &&
    r(2) > -tol &&
    r(0) + r(1) + r(2) < 1.0 + tol;
}



#ifdef LIBMESH_ENABLE_AMR
float Tet4::embedding_matrix (const unsigned int i,
                              const unsigned int j,
                              const unsigned int k) const
{
  // Choose an optimal diagonal, if one has not already been selected
  this->choose_diagonal();

  // Permuted j and k indices
  unsigned int
    jp=j,
    kp=k;

  if ((i>3) && (this->_diagonal_selection!=DIAG_02_13))
    {
      // Just the enum value cast to an unsigned int...
      const unsigned ds = static_cast<unsigned int>(this->_diagonal_selection);

      // Permute j, k:
      // ds==1      ds==2
      // 0 -> 1     0 -> 2
      // 1 -> 2     1 -> 0
      // 2 -> 0     2 -> 1
      if (jp != 3)
        jp = (jp+ds)%3;

      if (kp != 3)
        kp = (kp+ds)%3;
    }

  // Debugging
  // libMesh::err << "Selected diagonal " << _diagonal_selection << std::endl;
  // libMesh::err << "j=" << j << std::endl;
  // libMesh::err << "k=" << k << std::endl;
  // libMesh::err << "jp=" << jp << std::endl;
  // libMesh::err << "kp=" << kp << std::endl;

  // Call embedding matrix with permuted indices
  return this->_embedding_matrix[i][jp][kp];
}



// void Tet4::reselect_diagonal (const Diagonal diag)
// {
//   /* Make sure that the element has just been refined.  */
//   libmesh_assert(_children);
//   libmesh_assert_equal_to (n_children(), 8);
//   libmesh_assert_equal_to (_children[0]->refinement_flag(), JUST_REFINED);
//   libmesh_assert_equal_to (_children[1]->refinement_flag(), JUST_REFINED);
//   libmesh_assert_equal_to (_children[2]->refinement_flag(), JUST_REFINED);
//   libmesh_assert_equal_to (_children[3]->refinement_flag(), JUST_REFINED);
//   libmesh_assert_equal_to (_children[4]->refinement_flag(), JUST_REFINED);
//   libmesh_assert_equal_to (_children[5]->refinement_flag(), JUST_REFINED);
//   libmesh_assert_equal_to (_children[6]->refinement_flag(), JUST_REFINED);
//   libmesh_assert_equal_to (_children[7]->refinement_flag(), JUST_REFINED);
//
//   /* Check whether anything has to be changed.  */
//   if (_diagonal_selection!=diag)
//     {
//       /* Set new diagonal selection.  */
//       _diagonal_selection = diag;
//
//       /* The first four children do not have to be changed.  For the
//  others, only the nodes have to be changed.  Note that we have
//  to keep track of the nodes ourselves since there is no \p
//  MeshRefinement object with a valid \p _new_nodes_map
//  available.  */
//       for (unsigned int c=4; c<this->n_children(); c++)
// {
//   Elem * child = this->child_ptr(c);
//   for (unsigned int nc=0; nc<child->n_nodes(); nc++)
//     {
//       /* Unassign the current node.  */
//       child->set_node(nc) = nullptr;
//
//       /* We have to find the correct new node now.  We know
//  that it exists somewhere.  We make use of the fact
//  that the embedding matrix for these children consists
//  of entries 0.0 and 0.5 only.  Also, we make use of
//  the properties of the embedding matrix for the first
//  (unchanged) four children, which allow us to use a
//  simple mechanism to find the required node.  */
//
//
//       unsigned int first_05_in_embedding_matrix = libMesh::invalid_uint;
//       for (unsigned int n=0; n<this->n_nodes(); n++)
// {
//   if (this->embedding_matrix(c,nc,n) != 0.0)
//     {
//       /* It must be 0.5 then.  Check whether it's the
//  first or second time that we get a 0.5
//  value.  */
//       if (first_05_in_embedding_matrix==libMesh::invalid_uint)
// {
//   /* First time, so just remember this position.  */
//   first_05_in_embedding_matrix = n;
// }
//       else
// {
//   /* Second time, so we know now which node to
//      use.  */
//   child->set_node(nc) = this->child_ptr(n)->node_ptr(first_05_in_embedding_matrix);
// }
//
//     }
// }
//
//       /* Make sure that a node has been found.  */
//       libmesh_assert(child->node_ptr(nc));
//     }
// }
//     }
// }



// void Tet4::reselect_optimal_diagonal (const Diagonal exclude_this)
// {
//   Real diag_01_23 = (this->point(0)+this->point(1)-this->point(2)-this->point(3)).norm_sq();
//   Real diag_02_13 = (this->point(0)-this->point(1)+this->point(2)-this->point(3)).norm_sq();
//   Real diag_03_12 = (this->point(0)-this->point(1)-this->point(2)+this->point(3)).norm_sq();
//
//   Diagonal use_this = INVALID_DIAG;
//   switch (exclude_this)
//     {
//     case DIAG_01_23:
//       use_this = DIAG_02_13;
//       if (diag_03_12 < diag_02_13)
// {
//   use_this = DIAG_03_12;
// }
//       break;
//
//     case DIAG_02_13:
//       use_this = DIAG_03_12;
//       if (diag_01_23 < diag_03_12)
// {
//   use_this = DIAG_01_23;
// }
//       break;
//
//     case DIAG_03_12:
//       use_this = DIAG_02_13;
//       if (diag_01_23 < diag_02_13)
// {
//   use_this = DIAG_01_23;
// }
//       break;
//
//     default:
//       use_this = DIAG_02_13;
//       if (diag_01_23 < diag_02_13 || diag_03_12 < diag_02_13)
// {
//   if (diag_01_23 < diag_03_12)
//     {
//       use_this = DIAG_01_23;
//     }
//   else
//     {
//       use_this = DIAG_03_12;
//     }
// }
//       break;
//     }
//
//   reselect_diagonal (use_this);
// }
#endif // #ifdef LIBMESH_ENABLE_AMR

} // namespace libMesh
