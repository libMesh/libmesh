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
#include "libmesh/cell_tet4.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/face_tri3.h"

namespace libMesh
{



// ------------------------------------------------------------
// Tet4 class static member initializations
const unsigned int Tet4::side_nodes_map[4][3] =
  {
    {0, 2, 1}, // Side 0
    {0, 1, 3}, // Side 1
    {1, 2, 3}, // Side 2
    {2, 0, 3}  // Side 3
  };

const unsigned int Tet4::edge_nodes_map[6][2] =
  {
    {0, 1}, // Side 0
    {1, 2}, // Side 1
    {0, 2}, // Side 2
    {0, 3}, // Side 3
    {1, 3}, // Side 4
    {2, 3}  // Side 5
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
  for (unsigned int i = 0; i != 2; ++i)
    if (edge_nodes_map[e][i] == n)
      return true;
  return false;
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
  for (unsigned int i = 0; i != 3; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

AutoPtr<Elem> Tet4::build_side (const unsigned int i,
                                bool proxy) const
{
  libmesh_assert_less (i, this->n_sides());

  if (proxy)
    return AutoPtr<Elem>(new Side<Tri3,Tet4>(this,i));

  else
    {
      Elem* face = new Tri3;
      face->subdomain_id() = this->subdomain_id();

      switch (i)
        {
        case 0:
          {
            face->set_node(0) = this->get_node(0);
            face->set_node(1) = this->get_node(2);
            face->set_node(2) = this->get_node(1);
            break;
          }
        case 1:
          {
            face->set_node(0) = this->get_node(0);
            face->set_node(1) = this->get_node(1);
            face->set_node(2) = this->get_node(3);
            break;
          }
        case 2:
          {
            face->set_node(0) = this->get_node(1);
            face->set_node(1) = this->get_node(2);
            face->set_node(2) = this->get_node(3);
            break;
          }
        case 3:
          {
            face->set_node(0) = this->get_node(2);
            face->set_node(1) = this->get_node(0);
            face->set_node(2) = this->get_node(3);
            break;
          }
        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      return AutoPtr<Elem>(face);
    }

  libmesh_error_msg("We'll never get here!");
  return AutoPtr<Elem>();
}


AutoPtr<Elem> Tet4::build_edge (const unsigned int i) const
{
  libmesh_assert_less (i, this->n_edges());

  return AutoPtr<Elem>(new SideEdge<Edge2,Tet4>(this,i));
}


void Tet4::connectivity(const unsigned int libmesh_dbg_var(sc),
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
        conn[5] = this->node(3)+1;
        conn[6] = this->node(3)+1;
        conn[7] = this->node(3)+1;
        return;
      }

    case VTK:
      {
        conn.resize(4);
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



#ifdef LIBMESH_ENABLE_AMR

const float Tet4::_embedding_matrix[8][4][4] =
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
  Point a ( *this->get_node(3) - *this->get_node(0) );

  // b is the vector pointing from 0 to 1
  Point b ( *this->get_node(1) - *this->get_node(0) );

  // c is the vector pointing from 0 to 2
  Point c ( *this->get_node(2) - *this->get_node(0) );

  return (1.0 / 6.0) * (a * (b.cross(c)));
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
      dihedral_angles[k] = std::acos(n[i]*n[j] / n[i].size() / n[j].size()); // return value is between 0 and PI

  // Return max/min dihedral angles
  return std::make_pair(*std::min_element(dihedral_angles, dihedral_angles+6),
                        *std::max_element(dihedral_angles, dihedral_angles+6));

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

  // Call embedding matrx with permuted indices
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
//   Elem *child = this->child(c);
//   for (unsigned int nc=0; nc<child->n_nodes(); nc++)
//     {
//       /* Unassign the current node.  */
//       child->set_node(nc) = NULL;
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
//   /* First time, so just remeber this position.  */
//   first_05_in_embedding_matrix = n;
// }
//       else
// {
//   /* Second time, so we know now which node to
//      use.  */
//   child->set_node(nc) = this->child(n)->get_node(first_05_in_embedding_matrix);
// }
//
//     }
// }
//
//       /* Make sure that a node has been found.  */
//       libmesh_assert(child->get_node(nc));
//     }
// }
//     }
// }



// void Tet4::reselect_optimal_diagonal (const Diagonal exclude_this)
// {
//   Real diag_01_23 = (this->point(0)+this->point(1)-this->point(2)-this->point(3)).size_sq();
//   Real diag_02_13 = (this->point(0)-this->point(1)+this->point(2)-this->point(3)).size_sq();
//   Real diag_03_12 = (this->point(0)-this->point(1)-this->point(2)+this->point(3)).size_sq();
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
