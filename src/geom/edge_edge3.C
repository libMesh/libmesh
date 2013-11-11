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



// Local includes
#include "libmesh/edge_edge3.h"

namespace libMesh
{

#ifdef LIBMESH_ENABLE_AMR

const float Edge3::_embedding_matrix[2][3][3] =
{
  // embedding matrix for child 0
  {
    // 0    1    2
    {1.0, 0.0, 0.0}, // left
    {0.0, 0.0, 1.0}, // right
    {0.375,-0.125,0.75} // middle
  },

  // embedding matrix for child 1
  {
    // 0    1    2
    {0.0, 0.0, 1.0}, // left
    {0.0, 1.0, 0.0},  // right
    {-0.125,0.375,0.75} // middle
  }
};

#endif

bool Edge3::is_vertex(const unsigned int i) const
{
  if (i < 2)
    return true;
  return false;
}

bool Edge3::is_edge(const unsigned int i) const
{
  if (i < 2)
   return false;
  return true;
}

bool Edge3::is_face(const unsigned int ) const
{
  return false;
}

bool Edge3::is_node_on_side(const unsigned int n,
			    const unsigned int s) const
{
  libmesh_assert_less (s, 2);
  libmesh_assert_less (n, 3);
  return (s == n);
}

bool Edge3::is_node_on_edge(const unsigned int,
			    const unsigned int libmesh_dbg_var(e)) const
{
  libmesh_assert_equal_to (e, 0);
  return true;
}



bool Edge3::has_affine_map() const
{
  return (this->point(2).relative_fuzzy_equals
          ((this->point(0) + this->point(1))/2));
}



void Edge3::connectivity(const unsigned int sc,
			 const IOPackage iop,
			 std::vector<dof_id_type>& conn) const
{
  libmesh_assert_less_equal (sc, 1);
  libmesh_assert_less (sc, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  // Create storage
  conn.resize(2);

  switch (iop)
    {
    case TECPLOT:
      {
	switch (sc)
	  {
	  case 0:
	    conn[0] = this->node(0)+1;
	    conn[1] = this->node(2)+1;
	    return;

	  case 1:
	    conn[0] = this->node(2)+1;
	    conn[1] = this->node(1)+1;
	    return;

	  default:
	    libmesh_error();
	  }
      }


    case VTK:
      {
        conn.resize(3);
        conn[0] = this->node(0);
        conn[1] = this->node(1);
        conn[2] = this->node(2);
        return;

        /*
	switch (sc)
	  {
	  case 0:
	    conn[0] = this->node(0);
	    conn[1] = this->node(2);

	    return;

	  case 1:
	    conn[0] = this->node(2);
	    conn[1] = this->node(1);

	    return;

	  default:
	    libmesh_error();
          }
        */
      }

    default:
      {
	libmesh_error();
      }
    }
}



std::pair<unsigned short int, unsigned short int>
Edge3::second_order_child_vertex (const unsigned int) const
{
  return std::pair<unsigned short int, unsigned short int>(0,0);
}



Real Edge3::volume () const
{
  // Finding the (exact) length of a general quadratic element
  // is a surprisingly complicated formula.
  Point A = this->point(0) + this->point(1) - 2*this->point(2);
  Point B = (this->point(1) - this->point(0))/2;

  const Real a = A.size_sq();
  const Real b = 2.*(A*B);
  const Real c = B.size_sq();

  // Degenerate straight line case
  if (a < TOLERANCE*TOLERANCE*TOLERANCE)
    return (this->point(1) - this->point(0)).size();

  const Real ba=b/a;
  const Real ca=c/a;

  libmesh_assert (1.-ba+ca>0.);

  const Real s1 = std::sqrt(1. - ba + ca);
  const Real s2 = std::sqrt(1. + ba + ca);

  return 0.5*std::sqrt(a)*((1.-0.5*ba)*s1 +
			   (1.+0.5*ba)*s2 +
			   (ca - 0.25*ba*ba)*std::log( (1.-0.5*ba+s1)/(-1.-0.5*ba+s2) )
			   );
}

} // namespace libMesh
