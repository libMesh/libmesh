// $Id: cell_inf_prism6.C,v 1.27 2005-03-01 14:21:40 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include "libmesh_config.h"

#ifdef ENABLE_INFINITE_ELEMENTS

// C++ includes

// Local includes cont'd
#include "cell_inf_prism6.h"
#include "fe_interface.h"
#include "fe_type.h"



// ------------------------------------------------------------
// InfPrism6 class static member initializations
const unsigned int InfPrism6::side_nodes_map[4][4] =
{
  { 0, 1, 2, 99}, // Side 0
  { 0, 1, 3, 4},  // Side 1
  { 1, 2, 4, 5},  // Side 2
  { 2, 0, 5, 3}   // Side 3
};


// ------------------------------------------------------------
// InfPrism6 class member functions

bool InfPrism6::is_vertex(const unsigned int i) const
{
  if (i < 3)
    return true;
  return false;
}

bool InfPrism6::is_edge(const unsigned int i) const
{
  if (i < 3)
    return false;
  return true;
}

bool InfPrism6::is_face(const unsigned int) const
{
  return false;
}

bool InfPrism6::is_node_on_side(const unsigned int n,
				const unsigned int s) const
{
  assert(s < n_sides());
  for (unsigned int i = 0; i != 4; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool InfPrism6::contains_point (const Point& p) const
{
  /*
   * For infinite elements with linear base interpolation:
   *
   * make use of the fact that infinite elements do not
   * live inside the envelope.  Use a fast scheme to
   * check whether point \p p is inside or outside
   * our relevant part of the envelope.  Note that
   * this is not exclusive: only when the distance is less,
   * we are safe.  Otherwise, we cannot say anything. The 
   * envelope may be non-spherical, the physical point may lie
   * inside the envelope, outside the envelope, or even inside 
   * this infinite element.  Therefore if this fails,
   * fall back to the FEInterface::inverse_map()
   */
  const Point origin (this->origin());

  /*
   * determine the minimal distance of the base from the origin
   * use size_sq(), it is faster than size() and produces
   * the same behavior
   */
  const Real min_distance_sq = std::min((Point(this->point(0)-origin)).size_sq(),
				     std::min((Point(this->point(1)-origin)).size_sq(),
					      (Point(this->point(2)-origin)).size_sq()));

  /*
   * work with 1% allowable deviation.  We can still fall
   * back to the InfFE::inverse_map()
   */
  const Real conservative_p_dist_sq = 1.01 * (Point(p-origin).size_sq());


  if (conservative_p_dist_sq < min_distance_sq)
    {
      /*
       * the physical point is definitely not contained in the element
       */
      return false;
    }
  else
    {
      /*
       * Declare a basic FEType.  Will use default in the base,
       * and something else (not important) in radial direction.
       */
      FEType fe_type(default_order());
  
      const Point mapped_point = FEInterface::inverse_map(dim(),
							  fe_type,
							  this,
							  p,
							  1.e-4,
							  false);

      return FEInterface::on_reference_element(mapped_point, this->type());
    }
}



void InfPrism6::connectivity(const unsigned int sc,
			     const IOPackage iop,
			     std::vector<unsigned int>& conn) const
{
  assert (_nodes != NULL);
  assert (sc < this->n_sub_elem());
  assert (iop != INVALID_IO_PACKAGE);

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

    default:
      error();
    }

  error();
}





#ifdef ENABLE_AMR

const float InfPrism6::_embedding_matrix[4][6][6] =
{
  // embedding matrix for child 0
  {
    //          0           1           2           3           4           5 th parent Node
    {         1.0,        0.0,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 1
    {         0.5,        0.0,        0.5,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        1.0,        0.0,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.5,        0.5,        0.0}, // 4
    {         0.0,        0.0,        0.0,        0.5,        0.0,        0.5}  // 5
  },

  // embedding matrix for child 1
  {
    //          0           1           2           3           4           5 th parent Node
    {         0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.0,        1.0,        0.0,        0.0,        0.0,        0.0}, // 1
    {         0.0,        0.5,        0.5,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        0.5,        0.5,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.0,        1.0,        0.0}, // 4
    {         0.0,        0.0,        0.0,        0.0,        0.5,        0.5}  // 5
  },

  // embedding matrix for child 2
  {
    //          0           1           2           3           4           5 th parent Node
    {         0.5,        0.0,        0.5,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.0,        0.5,        0.5,        0.0,        0.0,        0.0}, // 1
    {         0.0,        0.0,        1.0,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        0.5,        0.0,        0.5}, // 3
    {         0.0,        0.0,        0.0,        0.0,        0.5,        0.5}, // 4
    {         0.0,        0.0,        0.0,        0.0,        0.0,        1.0}  // 5
  },

  // embedding matrix for child 3
  {
    //          0           1           2           3           4           5 th parent Node
    {         0.5,        0.5,        0.0,        0.0,        0.0,        0.0}, // 0th child N.
    {         0.0,        0.5,        0.5,        0.0,        0.0,        0.0}, // 1
    {         0.5,        0.0,        0.5,        0.0,        0.0,        0.0}, // 2
    {         0.0,        0.0,        0.0,        0.5,        0.5,        0.0}, // 3
    {         0.0,        0.0,        0.0,        0.0,        0.5,        0.5}, // 4
    {         0.0,        0.0,        0.0,        0.5,        0.0,        0.5}  // 5
  }
};



#endif

#endif  // ifdef ENABLE_INFINITE_ELEMENTS
