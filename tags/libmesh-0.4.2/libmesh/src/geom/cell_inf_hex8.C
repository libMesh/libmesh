// $Id: cell_inf_hex8.C,v 1.21 2004-01-03 15:37:43 benkirk Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2004  Benjamin S. Kirk, John W. Peterson
  
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
#include "cell_inf_hex8.h"
#include "fe_interface.h"
#include "fe_type.h"


// ------------------------------------------------------------
// InfHex8 class member functions

bool InfHex8::contains_point (const Point& p) const
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
   * Use size_sq() instead of size(), it is faster
   */
  const Real min_distance_sq = std::min((Point(this->point(0)-origin)).size_sq(),
					std::min((Point(this->point(1)-origin)).size_sq(),
						 std::min((Point(this->point(2)-origin)).size_sq(),
							  (Point(this->point(3)-origin)).size_sq())));

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






const std::vector<unsigned int> InfHex8::tecplot_connectivity(const unsigned int sc) const
{
  assert (_nodes != NULL);
  assert (sc < this->n_sub_elem());

  std::vector<unsigned int> conn(8);
  
  conn[0] = this->node(0)+1;
  conn[1] = this->node(1)+1;
  conn[2] = this->node(2)+1;
  conn[3] = this->node(3)+1;
  conn[4] = this->node(4)+1;
  conn[5] = this->node(5)+1;
  conn[6] = this->node(6)+1;
  conn[7] = this->node(7)+1;

  return conn;
}




#ifdef ENABLE_AMR

const float InfHex8::_embedding_matrix[4][8][8] =
{
  // embedding matrix for child 0
  {
    //     0      1      2      3      4      5      6      7 th parent N.(ode)
    {    1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}, // 0th child N.
    {    0.5,   0.5,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}, // 1
    {   0.25,  0.25,  0.25,  0.25,   0.0,   0.0,   0.0,   0.0}, // 2
    {    0.5,   0.0,   0.0,   0.5,   0.0,   0.0,   0.0,   0.0}, // 3
    {    0.0,   0.0,   0.0,   0.0,   1.0,   0.0,   0.0,   0.0}, // 4
    {    0.0,   0.0,   0.0,   0.0,   0.5,   0.5,   0.0,   0.0}, // 5
    {    0.0,   0.0,   0.0,   0.0,  0.25,  0.25,  0.25,  0.25}, // 6
    {    0.0,   0.0,   0.0,   0.0,   0.5,   0.0,   0.0,   0.5}  // 7
  },

  // embedding matrix for child 1
  {
    //     0      1      2      3      4      5      6      7 th parent N.(ode)
    {    0.5,   0.5,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}, // 0th child N.
    {    0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0}, // 1
    {    0.0,   0.5,   0.5,   0.0,   0.0,   0.0,   0.0,   0.0}, // 2
    {   0.25,  0.25,  0.25,  0.25,   0.0,   0.0,   0.0,   0.0}, // 3
    {    0.0,   0.0,   0.0,   0.0,   0.5,   0.5,   0.0,   0.0}, // 4
    {    0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0,   0.0}, // 5
    {    0.0,   0.0,   0.0,   0.0,   0.0,   0.5,   0.5,   0.0}, // 6
    {    0.0,   0.0,   0.0,   0.0,  0.25,  0.25,  0.25,  0.25}  // 7
  },

  // embedding matrix for child 2
  {
    //     0      1      2      3      4      5      6      7 th parent N.(ode)
    {    0.5,   0.0,   0.0,   0.5,   0.0,   0.0,   0.0,   0.0}, // 0th child N.
    {   0.25,  0.25,  0.25,  0.25,   0.0,   0.0,   0.0,   0.0}, // 1
    {    0.0,   0.0,   0.5,   0.5,   0.0,   0.0,   0.0,   0.0}, // 2
    {    0.0,   0.0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0}, // 3
    {    0.0,   0.0,   0.0,   0.0,   0.5,   0.0,   0.0,   0.5}, // 4
    {    0.0,   0.0,   0.0,   0.0,  0.25,  0.25,  0.25,  0.25}, // 5
    {    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.5,   0.5}, // 6
    {    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0}  // 7
  },

  // embedding matrix for child 3
  {
    //     0      1      2      3      4      5      6      7 th parent N.(ode)
    {   0.25,  0.25,  0.25,  0.25,   0.0,   0.0,   0.0,   0.0}, // 0th child N.
    {    0.0,   0.5,   0.5,   0.0,   0.0,   0.0,   0.0,   0.0}, // 1
    {    0.0,   0.0,   1.0,   0.0,   0.0,   0.0,   0.0,   0.0}, // 2
    {    0.0,   0.0,   0.5,   0.5,   0.0,   0.0,   0.0,   0.0}, // 3
    {    0.0,   0.0,   0.0,   0.0,  0.25,  0.25,  0.25,  0.25}, // 4
    {    0.0,   0.0,   0.0,   0.0,   0.0,   0.5,   0.5,   0.0}, // 5
    {    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   1.0,   0.0}, // 6
    {    0.0,   0.0,   0.0,   0.0,   0.0,   0.0,   0.5,   0.5}  // 7
  }
};



#endif

#endif // ifdef ENABLE_INFINITE_ELEMENTS
