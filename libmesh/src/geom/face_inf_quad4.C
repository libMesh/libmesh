// $Id: face_inf_quad4.C,v 1.16 2003-04-18 15:46:30 spetersen Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include "mesh_config.h"
#ifdef ENABLE_INFINITE_ELEMENTS


// Local includes cont'd
#include "face_inf_quad4.h"
#include "fe_interface.h"
#include "fe_type.h"



// ------------------------------------------------------------
// InfQuad4 class member functions
bool InfQuad4::contains_point (const Point& p) const
{
  /*
   * make use of the fact that infinite elements do not
   * live inside the envelope.  Use a fast scheme to
   * check whether point \p p is inside or outside
   * our relevant part of the envelope.  Note that
   * this is not exclusive: the point may be outside
   * the envelope, but contained in another infinite element.
   * Therefore, if the distance is greater, do fall back
   * to the scheme of using FEInterface::inverse_map().
   */
  const Point origin (this->origin());

  /*
   * determine the minimal distance of the base from the origin
   */
  const Real min_distance = std::min((Point(this->point(0)-origin)).size(),
				     (Point(this->point(1)-origin)).size());

  /*
   * work with 1% allowable deviation.  Can still fall
   * back to the InfFE::inverse_map()
   */
  const Real conservative_p_dist = 1.01 * (Point(p-origin).size());

  if (conservative_p_dist < min_distance)
    {
      /*
       * the physical point is definitely not contained
       * in the element, return false.
       */
      return false;
    }
  else
    {
      /*
       * cannot say anything, fall back to the FEInterface::inverse_map()
       *
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




#ifdef ENABLE_AMR

const float InfQuad4::_embedding_matrix[2][4][4] =
{
  // embedding matrix for child 0
  {
    // 0    1    2    3rd parent node
    {1.0, 0.0, 0.0, 0.0}, // 0th child node
    {0.5, 0.5, 0.0, 0.0}, // 1
    {0.0, 0.0, 1.0, 0.0}, // 2
    {0.0, 0.0, 0.5, 0.5}  // 3
  },

  // embedding matrix for child 1
  {
    // 0    1    2    3
    {0.5, 0.5, 0.0, 0.0}, // 0
    {0.0, 1.0, 0.0, 0.0}, // 1
    {0.0, 0.0, 0.5, 0.5}, // 2
    {0.0, 0.0, 0.0, 1.0}  // 3
  }
};

#endif



const std::vector<unsigned int> InfQuad4::tecplot_connectivity(const unsigned int sf) const
{
  assert (sf < this->n_sub_elem());
  
  std::vector<unsigned int> conn(4);

  conn[0] = this->node(0)+1;
  conn[1] = this->node(1)+1;
  conn[2] = this->node(3)+1;
  conn[3] = this->node(2)+1;

  return conn;
}



void InfQuad4::vtk_connectivity(const unsigned int sf,
				std::vector<unsigned int> *conn) const
{
  assert (_nodes != NULL);
  assert (sf < this->n_sub_elem());
  
  if (conn == NULL)
    conn = new std::vector<unsigned int>;

  conn->resize(4);

  (*conn)[0] = this->node(0);
  (*conn)[1] = this->node(1);
  (*conn)[2] = this->node(3);
  (*conn)[3] = this->node(2);

  return;
}


#endif // ifdef ENABLE_INFINITE_ELEMENTS
