// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS


// Local includes cont'd
#include "face_inf_quad4.h"
#include "fe_interface.h"
#include "fe_type.h"
#include "side.h"
#include "edge_edge2.h"
#include "edge_inf_edge2.h"

namespace libMesh
{



// ------------------------------------------------------------
// InfQuad4 class static member initializations
const unsigned int InfQuad4::side_nodes_map[3][2] =
{
  {0, 1}, // Side 0
  {1, 3}, // Side 1
  {0, 2}  // Side 2
};



#ifdef LIBMESH_ENABLE_AMR

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


// ------------------------------------------------------------
// InfQuad4 class member functions

bool InfQuad4::is_vertex(const unsigned int i) const
{
  if (i < 2)
    return true;
  return false;
}

bool InfQuad4::is_edge(const unsigned int i) const
{
  if (i < 2)
    return false;
  return true;
}

bool InfQuad4::is_face(const unsigned int) const
{
  return false;
}

bool InfQuad4::is_node_on_side(const unsigned int n,
			       const unsigned int s) const
{
  libmesh_assert(s < n_sides());
  for (unsigned int i = 0; i != 2; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool InfQuad4::contains_point (const Point& p, Real tol) const
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
   * use size_sq() instead of size(), it is slightly faster
   */
  const Real min_distance_sq = std::min((Point(this->point(0)-origin)).size_sq(),
					(Point(this->point(1)-origin)).size_sq());

  /*
   * work with 1% allowable deviation.  Can still fall
   * back to the InfFE::inverse_map()
   */
  const Real conservative_p_dist_sq = 1.01 * (Point(p-origin).size_sq());

  if (conservative_p_dist_sq < min_distance_sq)
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
							  tol,
							  false);

      return FEInterface::on_reference_element(mapped_point, this->type(), tol);
    }
}




AutoPtr<Elem> InfQuad4::build_side (const unsigned int i,
				    bool proxy) const
{
  // libmesh_assert (i < this->n_sides());

  if (proxy)
    {
      switch (i)
	{
	  // base
	case 0:
	  {
	    AutoPtr<Elem> ap(new Side<Edge2,InfQuad4>(this,i));
	    return ap;
	  }
	  // ifem edges
	case 1:
	case 2:
	  {
	    AutoPtr<Elem> ap(new Side<InfEdge2,InfQuad4>(this,i));
	    return ap;
	  }

	default:
	  libmesh_error();
	}
    }

  else
    {
      // FIXME: Find out how to return non-proxy side
      libmesh_error();
    }

  // How did we get here
  libmesh_error();
  AutoPtr<Elem> ap(NULL);  return ap;
}


void InfQuad4::connectivity(const unsigned int libmesh_dbg_var(sf),
			    const IOPackage iop,
			    std::vector<unsigned int>& conn) const
{
  libmesh_assert (sf < this->n_sub_elem());
  libmesh_assert (iop != INVALID_IO_PACKAGE);

  conn.resize(4);

  switch (iop)
    {
    case TECPLOT:
      {
	conn[0] = this->node(0)+1;
	conn[1] = this->node(1)+1;
	conn[2] = this->node(3)+1;
	conn[3] = this->node(2)+1;
	return;
      }
    case VTK:
      {
	conn[0] = this->node(0);
	conn[1] = this->node(1);
	conn[2] = this->node(3);
	conn[3] = this->node(2);
      }
    default:
      libmesh_error();
    }

  libmesh_error();
}

} // namespace libMesh


#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS


