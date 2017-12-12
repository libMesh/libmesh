// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS


// Local includes cont'd
#include "libmesh/face_inf_quad4.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_type.h"
#include "libmesh/side.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_inf_edge2.h"

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
  libmesh_assert_less (s, n_sides());
  for (unsigned int i = 0; i != 2; ++i)
    if (side_nodes_map[s][i] == n)
      return true;
  return false;
}

bool InfQuad4::contains_point (const Point & p, Real tol) const
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
  const Point my_origin (this->origin());

  /*
   * determine the minimal distance of the base from the origin
   * use norm_sq() instead of norm(), it is slightly faster
   */
  const Real min_distance_sq = std::min((Point(this->point(0)-my_origin)).norm_sq(),
                                        (Point(this->point(1)-my_origin)).norm_sq());

  /*
   * work with 1% allowable deviation.  Can still fall
   * back to the InfFE::inverse_map()
   */
  const Real conservative_p_dist_sq = 1.01 * (Point(p-my_origin).norm_sq());

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




std::unique_ptr<Elem> InfQuad4::build_side_ptr (const unsigned int i,
                                                bool proxy)
{
  // libmesh_assert_less (i, this->n_sides());

  if (proxy)
    {
      switch (i)
        {
          // base
        case 0:
          return libmesh_make_unique<Side<Edge2,InfQuad4>>(this,i);

          // ifem edges
        case 1:
        case 2:
          return libmesh_make_unique<Side<InfEdge2,InfQuad4>>(this,i);

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }
    }

  else
    {
      // Return value
      std::unique_ptr<Elem> edge;

      switch (i)
        {
        case 0:
          {
            edge = libmesh_make_unique<Edge2>();
            break;
          }

          // adjacent to another infinite element
        case 1:
        case 2:
          {
            edge = libmesh_make_unique<InfEdge2>();
            break;
          }

        default:
          libmesh_error_msg("Invalid side i = " << i);
        }

      edge->subdomain_id() = this->subdomain_id();

      // Set the nodes
      for (unsigned n=0; n<edge->n_nodes(); ++n)
        edge->set_node(n) = this->node_ptr(InfQuad4::side_nodes_map[i][n]);

      return edge;
    }
}


void InfQuad4::connectivity(const unsigned int libmesh_dbg_var(sf),
                            const IOPackage iop,
                            std::vector<dof_id_type> & conn) const
{
  libmesh_assert_less (sf, this->n_sub_elem());
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  conn.resize(4);

  switch (iop)
    {
    case TECPLOT:
      {
        conn[0] = this->node_id(0)+1;
        conn[1] = this->node_id(1)+1;
        conn[2] = this->node_id(3)+1;
        conn[3] = this->node_id(2)+1;
        return;
      }
    case VTK:
      {
        conn[0] = this->node_id(0);
        conn[1] = this->node_id(1);
        conn[2] = this->node_id(3);
        conn[3] = this->node_id(2);
      }
    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}

} // namespace libMesh


#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
