// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

// Local includes
#include "libmesh/face_inf_quad4.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_type.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_inf_edge2.h"
#include "libmesh/enum_io_package.h"
#include "libmesh/enum_order.h"
#include "libmesh/inf_fe_map.h"

namespace libMesh
{



// ------------------------------------------------------------
// InfQuad4 class static member initializations
const int InfQuad4::num_nodes;
const int InfQuad4::nodes_per_side;

const unsigned int InfQuad4::side_nodes_map[InfQuad4::num_sides][InfQuad4::nodes_per_side] =
  {
    {0, 1}, // Side 0
    {1, 3}, // Side 1
    {0, 2}  // Side 2
  };

#ifdef LIBMESH_ENABLE_AMR

const Real InfQuad4::_embedding_matrix[InfQuad4::num_children][InfQuad4::num_nodes][InfQuad4::num_nodes] =
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

bool InfQuad4::is_node_on_side(const unsigned int n,
                               const unsigned int s) const
{
  libmesh_assert_less (s, n_sides());
  return std::find(std::begin(side_nodes_map[s]),
                   std::end(side_nodes_map[s]),
                   n) != std::end(side_nodes_map[s]);
}

std::vector<unsigned>
InfQuad4::nodes_on_side(const unsigned int s) const
{
  libmesh_assert_less(s, n_sides());
  return {std::begin(side_nodes_map[s]), std::end(side_nodes_map[s])};
}

std::vector<unsigned>
InfQuad4::nodes_on_edge(const unsigned int e) const
{
  return nodes_on_side(e);
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
   * to the scheme of using inverse_map().
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
       * cannot say anything, fall back to the inverse_map()
       *
       * Declare a basic FEType.  Will use default in the base,
       * and something else (not important) in radial direction.
       */
      FEType fe_type(default_order());

      const Point mapped_point = InfFEMap::inverse_map(dim(), this, p,
                                                       tol, false);

      return this->on_reference_element(mapped_point, tol);
    }
}




Order InfQuad4::default_order() const
{
  return FIRST;
}



std::unique_ptr<Elem> InfQuad4::build_side_ptr (const unsigned int i)
{
  // libmesh_assert_less (i, this->n_sides());

  std::unique_ptr<Elem> edge;

  switch (i)
    {
    case 0:
      {
        edge = std::make_unique<Edge2>();
        break;
      }

      // adjacent to another infinite element
    case 1:
    case 2:
      {
        edge = std::make_unique<InfEdge2>();
        break;
      }

    default:
      libmesh_error_msg("Invalid side i = " << i);
    }

  // Set the nodes
  for (auto n : edge->node_index_range())
    edge->set_node(n, this->node_ptr(InfQuad4::side_nodes_map[i][n]));

  edge->set_interior_parent(this);
  edge->inherit_data_from(*this);

  return edge;
}



void InfQuad4::build_side_ptr (std::unique_ptr<Elem> & side,
                               const unsigned int i)
{
  this->side_ptr(side, i);
  side->set_interior_parent(this);
  side->inherit_data_from(*this);
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
        return;
      }
    default:
      libmesh_error_msg("Unsupported IO package " << iop);
    }
}


ElemType
InfQuad4::side_type (const unsigned int s) const
{
  libmesh_assert_less (s, 3);
  if (s == 0)
    return EDGE2;
  return INFEDGE2;
}


void InfQuad4::flip(BoundaryInfo * boundary_info)
{
  libmesh_assert(boundary_info);

  swap2nodes(0,1);
  swap2nodes(2,3);
  swap2neighbors(1,2);
  swap2boundarysides(1,2,boundary_info);
  swap2boundaryedges(1,2,boundary_info);
}

} // namespace libMesh


#endif // ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
