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



// C++ includes
#include <algorithm> // for std::sort
#include <iterator>  // for std::ostream_iterator
#include <sstream>
#include <limits>    // for std::numeric_limits<>
#include <cmath>     // for std::sqrt()

// Local includes
#include "libmesh/elem.h"
#include "libmesh/fe_type.h"
#include "libmesh/fe_interface.h"
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_edge3.h"
#include "libmesh/edge_edge4.h"
#include "libmesh/edge_inf_edge2.h"
#include "libmesh/face_tri3.h"
#include "libmesh/face_tri6.h"
#include "libmesh/face_quad4.h"
#include "libmesh/face_quad8.h"
#include "libmesh/face_quad9.h"
#include "libmesh/face_inf_quad4.h"
#include "libmesh/face_inf_quad6.h"
#include "libmesh/cell_tet4.h"
#include "libmesh/cell_tet10.h"
#include "libmesh/cell_hex8.h"
#include "libmesh/cell_hex20.h"
#include "libmesh/cell_hex27.h"
#include "libmesh/cell_inf_hex8.h"
#include "libmesh/cell_inf_hex16.h"
#include "libmesh/cell_inf_hex18.h"
#include "libmesh/cell_prism6.h"
#include "libmesh/cell_prism15.h"
#include "libmesh/cell_prism18.h"
#include "libmesh/cell_inf_prism6.h"
#include "libmesh/cell_inf_prism12.h"
#include "libmesh/cell_pyramid5.h"
#include "libmesh/fe_base.h"
#include "libmesh/mesh_base.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/remote_elem.h"
#include "libmesh/reference_elem.h"
#include "libmesh/string_to_enum.h"

#ifdef LIBMESH_ENABLE_PERIODIC
#include "libmesh/mesh.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/boundary_info.h"
#endif

namespace libMesh
{

// Initialize static member variables
const unsigned int Elem::type_to_n_nodes_map [] =
  {
    2,  // EDGE2
    3,  // EDGE3
    4,  // EDGE4

    3,  // TRI3
    6,  // TRI6

    4,  // QUAD4
    8,  // QUAD8
    9,  // QUAD9

    4,  // TET4
    10, // TET10

    8,  // HEX8
    20, // HEX20
    27, // HEX27

    6,  // PRISM6
    15, // PRISM15
    18, // PRISM18

    5,  // PYRAMID5

    2,  // INFEDGE2

    4,  // INFQUAD4
    6,  // INFQUAD6

    8,  // INFHEX8
    16, // INFHEX16
    18, // INFHEX18

    6,  // INFPRISM6
    16, // INFPRISM12

    1,  // NODEELEM
  };

const unsigned int Elem::type_to_n_sides_map [] =
  {
    2,  // EDGE2
    2,  // EDGE3
    2,  // EDGE4

    3,  // TRI3
    3,  // TRI6

    4,  // QUAD4
    4,  // QUAD8
    4,  // QUAD9

    4,  // TET4
    4,  // TET10

    6,  // HEX8
    6,  // HEX20
    6,  // HEX27

    5,  // PRISM6
    5,  // PRISM15
    5,  // PRISM18

    5,  // PYRAMID5

    2,  // INFEDGE2

    3,  // INFQUAD4
    3,  // INFQUAD6

    5,  // INFHEX8
    5,  // INFHEX16
    5,  // INFHEX18

    4,  // INFPRISM6
    4,  // INFPRISM12

    0,  // NODEELEM
  };


// ------------------------------------------------------------
// Elem class member funcions
AutoPtr<Elem> Elem::build(const ElemType type,
			  Elem* p)
{
  Elem* elem = NULL;

  switch (type)
    {
      // 1D elements
    case EDGE2:
      {
	elem = new Edge2(p);
	break;
      }
    case EDGE3:
      {
	elem = new Edge3(p);
	break;
      }
    case EDGE4:
      {
        elem = new Edge4(p);
        break;
      }



      // 2D elements
    case TRI3:
      {
	elem = new Tri3(p);
	break;
      }
    case TRI6:
      {
	elem = new Tri6(p);
	break;
      }
    case QUAD4:
      {
	elem = new Quad4(p);
	break;
      }
    case QUAD8:
      {
	elem = new Quad8(p);
	break;
      }
    case QUAD9:
      {
	elem = new Quad9(p);
	break;
      }


      // 3D elements
    case TET4:
      {
	elem = new Tet4(p);
	break;
      }
    case TET10:
      {
	elem = new Tet10(p);
	break;
      }
    case HEX8:
      {
	elem = new Hex8(p);
	break;
      }
    case HEX20:
      {
	elem = new Hex20(p);
	break;
      }
    case HEX27:
      {
	elem = new Hex27(p);
	break;
      }
    case PRISM6:
      {
	elem = new Prism6(p);
	break;
      }
    case PRISM15:
      {
	elem = new Prism15(p);
	break;
      }
    case PRISM18:
      {
	elem = new Prism18(p);
	break;
      }
    case PYRAMID5:
      {
	elem = new Pyramid5(p);
	break;
      }



#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

      // 1D infinite elements
    case INFEDGE2:
      {
	elem = new InfEdge2(p);
	break;
      }


      // 2D infinite elements
    case INFQUAD4:
      {
	elem = new InfQuad4(p);
	break;
      }
    case INFQUAD6:
      {
	elem = new InfQuad6(p);
	break;
      }


    // 3D infinite elements
    case INFHEX8:
      {
	elem = new InfHex8(p);
	break;
      }
    case INFHEX16:
      {
	elem = new InfHex16(p);
	break;
      }
    case INFHEX18:
      {
	elem = new InfHex18(p);
	break;
      }
    case INFPRISM6:
      {
	elem = new InfPrism6(p);
	break;
      }
    case INFPRISM12:
      {
	elem = new InfPrism12(p);
	break;
      }

#endif

    default:
      {
	libMesh::err << "ERROR: Undefined element type!." << std::endl;
	libmesh_error();
      }
    }


  AutoPtr<Elem> ap(elem);
  return ap;
}



const Elem* Elem::reference_elem () const
{
  return &(ReferenceElem::get(this->type()));
}



Point Elem::centroid() const
{
  Point cp;

  for (unsigned int n=0; n<this->n_vertices(); n++)
    cp.add (this->point(n));

  return (cp /= static_cast<Real>(this->n_vertices()));
}



Real Elem::hmin() const
{
  Real h_min=std::numeric_limits<Real>::max();

  for (unsigned int n_outer=0; n_outer<this->n_vertices(); n_outer++)
    for (unsigned int n_inner=n_outer+1; n_inner<this->n_vertices(); n_inner++)
      {
	const Point diff = (this->point(n_outer) - this->point(n_inner));

	h_min = std::min(h_min,diff.size_sq());
      }

  return std::sqrt(h_min);
}



Real Elem::hmax() const
{
  Real h_max=0;

  for (unsigned int n_outer=0; n_outer<this->n_vertices(); n_outer++)
    for (unsigned int n_inner=n_outer+1; n_inner<this->n_vertices(); n_inner++)
      {
	const Point diff = (this->point(n_outer) - this->point(n_inner));

	h_max = std::max(h_max,diff.size_sq());
      }

  return std::sqrt(h_max);
}



Real Elem::length(const unsigned int n1,
		  const unsigned int n2) const
{
  libmesh_assert_less ( n1, this->n_vertices() );
  libmesh_assert_less ( n2, this->n_vertices() );

  return (this->point(n1) - this->point(n2)).size();
}



bool Elem::operator == (const Elem& rhs) const
{

    // Cast rhs to an Elem*
//    const Elem* rhs_elem = dynamic_cast<const Elem*>(&rhs);
    const Elem* rhs_elem = &rhs;

    // If we cannot cast to an Elem*, rhs must be a Node
//    if(rhs_elem == static_cast<const Elem*>(NULL))
//        return false;

//   libmesh_assert (n_nodes());
//   libmesh_assert (rhs.n_nodes());

//   // Elements can only be equal if they
//   // contain the same number of nodes.
//   if (this->n_nodes() == rhs.n_nodes())
//     {
//       // Create a set that contains our global
//       // node numbers and those of our neighbor.
//       // If the set is the same size as the number
//       // of nodes in both elements then they must
//       // be connected to the same nodes.
//       std::set<unsigned int> nodes_set;

//       for (unsigned int n=0; n<this->n_nodes(); n++)
//         {
//           nodes_set.insert(this->node(n));
//           nodes_set.insert(rhs.node(n));
//         }

//       // If this passes the elements are connected
//       // to the same global nodes
//       if (nodes_set.size() == this->n_nodes())
//         return true;
//     }

//   // If we get here it is because the elements either
//   // do not have the same number of nodes or they are
//   // connected to different nodes.  Either way they
//   // are not the same element
//   return false;

  // Useful typedefs
  typedef std::vector<dof_id_type>::iterator iterator;


  // Elements can only be equal if they
  // contain the same number of nodes.
  // However, we will only test the vertices,
  // which is sufficient & cheaper
  if (this->n_nodes() == rhs_elem->n_nodes())
    {
      // The number of nodes in the element
      const unsigned int nn = this->n_nodes();

      // Create a vector that contains our global
      // node numbers and those of our neighbor.
      // If the sorted, unique vector is the same size
      // as the number of nodes in both elements then
      // they must be connected to the same nodes.
      //
      // The vector will be no larger than 2*n_nodes(),
      // so we might as well reserve the space.
      std::vector<dof_id_type> common_nodes;
      common_nodes.reserve (2*nn);

      // Add the global indices of the nodes
      for (unsigned int n=0; n<nn; n++)
	{
	  common_nodes.push_back (this->node(n));
	  common_nodes.push_back (rhs_elem->node(n));
	}

      // Sort the vector and find out how long
      // the sorted vector is.
      std::sort (common_nodes.begin(), common_nodes.end());

      iterator new_end = std::unique (common_nodes.begin(),
				      common_nodes.end());

      const int new_size = libmesh_cast_int<int>
        (std::distance (common_nodes.begin(), new_end));

      // If this passes the elements are connected
      // to the same global vertex nodes
      if (new_size == static_cast<int>(nn))
	return true;
    }

  // If we get here it is because the elements either
  // do not have the same number of nodes or they are
  // connected to different nodes.  Either way they
  // are not the same element
  return false;
}



bool Elem::is_semilocal(const processor_id_type my_pid) const
{
  std::set<const Elem *> point_neighbors;

  this->find_point_neighbors(point_neighbors);

  std::set<const Elem*>::const_iterator       it  = point_neighbors.begin();
  const std::set<const Elem*>::const_iterator end = point_neighbors.end();

  for (; it != end; ++it)
    {
      const Elem* elem = *it;
      if (elem->processor_id() == my_pid)
        return true;
    }

  return false;
}



bool Elem::contains_vertex_of(const Elem *e) const
{
  // Our vertices are the first numbered nodes
  for (unsigned int n = 0; n != e->n_vertices(); ++n)
    if (this->contains_point(e->point(n)))
      return true;
  return false;
}



bool Elem::contains_edge_of(const Elem *e) const
{
  unsigned int num_contained_edges = 0;

  // Our vertices are the first numbered nodes
  for (unsigned int n = 0; n != e->n_vertices(); ++n)
    {
      if (this->contains_point(e->point(n)))
	{
	  num_contained_edges++;
	  if(num_contained_edges>=2)
	    {
	      return true;
	    }
	}
    }
  return false;
}



void Elem::find_point_neighbors(const Point &p,
                                std::set<const Elem *> &neighbor_set) const
{
  libmesh_assert(this->contains_point(p));

  neighbor_set.clear();
  neighbor_set.insert(this);

  std::set<const Elem *> untested_set, next_untested_set;
  untested_set.insert(this);

  while (!untested_set.empty())
    {
      // Loop over all the elements in the patch that haven't already
      // been tested
      std::set<const Elem*>::const_iterator       it  = untested_set.begin();
      const std::set<const Elem*>::const_iterator end = untested_set.end();

      for (; it != end; ++it)
        {
          const Elem* elem = *it;

          for (unsigned int s=0; s<elem->n_sides(); s++)
            {
              const Elem* current_neighbor = elem->neighbor(s);
              if (current_neighbor &&
                  current_neighbor != remote_elem)    // we have a real neighbor on this side
                {
                  if (current_neighbor->active())                // ... if it is active
                    {
                      if (current_neighbor->contains_point(p))   // ... and touches p
                        {
                          // Make sure we'll test it
                          if (!neighbor_set.count(current_neighbor))
                            next_untested_set.insert (current_neighbor);

                          // And add it
                          neighbor_set.insert (current_neighbor);
                        }
                    }
#ifdef LIBMESH_ENABLE_AMR
                  else                                 // ... the neighbor is *not* active,
                    {                                  // ... so add *all* neighboring
                                                       // active children that touch p
                      std::vector<const Elem*> active_neighbor_children;

                      current_neighbor->active_family_tree_by_neighbor
                        (active_neighbor_children, elem);

                      std::vector<const Elem*>::const_iterator
                        child_it = active_neighbor_children.begin();
                      const std::vector<const Elem*>::const_iterator
                        child_end = active_neighbor_children.end();
                      for (; child_it != child_end; ++child_it)
                        {
                          const Elem *current_child = *child_it;
                          if (current_child->contains_point(p))
                            {
                              // Make sure we'll test it
                              if (!neighbor_set.count(current_child))
                                next_untested_set.insert (current_child);

                              neighbor_set.insert (current_child);
                            }
                        }
                    }
#endif // #ifdef LIBMESH_ENABLE_AMR
                }
            }
        }
      untested_set.swap(next_untested_set);
      next_untested_set.clear();
    }
}



void Elem::find_point_neighbors(std::set<const Elem *> &neighbor_set) const
{
  neighbor_set.clear();
  neighbor_set.insert(this);

  std::set<const Elem *> untested_set, next_untested_set;
  untested_set.insert(this);

  while (!untested_set.empty())
    {
      // Loop over all the elements in the patch that haven't already
      // been tested
      std::set<const Elem*>::const_iterator       it  = untested_set.begin();
      const std::set<const Elem*>::const_iterator end = untested_set.end();

      for (; it != end; ++it)
        {
          const Elem* elem = *it;

          for (unsigned int s=0; s<elem->n_sides(); s++)
            {
              const Elem* current_neighbor = elem->neighbor(s);
              if (current_neighbor &&
                  current_neighbor != remote_elem)    // we have a real neighbor on this side
                {
                  if (current_neighbor->active())                // ... if it is active
                    {
                      if (this->contains_vertex_of(current_neighbor) // ... and touches us
                          || current_neighbor->contains_vertex_of(this))
                        {
                          // Make sure we'll test it
                          if (!neighbor_set.count(current_neighbor))
                            next_untested_set.insert (current_neighbor);

                          // And add it
                          neighbor_set.insert (current_neighbor);
                        }
                    }
#ifdef LIBMESH_ENABLE_AMR
                  else                                 // ... the neighbor is *not* active,
                    {                                  // ... so add *all* neighboring
                                                       // active children
                      std::vector<const Elem*> active_neighbor_children;

                      current_neighbor->active_family_tree_by_neighbor
                        (active_neighbor_children, elem);

                      std::vector<const Elem*>::const_iterator
                        child_it = active_neighbor_children.begin();
                      const std::vector<const Elem*>::const_iterator
                        child_end = active_neighbor_children.end();
                      for (; child_it != child_end; ++child_it)
                        {
                          const Elem *current_child = *child_it;
                          if (this->contains_vertex_of(current_child) ||
                              (current_child)->contains_vertex_of(this))
                            {
                              // Make sure we'll test it
                              if (!neighbor_set.count(current_child))
                                next_untested_set.insert (current_child);

                              neighbor_set.insert (current_child);
                            }
                        }
                    }
#endif // #ifdef LIBMESH_ENABLE_AMR
                }
            }
        }
      untested_set.swap(next_untested_set);
      next_untested_set.clear();
    }
}



void Elem::find_edge_neighbors(const Point& p1,
                               const Point& p2,
			       std::set<const Elem *> &neighbor_set) const
{
  // Simple but perhaps suboptimal code: find elements containing the
  // first point, then winnow this set down by removing elements which
  // don't also contain the second point

  libmesh_assert(this->contains_point(p2));
  this->find_point_neighbors(p1, neighbor_set);

  std::set<const Elem*>::iterator        it = neighbor_set.begin();
  const std::set<const Elem*>::iterator end = neighbor_set.end();

  while(it != end) {
    std::set<const Elem*>::iterator current = it++;

    const Elem* elem = *current;
    // This won't invalidate iterator it, because it is already
    // pointing to the next element
    if (!elem->contains_point(p2))
      neighbor_set.erase(current);
  }
}



void Elem::find_edge_neighbors(std::set<const Elem *> &neighbor_set) const
{
  neighbor_set.clear();
  neighbor_set.insert(this);

  std::set<const Elem *> untested_set, next_untested_set;
  untested_set.insert(this);

  while (!untested_set.empty())
    {
      // Loop over all the elements in the patch that haven't already
      // been tested
      std::set<const Elem*>::const_iterator       it  = untested_set.begin();
      const std::set<const Elem*>::const_iterator end = untested_set.end();

      for (; it != end; ++it)
        {
          const Elem* elem = *it;

          for (unsigned int s=0; s<elem->n_sides(); s++)
            {
              const Elem* current_neighbor = elem->neighbor(s);
              if (current_neighbor &&
                  current_neighbor != remote_elem)    // we have a real neighbor on this side
                {
                  if (current_neighbor->active())                // ... if it is active
                    {
                      if (this->contains_edge_of(current_neighbor) // ... and touches us
                          || current_neighbor->contains_edge_of(this))
                        {
                          // Make sure we'll test it
                          if (!neighbor_set.count(current_neighbor))
                            next_untested_set.insert (current_neighbor);

                          // And add it
                          neighbor_set.insert (current_neighbor);
                        }
                    }
#ifdef LIBMESH_ENABLE_AMR
                  else                                 // ... the neighbor is *not* active,
                    {                                  // ... so add *all* neighboring
                                                       // active children
                      std::vector<const Elem*> active_neighbor_children;

                      current_neighbor->active_family_tree_by_neighbor
                        (active_neighbor_children, elem);

                      std::vector<const Elem*>::const_iterator
                        child_it = active_neighbor_children.begin();
                      const std::vector<const Elem*>::const_iterator
                        child_end = active_neighbor_children.end();
                      for (; child_it != child_end; ++child_it)
                        {
                          const Elem *current_child = *child_it;
                          if (this->contains_edge_of(*child_it) ||
                              (*child_it)->contains_edge_of(this))
                            {
                              // Make sure we'll test it
                              if (!neighbor_set.count(current_child))
                                next_untested_set.insert (current_child);

                              neighbor_set.insert (current_child);
                            }
                        }
                    }
#endif // #ifdef LIBMESH_ENABLE_AMR
                }
            }
        }
      untested_set.swap(next_untested_set);
      next_untested_set.clear();
    }
}

#ifdef LIBMESH_ENABLE_PERIODIC

Elem* Elem::topological_neighbor (const unsigned int i,
                                  MeshBase & mesh,
                                  const PointLocatorBase& point_locator,
                                  const PeriodicBoundaries * pb)
{
  libmesh_assert_less (i, this->n_neighbors());

  Elem * neighbor_i = this->neighbor(i);
  if (neighbor_i != NULL)
    return neighbor_i;

  if (pb)
  {
    // Since the neighbor is NULL it must be on a boundary. We need
    // see if this is a periodic boundary in which case it will have a
    // topological neighbor

    std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(this, i);
    for (std::vector<boundary_id_type>::iterator j = boundary_ids.begin(); j != boundary_ids.end(); ++j)
      if (pb->boundary(*j))
      {
        // Since the point locator inside of periodic boundaries
        // returns a const pointer we will retrieve the proper
        // pointer directly from the mesh object.  Also since coarse
        // elements do not have more refined neighbors we need to make
        // sure that we don't return one of these types of neighbors.
        neighbor_i = mesh.elem(pb->neighbor(*j, point_locator, this, i)->id());
        if (level() < neighbor_i->level())
          neighbor_i = neighbor_i->parent();
        return neighbor_i;
      }
  }

  return NULL;
}



const Elem* Elem::topological_neighbor (const unsigned int i,
                                        const MeshBase & mesh,
                                        const PointLocatorBase& point_locator,
                                        const PeriodicBoundaries * pb) const
{
  libmesh_assert_less (i, this->n_neighbors());

  const Elem * neighbor_i = this->neighbor(i);
  if (neighbor_i != NULL)
    return neighbor_i;

  if (pb)
  {
    // Since the neighbor is NULL it must be on a boundary. We need
    // see if this is a periodic boundary in which case it will have a
    // topological neighbor

    std::vector<boundary_id_type> boundary_ids = mesh.boundary_info->boundary_ids(this, i);
    for (std::vector<boundary_id_type>::iterator j = boundary_ids.begin(); j != boundary_ids.end(); ++j)
      if (pb->boundary(*j))
      {
        // Since the point locator inside of periodic boundaries
        // returns a const pointer we will retrieve the proper
        // pointer directly from the mesh object.  Also since coarse
        // elements do not have more refined neighbors we need to make
        // sure that we don't return one of these types of neighbors.
        neighbor_i = mesh.elem(pb->neighbor(*j, point_locator, this, i)->id());
        if (level() < neighbor_i->level())
          neighbor_i = neighbor_i->parent();
        return neighbor_i;
      }
  }

  return NULL;
}


bool Elem::has_topological_neighbor (const Elem* elem,
                                     const MeshBase & mesh,
                                     const PointLocatorBase& point_locator,
                                     PeriodicBoundaries * pb) const
{
  // First see if this is a normal "interior" neighbor
  if (has_neighbor(elem))
    return true;

  for (unsigned int n=0; n<this->n_neighbors(); n++)
    if (this->topological_neighbor(n, mesh, point_locator, pb))
      return true;

  return false;
}


#endif

#ifdef DEBUG

void Elem::libmesh_assert_valid_node_pointers() const
{
  libmesh_assert(this->valid_id());
  for (unsigned int n=0; n != this->n_nodes(); ++n)
    {
      libmesh_assert(this->get_node(n));
      libmesh_assert(this->get_node(n)->valid_id());
    }
}



void Elem::libmesh_assert_valid_neighbors() const
{
  for (unsigned int s=0; s<this->n_neighbors(); s++)
    {
      const Elem *neigh = this->neighbor(s);

      // Any element might have a remote neighbor; checking
      // to make sure that's not inaccurate is tough.
      if (neigh == remote_elem)
        continue;

      if (neigh)
        {
          // Only subactive elements have subactive neighbors
          libmesh_assert (this->subactive() || !neigh->subactive());

          const Elem *elem = this;

          // If we're subactive but our neighbor isn't, its
          // return neighbor link will be to our first active
          // ancestor OR to our inactive ancestor of the same
          // level as neigh,
          if (this->subactive() && !neigh->subactive())
            {
              for (elem = this; !elem->active();
                   elem = elem->parent())
                libmesh_assert(elem);
            }
          else
            {
              unsigned int rev = neigh->which_neighbor_am_i(elem);
              libmesh_assert_less (rev, neigh->n_neighbors());

              if (this->subactive() && !neigh->subactive())
                {
                  while (neigh->neighbor(rev) != elem)
                    {
                      libmesh_assert(elem->parent());
                      elem = elem->parent();
                    }
                }
              else
                {
                  Elem *nn = neigh->neighbor(rev);
                  libmesh_assert(nn);

                  for (; elem != nn; elem = elem->parent())
                    libmesh_assert(elem);
                }
            }
        }
      // If we don't have a neighbor and we're not subactive, our
      // ancestors shouldn't have any neighbors in this same
      // direction.
      else if (!this->subactive())
        {
          const Elem *my_parent = this->parent();
          if (my_parent &&
          // A parent with a different dimension isn't really one of
          // our ancestors, it means we're on a boundary mesh and this
          // is an interior mesh element for which we're on a side.
          // Nothing to test for in that case.
              (my_parent->dim() == this->dim()))
            libmesh_assert (!my_parent->neighbor(s));
        }
    }
}

#endif // DEBUG



void Elem::make_links_to_me_local(unsigned int n)
{
  Elem *neigh = this->neighbor(n);

  // Don't bother calling this function unless it's necessary
  libmesh_assert(neigh);
  libmesh_assert(!neigh->is_remote());

  // We never have neighbors more refined than us
  libmesh_assert_less_equal (neigh->level(), this->level());

  // We never have subactive neighbors of non subactive elements
  libmesh_assert(!neigh->subactive() || this->subactive());

  // If we have a neighbor less refined than us then it must not
  // have any more refined active descendants we could have
  // pointed to instead.
  libmesh_assert(neigh->level() == this->level() ||
                 neigh->active());

  // If neigh is at our level, then its family might have
  // remote_elem neighbor links which need to point to us
  // instead, but if not, then we're done.
  if (neigh->level() != this->level())
    return;

  // If neigh is subactive then we're not updating its neighbor links
  // FIXME - this needs to change when we start using subactive
  // elements for more than just the two-phase
  // restriction/prolongation projections.
  if (neigh->subactive())
    return;

  // What side of neigh are we on?  We can't use the usual Elem
  // method because we're in the middle of restoring topology
  const AutoPtr<Elem> my_side = this->side(n);
  unsigned int nn = 0;
  for (; nn != neigh->n_sides(); ++nn)
    {
      const AutoPtr<Elem> neigh_side = neigh->side(nn);
      if (*my_side == *neigh_side)
        break;
    }

  // we had better be on *some* side of neigh
  libmesh_assert_less (nn, neigh->n_sides());

  // Find any elements that ought to point to elem
  std::vector<const Elem*> neigh_family;
#ifdef LIBMESH_ENABLE_AMR
  if (this->active())
    neigh->family_tree_by_side(neigh_family, nn);
  else
#endif
    neigh_family.push_back(neigh);

  // And point them to elem
  for (unsigned int i = 0; i != neigh_family.size(); ++i)
    {
      Elem* neigh_family_member = const_cast<Elem*>(neigh_family[i]);

      // Ideally, the neighbor link ought to either be correct
      // already or ought to be to remote_elem.
      //
      // However, if we're redistributing a newly created elem,
      // after an AMR step but before find_neighbors has fixed up
      // neighbor links, we might have an out of date neighbor
      // link to elem's parent instead.
#ifdef LIBMESH_ENABLE_AMR
      libmesh_assert((neigh_family_member->neighbor(nn) == this) ||
                     (neigh_family_member->neighbor(nn) == remote_elem)
		     || ((this->refinement_flag() == JUST_REFINED) &&
                      (this->parent() != NULL) &&
                      (neigh_family_member->neighbor(nn) == this->parent())));
#else
      libmesh_assert((neigh_family_member->neighbor(nn) == this) ||
                     (neigh_family_member->neighbor(nn) == remote_elem));
#endif

      neigh_family_member->set_neighbor(nn, this);
    }
}


void Elem::make_links_to_me_remote()
{
  libmesh_assert_not_equal_to (this, remote_elem);

  // We need to have handled any children first
#if defined(LIBMESH_ENABLE_AMR) && defined(DEBUG)
  if (this->has_children())
    for (unsigned int c = 0; c != this->n_children(); ++c)
      {
        Elem *current_child = this->child(c);
        libmesh_assert_equal_to (current_child, remote_elem);
      }
#endif

  // Remotify any neighbor links to non-subactive elements
  if (!this->subactive())
    {
      for (unsigned int s = 0; s != this->n_sides(); ++s)
        {
          Elem *neigh = this->neighbor(s);
          if (neigh && neigh != remote_elem && !neigh->subactive())
            {
	      // My neighbor should never be more refined than me; my real
	      // neighbor would have been its parent in that case.
	      libmesh_assert_greater_equal (this->level(), neigh->level());

              if (this->level() == neigh->level() &&
                  neigh->has_neighbor(this))
                {
#ifdef LIBMESH_ENABLE_AMR
	          // My neighbor may have descendants which also consider me a
	          // neighbor
                  std::vector<const Elem*> family;
                  neigh->family_tree_by_neighbor (family, this);

                  // FIXME - There's a lot of ugly const_casts here; we
                  // may want to make remote_elem non-const and create
                  // non-const versions of the family_tree methods
                  for (unsigned int i=0; i != family.size(); ++i)
                    {
                      Elem *n = const_cast<Elem*>(family[i]);
                      libmesh_assert (n);
                      if (n == remote_elem)
                        continue;
                      unsigned int my_s = n->which_neighbor_am_i(this);
                      libmesh_assert_less (my_s, n->n_neighbors());
                      libmesh_assert_equal_to (n->neighbor(my_s), this);
                      n->set_neighbor(my_s, const_cast<RemoteElem*>(remote_elem));
                    }
#else
                  unsigned int my_s = neigh->which_neighbor_am_i(this);
                  libmesh_assert_less (my_s, neigh->n_neighbors());
                  libmesh_assert_equal_to (neigh->neighbor(my_s), this);
                  neigh->set_neighbor(my_s, const_cast<RemoteElem*>(remote_elem));
#endif
                }
#ifdef LIBMESH_ENABLE_AMR
              // Even if my neighbor doesn't link back to me, it might
	      // have subactive descendants which do
	      else if (neigh->has_children())
                {
                  // If my neighbor at the same level doesn't have me as a
	          // neighbor, I must be subactive
	          libmesh_assert(this->level() > neigh->level() ||
			         this->subactive());

                  // My neighbor must have some ancestor of mine as a
	          // neighbor
	          Elem *my_ancestor = this->parent();
	          libmesh_assert(my_ancestor);
                  while (!neigh->has_neighbor(my_ancestor))
                    {
                      my_ancestor = my_ancestor->parent();
	              libmesh_assert(my_ancestor);
                    }

	          // My neighbor may have descendants which consider me a
	          // neighbor
                  std::vector<const Elem*> family;
                  neigh->family_tree_by_subneighbor (family, my_ancestor, this);

                  // FIXME - There's a lot of ugly const_casts here; we
                  // may want to make remote_elem non-const and create
                  // non-const versions of the family_tree methods
                  for (unsigned int i=0; i != family.size(); ++i)
                    {
                      Elem *n = const_cast<Elem*>(family[i]);
                      libmesh_assert (n);
                      if (n == remote_elem)
                        continue;
                      unsigned int my_s = n->which_neighbor_am_i(this);
                      libmesh_assert_less (my_s, n->n_neighbors());
                      libmesh_assert_equal_to (n->neighbor(my_s), this);
                      n->set_neighbor(my_s, const_cast<RemoteElem*>(remote_elem));
                    }
                }
#endif
            }
        }
    }

#ifdef LIBMESH_ENABLE_AMR
  // Remotify parent's child link
  Elem *my_parent = this->parent();
  if (my_parent &&
      // As long as it's not already remote
      my_parent != remote_elem &&
      // And it's a real parent, not an interior parent
      this->dim() == my_parent->dim())
    {
      unsigned int me = my_parent->which_child_am_i(this);
      libmesh_assert_equal_to (my_parent->child(me), this);
      my_parent->set_child(me, const_cast<RemoteElem*>(remote_elem));
    }
#endif
}



void Elem::write_connectivity (std::ostream& out_stream,
			       const IOPackage iop) const
{
  libmesh_assert (out_stream.good());
  libmesh_assert(_nodes);
  libmesh_assert_not_equal_to (iop, INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
	// This connectivity vector will be used repeatedly instead
	// of being reconstructed inside the loop.
	std::vector<dof_id_type> conn;
	for (unsigned int sc=0; sc <this->n_sub_elem(); sc++)
	  {
	    this->connectivity(sc, TECPLOT, conn);

	    std::copy(conn.begin(),
		      conn.end(),
		      std::ostream_iterator<dof_id_type>(out_stream, " "));

	    out_stream << '\n';
	  }
	return;
      }

    case UCD:
      {
	for (unsigned int i=0; i<this->n_nodes(); i++)
	  out_stream << this->node(i)+1 << "\t";

	out_stream << '\n';
	return;
      }

    default:
      libmesh_error();
    }

  libmesh_error();
}


// void Elem::write_tecplot_connectivity(std::ostream& out) const
// {
//   libmesh_assert (!out.bad());
//   libmesh_assert(_nodes);

//   // This connectivity vector will be used repeatedly instead
//   // of being reconstructed inside the loop.
//   std::vector<dof_id_type> conn;
//   for (unsigned int sc=0; sc <this->n_sub_elem(); sc++)
//     {
//       this->connectivity(sc, TECPLOT, conn);

//       std::copy(conn.begin(),
//  		conn.end(),
//  		std::ostream_iterator<dof_id_type>(out, " "));

//       out << std::endl;
//     }
// }



// void Elem::write_ucd_connectivity(std::ostream &out) const
// {
//   libmesh_assert (out);
//   libmesh_assert(_nodes);

//   for (unsigned int i=0; i<this->n_nodes(); i++)
//     out << this->node(i)+1 << "\t";

//   out << std::endl;
// }



Real Elem::quality (const ElemQuality q) const
{
  switch (q)
    {
      /**
       * I don't know what to do for this metric.
       */
    default:
      {
	libmesh_here();

	libMesh::err << "ERROR:  unknown quality metric: "
		      << q
		      << std::endl
		      << "Cowardly returning 1."
		      << std::endl;

	return 1.;
      }
    }


    // Will never get here...
    libmesh_error();
    return 0.;
}



bool Elem::ancestor() const
{
#ifdef LIBMESH_ENABLE_AMR

  if (this->active())
    return false;

if (!this->has_children())
    return false;
  if (this->child(0)->active())
    return true;

  return this->child(0)->ancestor();
#else
  return false;
#endif
}



#ifdef LIBMESH_ENABLE_AMR

void Elem::add_child (Elem* elem)
{
  if(_children == NULL)
    {
      _children = new Elem*[this->n_children()];

      for (unsigned int c=0; c<this->n_children(); c++)
	this->set_child(c, NULL);
    }

  for (unsigned int c=0; c<this->n_children(); c++)
    {
      if(this->_children[c] == NULL || this->_children[c] == remote_elem)
	{
	  libmesh_assert_equal_to (this, elem->parent());
	  this->set_child(c, elem);
	  return;
	}
    }

  libMesh::err << "Error: Tried to add a child to an element with full children array"
                << std::endl;
  libmesh_error();
}



void Elem::add_child (Elem* elem, unsigned int c)
{
  if(!this->has_children())
    {
      _children = new Elem*[this->n_children()];

      for (unsigned int i=0; i<this->n_children(); i++)
	this->set_child(i, NULL);
    }

  libmesh_assert (this->_children[c] == NULL || this->child(c) == remote_elem);
  libmesh_assert (elem == remote_elem || this == elem->parent());

  this->set_child(c, elem);
}



void Elem::replace_child (Elem* elem, unsigned int c)
{
  libmesh_assert(this->has_children());

  libmesh_assert(this->child(c));

  this->set_child(c, elem);
}



bool Elem::is_child_on_edge(const unsigned int libmesh_dbg_var(c),
                            const unsigned int e) const
{
  libmesh_assert_less (c, this->n_children());
  libmesh_assert_less (e, this->n_edges());

  AutoPtr<Elem> my_edge = this->build_edge(e);
  AutoPtr<Elem> child_edge = this->build_edge(e);

  // We're assuming that an overlapping child edge has the same
  // number and orientation as its parent
  return (child_edge->node(0) == my_edge->node(0) ||
      child_edge->node(1) == my_edge->node(1));
}


void Elem::family_tree (std::vector<const Elem*>& family,
			const bool reset) const
{
  // The "family tree" doesn't include subactive elements
  libmesh_assert(!this->subactive());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // Add this element to the family tree.
  family.push_back(this);

  // Recurse into the elements children, if it has them.
  // Do not clear the vector any more.
  if (!this->active())
    for (unsigned int c=0; c<this->n_children(); c++)
      if (!this->child(c)->is_remote())
	this->child(c)->family_tree (family, false);
}



void Elem::total_family_tree (std::vector<const Elem*>& family,
			      const bool reset) const
{
  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // Add this element to the family tree.
  family.push_back(this);

  // Recurse into the elements children, if it has them.
  // Do not clear the vector any more.
  if (this->has_children())
    for (unsigned int c=0; c<this->n_children(); c++)
      if (!this->child(c)->is_remote())
	this->child(c)->total_family_tree (family, false);
}



void Elem::active_family_tree (std::vector<const Elem*>& active_family,
			       const bool reset) const
{
  // The "family tree" doesn't include subactive elements
  libmesh_assert(!this->subactive());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    active_family.clear();

  // Add this element to the family tree if it is active
  if (this->active())
    active_family.push_back(this);

  // Otherwise recurse into the element's children.
  // Do not clear the vector any more.
  else
    for (unsigned int c=0; c<this->n_children(); c++)
      if (!this->child(c)->is_remote())
	this->child(c)->active_family_tree (active_family, false);
}



void Elem::family_tree_by_side (std::vector<const Elem*>& family,
                                const unsigned int s,
                                const bool reset)  const
{
  // The "family tree" doesn't include subactive elements
  libmesh_assert(!this->subactive());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  libmesh_assert_less (s, this->n_sides());

  // Add this element to the family tree.
  family.push_back(this);

  // Recurse into the elements children, if it has them.
  // Do not clear the vector any more.
  if (!this->active())
    for (unsigned int c=0; c<this->n_children(); c++)
      if (!this->child(c)->is_remote() && this->is_child_on_side(c, s))
        this->child(c)->family_tree_by_side (family, s, false);
}



void Elem::active_family_tree_by_side (std::vector<const Elem*>& family,
                                       const unsigned int s,
                                       const bool reset) const
{
  // The "family tree" doesn't include subactive elements
  libmesh_assert(!this->subactive());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  libmesh_assert_less (s, this->n_sides());

  // Add an active element to the family tree.
  if (this->active())
    family.push_back(this);

  // Or recurse into an ancestor element's children.
  // Do not clear the vector any more.
  else
    for (unsigned int c=0; c<this->n_children(); c++)
      if (!this->child(c)->is_remote() && this->is_child_on_side(c, s))
        this->child(c)->active_family_tree_by_side (family, s, false);
}



void Elem::family_tree_by_neighbor (std::vector<const Elem*>& family,
                                    const Elem* neighbor_in,
			            const bool reset) const
{
  // The "family tree" doesn't include subactive elements
  libmesh_assert(!this->subactive());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // This only makes sense if we're already a neighbor
  libmesh_assert (this->has_neighbor(neighbor_in));

  // Add this element to the family tree.
  family.push_back(this);

  // Recurse into the elements children, if it's not active.
  // Do not clear the vector any more.
  if (!this->active())
    for (unsigned int c=0; c<this->n_children(); c++)
      {
        Elem *current_child = this->child(c);
        if (current_child != remote_elem && current_child->has_neighbor(neighbor_in))
          current_child->family_tree_by_neighbor (family, neighbor_in, false);
      }
}



void Elem::family_tree_by_subneighbor (std::vector<const Elem*>& family,
                                       const Elem* neighbor_in,
                                       const Elem* subneighbor,
			               const bool reset) const
{
  // The "family tree" doesn't include subactive elements
  libmesh_assert(!this->subactive());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // To simplifly this function we need an existing neighbor
  libmesh_assert (neighbor_in);
  libmesh_assert_not_equal_to (neighbor_in, remote_elem);
  libmesh_assert (this->has_neighbor(neighbor_in));

  // This only makes sense if subneighbor descends from neighbor
  libmesh_assert (subneighbor);
  libmesh_assert_not_equal_to (subneighbor, remote_elem);
  libmesh_assert (neighbor_in->is_ancestor_of(subneighbor));

  // Add this element to the family tree if applicable.
  if (neighbor_in == subneighbor)
    family.push_back(this);

  // Recurse into the elements children, if it's not active.
  // Do not clear the vector any more.
  if (!this->active())
    for (unsigned int c=0; c != this->n_children(); ++c)
      {
        Elem *current_child = this->child(c);
        if (current_child != remote_elem)
          for (unsigned int s=0; s != current_child->n_sides(); ++s)
            {
	      Elem *child_neigh = current_child->neighbor(s);
              if (child_neigh &&
		  (child_neigh == neighbor_in ||
		   (child_neigh->parent() == neighbor_in &&
		    child_neigh->is_ancestor_of(subneighbor))))
                current_child->family_tree_by_subneighbor (family, child_neigh,
                                                           subneighbor, false);
            }
      }
}



void Elem::active_family_tree_by_neighbor (std::vector<const Elem*>& family,
                                           const Elem* neighbor_in,
			                   const bool reset) const
{
  // The "family tree" doesn't include subactive elements
  libmesh_assert(!this->subactive());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // This only makes sense if we're already a neighbor
  if (this->level() >= neighbor_in->level())
    libmesh_assert (this->has_neighbor(neighbor_in));

  // Add an active element to the family tree.
  if (this->active())
    family.push_back(this);

  // Or recurse into an ancestor element's children.
  // Do not clear the vector any more.
  else if (!this->active())
    for (unsigned int c=0; c<this->n_children(); c++)
      {
        Elem *current_child = this->child(c);
        if (current_child != remote_elem && current_child->has_neighbor(neighbor_in))
          current_child->active_family_tree_by_neighbor (family, neighbor_in, false);
      }
}



unsigned int Elem::min_p_level_by_neighbor(const Elem* neighbor_in,
                                           unsigned int current_min) const
{
  libmesh_assert(!this->subactive());
  libmesh_assert(neighbor_in->active());

  // If we're an active element this is simple
  if (this->active())
    return std::min(current_min, this->p_level());

  libmesh_assert(has_neighbor(neighbor_in));

  // The p_level() of an ancestor element is already the minimum
  // p_level() of its children - so if that's high enough, we don't
  // need to examine any children.
  if (current_min <= this->p_level())
    return current_min;

  unsigned int min_p_level = current_min;

  for (unsigned int c=0; c<this->n_children(); c++)
    {
      const Elem* const current_child = this->child(c);
      if (current_child != remote_elem && current_child->has_neighbor(neighbor_in))
        min_p_level =
	  current_child->min_p_level_by_neighbor(neighbor_in,
                                                 min_p_level);
    }

  return min_p_level;
}


unsigned int Elem::min_new_p_level_by_neighbor(const Elem* neighbor_in,
                                               unsigned int current_min) const
{
  libmesh_assert(!this->subactive());
  libmesh_assert(neighbor_in->active());

  // If we're an active element this is simple
  if (this->active())
    {
      unsigned int new_p_level = this->p_level();
      if (this->p_refinement_flag() == Elem::REFINE)
        new_p_level += 1;
      if (this->p_refinement_flag() == Elem::COARSEN)
        {
          libmesh_assert_greater (new_p_level, 0);
          new_p_level -= 1;
        }
      return std::min(current_min, new_p_level);
    }

  libmesh_assert(has_neighbor(neighbor_in));

  unsigned int min_p_level = current_min;

  for (unsigned int c=0; c<this->n_children(); c++)
    {
      const Elem* const current_child = this->child(c);
      if (current_child && current_child != remote_elem)
        if (current_child->has_neighbor(neighbor_in))
          min_p_level =
	    current_child->min_new_p_level_by_neighbor(neighbor_in,
                                                       min_p_level);
    }

  return min_p_level;
}

#endif // #ifdef LIBMESH_ENABLE_AMR




bool Elem::contains_point (const Point& p, Real tol) const
{
  // We currently allow the user to enlarge the bounding box by
  // providing a tol > TOLERANCE (so this routine is identical to
  // Elem::close_to_point()), but print a warning so that the
  // user can eventually switch his code over to calling close_to_point()
  // instead, which is intended to be used for this purpose.
  if ( tol > TOLERANCE )
    {
      libmesh_do_once(libMesh::err
		      << "WARNING: Resizing bounding box to match user-specified tolerance!\n"
		      << "In the future, calls to Elem::contains_point() with tol > TOLERANCE\n"
		      << "will be more optimized, but should not be used\n"
		      << "to search for points 'close to' elements!\n"
		      << "Instead, use Elem::close_to_point() for this purpose.\n"
		      << std::endl;);
      return this->point_test(p, tol, tol);
    }
  else
    return this->point_test(p, TOLERANCE, tol);
}




bool Elem::close_to_point (const Point& p, Real tol) const
{
  // This test uses the user's passed-in tolerance for the
  // bounding box test as well, thereby allowing the routine to
  // find points which are not only "in" the element, but also
  // "nearby" to within some tolerance.
  return this->point_test(p, tol, tol);
}




bool Elem::point_test(const Point& p, Real box_tol, Real map_tol) const
{
  libmesh_assert_greater (box_tol, 0.);
  libmesh_assert_greater (map_tol, 0.);

  // This is a great optimization on first order elements, but it
  // could return false negatives on higher orders
  if (this->default_order() == FIRST)
    {
      // Check to make sure the element *could* contain this point, so we
      // can avoid an expensive inverse_map call if it doesn't.
      bool
#if LIBMESH_DIM > 2
           point_above_min_z = false,
           point_below_max_z = false,
#endif
#if LIBMESH_DIM > 1
           point_above_min_y = false,
           point_below_max_y = false,
#endif
           point_above_min_x = false,
           point_below_max_x = false;

      // For relative bounding box checks in physical space
      const Real my_hmax = this->hmax();

      for (unsigned int n=0; n != this->n_nodes(); ++n)
        {
          Point pe = this->point(n);
          point_above_min_x = point_above_min_x || (pe(0) - my_hmax*box_tol <= p(0));
          point_below_max_x = point_below_max_x || (pe(0) + my_hmax*box_tol >= p(0));
#if LIBMESH_DIM > 1
          point_above_min_y = point_above_min_y || (pe(1) - my_hmax*box_tol <= p(1));
          point_below_max_y = point_below_max_y || (pe(1) + my_hmax*box_tol >= p(1));
#endif
#if LIBMESH_DIM > 2
          point_above_min_z = point_above_min_z || (pe(2) - my_hmax*box_tol <= p(2));
          point_below_max_z = point_below_max_z || (pe(2) + my_hmax*box_tol >= p(2));
#endif
        }

      if (
#if LIBMESH_DIM > 2
          !point_above_min_z ||
          !point_below_max_z ||
#endif
#if LIBMESH_DIM > 1
          !point_above_min_y ||
          !point_below_max_y ||
#endif
          !point_above_min_x ||
          !point_below_max_x)
        return false;
    }

  // Declare a basic FEType.  Will be a Lagrange
  // element by default.
  FEType fe_type(this->default_order());

  // To be on the safe side, we converge the inverse_map() iteration
  // to a slightly tighter tolerance than that requested by the
  // user...
  const Point mapped_point = FEInterface::inverse_map(this->dim(),
						      fe_type,
						      this,
						      p,
						      0.1*map_tol, // <- this is |dx| tolerance, the Newton residual should be ~ |dx|^2
						      /*secure=*/ false);

  // Check that the refspace point maps back to p!  This is only necessary
  // for 1D and 2D elements, 3D elements always live in 3D.
  //
  // TODO: The contains_point() function could most likely be implemented
  // more efficiently in the element sub-classes themselves, at least for
  // the linear element types.
  if (this->dim() < 3)
    {
      Point xyz = FEInterface::map(this->dim(),
                                   fe_type,
                                   this,
                                   mapped_point);

      // Compute the distance between the original point and the re-mapped point.
      // They should be in the same place.
      Real dist = (xyz - p).size();


      // If dist is larger than some fraction of the tolerance, then return false.
      // This can happen when e.g. a 2D element is living in 3D, and
      // FEInterface::inverse_map() maps p onto the projection of the element,
      // effectively "tricking" FEInterface::on_reference_element().
      if (dist > this->hmax() * map_tol)
        return false;
    }



  return FEInterface::on_reference_element(mapped_point, this->type(), map_tol);
}




void Elem::print_info (std::ostream& os) const
{
  os << this->get_info()
     << std::endl;
}



std::string Elem::get_info () const
{
  std::ostringstream oss;

  oss << "  Elem Information"                                      << '\n'
      << "   id()=";

  if (this->valid_id())
    oss << this->id();
  else
    oss << "invalid";

  oss << ", processor_id()=" << this->processor_id()               << '\n';

  oss << "   type()="    << Utility::enum_to_string(this->type())  << '\n'
      << "   dim()="     << this->dim()                            << '\n'
      << "   n_nodes()=" << this->n_nodes()                        << '\n';

  for (unsigned int n=0; n != this->n_nodes(); ++n)
    oss << "    " << n << *this->get_node(n);

  oss << "   n_sides()=" << this->n_sides()                        << '\n';

  for (unsigned int s=0; s != this->n_sides(); ++s)
    {
      oss << "    neighbor(" << s << ")=";
      if (this->neighbor(s))
        oss << this->neighbor(s)->id() << '\n';
      else
        oss << "NULL\n";
    }

  oss << "   hmin()=" << this->hmin()
      << ", hmax()=" << this->hmax()                               << '\n'
      << "   volume()=" << this->volume()                          << '\n'
      << "   active()=" << this->active()
      << ", ancestor()=" << this->ancestor()
      << ", subactive()=" << this->subactive()
      << ", has_children()=" << this->has_children()               << '\n'
      << "   parent()=";
  if (this->parent())
    oss << this->parent()->id() << '\n';
  else
    oss << "NULL\n";
  oss << "   level()=" << this->level()
      << ", p_level()=" << this->p_level()                         << '\n'
#ifdef LIBMESH_ENABLE_AMR
      << "   refinement_flag()=" << Utility::enum_to_string(this->refinement_flag())        << '\n'
      << "   p_refinement_flag()=" << Utility::enum_to_string(this->p_refinement_flag())    << '\n'
#endif
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
      << "   infinite()=" << this->infinite()    << '\n';
      if (this->infinite())
        oss << "   origin()=" << this->origin()    << '\n'
#endif
      ;

  oss << "   DoFs=";
  for (unsigned int s=0; s != this->n_systems(); ++s)
    for (unsigned int v=0; v != this->n_vars(s); ++v)
      for (unsigned int c=0; c != this->n_comp(s,v); ++c)
        oss << '(' << s << '/' << v << '/' << this->dof_number(s,v,c) << ") ";


  return oss.str();
}



void Elem::nullify_neighbors ()
{
  // Tell any of my neighbors about my death...
  // Looks strange, huh?
  for (unsigned int n=0; n<this->n_neighbors(); n++)
    {
      Elem* current_neighbor = this->neighbor(n);
      if (current_neighbor && current_neighbor != remote_elem)
        {
	  // Note:  it is possible that I see the neighbor
	  // (which is coarser than me)
	  // but they don't see me, so avoid that case.
	  if (current_neighbor->level() == this->level())
	    {
	      const unsigned int w_n_a_i = current_neighbor->which_neighbor_am_i(this);
              libmesh_assert_less (w_n_a_i, current_neighbor->n_neighbors());
	      current_neighbor->set_neighbor(w_n_a_i, NULL);
	      this->set_neighbor(n, NULL);
	    }
        }
    }
}




unsigned int Elem::n_second_order_adjacent_vertices (const unsigned int) const
{
  // for linear elements, always return 0
  return 0;
}



unsigned short int Elem::second_order_adjacent_vertex (const unsigned int,
						       const unsigned int) const
{
  // for linear elements, always return 0
  return 0;
}



std::pair<unsigned short int, unsigned short int>
Elem::second_order_child_vertex (const unsigned int) const
{
  // for linear elements, always return 0
  return std::pair<unsigned short int, unsigned short int>(0,0);
}



ElemType Elem::first_order_equivalent_type (const ElemType et)
{
  switch (et)
    {
    case EDGE2:
    case EDGE3:
    case EDGE4:
      return EDGE2;
    case TRI3:
    case TRI6:
      return TRI3;
    case QUAD4:
    case QUAD8:
    case QUAD9:
      return QUAD4;
    case TET4:
    case TET10:
      return TET4;
    case HEX8:
    case HEX27:
    case HEX20:
      return HEX8;
    case PRISM6:
    case PRISM15:
    case PRISM18:
      return PRISM6;

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

    case INFQUAD4:
    case INFQUAD6:
      return INFQUAD4;
    case INFHEX8:
    case INFHEX16:
    case INFHEX18:
      return INFHEX8;
    case INFPRISM6:
    case INFPRISM12:
      return INFPRISM6;

#endif

    default:
      // unknown element
      return INVALID_ELEM;
    }
}



ElemType Elem::second_order_equivalent_type (const ElemType et,
					     const bool full_ordered)
{
  /* for second-order elements, always return \p INVALID_ELEM
   * since second-order elements should not be converted
   * into something else.  Only linear elements should
   * return something sensible here
   */
  switch (et)
    {
    case EDGE2:
      {
	// full_ordered not relevant
	return EDGE3;
      }

    case TRI3:
      {
	// full_ordered not relevant
	return TRI6;
      }

    case QUAD4:
      {
	if (full_ordered)
	  return QUAD9;
	else
	  return QUAD8;
      }

    case TET4:
      {
	// full_ordered not relevant
	return TET10;
      }

    case HEX8:
      {
	// see below how this correlates with INFHEX8
	if (full_ordered)
	  return HEX27;
	else
	  return HEX20;
      }

    case PRISM6:
      {
	if (full_ordered)
	  return PRISM18;
	else
	  return PRISM15;
      }

    case PYRAMID5:
      {
	// libmesh_error();
	return INVALID_ELEM;
      }



#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS

    // infinite elements
    case INFEDGE2:
      {
	// libmesh_error();
	return INVALID_ELEM;
      }

    case INFQUAD4:
      {
	// full_ordered not relevant
	return INFQUAD6;
      }

    case INFHEX8:
      {
	/*
	 * Note that this matches with \p Hex8:
	 * For full-ordered, \p InfHex18 and \p Hex27
	 * belong together, and for not full-ordered,
	 * \p InfHex16 and \p Hex20 belong together.
	 */
	if (full_ordered)
	  return INFHEX18;
	else
	  return INFHEX16;
      }

    case INFPRISM6:
      {
	// full_ordered not relevant
	return INFPRISM12;
      }

#endif


    default:
      {
	// second-order element
	return INVALID_ELEM;
      }
    }
}



Elem::side_iterator Elem::boundary_sides_begin()
{
  Predicates::BoundarySide<SideIter> bsp;
  return side_iterator(this->_first_side(), this->_last_side(), bsp);
}




Elem::side_iterator Elem::boundary_sides_end()
{
  Predicates::BoundarySide<SideIter> bsp;
  return side_iterator(this->_last_side(), this->_last_side(), bsp);
}




Real Elem::volume () const
{
  // The default implementation builds a finite element of the correct
  // order and sums up the JxW contributions.  This can be expensive,
  // so the various element types can overload this method and compute
  // the volume more efficiently.
  FEType fe_type (this->default_order() , LAGRANGE);

  AutoPtr<FEBase> fe (FEBase::build(this->dim(),
				    fe_type));

   const std::vector<Real>& JxW = fe->get_JxW();

  // The default quadrature rule should integrate the mass matrix,
  // thus it should be plenty to compute the area
  QGauss qrule (this->dim(), fe_type.default_quadrature_order());

  fe->attach_quadrature_rule(&qrule);

  fe->reinit(this);

  Real vol=0.;
  for (unsigned int qp=0; qp<qrule.n_points(); ++qp)
    vol += JxW[qp];

  return vol;

}


// ------------------------------------------------------------
// Elem::PackedElem static data
const unsigned int Elem::PackedElem::header_size = 10;


// Elem::PackedElem member funcions
void Elem::PackedElem::pack (std::vector<largest_id_type> &conn, const Elem* elem)
{
  libmesh_assert(elem);

  // we can do at least this good. note that hopefully in general
  // the user will already have reserved the full space, which will render
  // this redundant
  conn.reserve (conn.size() + elem->packed_size());

#ifdef LIBMESH_ENABLE_AMR
  conn.push_back (static_cast<largest_id_type>(elem->level()));
  conn.push_back (static_cast<largest_id_type>(elem->p_level()));
  conn.push_back (static_cast<largest_id_type>(elem->refinement_flag()));
  conn.push_back (static_cast<largest_id_type>(elem->p_refinement_flag()));
#else
  conn.push_back (0);
  conn.push_back (0);
  conn.push_back (0);
  conn.push_back (0);
#endif
  conn.push_back (static_cast<largest_id_type>(elem->type()));
  conn.push_back (elem->processor_id());
  conn.push_back (elem->subdomain_id());
  conn.push_back (elem->id());

#ifdef LIBMESH_ENABLE_AMR
  // use parent_ID of -1 to indicate a level 0 element
  if (elem->level() == 0)
    {
      conn.push_back(-1);
      conn.push_back(-1);
    }
  else
    {
      conn.push_back(elem->parent()->id());
      conn.push_back(elem->parent()->which_child_am_i(elem));
    }
#else
  conn.push_back (-1);
  conn.push_back (-1);
#endif

  for (unsigned int n=0; n<elem->n_nodes(); n++)
    conn.push_back (elem->node(n));

  for (unsigned int n=0; n<elem->n_neighbors(); n++)
    {
      Elem *neigh = elem->neighbor(n);
      if (neigh)
        conn.push_back (neigh->id());
      else
        conn.push_back (-1);
    }

  elem->pack_indexing(std::back_inserter(conn));
}



Elem * Elem::PackedElem::unpack (MeshBase &mesh, Elem *parent) const
{

  Elem *elem = Elem::build(this->type(),parent).release();
  libmesh_assert (elem);

#ifdef LIBMESH_ENABLE_AMR
  if (this->level() != 0)
    {
      libmesh_assert(parent);
      parent->add_child(elem, this->which_child_am_i());
      libmesh_assert_equal_to (parent->type(), elem->type());
      libmesh_assert_equal_to (parent->child(this->which_child_am_i()), elem);
    }
#endif

  // Assign the refinement flags and levels
#ifdef LIBMESH_ENABLE_AMR
  elem->set_p_level(this->p_level());
  elem->set_refinement_flag(this->refinement_flag());
  elem->set_p_refinement_flag(this->p_refinement_flag());
  libmesh_assert_equal_to (elem->level(), this->level());

  // If this element definitely should have children, assign
  // remote_elem for now; later unpacked elements may overwrite that.
  if (!elem->active())
    for (unsigned int c=0; c != elem->n_children(); ++c)
      elem->add_child(const_cast<RemoteElem*>(remote_elem), c);
#endif

  // Assign the IDs
  elem->subdomain_id() = this->subdomain_id();
  elem->processor_id() = this->processor_id();
  elem->set_id()       = this->id();

  // Assign the connectivity
  libmesh_assert_equal_to (elem->n_nodes(), this->n_nodes());

  for (unsigned int n=0; n<elem->n_nodes(); n++)
    elem->set_node(n) = mesh.node_ptr (this->node(n));

  // Assign the connectivity
  libmesh_assert_equal_to (elem->n_neighbors(), this->n_neighbors());

  for (unsigned int n=0; n<elem->n_neighbors(); n++)
    {
      dof_id_type neighbor_id = this->neighbor(n);

      // We should only be unpacking elements sent by their owners,
      // and their owners should know all their neighbors
      libmesh_assert_not_equal_to (neighbor_id, remote_elem->id());

      if (neighbor_id == DofObject::invalid_id)
	continue;

      Elem *neigh = mesh.query_elem(neighbor_id);
      if (!neigh)
        {
          elem->set_neighbor(n, const_cast<RemoteElem*>(remote_elem));
	  continue;
	}

      // We never have neighbors more refined than us
      libmesh_assert_less_equal (neigh->level(), elem->level());

      // We never have subactive neighbors of non subactive elements
      libmesh_assert(!neigh->subactive() || elem->subactive());

      // If we have a neighbor less refined than us then it must not
      // have any more refined active descendants we could have
      // pointed to instead.
      libmesh_assert(neigh->level() == elem->level() ||
                     neigh->active());

      elem->set_neighbor(n, neigh);

      // If neigh is at elem's level, then its family might have
      // remote_elem neighbor links which need to point to elem
      // instead, but if not, then we're done.
      if (neigh->level() != elem->level())
	continue;

      // What side of neigh is elem on?  We can't use the usual Elem
      // method because we haven't finished restoring topology
      const AutoPtr<Elem> my_side = elem->side(n);
      unsigned int nn = 0;
      for (; nn != neigh->n_sides(); ++nn)
	{
          const AutoPtr<Elem> neigh_side = neigh->side(nn);
          if (*my_side == *neigh_side)
	    break;
	}

      // elem had better be on *some* side of neigh
      libmesh_assert_less (nn, neigh->n_sides());

      // Find any elements that ought to point to elem
      std::vector<const Elem*> neigh_family;
#ifdef LIBMESH_ENABLE_AMR
      if (!neigh->subactive())
        neigh->family_tree_by_side(neigh_family, nn);
#else
        neigh_family.push_back(neigh);
#endif

      // And point them to elem
      for (unsigned int i = 0; i != neigh_family.size(); ++i)
        {
          Elem* neigh_family_member = const_cast<Elem*>(neigh_family[i]);

	  // The neighbor link ought to either be correct already or
	  // ought to be to remote_elem
          libmesh_assert(neigh_family_member->neighbor(nn) == elem ||
                         neigh_family_member->neighbor(nn) == remote_elem);

	  neigh_family_member->set_neighbor(nn, elem);
        }
    }

  elem->unpack_indexing(this->indices());

  return elem;
}



unsigned int Elem::opposite_side(const unsigned int /*s*/) const
{
  // If the subclass didn't rederive this, using it is an error
  libmesh_not_implemented();
}



unsigned int Elem::opposite_node(const unsigned int /*n*/,
                                 const unsigned int /*s*/) const
{
  // If the subclass didn't rederive this, using it is an error
  libmesh_not_implemented();
}

} // namespace libMesh
