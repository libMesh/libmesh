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



// C++ includes
#include <algorithm> // for std::sort
#include <iterator>  // for std::ostream_iterator

// Local includes
#include "elem.h"
#include "fe_type.h"
#include "fe_interface.h"
#include "edge_edge2.h"
#include "edge_edge3.h"
#include "edge_edge4.h"
#include "edge_inf_edge2.h"
#include "face_tri3.h"
#include "face_tri6.h"
#include "face_quad4.h"
#include "face_quad8.h"
#include "face_quad9.h"
#include "face_inf_quad4.h"
#include "face_inf_quad6.h"
#include "cell_tet4.h"
#include "cell_tet10.h"
#include "cell_hex8.h"
#include "cell_hex20.h"
#include "cell_hex27.h"
#include "cell_inf_hex8.h"
#include "cell_inf_hex16.h"
#include "cell_inf_hex18.h"
#include "cell_prism6.h"
#include "cell_prism15.h"
#include "cell_prism18.h"
#include "cell_inf_prism6.h"
#include "cell_inf_prism12.h"
#include "cell_pyramid5.h"
#include "fe_base.h"
#include "quadrature_gauss.h"
#include "remote_elem.h"
#include "mesh_base.h"

// Initialize static member variables
const unsigned int Elem::_bp1 = 65449;
const unsigned int Elem::_bp2 = 48661;
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
	std::cerr << "ERROR: Undefined element type!." << std::endl;
	libmesh_error();
      }
    }
  

  AutoPtr<Elem> ap(elem);
  return ap;
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
  Real h_min=1.e30;

  for (unsigned int n_outer=0; n_outer<this->n_vertices(); n_outer++)
    for (unsigned int n_inner=n_outer+1; n_inner<this->n_vertices(); n_inner++)
      {
	const Point diff = (this->point(n_outer) - this->point(n_inner));
	
	h_min = std::min(h_min,diff.size());
      }

  return h_min;
}



Real Elem::hmax() const
{
  Real h_max=0;

  for (unsigned int n_outer=0; n_outer<this->n_vertices(); n_outer++)
    for (unsigned int n_inner=n_outer+1; n_inner<this->n_vertices(); n_inner++)
      {
	const Point diff = (this->point(n_outer) - this->point(n_inner));
	
	h_max = std::max(h_max,diff.size());
      }

  return h_max;
}



Real Elem::length(const unsigned int n1, 
		  const unsigned int n2) const
{
  libmesh_assert ( n1 < this->n_vertices() );
  libmesh_assert ( n2 < this->n_vertices() );

  return (this->point(n1) - this->point(n2)).size();
}



bool Elem::operator == (const DofObject& rhs) const
{

    // Cast rhs to an Elem*
    const Elem* rhs_elem = dynamic_cast<const Elem*>(&rhs);

    // If we cannot cast to an Elem*, rhs must be a Node
    if(rhs_elem == static_cast<const Elem*>(NULL))
        return false;

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
  typedef std::vector<unsigned int>::iterator iterator;

  
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
      std::vector<unsigned int> common_nodes;
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
      
      const int new_size = std::distance (common_nodes.begin(),
					  new_end);
      
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



bool Elem::contains_vertex_of(const Elem *e) const
{
  // Our vertices are the first numbered nodes
  for (unsigned int n = 0; n != e->n_vertices(); ++n)
    if (this->contains_point(e->point(n)))
      return true;
  return false;
}



void Elem::find_point_neighbors(std::set<const Elem *> &neighbor_set) const
{
  neighbor_set.clear();
  neighbor_set.insert(this);

  unsigned int old_size;
  do
    {
      old_size = neighbor_set.size();

      // Loop over all the elements in the patch
      std::set<const Elem*>::const_iterator       it  = neighbor_set.begin();
      const std::set<const Elem*>::const_iterator end = neighbor_set.end();

      for (; it != end; ++it)
        {
          const Elem* elem = *it;

          for (unsigned int s=0; s<elem->n_sides(); s++)
            {
              const Elem* neighbor = elem->neighbor(s);
              if (neighbor &&
                  neighbor != remote_elem)    // we have a real neighbor on this side
                {
                  if (neighbor->active())                // ... if it is active
                    {
                      if (this->contains_vertex_of(neighbor) // ... and touches us
                          || neighbor->contains_vertex_of(this))  
                        neighbor_set.insert (neighbor);  // ... then add it
                    }
#ifdef LIBMESH_ENABLE_AMR
                  else                                 // ... the neighbor is *not* active,
                    {                                  // ... so add *all* neighboring
                                                       // active children
                      std::vector<const Elem*> active_neighbor_children;
  
                      neighbor->active_family_tree_by_neighbor
                        (active_neighbor_children, elem);

                      std::vector<const Elem*>::const_iterator
                        child_it = active_neighbor_children.begin();
                      const std::vector<const Elem*>::const_iterator
                        child_end = active_neighbor_children.end();
                      for (; child_it != child_end; ++child_it)
                        if (this->contains_vertex_of(*child_it) ||
                            (*child_it)->contains_vertex_of(this))
                          neighbor_set.insert (*child_it);
                    }
#endif // #ifdef LIBMESH_ENABLE_AMR
                }
            }
        }
    }
  while (old_size != neighbor_set.size());
}



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
          // ancestor
          if (this->subactive() && !neigh->subactive())
            {
              for (elem = this; !elem->active();
                   elem = elem->parent())
                libmesh_assert(elem);
            }

          unsigned int rev = neigh->which_neighbor_am_i(elem);
          libmesh_assert (rev < neigh->n_neighbors());

          if (this->subactive() && !neigh->subactive())
            {
              libmesh_assert (neigh->neighbor(rev) == elem);
            }
          else
            {
              Elem *nn = neigh->neighbor(rev);
              libmesh_assert(nn);

              for (; elem != nn; elem = elem->parent())
                libmesh_assert(elem);
            }
        }
      else
        {
          const Elem *parent = this->parent();
          if (parent)
            libmesh_assert (!parent->neighbor(s));
        }
    }
}

#endif // DEBUG



void Elem::make_links_to_me_remote()
{
  libmesh_assert (this != remote_elem);

  // We need to have handled any children first
#if defined(LIBMESH_ENABLE_AMR) && defined(DEBUG)
  if (this->has_children())
    for (unsigned int c = 0; c != this->n_children(); ++c)
      {
        Elem *child = this->child(c);
        libmesh_assert (child == remote_elem);
      }
#endif

  // Remotify any neighbor links
  for (unsigned int s = 0; s != this->n_sides(); ++s)
    {
      Elem *neigh = this->neighbor(s);
      if (neigh && neigh != remote_elem)
        {
	  // My neighbor should never be more refined than me; my real
	  // neighbor would have been its parent in that case.
	  libmesh_assert(this->level() >= neigh->level());
	  
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
                  libmesh_assert (my_s < n->n_neighbors());
                  libmesh_assert (n->neighbor(my_s) == this);
                  n->set_neighbor(my_s, const_cast<RemoteElem*>(remote_elem));
                }
#else
              unsigned int my_s = neigh->which_neighbor_am_i(this);
              libmesh_assert (my_s < neigh->n_neighbors());
              libmesh_assert (neigh->neighbor(my_s) == this);
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
	      Elem *ancestor = this->parent();
	      libmesh_assert(ancestor);
              while (!neigh->has_neighbor(ancestor))
                {
                  ancestor = ancestor->parent();
	          libmesh_assert(ancestor);
                }

	      // My neighbor may have descendants which consider me a
	      // neighbor
              std::vector<const Elem*> family;
              neigh->family_tree_by_subneighbor (family, ancestor, this);

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
                  libmesh_assert (my_s < n->n_neighbors());
                  libmesh_assert (n->neighbor(my_s) == this);
                  n->set_neighbor(my_s, const_cast<RemoteElem*>(remote_elem));
                }
            }
#endif
        }
    }

#ifdef LIBMESH_ENABLE_AMR
  // Remotify parent's child link
  Elem *parent = this->parent();
  if (parent && parent != remote_elem)
    {
      unsigned int me = parent->which_child_am_i(this);
      libmesh_assert (parent->_children[me] == this);
      parent->_children[me] = const_cast<RemoteElem*>(remote_elem);
    }
#endif
}



void Elem::write_connectivity (std::ostream& out,
			       const IOPackage iop) const
{
  libmesh_assert (out.good());
  libmesh_assert (_nodes != NULL);
  libmesh_assert (iop != INVALID_IO_PACKAGE);

  switch (iop)
    {
    case TECPLOT:
      {
	// This connectivity vector will be used repeatedly instead
	// of being reconstructed inside the loop.
	std::vector<unsigned int> conn;
	for (unsigned int sc=0; sc <this->n_sub_elem(); sc++)
	  {
	    this->connectivity(sc, TECPLOT, conn);
	    
	    std::copy(conn.begin(),
		      conn.end(),
		      std::ostream_iterator<unsigned int>(out, " "));
	    
	    out << '\n';
	  }
	return;
      }

    case UCD:
      {
	for (unsigned int i=0; i<this->n_nodes(); i++)
	  out << this->node(i)+1 << "\t";
	
	out << '\n';
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
//   libmesh_assert (_nodes != NULL);

//   // This connectivity vector will be used repeatedly instead
//   // of being reconstructed inside the loop.
//   std::vector<unsigned int> conn;
//   for (unsigned int sc=0; sc <this->n_sub_elem(); sc++)
//     {
//       this->connectivity(sc, TECPLOT, conn);
      
//       std::copy(conn.begin(),
//  		conn.end(),
//  		std::ostream_iterator<unsigned int>(out, " "));
      
//       out << std::endl;
//     }
// }



// void Elem::write_ucd_connectivity(std::ostream &out) const
// {
//   libmesh_assert (out);
//   libmesh_assert (_nodes != NULL);

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
	here();

	std::cerr << "ERROR:  unknown quality metric: "
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
	_children[c] = NULL;
    }
  
  for (unsigned int c=0; c<this->n_children(); c++)
    {
      if(_children[c] == NULL || _children[c] == remote_elem)
	{
	  libmesh_assert (this == elem->parent());
	  _children[c] = elem;
	  return;
	}
    }

  std::cerr << "Error: Tried to add a child to an element with full children array"
            << std::endl;
  libmesh_error();
}



void Elem::add_child (Elem* elem, unsigned int c)
{
  if(_children == NULL)
    {
      _children = new Elem*[this->n_children()];
      
      for (unsigned int i=0; i<this->n_children(); i++)
	_children[i] = NULL;
    }
  
  libmesh_assert (_children[c] == NULL || _children[c] == remote_elem);
  libmesh_assert (this == elem->parent());

  _children[c] = elem;
}



bool Elem::is_child_on_edge(const unsigned int c,
                            const unsigned int e) const
{
  libmesh_assert (c < this->n_children());
  libmesh_assert (e < this->n_edges());

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

  libmesh_assert(s < this->n_sides());

  // Add this element to the family tree.
  family.push_back(this);

  // Recurse into the elements children, if it has them.
  // Do not clear the vector any more.
  if (!this->active())
    for (unsigned int c=0; c<this->n_children(); c++)
      if (this->child(c)->is_child_on_side(c, s))
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

  libmesh_assert(s < this->n_sides());

  // Add an active element to the family tree.
  if (this->active())
    family.push_back(this);

  // Or recurse into an ancestor element's children.
  // Do not clear the vector any more.
  else
    for (unsigned int c=0; c<this->n_children(); c++)
      if (this->child(c)->is_child_on_side(c, s))
        this->child(c)->active_family_tree_by_side (family, s, false);
}



void Elem::family_tree_by_neighbor (std::vector<const Elem*>& family,
                                    const Elem* neighbor,
			            const bool reset) const
{
  // The "family tree" doesn't include subactive elements
  libmesh_assert(!this->subactive());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // This only makes sense if we're already a neighbor
  libmesh_assert (this->has_neighbor(neighbor));

  // Add this element to the family tree.
  family.push_back(this);

  // Recurse into the elements children, if it's not active.
  // Do not clear the vector any more.
  if (!this->active())
    for (unsigned int c=0; c<this->n_children(); c++)
      {
        Elem *child = this->child(c);
        if (child != remote_elem && child->has_neighbor(neighbor))
          child->family_tree_by_neighbor (family, neighbor, false);
      }
}



void Elem::family_tree_by_subneighbor (std::vector<const Elem*>& family,
                                       const Elem* neighbor,
                                       const Elem* subneighbor,
			               const bool reset) const
{
  // The "family tree" doesn't include subactive elements
  libmesh_assert(!this->subactive());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // To simplifly this function we need an existing neighbor
  libmesh_assert (neighbor);
  libmesh_assert (neighbor != remote_elem);
  libmesh_assert (this->has_neighbor(neighbor));

  // This only makes sense if subneighbor descends from neighbor
  libmesh_assert (subneighbor);
  libmesh_assert (subneighbor != remote_elem);
  libmesh_assert (neighbor->is_ancestor_of(subneighbor));

  // Add this element to the family tree if applicable.
  if (neighbor == subneighbor)
    family.push_back(this);

  // Recurse into the elements children, if it's not active.
  // Do not clear the vector any more.
  if (!this->active())
    for (unsigned int c=0; c != this->n_children(); ++c)
      {
        Elem *child = this->child(c);
        if (child != remote_elem)
          for (unsigned int s=0; s != child->n_sides(); ++s)
            {
	      Elem *child_neigh = child->neighbor(s);
              if (child_neigh &&
		  (child_neigh == neighbor ||
		   (child_neigh->parent() == neighbor &&
		    child_neigh->is_ancestor_of(subneighbor))))
                child->family_tree_by_subneighbor (family, child_neigh,
                                                   subneighbor, false);
            }
      }
}



void Elem::active_family_tree_by_neighbor (std::vector<const Elem*>& family,
                                           const Elem* neighbor,
			                   const bool reset) const
{
  // The "family tree" doesn't include subactive elements
  libmesh_assert(!this->subactive());

  // Clear the vector if the flag reset tells us to.
  if (reset)
    family.clear();

  // This only makes sense if we're already a neighbor
  if (this->level() >= neighbor->level())
    libmesh_assert (this->has_neighbor(neighbor));

  // Add an active element to the family tree.
  if (this->active())
    family.push_back(this);

  // Or recurse into an ancestor element's children.
  // Do not clear the vector any more.
  else if (!this->active())
    for (unsigned int c=0; c<this->n_children(); c++)
      {
        Elem *child = this->child(c);
        if (child != remote_elem && child->has_neighbor(neighbor))
          child->active_family_tree_by_neighbor (family, neighbor, false);
      }
}



unsigned int Elem::min_p_level_by_neighbor(const Elem* neighbor,
                                           unsigned int current_min) const
{
  libmesh_assert(!this->subactive());
  libmesh_assert(neighbor->active());

  // If we're an active element this is simple
  if (this->active())
    return std::min(current_min, this->p_level());

  libmesh_assert(has_neighbor(neighbor));

  // The p_level() of an ancestor element is already the minimum
  // p_level() of its children - so if that's high enough, we don't
  // need to examine any children.
  if (current_min <= this->p_level())
    return current_min;

  unsigned int min_p_level = current_min;

  for (unsigned int c=0; c<this->n_children(); c++)
    {
      const Elem* const child = this->child(c);
      if (child != remote_elem && child->has_neighbor(neighbor))
        min_p_level =
	  child->min_p_level_by_neighbor(neighbor,
                                         min_p_level);
    }

  return min_p_level;
}


unsigned int Elem::min_new_p_level_by_neighbor(const Elem* neighbor,
                                               unsigned int current_min) const
{
  libmesh_assert(!this->subactive());
  libmesh_assert(neighbor->active());

  // If we're an active element this is simple
  if (this->active())
    {
      unsigned int new_p_level = this->p_level();
      if (this->p_refinement_flag() == Elem::REFINE)
        new_p_level += 1;
      if (this->p_refinement_flag() == Elem::COARSEN)
        {
          libmesh_assert (new_p_level > 0);
          new_p_level -= 1;
        }
      return std::min(current_min, new_p_level);
    }

  libmesh_assert(has_neighbor(neighbor));

  unsigned int min_p_level = current_min;

  for (unsigned int c=0; c<this->n_children(); c++)
    {
      const Elem* const child = this->child(c);
      if (child && child != remote_elem)
        if (child->has_neighbor(neighbor))
          min_p_level =
	    child->min_new_p_level_by_neighbor(neighbor,
                                               min_p_level);
    }

  return min_p_level;
}

#endif // #ifdef LIBMESH_ENABLE_AMR



bool Elem::contains_point (const Point& p) const
{
  // Declare a basic FEType.  Will ue a Lagrange
  // element by default.
  FEType fe_type(this->default_order());
  
  const Point mapped_point = FEInterface::inverse_map(this->dim(),
						      fe_type,
						      this,
						      p,
						      1.e-4,
						      false);

  return FEInterface::on_reference_element(mapped_point, this->type());
}



void Elem::nullify_neighbors ()
{
  // Tell any of my neighbors about my death...
  // Looks strange, huh?
  for (unsigned int n=0; n<this->n_neighbors(); n++)
    {
      Elem* neighbor = this->neighbor(n);
      if (neighbor && neighbor != remote_elem)
        {
	  // Note:  it is possible that I see the neighbor
	  // (which is coarser than me)
	  // but they don't see me, so avoid that case.
	  if (neighbor->level() == this->level())
	    {	
	      const unsigned int w_n_a_i = neighbor->which_neighbor_am_i(this);
              libmesh_assert (w_n_a_i < neighbor->n_neighbors());
	      neighbor->set_neighbor(w_n_a_i, NULL);
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
void Elem::PackedElem::pack (std::vector<int> &conn, const Elem* elem)
{
  libmesh_assert (elem != NULL);
  
  // we can do at least this good. note that hopefully in general
  // the user will already have reserved the full space, which will render
  // this redundant
  conn.reserve (conn.size() + Elem::PackedElem::header_size + elem->n_nodes());

#ifdef LIBMESH_ENABLE_AMR
  conn.push_back (static_cast<int>(elem->level()));
  conn.push_back (static_cast<int>(elem->p_level()));
  conn.push_back (static_cast<int>(elem->refinement_flag()));
  conn.push_back (static_cast<int>(elem->p_refinement_flag()));
#else
  conn.push_back (0);
  conn.push_back (0);
  conn.push_back (0);
  conn.push_back (0);
#endif
  conn.push_back (static_cast<int>(elem->type()));
  conn.push_back (static_cast<int>(elem->processor_id()));
  conn.push_back (static_cast<int>(elem->subdomain_id()));
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
}



Elem * Elem::PackedElem::unpack (MeshBase &mesh, Elem *parent) const
{  
  
  Elem *elem = Elem::build(this->type(),parent).release();
  libmesh_assert (elem);

#ifdef LIBMESH_ENABLE_AMR
  if (this->level() != 0) 
    {
      libmesh_assert (parent != NULL);
      parent->add_child(elem, this->which_child_am_i());
      libmesh_assert (parent->type() == elem->type());
      libmesh_assert (parent->child(this->which_child_am_i()) == elem);
    }
#endif

  // Assign the IDs
#ifdef LIBMESH_ENABLE_AMR
  elem->set_p_level(this->p_level());
  elem->set_refinement_flag(this->refinement_flag());
  elem->set_p_refinement_flag(this->p_refinement_flag());
  libmesh_assert (elem->level() == this->level());
#endif
  elem->subdomain_id() = this->subdomain_id();
  elem->processor_id() = this->processor_id();
  elem->set_id()       = this->id();
  
  // Assign the connectivity
  libmesh_assert (elem->n_nodes() == this->n_nodes());

  for (unsigned int n=0; n<elem->n_nodes(); n++)
    elem->set_node(n) = mesh.node_ptr (this->node(n));
  
  return elem;
}
