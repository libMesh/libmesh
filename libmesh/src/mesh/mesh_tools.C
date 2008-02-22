// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2007  Benjamin S. Kirk, John W. Peterson
  
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
#include <limits>
#include <set>

// Local includes
#include "elem.h"
#include "mesh_tools.h"
#include "mesh_base.h"
#include "parallel.h"
#include "sphere.h"
#include "serial_mesh.h"
#include "parallel_mesh.h"
#include "mesh_communication.h"
#include "elem_range.h"
#include "node_range.h"
#include "threads.h"



// ------------------------------------------------------------
// anonymous namespace for helper classes
namespace {
  /**
   * SumElemWeight(Range) sums the number of nodes per element 
   * for each element in the provided range. The join() method
   * defines how to combine the reduction operation from two
   * distinct instances of this class which may be executed on 
   * separate threads.
   */
  class SumElemWeight
  {
  public:
    SumElemWeight () : 
      _weight(0) 
    {}

    SumElemWeight (SumElemWeight &, Threads::split) : 
      _weight(0)
    {}
    
    void operator()(const ConstElemRange &range)
    { 
      for (ConstElemRange::const_iterator it = range.begin(); it !=range.end(); ++it)
	_weight += (*it)->n_nodes();
    }

    unsigned int weight() const 
    { return _weight; }
    
    void join (const SumElemWeight &other)
    { _weight += other.weight(); }

  private:
    unsigned int _weight;
  };


  /**
   * FindBBox(Range) computes the bounding box for the objects
   * in the specified range.  This class may be split and subranges
   * can be executed on separate threads.  The join() method
   * defines how the results from two separate threads are combined.
   */
  class FindBBox
  {
  public:
    FindBBox () :
      _vmin(3,  std::numeric_limits<Real>::max()),
      _vmax(3, -std::numeric_limits<Real>::max())
    {}

    FindBBox (FindBBox &other, Threads::split) :
      _vmin(other._vmin),
      _vmax(other._vmax)
    {}

    std::vector<Real> & min() { return _vmin; }
    std::vector<Real> & max() { return _vmax; }
    
    void operator()(const ConstNodeRange &range)
    {
      for (ConstNodeRange::const_iterator it = range.begin(); it != range.end(); ++it)
	{
          const Node *node = *it;
	  assert (node != NULL);
	  
	  for (unsigned int i=0; i<3; i++)
	    {
	      _vmin[i] = std::min(_vmin[i], (*node)(i));
	      _vmax[i] = std::max(_vmax[i], (*node)(i));
	    }      
        }
    }

    void operator()(const ConstElemRange &range)
    {
      for (ConstElemRange::const_iterator it = range.begin(); it != range.end(); ++it)
	{
          const Elem *elem = *it;
	  assert (elem != NULL);

	  for (unsigned int n=0; n<elem->n_nodes(); n++)
	    for (unsigned int i=0; i<3; i++)
	      {
		_vmin[i] = std::min(_vmin[i], elem->point(n)(i));
		_vmax[i] = std::max(_vmax[i], elem->point(n)(i));
	      }      
        }
    }

    void join (const FindBBox &other)
    {
      for (unsigned int i=0; i<3; i++)
	{
	  _vmin[i] = std::min(_vmin[i], other._vmin[i]);
	  _vmax[i] = std::max(_vmax[i], other._vmax[i]);
	}      
    }

    MeshTools::BoundingBox bbox () const
    {
      Point pmin(_vmin[0], _vmin[1], _vmin[2]);
      Point pmax(_vmax[0], _vmax[1], _vmax[2]);

      const MeshTools::BoundingBox ret_val(pmin, pmax);

      return ret_val;        
    }
    
  private:
    std::vector<Real> _vmin;
    std::vector<Real> _vmax;
  };
}



// ------------------------------------------------------------
// MeshTools functions
unsigned int MeshTools::total_weight(const MeshBase& mesh)
{
  if (!mesh.is_serial())
    {
      parallel_only();
      unsigned int weight = MeshTools::weight (mesh, libMesh::processor_id());
      Parallel::sum(weight);
      unsigned int unpartitioned_weight =
        MeshTools::weight (mesh, DofObject::invalid_processor_id);
      return weight + unpartitioned_weight;
    }
  
  SumElemWeight sew;
  
  Threads::parallel_reduce (ConstElemRange (mesh.elements_begin(),
					    mesh.elements_end()),
			    sew);
  return sew.weight();

}



unsigned int MeshTools::weight(const MeshBase& mesh, const unsigned int pid)
{
  SumElemWeight sew;
  
  Threads::parallel_reduce (ConstElemRange (mesh.pid_elements_begin(pid),
					    mesh.pid_elements_end(pid)),
			    sew);
  return sew.weight();
}



void MeshTools::build_nodes_to_elem_map (const MeshBase& mesh,
					 std::vector<std::vector<unsigned int> >& nodes_to_elem_map)
{
  nodes_to_elem_map.resize (mesh.n_nodes());

  MeshBase::const_element_iterator       el  = mesh.elements_begin();
  const MeshBase::const_element_iterator end = mesh.elements_end();

  for (; el != end; ++el)
    for (unsigned int n=0; n<(*el)->n_nodes(); n++)
      {
	assert ((*el)->node(n) < nodes_to_elem_map.size());
	assert ((*el)->id()    < mesh.n_elem());
	
	nodes_to_elem_map[(*el)->node(n)].push_back((*el)->id());
      }
}



void MeshTools::build_nodes_to_elem_map (const MeshBase& mesh,
					 std::vector<std::vector<const Elem*> >& nodes_to_elem_map)
{
  nodes_to_elem_map.resize (mesh.n_nodes());

  MeshBase::const_element_iterator       el  = mesh.elements_begin();
  const MeshBase::const_element_iterator end = mesh.elements_end();

  for (; el != end; ++el)
    for (unsigned int n=0; n<(*el)->n_nodes(); n++)
      {
	assert ((*el)->node(n) < nodes_to_elem_map.size());
	
	nodes_to_elem_map[(*el)->node(n)].push_back(*el);
      }
}



void MeshTools::find_boundary_nodes (const MeshBase& mesh,
				     std::vector<bool>& on_boundary)
{
  // Resize the vector which holds boundary nodes and fill with false.
  on_boundary.resize(mesh.n_nodes());
  std::fill(on_boundary.begin(),
	    on_boundary.end(),
	    false);

  // Loop over elements, find those on boundary, and
  // mark them as true in on_boundary.
  MeshBase::const_element_iterator       el  = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_elements_end(); 

  for (; el != end; ++el)
    for (unsigned int s=0; s<(*el)->n_neighbors(); s++)
      if ((*el)->neighbor(s) == NULL) // on the boundary
	{
	  const AutoPtr<Elem> side((*el)->build_side(s));
	  
	  for (unsigned int n=0; n<side->n_nodes(); n++)
	    on_boundary[side->node(n)] = true;
	}
}



MeshTools::BoundingBox
MeshTools::bounding_box(const MeshBase& mesh)
{
  // This function must be run on all processors at once
  parallel_only();

  FindBBox find_bbox;
  
  Threads::parallel_reduce (ConstNodeRange (mesh.local_nodes_begin(),
					    mesh.local_nodes_end()),
			    find_bbox);

  // and the unpartitioned nodes
  Threads::parallel_reduce (ConstNodeRange (mesh.pid_nodes_begin(DofObject::invalid_processor_id),
					    mesh.pid_nodes_end(DofObject::invalid_processor_id)),
			    find_bbox);

  // Compare the bounding boxes across processors
  Parallel::min(find_bbox.min());
  Parallel::max(find_bbox.max());

  return find_bbox.bbox();
}



Sphere
MeshTools::bounding_sphere(const MeshBase& mesh)
{
  BoundingBox bbox = bounding_box(mesh);

  const Real  diag = (bbox.second - bbox.first).size();
  const Point cent = (bbox.second + bbox.first)/2.;

  return Sphere (cent, .5*diag);
}



MeshTools::BoundingBox
MeshTools::processor_bounding_box (const MeshBase& mesh,
				   const unsigned int pid)
{
  assert (pid < libMesh::n_processors());

  FindBBox find_bbox;
  
  Threads::parallel_reduce (ConstElemRange (mesh.pid_elements_begin(pid),
					    mesh.pid_elements_end(pid)),
			    find_bbox);
  
  return find_bbox.bbox();
}



Sphere
MeshTools::processor_bounding_sphere (const MeshBase& mesh,
				      const unsigned int pid)
{
  BoundingBox bbox = processor_bounding_box(mesh,pid);

  const Real  diag = (bbox.second - bbox.first).size();
  const Point cent = (bbox.second + bbox.first)/2.;

  return Sphere (cent, .5*diag);
}



MeshTools::BoundingBox
MeshTools::subdomain_bounding_box (const MeshBase& mesh,
				   const unsigned int sid)
{
  assert (mesh.n_nodes() != 0);

  Point min( 1.e30,  1.e30,  1.e30);
  Point max(-1.e30, -1.e30, -1.e30);

  for (unsigned int e=0; e<mesh.n_elem(); e++)
    if (mesh.elem(e)->subdomain_id() == sid)
      for (unsigned int n=0; n<mesh.elem(e)->n_nodes(); n++)
	for (unsigned int i=0; i<mesh.spatial_dimension(); i++)
	  {
	    min(i) = std::min(min(i), mesh.point(mesh.elem(e)->node(n))(i));
	    max(i) = std::max(max(i), mesh.point(mesh.elem(e)->node(n))(i));
	  }      
  
  const BoundingBox ret_val(min, max);
  
  return ret_val;  
}



Sphere
MeshTools::subdomain_bounding_sphere (const MeshBase& mesh,
				      const unsigned int sid)
{
  BoundingBox bbox = subdomain_bounding_box(mesh,sid);

  const Real  diag = (bbox.second - bbox.first).size();
  const Point cent = (bbox.second + bbox.first)/2.;

  return Sphere (cent, .5*diag);
}



void MeshTools::elem_types (const MeshBase& mesh,
			    std::vector<ElemType>& et)
{
  MeshBase::const_element_iterator       el  = mesh.elements_begin();
  const MeshBase::const_element_iterator end = mesh.elements_end(); 

  // Automatically get the first type
  et.push_back((*el)->type());  ++el;

  // Loop over the rest of the elements.
  // If the current element type isn't in the
  // vector, insert it.
  for (; el != end; ++el)
    if (!std::count(et.begin(), et.end(), (*el)->type()))
      et.push_back((*el)->type());
}



unsigned int MeshTools::n_elem_of_type (const MeshBase& mesh,
					const ElemType type)
{
  return static_cast<unsigned int>(std::distance(mesh.type_elements_begin(type),
						 mesh.type_elements_end  (type)));
}



unsigned int MeshTools::n_active_elem_of_type (const MeshBase& mesh,
					       const ElemType type)
{
  return static_cast<unsigned int>(std::distance(mesh.active_type_elements_begin(type),
						 mesh.active_type_elements_end  (type)));
}

unsigned int MeshTools::n_non_subactive_elem_of_type_at_level(const MeshBase& mesh,
                                                const ElemType type,
                                                const unsigned int level)
{
  unsigned int cnt = 0;
  // iterate over the elements of the specified type
  MeshBase::const_element_iterator el = mesh.type_elements_begin(type);
  const MeshBase::const_element_iterator end = mesh.type_elements_end(type);

  for(; el!=end; ++el)
    if( ((*el)->level() == level) && !(*el)->subactive())
      cnt++;

  return cnt;
}


unsigned int MeshTools::n_active_local_levels(const MeshBase& mesh)
{
  unsigned int max_level = 0;

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for( ; el != end_el; ++el)
    max_level = std::max((*el)->level(), max_level);

  return max_level + 1;
}



unsigned int MeshTools::n_active_levels(const MeshBase& mesh)
{
  parallel_only();

  unsigned int nl = MeshTools::n_active_local_levels(mesh);

  MeshBase::const_element_iterator el =
    mesh.active_pid_elements_begin(DofObject::invalid_processor_id);
  const MeshBase::const_element_iterator end_el =
    mesh.active_pid_elements_end(DofObject::invalid_processor_id);

  for( ; el != end_el; ++el)
    nl = std::max((*el)->level() + 1, nl);

  Parallel::max(nl);
  return nl;
}



unsigned int MeshTools::n_local_levels(const MeshBase& mesh)
{
  unsigned int max_level = 0;

  MeshBase::const_element_iterator el = mesh.local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.local_elements_end();

  for( ; el != end_el; ++el)
    max_level = std::max((*el)->level(), max_level);

  return max_level + 1;
}
   


unsigned int MeshTools::n_levels(const MeshBase& mesh)
{
  parallel_only();

  unsigned int nl = MeshTools::n_local_levels(mesh);

  MeshBase::const_element_iterator el =
    mesh.pid_elements_begin(DofObject::invalid_processor_id);
  const MeshBase::const_element_iterator end_el =
    mesh.pid_elements_end(DofObject::invalid_processor_id);

  for( ; el != end_el; ++el)
    nl = std::max((*el)->level() + 1, nl);

  Parallel::max(nl);
  return nl;
}



void MeshTools::get_not_subactive_node_ids(const MeshBase& mesh,
    std::set<unsigned int>& not_subactive_node_ids)
{
  MeshBase::const_element_iterator el           = mesh.elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.elements_end();
  for( ; el != end_el; ++el)
  {
    Elem* elem = (*el);
    if(!elem->subactive())
      for (unsigned int n=0; n<elem->n_nodes(); ++n)
        not_subactive_node_ids.insert(elem->node(n));
  }
}



unsigned int MeshTools::n_elem (MeshBase::const_element_iterator& begin,
                                MeshBase::const_element_iterator& end)
{
  return std::distance(begin, end);
}
   


unsigned int MeshTools::n_p_levels (const MeshBase& mesh)
{
  parallel_only();
  
  unsigned int max_p_level = 0;

  // first my local elements
  MeshBase::const_element_iterator
    el     = mesh.local_elements_begin(),
    end_el = mesh.local_elements_end();

  for( ; el != end_el; ++el)
    max_p_level = std::max((*el)->p_level(), max_p_level);

  // then any unpartitioned objects
  el     = mesh.pid_elements_begin(DofObject::invalid_processor_id);
  end_el = mesh.pid_elements_end(DofObject::invalid_processor_id);

  for( ; el != end_el; ++el)
    max_p_level = std::max((*el)->p_level(), max_p_level);

  Parallel::max(max_p_level);
  return max_p_level + 1;
}



void MeshTools::find_nodal_neighbors(const MeshBase&, const Node& n, 
                                     std::vector<std::vector<const Elem*> >& nodes_to_elem_map, 
                                     std::vector<const Node*>& neighbors)
{
  unsigned int global_id = n.id();
  
  //Iterators to iterate through the elements that include this node
  std::vector<const Elem*>::const_iterator el     = nodes_to_elem_map[global_id].begin();
  std::vector<const Elem*>::const_iterator end_el = nodes_to_elem_map[global_id].end();
  
  unsigned int n_ed=0; //Number of edges on the element
  unsigned int ed=0; //Current edge
  unsigned int l_n=0; //Local node number
  unsigned int o_n=0; //Other node on this edge
  
  //Assume we find a edge... then prove ourselves wrong...
  bool found_edge=true;
  
  Node * node_to_save = NULL;
  
  //Look through the elements that contain this node
  //find the local node id... then find the side that
  //node lives on in the element
  //next, look for the _other_ node on that side
  //That other node is a "nodal_neighbor"... save it
  for(;el != end_el;el++)
  {
    //We only care about active elements...
    if((*el)->active())
    {
      n_ed=(*el)->n_edges();
      
      //Find the local node id
      while(global_id != (*el)->node(l_n++)) { }
      l_n--; //Hmmm... take the last one back off
      
      while(ed<n_ed)
      {
        
        //Find the edge the node is on
        while(found_edge && !(*el)->is_node_on_edge(l_n,ed++))
        {
          //This only happens if all the edges have already been found
          if(ed>=n_ed)
            found_edge=false;
        }
        
        //Did we find one?
        if(found_edge)
        {
          ed--; //Take the last one back off again
          
          //Now find the other node on that edge
          while(!(*el)->is_node_on_edge(o_n++,ed) || global_id==(*el)->node(o_n-1)) { }
          o_n--;
          
          //We've found one!  Save it..
          node_to_save=(*el)->get_node(o_n);
          
          //Search to see if we've already found this one
          std::vector<const Node*>::const_iterator result = std::find(neighbors.begin(),neighbors.end(),node_to_save);
          
          //If we didn't find it and add it to the vector
          if(result == neighbors.end())
            neighbors.push_back(node_to_save);
        }
        
        //Reset to look for another
        o_n=0;
        
        //Keep looking for edges, node may be on more than one edge
        ed++;
      }
      
      //Reset to get ready for the next element
      l_n=ed=0;
      found_edge=true;
    }
  }
}

void MeshTools::find_hanging_nodes_and_parents(const MeshBase& mesh, std::map<unsigned int, std::vector<unsigned int> >& hanging_nodes)
{
  MeshBase::const_element_iterator it  = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end = mesh.active_local_elements_end();
  
  //Loop through all the elements
  for (; it != end; ++it)
  {
    //Save it off for easier access
    const Elem* elem = (*it);
    
    //Right now this only works for quad4's
    //assert(elem->type() == libMeshEnums::QUAD4);
    if(elem->type() == libMeshEnums::QUAD4)
    {
      //Loop over the sides looking for sides that have hanging nodes
      //This code is inspired by compute_proj_constraints()
      for (unsigned int s=0; s<elem->n_sides(); s++)
      {
        //If not a boundary node
        if (elem->neighbor(s) != NULL)
        {
          // Get pointers to the element's neighbor.
          const Elem* neigh = elem->neighbor(s);
          
          //Is there a coarser element next to this one?
          if (neigh->level() < elem->level()) 
          {
            const Elem *ancestor = elem;
            while (neigh->level() < ancestor->level())
              ancestor = ancestor->parent();
            unsigned int s_neigh = neigh->which_neighbor_am_i(ancestor);
            assert (s_neigh < neigh->n_neighbors());
            
            //Couple of helper uints...
            unsigned int node1=0;
            unsigned int node2=0;
            unsigned int hanging_node=0;
                
            bool found_in_neighbor = false;
                
            //Find the two vertices that make up this side
            while(!elem->is_node_on_side(node1++,s)) { }
            node1--;
                
                //Start looking for the second one with the next node
            node2=node1+1;
            
            //Find the other one
            while(!elem->is_node_on_side(node2++,s)) { }
            node2--;
  
                //Pull out their global ids:
            node1 = elem->node(node1);
            node2 = elem->node(node2);
            
            //Now find which node is present in the neighbor
            //FIXME This assumes a level one rule!
            //The _other_ one is the hanging node
            
            //First look for the first one
            //FIXME could be streamlined a bit
            for(unsigned int n=0;n<neigh->n_sides();n++)
            {
              if(neigh->node(n) == node1)
                found_in_neighbor=true;
            }
                
                
            if(!found_in_neighbor)
              hanging_node=node1;
            else //If it wasn't node1 then it must be node2!
              hanging_node=node2;
            
            //Reset these for reuse
            node1=0;
            node2=0;
            
            //Find the first node that makes up the side in the neighbor (these should be the parent nodes)
            while(!neigh->is_node_on_side(node1++,s_neigh)) { }
            node1--;
                
            node2=node1+1;
            
            //Find the second node...
            while(!neigh->is_node_on_side(node2++,s_neigh)) { }
            node2--;
            
            //Save them if we haven't already found the parents for this one
            if(hanging_nodes[hanging_node].size()<2)
            {
              hanging_nodes[hanging_node].push_back(neigh->node(node1));
              hanging_nodes[hanging_node].push_back(neigh->node(node2));
            }
          }
        }
      }
    }
  }
}



void MeshTools::Private::globally_renumber_nodes_and_elements (MeshBase& mesh)
{
  MeshCommunication().assign_global_indices(mesh);
}



void MeshTools::Private::fix_broken_node_and_element_numbering (SerialMesh &mesh)
{
   // Nodes first
  for (unsigned int n=0; n<mesh._nodes.size(); n++)
    if (mesh._nodes[n] != NULL)
      mesh._nodes[n]->set_id() = n;

  // Elements next
  for (unsigned int e=0; e<mesh._elements.size(); e++)
    if (mesh._elements[e] != NULL)
      mesh._elements[e]->set_id() = e; 
}



void MeshTools::Private::fix_broken_node_and_element_numbering (ParallelMesh &mesh)
{
  // We need access to iterators for the underlying containers,
  // not the mapvector<> reimplementations.
  mapvector<Node*>::maptype &nodes = mesh._nodes;
  mapvector<Elem*>::maptype &elem  = mesh._elements;
  
  // Nodes first
  {
    mapvector<Node*>::maptype::iterator
      it  = nodes.begin(),
      end = nodes.end();

    for (; it != end; ++it)
      if (it->second != NULL)
	it->second->set_id() = it->first;
  }

  // Elements next
  {
    mapvector<Elem*>::maptype::iterator
      it  = elem.begin(),
      end = elem.end();

    for (; it != end; ++it)
      if (it->second != NULL)
	it->second->set_id() = it->first;
  }
}
