// $Id: mesh_tools.C,v 1.2 2005-02-22 22:17:41 jwpeterson Exp $

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



// C++ includes

// Local includes
#include "mesh_tools.h"
#include "mesh_base.h"
#include "elem.h"
#include "sphere.h"



// ------------------------------------------------------------
// MeshTools functions
unsigned int MeshTools::total_weight(const MeshBase& mesh)
{
  unsigned int weight=0;

  MeshBase::const_element_iterator       el  = mesh.elements_begin();
  const MeshBase::const_element_iterator end = mesh.elements_end(); 

  for ( ; el != end; ++el)
    weight += (*el)->n_nodes();
  
  return weight;
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
  // processor bounding box with no arguments
  // computes the global bounding box
  return processor_bounding_box(mesh);
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
  assert (mesh.n_nodes() != 0);

  Point min(1.e30,   1.e30,  1.e30);
  Point max(-1.e30, -1.e30, -1.e30);

  // By default no processor is specified and we compute
  // the bounding box for the whole domain.
  if (pid == libMesh::invalid_uint)
    {
      for (unsigned int n=0; n<mesh.n_nodes(); n++)
	for (unsigned int i=0; i<mesh.spatial_dimension(); i++)
	  {
	    min(i) = std::min(min(i), mesh.point(n)(i));
	    max(i) = std::max(max(i), mesh.point(n)(i));
	  }      
    }
  // if a specific processor id is specified then we need
  // to only consider those elements living on that processor
  else
    {
      MeshBase::const_element_iterator       el  = mesh.pid_elements_begin(pid);
      const MeshBase::const_element_iterator end = mesh.pid_elements_end(pid);

      for (; el != end; ++el)
	for (unsigned int n=0; n<(*el)->n_nodes(); n++)
	  for (unsigned int i=0; i<mesh.spatial_dimension(); i++)
	    {
	      min(i) = std::min(min(i), mesh.point((*el)->node(n))(i));
	      max(i) = std::max(max(i), mesh.point((*el)->node(n))(i));
	    }      
    }

  const BoundingBox ret_val(min, max);

  return ret_val;  
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

  // By default no subdomain is specified and we compute
  // the bounding box for the whole domain.
  if (sid == libMesh::invalid_uint)
    {
      for (unsigned int n=0; n<mesh.n_nodes(); n++)
	for (unsigned int i=0; i<mesh.spatial_dimension(); i++)
	  {
	    min(i) = std::min(min(i), mesh.point(n)(i));
	    max(i) = std::max(max(i), mesh.point(n)(i));
	  }      
    }

  // if a specific subdomain id is specified then we need
  // to only consider those elements living on that subdomain
  else
    {
      for (unsigned int e=0; e<mesh.n_elem(); e++)
	if (mesh.elem(e)->subdomain_id() == sid)
	  for (unsigned int n=0; n<mesh.elem(e)->n_nodes(); n++)
	    for (unsigned int i=0; i<mesh.spatial_dimension(); i++)
	      {
		min(i) = std::min(min(i), mesh.point(mesh.elem(e)->node(n))(i));
		max(i) = std::max(max(i), mesh.point(mesh.elem(e)->node(n))(i));
	      }      
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
