// $Id$

// The libMesh Finite Element Library.
// Copyright (C) 2002-2008 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner
  
// This library is free software; you can redistribute and/or
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



// library configuration
#include "libmesh_config.h"

// C++ includes
#include <algorithm> // for std::min
#include <map>       // for std::multimap
#include <sstream>   // for std::ostringstream


// Local includes
#include "boundary_info.h"
#include "elem.h"
#include "mesh_base.h"
#include "mesh_tools.h"
#include "parallel.h"
#include "partitioner.h"
#include "point_locator_base.h"

namespace libMesh
{



// ------------------------------------------------------------
// MeshBase class member functions
MeshBase::MeshBase (unsigned int d) :
  boundary_info  (new BoundaryInfo(*this)),
  _n_parts       (1),
  _dim           (d),
  _is_prepared   (false),
  _point_locator (NULL),
  _partitioner   (NULL)
{
  libmesh_assert (LIBMESH_DIM <= 3);
  libmesh_assert (LIBMESH_DIM >= _dim);
  libmesh_assert (libMesh::initialized());
}



MeshBase::MeshBase (const MeshBase& other_mesh) :
  boundary_info  (new BoundaryInfo(*this)), // no copy constructor defined for BoundaryInfo?
  _n_parts       (other_mesh._n_parts),
  _dim           (other_mesh._dim),
  _is_prepared   (other_mesh._is_prepared),
  _point_locator (NULL),
  _partitioner   (other_mesh._partitioner->clone())
{
}



MeshBase::~MeshBase()
{
  this->clear();

  libmesh_assert (!libMesh::closed());
}



void MeshBase::prepare_for_use (const bool skip_renumber_nodes_and_elements)
{  
  // Renumber the nodes and elements so that they in contiguous
  // blocks.  By default, skip_renumber_nodes_and_elements is false,
  // however we may skip this step by passing
  // skip_renumber_nodes_and_elements==true to this function.
  //
  // Instances where you if prepare_for_use() should not renumber the nodes
  // and elements include reading in e.g. an xda/r or gmv file. In
  // this case, the ordering of the nodes may depend on an accompanying
  // solution, and the node ordering cannot be changed.
  if(!skip_renumber_nodes_and_elements)
    this->renumber_nodes_and_elements();
  
  // Let all the elements find their neighbors
  this->find_neighbors();

  // Partition the mesh.
  this->partition();
  
  // If we're using ParallelMesh, we'll want it parallelized.
  this->delete_remote_elements();

  if(!skip_renumber_nodes_and_elements)
    this->renumber_nodes_and_elements();

  // Reset our PointLocator.  This needs to happen any time the elements
  // in the underlying elements in the mesh have changed, so we do it here.
  this->clear_point_locator();

  // The mesh is now prepared for use.
  _is_prepared = true;
}



void MeshBase::clear ()
{
  // Reset the number of partitions
  _n_parts = 1;

  // Reset the _is_prepared flag
  _is_prepared = false;

  // Clear boundary information
  this->boundary_info->clear();

  // Clear our point locator.
  this->clear_point_locator();
}



unsigned int MeshBase::n_subdomains() const
{
  // This requires an inspection on every processor
  parallel_only();

  const_element_iterator       el  = this->active_elements_begin();
  const const_element_iterator end = this->active_elements_end(); 

  std::set<unsigned int> subdomain_ids;
  
  for (; el!=end; ++el)
    subdomain_ids.insert((*el)->subdomain_id());

  // Some subdomains may only live on other processors
  Parallel::set_union(subdomain_ids, 0);

  unsigned int n_sbd_ids = subdomain_ids.size();
  Parallel::broadcast(n_sbd_ids);

  return n_sbd_ids;
}


unsigned int MeshBase::n_nodes_on_proc (const unsigned int proc_id) const
{
  // We're either counting a processor's nodes or unpartitioned
  // nodes
  libmesh_assert (proc_id < libMesh::n_processors() ||
	  proc_id == DofObject::invalid_processor_id);
  
  return static_cast<unsigned int>(std::distance (this->pid_nodes_begin(proc_id),
						  this->pid_nodes_end  (proc_id)));
}



unsigned int MeshBase::n_elem_on_proc (const unsigned int proc_id) const
{
  // We're either counting a processor's elements or unpartitioned
  // elements
  libmesh_assert (proc_id < libMesh::n_processors() ||
	  proc_id == DofObject::invalid_processor_id);
  
  return static_cast<unsigned int>(std::distance (this->pid_elements_begin(proc_id),
						  this->pid_elements_end  (proc_id)));
}



unsigned int MeshBase::n_active_elem_on_proc (const unsigned int proc_id) const
{
  libmesh_assert (proc_id < libMesh::n_processors());
  return static_cast<unsigned int>(std::distance (this->active_pid_elements_begin(proc_id),
						  this->active_pid_elements_end  (proc_id)));
}
  


unsigned int MeshBase::n_sub_elem () const
{
  unsigned int ne=0;

  const_element_iterator       el  = this->elements_begin();
  const const_element_iterator end = this->elements_end(); 

  for (; el!=end; ++el)
    ne += (*el)->n_sub_elem(); 

  return ne;
}



unsigned int MeshBase::n_active_sub_elem () const
{
  unsigned int ne=0;

  const_element_iterator       el  = this->active_elements_begin();
  const const_element_iterator end = this->active_elements_end(); 

  for (; el!=end; ++el)
    ne += (*el)->n_sub_elem(); 

  return ne;
}



std::string MeshBase::get_info() const
{
  std::ostringstream out;

  out << " Mesh Information:"                                  << '\n'
      << "  mesh_dimension()="    << this->mesh_dimension()    << '\n'
      << "  spatial_dimension()=" << this->spatial_dimension() << '\n'
      << "  n_nodes()="           << this->n_nodes()           << '\n'
      << "    n_local_nodes()="   << this->n_local_nodes()     << '\n'
      << "  n_elem()="            << this->n_elem()            << '\n'
      << "    n_local_elem()="    << this->n_local_elem()      << '\n'
#ifdef LIBMESH_ENABLE_AMR
      << "    n_active_elem()="   << this->n_active_elem()     << '\n'
#endif
      << "  n_subdomains()="      << this->n_subdomains()      << '\n'
      << "  n_processors()="      << this->n_processors()      << '\n'
      << "  processor_id()="      << this->processor_id()      << '\n';

  return out.str();
}


void MeshBase::print_info(std::ostream& os) const
{
  os << this->get_info()
     << std::endl;
}


std::ostream& operator << (std::ostream& os, const MeshBase& m)
{
  m.print_info(os);
  return os;
}




void MeshBase::partition (const unsigned int n_parts)
{
  if (partitioner().get()) // "NULL" means don't partition
    partitioner()->partition (*this, n_parts);
}



unsigned int MeshBase::recalculate_n_partitions()
{
  const_element_iterator       el  = this->active_elements_begin();
  const const_element_iterator end = this->active_elements_end(); 

  unsigned int max_proc_id=0;
  
  for (; el!=end; ++el)
    max_proc_id = std::max(max_proc_id, static_cast<unsigned int>((*el)->processor_id()));

  // The number of partitions is one more than the max processor ID.
  _n_parts = max_proc_id+1;

  return _n_parts;
}



const PointLocatorBase & MeshBase::point_locator () const
{
  // Double checked lock pattern for efficiency
  if (_point_locator.get() == NULL)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    
    if (_point_locator.get() == NULL)
      _point_locator.reset (PointLocatorBase::build(TREE, *this).release());
  }

  return *_point_locator;
}



void MeshBase::clear_point_locator ()
{
  _point_locator.reset(NULL);
}

} // namespace libMesh




