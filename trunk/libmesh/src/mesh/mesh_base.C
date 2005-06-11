// $Id: mesh_base.C,v 1.95 2005-06-11 05:11:31 jwpeterson Exp $

// The libMesh Finite Element Library.
// Copyright (C) 2002-2005  Benjamin S. Kirk, John W. Peterson
  
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
#include <sstream>   // for std::ostringstream
#include <map>       // for std::multimap


// Local includes
#include "mesh_base.h"
#include "libmesh_logging.h"
#include "metis_partitioner.h" // for default partitioning
#include "elem.h"
#include "boundary_info.h"



// ------------------------------------------------------------
// MeshBase class member functions
MeshBase::MeshBase (unsigned int d) :
  boundary_info (new BoundaryInfo(*this)),
  _n_sbd        (1),
  _n_parts      (1),
  _dim          (d),
  _is_prepared  (false)
{
  assert (DIM <= 3);
  assert (DIM >= _dim);
  assert (libMesh::initialized());
}



MeshBase::MeshBase (const MeshBase& other_mesh) :
  boundary_info      (new BoundaryInfo(*this)), // no copy constructor defined for BoundaryInfo?
  _n_sbd             (other_mesh._n_sbd),
  _n_parts           (other_mesh._n_parts),
  _dim               (other_mesh._dim),
  _is_prepared       (other_mesh._is_prepared)

{
}



MeshBase::~MeshBase()
{
  this->clear();

  assert (!libMesh::closed());
}



void MeshBase::prepare_for_use ()
{
  // Renumber the nodes and elements so that they in
  // contiguous blocks.
  this->renumber_nodes_and_elements();
  
  // Let all the elements find their neighbors
  this->find_neighbors();

  // Partition the mesh.
  this->partition();
  
  // The mesh is now prepared for use.
  _is_prepared = true;
}



unsigned int MeshBase::n_active_elem () const
{
    return static_cast<unsigned int>(std::distance (this->active_elements_begin(),
						    this->active_elements_end()));
//   unsigned int num=0;

//   const_element_iterator       el  = this->active_elements_begin();
//   const const_element_iterator end = this->active_elements_end(); 
  
//   for (; el!=end; ++el)
//     num++;

//   return num;
}

 

void MeshBase::clear ()
{
  
  // Reset the number of subdomains
  _n_sbd  = 1;

  // Reset the number of partitions
  _n_parts = 1;

  // Reset the _is_prepared flag
  _is_prepared = false;

  // Clear boundary information
  this->boundary_info->clear();
}



unsigned int MeshBase::n_elem_on_proc (const unsigned int proc_id) const
{
    assert (proc_id < libMesh::n_processors());
    return static_cast<unsigned int>(std::distance (this->pid_elements_begin(proc_id),
						    this->pid_elements_end  (proc_id)));
//   assert (proc_id < libMesh::n_processors());

//   unsigned int ne=0;

//   const_element_iterator       el  = this->pid_elements_begin(proc_id);
//   const const_element_iterator end = this->pid_elements_end(proc_id);

//   for (; el!=end; ++el)
//     ne++;

//   return ne;
}



unsigned int MeshBase::n_active_elem_on_proc (const unsigned int proc_id) const
{
    assert (proc_id < libMesh::n_processors());
    return static_cast<unsigned int>(std::distance (this->active_pid_elements_begin(proc_id),
						    this->active_pid_elements_end  (proc_id)));
//   assert (proc_id < libMesh::n_processors());

//   unsigned int ne=0;

//   const_element_iterator       el  = this->active_pid_elements_begin(proc_id);
//   const const_element_iterator end = this->active_pid_elements_end(proc_id);

  
//   for (; el!=end; ++el)
//     ne++;

//   return ne;
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
      << "  n_elem()="            << this->n_elem()            << '\n'
      << "   n_local_elem()="     << this->n_local_elem()      << '\n'
#ifdef ENABLE_AMR
      << "   n_active_elem()="    << this->n_active_elem()     << '\n'
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
  MetisPartitioner partitioner;
  partitioner.partition (*this, n_parts); 
}









