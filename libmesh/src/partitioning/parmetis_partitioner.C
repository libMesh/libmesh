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



// C++ Includes   -----------------------------------

// Local Includes -----------------------------------
#include "libmesh_config.h"
#include "mesh_base.h"
#include "parallel.h"
#include "parmetis_partitioner.h"
#include "libmesh_logging.h"
#include "elem.h"

#ifdef HAVE_PARMETIS

// Include the MPI header files, which must be accessible for
// ParMETIS to work properly.
#include "mpi.h"

// Include the ParMETIS header files
namespace Parmetis {
  extern "C" {
#     include "parmetis.h"
  }
}

// What to do if ParMETIS is not available
#else
#  include "metis_partitioner.h"
#endif // #ifdef HAVE_PARMETIS ... else ...



// ------------------------------------------------------------
// ParmetisPartitioner implementation
void ParmetisPartitioner::_do_partition (MeshBase& mesh,
					 const unsigned int n_sbdmns)
{
  assert (n_sbdmns > 0);

  // Check for an easy return
  if (n_sbdmns == 1)
    {
      this->single_partition (mesh);
      return;
    }

  // This function must be run on all processors at once
  parallel_only();

// What to do if the Parmetis library IS NOT present
#ifndef HAVE_PARMETIS

  here();
  std::cerr << "ERROR: The library has been built without"  << std::endl
	    << "Parmetis support.  Using a Metis"           << std::endl
	    << "partitioner instead!"                       << std::endl;

  MetisPartitioner mp;

  mp.partition (mesh, n_sbdmns);
  
// What to do if the Parmetis library IS present
#else

  START_LOG("partition()", "ParmetisPartitioner");

  // Initialize the data structures required by ParMETIS
  this->initialize (mesh, n_sbdmns);

  // build the graph corresponding to the mesh
  this->build_graph (mesh);
  

  // Partition the graph
  MPI_Comm mpi_comm = libMesh::COMM_WORLD;
  
  // Call the ParMETIS k-way partitioning algorithm.
  Parmetis::ParMETIS_V3_PartKway(&_vtxdist[0], &_xadj[0], &_adjncy[0], &_vwgt[0], NULL,
				 &_wgtflag, &_numflag, &_ncon, &_nparts, &_tpwgts[0],
				 &_ubvec[0], &_options[0], &_edgecut,
				 &_part[_first_local_elem],
				 &mpi_comm);

  // Collect the partioning information from all the processors.
  Parallel::sum(_part);
  
  // Assign the returned processor ids
  this->assign_partitioning (mesh);
  

  STOP_LOG ("partition()", "ParmetisPartitioner");
  
#endif // #ifndef HAVE_PARMETIS ... else ...
  
}



void ParmetisPartitioner::_do_repartition (MeshBase& mesh,
					   const unsigned int n_sbdmns)
{
  assert (n_sbdmns > 0);

  // Check for an easy return
  if (n_sbdmns == 1)
    {
      MeshBase::element_iterator elem_it  = mesh.elements_begin();
      const MeshBase::element_iterator elem_end = mesh.elements_end(); 

      for ( ; elem_it != elem_end; ++elem_it)
	(*elem_it)->processor_id() = 0;
      
      return;
    }

  // This function must be run on all processors at once
  parallel_only();

// What to do if the Parmetis library IS NOT present
#ifndef HAVE_PARMETIS

  here();
  std::cerr << "ERROR: The library has been built without"  << std::endl
	    << "Parmetis support.  Using a Metis"           << std::endl
	    << "partitioner instead!"                       << std::endl;

  MetisPartitioner mp;

  mp.partition (mesh, n_sbdmns);
  
// What to do if the Parmetis library IS present
#else

  error();
  
  START_LOG("repartition()", "ParmetisPartitioner");

  // Initialize the data structures required by ParMETIS
  this->initialize (mesh, n_sbdmns);

  // build the graph corresponding to the mesh
  this->build_graph (mesh);
  

  // Partition the graph
  std::vector<int> vsize(_vwgt.size(), 1);
  float itr = 1000.0;
  MPI_Comm mpi_comm = libMesh::COMM_WORLD;

  // Call the ParMETIS adaptive repartitioning method.  This respects the
  // original partitioning when computing the new partitioning so as to
  // minimize the required data redistribution.
  Parmetis::ParMETIS_V3_AdaptiveRepart(&_vtxdist[0], &_xadj[0], &_adjncy[0],
				       &_vwgt[0], &vsize[0], NULL,
				       &_wgtflag, &_numflag, &_ncon,
				       &_nparts, &_tpwgts[0],
				       &_ubvec[0], & itr, &_options[0],
				       &_edgecut, &_part[_first_local_elem],
				       &mpi_comm);

  // Collect the partioning information from all the processors.
  Parallel::sum(_part);
  
  // Assign the returned processor ids
  this->assign_partitioning (mesh);
  

  STOP_LOG ("repartition()", "ParmetisPartitioner");
  
#endif // #ifndef HAVE_PARMETIS ... else ...
  
}



// Only need to compile these methods if ParMETIS is present
#ifdef HAVE_PARMETIS

void ParmetisPartitioner::initialize (const MeshBase& mesh,
				      const unsigned int n_sbdmns)
{
  const unsigned int n_elem               = mesh.n_elem();
  const unsigned int n_active_local_elem  = mesh.n_active_local_elem();
  const unsigned int n_active_elem        = mesh.n_active_elem();
  const unsigned int n_procs              = libMesh::n_processors();
  
  // Set parameters.
  _wgtflag = 2;                          // weights on vertices only
  _ncon    = 1;                          // one weight per vertex
  _numflag = 0;                          // C-style 0-based numbering
  _nparts  = static_cast<int>(n_sbdmns); // number of subdomains to create
  _edgecut = 0;                          // the numbers of edges cut by the
                                         //   partition

  // Initialize data structures for ParMETIS
  _vtxdist.resize (n_procs+1);     std::fill (_vtxdist.begin(), _vtxdist.end(), 0);
  _tpwgts.resize  (_nparts);       std::fill (_tpwgts.begin(),  _tpwgts.end(),  1./_nparts);
  _ubvec.resize   (_ncon);         std::fill (_ubvec.begin(),   _ubvec.end(),   1.);
  _part.resize    (n_active_elem); std::fill (_part.begin(),    _part.end(), 0);
  _options.resize (5);
  _vwgt.resize    (n_active_local_elem);
  

  // Set the options
  _options[0] = 0; // use default options


  // Set up the vtxdist array.  This will be the same on each processor.
  // Consult the Parmetis documentation.
  {
    assert (_vtxdist.size() == libMesh::n_processors()+1);
    assert (_vtxdist[0] == 0);
    
    for (unsigned int proc_id=0; proc_id<libMesh::n_processors(); proc_id++)
      _vtxdist[proc_id+1] = _vtxdist[proc_id] + mesh.n_active_elem_on_proc(proc_id);
    
    assert (_vtxdist[libMesh::n_processors()] == static_cast<int>(n_active_elem));
  }

  
  // Metis will only consider the active elements.
  // We need to map the active element ids into a
  // contiguous range.
  _forward_map.resize (n_elem); std::fill (_forward_map.begin(),
					   _forward_map.end(),
					   libMesh::invalid_uint);
  _first_local_elem = 0;
  unsigned int el_num = 0;
  unsigned int local_el_num = 0;
  
  for (unsigned int proc_id=0; proc_id<libMesh::n_processors(); proc_id++)
    {
      if (proc_id == libMesh::processor_id()) _first_local_elem = el_num;

      MeshBase::const_element_iterator       elem_it  = mesh.active_pid_elements_begin(proc_id);
      const MeshBase::const_element_iterator elem_end = mesh.active_pid_elements_end(proc_id);

      for (; elem_it != elem_end; ++elem_it)
	{
	  assert ((*elem_it)->id() < _forward_map.size());
	  assert ( _forward_map[(*elem_it)->id()] == libMesh::invalid_uint);
	  
	  _forward_map[(*elem_it)->id()] = el_num;
	  el_num++;
	  
	  // maybe there is a better weight?
	  if (proc_id == libMesh::processor_id())
	    _vwgt[local_el_num++] = (*elem_it)->n_nodes();
	}
    }

  assert (el_num       == n_active_elem);
  assert (local_el_num == n_active_local_elem);
}



void ParmetisPartitioner::build_graph (const MeshBase& mesh)
{
  // build the graph in distributed CSR format.  Note that
  // the edges in the graph will correspond to
  // face neighbors

  // Reserve space in the adjacency array
  const unsigned int n_active_local_elem  = mesh.n_active_local_elem();
  _xadj.reserve (n_active_local_elem + 1);
    
  std::vector<const Elem*> neighbors_offspring;
    
  MeshBase::const_element_iterator       elem_it  = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator elem_end = mesh.active_local_elements_end(); 

  
  for (; elem_it != elem_end; ++elem_it)
    {
      const Elem* elem = *elem_it;
      
      assert (elem->id() < _forward_map.size());
      assert (_forward_map[elem->id()] !=
	      libMesh::invalid_uint);
	
      // The beginning of the adjacency array for this elem
      _xadj.push_back (_adjncy.size());

      // Loop over the element's neighbors.  An element
      // adjacency corresponds to a face neighbor
      for (unsigned int ms=0; ms<elem->n_neighbors(); ms++)
	{
	  const Elem* neighbor = elem->neighbor(ms);
	  
	  if (neighbor != NULL)
	    {  
	      // If the neighbor is active treat it
	      // as a connection
	      if (neighbor->active())
		{
		  assert (neighbor->id() < _forward_map.size());
		  assert (_forward_map[neighbor->id()] !=
			  libMesh::invalid_uint);

		  _adjncy.push_back (_forward_map[neighbor->id()]);
		}
  
#ifdef ENABLE_AMR
	      
	      // Otherwise we need to find all of the
	      // neighbor's children that are connected to
	      // us and add them
	      else
		{
		  // The side of the neighbor to which
		  // we are connected
		  const unsigned int ns =
		    neighbor->which_neighbor_am_i (elem);
		  
		  // Get all the active children (& grandchildren, etc...)
		  // of the neighbor.
		  neighbor->active_family_tree (neighbors_offspring);
		  
		  // Get all the neighbor's children that
		  // live on that side and are thus connected
		  // to us
		  for (unsigned int nc=0; nc<neighbors_offspring.size(); nc++)
		    {
		      const Elem* child =
			neighbors_offspring[nc];
		      
		      // This does not assume a level-1 mesh.
		      // Note that since children have sides numbered
		      // coincident with the parent then this is a sufficient test.
		      if (child->neighbor(ns) == elem)
			{
			  assert (child->active());
			  assert (child->id() < _forward_map.size());
			  assert (_forward_map[child->id()] !=
				  libMesh::invalid_uint);
			
			  _adjncy.push_back (_forward_map[child->id()]);
			}
		    }
		}

#endif /* ifdef ENABLE_AMR */


	    }
	}
    }
    
  // The end of the adjacency array for this elem
  _xadj.push_back (_adjncy.size());    
}



void ParmetisPartitioner::assign_partitioning (MeshBase& mesh)
{
  // Assign the returned processor ids
  MeshBase::element_iterator       elem_it  = mesh.active_elements_begin();
  const MeshBase::element_iterator elem_end = mesh.active_elements_end();

  for (; elem_it != elem_end; ++elem_it)
    {
      Elem* elem = *elem_it;

      assert (elem->id() < _forward_map.size());
      assert (_forward_map[elem->id()] != libMesh::invalid_uint);
      assert (_forward_map[elem->id()] < _part.size());
      
      const unsigned int elem_procid =
        static_cast<short int>(_part[_forward_map[elem->id()]]);
      while (elem)
        {
          elem->processor_id() = elem_procid;
          elem = elem->parent();
        }
      
    }
}

#endif // #ifdef HAVE_PARMETIS
