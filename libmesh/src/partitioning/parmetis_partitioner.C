// $Id: parmetis_partitioner.C,v 1.3 2003-07-16 18:42:12 benkirk Exp $

// The Next Great Finite Element Library.
// Copyright (C) 2002  Benjamin S. Kirk, John W. Peterson
  
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
#include "mesh_config.h"
#include "mesh.h"
#include "parmetis_partitioner.h"
#include "mesh_logging.h"

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
void ParmetisPartitioner::partition (const unsigned int n_sbdmns)
{
  assert (n_sbdmns > 0);

  // Check for an easy return
  if (n_sbdmns == 1)
    {
      elem_iterator       elem_it (_mesh.elements_begin());
      const elem_iterator elem_end(_mesh.elements_end());
      
      for ( ; elem_it != elem_end; ++elem_it)
	(*elem_it)->set_subdomain_id() = 
	  (*elem_it)->set_processor_id() =
	  0;
      
      return;
    }

// What to do if the Parmetis library IS NOT present
#ifndef HAVE_PARMETIS

  here();
  std::cerr << "ERROR: The library has been built without"  << std::endl
	    << "Parmetis support.  Using a Metis"           << std::endl
	    << "partitioner instead!"                       << std::endl;

  MetisPartitioner mp(_mesh);

  mp.partition (n_sbdmns);
  
// What to do if the Parmetis library IS present
#else

  START_LOG("partition()", "ParmetisPartitioner");

  // Initialize the data structures required by ParMETIS
  this->initialize (n_sbdmns);

  // build the graph corresponding to the _mesh
  this->build_graph ();
  

  // Partition the graph
  std::vector<int>  local_part(_part);
  MPI_Comm mpi_comm = MPI_COMM_WORLD;
  
  Parmetis::ParMETIS_V3_PartKway(&_vtxdist[0], &_xadj[0], &_adjncy[0], &_vwgt[0], NULL,
				 &_wgtflag, &_numflag, &_ncon, &_nparts, &_tpwgts[0],
				 &_ubvec[0], &_options[0], &_edgecut,
				 &local_part[_first_local_elem],
				 &mpi_comm);

  // Collect the partioning information from all the processors.
  assert (_part.size() == local_part.size());
  MPI_Allreduce (&local_part[0], &_part[0], _part.size(), MPI_INT, MPI_SUM,
		 MPI_COMM_WORLD);
  
  // Assign the returned processor ids
  this->assign_partitioning ();
  

  STOP_LOG("partition()", "ParmetisPartitioner");
  
#endif // #ifndef HAVE_PARMETIS ... else ...
  
}



void ParmetisPartitioner::repartition (const unsigned int n_sbdmns)
{
  assert (n_sbdmns > 0);

  // Check for an easy return
  if (n_sbdmns == 1)
    {
      elem_iterator       elem_it (_mesh.elements_begin());
      const elem_iterator elem_end(_mesh.elements_end());
      
      for ( ; elem_it != elem_end; ++elem_it)
	(*elem_it)->set_subdomain_id() = 
	  (*elem_it)->set_processor_id() =
	  0;
      
      return;
    }

// What to do if the Parmetis library IS NOT present
#ifndef HAVE_PARMETIS

  here();
  std::cerr << "ERROR: The library has been built without"  << std::endl
	    << "Parmetis support.  Using a Metis"           << std::endl
	    << "partitioner instead!"                       << std::endl;

  MetisPartitioner mp(_mesh);

  mp.partition (n_sbdmns);
  
// What to do if the Parmetis library IS present
#else

  error();
  
  START_LOG("repartition()", "ParmetisPartitioner");

  STOP_LOG("repartition()", "ParmetisPartitioner");
  
#endif // #ifndef HAVE_PARMETIS ... else ...
  
}



#ifdef HAVE_PARMETIS

void ParmetisPartitioner::initialize (const unsigned int n_sbdmns)
{
  // Get a constant reference to the _mesh, use that to avoid
  // corrupting the non-const _mesh
  const MeshBase& mesh = _mesh;

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
					   static_cast<unsigned int>(-1));
  _first_local_elem = 0;
  unsigned int el_num = 0;
  unsigned int local_el_num = 0;
  
  for (unsigned int proc_id=0; proc_id<libMesh::n_processors(); proc_id++)
    {
      if (proc_id == libMesh::processor_id()) _first_local_elem = el_num;
      
      const_active_pid_elem_iterator       elem_it (mesh.elements_begin(), proc_id);	
      const const_active_pid_elem_iterator elem_end(mesh.elements_end(),   proc_id);
      
      for (; elem_it != elem_end; ++elem_it)
	{
	  assert ((*elem_it)->id() < _forward_map.size());
	  assert ( _forward_map[(*elem_it)->id()] == static_cast<unsigned int>(-1));
	  
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



void ParmetisPartitioner::build_graph ()
{
  // Get a constant reference to the _mesh, use that to avoid
  // corrupting the non-const _mesh
  const MeshBase& mesh = _mesh;
  
  // build the graph in distributed CSR format.  Note that
  // the edges in the graph will correspond to
  // face neighbors

  // Reserve space in the adjacency array
  const unsigned int n_active_local_elem  = mesh.n_active_local_elem();
  _xadj.reserve (n_active_local_elem + 1);
    
  std::vector<const Elem*> neighbors_offspring;
    
  const_active_local_elem_iterator       elem_it (mesh.elements_begin());
  const const_active_local_elem_iterator elem_end(mesh.elements_end());
  
  
  for (; elem_it != elem_end; ++elem_it)
    {
      const Elem* elem = *elem_it;
      
      assert (elem->id() < _forward_map.size());
      assert (_forward_map[elem->id()] !=
	      static_cast<unsigned int>(-1));
	
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
			  static_cast<unsigned int>(-1));

		  _adjncy.push_back (_forward_map[neighbor->id()]);
		}
	      
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
				  static_cast<unsigned int>(-1));
			
			  _adjncy.push_back (_forward_map[child->id()]);
			}
		    }
		}
	    }
	}
    }
    
  // The end of the adjacency array for this elem
  _xadj.push_back (_adjncy.size());    
}



void ParmetisPartitioner::assign_partitioning ()
{
  // Assign the returned processor ids
  active_elem_iterator       elem_it (_mesh.elements_begin());
  const active_elem_iterator elem_end(_mesh.elements_end());
  
  for (; elem_it != elem_end; ++elem_it)
    {
      Elem* elem = *elem_it;

      assert (elem->id() < _forward_map.size());
      assert (_forward_map[elem->id()] !=
	      static_cast<unsigned int>(-1));
      assert (_forward_map[elem->id()] < _part.size());
      
      elem->set_subdomain_id() =
	elem->set_processor_id() =
	static_cast<short int>(_part[_forward_map[elem->id()]]);
      
    }
}

#endif // #ifdef HAVE_PARMETIS
