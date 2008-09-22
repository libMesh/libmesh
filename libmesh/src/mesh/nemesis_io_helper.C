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


#include "nemesis_io_helper.h"

#if defined(LIBMESH_HAVE_NEMESIS_API) && defined(LIBMESH_HAVE_EXODUS_API)


// Initialize the various integer members to zero.  We can check
// these later to see if they've been properly initialized...
Nemesis_IO_Helper::Nemesis_IO_Helper(bool verbose) :
  ex2helper(verbose),
  nemesis_err_flag(0),
  num_nodes_global(0),
  num_elems_global(0),
  num_elem_blks_global(0),
  num_node_sets_global(0),
  num_side_sets_global(0),
  num_proc(0),
  num_proc_in_file(0),
  ftype('\0'),
  num_internal_nodes(0),
  num_border_nodes(0),
  num_external_nodes(0),
  num_internal_elems(0),
  num_border_elems(0),
  num_node_cmaps(0),
  num_elem_cmaps(0),
  _verbose(verbose)
{
  // Warn about using untested code!
  untested();
}


Nemesis_IO_Helper::~Nemesis_IO_Helper()
{
}



void Nemesis_IO_Helper::verbose (bool set_verbosity)
{
  _verbose = set_verbosity;

  // Also set verbosity in the exodus helper object
  ex2helper.verbose(_verbose);
}



void Nemesis_IO_Helper::get_init_global()
{
  nemesis_err_flag =
    Nemesis::ne_get_init_global(ex2helper.ex_id,
				&num_nodes_global,
				&num_elems_global,
				&num_elem_blks_global,
				&num_node_sets_global,
				&num_side_sets_global);
  ex2helper.check_err(nemesis_err_flag, "Error reading initial global data!");

  if (_verbose)
    {
      std::cout << "[" << libMesh::processor_id() << "] " << "num_nodes_global=" << num_nodes_global << std::endl;
      std::cout << "[" << libMesh::processor_id() << "] " << "num_elems_global=" << num_elems_global << std::endl;
      std::cout << "[" << libMesh::processor_id() << "] " << "num_elem_blks_global=" << num_elem_blks_global << std::endl;
      std::cout << "[" << libMesh::processor_id() << "] " << "num_node_sets_global=" << num_node_sets_global << std::endl;
      std::cout << "[" << libMesh::processor_id() << "] " << "num_side_sets_global=" << num_side_sets_global << std::endl;
    }
}



void Nemesis_IO_Helper::get_ss_param_global()
{
  if (num_side_sets_global > 0)
    {
      global_sideset_ids.resize(num_side_sets_global);
      num_global_side_counts.resize(num_side_sets_global);
      num_global_side_df_counts.resize(num_side_sets_global);
      
      nemesis_err_flag =
	Nemesis::ne_get_ss_param_global(ex2helper.ex_id,
					&global_sideset_ids[0],
					&num_global_side_counts[0],
					&num_global_side_df_counts[0]);
      ex2helper.check_err(nemesis_err_flag, "Error reading global sideset parameters!");

      if (_verbose)
	{
	  std::cout << "[" << libMesh::processor_id() << "] " << "Global Sideset IDs, Side Counts, and DF counts:" << std::endl;
	  for (unsigned int bn=0; bn<global_sideset_ids.size(); ++bn)
	    {
	      std::cout << "  [" << libMesh::processor_id() << "] "
			<< "global_sideset_ids["<<bn<<"]=" << global_sideset_ids[bn]
			<< ", num_global_side_counts["<<bn<<"]=" << num_global_side_counts[bn]
			<< ", num_global_side_df_counts["<<bn<<"]=" << num_global_side_df_counts[bn]
			<< std::endl;
	    }
	}
    }
}




void Nemesis_IO_Helper::get_ns_param_global()
{
  if (num_node_sets_global > 0)
    {
      global_nodeset_ids.resize(num_node_sets_global);
      num_global_node_counts.resize(num_node_sets_global);
      num_global_node_df_counts.resize(num_node_sets_global);
      
      nemesis_err_flag =
	Nemesis::ne_get_ns_param_global(ex2helper.ex_id,
					&global_nodeset_ids[0],
					&num_global_node_counts[0],
					&num_global_node_df_counts[0]);
      ex2helper.check_err(nemesis_err_flag, "Error reading global nodeset parameters!");

      if (_verbose)
	{
	  std::cout << "[" << libMesh::processor_id() << "] " << "Global Nodeset IDs, Node Counts, and DF counts:" << std::endl;
	  for (unsigned int bn=0; bn<global_nodeset_ids.size(); ++bn)
	    {
	      std::cout << "  [" << libMesh::processor_id() << "] "
			<< "global_nodeset_ids["<<bn<<"]=" << global_nodeset_ids[bn]
			<< ", num_global_node_counts["<<bn<<"]=" << num_global_node_counts[bn]
			<< ", num_global_node_df_counts["<<bn<<"]=" << num_global_node_df_counts[bn]
			<< std::endl;
	    }
	}
    }
}



void Nemesis_IO_Helper::get_eb_info_global()
{
  global_elem_blk_ids.resize(num_elem_blks_global);
  global_elem_blk_cnts.resize(num_elem_blks_global);
  
  nemesis_err_flag =
    Nemesis::ne_get_eb_info_global(ex2helper.ex_id,
				   &global_elem_blk_ids[0],
				   &global_elem_blk_cnts[0]);			  
  ex2helper.check_err(nemesis_err_flag, "Error reading global element block info!");

  if (_verbose)
    {
      std::cout << "[" << libMesh::processor_id() << "] " << "Global Element Block IDs and Counts:" << std::endl;
      for (unsigned int bn=0; bn<global_elem_blk_ids.size(); ++bn)
	{
	  std::cout << "  [" << libMesh::processor_id() << "] "
		    << "global_elem_blk_ids["<<bn<<"]=" << global_elem_blk_ids[bn]
		    << ", global_elem_blk_cnts["<<bn<<"]=" << global_elem_blk_cnts[bn]
		    << std::endl;
	}
    }
}



void Nemesis_IO_Helper::get_init_info()
{
  nemesis_err_flag =
    Nemesis::ne_get_init_info(ex2helper.ex_id,
			      &num_proc,
			      &num_proc_in_file,
			      &ftype);
  ex2helper.check_err(nemesis_err_flag, "Error reading initial info!");

  if (_verbose)
    {
      std::cout << "[" << libMesh::processor_id() << "] " << "num_proc=" << num_proc << std::endl;
      std::cout << "[" << libMesh::processor_id() << "] " << "num_proc_in_file=" << num_proc_in_file << std::endl;
      std::cout << "[" << libMesh::processor_id() << "] " << "ftype=" << ftype << std::endl;
    }
}



void Nemesis_IO_Helper::get_loadbal_param()
{
  nemesis_err_flag =
    Nemesis::ne_get_loadbal_param(ex2helper.ex_id,
				  &num_internal_nodes,
				  &num_border_nodes,
				  &num_external_nodes,
				  &num_internal_elems,
				  &num_border_elems,
				  &num_node_cmaps,
				  &num_elem_cmaps,
				  libMesh::processor_id() // The ID of the processor for which info is to be read
				  );
  ex2helper.check_err(nemesis_err_flag, "Error reading load balance parameters!");
	

  if (_verbose)
    {
      std::cout << "[" << libMesh::processor_id() << "] " << "num_internal_nodes=" << num_internal_nodes << std::endl;
      std::cout << "[" << libMesh::processor_id() << "] " << "num_border_nodes=" << num_border_nodes << std::endl;
      std::cout << "[" << libMesh::processor_id() << "] " << "num_external_nodes=" << num_external_nodes << std::endl;
      std::cout << "[" << libMesh::processor_id() << "] " << "num_internal_elems=" << num_internal_elems << std::endl;
      std::cout << "[" << libMesh::processor_id() << "] " << "num_border_elems=" << num_border_elems << std::endl;
      std::cout << "[" << libMesh::processor_id() << "] " << "num_node_cmaps=" << num_node_cmaps << std::endl;
      std::cout << "[" << libMesh::processor_id() << "] " << "num_elem_cmaps=" << num_elem_cmaps << std::endl;
    }
}



void Nemesis_IO_Helper::get_elem_map()
{
  elem_mapi.resize(num_internal_elems);
  elem_mapb.resize(num_border_elems);
  
  nemesis_err_flag =
    Nemesis::ne_get_elem_map(ex2helper.ex_id,
			     &elem_mapi[0],
			     &elem_mapb[0],
			     libMesh::processor_id()
			     );
  ex2helper.check_err(nemesis_err_flag, "Error reading element maps!");


  if (_verbose)
    {
      // These are not contiguous ranges....
      //std::cout << "[" << libMesh::processor_id() << "] " << "first interior elem id=" << elem_mapi[0] << std::endl;
      //std::cout << "[" << libMesh::processor_id() << "] " << "last interior elem id=" << elem_mapi.back() << std::endl;
      //std::cout << "[" << libMesh::processor_id() << "] " << "first boundary elem id=" << elem_mapb[0] << std::endl;
      //std::cout << "[" << libMesh::processor_id() << "] " << "last boundary elem id=" << elem_mapb.back() << std::endl;
      std::cout << "[" << libMesh::processor_id() << "] elem_mapi[i] = ";
      for (unsigned int i=0; i< static_cast<unsigned int>(std::min(10, num_internal_elems-1)); ++i)
	std::cout << elem_mapi[i] << ", ";
      std::cout << "... " << elem_mapi.back() << std::endl;

      std::cout << "[" << libMesh::processor_id() << "] elem_mapb[i] = ";
      for (unsigned int i=0; i< static_cast<unsigned int>(std::min(10, num_internal_elems-1)); ++i)
	std::cout << elem_mapb[i] << ", ";
      std::cout << "... " << elem_mapb.back() << std::endl;
    }
}




void Nemesis_IO_Helper::get_node_map()
{
  node_mapi.resize(num_internal_nodes);
  node_mapb.resize(num_border_nodes);
  node_mape.resize(num_external_nodes);
  
  nemesis_err_flag =
    Nemesis::ne_get_node_map(ex2helper.ex_id,
			     &node_mapi[0],
			     &node_mapb[0],
			     &node_mape[0], 
			     libMesh::processor_id()
			     );
  ex2helper.check_err(nemesis_err_flag, "Error reading node maps!");

  if (_verbose)
    {
      // Remark: The Exodus/Nemesis node numbring is always (?) 1-based!  This means the first interior node id will
      // always be == 1.
      std::cout << "[" << libMesh::processor_id() << "] " << "first interior node id=" << node_mapi[0] << std::endl;
      std::cout << "[" << libMesh::processor_id() << "] " << "last interior node id=" << node_mapi.back() << std::endl;

      std::cout << "[" << libMesh::processor_id() << "] " << "first boundary node id=" << node_mapb[0] << std::endl;
      std::cout << "[" << libMesh::processor_id() << "] " << "last boundary node id=" << node_mapb.back() << std::endl;

      // The number of external nodes is sometimes zero, don't try to access
      // node_mape.back() in this case!
      if (num_external_nodes > 0)
	{
	  std::cout << "[" << libMesh::processor_id() << "] " << "first external node id=" << node_mape[0] << std::endl;
	  std::cout << "[" << libMesh::processor_id() << "] " << "last external node id=" << node_mape.back() << std::endl;
	}
    }
}



void Nemesis_IO_Helper::get_cmap_params()
{
  node_cmap_ids.resize(num_node_cmaps);
  node_cmap_node_cnts.resize(num_node_cmaps);
  elem_cmap_ids.resize(num_elem_cmaps);
  elem_cmap_elem_cnts.resize(num_elem_cmaps);
  
  nemesis_err_flag =
    Nemesis::ne_get_cmap_params(ex2helper.ex_id,
				&node_cmap_ids[0],
				&node_cmap_node_cnts[0],
				&elem_cmap_ids[0],
				&elem_cmap_elem_cnts[0],
				libMesh::processor_id());
  ex2helper.check_err(nemesis_err_flag, "Error reading cmap parameters!");


  if (_verbose)
    {
      std::cout << "[" << libMesh::processor_id() << "] ";
      for (unsigned int i=0; i<node_cmap_ids.size(); ++i)
	std::cout << "node_cmap_ids["<<i<<"]=" << node_cmap_ids[i] << " ";
      std::cout << std::endl;
	
      std::cout << "[" << libMesh::processor_id() << "] ";
      for (unsigned int i=0; i<node_cmap_node_cnts.size(); ++i)
	std::cout << "node_cmap_node_cnts["<<i<<"]=" << node_cmap_node_cnts[i] << " ";
      std::cout << std::endl;

      std::cout << "[" << libMesh::processor_id() << "] ";
      for (unsigned int i=0; i<elem_cmap_ids.size(); ++i)
	std::cout << "elem_cmap_ids["<<i<<"]=" << elem_cmap_ids[i] << " ";
      std::cout << std::endl;

      std::cout << "[" << libMesh::processor_id() << "] ";
      for (unsigned int i=0; i<elem_cmap_elem_cnts.size(); ++i)
	std::cout << "elem_cmap_elem_cnts["<<i<<"]=" << elem_cmap_elem_cnts[i] << " ";
      std::cout << std::endl;
    }
}



void Nemesis_IO_Helper::get_node_cmap()
{
  node_cmap_node_ids.resize(num_node_cmaps);
  node_cmap_proc_ids.resize(num_node_cmaps);
  
  for (unsigned int i=0; i<node_cmap_node_ids.size(); ++i)
    {
      node_cmap_node_ids[i].resize(node_cmap_node_cnts[i]);
      node_cmap_proc_ids[i].resize(node_cmap_node_cnts[i]);

      nemesis_err_flag =
	Nemesis::ne_get_node_cmap(ex2helper.ex_id,
				  node_cmap_ids[i],
				  &node_cmap_node_ids[i][0],
				  &node_cmap_proc_ids[i][0],
				  libMesh::processor_id());
      ex2helper.check_err(nemesis_err_flag, "Error reading node cmap node and processor ids!");

      if (_verbose)
	{
	  std::cout << "[" << libMesh::processor_id() << "] node_cmap_node_ids["<<i<<"]=";
	  for (unsigned int j=0; j<node_cmap_node_ids[i].size(); ++j)
	    std::cout << node_cmap_node_ids[i][j] << " ";
	  std::cout << std::endl;

	  // This is basically a vector, all entries of which are = node_cmap_ids[i]
	  // Not sure if it's always guaranteed to be that or what...
	  std::cout << "[" << libMesh::processor_id() << "] node_cmap_proc_ids["<<i<<"]=";
	  for (unsigned int j=0; j<node_cmap_proc_ids[i].size(); ++j)
	    std::cout << node_cmap_proc_ids[i][j] << " ";
	  std::cout << std::endl;
	}
    }
}



void Nemesis_IO_Helper::get_elem_cmap()
{
  elem_cmap_elem_ids.resize(num_elem_cmaps);
  elem_cmap_side_ids.resize(num_elem_cmaps);
  elem_cmap_proc_ids.resize(num_elem_cmaps);
	  
  for (unsigned int i=0; i<elem_cmap_elem_ids.size(); ++i)
    {
      elem_cmap_elem_ids[i].resize(elem_cmap_elem_cnts[i]);
      elem_cmap_side_ids[i].resize(elem_cmap_elem_cnts[i]);
      elem_cmap_proc_ids[i].resize(elem_cmap_elem_cnts[i]);

      nemesis_err_flag =
	Nemesis::ne_get_elem_cmap(ex2helper.ex_id,
				  elem_cmap_ids[i],
				  &elem_cmap_elem_ids[i][0],
				  &elem_cmap_side_ids[i][0],
				  &elem_cmap_proc_ids[i][0],
				  libMesh::processor_id());
      ex2helper.check_err(nemesis_err_flag, "Error reading elem cmap elem, side, and processor ids!");


      if (_verbose)
	{
	  std::cout << "[" << libMesh::processor_id() << "] elem_cmap_elem_ids["<<i<<"]=";
	  for (unsigned int j=0; j<elem_cmap_elem_ids[i].size(); ++j)
	    std::cout << elem_cmap_elem_ids[i][j] << " ";
	  std::cout << std::endl;

	  // These must be the (local) side IDs (in the ExodusII face numbering scheme)
	  // of the sides shared across processors.
	  std::cout << "[" << libMesh::processor_id() << "] elem_cmap_side_ids["<<i<<"]=";
	  for (unsigned int j=0; j<elem_cmap_side_ids[i].size(); ++j)
	    std::cout << elem_cmap_side_ids[i][j] << " ";
	  std::cout << std::endl;

	  // This is basically a vector, all entries of which are = elem_cmap_ids[i]
	  // Not sure if it's always guaranteed to be that or what...
	  std::cout << "[" << libMesh::processor_id() << "] elem_cmap_proc_ids["<<i<<"]=";
	  for (unsigned int j=0; j<elem_cmap_proc_ids[i].size(); ++j)
	    std::cout << elem_cmap_proc_ids[i][j] << " ";
	  std::cout << std::endl;
	}
    }
}

#endif // #if defined(LIBMESH_HAVE_NEMESIS_API) && defined(LIBMESH_HAVE_EXODUS_API)
