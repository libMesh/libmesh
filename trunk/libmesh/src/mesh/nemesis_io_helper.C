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


// C++ headers
#include <set>

// Libmesh headers
#include "nemesis_io_helper.h"
#include "parallel_mesh.h"
#include "node.h"
#include "elem.h"
#include "boundary_info.h"
#include "parallel.h"

#if defined(LIBMESH_HAVE_NEMESIS_API) && defined(LIBMESH_HAVE_EXODUS_API)

namespace libMesh
{


// Initialize the various integer members to zero.  We can check
// these later to see if they've been properly initialized...
Nemesis_IO_Helper::Nemesis_IO_Helper(bool verbose) :
  ExodusII_IO_Helper(verbose),
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
  num_elem_cmaps(0)
{
  // Warn about using untested code!
  libmesh_experimental();
}


Nemesis_IO_Helper::~Nemesis_IO_Helper()
{
  // Our destructor is called from Nemesis_IO.  We close the Exodus file here since we have 
  // responsibility for managing the file's lifetime.
  this->ex_err = exII::ex_update(this->ex_id);
  this->check_err(ex_err, "Error flushing buffers to file.");
  this->close();
}



// void Nemesis_IO_Helper::verbose (bool set_verbosity)
// {
//   _verbose = set_verbosity;

//   // Also set verbosity in the exodus helper object
//   ex2helper.verbose(_verbose);
// }



void Nemesis_IO_Helper::get_init_global()
{
  nemesis_err_flag =
    Nemesis::ne_get_init_global(ex_id,
				&num_nodes_global,
				&num_elems_global,
				&num_elem_blks_global,
				&num_node_sets_global,
				&num_side_sets_global);
  this->check_err(nemesis_err_flag, "Error reading initial global data!");

  if (_verbose)
    {
      libMesh::out << "[" << libMesh::processor_id() << "] " << "num_nodes_global=" << num_nodes_global << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] " << "num_elems_global=" << num_elems_global << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] " << "num_elem_blks_global=" << num_elem_blks_global << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] " << "num_node_sets_global=" << num_node_sets_global << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] " << "num_side_sets_global=" << num_side_sets_global << std::endl;
    }
}



void Nemesis_IO_Helper::get_ss_param_global()
{
  if (num_side_sets_global > 0)
    {
      global_sideset_ids.resize(num_side_sets_global);
      num_global_side_counts.resize(num_side_sets_global);

      // df = "distribution factor", not really sure what that is.  I don't yet have a file
      // which has distribution factors so I guess we'll worry about it later...
      num_global_side_df_counts.resize(num_side_sets_global);
      
      nemesis_err_flag =
	Nemesis::ne_get_ss_param_global(ex_id,
					global_sideset_ids.empty()        ? NULL : &global_sideset_ids[0],
					num_global_side_counts.empty()    ? NULL : &num_global_side_counts[0],
					num_global_side_df_counts.empty() ? NULL : &num_global_side_df_counts[0]);
      this->check_err(nemesis_err_flag, "Error reading global sideset parameters!");

      if (_verbose)
	{
	  libMesh::out << "[" << libMesh::processor_id() << "] " << "Global Sideset IDs, Side Counts, and DF counts:" << std::endl;
	  for (unsigned int bn=0; bn<global_sideset_ids.size(); ++bn)
	    {
	      libMesh::out << "  [" << libMesh::processor_id() << "] "
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
	Nemesis::ne_get_ns_param_global(ex_id,
					global_nodeset_ids.empty()        ? NULL : &global_nodeset_ids[0],
					num_global_node_counts.empty()    ? NULL : &num_global_node_counts[0],
					num_global_node_df_counts.empty() ? NULL : &num_global_node_df_counts[0]);
      this->check_err(nemesis_err_flag, "Error reading global nodeset parameters!");

      if (_verbose)
	{
	  libMesh::out << "[" << libMesh::processor_id() << "] " << "Global Nodeset IDs, Node Counts, and DF counts:" << std::endl;
	  for (unsigned int bn=0; bn<global_nodeset_ids.size(); ++bn)
	    {
	      libMesh::out << "  [" << libMesh::processor_id() << "] "
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
    Nemesis::ne_get_eb_info_global(ex_id,
				   global_elem_blk_ids.empty()  ? NULL : &global_elem_blk_ids[0],
				   global_elem_blk_cnts.empty() ? NULL : &global_elem_blk_cnts[0]);			  
  this->check_err(nemesis_err_flag, "Error reading global element block info!");

  if (_verbose)
    {
      libMesh::out << "[" << libMesh::processor_id() << "] " << "Global Element Block IDs and Counts:" << std::endl;
      for (unsigned int bn=0; bn<global_elem_blk_ids.size(); ++bn)
	{
	  libMesh::out << "  [" << libMesh::processor_id() << "] "
		        << "global_elem_blk_ids["<<bn<<"]=" << global_elem_blk_ids[bn]
		        << ", global_elem_blk_cnts["<<bn<<"]=" << global_elem_blk_cnts[bn]
		        << std::endl;
	}
    }
}



void Nemesis_IO_Helper::get_init_info()
{
  nemesis_err_flag =
    Nemesis::ne_get_init_info(ex_id,
			      &num_proc,
			      &num_proc_in_file,
			      &ftype);
  this->check_err(nemesis_err_flag, "Error reading initial info!");

  if (_verbose)
    {
      libMesh::out << "[" << libMesh::processor_id() << "] " << "num_proc=" << num_proc << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] " << "num_proc_in_file=" << num_proc_in_file << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] " << "ftype=" << ftype << std::endl;
    }
}



void Nemesis_IO_Helper::get_loadbal_param()
{
  nemesis_err_flag =
    Nemesis::ne_get_loadbal_param(ex_id,
				  &num_internal_nodes,
				  &num_border_nodes,
				  &num_external_nodes,
				  &num_internal_elems,
				  &num_border_elems,
				  &num_node_cmaps,
				  &num_elem_cmaps,
				  libMesh::processor_id() // The ID of the processor for which info is to be read
				  );
  this->check_err(nemesis_err_flag, "Error reading load balance parameters!");
	

  if (_verbose)
    {
      libMesh::out << "[" << libMesh::processor_id() << "] " << "num_internal_nodes=" << num_internal_nodes << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] " << "num_border_nodes=" << num_border_nodes << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] " << "num_external_nodes=" << num_external_nodes << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] " << "num_internal_elems=" << num_internal_elems << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] " << "num_border_elems=" << num_border_elems << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] " << "num_node_cmaps=" << num_node_cmaps << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] " << "num_elem_cmaps=" << num_elem_cmaps << std::endl;
    }
}



void Nemesis_IO_Helper::get_elem_map()
{
  elem_mapi.resize(num_internal_elems);
  elem_mapb.resize(num_border_elems);
  
  nemesis_err_flag =
    Nemesis::ne_get_elem_map(ex_id,
			     elem_mapi.empty() ? NULL : &elem_mapi[0],
			     elem_mapb.empty() ? NULL : &elem_mapb[0],
			     libMesh::processor_id()
			     );
  this->check_err(nemesis_err_flag, "Error reading element maps!");


  if (_verbose)
    {
      // These are not contiguous ranges....
      //libMesh::out << "[" << libMesh::processor_id() << "] " << "first interior elem id=" << elem_mapi[0] << std::endl;
      //libMesh::out << "[" << libMesh::processor_id() << "] " << "last interior elem id=" << elem_mapi.back() << std::endl;
      //libMesh::out << "[" << libMesh::processor_id() << "] " << "first boundary elem id=" << elem_mapb[0] << std::endl;
      //libMesh::out << "[" << libMesh::processor_id() << "] " << "last boundary elem id=" << elem_mapb.back() << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] elem_mapi[i] = ";
      for (unsigned int i=0; i< static_cast<unsigned int>(num_internal_elems-1); ++i)
	libMesh::out << elem_mapi[i] << ", ";
      libMesh::out << "... " << elem_mapi.back() << std::endl;

      libMesh::out << "[" << libMesh::processor_id() << "] elem_mapb[i] = ";
      for (unsigned int i=0; i< static_cast<unsigned int>(std::min(10, num_border_elems-1)); ++i)
	libMesh::out << elem_mapb[i] << ", ";
      libMesh::out << "... " << elem_mapb.back() << std::endl;
    }
}




void Nemesis_IO_Helper::get_node_map()
{
  node_mapi.resize(num_internal_nodes);
  node_mapb.resize(num_border_nodes);
  node_mape.resize(num_external_nodes);
  
  nemesis_err_flag =
    Nemesis::ne_get_node_map(ex_id,
			     node_mapi.empty() ? NULL : &node_mapi[0],
			     node_mapb.empty() ? NULL : &node_mapb[0],
			     node_mape.empty() ? NULL : &node_mape[0], 
			     libMesh::processor_id()
			     );
  this->check_err(nemesis_err_flag, "Error reading node maps!");

  if (_verbose)
    {
      // Remark: The Exodus/Nemesis node numbring is always (?) 1-based!  This means the first interior node id will
      // always be == 1.
      libMesh::out << "[" << libMesh::processor_id() << "] " << "first interior node id=" << node_mapi[0] << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] " << "last interior node id=" << node_mapi.back() << std::endl;

      libMesh::out << "[" << libMesh::processor_id() << "] " << "first boundary node id=" << node_mapb[0] << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] " << "last boundary node id=" << node_mapb.back() << std::endl;

      // The number of external nodes is sometimes zero, don't try to access
      // node_mape.back() in this case!
      if (num_external_nodes > 0)
	{
	  libMesh::out << "[" << libMesh::processor_id() << "] " << "first external node id=" << node_mape[0] << std::endl;
	  libMesh::out << "[" << libMesh::processor_id() << "] " << "last external node id=" << node_mape.back() << std::endl;
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
    Nemesis::ne_get_cmap_params(ex_id,
				node_cmap_ids.empty()       ? NULL : &node_cmap_ids[0],
				node_cmap_node_cnts.empty() ? NULL : &node_cmap_node_cnts[0],
				elem_cmap_ids.empty()       ? NULL : &elem_cmap_ids[0],
				elem_cmap_elem_cnts.empty() ? NULL : &elem_cmap_elem_cnts[0],
				libMesh::processor_id());
  this->check_err(nemesis_err_flag, "Error reading cmap parameters!");


  if (_verbose)
    {
      libMesh::out << "[" << libMesh::processor_id() << "] ";
      for (unsigned int i=0; i<node_cmap_ids.size(); ++i)
	libMesh::out << "node_cmap_ids["<<i<<"]=" << node_cmap_ids[i] << " ";
      libMesh::out << std::endl;
	
      libMesh::out << "[" << libMesh::processor_id() << "] ";
      for (unsigned int i=0; i<node_cmap_node_cnts.size(); ++i)
	libMesh::out << "node_cmap_node_cnts["<<i<<"]=" << node_cmap_node_cnts[i] << " ";
      libMesh::out << std::endl;

      libMesh::out << "[" << libMesh::processor_id() << "] ";
      for (unsigned int i=0; i<elem_cmap_ids.size(); ++i)
	libMesh::out << "elem_cmap_ids["<<i<<"]=" << elem_cmap_ids[i] << " ";
      libMesh::out << std::endl;

      libMesh::out << "[" << libMesh::processor_id() << "] ";
      for (unsigned int i=0; i<elem_cmap_elem_cnts.size(); ++i)
	libMesh::out << "elem_cmap_elem_cnts["<<i<<"]=" << elem_cmap_elem_cnts[i] << " ";
      libMesh::out << std::endl;
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
	Nemesis::ne_get_node_cmap(ex_id,
				  node_cmap_ids.empty()         ? 0    : node_cmap_ids[i],
				  node_cmap_node_ids[i].empty() ? NULL : &node_cmap_node_ids[i][0],
				  node_cmap_proc_ids[i].empty() ? NULL : &node_cmap_proc_ids[i][0],
				  libMesh::processor_id());
      this->check_err(nemesis_err_flag, "Error reading node cmap node and processor ids!");

      if (_verbose)
	{
	  libMesh::out << "[" << libMesh::processor_id() << "] node_cmap_node_ids["<<i<<"]=";
	  for (unsigned int j=0; j<node_cmap_node_ids[i].size(); ++j)
	    libMesh::out << node_cmap_node_ids[i][j] << " ";
	  libMesh::out << std::endl;

	  // This is basically a vector, all entries of which are = node_cmap_ids[i]
	  // Not sure if it's always guaranteed to be that or what...
	  libMesh::out << "[" << libMesh::processor_id() << "] node_cmap_proc_ids["<<i<<"]=";
	  for (unsigned int j=0; j<node_cmap_proc_ids[i].size(); ++j)
	    libMesh::out << node_cmap_proc_ids[i][j] << " ";
	  libMesh::out << std::endl;
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
	Nemesis::ne_get_elem_cmap(ex_id,
				  elem_cmap_ids.empty()         ? 0    : elem_cmap_ids[i],
				  elem_cmap_elem_ids[i].empty() ? NULL : &elem_cmap_elem_ids[i][0],
				  elem_cmap_side_ids[i].empty() ? NULL : &elem_cmap_side_ids[i][0],
				  elem_cmap_proc_ids[i].empty() ? NULL : &elem_cmap_proc_ids[i][0],
				  libMesh::processor_id());
      this->check_err(nemesis_err_flag, "Error reading elem cmap elem, side, and processor ids!");


      if (_verbose)
	{
	  libMesh::out << "[" << libMesh::processor_id() << "] elem_cmap_elem_ids["<<i<<"]=";
	  for (unsigned int j=0; j<elem_cmap_elem_ids[i].size(); ++j)
	    libMesh::out << elem_cmap_elem_ids[i][j] << " ";
	  libMesh::out << std::endl;

	  // These must be the (local) side IDs (in the ExodusII face numbering scheme)
	  // of the sides shared across processors.
	  libMesh::out << "[" << libMesh::processor_id() << "] elem_cmap_side_ids["<<i<<"]=";
	  for (unsigned int j=0; j<elem_cmap_side_ids[i].size(); ++j)
	    libMesh::out << elem_cmap_side_ids[i][j] << " ";
	  libMesh::out << std::endl;

	  // This is basically a vector, all entries of which are = elem_cmap_ids[i]
	  // Not sure if it's always guaranteed to be that or what...
	  libMesh::out << "[" << libMesh::processor_id() << "] elem_cmap_proc_ids["<<i<<"]=";
	  for (unsigned int j=0; j<elem_cmap_proc_ids[i].size(); ++j)
	    libMesh::out << elem_cmap_proc_ids[i][j] << " ";
	  libMesh::out << std::endl;
	}
    }
}




void Nemesis_IO_Helper::put_init_info(unsigned num_proc,
				      unsigned num_proc_in_file,
				      const char* ftype)
{
  nemesis_err_flag =
    Nemesis::ne_put_init_info(ex_id,
			      num_proc,
			      num_proc_in_file,
			      const_cast<char*>(ftype));

  this->check_err(nemesis_err_flag, "Error writing initial information!");
}




void Nemesis_IO_Helper::put_init_global(unsigned num_nodes_global,
					unsigned num_elems_global,
					unsigned num_elem_blks_global,
					unsigned num_node_sets_global,
					unsigned num_side_sets_global)
{
  nemesis_err_flag = 
    Nemesis::ne_put_init_global(ex_id,
				num_nodes_global,
				num_elems_global,
				num_elem_blks_global,
				num_node_sets_global,
				num_side_sets_global);

  this->check_err(nemesis_err_flag, "Error writing initial global data!");
}



void Nemesis_IO_Helper::put_eb_info_global(std::vector<int>& global_elem_blk_ids,
					   std::vector<int>& global_elem_blk_cnts)
{
  nemesis_err_flag = 
    Nemesis::ne_put_eb_info_global(ex_id,
				   &global_elem_blk_ids[0],
				   &global_elem_blk_cnts[0]);

  this->check_err(nemesis_err_flag, "Error writing global element block information!");
}




void Nemesis_IO_Helper::put_ns_param_global(std::vector<int>& global_nodeset_ids,
					    std::vector<int>& num_global_node_counts,
					    std::vector<int>& num_global_node_df_counts)
{
  // Only add nodesets if there are some
  if(global_nodeset_ids.size())
  {
    nemesis_err_flag = 
      Nemesis::ne_put_ns_param_global(ex_id,
                                      &global_nodeset_ids[0],
                                      &num_global_node_counts[0],
                                      &num_global_node_df_counts[0]);
  }

  this->check_err(nemesis_err_flag, "Error writing global nodeset parameters!");
}




void Nemesis_IO_Helper::put_ss_param_global(std::vector<int>& global_sideset_ids,
					    std::vector<int>& num_global_side_counts,
					    std::vector<int>& num_global_side_df_counts)
{
  nemesis_err_flag = 
    Nemesis::ne_put_ss_param_global(ex_id,
				    &global_sideset_ids[0],
				    &num_global_side_counts[0],
				    &num_global_side_df_counts[0]);

  this->check_err(nemesis_err_flag, "Error writing global sideset parameters!");
}




void Nemesis_IO_Helper::put_loadbal_param(unsigned num_internal_nodes,
					  unsigned num_border_nodes,
					  unsigned num_external_nodes,
					  unsigned num_internal_elems,
					  unsigned num_border_elems,
					  unsigned num_node_cmaps,
					  unsigned num_elem_cmaps)
{
  nemesis_err_flag = 
    Nemesis::ne_put_loadbal_param(ex_id,
				  num_internal_nodes,
				  num_border_nodes,
				  num_external_nodes,
				  num_internal_elems,
				  num_border_elems,
				  num_node_cmaps,
				  num_elem_cmaps,
				  libMesh::processor_id());

  this->check_err(nemesis_err_flag, "Error writing loadbal parameters!");
}





void Nemesis_IO_Helper::put_cmap_params(std::vector<int>& node_cmap_ids,
					std::vector<int>& node_cmap_node_cnts,
					std::vector<int>& elem_cmap_ids,
					std::vector<int>& elem_cmap_elem_cnts)
{
  nemesis_err_flag = 
    Nemesis::ne_put_cmap_params(ex_id,
				&node_cmap_ids[0],
				&node_cmap_node_cnts[0],
				&elem_cmap_ids[0],
				&elem_cmap_elem_cnts[0],
				libMesh::processor_id());

  this->check_err(nemesis_err_flag, "Error writing cmap parameters!");
}




void Nemesis_IO_Helper::put_node_cmap(std::vector<std::vector<int> >& node_cmap_node_ids,
				      std::vector<std::vector<int> >& node_cmap_proc_ids)
{

  // Print to screen what we are about to print to Nemesis file
  if (_verbose)
    {
      for (unsigned i=0; i<node_cmap_node_ids.size(); ++i)
	{
	  libMesh::out << "[" << libMesh::processor_id() << "] put_node_cmap() : nodes communicated to proc " 
		       << this->node_cmap_ids[i]
		       << " = ";
	  for (unsigned j=0; j<node_cmap_node_ids[i].size(); ++j)
	    libMesh::out << node_cmap_node_ids[i][j] << " ";
	  libMesh::out << std::endl;
	}

      for (unsigned i=0; i<node_cmap_node_ids.size(); ++i)
	{
	  libMesh::out << "[" << libMesh::processor_id() << "] put_node_cmap() : processor IDs = ";
	  for (unsigned j=0; j<node_cmap_proc_ids[i].size(); ++j)
	    libMesh::out << node_cmap_proc_ids[i][j] << " ";
	  libMesh::out << std::endl;
	}
    }

  for (unsigned int i=0; i<node_cmap_node_ids.size(); ++i)
    {
      nemesis_err_flag = 
	Nemesis::ne_put_node_cmap(ex_id,
				  this->node_cmap_ids[i],
				  &node_cmap_node_ids[i][0],
				  &node_cmap_proc_ids[i][0],
				  libMesh::processor_id());
      
      this->check_err(nemesis_err_flag, "Error writing node communication map to file!");
    }
}




void Nemesis_IO_Helper::put_node_map(std::vector<int>& node_mapi,
				     std::vector<int>& node_mapb,
				     std::vector<int>& node_mape)
{
  nemesis_err_flag = 
    Nemesis::ne_put_node_map(ex_id,
			     &node_mapi[0],
			     &node_mapb[0],
			     node_mape.empty() ? NULL : &node_mape[0], // Don't take address of empty vector...
			     libMesh::processor_id());

  this->check_err(nemesis_err_flag, "Error writing Nemesis internal and border node maps to file!");
}




void Nemesis_IO_Helper::put_elem_cmap(std::vector<std::vector<int> >& elem_cmap_elem_ids,
				      std::vector<std::vector<int> >& elem_cmap_side_ids,
				      std::vector<std::vector<int> >& elem_cmap_proc_ids)
{
  for (unsigned int i=0; i<elem_cmap_ids.size(); ++i)
    {
      nemesis_err_flag = 
	Nemesis::ne_put_elem_cmap(ex_id,
				  this->elem_cmap_ids[i],
				  &elem_cmap_elem_ids[i][0],
				  &elem_cmap_side_ids[i][0],
				  &elem_cmap_proc_ids[i][0],
				  libMesh::processor_id());
      
      this->check_err(nemesis_err_flag, "Error writing elem communication map to file!");
    }  
}




void Nemesis_IO_Helper::put_elem_map(std::vector<int>& elem_mapi,
				     std::vector<int>& elem_mapb)
{
  nemesis_err_flag =
    Nemesis::ne_put_elem_map(ex_id,
			     &elem_mapi[0],
			     &elem_mapb[0],
			     libMesh::processor_id());
  
  this->check_err(nemesis_err_flag, "Error writing Nemesis internal and border element maps to file!");
}






void Nemesis_IO_Helper::put_n_coord(unsigned start_node_num,
				    unsigned num_nodes,
				    std::vector<Real>& x_coor,
				    std::vector<Real>& y_coor,
				    std::vector<Real>& z_coor)
{
  nemesis_err_flag = 
    Nemesis::ne_put_n_coord(ex_id,
			    start_node_num,
			    num_nodes,
			    &x_coor[0],
			    &y_coor[0],
			    &z_coor[0]);

  this->check_err(nemesis_err_flag, "Error writing coords to file!");
}








// Note: we can't reuse the ExodusII_IO_Helper code directly, since it only runs
// on processor 0.  TODO: We could have the body of this function as a separate 
// function and then ExodusII_IO_Helper would only call it if on processor 0...
void Nemesis_IO_Helper::create(std::string filename)
{
  // Fall back on double precision when necessary since ExodusII
  // doesn't seem to support long double
  comp_ws = std::min(sizeof(Real),sizeof(double));
  io_ws = std::min(sizeof(Real),sizeof(double));
    
  this->ex_id = exII::ex_create(filename.c_str(), EX_CLOBBER, &comp_ws, &io_ws);

  check_err(ex_id, "Error creating Nemesis mesh file.");

  if (_verbose)
    libMesh::out << "File created successfully." << std::endl;

  this->_created = true;
}




void Nemesis_IO_Helper::initialize(std::string title, const MeshBase & mesh)
{
  // Make sure that the reference passed in is really a ParallelMesh
  const ParallelMesh& pmesh = libmesh_cast_ref<const ParallelMesh&>(mesh);

  // According to Nemesis documentation, first call when writing should be to
  // ne_put_init_info().  Our reader doesn't actually call this, but we should
  // strive to be as close to a normal nemesis file as possible...
  this->put_init_info(libMesh::n_processors(), 1, "p");

 
  // Ready to put global initial information into the file.  This consists of
  // three parts labelled I, II, and III below...

  //
  // I.) Need to compute the number of global element blocks.  To be consistent with
  // Exodus, we also incorrectly associate the number of element blocks with the
  // number of libmesh subdomains...
  //
  {
    // 1.) Loop over active local elements, build up set of subdomain IDs.  FIXME: There 
    // should probably be one and only loop over active elements if we can help it.
    std::set<subdomain_id_type> global_subdomain_ids;

    // This map keeps track of the number of elements in each subdomain
    std::map<subdomain_id_type, unsigned> global_subdomain_counts;

    ParallelMesh::const_element_iterator elem_it = pmesh.active_local_elements_begin();
    ParallelMesh::const_element_iterator elem_end = pmesh.active_local_elements_end();
    
    for (; elem_it != elem_end; ++elem_it)
      {
	const Elem* elem = *elem_it;
	
	global_subdomain_ids.insert(elem->subdomain_id());

	// Increment the count of elements in this subdomain
	global_subdomain_counts[elem->subdomain_id()]++;
      }

    {
      // 2.) Copy local subdomain IDs into a vector for communication
      std::vector<subdomain_id_type> global_subdomain_ids_vector(global_subdomain_ids.begin(),
								 global_subdomain_ids.end());

      // 3.) Gather them into an enlarged vector
      Parallel::allgather(global_subdomain_ids_vector);

      // 4.) Insert any new IDs into the set (any duplicates will be dropped)
      global_subdomain_ids.insert(global_subdomain_ids_vector.begin(),
				  global_subdomain_ids_vector.end());
    }
  
    // 5.) Now global_subdomain_ids actually contains a global list of all subdomain IDs
    this->num_elem_blks_global = global_subdomain_ids.size();

    // Print the number of elements found locally in each subdomain
    if (_verbose)
      {
	libMesh::out << "[" << libMesh::processor_id() << "] ";
	for (std::map<subdomain_id_type, unsigned>::iterator it=global_subdomain_counts.begin();
	     it != global_subdomain_counts.end();
	     ++it)
	  {
	    libMesh::out << "ID: " 
			 << static_cast<unsigned>((*it).first) 
			 << ", Count: " << (*it).second << ", ";
	  }
	libMesh::out << std::endl;
      }

    // 6.) Parallel::sum up the number of elements in each block.  We know the global
    // subdomain IDs, so pack them into a vector one by one.  Use a vector of int since
    // that is what Nemesis wants
    this->global_elem_blk_ids.resize(global_subdomain_ids.size());

    unsigned cnt=0;
    for (std::set<subdomain_id_type>::iterator it=global_subdomain_ids.begin();
	 it != global_subdomain_ids.end(); ++it)
      {
	// Find the entry in the local map, note: if not found, will be created with 0 default value, which is OK...
	this->global_elem_blk_ids[cnt++] = global_subdomain_counts[*it];
      }

    // Sum up
    Parallel::sum(this->global_elem_blk_ids);

    if (_verbose)
      {
	libMesh::out << "[" << libMesh::processor_id() << "] global_elem_blk_ids = ";
	for (unsigned i=0; i<this->global_elem_blk_ids.size(); ++i)
	  libMesh::out << this->global_elem_blk_ids[i] << ", ";
	libMesh::out << std::endl;
      }

    // 7.) Create a vector<int> from the global_subdomain_ids set, for passing to Nemesis
    this->global_elem_blk_cnts.clear();
    this->global_elem_blk_cnts.insert(this->global_elem_blk_cnts.end(), // pos
				      global_subdomain_ids.begin(), 
				      global_subdomain_ids.end());


    // 8.) We will call put_eb_info_global later, it must be called after this->put_init_global().
  }



  //
  // II.) Determine the global number of nodesets by communication.
  // This code relies on BoundaryInfo storing side and node
  // boundary IDs separately at the time they are added to the
  // BoundaryInfo object.
  //
  {
    // 1.) Get reference to the set of node boundary IDs
    std::set<short> global_node_boundary_ids(mesh.boundary_info->get_node_boundary_ids().begin(),
					     mesh.boundary_info->get_node_boundary_ids().end());

    // 2.) Copy local IDs into a vector for communication until Roy fixes Parallel::set_union()
    std::vector<short> global_node_boundary_ids_vector(global_node_boundary_ids.begin(),
						       global_node_boundary_ids.end());
    
    // 3.) Gather them into an enlarged vector
    Parallel::allgather(global_node_boundary_ids_vector);

    // 4.) Insert any new IDs into the set (any duplicates will be dropped)
    global_node_boundary_ids.insert(global_node_boundary_ids_vector.begin(),
				    global_node_boundary_ids_vector.end());    

    // 5.) Now global_node_boundary_ids actually contains a global list of all node boundary IDs
    this->num_node_sets_global = global_node_boundary_ids.size();

    // 6.) Create a vector<int> from the global_node_boundary_ids set
    this->global_nodeset_ids.clear();
    this->global_nodeset_ids.insert(this->global_nodeset_ids.end(),
				    global_node_boundary_ids.begin(),
				    global_node_boundary_ids.end());

    if (_verbose)
      {
	libMesh::out << "[" << libMesh::processor_id() << "] global_nodeset_ids = ";
	for (unsigned i=0; i<global_nodeset_ids.size(); ++i) 
	  libMesh::out << global_nodeset_ids[i] << ", ";
	libMesh::out << std::endl;
      }

    // 7.) We also need to know the number of nodes which is in each of the nodesets, globally.
    // There is probably a better way to do this...
    std::vector<unsigned> node_list;
    std::vector<short> node_boundary_id_list;
    mesh.boundary_info->build_node_list(node_list, node_boundary_id_list);

    // Make sure we don't have any left over information
    this->num_global_node_counts.clear();
    this->num_global_node_counts.resize(this->global_nodeset_ids.size());

    // Unfortunately, we can't just count up all occurrences of a given id,
    // that would give us duplicate entries when we do the parallel summation.
    // So instead, only count entries for nodes owned by this processor.
    // Start by getting rid of all non-local node entries from the vectors.
    std::vector<unsigned>::iterator it_node=node_list.begin();
    std::vector<short>::iterator it_id=node_boundary_id_list.begin();
    for ( ; it_node != node_list.end(); )
      {
	if (mesh.node_ptr( *it_node )->processor_id() != libMesh::processor_id())
	  {
	    // Get rid of this node, but do it efficiently for a vector, by popping
	    // it off the back.
	    std::swap (*it_node, node_list.back() );
	    std::swap (*it_id, node_boundary_id_list.back() );

	    node_list.pop_back();
	    node_boundary_id_list.pop_back();
	  }
	else // node is local, go to next
	  {
	    ++it_node;
	    ++it_id;
	  }
      }

    // Now we can do the local count for each ID...
    for (unsigned i=0; i<global_nodeset_ids.size(); ++i)
      {
	this->num_global_node_counts[i] = std::count(node_boundary_id_list.begin(),
						     node_boundary_id_list.end(),
						     this->global_nodeset_ids[i]);
      }

    // And finally we can sum them up
    Parallel::sum(this->num_global_node_counts);

    if (_verbose)
      {
	libMesh::out << "[" << libMesh::processor_id() << "] num_global_node_counts = ";
	for (unsigned i=0; i<num_global_node_counts.size(); ++i) 
	  libMesh::out << num_global_node_counts[i] << ", ";
	libMesh::out << std::endl;
      }
  }

  

  //
  // III.) Need to compute the global number of sidesets by communication:
  // This code relies on BoundaryInfo storing side and node
  // boundary IDs separately at the time they are added to the
  // BoundaryInfo object.
  //
  {
    // 1.) Get reference to the set of side boundary IDs
    std::set<short> global_side_boundary_ids(mesh.boundary_info->get_side_boundary_ids().begin(),
					     mesh.boundary_info->get_side_boundary_ids().end());

    // 2.) Copy local IDs into a vector for communication until Roy fixes Parallel::set_union()
    std::vector<short> global_side_boundary_ids_vector(global_side_boundary_ids.begin(),
						       global_side_boundary_ids.end());
    
    // 3.) Gather them into an enlarged vector
    Parallel::allgather(global_side_boundary_ids_vector);

    // 4.) Insert any new IDs into the set (any duplicates will be dropped)
    global_side_boundary_ids.insert(global_side_boundary_ids_vector.begin(),
				    global_side_boundary_ids_vector.end());    

    // 5.) Now global_side_boundary_ids actually contains a global list of all side boundary IDs
    this->num_side_sets_global = global_side_boundary_ids.size();

    // 6.) Pack these sidesets into a vector so they can be written by Nemesis
    this->global_sideset_ids.clear(); // Make sure there is no leftover information
    this->global_sideset_ids.insert(this->global_sideset_ids.end(),
				    global_side_boundary_ids.begin(),
				    global_side_boundary_ids.end());

    if (_verbose)
      {
	libMesh::out << "[" << libMesh::processor_id() << "] global_sideset_ids = ";
	for (unsigned i=0; i<this->global_sideset_ids.size(); ++i) 
	  libMesh::out << this->global_sideset_ids[i] << ", ";
	libMesh::out << std::endl;
      }

    // We also need global counts of sides in each of the sidesets.  Again, there may be a 
    // better way to do this...
    std::vector<unsigned int> elem_list;
    std::vector<unsigned short int> side_list;
    std::vector<short int> id_list;
    mesh.boundary_info->build_side_list(elem_list, side_list, id_list);
    
    // Similarly to the nodes, we can't count any sides for elements which aren't local
    std::vector<unsigned>::iterator it_elem=elem_list.begin();
    std::vector<unsigned short>::iterator it_side=side_list.begin();
    std::vector<short>::iterator it_id=id_list.begin();
    for ( ; it_elem != elem_list.end(); )
      {
	if (mesh.elem( *it_elem )->processor_id() != libMesh::processor_id())
	  {
	    // Get rid of this elem, but do it efficiently for a vector, by popping
	    // it off the back.
	    std::swap (*it_elem, elem_list.back() );
	    std::swap (*it_side, side_list.back() );
	    std::swap (*it_id, id_list.back() );

	    elem_list.pop_back();
	    side_list.pop_back();
	    id_list.pop_back();
	  }
	else // elem is local, go to next
	  {
	    ++it_elem;
	    ++it_side;
	    ++it_id;
	  }
      }
    
    this->num_global_side_counts.clear(); // Make sure we don't have any leftover information
    this->num_global_side_counts.resize(this->global_sideset_ids.size());
    
    // Get the count for each global sideset ID
    for (unsigned i=0; i<global_sideset_ids.size(); ++i)
      {
	this->num_global_side_counts[i] = std::count(id_list.begin(),
						     id_list.end(),
						     this->global_sideset_ids[i]);
      }

    if (_verbose)
      {
	libMesh::out << "[" << libMesh::processor_id() << "] num_global_side_counts = ";
	for (unsigned i=0; i<this->num_global_side_counts.size(); ++i) 
	  libMesh::out << this->num_global_side_counts[i] << ", ";
	libMesh::out << std::endl;
      }

    // Finally sum up the result
    Parallel::sum(this->num_global_side_counts);

    if (_verbose)
      {
	libMesh::out << "[" << libMesh::processor_id() << "] num_global_side_counts = ";
	for (unsigned i=0; i<this->num_global_side_counts.size(); ++i) 
	  libMesh::out << this->num_global_side_counts[i] << ", ";
	libMesh::out << std::endl;
      }

  }
  

  // Write the global data to file
  this->put_init_global(pmesh.parallel_n_nodes(),
			pmesh.parallel_n_elem(),
			this->num_elem_blks_global, /* I.   */
			this->num_node_sets_global, /* II.  */
			this->num_side_sets_global  /* III. */
			); 

  // Next, we'll write global element block information to the file.  This was already
  // gathered in step I. above
  this->put_eb_info_global(this->global_elem_blk_cnts,
			   this->global_elem_blk_ids);


  // Next, write global nodeset information to the file.  This was already gathered in
  // step II. above.  
  this->num_global_node_df_counts.clear();
  this->num_global_node_df_counts.resize(this->global_nodeset_ids.size()); // distribution factors all zero...
  this->put_ns_param_global(this->global_nodeset_ids,
			    this->num_global_node_counts,
			    this->num_global_node_df_counts);

  
  // Next, write global sideset information to the file.  This was already gathered in
  // step III. above.  
  this->num_global_side_df_counts.clear();
  this->num_global_side_df_counts.resize(this->global_sideset_ids.size()); // distribution factors all zero...
  this->put_ss_param_global(this->global_sideset_ids,
			    this->num_global_side_counts,
			    this->num_global_side_df_counts);
  

  /**
   * Before we go any further we need to derive consistent node and element numbering schemes for all
   * local elems and nodes connected to local elements.
   *
   * Elements have to be numbered contiguously based on what block number they are in.  Therefore we have
   * to do a bit of work to get the block (ie subdomain) numbers first and store them off as block_ids.
   */
  {
    ParallelMesh::const_element_iterator elem_it = pmesh.active_local_elements_begin();
    ParallelMesh::const_element_iterator elem_end = pmesh.active_local_elements_end();

    // First loop over the elements to figure out which elements are in which subdomain
    for (; elem_it != elem_end; ++elem_it)
    {
      Elem * elem = *elem_it;

      // Grab the nodes while we're here
      for (unsigned int n=0; n<elem->n_nodes(); ++n)
        nodes_attached_to_local_elems.insert( elem->node(n) );

      unsigned int cur_subdomain = elem->subdomain_id();

      if(cur_subdomain == 0)
        cur_subdomain = std::numeric_limits<subdomain_id_type>::max();
     
      subdomain_map[cur_subdomain].push_back(elem->id());
    }

    // Set num_nodes which is used by exodusII_io_helper
    num_nodes = nodes_attached_to_local_elems.size();

    // Now come up with a 1 based numbering for these nodes
    exodus_node_num_to_libmesh.clear(); // Make sure it's empty

    // Set the map for nodes
    for(std::set<int>::iterator it = nodes_attached_to_local_elems.begin();
        it != nodes_attached_to_local_elems.end();
        ++it)
    {
      exodus_node_num_to_libmesh.push_back(*it);
      libmesh_node_num_to_exodus[*it] = exodus_node_num_to_libmesh.size();
    }

    /**
     * Now we're going to loop over the subdomain map and build a few things right
     * now that we'll use later.
     */
    exodus_elem_num_to_libmesh.clear(); // Make sure it's empty
    block_ids.clear(); // Make sure it's empty

    // Now loop over each subdomain and get a unique numbering for the elements
    for(std::map<unsigned int, std::vector<unsigned int>  >::iterator it = subdomain_map.begin();
        it != subdomain_map.end();
        ++it)
    {
      block_ids.push_back((*it).first);

      std::vector<unsigned int> & tmp_vec = (*it).second;

      ExodusII_IO_Helper::ElementMaps em;

      //Use the first element in this block to get representative information.
      //Note that Exodus assumes all elements in a block are of the same type!
      //We are using that same assumption here!
      const ExodusII_IO_Helper::Conversion conv = em.assign_conversion(mesh.elem(tmp_vec[0])->type());
      num_nodes_per_elem = mesh.elem(tmp_vec[0])->n_nodes();    

      std::vector<int> & this_block_connectivity = block_id_to_elem_connectivity[(*it).first];
      
      this_block_connectivity.resize(tmp_vec.size()*num_nodes_per_elem);

      for (unsigned int i=0; i<tmp_vec.size(); i++)
      {
        unsigned int elem_id = tmp_vec[i];

        // Set the number map for elements
        exodus_elem_num_to_libmesh.push_back(elem_id);
        libmesh_elem_num_to_exodus[elem_id] = exodus_elem_num_to_libmesh.size();

        Elem * elem = mesh.elem(elem_id);
          
        for (unsigned int j=0; j < static_cast<unsigned int>(num_nodes_per_elem); j++)
        {  
          const unsigned int connect_index   = (i*num_nodes_per_elem)+j;
          const unsigned int elem_node_index = conv.get_node_map(j);
          this_block_connectivity[connect_index] = libmesh_node_num_to_exodus[elem->node(elem_node_index)];
        }
      }
    }
  }

  
  // Next, we're going to write "load balance" parameters
  
  // First we'll collect IDs of border nodes.

  // The set which will eventually contain the IDs of "border nodes".  These are nodes
  // that lie on the boundary between one or more processors.
  std::set<unsigned> border_node_ids;

  // map from processor ID to set of nodes which elements from this processor "touch",
  // that is, 
  // proc_nodes_touched[p] = (set all node IDs found in elements owned by processor p)
  std::map<unsigned, std::set<unsigned> > proc_nodes_touched;
  typedef std::map<unsigned, std::set<unsigned> >::iterator proc_nodes_touched_iterator;

  // Another map to store sets of intersections with each processor (other than ourself, of course)
  std::map<unsigned, std::set<unsigned> > proc_nodes_touched_intersections;

  // We are going to create a lot of intermediate data structures here, so make sure
  // as many as possible all cleaned up by creating scope!
  {
    // Loop over active (not just active local) elements, make sets of node IDs for each
    // processor which has an element that "touches" a node.
    {
      ParallelMesh::const_element_iterator elem_it = pmesh.active_elements_begin();
      ParallelMesh::const_element_iterator elem_end = pmesh.active_elements_end();

      for (; elem_it != elem_end; ++elem_it)
	{
	  const Elem* elem = *elem_it;

	  // Get reference to the set for this processor.  If it does not exist
	  // it will be created.
	  std::set<unsigned>& set_p = proc_nodes_touched[ elem->processor_id() ];
	
	  // Insert all nodes touched by this element into the set
	  for (unsigned int node=0; node<elem->n_nodes(); ++node)
	    set_p.insert(elem->node(node));
	}
    }

    // The number of node communication maps is the number of other processors
    // with which we share nodes. (I think.) This is just the size of the map we just
    // created, minus 1.
    this->num_node_cmaps = proc_nodes_touched.size() - 1;
    libmesh_assert(static_cast<unsigned>(this->num_node_cmaps) < libMesh::n_processors());

    if (_verbose)
      {
	libMesh::out << "[" << libMesh::processor_id() 
		     << "] proc_nodes_touched contains " 
		     << proc_nodes_touched.size()       
		     << " sets of nodes."
		     << std::endl;

	for (proc_nodes_touched_iterator it = proc_nodes_touched.begin();
	     it != proc_nodes_touched.end();
	     ++it)
	  {
	    libMesh::out << "[" << libMesh::processor_id() 
			 << "] proc_nodes_touched[" << (*it).first << "] has " 
			 << (*it).second.size()
			 << " entries."
			 << std::endl;
	  }
      }
  

    // Loop over all the sets we just created and compute intersections with the
    // this processor's set.  Obviously, don't intersect with ourself.
    for (proc_nodes_touched_iterator it = proc_nodes_touched.begin();
	 it != proc_nodes_touched.end();
	 ++it)
      {
	// Don't compute intersections with ourself
	if ((*it).first == libMesh::processor_id())
	  continue;
      
	// Otherwise, compute intersection with other processor and ourself
	std::set<unsigned>& my_set = proc_nodes_touched[libMesh::processor_id()];
	std::set<unsigned>& other_set = (*it).second;
	std::set<unsigned>& result_set = proc_nodes_touched_intersections[ (*it).first ]; // created if does not exist
      
	std::set_intersection(my_set.begin(), my_set.end(),
			      other_set.begin(), other_set.end(),
			      std::inserter(result_set, result_set.end()));
      }

    if (_verbose)
      {
	for (proc_nodes_touched_iterator it = proc_nodes_touched_intersections.begin();
	     it != proc_nodes_touched_intersections.end();
	     ++it)
	  {
	    libMesh::out << "[" << libMesh::processor_id() 
			 << "] proc_nodes_touched_intersections[" << (*it).first << "] has " 
			 << (*it).second.size()
			 << " entries."
			 << std::endl;
	  }
      }

    // Compute the set_union of all the preceding intersections.  This will be the set of
    // border node IDs for this processor.  Remember border_node_ids was declared waaay up
    // there at the beginning of this scope...
    for (proc_nodes_touched_iterator it = proc_nodes_touched_intersections.begin();
	 it != proc_nodes_touched_intersections.end();
	 ++it)
      {
	std::set<unsigned>& other_set = (*it).second;
	std::set<unsigned> intermediate_result; // Don't think we can insert into one of the sets we're unioning...

	std::set_union(border_node_ids.begin(), border_node_ids.end(),
		       other_set.begin(), other_set.end(),
		       std::inserter(intermediate_result, intermediate_result.end()));
      
	// Swap our intermediate result into the final set
	border_node_ids.swap(intermediate_result);
      }

    if (_verbose)
      {
	libMesh::out << "[" << libMesh::processor_id() 
		     << "] border_node_ids.size()=" << border_node_ids.size() 
		     << std::endl;
      }
  } // end scope for border node ID creation
  
  // Store the number of border node IDs to be written to Nemesis file
  this->num_border_nodes = border_node_ids.size();

  

  // Next we'll collect numbers of internal and border elements, and internal nodes.
  // Note: "A border node does not a border element make...", that is, just because one
  // of an element's nodes has been identified as a border node, the element is not
  // necessarily a border element.  It must have a side on the boundary between processors,
  // i.e. have a face neighbor with a different processor id...
  std::set<unsigned> border_elem_ids;
  std::set<unsigned> internal_elem_ids;
  std::set<unsigned> internal_node_ids;

  // Map between processor ID and set of <element,side> pairs bordering that processor ID.
  // This object is used later, so we define it here outside the current local scope.
  std::map<unsigned, std::set<std::pair<unsigned,unsigned> > > proc_border_elem_sets;
  typedef std::map<unsigned, std::set<std::pair<unsigned,unsigned> > >::iterator proc_border_elem_sets_iterator;
   {
     ParallelMesh::const_element_iterator elem_it = pmesh.active_local_elements_begin();
     ParallelMesh::const_element_iterator elem_end = pmesh.active_local_elements_end();

     // A set of all node IDs encountered as we loop over active local elements.
     // IDs in this set which are not in the border_node_ids set are so-called
     // "internal nodes".
     std::set<unsigned> all_node_ids;

     // Set of all local, active element IDs.  After we have identified border element
     // IDs, the set_difference between this set and the border_elem_ids set will give us
     // the set of internal_elem_ids.
     std::set<unsigned> all_elem_ids;
     
     // A set of processor IDs which elements on this processor have as
     // neighbors.  The size of this set will determine the number of
     // element communication maps in Exodus.
     std::set<unsigned> neighboring_processor_ids;

     // Will be used to create conversion objects capable of mapping libmesh
     // element numberings into Nemesis numberings.
     ExodusII_IO_Helper::ElementMaps element_mapper;
     
     for (; elem_it != elem_end; ++elem_it)
       {
	 const Elem* elem = *elem_it;
	 
	 // Add this Elem's ID to all_elem_ids, later we will take the difference
	 // between this set and the set of border_elem_ids, to get the set of
	 // internal_elem_ids.
	 all_elem_ids.insert(elem->id());
	 
	 // Will be set to true if element is determined to be a border element
	 bool is_border_elem = false;

	 // Construct a conversion object for this Element.  This will help us map
	 // Libmesh numberings into Nemesis numberings for sides.
	 ExodusII_IO_Helper::Conversion conv = element_mapper.assign_conversion(elem->type());

	 // Add all this element's node IDs to the set of all node IDs.
	 // The set of internal_node_ids will be the set difference between
	 // the set of all nodes and the set of border nodes.
	 //
	 // In addition, if any node of a local node is listed in the
	 // border nodes list, then this element goes into the proc_border_elem_sets.
	 // Note that there is not a 1:1 correspondence between 
	 // border_elem_ids and the entries which go into proc_border_elem_sets.
	 // The latter is for communication purposes, ie determining which elements
	 // should be shared between processors.
	 for (unsigned int node=0; node<elem->n_nodes(); ++node)
	   {
	     all_node_ids.insert(elem->node(node));

	     //
	     // Uncomment below to add all elements which touch a border node
	     // to the set of communicated elements.  This seems to give too many
	     // elements for communication and is not needed...
	     //

// not needed	     // If this node has already been identified in the border_node_ids set,
// not needed	     // add this element to the appropriate set(s)? in proc_border_elem_sets
// not needed	     if (border_node_ids.find(elem->node(node)) != border_node_ids.end())
// not needed	       {
// not needed		 if (_verbose)
// not needed		   {
// not needed		     libMesh::out << "[" << libMesh::processor_id() << "] Identified Elem " 
// not needed				  << elem->id()
// not needed				  << " to be communicated because of border node "
// not needed				  << elem->node(node)
// not needed				  << std::endl;
// not needed		     
// not needed		   }
// not needed		 
// not needed		 // We need to communicate this element to any processor (other than ourself)
// not needed		 // which touches this border node.  Loop over the proc_nodes_touched_intersections
// not needed		 // sets, look for the current node there, and, if found, add this element to
// not needed		 // set of communicated elements for that processor.
// not needed		 for (proc_nodes_touched_iterator it = proc_nodes_touched_intersections.begin();
// not needed		      it != proc_nodes_touched_intersections.end();
// not needed		      ++it)
// not needed		   {
// not needed		     std::set<unsigned>& intersecting_nodes = (*it).second;
// not needed
// not needed		     if (intersecting_nodes.find(elem->node(node)) != intersecting_nodes.end())
// not needed		       {
// not needed			 if (_verbose)
// not needed			   {
// not needed			     libMesh::out << "[" << libMesh::processor_id() << "] Adding Elem " 
// not needed					  << elem->id()
// not needed					  << " to set of elements to be communicated to processor "
// not needed					  << (*it).first
// not needed					  << std::endl;
// not needed		     
// not needed			   }
// not needed			 
// not needed			 proc_border_elem_sets[ (*it).first ].insert( elem->id() );
// not needed		       }
// not needed		   }
// not needed
// not needed	       } // end if node found in set of border node IDs
	   } // end loop over element's nodes

	 // Loop over element's neighbors, see if it has a neighbor which is off-processor
	 for (unsigned int n=0; n<elem->n_neighbors(); ++n)
	   {
	     if (elem->neighbor(n) != NULL)
	       {
		 unsigned neighbor_proc_id = elem->neighbor(n)->processor_id();

		 // If my neighbor has a different processor ID, I must be a border element.
		 // Also track the neighboring processor ID if it is are different from our processor ID
		 if (neighbor_proc_id != libMesh::processor_id())
		   {
		     is_border_elem = true;
		     neighboring_processor_ids.insert(neighbor_proc_id);
		     
		     // Convert libmesh side(n) of this element into a side ID for Nemesis
		     unsigned nemesis_side_id = conv.get_inverse_side_map(n);
		     
		     if (_verbose)
		       libMesh::out << "[" << libMesh::processor_id() << "] LibMesh side " 
				    << n 
				    << " mapped to (1-based) Exodus side "
				    << nemesis_side_id
				    << std::endl;

		     // Add this element's ID and the ID of the side which is on the boundary
		     // to the set of border elements for this processor.
		     // Note: if the set does not already exist, this creates it.
		     // FIXME: Map Nemesis sides to libmesh sides?
		     proc_border_elem_sets[ neighbor_proc_id ].insert( std::make_pair(elem->id(), nemesis_side_id) );
		   }
	       }
	     // else 
	     //   // boundary element with NULL neighbor, assign invalid_uint
	     //   neighbor_proc_ids[n] = libMesh::invalid_uint;
	   } // end for loop over neighbors
	 
	 // If we're on a border element, add it to the set
	 if (is_border_elem)
	   border_elem_ids.insert( elem->id() );

       } // end for loop over active local elements
     
     // Take the set_difference between all elements and border elements to get internal
     // element IDs
     std::set_difference(all_elem_ids.begin(), all_elem_ids.end(),
			 border_elem_ids.begin(), border_elem_ids.end(),
			 std::inserter(internal_elem_ids, internal_elem_ids.end()));

     // Take the set_difference between all nodes and border nodes to get internal nodes
     std::set_difference(all_node_ids.begin(), all_node_ids.end(),
			 border_node_ids.begin(), border_node_ids.end(),
			 std::inserter(internal_node_ids, internal_node_ids.end()));

     if (_verbose)
       {
	 libMesh::out << "[" << libMesh::processor_id() << "] neighboring_processor_ids = ";
	 for (std::set<unsigned>::iterator it = neighboring_processor_ids.begin();
	      it != neighboring_processor_ids.end();
	      ++it)
	   {
	     libMesh::out << *it << " ";
	   }
	 libMesh::out << std::endl;
       }

     // The size of the neighboring_processor_ids set should be the number of element communication maps 
     this->num_elem_cmaps = neighboring_processor_ids.size();

     if (_verbose)
       libMesh::out << "[" << libMesh::processor_id() << "] " 
		    << "Number of neighboring processor IDs="
		    << this->num_elem_cmaps
		    << std::endl;

     if (_verbose)
       {
	 // Print out counts of border elements for each processor
	 for (proc_border_elem_sets_iterator it=proc_border_elem_sets.begin();
	      it != proc_border_elem_sets.end(); ++it)
	   {
	     libMesh::out << "[" << libMesh::processor_id() << "] "
			  << "Proc " 
			  << (*it).first << " communicates " 
			  << (*it).second.size() << " elements." << std::endl;
	   }
       }

// not needed     // In general, some elements may be communicated to multiple processors and so 
// not needed     // may appear in multiple sets in proc_border_elem_sets.  To determine the actual
// not needed     // count of elements which need to be "communicated" take the union of these sets.
// not needed     std::set<unsigned> union_of_communicated_elem_ids;
// not needed
// not needed     for (proc_border_elem_sets_iterator it=proc_border_elem_sets.begin();
// not needed	  it != proc_border_elem_sets.end(); ++it)
// not needed       {
// not needed	 std::set<unsigned>& other_set = (*it).second;
// not needed	 std::set<unsigned> intermediate_result;
// not needed	 
// not needed	 std::set_union(union_of_communicated_elem_ids.begin(), union_of_communicated_elem_ids.end(),
// not needed			other_set.begin(), other_set.end(),
// not needed			std::inserter(intermediate_result, intermediate_result.end()));
// not needed      
// not needed	// Swap our intermediate result into the final set
// not needed	union_of_communicated_elem_ids.swap(intermediate_result);
// not needed       }
// not needed
// not needed     if (_verbose)
// not needed       {
// not needed	 libMesh::out << "[" << libMesh::processor_id() << "] "
// not needed		      << "Size of unique set of element IDs requiring communication = " 
// not needed		      << union_of_communicated_elem_ids.size() << std::endl;
// not needed
// not needed	 // Print out actual IDs (1-based) of elements to be communicated
// not needed	 libMesh::out << "[" << libMesh::processor_id() << "] Communicating Elements: ";
// not needed	 for (std::set<unsigned>::iterator it=union_of_communicated_elem_ids.begin();
// not needed	      it != union_of_communicated_elem_ids.end();
// not needed	      ++it)
// not needed	   {
// not needed	     std::cout << (*it) + 1 << ", ";
// not needed	   }
// not needed	 std::cout << std::endl;
// not needed       }
     
     
   } // end scope for active+local element loop

   // Store the number of internal and border elements, and the number of internal nodes,
   // to be written to the Nemesis file.     
   this->num_internal_elems = internal_elem_ids.size();
   this->num_border_elems = border_elem_ids.size();
   this->num_internal_nodes = internal_node_ids.size();

   

  if (_verbose)
    {
      libMesh::out << "[" << libMesh::processor_id() << "] num_internal_nodes=" << this->num_internal_nodes << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] num_border_nodes=" << this->num_border_nodes << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] num_border_elems=" << this->num_border_elems << std::endl;
      libMesh::out << "[" << libMesh::processor_id() << "] num_internal_elems=" << this->num_internal_elems << std::endl;
    }

  
  // Write the loadbal information to the file
  this->put_loadbal_param(this->num_internal_nodes,
			  this->num_border_nodes,
			  this->num_external_nodes,
			  this->num_internal_elems,
			  this->num_border_elems,
			  this->num_node_cmaps,
			  this->num_elem_cmaps);
  


  
  // Next, write the communication map parameters to the file...

  // For the nodes, these are the number of entries in the sets in proc_nodes_touched_intersections
  // map computed above.  Note: this map does not contain self-intersections so we can loop over it
  // directly.
  this->node_cmap_node_cnts.clear(); // Make sure we don't have any leftover information...
  this->node_cmap_ids.clear();       // Make sure we don't have any leftover information...
  this->node_cmap_node_cnts.resize(this->num_node_cmaps);
  this->node_cmap_ids.resize(this->num_node_cmaps);

  {
    unsigned cnt=0; // Index into the vector
    for (proc_nodes_touched_iterator it = proc_nodes_touched_intersections.begin();
	 it != proc_nodes_touched_intersections.end();
	 ++it)
      {
	this->node_cmap_ids[cnt] = (*it).first; // The ID of the proc we communicate with 
	this->node_cmap_node_cnts[cnt] = (*it).second.size(); // The number of nodes we communicate
	cnt++; // increment vector index!
      }
  }

  // Print the packed vectors we just filled
  if (_verbose)
    {
      libMesh::out << "[" << libMesh::processor_id() << "] node_cmap_node_cnts = ";
      for (unsigned i=0; i<node_cmap_node_cnts.size(); ++i)
	libMesh::out << node_cmap_node_cnts[i] << ", ";
      libMesh::out << std::endl;

      libMesh::out << "[" << libMesh::processor_id() << "] node_cmap_ids = ";
      for (unsigned i=0; i<node_cmap_ids.size(); ++i)
	libMesh::out << node_cmap_ids[i] << ", ";
      libMesh::out << std::endl;
    }

  // For the elements, we have not yet computed all this information..
  this->elem_cmap_elem_cnts.clear(); // Make sure we don't have any leftover information...
  this->elem_cmap_ids.clear();       // Make sure we don't have any leftover information...
  this->elem_cmap_elem_cnts.resize(this->num_elem_cmaps);
  this->elem_cmap_ids.resize(this->num_elem_cmaps);

  // Pack the elem_cmap_ids and elem_cmap_elem_cnts vectors
  {
    unsigned cnt=0; // Index into the vectors we're filling
    for (proc_border_elem_sets_iterator it = proc_border_elem_sets.begin();
	 it != proc_border_elem_sets.end();
	 ++it)
      {
	this->elem_cmap_ids[cnt] = (*it).first; // The ID of the proc we communicate with 
	this->elem_cmap_elem_cnts[cnt] = (*it).second.size(); // The number of elems we communicate to/from that proc
	cnt++; // increment vector index!
      }
  }
  
  // Print the packed vectors we just filled
  if (_verbose)
    {
      libMesh::out << "[" << libMesh::processor_id() << "] elem_cmap_elem_cnts = ";
      for (unsigned i=0; i<elem_cmap_elem_cnts.size(); ++i)
	libMesh::out << elem_cmap_elem_cnts[i] << ", ";
      libMesh::out << std::endl;

      libMesh::out << "[" << libMesh::processor_id() << "] elem_cmap_ids = ";
      for (unsigned i=0; i<elem_cmap_ids.size(); ++i)
	libMesh::out << elem_cmap_ids[i] << ", ";
      libMesh::out << std::endl;
    }

  // Write communication map parameters to file.  Seems like this should write the
  // n_comm_ids and e_comm_ids entries in the exodus file, but they are still zero 
  // after calling this function...
  this->put_cmap_params(this->node_cmap_ids,
			this->node_cmap_node_cnts,
			this->elem_cmap_ids,
			this->elem_cmap_elem_cnts);


  // Write the node communication maps to file.  Be sure to use 1-based numbering
  // for Exodus/Nemesis...  The node IDs which are communicated are the ones currently stored in
  // proc_nodes_touched_intersections.
  
  // Make sure there's no left-over information
  this->node_cmap_node_ids.clear();
  this->node_cmap_proc_ids.clear();

  // Allocate enough space for all our node maps
  this->node_cmap_node_ids.resize(this->num_node_cmaps);
  this->node_cmap_proc_ids.resize(this->num_node_cmaps);
  {
    unsigned cnt=0; // Index into vectors
    for (proc_nodes_touched_iterator it = proc_nodes_touched_intersections.begin();
	 it != proc_nodes_touched_intersections.end();
	 ++it)
      {
	// Make sure the current node_cmap_id matches the index in our map of node intersections
	libmesh_assert( static_cast<unsigned>(this->node_cmap_ids[cnt]) == (*it).first );
	
	// Get reference to the set of IDs to be packed into the vector.
	std::set<unsigned>& node_set = (*it).second;
	
	//std::cout << "[" << libMesh::processor_id() << "] node_set.size()=" << node_set.size() << std::endl;

	// Resize the vectors to receive their payload
	this->node_cmap_node_ids[cnt].resize(node_set.size());
	this->node_cmap_proc_ids[cnt].resize(node_set.size());
	
	std::set<unsigned>::iterator node_set_iter = node_set.begin();
	
	// Pack the vectors with node IDs and processor IDs.
	for (unsigned j=0; j<this->node_cmap_node_ids[cnt].size(); ++j, ++node_set_iter)
	  {
	    this->node_cmap_node_ids[cnt][j] = libmesh_node_num_to_exodus[*node_set_iter];//(*node_set_iter) + 1; // Exodus is 1-based
	    this->node_cmap_proc_ids[cnt][j] = (*it).first;
	  }

	cnt++;// increment vector index to go to next processor
      }
  } // end scope for packing 

  // if (_verbose)
  //   libMesh::out << "Done packing." << std::endl;
  
  // Print out the vectors we just packed
  if (_verbose)
    {
      for (unsigned i=0; i<this->node_cmap_node_ids.size(); ++i)
	{
	  libMesh::out << "[" << libMesh::processor_id() << "] nodes communicated to proc " 
		       << this->node_cmap_ids[i]
		       << " = ";
	  for (unsigned j=0; j<this->node_cmap_node_ids[i].size(); ++j)
	    libMesh::out << this->node_cmap_node_ids[i][j] << " ";
	  libMesh::out << std::endl;
	}

      for (unsigned i=0; i<this->node_cmap_node_ids.size(); ++i)
	{
	  libMesh::out << "[" << libMesh::processor_id() << "] processor ID node communicated to = ";
	  for (unsigned j=0; j<this->node_cmap_proc_ids[i].size(); ++j)
	    libMesh::out << this->node_cmap_proc_ids[i][j] << " ";
	  libMesh::out << std::endl;
	}
    }
  
  // Write the packed node communication vectors to file.  
  this->put_node_cmap(this->node_cmap_node_ids,
		      this->node_cmap_proc_ids);
  

  // Write the Nemesis node maps (internal, border, and external nodes) to file.

  // Make sure we don't have any leftover information
  this->node_mapi.clear();
  this->node_mapb.clear();
  this->node_mape.clear();

  // Make sure there's enough space to hold all our node IDs
  this->node_mapi.resize(internal_node_ids.size());
  this->node_mapb.resize(border_node_ids.size());

  // Copy set contents into vectors
  //
  // Can't use insert, since we are copying unsigned's into vector<int>...
  // this->node_mapi.insert(internal_node_ids.begin(), internal_node_ids.end());
  // this->node_mapb.insert(boundary_node_ids.begin(), boundary_node_ids.end());
  {
    unsigned cnt = 0;
    for (std::set<unsigned>::iterator it=internal_node_ids.begin();
	 it != internal_node_ids.end();
	 ++it, ++cnt)
      this->node_mapi[cnt] = libmesh_node_num_to_exodus[*it];// + 1; // Exodus is 1-based!
  }

  {
    unsigned cnt=0;
    for (std::set<unsigned>::iterator it=border_node_ids.begin();
	 it != border_node_ids.end();
	 ++it, ++cnt)
      this->node_mapb[cnt] = libmesh_node_num_to_exodus[*it];// + 1; // Exodus is 1-based!
  }

  // Call the Nemesis API to write these arrays to file.  These node maps don't
  // quite match up to the reference files, similar to the entries in the node
  // communication maps.
  this->put_node_map(this->node_mapi,
		     this->node_mapb,
		     this->node_mape);



  // Write the Nemesis element communication maps, this includes border 
  // element IDs, sides which are on the border, and the processors to which
  // they are to be communicated...

  // Make sure there is no leftover information
  this->elem_cmap_elem_ids.clear();
  this->elem_cmap_side_ids.clear();
  this->elem_cmap_proc_ids.clear();
  
  // Allocate enough space for all our element maps
  this->elem_cmap_elem_ids.resize(this->num_elem_cmaps);
  this->elem_cmap_side_ids.resize(this->num_elem_cmaps);
  this->elem_cmap_proc_ids.resize(this->num_elem_cmaps);
  {
    unsigned cnt=0; // Index into vectors
    for (proc_border_elem_sets_iterator it=proc_border_elem_sets.begin();
	 it != proc_border_elem_sets.end();
	 ++it)
      {
	// Make sure the current elem_cmap_id matches the index in our map of node intersections
	libmesh_assert( static_cast<unsigned>(this->elem_cmap_ids[cnt]) == (*it).first );

	// Get reference to the set of IDs to be packed into the vector
	std::set<std::pair<unsigned,unsigned> >& elem_set = (*it).second;
	
	// Resize the vectors to receive their payload
	this->elem_cmap_elem_ids[cnt].resize(elem_set.size());
	this->elem_cmap_side_ids[cnt].resize(elem_set.size());
	this->elem_cmap_proc_ids[cnt].resize(elem_set.size());

	std::set<std::pair<unsigned,unsigned> >::iterator elem_set_iter = elem_set.begin();

	// Pack the vectors with elem IDs, side IDs, and processor IDs.
	for (unsigned j=0; j<this->elem_cmap_elem_ids[cnt].size(); ++j, ++elem_set_iter)
	  {
	    this->elem_cmap_elem_ids[cnt][j] = libmesh_elem_num_to_exodus[(*elem_set_iter).first];//  + 1; // Elem ID, Exodus is 1-based
	    this->elem_cmap_side_ids[cnt][j] = (*elem_set_iter).second;     // Side ID, this has already been converted above
	    this->elem_cmap_proc_ids[cnt][j] = (*it).first; // All have the same processor ID
	  }

	cnt++;// increment vector index to go to next processor
      }
  } // end scope for packing
  

  // Call the Nemesis API to write the packed vectors to file
  this->put_elem_cmap(this->elem_cmap_elem_ids,
		      this->elem_cmap_side_ids,
		      this->elem_cmap_proc_ids);
  





  // Write the Nemesis element maps (internal and border) to file.
  
  // Make sure we don't have any leftover info
  this->elem_mapi.clear();
  this->elem_mapb.clear();
  
  // Copy set contents into vectors
  this->elem_mapi.resize(internal_elem_ids.size());
  this->elem_mapb.resize(border_elem_ids.size());
  
  {
    unsigned cnt = 0;
    for (std::set<unsigned>::iterator it=internal_elem_ids.begin();
	 it != internal_elem_ids.end();
	 ++it, ++cnt)
      this->elem_mapi[cnt] = libmesh_elem_num_to_exodus[(*it)]; // + 1; // Exodus is 1-based!
  }

  {
    unsigned cnt = 0;
    for (std::set<unsigned>::iterator it=border_elem_ids.begin();
	 it != border_elem_ids.end();
	 ++it, ++cnt)
      this->elem_mapb[cnt] = libmesh_elem_num_to_exodus[(*it)]; // + 1; // Exodus is 1-based!
  }
  
  // Call the Nemesis API to write the internal and border element IDs.
  // FIXME: These ids don't exactly match, similar to the nodes...I assume
  // the numberings in the reference files must use the local node maps
  // that are also in the file?
  this->put_elem_map(this->elem_mapi,
		     this->elem_mapb);

  // FIXME: we also probably need to call essentially the initialize() function in
  // ExodusII_IO_Helper, minus the early return for not being on processor 0...
  // Note that the Exodus writer also incorrectly associates the number of "blocks"
  // written to the Exodus file with the subdomain IDs used in libmesh.
  // "Blocks" in exodus are just supposed to be composed of elements of a particular
  // geometric type, and be unrelated to libmesh subdomain IDs...
  
  // i.e. call something like:


  
///////  //ExodusII_IO_Helper::initialize(title, mesh); //////////////////
  {    
    num_dim = mesh.spatial_dimension();

    // Find the number of nodes... which are all the ones attached to local active elements
    
    num_elem = static_cast<unsigned int>(std::distance (pmesh.active_local_elements_begin(),
                                                        pmesh.active_local_elements_end()));

    /*
      std::vector<short int> unique_side_boundaries;
      std::vector<short int> unique_node_boundaries;

      mesh.boundary_info->build_side_boundary_ids(unique_side_boundaries);
      mesh.boundary_info->build_node_boundary_ids(unique_node_boundaries);
    */
  
//  num_side_sets = unique_side_boundaries.size();
//  num_node_sets = unique_node_boundaries.size();
    
    num_side_sets = 0;  
    num_node_sets = 0;
  
    //loop through element and map between block and element vector
    std::map<subdomain_id_type, std::vector<unsigned int>  > subdomain_map;

    MeshBase::const_element_iterator it = mesh.active_local_elements_begin(); 
    const MeshBase::const_element_iterator end = mesh.active_local_elements_end();
    for (; it != end; ++it)
    {
      Elem * elem = *it;
      subdomain_id_type cur_subdomain = elem->subdomain_id();
      
      if(cur_subdomain == 0)
        cur_subdomain = std::numeric_limits<subdomain_id_type>::max();

      subdomain_map[cur_subdomain].push_back(elem->id());
    }
    
    // For some reason it seems as if this is correct:
    num_elem_blk = num_elem_blks_global;//subdomain_map.size();

    ex_err = exII::ex_put_init(ex_id,
                               title.c_str(),
                               num_dim,
                               num_nodes,
                               num_elem,
                               num_elem_blk,
                               num_node_sets,
                               num_side_sets);
    
    check_err(ex_err, "Error initializing new Nemesis file.");
  }

  write_nodal_coordinates(mesh);
  write_elements(mesh);

  ex_err = exII::ex_update(ex_id);
}


void Nemesis_IO_Helper::write_nodal_coordinates(const MeshBase & mesh)
{
  // Make sure that the reference passed in is really a ParallelMesh
  const ParallelMesh& pmesh = libmesh_cast_ref<const ParallelMesh&>(mesh);
  
  int local_num_nodes = exodus_node_num_to_libmesh.size();
  
  x.resize(local_num_nodes);
  y.resize(local_num_nodes);
  z.resize(local_num_nodes);
  
  // Just loop over our list outputing the nodes the way we built the map
  for (unsigned int i=0; i<local_num_nodes; ++i)
  {
    const Node & node = *mesh.node_ptr(exodus_node_num_to_libmesh[i]);
    x[i]=node(0);
    y[i]=node(1);
    z[i]=node(2);
  }
  ex_err = exII::ex_put_coord(ex_id, &x[0], &y[0], &z[0]);
  check_err(ex_err, "Error writing node coordinates");

  ex_err = exII::ex_put_node_num_map(ex_id, &exodus_node_num_to_libmesh[0]);
  check_err(ex_err, "Error writing node num map");
}

void Nemesis_IO_Helper::write_elements(const MeshBase & mesh)
{
  // Iterate over the map we made earlier to write out the elements and their connectivity
  for(std::map<int, std::vector<int> >::iterator it = block_id_to_elem_connectivity.begin();
      it != block_id_to_elem_connectivity.end();
      ++it)
  {
    int block = (*it).first;
    std::vector<int> & this_block_connectivity = (*it).second;
    std::vector<unsigned int> & elements_in_this_block = subdomain_map[block];

    ExodusII_IO_Helper::ElementMaps em;

    //Use the first element in this block to get representative information.
    //Note that Exodus assumes all elements in a block are of the same type!
    //We are using that same assumption here!
    const ExodusII_IO_Helper::Conversion conv = em.assign_conversion(mesh.elem(elements_in_this_block[0])->type());
    num_nodes_per_elem = mesh.elem(elements_in_this_block[0])->n_nodes();

    ex_err = exII::ex_put_elem_block(ex_id, block, conv.exodus_elem_type().c_str(), elements_in_this_block.size(),num_nodes_per_elem,0);
    check_err(ex_err, "Error writing element block.");
  
    ex_err = exII::ex_put_elem_conn(ex_id, block, &this_block_connectivity[0]);
    check_err(ex_err, "Error writing element connectivities");
  }
  
  ex_err = exII::ex_put_elem_num_map(ex_id, &exodus_elem_num_to_libmesh[0]);
  check_err(ex_err, "Error writing element map");
}

void Nemesis_IO_Helper::write_nodal_solution(const std::vector<Number> & values, const std::vector<std::string> names, int timestep)
{
  int num_vars = names.size();
  int num_values = values.size();

  for (int c=0; c<num_vars; c++)
  {
    std::vector<Number> cur_soln(num_nodes);

    //Copy out this variable's solution
    for(int i=0; i<num_nodes; i++)
      cur_soln[i] = values[exodus_node_num_to_libmesh[i]*num_vars + c];
    
    write_nodal_values(c+1,cur_soln,timestep);
  }
}

std::string Nemesis_IO_Helper::construct_nemesis_filename(const std::string& base_filename)
{
  // Build a filename for this processor.  This code is cut-n-pasted from the read function
  // and should probably be put into a separate function...
  std::ostringstream file_oss;

  // We have to be a little careful here: Nemesis left pads its file
  // numbers based on the largest processor ID, so for example on 128
  // processors, we'd have:
  // mesh.e.128.001
  // mesh.e.128.002
  // ...
  // mesh.e.128.099
  // mesh.e.128.100
  // ...
  // mesh.e.128.127

  // Find the length of the highest processor ID 
  file_oss << (libMesh::n_processors()-1);
  unsigned field_width = file_oss.str().size();
  
  if (_verbose)
    libMesh::out << "field_width=" << field_width << std::endl;

  file_oss.str(""); // reset the string stream
  file_oss << base_filename  
	   << '.' << libMesh::n_processors() 
	   << '.' << std::setfill('0') << std::setw(field_width) << libMesh::processor_id();

  // Return the resulting string
  return file_oss.str();
}

} // namespace libMesh

#endif // #if defined(LIBMESH_HAVE_NEMESIS_API) && defined(LIBMESH_HAVE_EXODUS_API)
