// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include <iomanip>
#include <set>
#include <sstream>

// Libmesh headers
#include "libmesh/nemesis_io_helper.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/parallel.h"
#include "libmesh/mesh_base.h"
#include "libmesh/numeric_vector.h"

#if defined(LIBMESH_HAVE_NEMESIS_API) && defined(LIBMESH_HAVE_EXODUS_API)

namespace libMesh
{


// Initialize the various integer members to zero.  We can check
// these later to see if they've been properly initialized...
// The parent ExodusII_IO_Helper is created with the run_only_on_proc0
// flag set to false, so that we can make use of its functionality
// on multiple processors.
Nemesis_IO_Helper::Nemesis_IO_Helper(const ParallelObject & parent,
                                     bool verbose_in, bool single_precision) :
  ExodusII_IO_Helper(parent, verbose_in, /*run_only_on_proc0=*/false, /*single_precision=*/single_precision),
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
  // responsibility for managing the file's lifetime.  Only call ex_update() if the file was
  // opened for writing!
  if (this->opened_for_writing)
    {
      this->ex_err = exII::ex_update(this->ex_id);
      EX_EXCEPTIONLESS_CHECK_ERR(ex_err, "Error flushing buffers to file.");
    }
  this->close();
}



void Nemesis_IO_Helper::get_init_global()
{
  nemesis_err_flag =
    Nemesis::ne_get_init_global(ex_id,
                                &num_nodes_global,
                                &num_elems_global,
                                &num_elem_blks_global,
                                &num_node_sets_global,
                                &num_side_sets_global);
  EX_CHECK_ERR(nemesis_err_flag, "Error reading initial global data!");

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] " << "num_nodes_global=" << num_nodes_global << std::endl;
      libMesh::out << "[" << this->processor_id() << "] " << "num_elems_global=" << num_elems_global << std::endl;
      libMesh::out << "[" << this->processor_id() << "] " << "num_elem_blks_global=" << num_elem_blks_global << std::endl;
      libMesh::out << "[" << this->processor_id() << "] " << "num_node_sets_global=" << num_node_sets_global << std::endl;
      libMesh::out << "[" << this->processor_id() << "] " << "num_side_sets_global=" << num_side_sets_global << std::endl;
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
                                        &global_sideset_ids[0],
                                        &num_global_side_counts[0],
                                        &num_global_side_df_counts[0]);
      EX_CHECK_ERR(nemesis_err_flag, "Error reading global sideset parameters!");

      if (verbose)
        {
          libMesh::out << "[" << this->processor_id() << "] " << "Global Sideset IDs, Side Counts, and DF counts:" << std::endl;
          for (std::size_t bn=0; bn<global_sideset_ids.size(); ++bn)
            {
              libMesh::out << "  [" << this->processor_id() << "] "
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
                                        &global_nodeset_ids[0],
                                        &num_global_node_counts[0],
                                        &num_global_node_df_counts[0]);
      EX_CHECK_ERR(nemesis_err_flag, "Error reading global nodeset parameters!");

      if (verbose)
        {
          libMesh::out << "[" << this->processor_id() << "] " << "Global Nodeset IDs, Node Counts, and DF counts:" << std::endl;
          for (std::size_t bn=0; bn<global_nodeset_ids.size(); ++bn)
            {
              libMesh::out << "  [" << this->processor_id() << "] "
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

  if (num_elem_blks_global > 0)
    {
      nemesis_err_flag =
        Nemesis::ne_get_eb_info_global(ex_id,
                                       &global_elem_blk_ids[0],
                                       &global_elem_blk_cnts[0]);
      EX_CHECK_ERR(nemesis_err_flag, "Error reading global element block info!");
    }

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] " << "Global Element Block IDs and Counts:" << std::endl;
      for (std::size_t bn=0; bn<global_elem_blk_ids.size(); ++bn)
        {
          libMesh::out << "  [" << this->processor_id() << "] "
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
  EX_CHECK_ERR(nemesis_err_flag, "Error reading initial info!");

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] " << "num_proc=" << num_proc << std::endl;
      libMesh::out << "[" << this->processor_id() << "] " << "num_proc_in_file=" << num_proc_in_file << std::endl;
      libMesh::out << "[" << this->processor_id() << "] " << "ftype=" << ftype << std::endl;
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
                                  this->processor_id() // The ID of the processor for which info is to be read
                                  );
  EX_CHECK_ERR(nemesis_err_flag, "Error reading load balance parameters!");


  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] " << "num_internal_nodes=" << num_internal_nodes << std::endl;
      libMesh::out << "[" << this->processor_id() << "] " << "num_border_nodes=" << num_border_nodes << std::endl;
      libMesh::out << "[" << this->processor_id() << "] " << "num_external_nodes=" << num_external_nodes << std::endl;
      libMesh::out << "[" << this->processor_id() << "] " << "num_internal_elems=" << num_internal_elems << std::endl;
      libMesh::out << "[" << this->processor_id() << "] " << "num_border_elems=" << num_border_elems << std::endl;
      libMesh::out << "[" << this->processor_id() << "] " << "num_node_cmaps=" << num_node_cmaps << std::endl;
      libMesh::out << "[" << this->processor_id() << "] " << "num_elem_cmaps=" << num_elem_cmaps << std::endl;
    }
}



void Nemesis_IO_Helper::get_elem_map()
{
  elem_mapi.resize(num_internal_elems);
  elem_mapb.resize(num_border_elems);

  nemesis_err_flag =
    Nemesis::ne_get_elem_map(ex_id,
                             elem_mapi.empty() ? libmesh_nullptr : &elem_mapi[0],
                             elem_mapb.empty() ? libmesh_nullptr : &elem_mapb[0],
                             this->processor_id()
                             );
  EX_CHECK_ERR(nemesis_err_flag, "Error reading element maps!");


  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] elem_mapi[i] = ";
      for (unsigned int i=0; i< static_cast<unsigned int>(num_internal_elems-1); ++i)
        libMesh::out << elem_mapi[i] << ", ";
      libMesh::out << "... " << elem_mapi.back() << std::endl;

      libMesh::out << "[" << this->processor_id() << "] elem_mapb[i] = ";
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
                             node_mapi.empty() ? libmesh_nullptr : &node_mapi[0],
                             node_mapb.empty() ? libmesh_nullptr : &node_mapb[0],
                             node_mape.empty() ? libmesh_nullptr : &node_mape[0],
                             this->processor_id()
                             );
  EX_CHECK_ERR(nemesis_err_flag, "Error reading node maps!");

  if (verbose)
    {
      // Remark: The Exodus/Nemesis node numbring is always (?) 1-based!  This means the first interior node id will
      // always be == 1.
      libMesh::out << "[" << this->processor_id() << "] " << "first interior node id=" << node_mapi[0] << std::endl;
      libMesh::out << "[" << this->processor_id() << "] " << "last interior node id=" << node_mapi.back() << std::endl;

      libMesh::out << "[" << this->processor_id() << "] " << "first boundary node id=" << node_mapb[0] << std::endl;
      libMesh::out << "[" << this->processor_id() << "] " << "last boundary node id=" << node_mapb.back() << std::endl;

      // The number of external nodes is sometimes zero, don't try to access
      // node_mape.back() in this case!
      if (num_external_nodes > 0)
        {
          libMesh::out << "[" << this->processor_id() << "] " << "first external node id=" << node_mape[0] << std::endl;
          libMesh::out << "[" << this->processor_id() << "] " << "last external node id=" << node_mape.back() << std::endl;
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
                                node_cmap_ids.empty()       ? libmesh_nullptr : &node_cmap_ids[0],
                                node_cmap_node_cnts.empty() ? libmesh_nullptr : &node_cmap_node_cnts[0],
                                elem_cmap_ids.empty()       ? libmesh_nullptr : &elem_cmap_ids[0],
                                elem_cmap_elem_cnts.empty() ? libmesh_nullptr : &elem_cmap_elem_cnts[0],
                                this->processor_id());
  EX_CHECK_ERR(nemesis_err_flag, "Error reading cmap parameters!");


  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] ";
      for (std::size_t i=0; i<node_cmap_ids.size(); ++i)
        libMesh::out << "node_cmap_ids["<<i<<"]=" << node_cmap_ids[i] << " ";
      libMesh::out << std::endl;

      libMesh::out << "[" << this->processor_id() << "] ";
      for (std::size_t i=0; i<node_cmap_node_cnts.size(); ++i)
        libMesh::out << "node_cmap_node_cnts["<<i<<"]=" << node_cmap_node_cnts[i] << " ";
      libMesh::out << std::endl;

      libMesh::out << "[" << this->processor_id() << "] ";
      for (std::size_t i=0; i<elem_cmap_ids.size(); ++i)
        libMesh::out << "elem_cmap_ids["<<i<<"]=" << elem_cmap_ids[i] << " ";
      libMesh::out << std::endl;

      libMesh::out << "[" << this->processor_id() << "] ";
      for (std::size_t i=0; i<elem_cmap_elem_cnts.size(); ++i)
        libMesh::out << "elem_cmap_elem_cnts["<<i<<"]=" << elem_cmap_elem_cnts[i] << " ";
      libMesh::out << std::endl;
    }
}



void Nemesis_IO_Helper::get_node_cmap()
{
  node_cmap_node_ids.resize(num_node_cmaps);
  node_cmap_proc_ids.resize(num_node_cmaps);

  for (std::size_t i=0; i<node_cmap_node_ids.size(); ++i)
    {
      node_cmap_node_ids[i].resize(node_cmap_node_cnts[i]);
      node_cmap_proc_ids[i].resize(node_cmap_node_cnts[i]);

      // Don't call ne_get_node_cmap() if there is nothing there to
      // get, Nemesis throws an error in this case.
      if (node_cmap_node_cnts[i] > 0)
        {
          nemesis_err_flag =
            Nemesis::ne_get_node_cmap(ex_id,
                                      node_cmap_ids[i],
                                      &node_cmap_node_ids[i][0],
                                      &node_cmap_proc_ids[i][0],
                                      this->processor_id());
          EX_CHECK_ERR(nemesis_err_flag, "Error reading node cmap node and processor ids!");
        }

      if (verbose)
        {
          libMesh::out << "[" << this->processor_id() << "] node_cmap_node_ids["<<i<<"]=";
          for (std::size_t j=0; j<node_cmap_node_ids[i].size(); ++j)
            libMesh::out << node_cmap_node_ids[i][j] << " ";
          libMesh::out << std::endl;

          // This is basically a vector, all entries of which are = node_cmap_ids[i]
          // Not sure if it's always guaranteed to be that or what...
          libMesh::out << "[" << this->processor_id() << "] node_cmap_proc_ids["<<i<<"]=";
          for (std::size_t j=0; j<node_cmap_proc_ids[i].size(); ++j)
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

  for (std::size_t i=0; i<elem_cmap_elem_ids.size(); ++i)
    {
      elem_cmap_elem_ids[i].resize(elem_cmap_elem_cnts[i]);
      elem_cmap_side_ids[i].resize(elem_cmap_elem_cnts[i]);
      elem_cmap_proc_ids[i].resize(elem_cmap_elem_cnts[i]);

      if (elem_cmap_elem_cnts[i] > 0)
        {
          nemesis_err_flag =
            Nemesis::ne_get_elem_cmap(ex_id,
                                      elem_cmap_ids[i],
                                      &elem_cmap_elem_ids[i][0],
                                      &elem_cmap_side_ids[i][0],
                                      &elem_cmap_proc_ids[i][0],
                                      this->processor_id());
          EX_CHECK_ERR(nemesis_err_flag, "Error reading elem cmap elem, side, and processor ids!");
        }

      if (verbose)
        {
          libMesh::out << "[" << this->processor_id() << "] elem_cmap_elem_ids["<<i<<"]=";
          for (std::size_t j=0; j<elem_cmap_elem_ids[i].size(); ++j)
            libMesh::out << elem_cmap_elem_ids[i][j] << " ";
          libMesh::out << std::endl;

          // These must be the (local) side IDs (in the ExodusII face numbering scheme)
          // of the sides shared across processors.
          libMesh::out << "[" << this->processor_id() << "] elem_cmap_side_ids["<<i<<"]=";
          for (std::size_t j=0; j<elem_cmap_side_ids[i].size(); ++j)
            libMesh::out << elem_cmap_side_ids[i][j] << " ";
          libMesh::out << std::endl;

          // This is basically a vector, all entries of which are = elem_cmap_ids[i]
          // Not sure if it's always guaranteed to be that or what...
          libMesh::out << "[" << this->processor_id() << "] elem_cmap_proc_ids["<<i<<"]=";
          for (std::size_t j=0; j<elem_cmap_proc_ids[i].size(); ++j)
            libMesh::out << elem_cmap_proc_ids[i][j] << " ";
          libMesh::out << std::endl;
        }
    }
}




void Nemesis_IO_Helper::put_init_info(unsigned num_proc_in,
                                      unsigned num_proc_in_file_in,
                                      const char * ftype_in)
{
  nemesis_err_flag =
    Nemesis::ne_put_init_info(ex_id,
                              num_proc_in,
                              num_proc_in_file_in,
                              const_cast<char *>(ftype_in));

  EX_CHECK_ERR(nemesis_err_flag, "Error writing initial information!");
}




void Nemesis_IO_Helper::put_init_global(dof_id_type num_nodes_global_in,
                                        dof_id_type num_elems_global_in,
                                        unsigned num_elem_blks_global_in,
                                        unsigned num_node_sets_global_in,
                                        unsigned num_side_sets_global_in)
{
  nemesis_err_flag =
    Nemesis::ne_put_init_global(ex_id,
                                num_nodes_global_in,
                                num_elems_global_in,
                                num_elem_blks_global_in,
                                num_node_sets_global_in,
                                num_side_sets_global_in);

  EX_CHECK_ERR(nemesis_err_flag, "Error writing initial global data!");
}



void Nemesis_IO_Helper::put_eb_info_global(std::vector<int> & global_elem_blk_ids_in,
                                           std::vector<int> & global_elem_blk_cnts_in)
{
  nemesis_err_flag =
    Nemesis::ne_put_eb_info_global(ex_id,
                                   &global_elem_blk_ids_in[0],
                                   &global_elem_blk_cnts_in[0]);

  EX_CHECK_ERR(nemesis_err_flag, "Error writing global element block information!");
}




void Nemesis_IO_Helper::put_ns_param_global(std::vector<int> & global_nodeset_ids_in,
                                            std::vector<int> & num_global_node_counts_in,
                                            std::vector<int> & num_global_node_df_counts_in)
{
  // Only add nodesets if there are some
  if (global_nodeset_ids.size())
    {
      nemesis_err_flag =
        Nemesis::ne_put_ns_param_global(ex_id,
                                        &global_nodeset_ids_in[0],
                                        &num_global_node_counts_in[0],
                                        &num_global_node_df_counts_in[0]);
    }

  EX_CHECK_ERR(nemesis_err_flag, "Error writing global nodeset parameters!");
}




void Nemesis_IO_Helper::put_ss_param_global(std::vector<int> & global_sideset_ids_in,
                                            std::vector<int> & num_global_side_counts_in,
                                            std::vector<int> & num_global_side_df_counts_in)
{
  // Only add sidesets if there are some
  if (global_sideset_ids.size())
    {
      nemesis_err_flag =
        Nemesis::ne_put_ss_param_global(ex_id,
                                        &global_sideset_ids_in[0],
                                        &num_global_side_counts_in[0],
                                        &num_global_side_df_counts_in[0]);
    }

  EX_CHECK_ERR(nemesis_err_flag, "Error writing global sideset parameters!");
}




void Nemesis_IO_Helper::put_loadbal_param(unsigned num_internal_nodes_in,
                                          unsigned num_border_nodes_in,
                                          unsigned num_external_nodes_in,
                                          unsigned num_internal_elems_in,
                                          unsigned num_border_elems_in,
                                          unsigned num_node_cmaps_in,
                                          unsigned num_elem_cmaps_in)
{
  nemesis_err_flag =
    Nemesis::ne_put_loadbal_param(ex_id,
                                  num_internal_nodes_in,
                                  num_border_nodes_in,
                                  num_external_nodes_in,
                                  num_internal_elems_in,
                                  num_border_elems_in,
                                  num_node_cmaps_in,
                                  num_elem_cmaps_in,
                                  this->processor_id());

  EX_CHECK_ERR(nemesis_err_flag, "Error writing loadbal parameters!");
}





void Nemesis_IO_Helper::put_cmap_params(std::vector<int> & node_cmap_ids_in,
                                        std::vector<int> & node_cmap_node_cnts_in,
                                        std::vector<int> & elem_cmap_ids_in,
                                        std::vector<int> & elem_cmap_elem_cnts_in)
{
  // We might not have cmaps on every processor in some corner
  // cases
  if (node_cmap_ids.size())
    {
      nemesis_err_flag =
        Nemesis::ne_put_cmap_params(ex_id,
                                    &node_cmap_ids_in[0],
                                    &node_cmap_node_cnts_in[0],
                                    &elem_cmap_ids_in[0],
                                    &elem_cmap_elem_cnts_in[0],
                                    this->processor_id());
    }

  EX_CHECK_ERR(nemesis_err_flag, "Error writing cmap parameters!");
}




void Nemesis_IO_Helper::put_node_cmap(std::vector<std::vector<int> > & node_cmap_node_ids_in,
                                      std::vector<std::vector<int> > & node_cmap_proc_ids_in)
{

  // Print to screen what we are about to print to Nemesis file
  if (verbose)
    {
      for (std::size_t i=0; i<node_cmap_node_ids_in.size(); ++i)
        {
          libMesh::out << "[" << this->processor_id() << "] put_node_cmap() : nodes communicated to proc "
                       << this->node_cmap_ids[i]
                       << " = ";
          for (std::size_t j=0; j<node_cmap_node_ids_in[i].size(); ++j)
            libMesh::out << node_cmap_node_ids_in[i][j] << " ";
          libMesh::out << std::endl;
        }

      for (std::size_t i=0; i<node_cmap_node_ids_in.size(); ++i)
        {
          libMesh::out << "[" << this->processor_id() << "] put_node_cmap() : processor IDs = ";
          for (std::size_t j=0; j<node_cmap_proc_ids_in[i].size(); ++j)
            libMesh::out << node_cmap_proc_ids_in[i][j] << " ";
          libMesh::out << std::endl;
        }
    }

  for (std::size_t i=0; i<node_cmap_node_ids_in.size(); ++i)
    {
      nemesis_err_flag =
        Nemesis::ne_put_node_cmap(ex_id,
                                  this->node_cmap_ids[i],
                                  &node_cmap_node_ids_in[i][0],
                                  &node_cmap_proc_ids_in[i][0],
                                  this->processor_id());

      EX_CHECK_ERR(nemesis_err_flag, "Error writing node communication map to file!");
    }
}




void Nemesis_IO_Helper::put_node_map(std::vector<int> & node_mapi_in,
                                     std::vector<int> & node_mapb_in,
                                     std::vector<int> & node_mape_in)
{
  nemesis_err_flag =
    Nemesis::ne_put_node_map(ex_id,
                             node_mapi_in.empty() ? libmesh_nullptr : &node_mapi_in[0],
                             node_mapb_in.empty() ? libmesh_nullptr : &node_mapb_in[0],
                             node_mape_in.empty() ? libmesh_nullptr : &node_mape_in[0], // Don't take address of empty vector...
                             this->processor_id());

  EX_CHECK_ERR(nemesis_err_flag, "Error writing Nemesis internal and border node maps to file!");
}




void Nemesis_IO_Helper::put_elem_cmap(std::vector<std::vector<int> > & elem_cmap_elem_ids_in,
                                      std::vector<std::vector<int> > & elem_cmap_side_ids_in,
                                      std::vector<std::vector<int> > & elem_cmap_proc_ids_in)
{
  for (std::size_t i=0; i<elem_cmap_ids.size(); ++i)
    {
      nemesis_err_flag =
        Nemesis::ne_put_elem_cmap(ex_id,
                                  this->elem_cmap_ids[i],
                                  &elem_cmap_elem_ids_in[i][0],
                                  &elem_cmap_side_ids_in[i][0],
                                  &elem_cmap_proc_ids_in[i][0],
                                  this->processor_id());

      EX_CHECK_ERR(nemesis_err_flag, "Error writing elem communication map to file!");
    }
}




void Nemesis_IO_Helper::put_elem_map(std::vector<int> & elem_mapi_in,
                                     std::vector<int> & elem_mapb_in)
{
  nemesis_err_flag =
    Nemesis::ne_put_elem_map(ex_id,
                             elem_mapi_in.empty() ? libmesh_nullptr : &elem_mapi_in[0],
                             elem_mapb_in.empty() ? libmesh_nullptr : &elem_mapb_in[0],
                             this->processor_id());

  EX_CHECK_ERR(nemesis_err_flag, "Error writing Nemesis internal and border element maps to file!");
}






void Nemesis_IO_Helper::put_n_coord(unsigned start_node_num,
                                    unsigned num_nodes_in,
                                    std::vector<Real> & x_coor,
                                    std::vector<Real> & y_coor,
                                    std::vector<Real> & z_coor)
{
  nemesis_err_flag =
    Nemesis::ne_put_n_coord(ex_id,
                            start_node_num,
                            num_nodes_in,
                            &x_coor[0],
                            &y_coor[0],
                            &z_coor[0]);

  EX_CHECK_ERR(nemesis_err_flag, "Error writing coords to file!");
}








// Note: we can't reuse the ExodusII_IO_Helper code directly, since it only runs
// on processor 0.  TODO: We could have the body of this function as a separate
// function and then ExodusII_IO_Helper would only call it if on processor 0...
void Nemesis_IO_Helper::create(std::string filename)
{
  // Fall back on double precision when necessary since ExodusII
  // doesn't seem to support long double
  int
    comp_ws = 0,
    io_ws = 0;

  if (_single_precision)
    {
      comp_ws = sizeof(float);
      io_ws = sizeof(float);
    }
  else
    {
      comp_ws = cast_int<int>(std::min(sizeof(Real), sizeof(double)));
      io_ws = cast_int<int>(std::min(sizeof(Real), sizeof(double)));
    }

  this->ex_id = exII::ex_create(filename.c_str(), EX_CLOBBER, &comp_ws, &io_ws);

  EX_CHECK_ERR(ex_id, "Error creating Nemesis mesh file.");

  if (verbose)
    libMesh::out << "File created successfully." << std::endl;

  this->opened_for_writing = true;
}








void Nemesis_IO_Helper::initialize(std::string title_in, const MeshBase & mesh, bool /*use_discontinuous*/)
{
  // Make sure that the reference passed in is really a DistributedMesh
  // const DistributedMesh & pmesh = cast_ref<const DistributedMesh &>(mesh);
  const MeshBase & pmesh = mesh;

  // According to Nemesis documentation, first call when writing should be to
  // ne_put_init_info().  Our reader doesn't actually call this, but we should
  // strive to be as close to a normal nemesis file as possible...
  this->put_init_info(this->n_processors(), 1, "p");


  // Gather global "initial" information for Nemesis.  This consists of
  // three parts labelled I, II, and III below...

  //
  // I.) Need to compute the number of global element blocks.  To be consistent with
  // Exodus, we also incorrectly associate the number of element blocks with the
  // number of libmesh subdomains...
  //
  this->compute_num_global_elem_blocks(pmesh);

  //
  // II.) Determine the global number of nodesets by communication.
  // This code relies on BoundaryInfo storing side and node
  // boundary IDs separately at the time they are added to the
  // BoundaryInfo object.
  //
  this->compute_num_global_nodesets(pmesh);

  //
  // III.) Need to compute the global number of sidesets by communication:
  // This code relies on BoundaryInfo storing side and node
  // boundary IDs separately at the time they are added to the
  // BoundaryInfo object.
  //
  this->compute_num_global_sidesets(pmesh);

  // Now write the global data obtained in steps I, II, and III to the Nemesis file
  this->put_init_global(pmesh.parallel_n_nodes(),
                        pmesh.parallel_n_elem(),
                        this->num_elem_blks_global, /* I.   */
                        this->num_node_sets_global, /* II.  */
                        this->num_side_sets_global  /* III. */
                        );

  // Next, we'll write global element block information to the file.  This was already
  // gathered in step I. above
  this->put_eb_info_global(this->global_elem_blk_ids,
                           this->global_elem_blk_cnts);


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


  // Before we go any further we need to derive consistent node and
  // element numbering schemes for all local elems and nodes connected
  // to local elements.
  //
  // Must be called *after* the local_subdomain_counts map has been constructed
  // by the compute_num_global_elem_blocks() function!
  this->build_element_and_node_maps(pmesh);

  // Next step is to write "load balance" parameters.  Several things need to
  // be computed first though...

  // First we'll collect IDs of border nodes.
  this->compute_border_node_ids(pmesh);

  // Next we'll collect numbers of internal and border elements, and internal nodes.
  // Note: "A border node does not a border element make...", that is, just because one
  // of an element's nodes has been identified as a border node, the element is not
  // necessarily a border element.  It must have a side on the boundary between processors,
  // i.e. have a face neighbor with a different processor id...
  this->compute_internal_and_border_elems_and_internal_nodes(pmesh);

  // Finally we are ready to write the loadbal information to the file
  this->put_loadbal_param(this->num_internal_nodes,
                          this->num_border_nodes,
                          this->num_external_nodes,
                          this->num_internal_elems,
                          this->num_border_elems,
                          this->num_node_cmaps,
                          this->num_elem_cmaps);


  // Now we need to compute the "communication map" parameters.  These are basically
  // lists of nodes and elements which need to be communicated between different processors
  // when the mesh file is read back in.
  this->compute_communication_map_parameters();

  // Write communication map parameters to file.
  this->put_cmap_params(this->node_cmap_ids,
                        this->node_cmap_node_cnts,
                        this->elem_cmap_ids,
                        this->elem_cmap_elem_cnts);


  // Ready the node communication maps.  The node IDs which
  // are communicated are the ones currently stored in
  // proc_nodes_touched_intersections.
  this->compute_node_communication_maps();

  // Write the packed node communication vectors to file.
  this->put_node_cmap(this->node_cmap_node_ids,
                      this->node_cmap_proc_ids);


  // Ready the node maps.  These have nothing to do with communiction, they map
  // the nodes to internal, border, and external nodes in the file.
  this->compute_node_maps();

  // Call the Nemesis API to write the node maps to file.
  this->put_node_map(this->node_mapi,
                     this->node_mapb,
                     this->node_mape);



  // Ready the element communication maps.  This includes border
  // element IDs, sides which are on the border, and the processors to which
  // they are to be communicated...
  this->compute_elem_communication_maps();



  // Call the Nemesis API to write the packed element communication maps vectors to file
  this->put_elem_cmap(this->elem_cmap_elem_ids,
                      this->elem_cmap_side_ids,
                      this->elem_cmap_proc_ids);






  // Ready the Nemesis element maps (internal and border) for writing to file.
  this->compute_element_maps();

  // Call the Nemesis API to write the internal and border element IDs.
  this->put_elem_map(this->elem_mapi,
                     this->elem_mapb);


  // Now write Exodus-specific initialization information, some of which is
  // different when you are using Nemesis.
  this->write_exodus_initialization_info(pmesh, title_in);
} // end initialize()






void Nemesis_IO_Helper::write_exodus_initialization_info(const MeshBase & pmesh,
                                                         const std::string & title_in)
{
  // This follows the convention of Exodus: we always write out the mesh as LIBMESH_DIM-dimensional,
  // even if it is 2D...
  this->num_dim = LIBMESH_DIM;

  this->num_elem = static_cast<unsigned int>(std::distance (pmesh.active_local_elements_begin(),
                                                            pmesh.active_local_elements_end()));

  // Exodus will also use *global* number of side and node sets,
  // though it will not write out entries for all of them...
  this->num_side_sets =
    cast_int<int>(this->global_sideset_ids.size());
  this->num_node_sets =
    cast_int<int>(this->global_nodeset_ids.size());

  // We need to write the global number of blocks, even though this processor might not have
  // elements in some of them!
  this->num_elem_blk = this->num_elem_blks_global;

  ex_err = exII::ex_put_init(ex_id,
                             title_in.c_str(),
                             this->num_dim,
                             this->num_nodes,
                             this->num_elem,
                             this->num_elem_blk,
                             this->num_node_sets,
                             this->num_side_sets);

  EX_CHECK_ERR(ex_err, "Error initializing new Nemesis file.");
}





void Nemesis_IO_Helper::compute_element_maps()
{
  // Make sure we don't have any leftover info
  this->elem_mapi.clear();
  this->elem_mapb.clear();

  // Copy set contents into vectors
  this->elem_mapi.resize(this->internal_elem_ids.size());
  this->elem_mapb.resize(this->border_elem_ids.size());

  {
    unsigned cnt = 0;
    std::set<unsigned>::iterator
      it = this->internal_elem_ids.begin(),
      end = this->internal_elem_ids.end();

    for (; it != end; ++it, ++cnt)
      {
        std::map<int, int>::iterator elem_it = libmesh_elem_num_to_exodus.find(*it);
        if (elem_it == libmesh_elem_num_to_exodus.end())
          libmesh_error_msg("Elem number " << *it << " not found in libmesh_elem_num_to_exodus map.");
        this->elem_mapi[cnt] = elem_it->second;
      }
  }

  {
    unsigned cnt = 0;
    std::set<unsigned>::iterator
      it = this->border_elem_ids.begin(),
      end = this->border_elem_ids.end();

    for (; it != end; ++it, ++cnt)
      {
        std::map<int, int>::iterator elem_it = libmesh_elem_num_to_exodus.find(*it);
        if (elem_it == libmesh_elem_num_to_exodus.end())
          libmesh_error_msg("Elem number " << *it << " not found in libmesh_elem_num_to_exodus map.");
        this->elem_mapb[cnt] = elem_it->second;
      }
  }
}



void Nemesis_IO_Helper::compute_elem_communication_maps()
{
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
    proc_border_elem_sets_iterator
      it = this->proc_border_elem_sets.begin(),
      end = this->proc_border_elem_sets.end();

    for (; it != end; ++it)
      {
        // Make sure the current elem_cmap_id matches the index in our map of node intersections
        libmesh_assert_equal_to (static_cast<unsigned>(this->elem_cmap_ids[cnt]), it->first);

        // Get reference to the set of IDs to be packed into the vector
        std::set<std::pair<unsigned,unsigned> > & elem_set = it->second;

        // Resize the vectors to receive their payload
        this->elem_cmap_elem_ids[cnt].resize(elem_set.size());
        this->elem_cmap_side_ids[cnt].resize(elem_set.size());
        this->elem_cmap_proc_ids[cnt].resize(elem_set.size());

        std::set<std::pair<unsigned,unsigned> >::iterator elem_set_iter = elem_set.begin();

        // Pack the vectors with elem IDs, side IDs, and processor IDs.
        for (std::size_t j=0; j<this->elem_cmap_elem_ids[cnt].size(); ++j, ++elem_set_iter)
          {
            std::map<int, int>::iterator elem_it = libmesh_elem_num_to_exodus.find((*elem_set_iter).first);

            if (elem_it == libmesh_elem_num_to_exodus.end())
              libmesh_error_msg("Elem number " << (*elem_set_iter).first << " not found in libmesh_elem_num_to_exodus map.");

            this->elem_cmap_elem_ids[cnt][j] = elem_it->second;
            this->elem_cmap_side_ids[cnt][j] = (*elem_set_iter).second;     // Side ID, this has already been converted above
            this->elem_cmap_proc_ids[cnt][j] = it->first; // All have the same processor ID
          }

        // increment vector index to go to next processor
        cnt++;
      }
  } // end scope for packing
}





void Nemesis_IO_Helper::compute_node_maps()
{
  // Make sure we don't have any leftover information
  this->node_mapi.clear();
  this->node_mapb.clear();
  this->node_mape.clear();

  // Make sure there's enough space to hold all our node IDs
  this->node_mapi.resize(this->internal_node_ids.size());
  this->node_mapb.resize(this->border_node_ids.size());

  // Copy set contents into vectors
  {
    unsigned cnt = 0;
    std::set<unsigned>::iterator
      it = this->internal_node_ids.begin(),
      end = this->internal_node_ids.end();

    for (; it != end; ++it, ++cnt)
      {
        std::map<int, int>::iterator node_it = libmesh_node_num_to_exodus.find(*it);
        if (node_it == libmesh_node_num_to_exodus.end())
          libmesh_error_msg("Node number " << *it << " not found in libmesh_node_num_to_exodus map.");
        this->node_mapi[cnt] = node_it->second;
      }
  }

  {
    unsigned cnt=0;
    std::set<unsigned>::iterator
      it = this->border_node_ids.begin(),
      end = this->border_node_ids.end();

    for (; it != end; ++it, ++cnt)
      {
        std::map<int, int>::iterator node_it = libmesh_node_num_to_exodus.find(*it);
        if (node_it == libmesh_node_num_to_exodus.end())
          libmesh_error_msg("Node number " << *it << " not found in libmesh_node_num_to_exodus map.");
        this->node_mapb[cnt] = node_it->second;
      }
  }
}





void Nemesis_IO_Helper::compute_node_communication_maps()
{
  // Make sure there's no left-over information
  this->node_cmap_node_ids.clear();
  this->node_cmap_proc_ids.clear();

  // Allocate enough space for all our node maps
  this->node_cmap_node_ids.resize(this->num_node_cmaps);
  this->node_cmap_proc_ids.resize(this->num_node_cmaps);
  {
    unsigned cnt=0; // Index into vectors
    proc_nodes_touched_iterator
      it = this->proc_nodes_touched_intersections.begin(),
      end = this->proc_nodes_touched_intersections.end();

    for (; it != end; ++it)
      {
        // Make sure the current node_cmap_id matches the index in our map of node intersections
        libmesh_assert_equal_to (static_cast<unsigned>(this->node_cmap_ids[cnt]), it->first);

        // Get reference to the set of IDs to be packed into the vector.
        std::set<unsigned> & node_set = it->second;

        // Resize the vectors to receive their payload
        this->node_cmap_node_ids[cnt].resize(node_set.size());
        this->node_cmap_proc_ids[cnt].resize(node_set.size());

        std::set<unsigned>::iterator node_set_iter = node_set.begin();

        // Pack the vectors with node IDs and processor IDs.
        for (std::size_t j=0; j<this->node_cmap_node_ids[cnt].size(); ++j, ++node_set_iter)
          {
            std::map<int, int>::iterator node_it = libmesh_node_num_to_exodus.find(*node_set_iter);
            if (node_it == libmesh_node_num_to_exodus.end())
              libmesh_error_msg("Node number " << *node_set_iter << " not found in libmesh_node_num_to_exodus map.");

            this->node_cmap_node_ids[cnt][j] = node_it->second;
            this->node_cmap_proc_ids[cnt][j] = it->first;
          }

        // increment vector index to go to next processor
        cnt++;
      }
  } // end scope for packing

  // Print out the vectors we just packed
  if (verbose)
    {
      for (std::size_t i=0; i<this->node_cmap_node_ids.size(); ++i)
        {
          libMesh::out << "[" << this->processor_id() << "] nodes communicated to proc "
                       << this->node_cmap_ids[i]
                       << " = ";
          for (std::size_t j=0; j<this->node_cmap_node_ids[i].size(); ++j)
            libMesh::out << this->node_cmap_node_ids[i][j] << " ";
          libMesh::out << std::endl;
        }

      for (std::size_t i=0; i<this->node_cmap_node_ids.size(); ++i)
        {
          libMesh::out << "[" << this->processor_id() << "] processor ID node communicated to = ";
          for (std::size_t j=0; j<this->node_cmap_proc_ids[i].size(); ++j)
            libMesh::out << this->node_cmap_proc_ids[i][j] << " ";
          libMesh::out << std::endl;
        }
    }
}




void Nemesis_IO_Helper::compute_communication_map_parameters()
{
  // For the nodes, these are the number of entries in the sets in proc_nodes_touched_intersections
  // map computed above.  Note: this map does not contain self-intersections so we can loop over it
  // directly.
  this->node_cmap_node_cnts.clear(); // Make sure we don't have any leftover information...
  this->node_cmap_ids.clear();       // Make sure we don't have any leftover information...
  this->node_cmap_node_cnts.resize(this->num_node_cmaps);
  this->node_cmap_ids.resize(this->num_node_cmaps);

  {
    unsigned cnt=0; // Index into the vector
    proc_nodes_touched_iterator
      it = this->proc_nodes_touched_intersections.begin(),
      end = this->proc_nodes_touched_intersections.end();

    for (; it != end; ++it)
      {
        this->node_cmap_ids[cnt] = it->first; // The ID of the proc we communicate with
        this->node_cmap_node_cnts[cnt] = cast_int<int>(it->second.size()); // The number of nodes we communicate
        cnt++; // increment vector index!
      }
  }

  // Print the packed vectors we just filled
  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] node_cmap_node_cnts = ";
      for (std::size_t i=0; i<node_cmap_node_cnts.size(); ++i)
        libMesh::out << node_cmap_node_cnts[i] << ", ";
      libMesh::out << std::endl;

      libMesh::out << "[" << this->processor_id() << "] node_cmap_ids = ";
      for (std::size_t i=0; i<node_cmap_ids.size(); ++i)
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
    proc_border_elem_sets_iterator
      it = this->proc_border_elem_sets.begin(),
      end = this->proc_border_elem_sets.end();

    for (; it != end; ++it)
      {
        this->elem_cmap_ids[cnt] = it->first; // The ID of the proc we communicate with
        this->elem_cmap_elem_cnts[cnt] = cast_int<int>(it->second.size()); // The number of elems we communicate to/from that proc
        cnt++; // increment vector index!
      }
  }

  // Print the packed vectors we just filled
  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] elem_cmap_elem_cnts = ";
      for (std::size_t i=0; i<elem_cmap_elem_cnts.size(); ++i)
        libMesh::out << elem_cmap_elem_cnts[i] << ", ";
      libMesh::out << std::endl;

      libMesh::out << "[" << this->processor_id() << "] elem_cmap_ids = ";
      for (std::size_t i=0; i<elem_cmap_ids.size(); ++i)
        libMesh::out << elem_cmap_ids[i] << ", ";
      libMesh::out << std::endl;
    }
}




void
Nemesis_IO_Helper::compute_internal_and_border_elems_and_internal_nodes(const MeshBase & pmesh)
{
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

  MeshBase::const_element_iterator elem_it = pmesh.active_local_elements_begin();
  MeshBase::const_element_iterator elem_end = pmesh.active_local_elements_end();

  for (; elem_it != elem_end; ++elem_it)
    {
      const Elem * elem = *elem_it;

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
          this->nodes_attached_to_local_elems.insert(elem->node_id(node));
        } // end loop over element's nodes

      // Loop over element's neighbors, see if it has a neighbor which is off-processor
      for (unsigned int n=0; n<elem->n_neighbors(); ++n)
        {
          if (elem->neighbor_ptr(n) != libmesh_nullptr)
            {
              unsigned neighbor_proc_id = elem->neighbor_ptr(n)->processor_id();

              // If my neighbor has a different processor ID, I must be a border element.
              // Also track the neighboring processor ID if it is are different from our processor ID
              if (neighbor_proc_id != this->processor_id())
                {
                  is_border_elem = true;
                  neighboring_processor_ids.insert(neighbor_proc_id);

                  // Convert libmesh side(n) of this element into a side ID for Nemesis
                  unsigned nemesis_side_id = conv.get_inverse_side_map(n);

                  if (verbose)
                    libMesh::out << "[" << this->processor_id() << "] LibMesh side "
                                 << n
                                 << " mapped to (1-based) Exodus side "
                                 << nemesis_side_id
                                 << std::endl;

                  // Add this element's ID and the ID of the side which is on the boundary
                  // to the set of border elements for this processor.
                  // Note: if the set does not already exist, this creates it.
                  this->proc_border_elem_sets[ neighbor_proc_id ].insert( std::make_pair(elem->id(), nemesis_side_id) );
                }
            }
        } // end for loop over neighbors

      // If we're on a border element, add it to the set
      if (is_border_elem)
        this->border_elem_ids.insert( elem->id() );

    } // end for loop over active local elements

  // Take the set_difference between all elements and border elements to get internal
  // element IDs
  std::set_difference(all_elem_ids.begin(), all_elem_ids.end(),
                      this->border_elem_ids.begin(), this->border_elem_ids.end(),
                      std::inserter(this->internal_elem_ids, this->internal_elem_ids.end()));

  // Take the set_difference between all nodes and border nodes to get internal nodes
  std::set_difference(this->nodes_attached_to_local_elems.begin(), this->nodes_attached_to_local_elems.end(),
                      this->border_node_ids.begin(), this->border_node_ids.end(),
                      std::inserter(this->internal_node_ids, this->internal_node_ids.end()));

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] neighboring_processor_ids = ";
      for (std::set<unsigned>::iterator it = neighboring_processor_ids.begin();
           it != neighboring_processor_ids.end();
           ++it)
        {
          libMesh::out << *it << " ";
        }
      libMesh::out << std::endl;
    }

  // The size of the neighboring_processor_ids set should be the number of element communication maps
  this->num_elem_cmaps =
    cast_int<int>(neighboring_processor_ids.size());

  if (verbose)
    libMesh::out << "[" << this->processor_id() << "] "
                 << "Number of neighboring processor IDs="
                 << this->num_elem_cmaps
                 << std::endl;

  if (verbose)
    {
      // Print out counts of border elements for each processor
      for (proc_border_elem_sets_iterator it=this->proc_border_elem_sets.begin();
           it != this->proc_border_elem_sets.end(); ++it)
        {
          libMesh::out << "[" << this->processor_id() << "] "
                       << "Proc "
                       << it->first << " communicates "
                       << it->second.size() << " elements." << std::endl;
        }
    }

  // Store the number of internal and border elements, and the number of internal nodes,
  // to be written to the Nemesis file.
  this->num_internal_elems =
    cast_int<int>(this->internal_elem_ids.size());
  this->num_border_elems   =
    cast_int<int>(this->border_elem_ids.size());
  this->num_internal_nodes =
    cast_int<int>(this->internal_node_ids.size());

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] num_internal_nodes=" << this->num_internal_nodes << std::endl;
      libMesh::out << "[" << this->processor_id() << "] num_border_nodes=" << this->num_border_nodes << std::endl;
      libMesh::out << "[" << this->processor_id() << "] num_border_elems=" << this->num_border_elems << std::endl;
      libMesh::out << "[" << this->processor_id() << "] num_internal_elems=" << this->num_internal_elems << std::endl;
    }
}



void Nemesis_IO_Helper::compute_num_global_sidesets(const MeshBase & pmesh)
{
  // 1.) Get reference to the set of side boundary IDs
  std::set<boundary_id_type> global_side_boundary_ids
    (pmesh.get_boundary_info().get_side_boundary_ids().begin(),
     pmesh.get_boundary_info().get_side_boundary_ids().end());

  // 2.) Gather boundary side IDs from other processors
  this->comm().set_union(global_side_boundary_ids);

  // 3.) Now global_side_boundary_ids actually contains a global list of all side boundary IDs
  this->num_side_sets_global =
    cast_int<int>(global_side_boundary_ids.size());

  // 4.) Pack these sidesets into a vector so they can be written by Nemesis
  this->global_sideset_ids.clear(); // Make sure there is no leftover information
  this->global_sideset_ids.insert(this->global_sideset_ids.end(),
                                  global_side_boundary_ids.begin(),
                                  global_side_boundary_ids.end());

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] global_sideset_ids = ";
      for (std::size_t i=0; i<this->global_sideset_ids.size(); ++i)
        libMesh::out << this->global_sideset_ids[i] << ", ";
      libMesh::out << std::endl;
    }

  // We also need global counts of sides in each of the sidesets.  Again, there may be a
  // better way to do this...
  std::vector<dof_id_type> bndry_elem_list;
  std::vector<unsigned short int> bndry_side_list;
  std::vector<boundary_id_type> bndry_id_list;
  pmesh.get_boundary_info().build_side_list(bndry_elem_list, bndry_side_list, bndry_id_list);

  // Similarly to the nodes, we can't count any sides for elements which aren't local
  std::vector<dof_id_type>::iterator it_elem=bndry_elem_list.begin();
  std::vector<unsigned short>::iterator it_side=bndry_side_list.begin();
  std::vector<boundary_id_type>::iterator it_id=bndry_id_list.begin();

  // New end iterators, to be updated as we find non-local IDs
  std::vector<dof_id_type>::iterator new_bndry_elem_list_end = bndry_elem_list.end();
  std::vector<unsigned short>::iterator new_bndry_side_list_end = bndry_side_list.end();
  std::vector<boundary_id_type>::iterator new_bndry_id_list_end = bndry_id_list.end();

  for ( ; it_elem != new_bndry_elem_list_end; )
    {
      if (pmesh.elem_ref(*it_elem).processor_id() != this->processor_id())
        {
          // Back up the new end iterators to prepare for swap
          --new_bndry_elem_list_end;
          --new_bndry_side_list_end;
          --new_bndry_id_list_end;

          // Swap places, the non-local elem will now be "past-the-end"
          std::swap (*it_elem, *new_bndry_elem_list_end);
          std::swap (*it_side, *new_bndry_side_list_end);
          std::swap (*it_id, *new_bndry_id_list_end);
        }
      else // elem is local, go to next
        {
          ++it_elem;
          ++it_side;
          ++it_id;
        }
    }

  // Erase from "new" end to old end on each vector.
  bndry_elem_list.erase(new_bndry_elem_list_end, bndry_elem_list.end());
  bndry_side_list.erase(new_bndry_side_list_end, bndry_side_list.end());
  bndry_id_list.erase(new_bndry_id_list_end, bndry_id_list.end());

  this->num_global_side_counts.clear(); // Make sure we don't have any leftover information
  this->num_global_side_counts.resize(this->global_sideset_ids.size());

  // Get the count for each global sideset ID
  for (std::size_t i=0; i<global_sideset_ids.size(); ++i)
    {
      this->num_global_side_counts[i] = cast_int<int>
        (std::count(bndry_id_list.begin(),
                    bndry_id_list.end(),
                    cast_int<boundary_id_type>(this->global_sideset_ids[i])));
    }

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] num_global_side_counts = ";
      for (std::size_t i=0; i<this->num_global_side_counts.size(); ++i)
        libMesh::out << this->num_global_side_counts[i] << ", ";
      libMesh::out << std::endl;
    }

  // Finally sum up the result
  this->comm().sum(this->num_global_side_counts);

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] num_global_side_counts = ";
      for (std::size_t i=0; i<this->num_global_side_counts.size(); ++i)
        libMesh::out << this->num_global_side_counts[i] << ", ";
      libMesh::out << std::endl;
    }
}






void Nemesis_IO_Helper::compute_num_global_nodesets(const MeshBase & pmesh)
{
  std::set<boundary_id_type> local_node_boundary_ids;

  // 1.) Get reference to the set of node boundary IDs *for this processor*
  std::set<boundary_id_type> global_node_boundary_ids
    (pmesh.get_boundary_info().get_node_boundary_ids().begin(),
     pmesh.get_boundary_info().get_node_boundary_ids().end());

  // Save a copy of the local_node_boundary_ids...
  local_node_boundary_ids = global_node_boundary_ids;

  // 2.) Gather boundary node IDs from other processors
  this->comm().set_union(global_node_boundary_ids);

  // 3.) Now global_node_boundary_ids actually contains a global list of all node boundary IDs
  this->num_node_sets_global =
    cast_int<int>(global_node_boundary_ids.size());

  // 4.) Create a vector<int> from the global_node_boundary_ids set
  this->global_nodeset_ids.clear();
  this->global_nodeset_ids.insert(this->global_nodeset_ids.end(),
                                  global_node_boundary_ids.begin(),
                                  global_node_boundary_ids.end());

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] global_nodeset_ids = ";
      for (std::size_t i=0; i<global_nodeset_ids.size(); ++i)
        libMesh::out << global_nodeset_ids[i] << ", ";
      libMesh::out << std::endl;

      libMesh::out << "[" << this->processor_id() << "] local_node_boundary_ids = ";
      for (std::set<boundary_id_type>::iterator it = local_node_boundary_ids.begin();
           it != local_node_boundary_ids.end();
           ++it)
        libMesh::out << *it << ", ";
      libMesh::out << std::endl;
    }

  // 7.) We also need to know the number of nodes which is in each of the nodesets, globally.
  // There is probably a better way to do this...
  std::vector<dof_id_type> boundary_node_list;
  std::vector<boundary_id_type> boundary_node_boundary_id_list;
  pmesh.get_boundary_info().build_node_list
    (boundary_node_list, boundary_node_boundary_id_list);

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] boundary_node_list.size()="
                   << boundary_node_list.size() << std::endl;
      libMesh::out << "[" << this->processor_id() << "] (boundary_node_id, boundary_id) = ";
      for (std::size_t i=0; i<boundary_node_list.size(); ++i)
        {
          libMesh::out << "(" << boundary_node_list[i] << ", " << boundary_node_boundary_id_list[i] << ") ";
        }
      libMesh::out << std::endl;
    }

  // Now get the global information.  In this case, we only want to count boundary
  // information for nodes *owned* by this processor, so there are no duplicates.

  // Make sure we don't have any left over information
  this->num_global_node_counts.clear();
  this->num_global_node_counts.resize(this->global_nodeset_ids.size());

  // Unfortunately, we can't just count up all occurrences of a given id,
  // that would give us duplicate entries when we do the parallel summation.
  // So instead, only count entries for nodes owned by this processor.
  // Start by getting rid of all non-local node entries from the vectors.
  std::vector<dof_id_type>::iterator it_node=boundary_node_list.begin();
  std::vector<boundary_id_type>::iterator it_id=boundary_node_boundary_id_list.begin();

  // New end iterators, to be updated as we find non-local IDs
  std::vector<dof_id_type>::iterator new_node_list_end = boundary_node_list.end();
  std::vector<boundary_id_type>::iterator new_boundary_id_list_end = boundary_node_boundary_id_list.end();
  for ( ; it_node != new_node_list_end; )
    {
      if (pmesh.node_ptr( *it_node )->processor_id() != this->processor_id())
        {
          // Back up the new end iterators to prepare for swap
          --new_node_list_end;
          --new_boundary_id_list_end;

          // Swap places, the non-local node will now be "past-the-end"
          std::swap (*it_node, *new_node_list_end);
          std::swap (*it_id, *new_boundary_id_list_end);
        }
      else // node is local, go to next
        {
          ++it_node;
          ++it_id;
        }
    }

  // Erase from "new" end to old end on each vector.
  boundary_node_list.erase(new_node_list_end, boundary_node_list.end());
  boundary_node_boundary_id_list.erase(new_boundary_id_list_end, boundary_node_boundary_id_list.end());

  // Now we can do the local count for each ID...
  for (std::size_t i=0; i<global_nodeset_ids.size(); ++i)
    {
      this->num_global_node_counts[i] = cast_int<int>
        (std::count(boundary_node_boundary_id_list.begin(),
                    boundary_node_boundary_id_list.end(),
                    cast_int<boundary_id_type>(this->global_nodeset_ids[i])));
    }

  // And finally we can sum them up
  this->comm().sum(this->num_global_node_counts);

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] num_global_node_counts = ";
      for (std::size_t i=0; i<num_global_node_counts.size(); ++i)
        libMesh::out << num_global_node_counts[i] << ", ";
      libMesh::out << std::endl;
    }
}




void Nemesis_IO_Helper::compute_num_global_elem_blocks(const MeshBase & pmesh)
{
  // 1.) Loop over active local elements, build up set of subdomain IDs.
  std::set<subdomain_id_type> global_subdomain_ids;

  // This map keeps track of the number of elements in each subdomain over all processors
  std::map<subdomain_id_type, unsigned> global_subdomain_counts;

  MeshBase::const_element_iterator elem_it = pmesh.active_local_elements_begin();
  MeshBase::const_element_iterator elem_end = pmesh.active_local_elements_end();

  for (; elem_it != elem_end; ++elem_it)
    {
      const Elem * elem = *elem_it;

      subdomain_id_type cur_subdomain = elem->subdomain_id();

      /*
      // We can't have a zero subdomain ID in Exodus (for some reason?)
      // so map zero subdomains to a max value...
      if (cur_subdomain == 0)
      cur_subdomain = std::numeric_limits<subdomain_id_type>::max();
      */

      global_subdomain_ids.insert(cur_subdomain);

      // Increment the count of elements in this subdomain
      global_subdomain_counts[cur_subdomain]++;
    }

  // We're next going to this->comm().sum the subdomain counts, so save the local counts
  this->local_subdomain_counts = global_subdomain_counts;

  {
    // 2.) Copy local subdomain IDs into a vector for communication
    std::vector<subdomain_id_type> global_subdomain_ids_vector(global_subdomain_ids.begin(),
                                                               global_subdomain_ids.end());

    // 3.) Gather them into an enlarged vector
    this->comm().allgather(global_subdomain_ids_vector);

    // 4.) Insert any new IDs into the set (any duplicates will be dropped)
    global_subdomain_ids.insert(global_subdomain_ids_vector.begin(),
                                global_subdomain_ids_vector.end());
  }

  // 5.) Now global_subdomain_ids actually contains a global list of all subdomain IDs
  this->num_elem_blks_global =
    cast_int<int>(global_subdomain_ids.size());

  // Print the number of elements found locally in each subdomain
  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] ";
      for (std::map<subdomain_id_type, unsigned>::iterator it=global_subdomain_counts.begin();
           it != global_subdomain_counts.end();
           ++it)
        {
          libMesh::out << "ID: "
                       << static_cast<unsigned>(it->first)
                       << ", Count: " << it->second << ", ";
        }
      libMesh::out << std::endl;
    }

  // 6.) this->comm().sum up the number of elements in each block.  We know the global
  // subdomain IDs, so pack them into a vector one by one.  Use a vector of int since
  // that is what Nemesis wants
  this->global_elem_blk_cnts.resize(global_subdomain_ids.size());

  unsigned cnt=0;
  for (std::set<subdomain_id_type>::iterator it=global_subdomain_ids.begin();
       it != global_subdomain_ids.end(); ++it)
    {
      // Find the entry in the local map, note: if not found, will be created with 0 default value, which is OK...
      this->global_elem_blk_cnts[cnt++] = global_subdomain_counts[*it];
    }

  // Sum up subdomain counts from all processors
  this->comm().sum(this->global_elem_blk_cnts);

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] global_elem_blk_cnts = ";
      for (std::size_t i=0; i<this->global_elem_blk_cnts.size(); ++i)
        libMesh::out << this->global_elem_blk_cnts[i] << ", ";
      libMesh::out << std::endl;
    }

  // 7.) Create a vector<int> from the global_subdomain_ids set, for passing to Nemesis
  this->global_elem_blk_ids.clear();
  this->global_elem_blk_ids.insert(this->global_elem_blk_ids.end(), // pos
                                   global_subdomain_ids.begin(),
                                   global_subdomain_ids.end());

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] global_elem_blk_ids = ";
      for (std::size_t i=0; i<this->global_elem_blk_ids.size(); ++i)
        libMesh::out << this->global_elem_blk_ids[i] << ", ";
      libMesh::out << std::endl;
    }


  // 8.) We will call put_eb_info_global later, it must be called after this->put_init_global().
}




void Nemesis_IO_Helper::build_element_and_node_maps(const MeshBase & pmesh)
{
  // If we don't have any local subdomains, it had better be because
  // we don't have any local elements
#ifdef DEBUG
  if (local_subdomain_counts.empty())
    {
      libmesh_assert(pmesh.active_local_elements_begin() ==
                     pmesh.active_local_elements_end());
      libmesh_assert(this->nodes_attached_to_local_elems.empty());
    }
#endif

  // Elements have to be numbered contiguously based on what block
  // number they are in.  Therefore we have to do a bit of work to get
  // the block (ie subdomain) numbers first and store them off as
  // block_ids.

  // Make sure there is no leftover information in the subdomain_map, and reserve
  // enough space to store the elements we need.
  this->subdomain_map.clear();
  for (std::map<subdomain_id_type, unsigned>::iterator it=this->local_subdomain_counts.begin();
       it != this->local_subdomain_counts.end();
       ++it)
    {
      subdomain_id_type cur_subdomain = it->first;

      /*
      // We can't have a zero subdomain ID in Exodus (for some reason?)
      // so map zero subdomains to a max value...
      if (cur_subdomain == 0)
      cur_subdomain = std::numeric_limits<subdomain_id_type>::max();
      */

      if (verbose)
        {
          libMesh::out << "[" << this->processor_id() << "] "
                       << "local_subdomain_counts [" << static_cast<unsigned>(cur_subdomain) << "]= "
                       << it->second
                       << std::endl;
        }

      // *it.first is the subodmain ID, *it.second is the number of elements it contains
      this->subdomain_map[ cur_subdomain ].reserve( it->second );
    }


  // First loop over the elements to figure out which elements are in which subdomain
  MeshBase::const_element_iterator elem_it = pmesh.active_local_elements_begin();
  MeshBase::const_element_iterator elem_end = pmesh.active_local_elements_end();

  for (; elem_it != elem_end; ++elem_it)
    {
      const Elem * elem = *elem_it;

      // Grab the nodes while we're here.
      for (unsigned int n=0; n<elem->n_nodes(); ++n)
        this->nodes_attached_to_local_elems.insert( elem->node_id(n) );

      subdomain_id_type cur_subdomain = elem->subdomain_id();

      this->subdomain_map[cur_subdomain].push_back
        (cast_int<unsigned>(elem->id()));
    }

  // Set num_nodes which is used by exodusII_io_helper
  this->num_nodes =
    cast_int<int>(this->nodes_attached_to_local_elems.size());

  // Now come up with a 1-based numbering for these nodes
  this->exodus_node_num_to_libmesh.clear(); // Make sure it's empty
  this->exodus_node_num_to_libmesh.reserve(this->nodes_attached_to_local_elems.size());

  // Also make sure there's no leftover information in the map which goes the
  // other direction.
  this->libmesh_node_num_to_exodus.clear();

  // Set the map for nodes
  for (std::set<int>::iterator it = this->nodes_attached_to_local_elems.begin();
       it != this->nodes_attached_to_local_elems.end();
       ++it)
    {
      // I.e. given exodus_node_id,
      // exodus_node_num_to_libmesh[ exodus_node_id ] returns the libmesh ID for that node.
      // Note that even though most of Exodus is 1-based, this code will map an Exodus ID of
      // zero to some libmesh node ID.  Is that a problem?
      this->exodus_node_num_to_libmesh.push_back(*it);

      // Likewise, given libmesh_node_id,
      // libmesh_node_num_to_exodus[ libmesh_node_id ] returns the *Exodus* ID for that node.
      // Unlike the exodus_node_num_to_libmesh vector above, this one is a std::map
      this->libmesh_node_num_to_exodus[*it] =
        cast_int<int>(this->exodus_node_num_to_libmesh.size()); // should never be zero...
    }

  // Now we're going to loop over the subdomain map and build a few things right
  // now that we'll use later.

  // First make sure our data structures don't have any leftover data...
  this->exodus_elem_num_to_libmesh.clear();
  this->block_ids.clear();
  this->libmesh_elem_num_to_exodus.clear();

  // Now loop over each subdomain and get a unique numbering for the elements
  for (std::map<subdomain_id_type, std::vector<unsigned int> >::iterator it = this->subdomain_map.begin();
       it != this->subdomain_map.end();
       ++it)
    {
      block_ids.push_back(it->first);

      // Vector of element IDs for this subdomain
      std::vector<unsigned int> & elem_ids_this_subdomain = it->second;

      // The code below assumes this subdomain block is not empty, make sure that's the case!
      if (elem_ids_this_subdomain.size() == 0)
        libmesh_error_msg("Error, no element IDs found in subdomain " << it->first);

      ExodusII_IO_Helper::ElementMaps em;

      // Use the first element in this block to get representative information.
      // Note that Exodus assumes all elements in a block are of the same type!
      // We are using that same assumption here!
      const ExodusII_IO_Helper::Conversion conv = em.assign_conversion
        (pmesh.elem_ref(elem_ids_this_subdomain[0]).type());
      this->num_nodes_per_elem =
        pmesh.elem_ref(elem_ids_this_subdomain[0]).n_nodes();

      // Get a reference to the connectivity vector for this subdomain.  This vector
      // is most likely empty, we are going to fill it up now.
      std::vector<int> & current_block_connectivity = this->block_id_to_elem_connectivity[it->first];

      // Just in case it's not already empty...
      current_block_connectivity.clear();
      current_block_connectivity.resize(elem_ids_this_subdomain.size() * this->num_nodes_per_elem);

      for (std::size_t i=0; i<elem_ids_this_subdomain.size(); i++)
        {
          unsigned int elem_id = elem_ids_this_subdomain[i];

          // Set the number map for elements
          this->exodus_elem_num_to_libmesh.push_back(elem_id);
          this->libmesh_elem_num_to_exodus[elem_id] =
            cast_int<int>(this->exodus_elem_num_to_libmesh.size());

          const Elem & elem = pmesh.elem_ref(elem_id);

          // Exodus/Nemesis want every block to have the same element type
          // libmesh_assert_equal_to (elem->type(), conv.get_canonical_type());

          // But we can get away with writing e.g. HEX8 and INFHEX8 in
          // the same block...
          libmesh_assert_equal_to (elem.n_nodes(), Elem::build(conv.get_canonical_type(), libmesh_nullptr)->n_nodes());

          for (unsigned int j=0; j < static_cast<unsigned int>(this->num_nodes_per_elem); j++)
            {
              const unsigned int connect_index   = (i*this->num_nodes_per_elem)+j;
              const unsigned int elem_node_index = conv.get_node_map(j);

              std::map<int, int>::iterator node_it = libmesh_node_num_to_exodus.find(elem.node_id(elem_node_index));
              if (node_it == libmesh_node_num_to_exodus.end())
                libmesh_error_msg("Node number " << elem.node_id(elem_node_index) << " not found in libmesh_node_num_to_exodus map.");

              current_block_connectivity[connect_index] = node_it->second;
            }
        } // End loop over elems in this subdomain
    } // end loop over subdomain_map
}





void Nemesis_IO_Helper::compute_border_node_ids(const MeshBase & pmesh)
{
  // The set which will eventually contain the IDs of "border nodes".  These are nodes
  // that lie on the boundary between one or more processors.
  //std::set<unsigned> border_node_ids;

  // map from processor ID to set of nodes which elements from this processor "touch",
  // that is,
  // proc_nodes_touched[p] = (set all node IDs found in elements owned by processor p)
  std::map<unsigned, std::set<unsigned> > proc_nodes_touched;


  // We are going to create a lot of intermediate data structures here, so make sure
  // as many as possible all cleaned up by creating scope!
  {
    // Loop over active (not just active local) elements, make sets of node IDs for each
    // processor which has an element that "touches" a node.
    {
      MeshBase::const_element_iterator elem_it = pmesh.active_elements_begin();
      MeshBase::const_element_iterator elem_end = pmesh.active_elements_end();

      for (; elem_it != elem_end; ++elem_it)
        {
          const Elem * elem = *elem_it;

          // Get reference to the set for this processor.  If it does not exist
          // it will be created.
          std::set<unsigned> & set_p = proc_nodes_touched[ elem->processor_id() ];

          // Insert all nodes touched by this element into the set
          for (unsigned int node=0; node<elem->n_nodes(); ++node)
            set_p.insert(elem->node_id(node));
        }
    }

    // The number of node communication maps is the number of other processors
    // with which we share nodes. (I think.) This is just the size of the map we just
    // created, minus 1.
    this->num_node_cmaps =
      cast_int<int>(proc_nodes_touched.size() - 1);

    // If we've got no elements on this processor and haven't touched
    // any nodes, however, then that's 0 other processors with which
    // we share nodes, not -1.
    if (this->num_node_cmaps == -1)
      {
        libmesh_assert (pmesh.active_elements_begin() == pmesh.active_elements_end());
        this->num_node_cmaps = 0;
      }

    // We can't be connecting to more processors than exist outside
    // ourselves
    libmesh_assert_less (static_cast<unsigned>(this->num_node_cmaps), this->n_processors());

    if (verbose)
      {
        libMesh::out << "[" << this->processor_id()
                     << "] proc_nodes_touched contains "
                     << proc_nodes_touched.size()
                     << " sets of nodes."
                     << std::endl;

        for (proc_nodes_touched_iterator it = proc_nodes_touched.begin();
             it != proc_nodes_touched.end();
             ++it)
          {
            libMesh::out << "[" << this->processor_id()
                         << "] proc_nodes_touched[" << it->first << "] has "
                         << it->second.size()
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
        if (it->first == this->processor_id())
          continue;

        // Otherwise, compute intersection with other processor and ourself
        std::set<unsigned> & my_set = proc_nodes_touched[this->processor_id()];
        std::set<unsigned> & other_set = it->second;
        std::set<unsigned> & result_set = this->proc_nodes_touched_intersections[ it->first ]; // created if does not exist

        std::set_intersection(my_set.begin(), my_set.end(),
                              other_set.begin(), other_set.end(),
                              std::inserter(result_set, result_set.end()));
      }

    if (verbose)
      {
        for (proc_nodes_touched_iterator it = this->proc_nodes_touched_intersections.begin();
             it != this->proc_nodes_touched_intersections.end();
             ++it)
          {
            libMesh::out << "[" << this->processor_id()
                         << "] this->proc_nodes_touched_intersections[" << it->first << "] has "
                         << it->second.size()
                         << " entries."
                         << std::endl;
          }
      }

    // Compute the set_union of all the preceding intersections.  This will be the set of
    // border node IDs for this processor.
    for (proc_nodes_touched_iterator it = this->proc_nodes_touched_intersections.begin();
         it != this->proc_nodes_touched_intersections.end();
         ++it)
      {
        std::set<unsigned> & other_set = it->second;
        std::set<unsigned> intermediate_result; // Don't think we can insert into one of the sets we're unioning...

        std::set_union(this->border_node_ids.begin(), this->border_node_ids.end(),
                       other_set.begin(), other_set.end(),
                       std::inserter(intermediate_result, intermediate_result.end()));

        // Swap our intermediate result into the final set
        this->border_node_ids.swap(intermediate_result);
      }

    if (verbose)
      {
        libMesh::out << "[" << this->processor_id()
                     << "] border_node_ids.size()=" << this->border_node_ids.size()
                     << std::endl;
      }
  } // end scope for border node ID creation

  // Store the number of border node IDs to be written to Nemesis file
  this->num_border_nodes = cast_int<int>(this->border_node_ids.size());
}





void Nemesis_IO_Helper::write_nodesets(const MeshBase & mesh)
{
  // Write the nodesets.  In Nemesis, the idea is to "create space" for the global
  // set of boundary nodesets, but to only write node IDs which are local to the current
  // processor.  This is what is done in Nemesis files created by the "loadbal" script.

  // Store a map of vectors for boundary node IDs on this processor.
  // Use a vector of int here so it can be passed directly to Exodus.
  std::map<boundary_id_type, std::vector<int> > local_node_boundary_id_lists;
  typedef std::map<boundary_id_type, std::vector<int> >::iterator local_node_boundary_id_lists_iterator;

  // FIXME: We should build this list only one time!!  We already built it above, but we
  // did not have the libmesh to exodus node mapping at that time... for now we'll just
  // build it here again, hopefully it's small relative to the size of the entire mesh.
  std::vector<dof_id_type> boundary_node_list;
  std::vector<boundary_id_type> boundary_node_boundary_id_list;
  mesh.get_boundary_info().build_node_list
    (boundary_node_list, boundary_node_boundary_id_list);

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] boundary_node_list.size()="
                   << boundary_node_list.size() << std::endl;
      libMesh::out << "[" << this->processor_id() << "] (boundary_node_id, boundary_id) = ";
      for (std::size_t i=0; i<boundary_node_list.size(); ++i)
        {
          libMesh::out << "(" << boundary_node_list[i] << ", " << boundary_node_boundary_id_list[i] << ") ";
        }
      libMesh::out << std::endl;
    }

  // For each node in the node list, add it to the vector of node IDs for that
  // set for the local processor.  This will be used later when writing Exodus
  // nodesets.
  for (std::size_t i=0; i<boundary_node_list.size(); ++i)
    {
      // Don't try to grab a reference to the vector unless the current node is attached
      // to a local element.  Otherwise, another processor will be responsible for writing it in its nodeset.
      std::map<int, int>::iterator it = this->libmesh_node_num_to_exodus.find( boundary_node_list[i] );

      if ( it != this->libmesh_node_num_to_exodus.end() )
        {
          // Get reference to the vector where this node ID will be inserted.  If it
          // doesn't yet exist, this will create it.
          std::vector<int> & current_id_set = local_node_boundary_id_lists[ boundary_node_boundary_id_list[i] ];

          // Push back Exodus-mapped node ID for this set
          // TODO: reserve space in these vectors somehow.
          current_id_set.push_back( it->second );
        }
    }

  // See what we got
  if (verbose)
    {
      for (std::map<boundary_id_type, std::vector<int> >::iterator it = local_node_boundary_id_lists.begin();
           it != local_node_boundary_id_lists.end();
           ++it)
        {
          libMesh::out << "[" << this->processor_id() << "] ID: " << it->first << ", ";

          std::vector<int> & current_id_set = it->second;

          // Libmesh node ID (Exodus Node ID)
          for (std::size_t j=0; j<current_id_set.size(); ++j)
            libMesh::out << current_id_set[j]
                         << ", ";

          libMesh::out << std::endl;
        }
    }

  // Loop over *global* nodeset IDs, call the Exodus API.  Note that some nodesets may be empty
  // for a given processor.
  if (global_nodeset_ids.size() > 0) {
  NamesData names_table(global_nodeset_ids.size(), MAX_STR_LENGTH);

  for (std::size_t i=0; i<this->global_nodeset_ids.size(); ++i)
    {
      const std::string & current_ns_name =
        mesh.get_boundary_info().get_nodeset_name(global_nodeset_ids[i]);

      // Store this name in a data structure that will be used to
      // write sideset names to file.
      names_table.push_back_entry(current_ns_name);

      if (verbose)
        {
          libMesh::out << "[" << this->processor_id()
                       << "] Writing out Exodus nodeset info for ID: " << global_nodeset_ids[i]
                       << ", Name: " << current_ns_name
                       << std::endl;
        }

      // Convert current global_nodeset_id into an exodus ID, which can't be zero...
      int exodus_id = global_nodeset_ids[i];

      /*
      // Exodus can't handle zero nodeset IDs (?)  Use max short here since
      // when libmesh reads it back in, it will want to store it as a short...
      if (exodus_id==0)
      exodus_id = std::numeric_limits<short>::max();
      */

      // Try to find this boundary ID in the local list we created
      local_node_boundary_id_lists_iterator it =
        local_node_boundary_id_lists.find (cast_int<boundary_id_type>(this->global_nodeset_ids[i]));

      // No nodes found for this boundary ID on this processor
      if (it == local_node_boundary_id_lists.end())
        {
          if (verbose)
            libMesh::out << "[" << this->processor_id()
                         << "] No nodeset data for ID: " << global_nodeset_ids[i]
                         << " on this processor." << std::endl;

          // Call the Exodus interface to write the parameters of this node set
          this->ex_err = exII::ex_put_node_set_param(this->ex_id,
                                                     exodus_id,
                                                     0, /* No nodes for this ID */
                                                     0  /* No distribution factors */);
          EX_CHECK_ERR(this->ex_err, "Error writing nodeset parameters in Nemesis");

        }
      else // Boundary ID *was* found in list
        {
          // Get reference to the vector of node IDs
          std::vector<int> & current_nodeset_ids = it->second;

          // Call the Exodus interface to write the parameters of this node set
          this->ex_err = exII::ex_put_node_set_param(this->ex_id,
                                                     exodus_id,
                                                     current_nodeset_ids.size(),
                                                     0  /* No distribution factors */);

          EX_CHECK_ERR(this->ex_err, "Error writing nodeset parameters in Nemesis");

          // Call Exodus interface to write the actual node IDs for this boundary ID
          this->ex_err = exII::ex_put_node_set(this->ex_id,
                                               exodus_id,
                                               &current_nodeset_ids[0]);

          EX_CHECK_ERR(this->ex_err, "Error writing nodesets in Nemesis");

        }
    } // end loop over global nodeset IDs

  // Write out the nodeset names
  ex_err = exII::ex_put_names(ex_id,
                              exII::EX_NODE_SET,
                              names_table.get_char_star_star());
  EX_CHECK_ERR(ex_err, "Error writing nodeset names");
  } // end for loop over global nodeset IDs
}




void Nemesis_IO_Helper::write_sidesets(const MeshBase & mesh)
{
  // Write the sidesets.  In Nemesis, the idea is to "create space" for the global
  // set of boundary sidesets, but to only write sideset IDs which are local to the current
  // processor.  This is what is done in Nemesis files created by the "loadbal" script.
  // See also: ExodusII_IO_Helper::write_sidesets()...


  // Store a map of vectors for boundary side IDs on this processor.
  // Use a vector of int here so it can be passed directly to Exodus.
  std::map<boundary_id_type, std::vector<int> > local_elem_boundary_id_lists;
  std::map<boundary_id_type, std::vector<int> > local_elem_boundary_id_side_lists;
  typedef std::map<boundary_id_type, std::vector<int> >::iterator local_elem_boundary_id_lists_iterator;

  ExodusII_IO_Helper::ElementMaps em;

  // FIXME: We already built this list once, we should reuse that information!
  std::vector< dof_id_type > bndry_elem_list;
  std::vector< unsigned short int > bndry_side_list;
  std::vector< boundary_id_type > bndry_id_list;

  mesh.get_boundary_info().build_side_list
    (bndry_elem_list, bndry_side_list, bndry_id_list);

  // Integer looping, skipping non-local elements
  for (std::size_t i=0; i<bndry_elem_list.size(); ++i)
    {
      // Get pointer to current Elem
      const Elem * elem = mesh.elem_ptr(bndry_elem_list[i]);

      std::vector<const Elem *> family;
#ifdef LIBMESH_ENABLE_AMR
      // We need to build up active elements if AMR is enabled and add
      // them to the exodus sidesets instead of the potentially inactive "parent" elements
      // Technically we don't need to "reset" the tree since the vector was just created.
      elem->active_family_tree_by_side(family, bndry_side_list[i], /*reset tree=*/false);
#else
      // If AMR is not even enabled, just push back the element itself
      family.push_back( elem );
#endif

      // Loop over all the elements in the family tree, store their converted IDs
      // and side IDs to the map's vectors.  TODO: Somehow reserve enough space for these
      // push_back's...
      for (std::size_t j=0; j<family.size(); ++j)
        {
          const dof_id_type f_id = family[j]->id();
          const Elem & f = mesh.elem_ref(f_id);

          // If element is local, process it
          if (f.processor_id() == this->processor_id())
            {
              const ExodusII_IO_Helper::Conversion conv = em.assign_conversion(f.type());

              // Use the libmesh to exodus datastructure map to get the proper sideset IDs
              // The datastructure contains the "collapsed" contiguous ids.
              //
              // We know the parent element is local, but let's be absolutely sure that all the children have been
              // actually mapped to Exodus IDs before we blindly try to add them...
              std::map<int,int>::iterator it = this->libmesh_elem_num_to_exodus.find( f_id );
              if (it != this->libmesh_elem_num_to_exodus.end())
                {
                  local_elem_boundary_id_lists[ bndry_id_list[i] ].push_back( it->second );
                  local_elem_boundary_id_side_lists[ bndry_id_list[i] ].push_back(conv.get_inverse_side_map( bndry_side_list[i] ));
                }
              else
                libmesh_error_msg("Error, no Exodus mapping for Elem " \
                                  << f_id                              \
                                  << " on processor "                 \
                                  << this->processor_id());
            }
        }
    }


  // Loop over *global* sideset IDs, call the Exodus API.  Note that some sidesets may be empty
  // for a given processor.
  if (global_sideset_ids.size() > 0) {
  NamesData names_table(global_sideset_ids.size(), MAX_STR_LENGTH);

  for (std::size_t i=0; i<this->global_sideset_ids.size(); ++i)
    {
      const std::string & current_ss_name =
        mesh.get_boundary_info().get_sideset_name(global_sideset_ids[i]);

      // Store this name in a data structure that will be used to
      // write sideset names to file.
      names_table.push_back_entry(current_ss_name);

      if (verbose)
        {
          libMesh::out << "[" << this->processor_id()
                       << "] Writing out Exodus sideset info for ID: " << global_sideset_ids[i]
                       << ", Name: " << current_ss_name
                       << std::endl;
        }

      // Convert current global_sideset_id into an exodus ID, which can't be zero...
      int exodus_id = global_sideset_ids[i];

      /*
      // Exodus can't handle zero sideset IDs (?)  Use max short here since
      // when libmesh reads it back in, it will want to store it as a short...
      if (exodus_id==0)
      exodus_id = std::numeric_limits<short>::max();
      */

      // Try to find this boundary ID in the local list we created
      local_elem_boundary_id_lists_iterator it =
        local_elem_boundary_id_lists.find (cast_int<boundary_id_type>(this->global_sideset_ids[i]));

      // No sides found for this boundary ID on this processor
      if (it == local_elem_boundary_id_lists.end())
        {
          if (verbose)
            libMesh::out << "[" << this->processor_id()
                         << "] No sideset data for ID: " << global_sideset_ids[i]
                         << " on this processor." << std::endl;

          // Call the Exodus interface to write the parameters of this side set
          this->ex_err = exII::ex_put_side_set_param(this->ex_id,
                                                     exodus_id,
                                                     0, /* No sides for this ID */
                                                     0  /* No distribution factors */);
          EX_CHECK_ERR(this->ex_err, "Error writing sideset parameters in Nemesis");

        }
      else // Boundary ID *was* found in list
        {
          // Get iterator to sides vector as well
          local_elem_boundary_id_lists_iterator it_sides =
            local_elem_boundary_id_side_lists.find (cast_int<boundary_id_type>(this->global_sideset_ids[i]));

          libmesh_assert (it_sides != local_elem_boundary_id_side_lists.end());

          // Get reference to the vector of elem IDs
          std::vector<int> & current_sideset_elem_ids = it->second;

          // Get reference to the vector of side IDs
          std::vector<int> & current_sideset_side_ids = (*it_sides).second;

          // Call the Exodus interface to write the parameters of this side set
          this->ex_err = exII::ex_put_side_set_param(this->ex_id,
                                                     exodus_id,
                                                     current_sideset_elem_ids.size(),
                                                     0  /* No distribution factors */);

          EX_CHECK_ERR(this->ex_err, "Error writing sideset parameters in Nemesis");

          // Call Exodus interface to write the actual side IDs for this boundary ID
          this->ex_err = exII::ex_put_side_set(this->ex_id,
                                               exodus_id,
                                               &current_sideset_elem_ids[0],
                                               &current_sideset_side_ids[0]);

          EX_CHECK_ERR(this->ex_err, "Error writing sidesets in Nemesis");
        }
    } // end for loop over global sideset IDs

  // Write sideset names to file.  Some of these may be blank strings
  // if the current processor didn't have all the sideset names for
  // any reason...
  ex_err = exII::ex_put_names(this->ex_id,
                              exII::EX_SIDE_SET,
                              names_table.get_char_star_star());
  EX_CHECK_ERR(ex_err, "Error writing sideset names");

  } // end if (global_sideset_ids.size() > 0)
}



void Nemesis_IO_Helper::write_nodal_coordinates(const MeshBase & mesh, bool /*use_discontinuous*/)
{
  // Make sure that the reference passed in is really a DistributedMesh
  // const DistributedMesh & pmesh = cast_ref<const DistributedMesh &>(mesh);

  unsigned local_num_nodes =
    cast_int<unsigned int>(this->exodus_node_num_to_libmesh.size());

  x.resize(local_num_nodes);
  y.resize(local_num_nodes);
  z.resize(local_num_nodes);

  // Just loop over our list outputing the nodes the way we built the map
  for (unsigned int i=0; i<local_num_nodes; ++i)
    {
      const Point & pt = mesh.point(this->exodus_node_num_to_libmesh[i]);
      x[i]=pt(0);
      y[i]=pt(1);
      z[i]=pt(2);
    }

  if (local_num_nodes)
    {
      if (_single_precision)
        {
          std::vector<float>
            x_single(x.begin(), x.end()),
            y_single(y.begin(), y.end()),
            z_single(z.begin(), z.end());

          ex_err = exII::ex_put_coord(ex_id, &x_single[0], &y_single[0], &z_single[0]);
        }
      else
        {
          // Call Exodus API to write nodal coordinates...
          ex_err = exII::ex_put_coord(ex_id, &x[0], &y[0], &z[0]);
        }
      EX_CHECK_ERR(ex_err, "Error writing node coordinates");

      // And write the nodal map we created for them
      ex_err = exII::ex_put_node_num_map(ex_id, &(this->exodus_node_num_to_libmesh[0]));
      EX_CHECK_ERR(ex_err, "Error writing node num map");
    }
  else // Does the Exodus API want us to write empty nodal coordinates?
    {
      ex_err = exII::ex_put_coord(ex_id, libmesh_nullptr, libmesh_nullptr, libmesh_nullptr);
      EX_CHECK_ERR(ex_err, "Error writing empty node coordinates");

      ex_err = exII::ex_put_node_num_map(ex_id, libmesh_nullptr);
      EX_CHECK_ERR(ex_err, "Error writing empty node num map");
    }
}





void Nemesis_IO_Helper::write_elements(const MeshBase & mesh, bool /*use_discontinuous*/)
{
  // Only write elements if there are elements blocks available.
  if (this->num_elem_blks_global > 0) {

  // Data structure to store element block names that will be used to
  // write the element block names to file.
  NamesData names_table(this->num_elem_blks_global, MAX_STR_LENGTH);

  // Loop over all blocks, even if we don't have elements in each block.
  // If we don't have elements we need to write out a 0 for that block...
  for (unsigned int i=0; i<static_cast<unsigned>(this->num_elem_blks_global); ++i)
    {
      // Even if there are no elements for this block on the current
      // processor, we still want to write its name to file, if
      // possible. MeshBase::subdomain_name() will just return an
      // empty string if there is no name associated with the current
      // block.
      names_table.push_back_entry(mesh.subdomain_name(this->global_elem_blk_ids[i]));

      // Search for the current global block ID in the map
      std::map<int, std::vector<int> >::iterator it =
        this->block_id_to_elem_connectivity.find( this->global_elem_blk_ids[i] );

      // If not found, write a zero to file....
      if (it == this->block_id_to_elem_connectivity.end())
        {
          this->ex_err = exII::ex_put_elem_block(this->ex_id,
                                                 this->global_elem_blk_ids[i],
                                                 "Empty",
                                                 0, /* n. elements in this block */
                                                 0, /* n. nodes per element */
                                                 0);  /* number of attributes per element */

          EX_CHECK_ERR(this->ex_err, "Error writing element block from Nemesis.");
        }

      // Otherwise, write the actual block information and connectivity to file
      else
        {
          subdomain_id_type block =
            cast_int<subdomain_id_type>(it->first);
          std::vector<int> & this_block_connectivity = it->second;
          std::vector<unsigned int> & elements_in_this_block = subdomain_map[block];

          ExodusII_IO_Helper::ElementMaps em;

          //Use the first element in this block to get representative information.
          //Note that Exodus assumes all elements in a block are of the same type!
          //We are using that same assumption here!
          const ExodusII_IO_Helper::Conversion conv =
            em.assign_conversion(mesh.elem_ref(elements_in_this_block[0]).type());

          this->num_nodes_per_elem =
            mesh.elem_ref(elements_in_this_block[0]).n_nodes();

          ex_err = exII::ex_put_elem_block(ex_id,
                                           block,
                                           conv.exodus_elem_type().c_str(),
                                           elements_in_this_block.size(),
                                           num_nodes_per_elem,
                                           0);
          EX_CHECK_ERR(ex_err, "Error writing element block from Nemesis.");

          ex_err = exII::ex_put_elem_conn(ex_id,
                                          block,
                                          &this_block_connectivity[0]);
          EX_CHECK_ERR(ex_err, "Error writing element connectivities from Nemesis.");
        }
    } // end loop over global block IDs

  // Only call this once, not in the loop above!
  ex_err = exII::ex_put_elem_num_map(ex_id,
                                     exodus_elem_num_to_libmesh.empty() ? libmesh_nullptr : &exodus_elem_num_to_libmesh[0]);
  EX_CHECK_ERR(ex_err, "Error writing element map");

  // Write the element block names to file.
  ex_err = exII::ex_put_names(ex_id, exII::EX_ELEM_BLOCK, names_table.get_char_star_star());
  EX_CHECK_ERR(ex_err, "Error writing element block names");
  } // end if (this->num_elem_blks_global > 0)
}





void Nemesis_IO_Helper::write_nodal_solution(const std::vector<Number> & values,
                                             const std::vector<std::string> & names,
                                             int timestep)
{
  int num_vars = cast_int<int>(names.size());
  //int num_values = values.size(); // Not used?

  for (int c=0; c<num_vars; c++)
    {
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
      std::vector<Real> real_parts(num_nodes);
      std::vector<Real> imag_parts(num_nodes);
      std::vector<Real> magnitudes(num_nodes);

      for (int i=0; i<num_nodes; ++i)
        {
          Number value = values[this->exodus_node_num_to_libmesh[i]*num_vars + c];
          real_parts[i] = value.real();
          imag_parts[i] = value.imag();
          magnitudes[i] = std::abs(value);
        }
      write_nodal_values(3*c+1,real_parts,timestep);
      write_nodal_values(3*c+2,imag_parts,timestep);
      write_nodal_values(3*c+3,magnitudes,timestep);
#else
      std::vector<Number> cur_soln(num_nodes);

      // Copy out this variable's solution
      for (int i=0; i<num_nodes; i++)
        cur_soln[i] = values[this->exodus_node_num_to_libmesh[i]*num_vars + c];

      write_nodal_values(c+1,cur_soln,timestep);
#endif
    }
}



void Nemesis_IO_Helper::write_nodal_solution(const NumericVector<Number> & parallel_soln,
                                             const std::vector<std::string> & names,
                                             int timestep)
{
  int num_vars = cast_int<int>(names.size());

#ifndef LIBMESH_USE_COMPLEX_NUMBERS

  for (int c=0; c<num_vars; c++)
    {
      // Fill up a std::vector with the dofs for the current variable
      std::vector<numeric_index_type> required_indices(num_nodes);

      for (int i=0; i<num_nodes; i++)
        required_indices[i] = this->exodus_node_num_to_libmesh[i]*num_vars + c;

      // Get the dof values required to write just our local part of
      // the solution vector.
      std::vector<Number> local_soln;
      parallel_soln.localize(local_soln, required_indices);

      // Call the ExodusII_IO_Helper function to write the data.
      write_nodal_values(c+1, local_soln, timestep);
    }

#else // LIBMESH_USE_COMPLEX_NUMBERS

  for (int c=0; c<num_vars; c++)
    {
      // Fill up a std::vector with the dofs for the current variable
      std::vector<numeric_index_type> required_indices(num_nodes);

      for (int i=0; i<num_nodes; i++)
        required_indices[i] = this->exodus_node_num_to_libmesh[i]*num_vars + c;

      // Get the dof values required to write just our local part of
      // the solution vector.
      std::vector<Number> local_soln;
      parallel_soln.localize(local_soln, required_indices);

      // We have the local (complex) values. Now extract the real,
      // imaginary, and magnitude values from them.
      std::vector<Real> real_parts(num_nodes);
      std::vector<Real> imag_parts(num_nodes);
      std::vector<Real> magnitudes(num_nodes);

      for (int i=0; i<num_nodes; ++i)
        {
          real_parts[i] = local_soln[i].real();
          imag_parts[i] = local_soln[i].imag();
          magnitudes[i] = std::abs(local_soln[i]);
        }

      // Write the real, imaginary, and magnitude values to file.
      write_nodal_values(3*c+1, real_parts, timestep);
      write_nodal_values(3*c+2, imag_parts, timestep);
      write_nodal_values(3*c+3, magnitudes, timestep);
    }

#endif


}




std::string Nemesis_IO_Helper::construct_nemesis_filename(const std::string & base_filename)
{
  // Build a filename for this processor.  This code is cut-n-pasted from the read function
  // and should probably be put into a separate function...
  std::ostringstream file_oss;

  // We have to be a little careful here: Nemesis left pads its file
  // numbers based on the number of processors, so for example on 10
  // processors, we'd have:
  // mesh.e.10.00
  // mesh.e.10.01
  // mesh.e.10.02
  // ...
  // mesh.e.10.09

  // And on 100 you'd have:
  // mesh.e.100.000
  // mesh.e.100.001
  // ...
  // mesh.e.128.099

  // Find the length of the highest processor ID
  file_oss << (this->n_processors());
  unsigned int field_width = cast_int<unsigned int>(file_oss.str().size());

  if (verbose)
    libMesh::out << "field_width=" << field_width << std::endl;

  file_oss.str(""); // reset the string stream
  file_oss << base_filename
           << '.' << this->n_processors()
           << '.' << std::setfill('0') << std::setw(field_width) << this->processor_id();

  // Return the resulting string
  return file_oss.str();
}

} // namespace libMesh

#endif // #if defined(LIBMESH_HAVE_NEMESIS_API) && defined(LIBMESH_HAVE_EXODUS_API)
