// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/boundary_info.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_interface.h"
#include "libmesh/int_range.h"
#include "libmesh/mesh_base.h"
#include "libmesh/node.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel.h"
#include "libmesh/utility.h"

#if defined(LIBMESH_HAVE_NEMESIS_API) && defined(LIBMESH_HAVE_EXODUS_API)

#include <libmesh/ignore_warnings.h>
namespace exII {
extern "C" {
#include "exodusII.h" // defines MAX_LINE_LENGTH, MAX_STR_LENGTH used later
}
}

// The Nemesis API header file.  Should already be
// correctly extern C'd but it doesn't hurt :)
namespace Nemesis {
extern "C" {
  // this include guard gets set by exodus, but we included it
  // in a namespace, so nemesis will not properly resolve e.g.
  // ex_entity_id in the global namespace.  undefine the guard
  // to get ne_nemesisI.h to properly include the typedefs
#  ifdef EXODUS_II_HDR
#    undef EXODUS_II_HDR
#  endif
#  ifdef EXODUSII_H
#    undef EXODUSII_H
#  endif
#  include "ne_nemesisI.h"
}
}

#include <libmesh/restore_warnings.h>

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
  num_elem_cmaps(0),
  write_complex_abs(true)
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



void Nemesis_IO_Helper::read_nodeset(int id)
{
  libmesh_assert_less (id, nodeset_ids.size());
  libmesh_assert_less (id, num_nodes_per_set.size());
  libmesh_assert_less (id, num_node_df_per_set.size());

  ex_err = exII::ex_get_set_param(ex_id,
                                  exII::EX_NODE_SET,
                                  nodeset_ids[id],
                                  &num_nodes_per_set[id],
                                  &num_node_df_per_set[id]);
  EX_CHECK_ERR(ex_err, "Error retrieving nodeset parameters.");
  message("Parameters retrieved successfully for nodeset: ", id);

  node_list.resize(num_nodes_per_set[id]);

  // Don't call ex_get_set unless there are actually nodes there to get.
  // Exodus prints an annoying warning message in DEBUG mode otherwise...
  if (num_nodes_per_set[id] > 0)
    {
      ex_err = exII::ex_get_set(ex_id,
                                exII::EX_NODE_SET,
                                nodeset_ids[id],
                                node_list.data(),
                                nullptr); // set_extra_list, ignored for node sets

      EX_CHECK_ERR(ex_err, "Error retrieving nodeset data.");
      message("Data retrieved successfully for nodeset: ", id);
    }
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
                                        global_sideset_ids.data(),
                                        num_global_side_counts.data(),
                                        num_global_side_df_counts.data());
      EX_CHECK_ERR(nemesis_err_flag, "Error reading global sideset parameters!");

      if (verbose)
        {
          libMesh::out << "[" << this->processor_id() << "] " << "Global Sideset IDs, Side Counts, and DF counts:" << std::endl;
          for (auto bn : index_range(global_sideset_ids))
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
                                        global_nodeset_ids.data(),
                                        num_global_node_counts.data(),
                                        num_global_node_df_counts.data());
      EX_CHECK_ERR(nemesis_err_flag, "Error reading global nodeset parameters!");

      if (verbose)
        {
          libMesh::out << "[" << this->processor_id() << "] " << "Global Nodeset IDs, Node Counts, and DF counts:" << std::endl;
          for (auto bn : index_range(global_nodeset_ids))
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
                                       global_elem_blk_ids.data(),
                                       global_elem_blk_cnts.data());
      EX_CHECK_ERR(nemesis_err_flag, "Error reading global element block info!");
    }

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] " << "Global Element Block IDs and Counts:" << std::endl;
      for (auto bn : index_range(global_elem_blk_ids))
        {
          libMesh::out << "  [" << this->processor_id() << "] "
                       << "global_elem_blk_ids[" << bn << "]=" << global_elem_blk_ids[bn]
                       << ", global_elem_blk_cnts[" << bn << "]=" << global_elem_blk_cnts[bn]
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
                             elem_mapi.empty() ? nullptr : elem_mapi.data(),
                             elem_mapb.empty() ? nullptr : elem_mapb.data(),
                             this->processor_id()
                             );
  EX_CHECK_ERR(nemesis_err_flag, "Error reading element maps!");


  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] elem_mapi[i] = ";
      for (auto i : make_range(num_internal_elems-1))
        libMesh::out << elem_mapi[i] << ", ";
      libMesh::out << "... " << elem_mapi.back() << std::endl;

      libMesh::out << "[" << this->processor_id() << "] elem_mapb[i] = ";
      for (auto i : make_range(std::min(10, num_border_elems-1)))
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
                             node_mapi.empty() ? nullptr : node_mapi.data(),
                             node_mapb.empty() ? nullptr : node_mapb.data(),
                             node_mape.empty() ? nullptr : node_mape.data(),
                             this->processor_id()
                             );
  EX_CHECK_ERR(nemesis_err_flag, "Error reading node maps!");

  if (verbose)
    {
      // Remark: The Exodus/Nemesis node numbering is always (?) 1-based!  This means the first interior node id will
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
                                node_cmap_ids.empty()       ? nullptr : node_cmap_ids.data(),
                                node_cmap_node_cnts.empty() ? nullptr : node_cmap_node_cnts.data(),
                                elem_cmap_ids.empty()       ? nullptr : elem_cmap_ids.data(),
                                elem_cmap_elem_cnts.empty() ? nullptr : elem_cmap_elem_cnts.data(),
                                this->processor_id());
  EX_CHECK_ERR(nemesis_err_flag, "Error reading cmap parameters!");


  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] ";
      for (auto i : index_range(node_cmap_ids))
        libMesh::out << "node_cmap_ids[" << i << "]=" << node_cmap_ids[i] << " ";
      libMesh::out << std::endl;

      libMesh::out << "[" << this->processor_id() << "] ";
      for (auto i : index_range(node_cmap_node_cnts))
        libMesh::out << "node_cmap_node_cnts[" << i << "]=" << node_cmap_node_cnts[i] << " ";
      libMesh::out << std::endl;

      libMesh::out << "[" << this->processor_id() << "] ";
      for (auto i : index_range(elem_cmap_ids))
        libMesh::out << "elem_cmap_ids[" << i << "]=" << elem_cmap_ids[i] << " ";
      libMesh::out << std::endl;

      libMesh::out << "[" << this->processor_id() << "] ";
      for (auto i : index_range(elem_cmap_elem_cnts))
        libMesh::out << "elem_cmap_elem_cnts[" << i << "]=" << elem_cmap_elem_cnts[i] << " ";
      libMesh::out << std::endl;
    }
}



void Nemesis_IO_Helper::get_node_cmap()
{
  node_cmap_node_ids.resize(num_node_cmaps);
  node_cmap_proc_ids.resize(num_node_cmaps);

  for (auto i : index_range(node_cmap_node_ids))
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
                                      node_cmap_node_ids[i].data(),
                                      node_cmap_proc_ids[i].data(),
                                      this->processor_id());
          EX_CHECK_ERR(nemesis_err_flag, "Error reading node cmap node and processor ids!");
        }

      if (verbose)
        {
          libMesh::out << "[" << this->processor_id() << "] node_cmap_node_ids[" << i << "]=";
          for (const auto & dof : node_cmap_node_ids[i])
            libMesh::out << dof << " ";
          libMesh::out << std::endl;

          // This is basically a vector, all entries of which are = node_cmap_ids[i]
          // Not sure if it's always guaranteed to be that or what...
          libMesh::out << "[" << this->processor_id() << "] node_cmap_proc_ids[" << i << "]=";
          for (const auto & dof : node_cmap_proc_ids[i])
            libMesh::out << dof << " ";
          libMesh::out << std::endl;
        }
    }
}



void Nemesis_IO_Helper::get_elem_cmap()
{
  elem_cmap_elem_ids.resize(num_elem_cmaps);
  elem_cmap_side_ids.resize(num_elem_cmaps);
  elem_cmap_proc_ids.resize(num_elem_cmaps);

  for (auto i : index_range(elem_cmap_elem_ids))
    {
      elem_cmap_elem_ids[i].resize(elem_cmap_elem_cnts[i]);
      elem_cmap_side_ids[i].resize(elem_cmap_elem_cnts[i]);
      elem_cmap_proc_ids[i].resize(elem_cmap_elem_cnts[i]);

      if (elem_cmap_elem_cnts[i] > 0)
        {
          nemesis_err_flag =
            Nemesis::ne_get_elem_cmap(ex_id,
                                      elem_cmap_ids[i],
                                      elem_cmap_elem_ids[i].data(),
                                      elem_cmap_side_ids[i].data(),
                                      elem_cmap_proc_ids[i].data(),
                                      this->processor_id());
          EX_CHECK_ERR(nemesis_err_flag, "Error reading elem cmap elem, side, and processor ids!");
        }

      if (verbose)
        {
          libMesh::out << "[" << this->processor_id() << "] elem_cmap_elem_ids[" << i << "]=";
          for (const auto & dof : elem_cmap_elem_ids[i])
            libMesh::out << dof << " ";
          libMesh::out << std::endl;

          // These must be the (local) side IDs (in the ExodusII face numbering scheme)
          // of the sides shared across processors.
          libMesh::out << "[" << this->processor_id() << "] elem_cmap_side_ids[" << i << "]=";
          for (const auto & dof : elem_cmap_side_ids[i])
            libMesh::out << dof << " ";
          libMesh::out << std::endl;

          // This is basically a vector, all entries of which are = elem_cmap_ids[i]
          // Not sure if it's always guaranteed to be that or what...
          libMesh::out << "[" << this->processor_id() << "] elem_cmap_proc_ids[" << i << "]=";
          for (const auto & dof : elem_cmap_proc_ids[i])
            libMesh::out << dof << " ";
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
                                   global_elem_blk_ids_in.data(),
                                   global_elem_blk_cnts_in.data());

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
                                        global_nodeset_ids_in.data(),
                                        num_global_node_counts_in.data(),
                                        num_global_node_df_counts_in.data());
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
                                        global_sideset_ids_in.data(),
                                        num_global_side_counts_in.data(),
                                        num_global_side_df_counts_in.data());
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
  nemesis_err_flag =
    Nemesis::ne_put_cmap_params(ex_id,
                                node_cmap_ids_in.empty() ? nullptr : node_cmap_ids_in.data(),
                                node_cmap_node_cnts_in.empty() ? nullptr : node_cmap_node_cnts_in.data(),
                                elem_cmap_ids_in.empty() ? nullptr : elem_cmap_ids_in.data(),
                                elem_cmap_elem_cnts_in.empty() ? nullptr : elem_cmap_elem_cnts_in.data(),
                                this->processor_id());

  EX_CHECK_ERR(nemesis_err_flag, "Error writing cmap parameters!");
}




void Nemesis_IO_Helper::put_node_cmap(std::vector<std::vector<int>> & node_cmap_node_ids_in,
                                      std::vector<std::vector<int>> & node_cmap_proc_ids_in)
{
  // Print to screen what we are about to print to Nemesis file
  if (verbose)
    {
      for (auto i : index_range(node_cmap_node_ids_in))
        {
          libMesh::out << "[" << this->processor_id() << "] put_node_cmap() : nodes communicated to proc "
                       << this->node_cmap_ids[i]
                       << " = ";
          for (const auto & node_id : node_cmap_node_ids_in[i])
            libMesh::out << node_id << " ";
          libMesh::out << std::endl;
        }

      for (auto i : index_range(node_cmap_node_ids_in))
        {
          libMesh::out << "[" << this->processor_id() << "] put_node_cmap() : processor IDs = ";
          for (const auto & proc_id : node_cmap_proc_ids_in[i])
            libMesh::out << proc_id << " ";
          libMesh::out << std::endl;
        }
    }

  for (auto i : index_range(node_cmap_node_ids_in))
    {
      int * node_ids_ptr = node_cmap_node_ids_in[i].empty() ?
        nullptr : node_cmap_node_ids_in[i].data();
      int * proc_ids_ptr = node_cmap_proc_ids_in[i].empty() ?
        nullptr : node_cmap_proc_ids_in[i].data();

      nemesis_err_flag =
        Nemesis::ne_put_node_cmap(ex_id, this->node_cmap_ids[i],
                                  node_ids_ptr, proc_ids_ptr,
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
                             node_mapi_in.empty() ? nullptr : node_mapi_in.data(),
                             node_mapb_in.empty() ? nullptr : node_mapb_in.data(),
                             node_mape_in.empty() ? nullptr : node_mape_in.data(),
                             this->processor_id());

  EX_CHECK_ERR(nemesis_err_flag, "Error writing Nemesis internal and border node maps to file!");
}




void Nemesis_IO_Helper::put_elem_cmap(std::vector<std::vector<int>> & elem_cmap_elem_ids_in,
                                      std::vector<std::vector<int>> & elem_cmap_side_ids_in,
                                      std::vector<std::vector<int>> & elem_cmap_proc_ids_in)
{
  for (auto i : index_range(elem_cmap_ids))
    {
      nemesis_err_flag =
        Nemesis::ne_put_elem_cmap(ex_id,
                                  this->elem_cmap_ids[i],
                                  elem_cmap_elem_ids_in[i].data(),
                                  elem_cmap_side_ids_in[i].data(),
                                  elem_cmap_proc_ids_in[i].data(),
                                  this->processor_id());

      EX_CHECK_ERR(nemesis_err_flag, "Error writing elem communication map to file!");
    }
}




void Nemesis_IO_Helper::put_elem_map(std::vector<int> & elem_mapi_in,
                                     std::vector<int> & elem_mapb_in)
{
  nemesis_err_flag =
    Nemesis::ne_put_elem_map(ex_id,
                             elem_mapi_in.empty() ? nullptr : elem_mapi_in.data(),
                             elem_mapb_in.empty() ? nullptr : elem_mapb_in.data(),
                             this->processor_id());

  EX_CHECK_ERR(nemesis_err_flag, "Error writing Nemesis internal and border element maps to file!");
}



void Nemesis_IO_Helper::initialize(std::string title_in, const MeshBase & mesh, bool /*use_discontinuous*/)
{
  // Make sure that the reference passed in is really a DistributedMesh
  // const DistributedMesh & pmesh = cast_ref<const DistributedMesh &>(mesh);
  const MeshBase & pmesh = mesh;

  // If _write_as_dimension is nonzero, use it to set num_dim later in the Exodus file.
  if (_write_as_dimension)
    num_dim = _write_as_dimension;
  else if (_use_mesh_dimension_instead_of_spatial_dimension)
    num_dim = mesh.mesh_dimension();
  else
    num_dim = mesh.spatial_dimension();

  // According to Nemesis documentation, first call when writing should be to
  // ne_put_init_info().  Our reader doesn't actually call this, but we should
  // strive to be as close to a normal nemesis file as possible...
  this->put_init_info(this->n_processors(), 1, "p");


  // Gather global "initial" information for Nemesis.  This consists of
  // three parts labeled I, II, and III below...

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

  // Ready the node maps.  These have nothing to do with communication, they map
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
    for (const auto & id : this->internal_elem_ids)
      this->elem_mapi[cnt++] = libmesh_map_find(libmesh_elem_num_to_exodus, id);
  }

  {
    unsigned cnt = 0;
    for (const auto & id : this->border_elem_ids)
      this->elem_mapb[cnt++] = libmesh_map_find(libmesh_elem_num_to_exodus, id);
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
        std::set<std::pair<unsigned,unsigned>> & elem_set = it->second;

        // Resize the vectors to receive their payload
        this->elem_cmap_elem_ids[cnt].resize(elem_set.size());
        this->elem_cmap_side_ids[cnt].resize(elem_set.size());
        this->elem_cmap_proc_ids[cnt].resize(elem_set.size());

        std::set<std::pair<unsigned,unsigned>>::iterator elem_set_iter = elem_set.begin();

        // Pack the vectors with elem IDs, side IDs, and processor IDs.
        for (std::size_t j=0, eceis=this->elem_cmap_elem_ids[cnt].size(); j<eceis; ++j, ++elem_set_iter)
          {
            this->elem_cmap_elem_ids[cnt][j] =
              libmesh_map_find(libmesh_elem_num_to_exodus, elem_set_iter->first);
            this->elem_cmap_side_ids[cnt][j] = elem_set_iter->second;     // Side ID, this has already been converted above
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
    for (const auto & id : this->internal_node_ids)
      this->node_mapi[cnt++] = libmesh_map_find(libmesh_node_num_to_exodus, id);
  }

  {
    unsigned cnt=0;
    for (const auto & id : this->border_node_ids)
      this->node_mapb[cnt++] = libmesh_map_find(libmesh_node_num_to_exodus, id);
  }
}





void Nemesis_IO_Helper::compute_node_communication_maps()
{
  // Make sure there's no left-over information
  this->node_cmap_node_ids.clear();
  this->node_cmap_proc_ids.clear();

  libmesh_assert_less_equal
    (this->proc_nodes_touched_intersections.size(),
     std::size_t(this->num_node_cmaps));

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
        for (std::size_t j=0, nceis=this->node_cmap_node_ids[cnt].size(); j<nceis; ++j, ++node_set_iter)
          {
            this->node_cmap_node_ids[cnt][j] =
              libmesh_map_find(libmesh_node_num_to_exodus, *node_set_iter);
            this->node_cmap_proc_ids[cnt][j] = it->first;
          }

        // increment vector index to go to next processor
        cnt++;
      }
  } // end scope for packing

  // Print out the vectors we just packed
  if (verbose)
    {
      for (auto i : index_range(this->node_cmap_node_ids))
        {
          libMesh::out << "[" << this->processor_id() << "] nodes communicated to proc "
                       << this->node_cmap_ids[i]
                       << " = ";
          for (const auto & node_id : this->node_cmap_node_ids[i])
            libMesh::out << node_id << " ";
          libMesh::out << std::endl;
        }

      for (const auto & id_vec : this->node_cmap_node_ids)
        {
          libMesh::out << "[" << this->processor_id() << "] processor ID node communicated to = ";
          for (const auto & proc_id : id_vec)
            libMesh::out << proc_id << " ";
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
      for (const auto & node_cnt : node_cmap_node_cnts)
        libMesh::out << node_cnt << ", ";
      libMesh::out << std::endl;

      libMesh::out << "[" << this->processor_id() << "] node_cmap_ids = ";
      for (const auto & node_id : node_cmap_ids)
        libMesh::out << node_id << ", ";
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
      for (const auto & elem_cnt : elem_cmap_elem_cnts)
        libMesh::out << elem_cnt << ", ";
      libMesh::out << std::endl;

      libMesh::out << "[" << this->processor_id() << "] elem_cmap_ids = ";
      for (const auto & elem_id : elem_cmap_ids)
        libMesh::out << elem_id << ", ";
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

  for (const auto & elem : pmesh.active_local_element_ptr_range())
    {
      // Add this Elem's ID to all_elem_ids, later we will take the difference
      // between this set and the set of border_elem_ids, to get the set of
      // internal_elem_ids.
      all_elem_ids.insert(elem->id());

      // Will be set to true if element is determined to be a border element
      bool is_border_elem = false;

      // Construct a conversion object for this Element.  This will help us map
      // Libmesh numberings into Nemesis numberings for sides.
      const auto & conv = get_conversion(elem->type());

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
      for (auto node : elem->node_index_range())
        this->nodes_attached_to_local_elems.insert(elem->node_id(node));

      // Loop over element's neighbors, see if it has a neighbor which is off-processor
      for (auto n : elem->side_index_range())
        {
          if (elem->neighbor_ptr(n) != nullptr)
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
                  this->proc_border_elem_sets[ neighbor_proc_id ].emplace(elem->id(), nemesis_side_id);
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
      for (const auto & id : neighboring_processor_ids)
        libMesh::out << id << " ";
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
      for (const auto & [proc_id, set] : proc_border_elem_sets)
        {
          libMesh::out << "[" << this->processor_id() << "] "
                       << "Proc "
                       << proc_id << " communicates "
                       << set.size() << " elements." << std::endl;
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
      for (const auto & id : this->global_sideset_ids)
        libMesh::out << id << ", ";
      libMesh::out << std::endl;
    }

  // We also need global counts of sides in each of the sidesets.
  // Build a list of (elem, side, bc) tuples.
  typedef std::tuple<dof_id_type, unsigned short int, boundary_id_type> Tuple;
  std::vector<Tuple> bc_triples = pmesh.get_boundary_info().build_side_list();

  // Iterators to the beginning and end of the current range.
  std::vector<Tuple>::iterator
    it = bc_triples.begin(),
    new_end = bc_triples.end();

  while (it != new_end)
    {
      if (pmesh.elem_ref(std::get<0>(*it)).processor_id() != this->processor_id())
        {
          // Back up the new end iterators to prepare for swap
          --new_end;

          // Swap places, the non-local elem will now be "past-the-end"
          std::swap (*it, *new_end);
        }
      else // elem is local, go to next
        ++it;
    }

  // Erase from "new" end to old.
  bc_triples.erase(new_end, bc_triples.end());

  this->num_global_side_counts.clear(); // Make sure we don't have any leftover information
  this->num_global_side_counts.resize(this->global_sideset_ids.size());

  // Get the count for each global sideset ID
  for (auto i : index_range(global_sideset_ids))
    {
      int id = global_sideset_ids[i];
      this->num_global_side_counts[i] =
        cast_int<int>(std::count_if(bc_triples.begin(),
                                    bc_triples.end(),
                                    [id](const Tuple & t)->bool { return std::get<2>(t) == id; }));
    }

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] num_global_side_counts = ";
      for (const auto & cnt : this->num_global_side_counts)
        libMesh::out << cnt << ", ";
      libMesh::out << std::endl;
    }

  // Finally sum up the result
  this->comm().sum(this->num_global_side_counts);

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] num_global_side_counts = ";
      for (const auto & cnt : this->num_global_side_counts)
        libMesh::out << cnt << ", ";
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
      for (const auto & id : global_nodeset_ids)
        libMesh::out << id << ", ";
      libMesh::out << std::endl;

      libMesh::out << "[" << this->processor_id() << "] local_node_boundary_ids = ";
      for (const auto & id : local_node_boundary_ids)
        libMesh::out << id << ", ";
      libMesh::out << std::endl;
    }

  // 7.) We also need to know the number of nodes which is in each of the nodesets, globally.

  // Build list of (node-id, bc-id) tuples.
  typedef std::tuple<dof_id_type, boundary_id_type> Tuple;
  std::vector<Tuple> bc_tuples = pmesh.get_boundary_info().build_node_list();

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] boundary_node_list.size()="
                   << bc_tuples.size() << std::endl;
      libMesh::out << "[" << this->processor_id() << "] (boundary_node_id, boundary_id) = ";
      for (const auto & t : bc_tuples)
        libMesh::out << "(" << std::get<0>(t) << ", " << std::get<1>(t) << ") ";
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
  std::vector<Tuple>::iterator
    it = bc_tuples.begin(),
    new_end = bc_tuples.end();

  while (it != new_end)
    {
      if (pmesh.node_ptr(std::get<0>(*it))->processor_id() != this->processor_id())
        {
          // Back up the new end iterators to prepare for swap
          --new_end;

          // Swap places, the non-local node will now be "past-the-end"
          std::swap(*it, *new_end);
        }
      else // node is local, go to next
        ++it;
    }

  // Erase from "new" end to old end.
  bc_tuples.erase(new_end, bc_tuples.end());

  // Now we can do the local count for each ID...
  for (auto i : index_range(global_nodeset_ids))
    {
      int id = this->global_nodeset_ids[i];
      this->num_global_node_counts[i] =
        cast_int<int>(std::count_if(bc_tuples.begin(),
                                    bc_tuples.end(),
                                    [id](const Tuple & t)->bool { return std::get<1>(t) == id; }));
    }

  // And finally we can sum them up
  this->comm().sum(this->num_global_node_counts);

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] num_global_node_counts = ";
      for (const auto & cnt : num_global_node_counts)
        libMesh::out << cnt << ", ";
      libMesh::out << std::endl;
    }
}




void Nemesis_IO_Helper::compute_num_global_elem_blocks(const MeshBase & pmesh)
{
  // 1.) Loop over active local elements, build up set of subdomain IDs.
  std::set<subdomain_id_type> global_subdomain_ids;

  // This map keeps track of the number of elements in each subdomain over all processors
  std::map<subdomain_id_type, unsigned> global_subdomain_counts;

  for (const auto & elem : pmesh.active_local_element_ptr_range())
    {
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
      for (const auto & [subdomain_id, cnt] : global_subdomain_counts)
        {
          libMesh::out << "ID: "
                       << static_cast<unsigned>(subdomain_id)
                       << ", Count: " << cnt << ", ";
        }
      libMesh::out << std::endl;
    }

  // 6.) this->comm().sum up the number of elements in each block.  We know the global
  // subdomain IDs, so pack them into a vector one by one.  Use a vector of int since
  // that is what Nemesis wants
  this->global_elem_blk_cnts.resize(global_subdomain_ids.size());

  unsigned cnt=0;
  // Find the entry in the local map, note: if not found, will be created with 0 default value, which is OK...
  for (const auto & id : global_subdomain_ids)
    this->global_elem_blk_cnts[cnt++] = global_subdomain_counts[id];

  // Sum up subdomain counts from all processors
  this->comm().sum(this->global_elem_blk_cnts);

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] global_elem_blk_cnts = ";
      for (const auto & bc : this->global_elem_blk_cnts)
        libMesh::out << bc << ", ";
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
      for (const auto & id : this->global_elem_blk_ids)
        libMesh::out << id << ", ";
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
  this->build_subdomain_map(pmesh, true);

  // First loop over the elements to figure out which elements are in which subdomain
  for (const auto & elem : pmesh.active_local_element_ptr_range())
    {
      // Grab the nodes while we're here.
      for (auto n : elem->node_index_range())
        this->nodes_attached_to_local_elems.insert( elem->node_id(n) );
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
  for (const auto & id : nodes_attached_to_local_elems)
    {
      // I.e. given exodus_node_id,
      // exodus_node_num_to_libmesh[ exodus_node_id ] returns the
      // libmesh ID for that node, plus one.
      // Here we index libMesh IDs with an offset of 1 because they're
      // the literal numbers that get written to the exodus file, but
      // we index Exodus IDs with an offset of 0 because we read that
      // exodus data into a C++ vector.  Confused yet?
      this->exodus_node_num_to_libmesh.push_back(id+1);

      // Likewise, given libmesh_node_id,
      // libmesh_node_num_to_exodus[ libmesh_node_id] returns the
      // *Exodus* ID for that node.  Unlike the
      // exodus_node_num_to_libmesh vector above, this one is a
      // std::map.  We're never handing a data buffer from it over to
      // another API so we don't need to do any weird offsets with it.
      this->libmesh_node_num_to_exodus[id] =
        this->exodus_node_num_to_libmesh.size(); // should never be zero...
    }

  // Now we're going to loop over the subdomain map and build a few things right
  // now that we'll use later.

  // First make sure our data structures don't have any leftover data...
  this->exodus_elem_num_to_libmesh.clear();
  this->block_ids.clear();
  this->libmesh_elem_num_to_exodus.clear();

  // Now loop over each subdomain and get a unique numbering for the elements
  for (auto & [block_id, elem_ids_this_subdomain] : this->_subdomain_map)
    {
      block_ids.push_back(block_id);

      // The code below assumes this subdomain block is not empty, make sure that's the case!
      libmesh_error_msg_if(elem_ids_this_subdomain.size() == 0,
                           "Error, no element IDs found in subdomain " << block_id);

      // Use the first element in this block to get representative information.
      // Note that Exodus assumes all elements in a block are of the same type!
      // We are using that same assumption here!
      const auto & conv = get_conversion
        (pmesh.elem_ref(elem_ids_this_subdomain[0]).type());
      this->num_nodes_per_elem =
        pmesh.elem_ref(elem_ids_this_subdomain[0]).n_nodes();

      // Get a reference to the connectivity vector for this subdomain.  This vector
      // is most likely empty, we are going to fill it up now.
      std::vector<int> & current_block_connectivity = this->block_id_to_elem_connectivity[block_id];

      // Just in case it's not already empty...
      current_block_connectivity.clear();
      current_block_connectivity.resize(elem_ids_this_subdomain.size() * this->num_nodes_per_elem);

      for (auto i : index_range(elem_ids_this_subdomain))
        {
          auto elem_id = elem_ids_this_subdomain[i];

          // Set the number map for elements
          // exodus_elem_num_to_libmesh[ exodus_node_id ] returns the
          // libmesh ID for that element, plus one.
          // Like with nodes above, we index libMesh IDs with an
          // offset of 1 because they're the literal numbers that get
          // written to the exodus file, but we index Exodus IDs with
          // an offset of 0 because we read that exodus data into a
          // C++ vector.
          this->exodus_elem_num_to_libmesh.push_back(elem_id+1);

          // Likewise, given libmesh elem_id,
          // libmesh_elem_num_to_exodus[ elem_id ] returns the
          // *Exodus* ID for that node.  Unlike the
          // exodus_elem_num_to_libmesh vector above, this one is a
          // std::map.  We're never handing a data buffer from it over to
          // another API so we don't need to do any weird offsets with it.
          this->libmesh_elem_num_to_exodus[elem_id] =
            this->exodus_elem_num_to_libmesh.size();

          const Elem & elem = pmesh.elem_ref(elem_id);

          // Exodus/Nemesis want every block to have the same element type
          // libmesh_assert_equal_to (elem->type(), conv.libmesh_elem_type());

          // But we can get away with writing e.g. HEX8 and INFHEX8 in
          // the same block...
          libmesh_assert_equal_to (elem.n_nodes(), Elem::build(conv.libmesh_elem_type(), nullptr)->n_nodes());

          for (auto j : make_range(this->num_nodes_per_elem))
            {
              const unsigned int connect_index   = (i*this->num_nodes_per_elem)+j;
              const unsigned int elem_node_index = conv.get_inverse_node_map(j); // inverse node map is used for writing

              current_block_connectivity[connect_index] =
                libmesh_map_find(libmesh_node_num_to_exodus,
                                 elem.node_id(elem_node_index));
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
  std::map<unsigned, std::set<unsigned>> proc_nodes_touched;


  // We are going to create a lot of intermediate data structures here, so make sure
  // as many as possible all cleaned up by creating scope!
  {
    // Loop over active (not just active local) elements, make sets of node IDs for each
    // processor which has an element that "touches" a node.
    for (const auto & elem : pmesh.active_element_ptr_range())
      {
        // Get reference to the set for this processor.  If it does not exist
        // it will be created.
        std::set<unsigned> & set_p = proc_nodes_touched[ elem->processor_id() ];

        // Insert all nodes touched by this element into the set
        for (auto node : elem->node_index_range())
          set_p.insert(elem->node_id(node));
      }

    if (verbose)
      {
        libMesh::out << "[" << this->processor_id()
                     << "] proc_nodes_touched contains "
                     << proc_nodes_touched.size()
                     << " sets of nodes."
                     << std::endl;

        for (const auto & [proc_id, set] : proc_nodes_touched)
          libMesh::out << "[" << this->processor_id()
                       << "] proc_nodes_touched[" << proc_id << "] has "
                       << set.size()
                       << " entries."
                       << std::endl;
      }


    // Loop over all the sets we just created and compute intersections with the
    // this processor's set.  Obviously, don't intersect with ourself.
    this->proc_nodes_touched_intersections.clear();
    for (auto & [proc_id, other_set] : proc_nodes_touched)
      {
        // Don't compute intersections with ourself
        if (proc_id == this->processor_id())
          continue;

        std::set<unsigned int> this_intersection;

        // Otherwise, compute intersection with other processor and ourself
        std::set<unsigned> & my_set = proc_nodes_touched[this->processor_id()];

        std::set_intersection(my_set.begin(), my_set.end(),
                              other_set.begin(), other_set.end(),
                              std::inserter(this_intersection, this_intersection.end()));

        if (!this_intersection.empty())
          this->proc_nodes_touched_intersections.emplace
            (proc_id, std::move(this_intersection));
      }

    if (verbose)
      {
        for (const auto & [proc_id, set] : proc_nodes_touched_intersections)
          libMesh::out << "[" << this->processor_id()
                       << "] this->proc_nodes_touched_intersections[" << proc_id << "] has "
                       << set.size()
                       << " entries."
                       << std::endl;
      }

    // The number of node communication maps is the number of other processors
    // with which we share nodes.
    this->num_node_cmaps =
      cast_int<int>(proc_nodes_touched_intersections.size());

    // We can't be connecting to more processors than exist outside
    // ourselves
    libmesh_assert_less (this->num_node_cmaps, this->n_processors());

    // Compute the set_union of all the preceding intersections.  This will be the set of
    // border node IDs for this processor.
    for (auto & pr : proc_nodes_touched_intersections)
      {
        std::set<unsigned> & other_set = pr.second;
        std::set<unsigned> intermediate_result; // Don't think we can insert into one of the sets we're unioning...

        std::set_union(this->border_node_ids.begin(), this->border_node_ids.end(),
                       other_set.begin(), other_set.end(),
                       std::inserter(intermediate_result, intermediate_result.end()));

        // Swap our intermediate result into the final set
        this->border_node_ids.swap(intermediate_result);
      }

    libmesh_assert_less_equal
      (this->proc_nodes_touched_intersections.size(),
       std::size_t(this->num_node_cmaps));

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
  std::map<boundary_id_type, std::vector<int>> local_node_boundary_id_lists;

  // FIXME: We should build this list only one time!!  We already built it above, but we
  // did not have the libmesh to exodus node mapping at that time... for now we'll just
  // build it here again, hopefully it's small relative to the size of the entire mesh.

  // Build list of (node-id, bc-id) tuples.
  typedef std::tuple<dof_id_type, boundary_id_type> Tuple;
  std::vector<Tuple> bc_tuples = mesh.get_boundary_info().build_node_list();

  if (verbose)
    {
      libMesh::out << "[" << this->processor_id() << "] boundary_node_list.size()="
                   << bc_tuples.size() << std::endl;
      libMesh::out << "[" << this->processor_id() << "] (boundary_node_id, boundary_id) = ";
      for (const auto & t : bc_tuples)
        libMesh::out << "(" << std::get<0>(t) << ", " << std::get<1>(t) << ") ";
      libMesh::out << std::endl;
    }

  // For each node in the node list, add it to the vector of node IDs for that
  // set for the local processor.  This will be used later when writing Exodus
  // nodesets.
  for (const auto & t : bc_tuples)
    {
      // Don't try to grab a reference to the vector unless the current node is attached
      // to a local element.  Otherwise, another processor will be responsible for writing it in its nodeset.
      if (const auto it = this->libmesh_node_num_to_exodus.find(std::get<0>(t));
          it != this->libmesh_node_num_to_exodus.end())
        {
          // Get reference to the vector where this node ID will be inserted.  If it
          // doesn't yet exist, this will create it.
          std::vector<int> & current_id_set = local_node_boundary_id_lists[std::get<1>(t)];

          // Push back Exodus-mapped node ID for this set
          // TODO: reserve space in these vectors somehow.
          current_id_set.push_back( it->second );
        }
    }

  // See what we got
  if (verbose)
    {
      for (const auto & [bndry_id, set] : local_node_boundary_id_lists)
        {
          libMesh::out << "[" << this->processor_id() << "] ID: " << bndry_id << ", ";

          // Libmesh node ID (Exodus Node ID)
          for (const auto & id : set)
            libMesh::out << id << ", ";
          libMesh::out << std::endl;
        }
    }

  // Loop over *global* nodeset IDs, call the Exodus API.  Note that some nodesets may be empty
  // for a given processor.
  if (global_nodeset_ids.size() > 0)
    {
      NamesData names_table(global_nodeset_ids.size(), MAX_STR_LENGTH);

      for (const auto & nodeset_id : this->global_nodeset_ids)
        {
          const std::string & current_ns_name =
            mesh.get_boundary_info().get_nodeset_name
              (cast_int<boundary_id_type>(nodeset_id));

          // Store this name in a data structure that will be used to
          // write sideset names to file.
          names_table.push_back_entry(current_ns_name);

          if (verbose)
            {
              libMesh::out << "[" << this->processor_id()
                           << "] Writing out Exodus nodeset info for ID: " << nodeset_id
                           << ", Name: " << current_ns_name
                           << std::endl;
            }

          // Convert current global_nodeset_id into an exodus ID, which can't be zero...
          int exodus_id = nodeset_id;

          /*
          // Exodus can't handle zero nodeset IDs (?)  Use max short here since
          // when libmesh reads it back in, it will want to store it as a short...
          if (exodus_id==0)
          exodus_id = std::numeric_limits<short>::max();
          */

          // Try to find this boundary ID in the local list we created
          if (const auto it = local_node_boundary_id_lists.find (cast_int<boundary_id_type>(nodeset_id));
              it == local_node_boundary_id_lists.end())
            {
              // No nodes found for this boundary ID on this processor
              if (verbose)
                libMesh::out << "[" << this->processor_id()
                             << "] No nodeset data for ID: " << nodeset_id
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
              const std::vector<int> & current_nodeset_ids = it->second;

              // Call the Exodus interface to write the parameters of this node set
              this->ex_err = exII::ex_put_node_set_param(this->ex_id,
                                                         exodus_id,
                                                         current_nodeset_ids.size(),
                                                         0  /* No distribution factors */);

              EX_CHECK_ERR(this->ex_err, "Error writing nodeset parameters in Nemesis");

              // Call Exodus interface to write the actual node IDs for this boundary ID
              this->ex_err = exII::ex_put_node_set(this->ex_id,
                                                   exodus_id,
                                                   current_nodeset_ids.data());

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
  std::map<boundary_id_type, std::vector<int>> local_elem_boundary_id_lists;
  std::map<boundary_id_type, std::vector<int>> local_elem_boundary_id_side_lists;

  // FIXME: We already built this list once, we should reuse that information!
  std::vector<std::tuple<dof_id_type, unsigned short int, boundary_id_type>> bndry_elem_side_id_list =
    mesh.get_boundary_info().build_side_list();

  // Integer looping, skipping non-local elements
  for (const auto & t : bndry_elem_side_id_list)
    {
      // Get pointer to current Elem
      const Elem * elem = mesh.elem_ptr(std::get<0>(t));

      std::vector<const Elem *> family;
#ifdef LIBMESH_ENABLE_AMR
      // We need to build up active elements if AMR is enabled and add
      // them to the exodus sidesets instead of the potentially inactive "parent" elements
      // Technically we don't need to "reset" the tree since the vector was just created.
      elem->active_family_tree_by_side(family, std::get<1>(t), /*reset tree=*/false);
#else
      // If AMR is not even enabled, just push back the element itself
      family.push_back( elem );
#endif

      // Loop over all the elements in the family tree, store their converted IDs
      // and side IDs to the map's vectors.  TODO: Somehow reserve enough space for these
      // push_back's...
      for (const auto & tree_elem : family)
        {
          const dof_id_type f_id = tree_elem->id();
          const Elem & f = mesh.elem_ref(f_id);

          // If element is local, process it
          if (f.processor_id() == this->processor_id())
            {
              const auto & conv = get_conversion(f.type());

              // Use the libmesh to exodus data structure map to get the proper sideset IDs
              // The data structure contains the "collapsed" contiguous ids.
              //
              // We know the parent element is local, but let's be absolutely sure that all the children have been
              // actually mapped to Exodus IDs before we blindly try to add them...
              local_elem_boundary_id_lists[ std::get<2>(t) ].push_back( libmesh_map_find(libmesh_elem_num_to_exodus, f_id) );
              local_elem_boundary_id_side_lists[ std::get<2>(t) ].push_back(conv.get_inverse_side_map( std::get<1>(t) ));
            }
        }
    }


  // Loop over *global* sideset IDs, call the Exodus API.  Note that some sidesets may be empty
  // for a given processor.
  if (global_sideset_ids.size() > 0)
    {
      NamesData names_table(global_sideset_ids.size(), MAX_STR_LENGTH);

      for (const auto & exodus_id : this->global_sideset_ids)
        {
          const std::string & current_ss_name =
            mesh.get_boundary_info().get_sideset_name
              (cast_int<boundary_id_type>(exodus_id));

          // Store this name in a data structure that will be used to
          // write sideset names to file.
          names_table.push_back_entry(current_ss_name);

          if (verbose)
            {
              libMesh::out << "[" << this->processor_id()
                           << "] Writing out Exodus sideset info for ID: " << exodus_id
                           << ", Name: " << current_ss_name
                           << std::endl;
            }

          // Try to find this boundary ID in the local list we created
          if (const auto it = local_elem_boundary_id_lists.find (cast_int<boundary_id_type>(exodus_id));
              it == local_elem_boundary_id_lists.end())
            {
              // No sides found for this boundary ID on this processor
              if (verbose)
                libMesh::out << "[" << this->processor_id()
                             << "] No sideset data for ID: " << exodus_id
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
              // Get reference to the vector of elem IDs
              const std::vector<int> & current_sideset_elem_ids = it->second;

              // Get reference to the vector of side IDs
              std::vector<int> & current_sideset_side_ids =
                libmesh_map_find(local_elem_boundary_id_side_lists,
                                 cast_int<boundary_id_type>(exodus_id));

              // Call the Exodus interface to write the parameters of this side set
              this->ex_err = exII::ex_put_side_set_param(this->ex_id,
                                                         exodus_id,
                                                         current_sideset_elem_ids.size(),
                                                         0  /* No distribution factors */);

              EX_CHECK_ERR(this->ex_err, "Error writing sideset parameters in Nemesis");

              // Call Exodus interface to write the actual side IDs for this boundary ID
              this->ex_err = exII::ex_put_side_set(this->ex_id,
                                                   exodus_id,
                                                   current_sideset_elem_ids.data(),
                                                   current_sideset_side_ids.data());

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
  auto local_num_nodes = this->exodus_node_num_to_libmesh.size();

  x.resize(local_num_nodes);
  y.resize(local_num_nodes);
  z.resize(local_num_nodes);

  // Just loop over our list outputting the nodes the way we built the map
  for (auto i : make_range(local_num_nodes))
    {
      const Point & pt = mesh.point(this->exodus_node_num_to_libmesh[i]-1);
      x[i]=pt(0);
      y[i]=pt(1);
      z[i]=pt(2);
    }

  if (local_num_nodes)
    {
      // Call Exodus API to write nodal coordinates...
      ex_err = exII::ex_put_coord
        (ex_id,
         x.empty() ? nullptr : MappedOutputVector(x, _single_precision).data(),
         y.empty() ? nullptr : MappedOutputVector(y, _single_precision).data(),
         z.empty() ? nullptr : MappedOutputVector(z, _single_precision).data());
      EX_CHECK_ERR(ex_err, "Error writing node coordinates");

      // And write the nodal map we created for them
      ex_err = exII::ex_put_node_num_map(ex_id, this->exodus_node_num_to_libmesh.data());
      EX_CHECK_ERR(ex_err, "Error writing node num map");
    }
  else // Does the Exodus API want us to write empty nodal coordinates?
    {
      ex_err = exII::ex_put_coord(ex_id, nullptr, nullptr, nullptr);
      EX_CHECK_ERR(ex_err, "Error writing empty node coordinates");

      ex_err = exII::ex_put_node_num_map(ex_id, nullptr);
      EX_CHECK_ERR(ex_err, "Error writing empty node num map");
    }
}





void Nemesis_IO_Helper::write_elements(const MeshBase & mesh, bool /*use_discontinuous*/)
{
  // Only write elements if there are elements blocks available.
  if (this->num_elem_blks_global > 0)
    {
      // Data structure to store element block names that will be used to
      // write the element block names to file.
      NamesData names_table(this->num_elem_blks_global, MAX_STR_LENGTH);

      // Loop over all blocks, even if we don't have elements in each block.
      // If we don't have elements we need to write out a 0 for that block...
      for (auto i : make_range(this->num_elem_blks_global))
        {
          // Even if there are no elements for this block on the current
          // processor, we still want to write its name to file, if
          // possible. MeshBase::subdomain_name() will just return an
          // empty string if there is no name associated with the current
          // block.
          names_table.push_back_entry
            (mesh.subdomain_name(cast_int<subdomain_id_type>(this->global_elem_blk_ids[i])));

          // Search for the current global block ID in the map
          if (const auto it = this->block_id_to_elem_connectivity.find( this->global_elem_blk_ids[i] );
              it == this->block_id_to_elem_connectivity.end())
            {
              // If not found, write a zero to file....
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
              const std::vector<int> & this_block_connectivity = it->second;
              std::vector<dof_id_type> & elements_in_this_block = this->_subdomain_map[block];

              // Use the first element in this block to get representative information.
              // Note that Exodus assumes all elements in a block are of the same type!
              // We are using that same assumption here!
              const auto & conv =
                get_conversion(mesh.elem_ref(elements_in_this_block[0]).type());

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
                                              this_block_connectivity.data());
              EX_CHECK_ERR(ex_err, "Error writing element connectivities from Nemesis.");
            }
        } // end loop over global block IDs

      // Only call this once, not in the loop above!
      ex_err = exII::ex_put_elem_num_map(ex_id,
                                         exodus_elem_num_to_libmesh.empty() ? nullptr : exodus_elem_num_to_libmesh.data());
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
          Number value = values[(this->exodus_node_num_to_libmesh[i]-1)*num_vars + c];
          real_parts[i] = value.real();
          imag_parts[i] = value.imag();
          magnitudes[i] = std::abs(value);
        }
      write_nodal_values(3*c+1,real_parts,timestep);
      write_nodal_values(3*c+2,imag_parts,timestep);
      write_nodal_values(3*c+3,magnitudes,timestep);
#else
      std::vector<Number> cur_soln(this->num_nodes);

      // Copy out this variable's solution
      for (int i=0; i<this->num_nodes; i++)
        cur_soln[i] = values[(this->exodus_node_num_to_libmesh[i]-1)*num_vars + c];

      write_nodal_values(c+1,cur_soln,timestep);
#endif
    }
}



void Nemesis_IO_Helper::write_nodal_solution(const NumericVector<Number> & parallel_soln,
                                             const std::vector<std::string> & names,
                                             int timestep,
                                             const std::vector<std::string> & output_names)
{
  int num_vars = cast_int<int>(names.size());

  for (int c=0; c<num_vars; c++)
    {
      // Find the position of names[c] in the output_names vector, if it exists.
      auto pos = std::find(output_names.begin(), output_names.end(), names[c]);

      // Skip names[c] if it's not supposed to be output.
      if (pos == output_names.end())
        continue;

      // Compute the (zero-based) index which determines which
      // variable this will be as far as Nemesis is concerned.  This
      // will be used below in the write_nodal_values() call.
      int variable_name_position =
        cast_int<int>(std::distance(output_names.begin(), pos));

      // Fill up a std::vector with the dofs for the current variable
      std::vector<numeric_index_type> required_indices(this->num_nodes);

      for (int i=0; i<this->num_nodes; i++)
        required_indices[i] = static_cast<dof_id_type>(this->exodus_node_num_to_libmesh[i]-1) * num_vars + c;

      // Get the dof values required to write just our local part of
      // the solution vector.
      std::vector<Number> local_soln;
      parallel_soln.localize(local_soln, required_indices);

#ifndef LIBMESH_USE_COMPLEX_NUMBERS
      // Call the ExodusII_IO_Helper function to write the data.
      write_nodal_values(variable_name_position + 1, local_soln, timestep);
#else
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
      write_nodal_values(3 * variable_name_position + 1, real_parts, timestep);
      write_nodal_values(3 * variable_name_position + 2, imag_parts, timestep);
      write_nodal_values(3 * variable_name_position + 3, magnitudes, timestep);
#endif
    }
}



void Nemesis_IO_Helper::write_nodal_solution(const EquationSystems & es,
                                             const std::vector<std::pair<unsigned int, unsigned int>> & var_nums,
                                             int timestep,
                                             const std::vector<std::string> & output_names)
{
  const MeshBase & mesh = es.get_mesh();

  // FIXME - half this code might be replaceable with a call to
  // EquationSystems::build_parallel_solution_vector()...

  for (auto [sys_num, var] : var_nums)
    {
      const System & sys = es.get_system(sys_num);
      const std::string & name = sys.variable_name(var);

      auto pos = std::find(output_names.begin(), output_names.end(), name);

      // Skip this name if it's not supposed to be output.
      if (pos == output_names.end())
        continue;

      // Compute the (zero-based) index which determines which
      // variable this will be as far as Nemesis is concerned.  This
      // will be used below in the write_nodal_values() call.
      int variable_name_position =
        cast_int<int>(std::distance(output_names.begin(), pos));

      // Fill up a std::vector with the dofs for the current variable
      std::vector<numeric_index_type> required_indices(this->num_nodes);

      // Get the dof values required to write just our local part of
      // the solution vector.
      std::vector<Number> local_soln;

      const FEType type = sys.variable_type(var);
      if (type.family == SCALAR)
        {
          std::vector<numeric_index_type> scalar_indices;
          sys.get_dof_map().SCALAR_dof_indices(scalar_indices, var);
          for (int i=0; i<this->num_nodes; i++)
            required_indices[i] = scalar_indices[0];
          sys.current_local_solution->get(required_indices, local_soln);
        }
      else
        {
          // If we have DoFs at all nodes, e.g. for isoparametric
          // elements, this is easy:
          bool found_all_indices = true;
          for (int i=0; i<this->num_nodes; i++)
            {
              const Node & node = mesh.node_ref(this->exodus_node_num_to_libmesh[i]-1);
              if (node.n_comp(sys_num, var))
                required_indices[i] = node.dof_number(sys_num, var, 0);
              else
                {
                  found_all_indices = false;
                  break;
                }
            }

          if (found_all_indices)
            sys.current_local_solution->get(required_indices, local_soln);
          // Fine, we'll do it the hard way
          if (!found_all_indices)
            {
              local_soln.resize(num_nodes);

              const Variable & var_description = sys.variable(var);
              const DofMap & dof_map           = sys.get_dof_map();

              NumericVector<Number> & sys_soln(*sys.current_local_solution);
              std::vector<Number>      elem_soln;   // The finite element solution
              std::vector<Number>      nodal_soln;  // The FE solution interpolated to the nodes
              std::vector<dof_id_type> dof_indices; // The DOF indices for the finite element

              for (const auto & elem : mesh.active_local_element_ptr_range())
                if (var_description.active_on_subdomain(elem->subdomain_id()))
                  {
                    dof_map.dof_indices (elem, dof_indices, var);
                    elem_soln.resize(dof_indices.size());

                    for (auto i : index_range(dof_indices))
                      elem_soln[i] = sys_soln(dof_indices[i]);

                    FEInterface::nodal_soln (elem->dim(),
                                             type,
                                             elem,
                                             elem_soln,
                                             nodal_soln);

                    // infinite elements should be skipped...
                    if (!elem->infinite())
                      for (auto n : elem->node_index_range())
                        {
                          const std::size_t exodus_num =
                            libmesh_node_num_to_exodus[elem->node_id(n)];
                          libmesh_assert_greater(exodus_num, 0);
                          libmesh_assert_less(exodus_num-1, local_soln.size());
                          local_soln[exodus_num-1] = nodal_soln[n];
                        }
                  }
            }
        }

#ifndef LIBMESH_USE_COMPLEX_NUMBERS
      // Call the ExodusII_IO_Helper function to write the data.
      write_nodal_values(variable_name_position + 1, local_soln, timestep);
#else
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
      write_nodal_values(3 * variable_name_position + 1, real_parts, timestep);
      write_nodal_values(3 * variable_name_position + 2, imag_parts, timestep);
      write_nodal_values(3 * variable_name_position + 3, magnitudes, timestep);
#endif
    }
}



void
Nemesis_IO_Helper::initialize_element_variables(std::vector<std::string> names,
                                                const std::vector<std::set<subdomain_id_type>> & vars_active_subdomains)
{
  // Quick return if there are no element variables to write
  if (names.size() == 0)
    return;

  // Quick return if we have already called this function
  if (_elem_vars_initialized)
    return;

  // Be sure that variables in the file match what we are asking for
  if (num_elem_vars > 0)
    {
      this->check_existing_vars(ELEMENTAL, names, this->elem_var_names);
      return;
    }

  // Set the flag so we can skip this stuff on subsequent calls to
  // initialize_element_variables()
  _elem_vars_initialized = true;

  this->write_var_names(ELEMENTAL, names);

  // Create a truth table from global_elem_blk_ids and the information
  // in vars_active_subdomains. Create a truth table of
  // size global_elem_blk_ids.size() * names.size().
  std::vector<int> truth_tab(global_elem_blk_ids.size() * names.size());
  for (auto blk : index_range(global_elem_blk_ids))
    for (auto var : index_range(names))
      if (vars_active_subdomains[var].empty() ||
          vars_active_subdomains[var].count(cast_int<subdomain_id_type>(global_elem_blk_ids[blk])))
        truth_tab[names.size() * blk + var] = 1;

  // Write truth table to file.
  if (truth_tab.size())
    {
      ex_err = exII::ex_put_elem_var_tab(ex_id,
                                         cast_int<int>(global_elem_blk_ids.size()),
                                         cast_int<int>(names.size()),
                                         truth_tab.data());
      EX_CHECK_ERR(ex_err, "Error writing element truth table.");
    }
}



void
Nemesis_IO_Helper::write_element_values(const MeshBase & mesh,
                                        const EquationSystems & es,
                                        const std::vector<std::pair<unsigned int, unsigned int>> &var_nums,
                                        int timestep,
                                        const std::vector<std::set<subdomain_id_type>> & vars_active_subdomains)
{
  // For each variable in names,
  //   For each subdomain in subdomain_map,
  //     If this (subdomain, variable) combination is active
  //       For each component in variable
  //         Extract element values into local_soln (localize is a collective)
  //         Write local_soln to file
  //   Update var_ctr with number of vector components for variable
  //
  unsigned int var_ctr = 0;
  for (auto v : index_range(var_nums))
    {
      const unsigned int sys_num = var_nums[v].first;
      const unsigned int var = var_nums[v].second;
      const System & system = es.get_system(sys_num);

      // We need to check if the constant monomial is a scalar or a vector and set the number of
      // components as the mesh spatial dimension for the latter as per es.find_variable_numbers().
      // Even for the case where a variable is not active on any subdomain belonging to the
      // processor, we still need to know this number to update 'var_ctr'.
      const unsigned int n_comps =
        (system.variable_type(var) == FEType(CONSTANT, MONOMIAL_VEC)) ? mesh.spatial_dimension() : 1;

      // Get list of active subdomains for variable v
      const auto & active_subdomains = vars_active_subdomains[v];

      for (const int sbd_id_int : global_elem_blk_ids)
        {
          const subdomain_id_type sbd_id =
            cast_int<subdomain_id_type>(sbd_id_int);
          auto it = this->_subdomain_map.find(sbd_id);
          const std::vector<dof_id_type> empty_vec;
          const std::vector<dof_id_type> & elem_ids =
            (it == this->_subdomain_map.end()) ? empty_vec : it->second;

          // Possibly skip this (variable, subdomain) combination. Also, check that there is at
          // least one element on the subdomain... Indeed, it is possible to have zero elements,
          // e.g., when running "adaptivity_ex3" in parallel with the 'dimension=1' argument.
          if ((active_subdomains.empty() || active_subdomains.count(sbd_id)) && elem_ids.size())
            {
              std::vector<numeric_index_type> required_indices;
              required_indices.reserve(elem_ids.size());

              // The number of DOF components needs to be equal to the expected number so that we
              // know where to store data to correctly correspond to variable names - verify this by
              // accessing the n_comp method for the last element ID, which should return the same
              // value for all elements on a given subdomain, so we only need to check this once.
              libmesh_assert_equal_to(n_comps, mesh.elem_ref(elem_ids.back()).n_comp(sys_num, var));

              // Loop through the DOFs of the variable and write the values for it on each element.
              // The variable name should have been decomposed by es.find_variable_numbers().
              for (unsigned int comp = 0; comp < n_comps; ++comp)
              {
                for (const auto & id : elem_ids)
                    required_indices.push_back(mesh.elem_ref(id).dof_number(sys_num, var, comp));

                std::vector<Number> local_soln;
                system.current_local_solution->get(required_indices, local_soln);

                // reset for the next component
                required_indices.clear();

                // It's possible that there's nothing for us to write:
                // we may not be responsible for any elements on the
                // current subdomain.  We did still have to participate
                // in the localize() call above, however, since it is a
                // collective.
                if (local_soln.size())
                  {
#ifdef LIBMESH_USE_COMPLEX_NUMBERS
                    int stride = write_complex_abs ? 3 : 2;
                    std::vector<Real> local_soln_buffer(local_soln.size());
                    std::transform(local_soln.begin(), local_soln.end(),
                                   local_soln_buffer.begin(), [](Number n) { return n.real(); });
                    ex_err = exII::ex_put_elem_var(ex_id,
                                                   timestep,
                                                   static_cast<int>(stride*(var_ctr+comp)+1),
                                                   static_cast<int>(sbd_id),
                                                   static_cast<int>(local_soln.size()),
                                                   MappedOutputVector(local_soln_buffer, _single_precision).data());
                    EX_CHECK_ERR(ex_err, "Error writing element real values.");

                    std::transform(local_soln.begin(), local_soln.end(),
                                   local_soln_buffer.begin(), [](Number n) { return n.imag(); });
                    ex_err = exII::ex_put_elem_var(ex_id,
                                                   timestep,
                                                   static_cast<int>(stride*(var_ctr+comp)+2),
                                                   static_cast<int>(sbd_id),
                                                   static_cast<int>(local_soln.size()),
                                                   MappedOutputVector(local_soln_buffer, _single_precision).data());
                    EX_CHECK_ERR(ex_err, "Error writing element imaginary values.");

                    if (write_complex_abs)
                      {
                        std::transform(local_soln.begin(), local_soln.end(),
                                       local_soln_buffer.begin(), [](Number n) { return std::abs(n); });
                        ex_err = exII::ex_put_elem_var(ex_id,
                                                       timestep,
                                                       static_cast<int>(stride*(var_ctr+comp)+2),
                                                       static_cast<int>(sbd_id),
                                                       static_cast<int>(local_soln.size()),
                                                       MappedOutputVector(local_soln_buffer, _single_precision).data());
                        EX_CHECK_ERR(ex_err, "Error writing element magnitudes.");
                      }
#else // LIBMESH_USE_COMPLEX_NUMBERS
                    ex_err = exII::ex_put_elem_var(ex_id,
                                                   timestep,
                                                   static_cast<int>(var_ctr+comp+1),
                                                   static_cast<int>(sbd_id),
                                                   static_cast<int>(local_soln.size()),
                                                   MappedOutputVector(local_soln, _single_precision).data());
                    EX_CHECK_ERR(ex_err, "Error writing element values.");
#endif // LIBMESH_USE_COMPLEX_NUMBERS
                  }
              } // end loop over vector components
            }
        } // end loop over active subdomains

        var_ctr += n_comps;
    } // end loop over vars

  this->update();
}



void Nemesis_IO_Helper::read_var_names_impl(const char * var_type,
                                            int & count,
                                            std::vector<std::string> & result)
{
  // Most of what we need to do is the same as for Exodus
  this->ExodusII_IO_Helper::read_var_names_impl(var_type, count, result);

  // But with tests where we have more processors than elements,
  // Nemesis doesn't let us put variable names in files written by
  // processors owning nothing, but we may still *need* those
  // variable names on every processor, so let's sync them up...

  processor_id_type pid_broadcasting_names = this->processor_id();
  const std::size_t n_names = result.size();
  if (!n_names)
    pid_broadcasting_names = DofObject::invalid_processor_id;

  libmesh_assert(this->comm().semiverify
                 (n_names ? nullptr : &n_names));

  this->comm().min(pid_broadcasting_names);

  if (pid_broadcasting_names != DofObject::invalid_processor_id)
    this->comm().broadcast(result, pid_broadcasting_names);
}



std::string Nemesis_IO_Helper::construct_nemesis_filename(std::string_view base_filename)
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
