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

#ifndef __nemesis_io_helper_h__
#define __nemesis_io_helper_h__

#include "libmesh_config.h"

#if defined(HAVE_NEMESIS_API) && defined(HAVE_EXODUS_API)

#include <vector>
#include "exodusII_io_helper.h"

// The Nemesis API header file.  Should already be
// correctly extern C'd but it doesn't hurt :)
namespace Nemesis {
  extern "C" {
#include "ne_nemesisI.h"
  }
}


/**
 * This is the \p Nemesis_IO_Helper class.  Think of it as
 * a big struct with storage for all the stuff one might
 * want to pull from a Nemesis file.  Also contains an
 * ExodusII_IO_Helper object, since Nemesis is based on
 * the same file format.
 *
 * @author Johw W. Peterson, 2008.
 */
class Nemesis_IO_Helper
{
public:
  /**
   * Constructor.
   */
  Nemesis_IO_Helper(bool verbose);

  /**
   * Destructor.
   */
  ~Nemesis_IO_Helper();
  
  // Member functions.  These just allocate memory for you and call the Nemesis
  // routines of the same name.  They also handle error checking for the Nemesis
  // return value.  Be careful calling these at random, some depend on others
  // being called first...
  void get_init_global();
  void get_ss_param_global();
  void get_ns_param_global();
  void get_eb_info_global();
  void get_init_info();
  void get_loadbal_param();
  void get_elem_map();
  void get_node_map();
  void get_cmap_params();
  void get_node_cmap();
  void get_elem_cmap();
  
  // Member data

  /**
   * Instance of the Exodus IO Helper.  We call the Exodus API through
   * this object.  Instead of creating forwarding functions for
   * everything in the ExodusII_IO_Helper class, just call them
   * directly through this object!
   */
  ExodusII_IO_Helper ex2helper;

  /**
   * All (?) Nemesis functions return an int.  If it's negative that signals an error!
   * Internally, we use the ExodusII_IO_Helper::check_err() function to check for errors.
   */
  int nemesis_err_flag;
    
  /**
   * Global initial information.  The names are self-explanatory
   * for the most part.  Used with Nemesis::ne_get_init_global().
   */
  int num_nodes_global;
  int num_elems_global;
  int num_elem_blks_global;
  int num_node_sets_global;
  int num_side_sets_global;

  /**
   * The number of processors for which the NEMESIS I file was created.
   * To be used with Nemesis::ne_get_init_info().
   */
  int num_proc;

  /**
   * The number of processors for which the NEMESIS I file stores information.
   * This is generally equal to 1 (1 CPU/file) at least for the splitting Derek gave us.
   * To be used with Nemesis::ne_get_init_info().   
   */
  int num_proc_in_file;

  /**
   * The type of file to be written. Either 's', for a scalar
   * load-balance file, or 'p' for a parallel file.
   * To be used with Nemesis::ne_get_init_info().   
   */
  char ftype;
  
  /**
   * Containers for reading global sideset (boundary conditions) information.  Each vector will
   * eventually have num_side_sets_global entries, and be used in calls to
   * Nemesis::ne_get_ss_param_global().
   *
   * It's an error to call ne_get_ss_param_global when num_side_sets_global==0
   */
  std::vector<int> global_sideset_ids;
  std::vector<int> num_global_side_counts;
  std::vector<int> num_global_side_df_counts;


  /**
   * Containers for reading global nodeset information.  One vector entry per nodeset.
   * Each vector will eventually have num_node_sets_global entries, and
   * will be used in calls to Nemesis::ne_get_ns_param_global().
   *
   * It's an error to call ne_get_ns_param_global when num_node_sets_global==0
   */
  std::vector<int> global_nodeset_ids;
  std::vector<int> num_global_node_counts;
  std::vector<int> num_global_node_df_counts;


  /**
   * Read the global element block IDs and counts.  These vectors will
   * eventually have num_elem_blks_global entries.  To be used with
   * Nemesis::ne_get_eb_info_global().
   */
  std::vector<int> global_elem_blk_ids;
  std::vector<int> global_elem_blk_cnts;
  

  /**
   * To be used with the Nemesis::ne_get_loadbal_param() routine.
   */

  /**
   * The number of FEM nodes contained in FEM elements wholly owned by the current processor.
   * To be used with the Nemesis::ne_get_loadbal_param() routine.
   */
  int num_internal_nodes;

  /**
   * The number of FEM nodes local to a processor but residing in an
   * element which also has FEM nodes on other processors.
   * To be used with the Nemesis::ne_get_loadbal_param() routine.
   */
  int num_border_nodes;
  
  /**
   * The number of FEM nodes that reside on another processor but
   * whose element partially resides on the current processor.
   * To be used with the Nemesis::ne_get_loadbal_param() routine.
   */
  int num_external_nodes;
  
  /**
   * The number of internal FEM elements. Elements local to this processor.
   * To be used with the Nemesis::ne_get_loadbal_param() routine.
   */
  int num_internal_elems;

  /**
   * The number of border FEM elements. Elements local to this
   * processor but whose FEM nodes reside on other processors as well.
   * To be used with the Nemesis::ne_get_loadbal_param() routine.
   */
  int num_border_elems;
  
  /**
   * The number of nodal communication maps for this processor. (One
   * per neighboring proc?)
   * To be used with the Nemesis::ne_get_loadbal_param() routine.
   */
  int num_node_cmaps;
  
  /**
   * The number of elemental communication maps for this
   * processor. (One per neighboring proc?)
   * To be used with the Nemesis::ne_get_loadbal_param() routine.
   */
  int num_elem_cmaps;


  
  /**
   * Vector which stores internal element IDs.  Will have length
   * num_internal_elems.
   * To be used with Nemesis::ne_get_elem_map().
   */
  std::vector<int> elem_mapi;
  
  /**
   * Vector which stores border element IDs.  Will have length
   * num_border_elems.
   * To be used with Nemesis::ne_get_elem_map().
   */
  std::vector<int> elem_mapb;


  
  /**
   * Vector which stores internal node IDs.  Will have length
   * num_internal_nodes.
   * To be used with Nemesis::ne_get_node_map().
   */
  std::vector<int> node_mapi;
  
  /**
   * Vector which stores border node IDs.  Will have length
   * num_border_nodes.
   * To be used with Nemesis::ne_get_node_map().
   */
  std::vector<int> node_mapb;

  /**
   * Vector which stores external node IDs.  Will have length
   * num_external_nodes.
   * To be used with Nemesis::ne_get_node_map().
   */
  std::vector<int> node_mape;
  
  
  /**
   * Vectors for storing the communication map parameters.
   * Each will eventually have length num_node_cmaps OR
   * num_elem_cmaps as appropriate.
   * For use with Nemesis::ne_get_cmap_params().
   */
  std::vector<int> node_cmap_ids;
  std::vector<int> node_cmap_node_cnts;
  std::vector<int> elem_cmap_ids;
  std::vector<int> elem_cmap_elem_cnts;


  /**
   * 2 vectors of vectors for storing the node communication IDs for this processor.
   * There will be num_node_cmaps rows, row i will have node_cmap_node_cnts[i] entries.
   * To be used with Nemesis::ne_get_node_cmap().
   *
   * Remark: node_cmap_proc_ids is a vector, all entries of which are = node_cmap_ids[i]
   * Not sure what the point of that is...
   */
  std::vector<std::vector<int> > node_cmap_node_ids;
  std::vector<std::vector<int> > node_cmap_proc_ids;


  /**
   * 3 vectors of vectors for storing element communication IDs for this processor.
   * There will be num_elem_cmaps rows, row i will have elem_cmap_elem_cnts[i] entries.
   * To be used with Nemesis::ne_get_elem_cmap().
   */
  std::vector<std::vector<int> > elem_cmap_elem_ids;
  std::vector<std::vector<int> > elem_cmap_side_ids;
  std::vector<std::vector<int> > elem_cmap_proc_ids;


private:
  bool _verbose;
  
};

#endif // #if defined(HAVE_NEMESIS_API) && defined(HAVE_EXODUS_API)
#endif // #ifndef __nemesis_io_helper_h__
