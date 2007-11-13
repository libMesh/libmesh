/*
 * Copyright (c) 2005 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
 * retains certain rights in this software.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.  
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

/*****************************************************************************
 *
 * exodusII.h - Exodus II include file, for general use
 *
 * author - Sandia National Laboratories
 *          
 * environment - UNIX
 *
 * exit conditions - 
 *
 * revision history - 
 *
 *  $Id: exodusII.h,v 1.5 2007/10/08 15:01:16 gdsjaar Exp $
 *****************************************************************************/

#include "netcdf.h"
#include "stddef.h"

#ifndef TRUE
#define TRUE -1
#endif

#ifndef FALSE
#define FALSE 0 
#endif

#ifndef EXODUS_II_HDR
#define EXODUS_II_HDR

/*
 * need following extern if this include file is used in a C++ program, to
 * keep the C++ compiler from mangling the function names.
 */
#ifdef __cplusplus
extern "C" {
#endif

  /*
   * The following are miscellaneous constants used in the EXODUS II API.
   */

#define EX_NOCLOBBER            0 /* Don't overwrite existing database, default */
#define EX_CLOBBER              1
#define EX_NORMAL_MODEL         2 /* disable mods that permit storage of larger models */
#define EX_LARGE_MODEL          4 /* enable mods that permit storage of larger models */
#define EX_NETCDF4              8 /* use the hdf5-based netcdf4 output */
#define EX_NOSHARE             16 /* Do not open netcdf file in "share" mode */
#define EX_SHARE               32 /* Do open netcdf file in "share" mode */

#define EX_READ                 0
#define EX_WRITE                1

#define EX_INQ_FILE_TYPE        1               /* inquire EXODUS II file type*/
#define EX_INQ_API_VERS         2               /* inquire API version number */
#define EX_INQ_DB_VERS          3               /* inquire database version   */
                                                /*   number                   */
#define EX_INQ_TITLE            4               /* inquire database title     */
#define EX_INQ_DIM              5               /* inquire number of          */
                                                /*   dimensions               */
#define EX_INQ_NODES            6               /* inquire number of nodes    */
#define EX_INQ_ELEM             7               /* inquire number of elements */
#define EX_INQ_ELEM_BLK         8               /* inquire number of element  */
                                                /*   blocks                   */
#define EX_INQ_NODE_SETS        9               /* inquire number of node sets*/
#define EX_INQ_NS_NODE_LEN      10              /* inquire length of node set */
                                                /*   node list                */
#define EX_INQ_SIDE_SETS        11              /* inquire number of side sets*/
#define EX_INQ_SS_NODE_LEN      12              /* inquire length of side set */
                                                /*   node list                */
#define EX_INQ_SS_ELEM_LEN      13              /* inquire length of side set */
                                                /*   element list             */
#define EX_INQ_QA               14              /* inquire number of QA       */
                                                /*   records                  */
#define EX_INQ_INFO             15              /* inquire number of info     */
                                                /*   records                  */
#define EX_INQ_TIME             16              /* inquire number of time     */
                                                /*   steps in the database    */
#define EX_INQ_EB_PROP          17              /* inquire number of element  */
                                                /*   block properties         */
#define EX_INQ_NS_PROP          18              /* inquire number of node set */
                                                /*   properties               */
#define EX_INQ_SS_PROP          19              /* inquire number of side set */
#define EX_INQ_NS_DF_LEN        20              /* inquire length of node set */
                                                /*   distribution factor  list*/
#define EX_INQ_SS_DF_LEN        21              /* inquire length of node set */
                                                /*   distribution factor  list*/
#define EX_INQ_LIB_VERS         22              /* inquire API Lib vers number*/
#define EX_INQ_EM_PROP          23              /* inquire number of element  */
                                                /*   map properties           */
#define EX_INQ_NM_PROP          24              /* inquire number of node     */
                                                /*   map properties           */
#define EX_INQ_ELEM_MAP         25              /* inquire number of element  */
                                                /*   maps                     */
#define EX_INQ_NODE_MAP         26              /* inquire number of node     */
                                                /*   maps                     */

  /*   properties               */
#define EX_ELEM_BLOCK           1               /* element block property code*/
#define EX_NODE_SET             2               /* node set property code     */
#define EX_SIDE_SET             3               /* side set property code     */
#define EX_ELEM_MAP             4               /* element map property code  */
#define EX_NODE_MAP             5               /* node map property code     */

  /*   max string lengths; constants that are used as netcdf dimensions must be
       of type long       */
#define MAX_STR_LENGTH          32L
#define MAX_VAR_NAME_LENGTH     20
#define MAX_LINE_LENGTH         80L
#define MAX_ERR_LENGTH          256

  /*   for netCDF 3.4, we estimate the size of the header; 
       if estimate is larger than this max, set the estimate to this max;
       I've never measured a header larger than 20K   */
#define MAX_HEADER_SIZE         30000

  /* routines for file initialization i/o */
  extern int ex_close (int exoid);
  extern int ex_cvt_nodes_to_sides(int exoid, int *num_elem_per_set,
				   int *num_nodes_per_set, int *side_sets_elem_index,
				   int *side_sets_node_index, int *side_sets_elem_list,
				   int *side_sets_node_list, int *side_sets_side_list);
  extern int ex_copy (int in_exoid, int out_exoid);
  extern int ex_create (const char *path, int cmode, int *comp_ws, int *io_ws);
  extern int ex_get_all_times (int   exoid, void *time_values);
  extern int ex_get_concat_node_sets (int   exoid,
				      int  *node_set_ids,
				      int  *num_nodes_per_set, 
				      int  *num_df_per_set, 
				      int  *node_sets_node_index,
				      int  *node_sets_df_index,
				      int  *node_sets_node_list, 
				      void *node_sets_dist_fact);
  extern int ex_get_coord_names (int    exoid,
				 char **coord_names);
  extern int ex_get_coord (int exoid,
			   void *x_coor,
			   void *y_coor,
			   void *z_coor);
  extern int ex_get_concat_side_sets (int   exoid,
				      int  *side_set_ids,
				      int  *num_elem_per_set,
				      int  *num_dist_per_set,
				      int  *side_sets_elem_index,
				      int  *side_sets_dist_index,
				      int  *side_sets_elem_list,
				      int  *side_sets_side_list,
				      void *side_sets_dist_fact);
  extern int ex_get_elem_attr_names (int   exoid,
				     int   elem_blk_id,
				     char **names);
  extern int ex_get_elem_attr (int   exoid,
			       int   elem_blk_id,
			       void *attrib);
  extern int ex_get_ids (int  exoid, int obj_type, int *ids);
  extern int ex_get_elem_blk_ids (int  exoid, int *ids);
  extern int ex_get_elem_block (int   exoid,
				int   elem_blk_id,
				char *elem_type,
				int  *num_elem_this_blk, 
				int  *num_nodes_per_elem,
				int  *num_attr);

  extern int ex_get_elem_conn (int   exoid,
			       int   elem_blk_id,
			       int  *connect);

  extern int ex_get_elem_map (int   exoid,
			      int   map_id,
			      int  *elem_map);
  extern int ex_get_elem_num_map (int  exoid,
				  int *elem_map);
  extern int ex_get_elem_var (int   exoid,
			      int   time_step,
			      int   elem_var_index,
			      int   elem_blk_id, 
			      int   num_elem_this_blk,
			      void *elem_var_vals);
  extern int ex_get_elem_varid (int  exoid,
				int *varid);
  extern int ex_get_elem_var_time (int   exoid,
				   int   elem_var_index,
				   int   elem_number,
				   int   beg_time_step, 
				   int   end_time_step,
				   void *elem_var_vals);
  extern int ex_get_coordinate_frames(int exoid, int *nframes, int *cf_ids,
				      void* pt_coordinates, char* tags);

  extern int ex_get_glob_vars (int   exoid,
			       int   time_step,
			       int   num_glob_vars,
			       void *glob_var_vals);

  extern int ex_get_glob_var_time (int   exoid,
				   int   glob_var_index,
				   int   beg_time_step,
				   int   end_time_step,
				   void *glob_var_vals);

  extern int ex_get_info (int exoid, char **info);

  extern int ex_get_init (int   exoid,
			  char *title,
			  int  *num_dim,
			  int  *num_nodes,
			  int  *num_elem, 
			  int  *num_elem_blk,
			  int  *num_node_sets,
			  int  *num_side_sets);

  extern int ex_get_map (int  exoid, int *elem_map);

  extern int ex_get_map_param (int   exoid,
			       int  *num_node_maps,
			       int  *num_elem_maps);

  extern int ex_get_name (int   exoid,
			  int   obj_type,
			  int   entity_id, 
			  char *name);

  extern int ex_get_names (int exoid,
			   int obj_type,
			   char **names);

  extern int ex_get_node_map (int   exoid,
			      int   map_id,
			      int  *node_map);

  extern int ex_get_node_num_map (int  exoid,
				  int *node_map);

  extern int ex_get_node_set_param (int  exoid,
				    int  node_set_id,
				    int *num_nodes_in_set,
				    int *num_df_in_set);

  extern int ex_get_node_set (int   exoid,
			      int   node_set_id,
			      int  *node_set_node_list);

  extern int ex_get_node_set_dist_fact  (int   exoid,
					 int   node_set_id,
					 void *node_set_dist_fact);

  extern int ex_get_node_set_ids (int  exoid,
				  int *ids);

  extern int ex_get_nset_var_tab (int  exoid,
				  int  num_nodesets,
				  int  num_nset_var,
				  int *nset_var_tab);

  extern int ex_get_nset_var (int   exoid,
			      int   time_step,
			      int   nset_var_index,
			      int   nset_id, 
			      int   num_node_this_nset,
			      void *nset_var_vals);

  extern int ex_get_nset_varid (int  exoid,
				int *varid);

  extern int ex_get_nodal_var (int   exoid,
			       int   time_step,
			       int   nodal_var_index,
			       int   num_nodes, 
			       void *nodal_var_vals);

  extern int ex_get_nodal_varid(int exoid, int *varid);

  extern int ex_get_nodal_var_time (int   exoid,
				    int   nodal_var_index,
				    int   node_number,
				    int   beg_time_step, 
				    int   end_time_step,
				    void *nodal_var_vals);

  extern int ex_get_nodal_varid_var(int   exoid,
				    int   time_step,
				    int   nodal_var_index,
				    int   num_nodes, 
				    int   varid,
				    void *nodal_var_vals);

  extern int ex_get_one_elem_attr (int   exoid,
				   int   elem_blk_id,
				   int   attrib_index,
				   void *attrib);

  extern int ex_get_prop_array (int   exoid,
				int   obj_type,
				const char *prop_name,
				int  *values);

  extern int ex_get_prop (int   exoid,
			  int   obj_type,
			  int   obj_id,
			  const char *prop_name,
			  int  *value);

  extern int ex_get_partial_elem_map (int   exoid,
				      int   map_id,
				      int ent_start,
				      int ent_count, 
				      int  *elem_map);

  extern int ex_get_prop_names (int    exoid,
				int    obj_type,
				char **prop_names);

  extern int ex_get_qa (int exoid,
			char *qa_record[][4]);
  extern int ex_get_side_set_node_list_len(int exoid,
					   int side_set_id,
					   int *side_set_node_list_len);
  extern int ex_get_side_set_param (int  exoid,
				    int  side_set_id,
				    int *num_side_in_set, 
				    int *num_dist_fact_in_set);
  extern int ex_get_side_set (int   exoid,
			      int   side_set_id,
			      int  *side_set_elem_list, 
			      int  *side_set_side_list);
  extern int ex_get_side_set_node_count(int exoid,
					int side_set_id,
					int *side_set_node_cnt_list);
  extern int ex_get_side_set_dist_fact (int   exoid,
					int   side_set_id,
					void *side_set_dist_fact);
  extern int ex_get_side_set_ids (int  exoid,
				  int *ids);
  extern int ex_get_side_set_node_list(int exoid,
				       int side_set_id,
				       int *side_set_node_cnt_list,
				       int *side_set_node_list);
  extern int ex_get_sset_var (int   exoid,
			      int   time_step,
			      int   sset_var_index,
			      int   sset_id, 
			      int   num_side_this_sset,
			      void *sset_var_vals);

  extern int ex_get_sset_var_tab (int  exoid,
				  int  num_sidesets,
				  int  num_sset_var,
				  int *sset_var_tab);
  extern int ex_get_sset_varid (int  exoid,
				int *varid);
  extern int ex_get_time (int   exoid,
			  int   time_step,
			  void *time_value);
  extern int ex_get_var_names (int   exoid,
			       const char *var_type,
			       int   num_vars,
			       char *var_names[]);
  extern int ex_get_varid (int  exoid, const char *var_type,
			   int *varid_arr);
  extern int ex_get_var_name (int   exoid,
			      const char *var_type,
			      int   var_num,
			      char *var_name);
  extern int ex_get_var_param (int   exoid,
			       const char *var_type,
			       int  *num_vars);

  extern int ex_get_object_truth_vector (int  exoid,
					 const char *var_type,
					 int  object_id,
					 int  num_var,
					 int *var_vector);
  
  extern int ex_get_var_tab (int  exoid,
			     const char *var_type,
			     int  num_blk,
			     int  num_var,
			     int *var_tab);
  
  extern int ex_get_elem_var_tab (int  exoid,
				  int  num_elem_blk,
				  int  num_elem_var,
				  int *elem_var_tab);
  extern int ex_open (const char  *path,
		      int    mode,
		      int   *comp_ws,
		      int   *io_ws,
		      float *version);

  extern int ex_put_attr_param (int   exoid,
				int   obj_type,
				int   obj_id,
				int   num_attrs);

  extern int ex_get_attr_param (int   exoid,
				int   obj_type,
				int   obj_id,
				int   *num_attrs);

  extern int ex_put_all_var_param (int exoid,
				   int num_g, int num_n,
				   int num_e, int *elem_var_tab,
				   int num_m, int *nset_var_tab,
				   int num_s, int *sset_var_tab);

  extern int ex_put_concat_elem_block (int    exoid,
				       const int*   elem_blk_id,
				       char *elem_type[],
				       const int*   num_elem_this_blk,
				       const int*   num_nodes_per_elem,
				       const int*   num_attr,
				       int    define_maps);

  extern int ex_put_concat_node_sets (int   exoid,
				      int  *node_set_ids,
				      int  *num_nodes_per_set,
				      int  *num_dist_per_set,
				      int  *node_sets_node_index,
				      int  *node_sets_df_index,
				      int  *node_sets_node_list,
				      void *node_sets_dist_fact);

  extern int ex_put_concat_side_sets (int   exoid,
				      int  *side_set_ids,
				      int  *num_elem_per_set,
				      int  *num_dist_per_set,
				      int  *side_sets_elem_index,
				      int  *side_sets_dist_index,
				      int  *side_sets_elem_list,
				      int  *side_sets_side_list,
				      void *side_sets_dist_fact);

  extern int ex_put_concat_var_param (int exoid, int num_g, int num_n,
				      int num_e, int num_elem_blk, int  *elem_var_tab);
  
  extern int ex_put_coord_names (int   exoid,
				 char *coord_names[]);
  extern int ex_put_coord (int   exoid,
			   const void *x_coor,
			   const void *y_coor,
			   const void *z_coor);
  extern int ex_put_elem_attr_names(int   exoid,
				    int   elem_blk_id,
				    char *names[]);
  extern int ex_put_elem_attr (int   exoid,
			       int   elem_blk_id,
			       const void *attrib);
  extern int ex_put_elem_block (int   exoid,
				int   elem_blk_id,
				const char *elem_type,
				int   num_elem_this_blk,
				int   num_nodes_per_elem,
				int   num_attr);

  extern int ex_put_elem_conn (int   exoid,
			       int   elem_blk_id,
			       const int  *connect);
  extern int ex_put_elem_map (int exoid,
			      int map_id,
			      const int *elem_map);
  extern int ex_put_elem_num_map (int  exoid,
				  const int *elem_map);
  extern int ex_put_elem_var (int   exoid,
			      int   time_step,
			      int   elem_var_index,
			      int   elem_blk_id,
			      int   num_elem_this_blk,
			      const void *elem_var_vals);

  extern int ex_put_coordinate_frames(int exoid, int nframes, const int cf_ids[], 
				      void* pt_coordinates, const char* tags);
  extern int ex_put_glob_vars (int   exoid,
			       int   time_step,
			       int   num_glob_vars,
			       const void *glob_var_vals);
  extern int ex_put_info (int   exoid, 
			  int   num_info,
			  char *info[]);
  extern int ex_put_init (int   exoid,
			  const char *title,
			  int   num_dim,
			  int   num_nodes,
			  int   num_elem,
			  int   num_elem_blk,
			  int   num_node_sets,
			  int   num_side_sets);

  extern int ex_put_map (int  exoid,
			 const int *elem_map);
  extern int ex_put_map_param (int   exoid,
			       int   num_node_maps,
			       int   num_elem_maps);
  extern int ex_put_name (int   exoid,
			  int   obj_type,
			  int   entity_id,
			  const char *name);
  extern int ex_put_names (int   exoid,
			   int   obj_type,
			   char *names[]);
  extern int ex_put_nodal_var (int   exoid,
			       int   time_step,
			       int   nodal_var_index,
			       int   num_nodes, 
			       const void *nodal_var_vals);

  extern int ex_put_nodal_varid_var(int   exoid,
				    int   time_step,
				    int   nodal_var_index,
				    int   num_nodes, 
				    int   varid,
				    const void *nodal_var_vals);

  extern int ex_put_node_map (int exoid,
			      int map_id,
			      const int *node_map);
  extern int ex_put_node_num_map (int  exoid,
				  const int *node_map);
  extern int ex_put_node_set_param (int exoid,
				    int node_set_id,
				    int num_nodes_in_set,
				    int num_dist_in_set);
  extern int ex_put_node_set (int   exoid,
			      int   node_set_id,
			      const int  *node_set_node_list);
  extern int ex_put_node_set_dist_fact  (int   exoid,
					 int   node_set_id,
					 const void *node_set_dist_fact);
  extern int ex_put_nset_var (int   exoid,
			      int   time_step,
			      int   nset_var_index,
			      int   nset_id,
			      int   num_nodes_this_nset,
			      const void *nset_var_vals);

  extern int ex_put_nset_var_tab (int  exoid,
				  int  num_nset,
				  int  num_nset_var,
				  int *nset_var_tab);
  extern int ex_put_one_elem_attr (int   exoid,
				   int   elem_blk_id,
				   int   attrib_index,
				   const void *attrib);
  extern int ex_put_partial_elem_map (int   exoid,
				      int   map_id,
				      int ent_start,
				      int ent_count, 
				      const int  *elem_map);

  extern int ex_put_prop (int   exoid,
			  int   obj_type,
			  int   obj_id,
			  const char *prop_name,
			  int   value);

  extern int ex_put_prop_array (int   exoid,
				int   obj_type,
				const char *prop_name,
				const int  *values);
  extern int ex_put_prop_names (int   exoid,
				int   obj_type,
				int   num_props,
				char **prop_names);
  extern int ex_put_qa (int   exoid,
			int   num_qa_records,
			char* qa_record[][4]);
  extern int ex_put_side_set_param (int exoid,
				    int side_set_id,
				    int num_side_in_set,
				    int num_dist_fact_in_set);
  extern int ex_put_side_set (int   exoid,
			      int   side_set_id,
			      const int  *side_set_elem_list,
			      const int  *side_set_side_list);
  extern int ex_put_side_set_dist_fact (int   exoid,
					int   side_set_id,
					const void *side_set_dist_fact);
  extern int ex_put_sset_var (int   exoid,
			      int   time_step,
			      int   sset_var_index,
			      int   sset_id,
			      int   num_faces_this_sset,
			      const void *sset_var_vals);

  extern int ex_put_sset_var_tab (int  exoid,
				  int  num_sset,
				  int  num_sset_var,
				  int *sset_var_tab);
  extern int ex_put_time (int   exoid,
			  int   time_step,
			  const void *time_value);
  extern int ex_put_varid_var(int   exoid,
			      int   time_step,
			      int   varid,
			      int   num_entity,
			      const void *var_vals);

  extern int ex_put_var_names (int   exoid,
			       const char *var_type,
			       int   num_vars,
			       char *var_names[]);
  extern int ex_put_var_name (int   exoid,
			      const char *var_type,
			      int   var_num,
			      const char *var_name);
  extern int ex_put_var_param (int   exoid,
			       const char *var_type,
			       int   num_vars);
  extern int ex_put_var_tab (int  exoid,
			     const char *var_type,
			     int  num_blk,
			     int  num_var,
			     int *var_tab);
  
  extern int ex_put_elem_var_tab (int  exoid,
				  int  num_elem_blk,
				  int  num_elem_var,
				  int *elem_var_tab);
  extern int ex_update (int exoid);
  extern int ex_get_num_props (int exoid, int obj_type);
  extern int ex_large_model(int exoid);
  extern size_t ex_header_size(int exoid);

  extern void ex_err(const char*, const char*, int);
  extern void ex_opts(int);
  extern int ex_inquire(int, int, int*, void*, char*);

  extern int ex_get_varid_var(int   exoid,
			      int   time_step,
			      int   varid,
			      int   num_entity,
			      void *var_vals);
  
  /* ERROR CODE DEFINITIONS AND STORAGE                                       */
  extern int exerrval;            /* shared error return value                */
  extern int exoptval;            /* error reporting flag (default is quiet)  */
  
#ifdef __cplusplus
}                               /* close brackets on extern "C" declaration */
#endif

#endif

/* ex_opts function codes - codes are OR'ed into exopts                     */
#define EX_VERBOSE      1       /* verbose mode message flag                */
#define EX_DEBUG        2       /* debug mode def                           */
#define EX_ABORT        4       /* abort mode flag def                      */

/* Exodus error return codes - exerrval return values:                      */
#define EX_MEMFAIL       1000   /* memory allocation failure flag def       */
#define EX_BADFILEMODE   1001   /* bad file mode def                        */
#define EX_BADFILEID     1002   /* bad file id def                          */
#define EX_WRONGFILETYPE 1003   /* wrong file type for function             */
#define EX_LOOKUPFAIL    1004   /* id table lookup failed                   */
#define EX_BADPARAM      1005   /* bad parameter passed                     */
#define EX_NULLENTITY   -1006   /* null entity found                        */
#define EX_MSG          -1000   /* message print code - no error implied    */
#define EX_PRTLASTMSG   -1001   /* print last error message msg code        */

#include "exodusII_ext.h"
