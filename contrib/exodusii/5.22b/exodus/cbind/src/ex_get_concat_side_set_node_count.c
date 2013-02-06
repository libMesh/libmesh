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
* exgcssc - ex_get_concat_side_set_node_count
*
* entry conditions - 
*   input parameters:
*       int     exoid                   exodus file id
*
* exit conditions - 
*       int     *side_set_node_cnt_list returned array of number of nodes for
*                                       side or face for all sidesets
* revision history - 
*
*****************************************************************************/

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "exodusII.h"
#include "exodusII_int.h"

/*! \cond INTERNAL */
static void *safe_free(void *array)
{
  if (array != 0) free(array);
  return 0;
}

/* Generic error message for element type/node count mapping...*/
#define EL_NODE_COUNT_ERROR sprintf(errmsg, \
      "Error: An element of type '%s' with %d nodes is not valid.",\
                      elem_blk_parms[i].elem_type,\
                      elem_blk_parms[i].num_nodes_per_elem);\
              ex_err("ex_get_side_set_node_count",errmsg,EX_MSG);\
              return(EX_FATAL);
/*! \endcond */

/*! \undoc */
int ex_get_concat_side_set_node_count(int exoid,
				      int *side_set_node_cnt_list)
{
  size_t m;
  int ii, i, j, iss, ioff; 
  ex_entity_id side_set_id;
  int  num_side_sets, num_elem_blks, num_df, ndim;
  int64_t tot_num_elem = 0, tot_num_ss_elem = 0, side, elem;
  void_int *elem_blk_ids       = NULL;
  void_int *side_set_ids       = NULL;
  void_int *ss_elem_ndx        = NULL;
  void_int *side_set_elem_list = NULL;
  void_int *side_set_side_list = NULL;
  size_t elem_ctr;
  int int_size, ids_size;
  int status;
  
  struct elem_blk_parm  *elem_blk_parms;

  char errmsg[MAX_ERR_LENGTH];

  exerrval = 0; /* clear error code */

  /* first check if any side sets are specified */
  /* inquire how many side sets have been stored */
  num_side_sets = ex_inquire_int(exoid, EX_INQ_SIDE_SETS);
  if (num_side_sets < 0) {
    sprintf(errmsg,
           "Error: failed to get number of side sets in file id %d",exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
    return(EX_FATAL);
  }

  if (num_side_sets == 0) {
    sprintf(errmsg,
           "Warning: no side sets defined in file id %d",exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,EX_WARN);
    return(EX_WARN);
  }

  num_elem_blks = ex_inquire_int(exoid, EX_INQ_ELEM_BLK);
  if (num_elem_blks < 0) {
    sprintf(errmsg,
           "Error: failed to get number of element blocks in file id %d",exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
    return(EX_FATAL);
  }

  tot_num_elem = ex_inquire_int(exoid, EX_INQ_ELEM);
  if (tot_num_elem < 0) {
    sprintf(errmsg,
           "Error: failed to get total number of elements in file id %d",exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
    return(EX_FATAL);
  }

  /* get the dimensionality of the coordinates;  this is necessary to
     distinguish between 2d TRIs and 3d TRIs */
  ndim = ex_inquire_int(exoid, EX_INQ_DIM);
  if (ndim < 0) {
    sprintf(errmsg,
           "Error: failed to get dimensionality in file id %d",exoid);
    ex_err("ex_cvt_nodes_to_sides",errmsg,exerrval);
    return(EX_FATAL);
  }

  int_size = sizeof(int);
  if (ex_int64_status(exoid) & EX_BULK_INT64_API)
    int_size = sizeof(int64_t);
  
  /* Allocate space for the element block ids */
  ids_size = sizeof(int);
  if (ex_int64_status(exoid) & EX_IDS_INT64_API) {
    ids_size = sizeof(int64_t);
  }

  if (!(elem_blk_ids=malloc(num_elem_blks*ids_size))) {
    exerrval = EX_MEMFAIL;
    sprintf(errmsg,
            "Error: failed to allocate space for element block ids for file id %d",
            exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
    goto error_ret;
  }

  if (ex_get_ids(exoid, EX_ELEM_BLOCK, elem_blk_ids) == -1) {
    sprintf(errmsg,
            "Error: failed to get element block ids in file id %d",
            exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,EX_MSG);
    return(EX_FATAL);
  } 

  /* Allocate space for the element block params */
  if (!(elem_blk_parms=malloc(num_elem_blks*sizeof(struct elem_blk_parm)))) {
    exerrval = EX_MEMFAIL;
    sprintf(errmsg,
      "Error: failed to allocate space for element block params for file id %d",
            exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
    goto error_ret;
  }

  elem_ctr = 0;
  for (i=0; i<num_elem_blks; i++)
  {
    ex_entity_id id;
    ex_block block;
    if (ex_int64_status(exoid) & EX_IDS_INT64_API) {
      id = ((int64_t*)elem_blk_ids)[i];
    } else {
      id = ((int*)elem_blk_ids)[i];
    }

    /* read in an element block parameter */
    block.type = EX_ELEM_BLOCK;
    block.id = id;

    /* read in an element block parameter */
    if ((ex_get_block_param (exoid, &block)) == -1) {
      sprintf(errmsg,
             "Error: failed to get element block  %"PRId64" parameters in file id %d",
              block.id, exoid);
      ex_err("ex_get_concat_side_set_node_count",errmsg,EX_MSG);
      return(EX_FATAL);
    }

    elem_blk_parms[i].num_elem_in_blk = block.num_entry;
    elem_blk_parms[i].num_nodes_per_elem = block.num_nodes_per_entry;
    elem_blk_parms[i].num_attr = block.num_attribute;

    for (m=0; m < strlen(block.topology); m++) {
      elem_blk_parms[i].elem_type[m] = toupper(block.topology[m]);
    }
    elem_blk_parms[i].elem_type[m] = '\0';

    if (strncmp(elem_blk_parms[i].elem_type,"CIRCLE",3) == 0)
    {
      elem_blk_parms[i].elem_type_val = EX_EL_CIRCLE;
      elem_blk_parms[i].num_sides = 1;
      elem_blk_parms[i].num_nodes_per_side[0] = 1;
    }
    else if (strncmp(elem_blk_parms[i].elem_type,"SPHERE",3) == 0)
    {
      elem_blk_parms[i].elem_type_val = EX_EL_SPHERE;
      elem_blk_parms[i].num_sides = 1;
        elem_blk_parms[i].num_nodes_per_side[0] = 1;
    }
    else if (strncmp(elem_blk_parms[i].elem_type,"QUAD",3) == 0)
    {
      elem_blk_parms[i].elem_type_val = EX_EL_QUAD;
      elem_blk_parms[i].num_sides = 4;
      if (elem_blk_parms[i].num_nodes_per_elem == 4) {
        elem_blk_parms[i].num_nodes_per_side[0] = 2;
        elem_blk_parms[i].num_nodes_per_side[1] = 2;
        elem_blk_parms[i].num_nodes_per_side[2] = 2;
        elem_blk_parms[i].num_nodes_per_side[3] = 2;
      } else if (elem_blk_parms[i].num_nodes_per_elem == 5) {
        elem_blk_parms[i].num_nodes_per_side[0] = 2;
        elem_blk_parms[i].num_nodes_per_side[1] = 2;
        elem_blk_parms[i].num_nodes_per_side[2] = 2;
        elem_blk_parms[i].num_nodes_per_side[3] = 2;
      } else if (elem_blk_parms[i].num_nodes_per_elem == 9 ||
                 elem_blk_parms[i].num_nodes_per_elem == 8) {
        elem_blk_parms[i].num_nodes_per_side[0] = 3;
        elem_blk_parms[i].num_nodes_per_side[1] = 3;
        elem_blk_parms[i].num_nodes_per_side[2] = 3;
        elem_blk_parms[i].num_nodes_per_side[3] = 3;
      } else {
        EL_NODE_COUNT_ERROR;
      }
    }
    else if (strncmp(elem_blk_parms[i].elem_type,"TRIANGLE",3) == 0)
    {
      elem_blk_parms[i].elem_type_val = EX_EL_TRIANGLE;
      if (ndim == 2) { /* 2d TRIs */
        elem_blk_parms[i].num_sides = 3;
        if (elem_blk_parms[i].num_nodes_per_elem == 3) {
          elem_blk_parms[i].num_nodes_per_side[0] = 2;
          elem_blk_parms[i].num_nodes_per_side[1] = 2;
          elem_blk_parms[i].num_nodes_per_side[2] = 2;
        } else if (elem_blk_parms[i].num_nodes_per_elem == 6) {
          elem_blk_parms[i].num_nodes_per_side[0] = 3;
          elem_blk_parms[i].num_nodes_per_side[1] = 3;
          elem_blk_parms[i].num_nodes_per_side[2] = 3;
        }
      } else if (ndim == 3) { /* 3d TRIs -- triangular shell*/
        elem_blk_parms[i].num_sides = 5; /* 2 Faces and 3 Edges */
        if (elem_blk_parms[i].num_nodes_per_elem == 3) {
          elem_blk_parms[i].num_nodes_per_side[0] = 3;
          elem_blk_parms[i].num_nodes_per_side[1] = 3;
          elem_blk_parms[i].num_nodes_per_side[2] = 2;
          elem_blk_parms[i].num_nodes_per_side[3] = 2;
          elem_blk_parms[i].num_nodes_per_side[4] = 2;
        } else if (elem_blk_parms[i].num_nodes_per_elem == 6) {
          elem_blk_parms[i].num_nodes_per_side[0] = 6;
          elem_blk_parms[i].num_nodes_per_side[1] = 6;
          elem_blk_parms[i].num_nodes_per_side[2] = 3;
          elem_blk_parms[i].num_nodes_per_side[3] = 3;
          elem_blk_parms[i].num_nodes_per_side[4] = 3;
        } else {
          EL_NODE_COUNT_ERROR;
        }
      }
    }
    else if (strncmp(elem_blk_parms[i].elem_type,"SHELL",3) == 0)
    {
      elem_blk_parms[i].elem_type_val = EX_EL_SHELL;

      if (elem_blk_parms[i].num_nodes_per_elem == 2) {/* KLUDGE for 2D Shells*/
        elem_blk_parms[i].num_sides = 2; 
        elem_blk_parms[i].num_nodes_per_side[0] = 2;
        elem_blk_parms[i].num_nodes_per_side[1] = 2;
      } else if (elem_blk_parms[i].num_nodes_per_elem == 4) {
        elem_blk_parms[i].num_sides = 6;  /* 2 Faces, 4 Edges */
        elem_blk_parms[i].num_nodes_per_side[0] = 4;
        elem_blk_parms[i].num_nodes_per_side[1] = 4;
        elem_blk_parms[i].num_nodes_per_side[2] = 2;
        elem_blk_parms[i].num_nodes_per_side[3] = 2;
        elem_blk_parms[i].num_nodes_per_side[4] = 2;
        elem_blk_parms[i].num_nodes_per_side[5] = 2;
      } else if (elem_blk_parms[i].num_nodes_per_elem == 8 ||
                 elem_blk_parms[i].num_nodes_per_elem == 9) {
        elem_blk_parms[i].num_sides = 6;  /* 2 Faces, 4 Edges */
        elem_blk_parms[i].num_nodes_per_side[0] =
          elem_blk_parms[i].num_nodes_per_elem; /* 8 or 9 */
        elem_blk_parms[i].num_nodes_per_side[1] =
          elem_blk_parms[i].num_nodes_per_elem; /* 8 or 9 */
        elem_blk_parms[i].num_nodes_per_side[2] = 3;
        elem_blk_parms[i].num_nodes_per_side[3] = 3;
        elem_blk_parms[i].num_nodes_per_side[4] = 3;
        elem_blk_parms[i].num_nodes_per_side[5] = 3;
      } else {
        EL_NODE_COUNT_ERROR;
      }
    }
    else if (strncmp(elem_blk_parms[i].elem_type,"HEX",3) == 0)
    {
      elem_blk_parms[i].elem_type_val = EX_EL_HEX;
      elem_blk_parms[i].num_sides = 6;  
      /* determine side set node stride */
      if (elem_blk_parms[i].num_nodes_per_elem == 8) {  /* 8-node bricks */
        elem_blk_parms[i].num_nodes_per_side[0] = 4;
        elem_blk_parms[i].num_nodes_per_side[1] = 4;
        elem_blk_parms[i].num_nodes_per_side[2] = 4;
        elem_blk_parms[i].num_nodes_per_side[3] = 4;
        elem_blk_parms[i].num_nodes_per_side[4] = 4;
        elem_blk_parms[i].num_nodes_per_side[5] = 4;
      } else if (elem_blk_parms[i].num_nodes_per_elem == 9) { /* 9-node bricks */
        elem_blk_parms[i].num_nodes_per_side[0] = 4;
        elem_blk_parms[i].num_nodes_per_side[1] = 4;
        elem_blk_parms[i].num_nodes_per_side[2] = 4;
        elem_blk_parms[i].num_nodes_per_side[3] = 4;
        elem_blk_parms[i].num_nodes_per_side[4] = 4;
        elem_blk_parms[i].num_nodes_per_side[5] = 4;
      } else if (elem_blk_parms[i].num_nodes_per_elem == 12) { /* HEXSHELLS */
        elem_blk_parms[i].num_nodes_per_side[0] = 6;
        elem_blk_parms[i].num_nodes_per_side[1] = 6;
        elem_blk_parms[i].num_nodes_per_side[2] = 6;
        elem_blk_parms[i].num_nodes_per_side[3] = 6;
        elem_blk_parms[i].num_nodes_per_side[4] = 4;
        elem_blk_parms[i].num_nodes_per_side[5] = 4;
      } else if (elem_blk_parms[i].num_nodes_per_elem == 20) { /* 20-node bricks */
        elem_blk_parms[i].num_nodes_per_side[0] = 8;
        elem_blk_parms[i].num_nodes_per_side[1] = 8;
        elem_blk_parms[i].num_nodes_per_side[2] = 8;
        elem_blk_parms[i].num_nodes_per_side[3] = 8;
        elem_blk_parms[i].num_nodes_per_side[4] = 8;
        elem_blk_parms[i].num_nodes_per_side[5] = 8;
      } else if (elem_blk_parms[i].num_nodes_per_elem == 27) { /* 27-node bricks */
        elem_blk_parms[i].num_nodes_per_side[0] = 9;
        elem_blk_parms[i].num_nodes_per_side[1] = 9;
        elem_blk_parms[i].num_nodes_per_side[2] = 9;
        elem_blk_parms[i].num_nodes_per_side[3] = 9;
        elem_blk_parms[i].num_nodes_per_side[4] = 9;
        elem_blk_parms[i].num_nodes_per_side[5] = 9;
      } else {
        EL_NODE_COUNT_ERROR;
      }
    }
    else if (strncmp(elem_blk_parms[i].elem_type,"TETRA",3) == 0)
    {
      elem_blk_parms[i].elem_type_val = EX_EL_TETRA;
      elem_blk_parms[i].num_sides = 4;  
      /* determine side set node stride */
      if (elem_blk_parms[i].num_nodes_per_elem == 4) {
        elem_blk_parms[i].num_nodes_per_side[0] = 3;
        elem_blk_parms[i].num_nodes_per_side[1] = 3;
        elem_blk_parms[i].num_nodes_per_side[2] = 3;
        elem_blk_parms[i].num_nodes_per_side[3] = 3;
      } else if (elem_blk_parms[i].num_nodes_per_elem == 8) {
        elem_blk_parms[i].num_nodes_per_side[0] = 4;
        elem_blk_parms[i].num_nodes_per_side[1] = 4;
        elem_blk_parms[i].num_nodes_per_side[2] = 4;
        elem_blk_parms[i].num_nodes_per_side[3] = 4;
      } else if (elem_blk_parms[i].num_nodes_per_elem == 10) {
        elem_blk_parms[i].num_nodes_per_side[0] = 6;
        elem_blk_parms[i].num_nodes_per_side[1] = 6;
        elem_blk_parms[i].num_nodes_per_side[2] = 6;
        elem_blk_parms[i].num_nodes_per_side[3] = 6;
      } else {
        EL_NODE_COUNT_ERROR;
      }
    }
    else if (strncmp(elem_blk_parms[i].elem_type,"WEDGE",3) == 0)
    {
      elem_blk_parms[i].elem_type_val = EX_EL_WEDGE;
      elem_blk_parms[i].num_sides = 5;  
      if (elem_blk_parms[i].num_nodes_per_elem == 6) {
        elem_blk_parms[i].num_nodes_per_side[0] = 4;
        elem_blk_parms[i].num_nodes_per_side[1] = 4;
        elem_blk_parms[i].num_nodes_per_side[2] = 4;
        elem_blk_parms[i].num_nodes_per_side[3] = 3;
        elem_blk_parms[i].num_nodes_per_side[4] = 3;
      } else if (elem_blk_parms[i].num_nodes_per_elem == 15){
        elem_blk_parms[i].num_nodes_per_side[0] = 8;
        elem_blk_parms[i].num_nodes_per_side[1] = 8;
        elem_blk_parms[i].num_nodes_per_side[2] = 8;
        elem_blk_parms[i].num_nodes_per_side[3] = 6;
        elem_blk_parms[i].num_nodes_per_side[4] = 6;
      } else {
        EL_NODE_COUNT_ERROR;
      }
    }
    else if (strncmp(elem_blk_parms[i].elem_type,"PYRAMID",3) == 0)
    {
      elem_blk_parms[i].elem_type_val = EX_EL_PYRAMID;
      elem_blk_parms[i].num_sides = 5;  
      if (elem_blk_parms[i].num_nodes_per_elem == 5) {
        elem_blk_parms[i].num_nodes_per_side[0] = 3;
        elem_blk_parms[i].num_nodes_per_side[1] = 3;
        elem_blk_parms[i].num_nodes_per_side[2] = 3;
        elem_blk_parms[i].num_nodes_per_side[3] = 3;
        elem_blk_parms[i].num_nodes_per_side[4] = 4;
      } else if (elem_blk_parms[i].num_nodes_per_elem == 13){
        elem_blk_parms[i].num_nodes_per_side[0] = 6;
        elem_blk_parms[i].num_nodes_per_side[1] = 6;
        elem_blk_parms[i].num_nodes_per_side[2] = 6;
        elem_blk_parms[i].num_nodes_per_side[3] = 6;
        elem_blk_parms[i].num_nodes_per_side[4] = 8;
      } else {
        EL_NODE_COUNT_ERROR;
      }
    }
    else if (strncmp(elem_blk_parms[i].elem_type,"BEAM",3) == 0)
    {
      elem_blk_parms[i].elem_type_val = EX_EL_BEAM;
      elem_blk_parms[i].num_sides = 2;  

      if (elem_blk_parms[i].num_nodes_per_elem == 2) {
        elem_blk_parms[i].num_nodes_per_side[0] = 2;
        elem_blk_parms[i].num_nodes_per_side[1] = 2;
      } else if (elem_blk_parms[i].num_nodes_per_elem == 3){
        elem_blk_parms[i].num_nodes_per_side[0] = 3;
        elem_blk_parms[i].num_nodes_per_side[1] = 3;
      } else {
        EL_NODE_COUNT_ERROR;
      }
    }
    else if ( (strncmp(elem_blk_parms[i].elem_type,"TRUSS",3) == 0) ||
              (strncmp(elem_blk_parms[i].elem_type,"BAR",3) == 0) ||
              (strncmp(elem_blk_parms[i].elem_type,"EDGE",3) == 0) )
    {
      elem_blk_parms[i].elem_type_val = EX_EL_TRUSS;
      elem_blk_parms[i].num_sides = 2;  

      if (elem_blk_parms[i].num_nodes_per_elem == 2) {
        elem_blk_parms[i].num_nodes_per_side[0] = 2;
        elem_blk_parms[i].num_nodes_per_side[1] = 2;
      } else if (elem_blk_parms[i].num_nodes_per_elem == 3) {
        elem_blk_parms[i].num_nodes_per_side[0] = 3;
        elem_blk_parms[i].num_nodes_per_side[1] = 3;
      } else {
        EL_NODE_COUNT_ERROR;
      }
    }
    /* Used for an empty block in a parallel decomposition */
    else if (strncmp(elem_blk_parms[i].elem_type,"NULL",3) == 0) {
      elem_blk_parms[i].elem_type_val = EX_EL_NULL_ELEMENT;
      elem_blk_parms[i].num_sides = 0;  
      elem_blk_parms[i].num_nodes_per_side[0] = 0;
      elem_blk_parms[i].num_elem_in_blk = 0;
    } else {
      /* unsupported element type; no problem if no sides specified for
         this element block */
      elem_blk_parms[i].elem_type_val = EX_EL_UNK;
      elem_blk_parms[i].num_sides = 0;  
      elem_blk_parms[i].num_nodes_per_side[0] = 0;
    }

    elem_blk_parms[i].elem_blk_id = id;    /* save id */
    elem_ctr += elem_blk_parms[i].num_elem_in_blk;
    elem_blk_parms[i].elem_ctr = elem_ctr;      /* save elem number max */
  }


  /* Finally... Create the list of node counts for each face in the
   * side set.
   */
  /* Allocate space for the sideset ids */
  if (!(side_set_ids=malloc(num_side_sets*ids_size))) {
    exerrval = EX_MEMFAIL;
    sprintf(errmsg,
            "Error: failed to allocate space for side set ids for file id %d",
            exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
    goto error_ret;
  }

  if (ex_get_ids(exoid, EX_SIDE_SET, side_set_ids) == -1) {
    sprintf(errmsg,
            "Error: failed to get side set ids in file id %d",
            exoid);
    ex_err("ex_get_concat_side_set_node_count",errmsg,EX_MSG);
    goto error_ret;
  } 

  /* Lookup index of side set id in VAR_SS_IDS array */
  ioff = 0;
  for (iss=0; iss<num_side_sets; iss++) {
    if (ex_int64_status(exoid) & EX_IDS_INT64_API) {
      side_set_id = ((int64_t*)side_set_ids)[iss];
    } else {
      side_set_id = ((int*)side_set_ids)[iss];
    }
    
    /* First determine the  # of elements in the side set*/
    if (int_size == sizeof(int64_t)) {
      status = ex_get_set_param(exoid,EX_SIDE_SET, side_set_id,&tot_num_ss_elem,&num_df);
    } else {
      int tot, df;
      status = ex_get_set_param(exoid,EX_SIDE_SET, side_set_id,&tot,&df);
      tot_num_ss_elem = tot;
      num_df = df;
    }

    if (status != EX_NOERR) {
      sprintf(errmsg,
	      "Error: failed to get number of elements in side set  %"PRId64" in file id %d",
	      side_set_id, exoid);
      ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
      goto error_ret;
    }
    
    if (tot_num_ss_elem == 0)
      continue;
   
    /* Allocate space for the side set element list */
    if (!(side_set_elem_list=malloc(tot_num_ss_elem*int_size))) {
      exerrval = EX_MEMFAIL;
      sprintf(errmsg,
	      "Error: failed to allocate space for side set element list for file id %d",
	      exoid);
      ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
      goto error_ret;
    }
    
    /* Allocate space for the side set side list */
    if (!(side_set_side_list=malloc(tot_num_ss_elem*int_size))) {
      exerrval = EX_MEMFAIL;
      sprintf(errmsg,
	      "Error: failed to allocate space for side set side list for file id %d",
	      exoid);
      ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
      goto error_ret;
    }
    
    if (ex_get_set(exoid, EX_SIDE_SET, side_set_id, 
		   side_set_elem_list, side_set_side_list) == -1) {
      sprintf(errmsg,
	      "Error: failed to get side set  %"PRId64" in file id %d",
	      side_set_id, exoid);
      ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
      goto error_ret;
    }
    
    /* Allocate space for the ss element index array */
    if (!(ss_elem_ndx=malloc(tot_num_ss_elem*int_size))) {
      exerrval = EX_MEMFAIL;
      sprintf(errmsg,
	      "Error: failed to allocate space for side set elem sort array for file id %d",
	      exoid);
      ex_err("ex_get_concat_side_set_node_count",errmsg,exerrval);
      goto error_ret;
    }
    
    /* Sort side set element list into index array  - non-destructive */
    if (int_size == sizeof(int64_t)) {
      /* Sort side set element list into index array  - non-destructive */
      int64_t *elems = (int64_t*)ss_elem_ndx;
      for (i=0;i<tot_num_ss_elem;i++) {
	elems[i] = i; /* init index array to current position */
      }
      ex_iqsort64(side_set_elem_list, ss_elem_ndx,tot_num_ss_elem);
    } else {
      /* Sort side set element list into index array  - non-destructive */
      int *elems = (int*)ss_elem_ndx;
      for (i=0;i<tot_num_ss_elem;i++) {
	elems[i] = i; /* init index array to current position */
      }
      ex_iqsort(side_set_elem_list, ss_elem_ndx,tot_num_ss_elem);
    }
    
    j = 0; /* The current element block... */
    for (ii=0;ii<tot_num_ss_elem;ii++) {
      
      int64_t elem_ndx;
      if (ex_int64_status(exoid) & EX_BULK_INT64_API) {
	elem_ndx = ((int64_t*)ss_elem_ndx)[ii];
	elem     = ((int64_t*)side_set_elem_list)[elem_ndx];
	side     = ((int64_t*)side_set_side_list)[elem_ndx]-1;
      } else {
	elem_ndx = ((int*)ss_elem_ndx)[ii];
	elem     = ((int*)side_set_elem_list)[elem_ndx];
	side     = ((int*)side_set_side_list)[elem_ndx]-1;
      }
      
      /*
       * Since the elements are being accessed in sorted, order, the
       * block that contains the elements must progress sequentially
       * from block 0 to block[num_elem_blks-1]. Once we find an element
       * not in this block, find a following block that contains it...
       */       
      for ( ; j<num_elem_blks; j++) {
	if (elem <= elem_blk_parms[j].elem_ctr) {
	  break;
	}
      }
      
      if (j < num_elem_blks) {
	assert(side < elem_blk_parms[j].num_sides);  
	side_set_node_cnt_list[elem_ndx+ioff] = elem_blk_parms[j].num_nodes_per_side[side];
      } else {
	exerrval = EX_BADPARAM;
	sprintf(errmsg,
		"Error: Invalid element number  %"PRId64" found in side set  %"PRId64" in file %d",
		elem, side_set_id, exoid);
	ex_err("ex_get_concat_side_set_node_count",errmsg,EX_MSG);
	goto error_ret;
      }
    }
    ss_elem_ndx        = safe_free(ss_elem_ndx);
    side_set_elem_list = safe_free(side_set_elem_list);
    side_set_side_list = safe_free(side_set_side_list);
    ioff += tot_num_ss_elem;
  }
    
  /* All done: release allocated memory */
  safe_free(elem_blk_ids);
  safe_free(side_set_ids);

  return(EX_NOERR);

 error_ret:
  safe_free(elem_blk_ids);
  safe_free(side_set_ids);
  safe_free(ss_elem_ndx);
  safe_free(side_set_elem_list);
  safe_free(side_set_side_list);
  return (EX_FATAL);
}
