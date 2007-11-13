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
* expmp - ex_get_attr_param
*
* entry conditions - 
*   input parameters:
*       int     exoid           exodus file id
*       int     obj_type        block/set type (node, edge, face, elem)
*       int     obj_id          block/set id (ignored for NODAL)       
*       int     num_attrs       number of attributes
*
* exit conditions - 
*
*
*****************************************************************************/

#include "exodusII.h"
#include "exodusII_int.h"

/*!
 * retrieves the number of attributes.
 */

int ex_get_attr_param (int   exoid,
		       int   obj_type,
		       int   obj_id,
		       int*  num_attrs)
{
  int dimid;
  
  char errmsg[MAX_ERR_LENGTH];
  const char *tname;
  const char *vobjids;
  const char *dnumobjatt;

  int obj_id_ndx;
  long lnum_attr_per_entry;
  
  switch (obj_type) {
  case EX_NODE_SET:
    tname = "node set";
    vobjids = VAR_NS_IDS;
    break;
  case EX_EDGE_SET:
    tname = "edge set";
    vobjids = VAR_ES_IDS;
    break;
  case EX_FACE_SET:
    tname = "face set";
    vobjids = VAR_FS_IDS;
    break;
  case EX_ELEM_SET:
    tname = "element set";
    vobjids = VAR_ELS_IDS;
    break;
  case EX_NODAL:
    tname = "node block";
    break;
  case EX_EDGE_BLOCK:
    tname = "edge block";
    vobjids = VAR_ID_ED_BLK;
    break;
  case EX_FACE_BLOCK:
    tname = "face block";
    vobjids = VAR_ID_FA_BLK;
    break;
  case EX_ELEM_BLOCK:
    tname = "element block";
    vobjids = VAR_ID_EL_BLK;
    break;
  default:
    exerrval = EX_BADPARAM;
    sprintf( errmsg, "Error: Invalid object type (%d) specified for file id %d",
	     obj_type, exoid );
    ex_err( "ex_get_attr", errmsg, exerrval );
    return (EX_FATAL);
  }

  /* Determine index of obj_id in vobjids array */
  if (obj_type == EX_NODAL)
    obj_id_ndx = 0;
  else {
    obj_id_ndx = ex_id_lkup(exoid,vobjids,obj_id);
    
    if (exerrval != 0) {
      if (exerrval == EX_NULLENTITY) {
	*num_attrs = 0;
	return (EX_NOERR);
      } else {
	sprintf(errmsg,
		"Warning: failed to locate %s id %d in %s array in file id %d",
		tname,obj_id,vobjids, exoid);
	ex_err("ex_get_attr_param",errmsg,exerrval);
	return (EX_WARN);
      }
    }
  }

  switch (obj_type) {
  case EX_NODE_SET:
    tname = "node set";
    dnumobjatt = DIM_NUM_ATT_IN_NS(obj_id_ndx);
    break;
  case EX_EDGE_SET:
    tname = "edge set";
    dnumobjatt = DIM_NUM_ATT_IN_ES(obj_id_ndx);
    break;
  case EX_FACE_SET:
    tname = "face set";
    dnumobjatt = DIM_NUM_ATT_IN_FS(obj_id_ndx);
    break;
  case EX_ELEM_SET:
    tname = "elem set";
    dnumobjatt = DIM_NUM_ATT_IN_ELS(obj_id_ndx);
    break;
  case EX_NODAL:
    tname = "node block";
    dnumobjatt = DIM_NUM_ATT_IN_NBLK;
    break;
  case EX_EDGE_BLOCK:
    tname = "edge block";
    dnumobjatt = DIM_NUM_ATT_IN_EBLK(obj_id_ndx);
    break;
  case EX_FACE_BLOCK:
    tname = "face block";
    dnumobjatt = DIM_NUM_ATT_IN_FBLK(obj_id_ndx);
    break;
  case EX_ELEM_BLOCK:
    tname = "element block";
    dnumobjatt = DIM_NUM_ATT_IN_BLK(obj_id_ndx);
    break;
  default:
    exerrval = EX_BADPARAM;
    sprintf(errmsg, "Error: Bad block type (%d) specified for file id %d",
	    obj_type, exoid );
    ex_err("ex_get_attr_param",errmsg,exerrval);
    return (EX_FATAL);
  }

  exerrval = 0; /* clear error code */

  if ((dimid = ncdimid (exoid, dnumobjatt)) == -1) {
    /* dimension is undefined */
    *num_attrs = 0;
  } else {
    if (ncdiminq (exoid, dimid, (char *)0, &lnum_attr_per_entry) == -1)
      {
	exerrval = ncerr;
	sprintf(errmsg,
		"Error: failed to get number of attributes in %s %d in file id %d",
		tname,obj_id, exoid);
	ex_err("ex_get_attr_param",errmsg, exerrval);
	return(EX_FATAL);
      }
    *num_attrs = lnum_attr_per_entry;
  }
  return (EX_NOERR);
}

