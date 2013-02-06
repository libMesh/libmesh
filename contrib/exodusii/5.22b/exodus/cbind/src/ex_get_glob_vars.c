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
* exggv - ex_get_glob_vars
*
* entry conditions - 
*   input parameters:
*       int     exoid                   exodus file id
*       int     time_step               time step number
*       int     num_glob_vars           number of global vars in file
*
* exit conditions - 
*       float*  glob_var_vals           array of global variable values
*
* revision history - 
*
*
*****************************************************************************/

#include "exodusII.h"
#include "exodusII_int.h"

/*!
\ingroup ResultsData 

The function ex_get_glob_vars() reads the values of all the
global variables for a single time step. Memory must be allocated for
the global variables values array before this function is invoked.

Because global variables are floating point values, the application
code must declare the array passed to be the appropriate type
(float or double) to match the compute word size passed in
ex_create() or ex_open().

In case of an error, ex_get_glob_vars() returns a negative
number; a warning will return a positive number. Possible causes of
errors include:

 - data file not properly opened with call to ex_create() or ex_open()
 - no global variables stored in the file.
 - a warning value is returned if no global variables are stored in the file.

\param[in] exoid        exodus file ID returned from a previous call to ex_create() or ex_open().

\param[in] time_step    The time step, as described under ex_put_time(), at 
                        which the global variable values are desired. This is essentially 
                        an index (in the time dimension) into the global variable values 
                        array stored in the database. The first time step is 1.

\param[in]  num_glob_vars The number of global variables stored in the database.
\param[out] glob_var_vals Returned array of num_glob_vars global variable values 
                          for the time_step'th time step.

The following is an example code segment that reads all the global
variables at one time step:

\verbatim
int num_glo_vars, error, time_step;
float *var_values;

error = ex_get_variable_param (idexo, EX_GLOBAL, &num_glo_vars);
var_values = (float *) calloc (num_glo_vars, sizeof(float));
error = ex_get_glob_vars (idexo, time_step, num_glo_vars, 
                          var_values);
\endverbatim
 */

int ex_get_glob_vars (int   exoid,
                      int   time_step,
                      int   num_glob_vars,
                      void *glob_var_vals)
{
   int varid;
   int status;
   size_t start[2], count[2];
   char errmsg[MAX_ERR_LENGTH];

   exerrval = 0; /* clear error code */

   /* inquire previously defined variable */
   if ((status = nc_inq_varid (exoid, VAR_GLO_VAR, &varid)) != NC_NOERR) {
     exerrval = status;
     sprintf(errmsg,
            "Warning: failed to locate global variables in file id %d",
            exoid);
     ex_err("ex_get_glob_vars",errmsg,exerrval);
     return (EX_WARN);
   }

   /* read values of global variables */
   start[0] = --time_step;
   start[1] = 0;

   count[0] = 1;
   count[1] = num_glob_vars;

   if (ex_comp_ws(exoid) == 4) {
     status = nc_get_vara_float(exoid, varid, start, count, glob_var_vals);
   } else {
     status = nc_get_vara_double(exoid, varid, start, count, glob_var_vals);
   }

   if (status != NC_NOERR) {
     exerrval = status;
     sprintf(errmsg,
            "Error: failed to get global variable values from file id %d",
            exoid);
     ex_err("ex_get_glob_vars",errmsg,exerrval);
     return (EX_FATAL);
   }
   return (EX_NOERR);
}
