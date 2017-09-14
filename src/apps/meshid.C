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


// Open the mesh named in command line arguments,
// update ids of one or more of the following entries
// as prescribed on the command line:
// blocks, sidesets and nodesets

#include <iostream>
#include <time.h> // time
#include <stdlib.h> // rand, srand

#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh_nullptr.h"

#ifdef LIBMESH_HAVE_EXODUS_API

#include "exodusII.h"
#include "exodusII_int.h"
#include "libmesh/getpot.h"

#define EXODUS_DIM 0x8
#define BLOCKS     0x4
#define SIDES      0x2
#define NODES      0x1

// Report error from NetCDF and exit
void handle_error(int error, std::string message);

// Report error in input flags and exit
void usage_error(const char * progname);

// Generate a random string, seed RNG with current time.
void gen_random_string(std::string & s, const int len);



int main(int argc, char ** argv)
{
  GetPot cl(argc, argv);

  // Command line parsing
  if (!cl.search("--input"))
    {
      std::cerr << "No --input argument found!" << std::endl;
      usage_error(argv[0]);
    }
  const char * meshname = cl.next("");

  if (!cl.search("--oldid"))
    {
      std::cerr << "No --oldid argument found!" << std::endl;
      usage_error(argv[0]);
    }
  long oldid = cl.next(0);

  if (!cl.search("--newid"))
    {
      std::cerr << "No --newid argument found!" << std::endl;
      usage_error(argv[0]);
    }
  long newid = cl.next(0);

  unsigned char flags = 0;
  if (cl.search("--nodesetonly"))
    flags |= NODES;
  if (cl.search("--sidesetonly"))
    flags |= SIDES;
  if (cl.search("--blockonly"))
    flags |= BLOCKS;
  if (cl.search("--dim"))
    flags |= EXODUS_DIM;

  // No command line flags were set, turn on NODES, SIDES, and BLOCKS
  if (!flags)
    flags = NODES | SIDES | BLOCKS; // ALL except EXODUS_DIM on

  // flags are exclusive
  if (flags != NODES &&
      flags != SIDES &&
      flags != BLOCKS &&
      flags != EXODUS_DIM &&
      flags != (NODES | SIDES | BLOCKS))
    {
      std::cerr << "Only one of the following options may be active [--nodesetonly | --sidesetonly | --blockonly | --dim]!" << std::endl;
      usage_error(argv[0]);
    }

  // Processing
  std::string var_name, dim_name;
  int status;
  int nc_id, var_id, dim_id;
  size_t dim_len;

  status = nc_open (meshname, NC_WRITE, &nc_id);
  if (status != NC_NOERR) handle_error(status, "Error while opening file.");

  for (unsigned char mask = 8; mask; mask/=2)
    {
      // These are char *'s #defined in exodusII_int.h
      switch (flags & mask)
        {
        case BLOCKS:
          dim_name = DIM_NUM_EL_BLK;
          var_name = VAR_ID_EL_BLK;
          break;
        case SIDES:
          dim_name = DIM_NUM_SS;
          var_name = VAR_SS_IDS;
          break;
        case NODES:
          dim_name = DIM_NUM_NS;
          var_name = VAR_NS_IDS;
          break;
        case EXODUS_DIM:
          dim_name = DIM_NUM_DIM;
          // var_name not used for setting dimension
          break;
        default:
          // We don't match this flag, so go to the next mask
          continue;
        }

      // Get the ID and length of the variable in question - stored in a dimension field
      status = nc_inq_dimid (nc_id, dim_name.c_str(), &dim_id);
      if (status != NC_NOERR) handle_error(status, "Error while inquiring about a dimension's ID.");

      status = nc_inq_dimlen (nc_id, dim_id, &dim_len);
      if (status != NC_NOERR) handle_error(status, "Error while inquiring about a dimension's length.");

      if ((flags & mask) != EXODUS_DIM)
        {
          // Now get the variable values themselves
          std::vector<long> var_vals(dim_len);

          status = nc_inq_varid (nc_id, var_name.c_str(), &var_id);
          if (status != NC_NOERR) handle_error(status, "Error while inquiring about a variable's ID.");

          status = nc_get_var_long (nc_id, var_id, &var_vals[0]);
          if (status != NC_NOERR) handle_error(status, "Error while retrieving a variable's values.");

          // Update the variable value specified on the command line
          for (unsigned int i=0; i<dim_len; ++i)
            if (var_vals[i] == oldid)
              var_vals[i] = newid;

          // Save that value back to the NetCDF database
          status = nc_put_var_long (nc_id, var_id, &var_vals[0]);
          if (status != NC_NOERR) handle_error(status, "Error while writing a variable's values.");
        }

      // Redefine the dimension
      else
        {
          // The value stored in dim_len is actually the dimension?
          if (dim_len == (size_t)oldid)
            {
              // Trying to change def's always raises
              // Error -38: /* Operation not allowed in data mode */
              // unless you are in "define" mode.  So let's go there now.

              // Try to put the file into define mode
              status = nc_redef(nc_id);
              if (status != NC_NOERR) handle_error(status, "Error while putting file into define mode.");

              // Rename the "num_dim" dimension.  Note: this will fail if there is already a dimension
              // which has the name you are trying to use.  This can happen, for example if you have already
              // changed the dimension of this exodus file once using this very script.  There appears
              // to be no way to delete a dimension using basic NetCDF interfaces, so to workaround this
              // we just rename it to an arbitrary unique string that Exodus will ignore.

              // Construct a string with 6 random alpha-numeric characters at the end.
              std::string random_dim_name;
              gen_random_string(random_dim_name, 6);
              random_dim_name = std::string("ignored_num_dim_") + random_dim_name;

              // Rename the old dimension variable to our randomly-chosen name
              status = nc_rename_dim(nc_id, dim_id, random_dim_name.c_str());
              if (status != NC_NOERR) handle_error(status, "Error while trying to rename num_dim.");

              // Now define a new "num_dim" value of newid
              int dummy=0;
              status = nc_def_dim (nc_id, dim_name.c_str(), newid, &dummy);
              if (status != NC_NOERR) handle_error(status, "Error while trying to define num_dim.");

              // Leave define mode
              status = nc_enddef(nc_id);
              if (status != NC_NOERR) handle_error(status, "Error while leaving define mode.");
            }
        }
    } // end for

  // Write out the dataset
  status = nc_close(nc_id);

  return (status != NC_NOERR);
}

#else // LIBMESH_HAVE_EXODUS_API

int main(int, char **)
{
  std::cerr << "Error: meshid requires libMesh configured with --enable-exodus" << std::endl;
}
#endif // LIBMESH_HAVE_EXODUS_API



void handle_error(int error, std::string message)
{
  std::cout << "Error " << error << " occurred while working with the netCDF API" << std::endl;
  std::cout << message << std::endl;

  exit(1);
}



void usage_error(const char * progname)
{
  std::cout << "Usage: " << progname
            << " --input inputmesh --oldid <n> --newid <n> [--nodesetonly | --sidesetonly | --blockonly]"
            << std::endl;
  exit(1);
}



void gen_random_string(std::string & s, const int len)
{
  static const char alphanum[] =
    "0123456789"
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz";

  // Seed the random number generator with the current time
  srand( static_cast<unsigned>(time(libmesh_nullptr)) );

  s.resize(len);
  for (int i = 0; i < len; ++i)
    {
      unsigned int alphai = static_cast<unsigned int>
        (rand() / (RAND_MAX+1.0) * (sizeof(alphanum)-1));
      s[i] = alphanum[alphai];
    }
}
