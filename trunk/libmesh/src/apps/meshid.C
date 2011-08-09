// Open the mesh named in command line arguments,
// update ids of one or more of the following enties
// as perscribed on the command line:
// blocks, sidesets and nodesets

#include <iostream>

#include "exodusII.h"
#include "exodusII_int.h"
#include "getpot.h"

#define BLOCKS  0x4
#define SIDES   0x2
#define NODES   0x1

void handle_error(int error)
{
  std::cout << "An error occured while working with the netCDF API" << std::endl;

  exit(1);
}

void usage_error(const char *progname)
{
  std::cout << "Usage: " << progname
            << " --input inputmesh --oldid <n> --newid <n> [--nodesetonly | --sidesetonly | --blockonly]" 
            << std::endl;
  exit(1);
}

int main(int argc, char** argv)
{
  GetPot cl(argc, argv);
  
  // Command line parsing
  if(!cl.search("--input"))
  {
    std::cerr << "No --input argument found!" << std::endl;
    usage_error(argv[0]);
  } 
  const char* meshname = cl.next("");

  if(!cl.search("--oldid"))
  {
    std::cerr << "No --oldid argument found!" << std::endl;
    usage_error(argv[0]);
  }
  long oldid = cl.next(0);

  if(!cl.search("--newid"))
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
  if (!flags)
    flags = NODES | SIDES | BLOCKS; // ALL ON

  // flags are exclusive
  if (flags != NODES && flags != SIDES && flags != BLOCKS  && flags != (NODES | SIDES | BLOCKS))
  {
    std::cerr << "Only one of the following options may be active [--nodesetonly | --sidesetonly | --blockonly]!" << std::endl;
    usage_error(argv[0]);
  }
  
  // Processing
  std::string var_name, dim_name;
  int status;
  int nc_id, var_id, dim_id, var_len;
  size_t dim_len;

  status = nc_open (meshname, NC_WRITE, &nc_id);
  if (status != NC_NOERR) handle_error(status);

  for (unsigned char mask = 1<<2; mask; mask>>=1)
  {
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
    default:
      continue;
    }

    // Get the length of the "variable" in question - stored in a dimension field
    status = nc_inq_dimid (nc_id, dim_name.c_str(), &dim_id);
    if (status != NC_NOERR) handle_error(status);

    status = nc_inq_dimlen (nc_id, dim_id, &dim_len);
    if (status != NC_NOERR) handle_error(status);
    
    // Now get the variable values themselves
    std::vector<long> var_vals(dim_len);
    
    status = nc_inq_varid (nc_id, var_name.c_str(), &var_id);
    if (status != NC_NOERR) handle_error(status);

    status = nc_get_var_long (nc_id, var_id, &var_vals[0]);
    if (status != NC_NOERR) handle_error(status);

    // Update the variable value specified on the command line
    for (unsigned int i=0; i<dim_len; ++i)
      if (var_vals[i] == oldid)
        var_vals[i] = newid;

    // Save that value back to the NetCDF database
    status = nc_put_var_long (nc_id, var_id, &var_vals[0]);
    if (status != NC_NOERR) handle_error(status);
  }

  // Write out the dataset
  status = nc_close(nc_id);
  
  exit(0);
}
