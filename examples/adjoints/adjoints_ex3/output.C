#include <fstream>
#include <unistd.h>

#include "libmesh/libmesh_common.h"

#include "output.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

void start_output(unsigned int timesteps,
                  std::string filename,
                  std::string varname)
{
  // Truncate then append if we're doing a restart
  if (timesteps)
    {
      std::ifstream infile (filename.c_str());
      char buf[1025];
      if (!infile.is_open())
        libmesh_error_msg("Error opening file " << filename);

      // Count off the lines we want to save, including
      // the original header
      for (unsigned int i=0; i != timesteps+1; ++i)
        infile.getline(buf, 1024);
      if (!infile.good())
        libmesh_error_msg("Error reading line from file " << filename);
      unsigned int length = cast_int<unsigned int>(infile.tellg());
      infile.close();

      // Then throw away the rest
      int err = truncate(filename.c_str(), length);
      if (err != 0)
        libmesh_error_msg("Error truncating file " << filename);
    }
  // Write new headers if we're starting a new simulation
  else
    {
      std::ofstream outfile (filename.c_str(), std::ios_base::out);
      outfile << varname << " = [" << std::endl;
    }
}
