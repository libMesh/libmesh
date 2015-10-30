// Open the getpot input file given by the input file name; write out
// all GetPot object data to the output file name

#include "libmesh/libmesh_config.h"
#include "libmesh/getpot.h"

#include <fstream>

int main(int argc, char** argv)
{
  using namespace libMesh;

  if (argc < 2)
    libmesh_error_msg("Usage: " << argv[0] << " inputconfigfile [outputconfigfile]");

  GetPot gp(argv[1]);

  std::ostream *my_out;
  std::ofstream fout;
  fout.exceptions ( std::ofstream::failbit | std::ofstream::badbit );

  if (argc < 3)
    my_out = &std::cout;
  else
    {
      fout.open(argv[2]);
      my_out = &fout;
    }

  gp.print("", *my_out, 1);
}
