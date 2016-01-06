// Open the getpot input file given by the input file name; write out
// all GetPot object data to the output file name

#include "libmesh/libmesh_config.h"
#include "libmesh/parsed_function.h"
#include "libmesh/point.h"

#include <cstdlib>

int main(int argc, char ** argv)
{
  using namespace libMesh;

  if (argc < 2)
    libmesh_error_msg("Usage: " << argv[0] << " function_to_eval [x] [y] [z] [t]");

  std::string function_string = argv[1];

  ParsedFunction<> func(function_string);

  const Point p ( (argc > 2) ? std::atof(argv[2]) : 0.0,
                  (argc > 3) ? std::atof(argv[3]) : 0.0,
                  (argc > 4) ? std::atof(argv[4]) : 0.0 );

  const libMesh::Real t = (argc > 5) ? std::atof(argv[5]) : 0.0;

  const libMesh::Number out = func(p,t);

  libMesh::out << "out = " << out << std::endl;

  return 0;
}
