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
