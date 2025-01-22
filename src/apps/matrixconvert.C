// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// libmesh includes
#include "libmesh/libmesh.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/sparse_matrix.h"

using namespace libMesh;

// Open the matrix named in a command line argument,
// and rewrite it out in matlab format (our most portable currently)
// to the given output filename.
int main(int argc, char ** argv)
{
  LibMeshInit init(argc, argv);

  if (argc != 3)
    libmesh_error_msg
      ("Usage: " << argv[0] <<
       " inputmatrix outputmatrix");

  auto mat = SparseMatrix<Number>::build(init.comm());

  LOG_CALL("mat.read()", "main", mat->read(argv[1]));

  libMesh::out << "Loaded " << mat->m() << " by " << mat->n() <<
    " matrix " << argv[1] << std::endl;

  LOG_CALL("mat.print_matlab()", "main", mat->print_matlab(argv[2]));

  libMesh::out << "Wrote output " << argv[2] << std::endl;

  return 0;
}
