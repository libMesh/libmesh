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
#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/linear_solver.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"

using namespace libMesh;

void usage(const std::string & prog_name);

// Open the matrix and the right-hand-side vector named in command
// line arguments, solve them with a linear solver, and optionally
// write the solution to the given output filename.
int main(int argc, char ** argv)
{
  LibMeshInit init(argc, argv);

  if (argc < 3)
    usage(argv[0]);

  // Create a GetPot object to parse the command line
  GetPot command_line (argc, argv);

  std::string matfile, rhsfile;
  if (command_line.search(1, "-m"))
    matfile = command_line.next(matfile);
  else
    usage(argv[0]);

  if (command_line.search(1, "-r"))
    rhsfile = command_line.next(rhsfile);
  else
    usage(argv[0]);

  auto mat = SparseMatrix<Number>::build(init.comm());

  LOG_CALL("mat.read()", "main", mat->read(matfile));

  libMesh::out << "Loaded " << mat->m() << " by " << mat->n() <<
    " matrix " << matfile << std::endl;

  auto rhs = NumericVector<Number>::build(init.comm()),
       solution = NumericVector<Number>::build(init.comm());

  LOG_CALL("rhs.read()", "main", rhs->read_matlab(rhsfile));

  libMesh::out << "Loaded length-" << rhs->size() <<
    " vector " << rhsfile << std::endl;

  solution->init(rhs->size(), rhs->local_size());

  auto solver = LinearSolver<Number>::build(init.comm());

  Real tol = 1e-9;
  if (command_line.search(1, "-t"))
    tol = command_line.next(tol);

  unsigned int n_iters = 1e3;
  if (command_line.search(1, "-n"))
    n_iters = command_line.next(n_iters);

  auto [iters, res] = solver->solve(*mat, *solution, *rhs, tol, n_iters);

  libMesh::out << "Finished in " << iters << " iterations with " << res
      << " residual " << std::endl;

  return 0;
}


void usage(const std::string & prog_name)
{
  std::ostringstream helpList;
  helpList <<
    "Usage: " << prog_name <<
    " [options]\n"
    "\n"
    "options:\n"
    "    -m <string> Matrix filename (required)\n"
    "    -r <string> RHS Vector filename (required)\n"
    "    -t <tol>    Linear Solver tolerance (default 1e-9)\n"
    "    -n <n_iter> Linear Solver max iterations (default 1e3)\n";

  libMesh::out << helpList.str();
  exit(0);
}
