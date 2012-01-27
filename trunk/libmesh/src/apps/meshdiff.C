// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// Open the two solution files named on standard input, find all
// variables they have in common, and output the Hilbert norms of the
// differences between them.

#include "libmesh.h"

#include "mesh.h"
#include "equation_systems.h"
#include "exact_solution.h"

using namespace libMesh;

unsigned int dim = 2; // This gets overridden by most mesh formats

int main(int argc, char** argv)
{
  LibMeshInit init(argc, argv);

  Mesh mesh1(dim), mesh2(dim);
  EquationSystems es1(mesh1), es2(mesh2);

  std::cout << "Usage: " << argv[0]
            << " coarsemesh coarsesolution finemesh finesolution" << std::endl;

  mesh1.read(argv[1]);
  std::cout << "Loaded coarse mesh " << argv[1] << std::endl;
  es1.read(argv[2]);
  std::cout << "Loaded coarse solution " << argv[2] << std::endl;
  mesh2.read(argv[3]);
  std::cout << "Loaded fine mesh " << argv[3] << std::endl;
  es2.read(argv[4]);
  std::cout << "Loaded fine solution " << argv[4] << std::endl;

  ExactSolution exact_sol(es1);
  exact_sol.attach_reference_solution(&es2);

  std::vector<std::string> sysnames;
  sysnames.reserve(es1.n_systems());

  for (unsigned int i = 0; i != es1.n_systems(); ++i)
    if (es2.has_system(es1.get_system(i).name()))
      sysnames.push_back(es1.get_system(i).name());
    else
      std::cout << "Coarse solution system "
                << es1.get_system(i).name()
		<< " not found in fine solution!" << std::endl;

  for (unsigned int i = 0; i != es2.n_systems(); ++i)
    if (!es1.has_system(es2.get_system(i).name()))
      std::cout << "Fine solution system "
                << es2.get_system(i).name()
		<< " not found in coarse solution!" << std::endl;

  if (!es1.n_systems() && !es2.n_systems())
    std::cout << "No systems found in fine or coarse solution!" 
              << std::endl;

  for (unsigned int i = 0; i != sysnames.size(); ++i)
    {
      const std::string sysname = sysnames[i];
      const System &sys1 = es1.get_system(sysname);
      const System &sys2 = es2.get_system(sysname);

      for (unsigned int j = 0; j != sys1.n_vars(); ++j)
        {
          const std::string varname = sys1.variable_name(j);

          if (!sys2.has_variable(varname))
            {
              std::cout << "Fine solution system " << sysname 
                        << " variable " << varname
		        << " not found in coarse solution!" << std::endl;
              continue;
            }
          const unsigned int j2 = sys2.variable_number(varname);

          exact_sol.compute_error(sysname, varname);

          std::cout << "Errors in system " << sysname << ", variable " << varname << ":" << std::endl;
          std::cout << "L2 error: " << exact_sol.l2_error(sysname, varname) 
                    << ", fine norm: " << sys1.calculate_norm(*sys1.solution, j, L2)
                    << ", coarse norm: " << sys2.calculate_norm(*sys2.solution, j2, L2) << std::endl;
          std::cout << "H1 error: " << exact_sol.h1_error(sysname, varname)
                    << ", fine norm: " << sys1.calculate_norm(*sys1.solution, j, H1)
                    << ", coarse norm: " << sys2.calculate_norm(*sys2.solution, j2, H1) << std::endl;
          std::cout << "H2 error: " << exact_sol.h2_error(sysname, varname)
                    << ", fine norm: " << sys1.calculate_norm(*sys1.solution, j, H2)
                    << ", coarse norm: " << sys2.calculate_norm(*sys2.solution, j2, H2) << std::endl;
        }

      for (unsigned int j = 0; j != sys2.n_vars(); ++j)
        {
          const std::string varname = sys2.variable_name(j);

          if (!sys1.has_variable(varname))
            {
              std::cout << "Coarse solution system " << sysname 
                        << " variable " << varname
		        << " not found in fine solution!" << std::endl;
              continue;
            }
        }

      if (!sys1.n_vars() && !sys2.n_vars())
        std::cout << "No variables found in fine or coarse solution system "
                  << sysname << '!' << std::endl;
    }

  return 0;
}
