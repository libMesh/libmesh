// Open the two solution files named on standard input, find all
// variables they have in common, and output the Hilbert norms of the
// differences between them.

#include "libmesh.h"

#include "mesh.h"
#include "equation_systems.h"
#include "exact_solution.h"

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
  es1.read(argv[2], libMeshEnums::READ);
  std::cout << "Loaded coarse solution " << argv[2] << std::endl;
  mesh2.read(argv[3]);
  std::cout << "Loaded fine mesh " << argv[3] << std::endl;
  es2.read(argv[4], libMeshEnums::READ);
  std::cout << "Loaded fine solution " << argv[4] << std::endl;

  ExactSolution exact_sol(es1);
  exact_sol.attach_reference_solution(&es2);

  std::vector<std::string> sysnames;
  sysnames.reserve(es1.n_systems());

  for (unsigned int i = 0; i != es1.n_systems(); ++i)
    if (es2.has_system(es1.get_system(i).name()))
      sysnames.push_back(es1.get_system(i).name());

  for (unsigned int i = 0; i != sysnames.size(); ++i)
    {
      std::string sysname = sysnames[i];
      System &sys1 = es1.get_system(sysname);
      System &sys2 = es2.get_system(sysname);

      for (unsigned int j = 0; j != sys1.n_vars(); ++j)
        {
          std::string varname = sys1.variable_name(j);

          if (!sys2.has_variable(varname))
            continue;

          exact_sol.compute_error(sysname, varname);

          std::cout << "Errors in system " << sysname << ", variable " << varname << ":" << std::endl;
          std::cout << "L2 error: " << exact_sol.l2_error(sysname, varname) << std::endl;
          std::cout << "H1 error: " << exact_sol.h1_error(sysname, varname) << std::endl;
          std::cout << "H2 error: " << exact_sol.h2_error(sysname, varname) << std::endl;
        }
    }

  return 0;
}
