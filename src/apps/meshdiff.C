// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exact_solution.h"
#include "libmesh/mesh_function.h"
#include "libmesh/namebased_io.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/enum_norm_type.h"

using namespace libMesh;

unsigned char dim = 2; // This gets overridden by most mesh formats

int main(int argc, char ** argv)
{
  LibMeshInit init(argc, argv);

  Mesh coarse_mesh(init.comm(), dim), fine_mesh(init.comm(), dim);
  EquationSystems coarse_es(coarse_mesh), fine_es(fine_mesh);

  libMesh::out << "Usage: " << argv[0]
               << " coarsemesh coarsesolution finemesh finesolution [outputdiff]" << std::endl;

  if (argc < 5)
    libmesh_error();

  coarse_mesh.read(argv[1]);
  libMesh::out << "Loaded coarse mesh " << argv[1] << std::endl;
  coarse_es.read(argv[2]);
  libMesh::out << "Loaded coarse solution " << argv[2] << std::endl;
  fine_mesh.read(argv[3]);
  libMesh::out << "Loaded fine mesh " << argv[3] << std::endl;
  fine_es.read(argv[4]);
  libMesh::out << "Loaded fine solution " << argv[4] << std::endl;

  ExactSolution exact_sol(coarse_es);
  exact_sol.attach_reference_solution(&fine_es);

  std::vector<std::string> sysnames;
  sysnames.reserve(coarse_es.n_systems());

  for (unsigned int i = 0; i != coarse_es.n_systems(); ++i)
    if (fine_es.has_system(coarse_es.get_system(i).name()))
      sysnames.push_back(coarse_es.get_system(i).name());
    else
      libMesh::out << "Coarse solution system "
                   << coarse_es.get_system(i).name()
                   << " not found in fine solution!" << std::endl;

  for (unsigned int i = 0; i != fine_es.n_systems(); ++i)
    if (!coarse_es.has_system(fine_es.get_system(i).name()))
      libMesh::out << "Fine solution system "
                   << fine_es.get_system(i).name()
                   << " not found in coarse solution!" << std::endl;

  if (!coarse_es.n_systems() && !fine_es.n_systems())
    libMesh::out << "No systems found in fine or coarse solution!"
                 << std::endl;

  for (std::size_t i = 0; i != sysnames.size(); ++i)
    {
      const std::string sysname = sysnames[i];
      const System & coarse_sys = coarse_es.get_system(sysname);
      const System & fine_sys = fine_es.get_system(sysname);

      for (unsigned int j = 0; j != coarse_sys.n_vars(); ++j)
        {
          const std::string varname = coarse_sys.variable_name(j);

          if (!fine_sys.has_variable(varname))
            {
              libMesh::out << "Fine solution system " << sysname
                           << " variable " << varname
                           << " not found in coarse solution!" << std::endl;
              continue;
            }
          const unsigned int j2 = fine_sys.variable_number(varname);

          exact_sol.compute_error(sysname, varname);

          libMesh::out << "Errors in system " << sysname << ", variable " << varname << ":" << std::endl;
          libMesh::out << "L2 error: " << exact_sol.l2_error(sysname, varname)
                       << ", fine norm: " << coarse_sys.calculate_norm(*coarse_sys.solution, j, L2)
                       << ", coarse norm: " << fine_sys.calculate_norm(*fine_sys.solution, j2, L2) << std::endl;
          libMesh::out << "H1 error: " << exact_sol.h1_error(sysname, varname)
                       << ", fine norm: " << coarse_sys.calculate_norm(*coarse_sys.solution, j, H1)
                       << ", coarse norm: " << fine_sys.calculate_norm(*fine_sys.solution, j2, H1) << std::endl;
          libMesh::out << "H2 error: " << exact_sol.h2_error(sysname, varname)
                       << ", fine norm: " << coarse_sys.calculate_norm(*coarse_sys.solution, j, H2)
                       << ", coarse norm: " << fine_sys.calculate_norm(*fine_sys.solution, j2, H2) << std::endl;
        }

      for (unsigned int j = 0; j != fine_sys.n_vars(); ++j)
        {
          const std::string varname = fine_sys.variable_name(j);

          if (!coarse_sys.has_variable(varname))
            {
              libMesh::out << "Coarse solution system " << sysname
                           << " variable " << varname
                           << " not found in fine solution!" << std::endl;
              continue;
            }
        }

      if (!coarse_sys.n_vars() && !fine_sys.n_vars())
        libMesh::out << "No variables found in fine or coarse solution system "
                     << sysname << '!' << std::endl;
    }

  if (argc > 5)
    {
      libMesh::out << "Writing diff solution " << argv[5] << std::endl;

      for (std::size_t i = 0; i != sysnames.size(); ++i)
        {
          const std::string sysname = sysnames[i];
          const System & coarse_sys = coarse_es.get_system(sysname);
          const System & fine_sys = fine_es.get_system(sysname);

          std::unique_ptr<NumericVector<Number>> fine_solution = fine_sys.solution->clone();
          fine_sys.solution->zero();

          std::vector<unsigned int>
            var_remapping(fine_sys.n_vars(), libMesh::invalid_uint);

          for (unsigned int j = 0; j != fine_sys.n_vars(); ++j)
            {
              const std::string varname = fine_sys.variable_name(j);

              if (coarse_sys.has_variable(varname))
                var_remapping[j] = coarse_sys.variable_number(varname);
            }

          MeshFunction coarse_solution
            (coarse_es, *coarse_sys.solution,
             coarse_sys.get_dof_map(), var_remapping);
          coarse_solution.init();
          const DenseVector<Number> oom_value(fine_sys.n_vars(), 0);
          coarse_solution.enable_out_of_mesh_mode(oom_value);

          fine_sys.project_solution(&coarse_solution);
          *fine_sys.solution -= *fine_solution;
        }

      NameBasedIO(fine_mesh).write_equation_systems (argv[5], fine_es);
    }

  return 0;
}
