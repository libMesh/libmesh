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


// Open the solution files named on standard input, take the average
// of their fields' values, and output to a new solution file.

#include "libmesh/libmesh.h"

#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/numeric_vector.h"

#include "libmesh/exodusII_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/gmv_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/medit_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/vtk_io.h"

using namespace libMesh;

unsigned char dim = 2; // This gets overridden by most mesh formats

int main(int argc, char ** argv)
{
  LibMeshInit init(argc, argv);

  Mesh mesh1(init.comm(), dim);
  EquationSystems es1(mesh1);

  libMesh::out << "Usage: " << argv[0]
               << " outputsolution mesh firstsolution [secondsolution ...]" << std::endl;

  mesh1.read(argv[2]);
  libMesh::out << "Loaded mesh " << argv[2] << std::endl;

  mesh1.print_info();

  es1.read(argv[3],
           EquationSystems::READ_HEADER |
           EquationSystems::READ_DATA |
           EquationSystems::READ_BASIC_ONLY);
  libMesh::out << "Loaded first solution " << argv[3] << std::endl;

  es1.print_info();

  std::vector<std::string> sysnames;
  std::vector<NumericVector<libMesh::Number> *> summed_solutions;

  for (unsigned int s = 0; s != es1.n_systems(); ++s)
    {
      sysnames.push_back(es1.get_system(s).name());
      summed_solutions.push_back(es1.get_system(s).solution->clone().release());
    }

  for (int i=4; i < argc; ++i)
    {
      Mesh mesh2(init.comm(), dim);
      EquationSystems es2(mesh2);
      mesh2.read(argv[2]);
      es2.read(argv[i],
               EquationSystems::READ_HEADER |
               EquationSystems::READ_DATA |
               EquationSystems::READ_BASIC_ONLY);
      libMesh::out << "Loaded next solution " << argv[i] << std::endl;

      for (std::size_t s = 0; s != sysnames.size(); ++s)
        {
          if (!es2.has_system(sysnames[s]))
            libmesh_error_msg("EquationSystems object does not have " << sysnames[s]);

          (*summed_solutions[s]) += *es2.get_system(s).solution;
        }
    }

  int n_solutions = argc - 3;

  for (std::size_t s = 0; s != sysnames.size(); ++s)
    {
      (*summed_solutions[s]) /= n_solutions;
      es1.get_system(s).solution->swap(*summed_solutions[s]);
      es1.get_system(s).solution->close();

      delete summed_solutions[s];
    }

  std::string outputname(argv[1]);

  if (outputname.find(".xda") != std::string::npos ||
      outputname.find(".xdr") != std::string::npos)
    es1.write(outputname.c_str(),
              EquationSystems::WRITE_DATA |
              EquationSystems::WRITE_ADDITIONAL_DATA);
  else if (outputname.find(".gnuplot") != std::string::npos)
    GnuPlotIO(mesh1).write_equation_systems (outputname, es1);
  else if (outputname.find(".gmv") != std::string::npos)
    GMVIO(mesh1).write_equation_systems (outputname, es1);
  else if (outputname.find(".mesh") != std::string::npos)
    MEDITIO(mesh1).write_equation_systems (outputname, es1);
  else if (outputname.find(".msh") != std::string::npos)
    GmshIO(mesh1).write_equation_systems (outputname, es1);
  else if (outputname.find(".plt") != std::string::npos)
    TecplotIO(mesh1).write_equation_systems (outputname, es1);
  else if (outputname.find(".vtu") != std::string::npos)
    VTKIO(mesh1).write_equation_systems (outputname, es1);
  else if (outputname.find(".e") != std::string::npos)
    ExodusII_IO(mesh1).write_equation_systems (outputname, es1);
  else if (outputname.find(".n") != std::string::npos)
    Nemesis_IO(mesh1).write_equation_systems (outputname, es1);
  libMesh::out << "Wrote solution " << outputname << std::endl;

  return 0;
}
