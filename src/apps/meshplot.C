// Open the mesh and solution file named in command line arguments,
// and plot them out to the given filename.

#include "libmesh/libmesh.h"

#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"

#include "libmesh/exodusII_io.h"
#include "libmesh/gmsh_io.h"
#include "libmesh/gmv_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/medit_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/tecplot_io.h"
#include "libmesh/vtk_io.h"

int main(int argc, char** argv)
{
  using namespace libMesh;

  LibMeshInit init(argc, argv);

  Mesh mesh(init.comm());
  EquationSystems es(mesh);

  libMesh::out << "Usage: " << argv[0]
            << " inputmesh inputsolution outputplot" << std::endl;

  START_LOG("mesh.read()", "main");
  mesh.read(argv[1]);
  STOP_LOG("mesh.read()", "main");
  libMesh::out << "Loaded mesh " << argv[1] << std::endl;
  mesh.print_info();

  START_LOG("es.read()", "main");
  std::string solnname = argv[2];

  es.read(solnname,
          EquationSystems::READ_HEADER |
          EquationSystems::READ_DATA |
          EquationSystems::READ_ADDITIONAL_DATA |
          EquationSystems::READ_BASIC_ONLY);
  STOP_LOG("es.read()", "main");
  libMesh::out << "Loaded solution " << solnname << std::endl;
  es.print_info();

  START_LOG("write_equation_systems()", "main");
  std::string outputname(argv[3]);

  if (outputname.find(".gnuplot") != std::string::npos)
    GnuPlotIO(mesh).write_equation_systems (outputname, es);
  else if (outputname.find(".gmv") != std::string::npos)
    GMVIO(mesh).write_equation_systems (outputname, es);
  else if (outputname.find(".mesh") != std::string::npos)
    MEDITIO(mesh).write_equation_systems (outputname, es);
  else if (outputname.find(".msh") != std::string::npos)
    GmshIO(mesh).write_equation_systems (outputname, es);
  else if (outputname.find(".plt") != std::string::npos)
    TecplotIO(mesh).write_equation_systems (outputname, es);
  else if (outputname.find(".vtu") != std::string::npos)
    VTKIO(mesh).write_equation_systems (outputname, es);
  else if (outputname.find(".e") != std::string::npos)
    ExodusII_IO(mesh).write_equation_systems (outputname, es);
  else if (outputname.find(".n") != std::string::npos)
    Nemesis_IO(mesh).write_equation_systems (outputname, es);

  STOP_LOG("write_equation_systems()", "main");
  libMesh::out << "Wrote output " << outputname << std::endl;
}
