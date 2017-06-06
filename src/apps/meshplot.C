// Open the mesh and solution file named in command line arguments,
// and plot them out to the given filename.

#include "libmesh/libmesh.h"

#include "libmesh/equation_systems.h"
#include "libmesh/mesh.h"
#include "libmesh/namebased_io.h"

int main(int argc, char ** argv)
{
  using namespace libMesh;

  LibMeshInit init(argc, argv);

  Mesh mesh(init.comm());
  EquationSystems es(mesh);

  if (argc < 3 || argc > 4)
    libmesh_error_msg
      ("Usage: " << argv[0] <<
       " inputmesh [inputsolution] outputplot");


  START_LOG("mesh.read()", "main");
  mesh.read(argv[1]);
  STOP_LOG("mesh.read()", "main");
  libMesh::out << "Loaded mesh " << argv[1] << std::endl;
  mesh.print_info();

  if (argc > 3)
    {
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
    }

  START_LOG("write_equation_systems()", "main");
  std::string outputname(argv[argc-1]);
  if ((outputname.find(".xda") != std::string::npos) ||
      (outputname.find(".xdr") != std::string::npos))
    {
      mesh.write("mesh-"+outputname);
      es.write("soln-"+outputname);
    }
  else
    NameBasedIO(mesh).write_equation_systems (outputname, es);
  STOP_LOG("write_equation_systems()", "main");

  libMesh::out << "Wrote output " << outputname << std::endl;
}
