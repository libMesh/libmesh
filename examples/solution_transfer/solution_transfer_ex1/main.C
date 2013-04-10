// Basic include files needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/explicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_DTK
#include "libmesh/dtk_solution_transfer.h"
#endif

// Bring in everything from the libMesh namespace
using namespace libMesh;


Number initial_value(const Point& p,
                     const Parameters& parameters,
                     const std::string&,
                     const std::string&)
{
  return p(0)*p(0) + 1; // x^2 + 1
}

void initialize(EquationSystems& es,
                const std::string& system_name)
{
  ExplicitSystem & system = es.get_system<ExplicitSystem>(system_name);
  es.parameters.set<Real> ("time") = system.time = 0;
  system.project_solution(initial_value, NULL, es.parameters);
}

int main(int argc, char* argv[])
{
  LibMeshInit init (argc, argv);

#ifdef LIBMESH_HAVE_DTK

  Mesh from_mesh(init.comm());
  MeshTools::Generation::build_cube(from_mesh, 4, 4, 4, 0, 1, 0, 1, 0, 1, HEX8);
  from_mesh.print_info();
  EquationSystems from_es(from_mesh);
  System & from_sys = from_es.add_system<ExplicitSystem>("From");
  unsigned int from_var = from_sys.add_variable("from");
  from_sys.attach_init_function(initialize);
  from_es.init();

  ExodusII_IO(from_mesh).write_equation_systems("from.e", from_es);

  Mesh to_mesh;
  MeshTools::Generation::build_cube(to_mesh, 5, 5, 5, 0, 1, 0, 1, 0, 1, TET4);
  to_mesh.print_info();
  EquationSystems to_es(to_mesh);
  System & to_sys = to_es.add_system<ExplicitSystem>("To");
  unsigned int to_var = to_sys.add_variable("to");
  to_es.init();

  DTKSolutionTransfer dtk_transfer;

  dtk_transfer.transfer(from_sys.variable(from_var), to_sys.variable(to_var));

  to_es.update();
  ExodusII_IO(to_mesh).write_equation_systems("to.e", to_es);

#endif

  return 0;
}
