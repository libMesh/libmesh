// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/getpot.h"

// Application includes
#include "advection_system.h"

// C++ include files that we need
#include <iostream>

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function to set the initial condition
void set_initial_condition (EquationSystems& es,
                            const std::string& system_name);
Number initial_condition (const Point& p,
                          const Parameters& parameters,
                          const std::string&,
                          const std::string&);

// The main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // Create GetPot object to parse the command line
  GetPot command_line (argc, argv);

  // Allow the user specify a custom input file on the command line
  // using the "-i" flag, and provide a sensible default value in case
  // they don't.
  std::string parameters_filename = "advection_2D.in";
  if ( command_line.search(2, "-i", "--input_file") )
    parameters_filename = command_line.next(parameters_filename);

  // Read mesh options from the input file
  GetPot infile(parameters_filename);
  unsigned int n_elem = infile("n_elem", 1);       // Determines the number of elements in the "truth" mesh
  const unsigned int dim = 2;                      // The number of spatial dimensions

  // Build a mesh
  Mesh mesh (init.comm());
  MeshTools::Generation::build_square (mesh,
                                       n_elem, n_elem,
                                       0., 1.,
                                       0., 1.,
                                       QUAD9);

  // Create an equation systems object
  EquationSystems equation_systems (mesh);

  // Create an AdvectionSystem within the EquationSystems
  AdvectionSystem & advection_system =
    equation_systems.add_system<AdvectionSystem> ("Advection2D");

  // Process the input file values
  advection_system.process_parameters_file(parameters_filename);

  // Set the initial condition function
  advection_system.attach_init_function (set_initial_condition);

  equation_systems.init ();

  equation_systems.print_info();
  mesh.print_info();

  advection_system.print_info();

  advection_system.assemble_all_matrices();
  // advection_system.write_out_discretization_matrices();

#ifdef LIBMESH_HAVE_EXODUS_API
  ExodusII_IO(mesh).write_equation_systems ("claw_solution.0000.e", equation_systems);
#endif

  advection_system.solve_conservation_law();

  return 0;
}

void set_initial_condition (EquationSystems& es,
                            const std::string&)
{
  LinearImplicitSystem & system =
    es.get_system<LinearImplicitSystem>("Advection2D");

  system.project_solution(initial_condition, /*gradient=*/nullptr, es.parameters);
}

Number initial_condition (const Point& p,
                          const Parameters& ,
                          const std::string&,
                          const std::string&)
{
  Real x = p(0);
  Real y = p(1);

  return std::exp( -100.*(std::pow(x-0.5,2.) + std::pow(y-0.5,2.)) );
}

