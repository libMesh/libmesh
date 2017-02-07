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



// <h1>Adaptivity Example 5 - Periodic Boundary Conditions with Adaptive Mesh Refinement</h1>
// \author Benjamin S. Kirk
// \date 2003
//
// This example uses the same simple, linear transient
// system as in example 10; however in this case periodic boundary
// conditions are applied at the sides of the domain.
//
// This code also contains an example use of ParsedFunction, to
// allow users to specify an exact solution on the command line.

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>

// Basic include file needed for the mesh functionality.
#include "libmesh/libmesh.h"
#include "libmesh/replicated_mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/gmv_io.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"

#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/parsed_function.h"

#include "libmesh/getpot.h"

// This example will solve a linear transient system,
// so we need to include the TransientLinearImplicitSystem definition.
#include "libmesh/transient_system.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/vector_value.h"

// To refine the mesh we need an ErrorEstimator
// object to figure out which elements to refine.
#include "libmesh/error_vector.h"
#include "libmesh/kelly_error_estimator.h"

// The definition of a geometric element
#include "libmesh/elem.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This function will assemble the system
// matrix and right-hand-side at each time step.  Note that
// since the system is linear we technically do not need to
// assmeble the matrix at each time step, but we will anyway.
// In subsequent examples we will employ adaptive mesh refinement,
// and with a changing mesh it will be necessary to rebuild the
// system matrix.
#ifdef LIBMESH_ENABLE_AMR
void assemble_cd (EquationSystems & es,
                  const std::string & system_name);
#endif

// Function prototype.  This function will initialize the system.
// Initialization functions are optional for systems.  They allow
// you to specify the initial values of the solution.  If an
// initialization function is not provided then the default (0)
// solution is provided.
void init_cd (EquationSystems & es,
              const std::string & system_name);

// Exact solution function prototype.  This gives the exact
// solution as a function of space and time.  In this case the
// initial condition will be taken as the exact solution at time 0,
// as will the Dirichlet boundary conditions at time t.
Real exact_solution (const Real x,
                     const Real y,
                     const Real t);

Number exact_value (const Point & p,
                    const Parameters & parameters,
                    const std::string &,
                    const std::string &)
{
  return exact_solution(p(0), p(1), parameters.get<Real> ("time"));
}

// With --enable-fparser, the user can also optionally set their own
// exact solution equations.
UniquePtr<FunctionBase<Number> > parsed_solution;


// Returns a string with 'number' formatted and placed directly
// into the string in some way
std::string exodus_filename(unsigned number);


// Begin the main program.  Note that the first
// statement in the program throws an error if
// you are in complex number mode, since this
// example is only intended to work with real
// numbers.
int main (int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

#if !defined(LIBMESH_ENABLE_AMR)
  libmesh_example_requires(false, "--enable-amr");
#elif !defined(LIBMESH_HAVE_XDR)
  // We use XDR support in our output here
  libmesh_example_requires(false, "--enable-xdr");
#elif !defined(LIBMESH_ENABLE_PERIODIC)
  libmesh_example_requires(false, "--enable-periodic");
#elif (LIBMESH_DOF_ID_BYTES == 8)
  libmesh_example_requires(false, "--with-dof-id-bytes=4");
#else

  // Our Trilinos interface does not yet support adaptive transient
  // problems
  libmesh_example_requires(libMesh::default_solver_package() != TRILINOS_SOLVERS, "--enable-petsc");

  // Brief message to the user regarding the program name
  // and command line arguments.

  // Use commandline parameter to specify if we are to
  // read in an initial solution or generate it ourself
  libMesh::out << "Usage:\n"
               <<"\t " << argv[0] << " -init_timestep 0\n"
               << "OR\n"
               <<"\t " << argv[0] << " -read_solution -init_timestep 26\n"
               << std::endl;

  libMesh::out << "Running: " << argv[0];

  for (int i=1; i<argc; i++)
    libMesh::out << " " << argv[i];

  libMesh::out << std::endl << std::endl;

  // Create a GetPot object to parse the command line
  GetPot command_line (argc, argv);


  // This boolean value is obtained from the command line, it is true
  // if the flag "-read_solution" is present, false otherwise.
  // It indicates whether we are going to read in
  // the mesh and solution files "saved_mesh.xda" and "saved_solution.xda"
  // or whether we are going to start from scratch by just reading
  // "mesh.xda"
  const bool read_solution = command_line.search("-read_solution");

  // This value is also obtained from the commandline and it specifies the
  // initial value for the t_step looping variable. We must
  // distinguish between the two cases here, whether we read in the
  // solution or we started from scratch, so that we do not overwrite the
  // gmv output files.
  unsigned int init_timestep = 0;

  // Search the command line for the "init_timestep" flag and if it is
  // present, set init_timestep accordingly.
  if (command_line.search("-init_timestep"))
    init_timestep = command_line.next(0);
  else
    {
      // This handy function will print the file name, line number,
      // specified message, and then throw an exception.
      libmesh_error_msg("ERROR: Initial timestep not specified!");
    }

  // This value is also obtained from the command line, and specifies
  // the number of time steps to take.
  unsigned int n_timesteps = 0;

  // Again do a search on the command line for the argument
  if (command_line.search("-n_timesteps"))
    n_timesteps = command_line.next(0);
  else
    libmesh_error_msg("ERROR: Number of timesteps not specified");

  // The user can specify a different exact solution on the command
  // line, if we have an expression parser compiled in
#ifdef LIBMESH_HAVE_FPARSER
  const bool have_expression = command_line.search("-exact_solution");
#else
  const bool have_expression = false;
#endif
  if (have_expression)
    parsed_solution.reset
      (new ParsedFunction<Number>(command_line.next(std::string())));

  // Skip this 2D example if libMesh was compiled as 1D-only.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Create a new mesh on the default MPI communicator.
  // DistributedMesh currently has a bug which is triggered by this
  // example.
  ReplicatedMesh mesh(init.comm());

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);
  MeshRefinement mesh_refinement (mesh);

  // Declare the system and its variables.
  // Begin by creating a transient system
  // named "Convection-Diffusion".
  TransientLinearImplicitSystem & system =
    equation_systems.add_system<TransientLinearImplicitSystem>
    ("Convection-Diffusion");

  // Give the system a pointer to the assembly function.
  system.attach_assemble_function (assemble_cd);

  // Creating and attaching Periodic Boundaries
  DofMap & dof_map = system.get_dof_map();

  // Create a boundary periodic with one displaced 2.0 in the x
  // direction
  PeriodicBoundary horz(RealVectorValue(2.0, 0., 0.));

  // Connect boundary ids 3 and 1 with it
  horz.myboundary = 3;
  horz.pairedboundary = 1;

  // Add it to the PeriodicBoundaries
  dof_map.add_periodic_boundary(horz);

  // Create a boundary periodic with one displaced 2.0 in the y
  // direction
  PeriodicBoundary vert(RealVectorValue(0., 2.0, 0.));

  // Connect boundary ids 0 and 2 with it
  vert.myboundary = 0;
  vert.pairedboundary = 2;

  // Add it to the PeriodicBoundaries
  dof_map.add_periodic_boundary(vert);

  // Next build or read the mesh.  We do this only *after* generating
  // periodic boundaries; otherwise a DistributedMesh won't know to
  // retain periodic neighbor elements.

  // First we process the case where we do not read in the solution
  if (!read_solution)
    {
      MeshTools::Generation::build_square(mesh, 2, 2, 0., 2., 0., 2., QUAD4);

      // Again do a search on the command line for an argument
      unsigned int n_refinements = 5;
      if (command_line.search("-n_refinements"))
        n_refinements = command_line.next(0);

      // Uniformly refine the mesh 5 times
      if (!read_solution)
        mesh_refinement.uniformly_refine (n_refinements);

      // Print information about the mesh to the screen.
      mesh.print_info();

      // Adds the variable "u" to "Convection-Diffusion".  "u"
      // will be approximated using first-order approximation.
      system.add_variable ("u", FIRST);

      // Give the system a pointer to the initialization function.
      system.attach_init_function (init_cd);
    }
  // Otherwise we read in the solution and mesh
  else
    {
      // Read in the mesh stored in "saved_mesh.xda"
      mesh.read("saved_mesh.xdr");

      // Print information about the mesh to the screen.
      mesh.print_info();

      // Read in the solution stored in "saved_solution.xda"
      equation_systems.read("saved_solution.xdr", DECODE);
    }

  // Initialize the data structures for the equation system.
  if (!read_solution)
    equation_systems.init ();
  else
    equation_systems.reinit ();

  // Print out the H1 norm of the initialized or saved solution, for
  // verification purposes:
  Real H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));

  libMesh::out << "Initial H1 norm = " << H1norm << std::endl << std::endl;

  // Prints information about the system to the screen.
  equation_systems.print_info();

  equation_systems.parameters.set<unsigned int>
    ("linear solver maximum iterations") = 250;
  equation_systems.parameters.set<Real>
    ("linear solver tolerance") = TOLERANCE;

  if (!read_solution)
    {
      // Write out the initial condition
#ifdef LIBMESH_HAVE_GMV
      GMVIO(mesh).write_equation_systems ("out.gmv.000",
                                          equation_systems);
#endif
#ifdef LIBMESH_HAVE_EXODUS_API
      ExodusII_IO(mesh).write_equation_systems (exodus_filename(0),
                                                equation_systems);
#endif
    }
  else
    {
      // Write out the solution that was read in
#ifdef LIBMESH_HAVE_GMV
      GMVIO(mesh).write_equation_systems ("solution_read_in.gmv",
                                          equation_systems);
#endif
#ifdef LIBMESH_HAVE_EXODUS_API
      ExodusII_IO(mesh).write_equation_systems ("solution_read_in.e",
                                                equation_systems);
#endif
    }


  // The Convection-Diffusion system requires that we specify
  // the flow velocity.  We will specify it as a RealVectorValue
  // data type and then use the Parameters object to pass it to
  // the assemble function.
  equation_systems.parameters.set<RealVectorValue>("velocity") =
    RealVectorValue (0.8, 0.8);

  // The Convection-Diffusion system also requires a specified
  // diffusivity.  We use an isotropic (hence Real) value.
  equation_systems.parameters.set<Real>("diffusivity") = 0.01;

  // Solve the system "Convection-Diffusion".  This will be done by
  // looping over the specified time interval and calling the
  // solve() member at each time step.  This will assemble the
  // system and call the linear solver.
  const Real dt = 0.025;
  system.time = init_timestep*dt;

  // Tell the MeshRefinement object about the periodic boundaries so
  // that it can get heuristics like level-one conformity and
  // unrefined island elimination right.
  mesh_refinement.set_periodic_boundaries_ptr(dof_map.get_periodic_boundaries());

  // We do 25 timesteps both before and after writing out the
  // intermediate solution
  for (unsigned int t_step=init_timestep; t_step<(init_timestep+n_timesteps); t_step++)
    {
      // Increment the time counter, set the time and the
      // time step size as parameters in the EquationSystem.
      system.time += dt;

      equation_systems.parameters.set<Real> ("time") = system.time;
      equation_systems.parameters.set<Real> ("dt")   = dt;

      // A pretty update message
      libMesh::out << " Solving time step ";

      {
        // Save flags to avoid polluting cout with custom precision values, etc.
        std::ios_base::fmtflags os_flags = libMesh::out.flags();

        libMesh::out << t_step
                     << ", time="
                     << std::setw(6)
                     << std::setprecision(3)
                     << std::setfill('0')
                     << std::left
                     << system.time
                     << "..."
                     << std::endl;

        // Restore flags
        libMesh::out.flags(os_flags);
      }

      // At this point we need to update the old
      // solution vector.  The old solution vector
      // will be the current solution vector from the
      // previous time step.  We will do this by extracting the
      // system from the EquationSystems object and using
      // vector assignment.  Since only TransientLinearImplicitSystems
      // (and systems derived from them) contain old solutions
      // we need to specify the system type when we ask for it.
      *system.old_local_solution = *system.current_local_solution;

      // The number of refinement steps per time step.
      unsigned int max_r_steps = 1;
      if (command_line.search("-max_r_steps"))
        max_r_steps = command_line.next(0);

      // A refinement loop.
      for (unsigned int r_step=0; r_step<max_r_steps+1; r_step++)
        {
          // Assemble & solve the linear system
          system.solve();

          // Print out the H1 norm, for verification purposes:
          H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));
          libMesh::out << "H1 norm = " << H1norm << std::endl;

          // Possibly refine the mesh
          if (r_step+1 <= max_r_steps)
            {
              libMesh::out << "  Refining the mesh..." << std::endl;

              // The ErrorVector is a particular StatisticsVector
              // for computing error information on a finite element mesh.
              ErrorVector error;

              // The ErrorEstimator class interrogates a finite element
              // solution and assigns to each element a positive error value.
              // This value is used for deciding which elements to refine
              // and which to coarsen.
              //ErrorEstimator* error_estimator = new KellyErrorEstimator;
              KellyErrorEstimator error_estimator;

              // Compute the error for each active element using the provided
              // flux_jump indicator.  Note in general you will need to
              // provide an error estimator specifically designed for your
              // application.
              error_estimator.estimate_error (system,
                                              error);

              // This takes the error in error and decides which elements
              // will be coarsened or refined.  Any element within 20% of the
              // maximum error on any element will be refined, and any
              // element within 7% of the minimum error on any element might
              // be coarsened. Note that the elements flagged for refinement
              // will be refined, but those flagged for coarsening _might_ be
              // coarsened.
              mesh_refinement.refine_fraction() = 0.80;
              mesh_refinement.coarsen_fraction() = 0.07;
              mesh_refinement.max_h_level() = 5;
              mesh_refinement.flag_elements_by_error_fraction (error);

              // This call actually refines and coarsens the flagged
              // elements.
              mesh_refinement.refine_and_coarsen_elements();

              // This call reinitializes the EquationSystems object for
              // the newly refined mesh.  One of the steps in the
              // reinitialization is projecting the solution,
              // old_solution, etc... vectors from the old mesh to
              // the current one.
              equation_systems.reinit ();
            }
        }

      // Again do a search on the command line for an argument
      unsigned int output_freq = 10;
      if (command_line.search("-output_freq"))
        output_freq = command_line.next(0);

      // Output every 10 timesteps to file.
      if ((t_step+1)%output_freq == 0)
        {
#ifdef LIBMESH_HAVE_GMV
          // std::ostringstream file_name;
          // out << "out.gmv."
          //     << std::setw(3)
          //     << std::setfill('0')
          //     << std::right
          //     << t_step+1;
          // GMVIO(mesh).write_equation_systems (file_name.str(),
          //                                     equation_systems);
#endif
#ifdef LIBMESH_HAVE_EXODUS_API
          // So... if paraview is told to open a file called out.e.{N}, it automatically tries to
          // open out.e.{N-1}, out.e.{N-2}, etc.  If we name the file something else, we can work
          // around that issue, but the right thing to do (for adaptive meshes) is to write a filename
          // with the adaptation step into a separate file.
          ExodusII_IO(mesh).write_equation_systems (exodus_filename(t_step+1),
                                                    equation_systems);
#endif
        }
    }

  if (!read_solution)
    {
      // Print out the H1 norm of the saved solution, for verification purposes:
      H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));

      libMesh::out << "Final H1 norm = " << H1norm << std::endl << std::endl;

      mesh.write("saved_mesh.xdr");
      equation_systems.write("saved_solution.xdr", ENCODE);
#ifdef LIBMESH_HAVE_GMV
      GMVIO(mesh).write_equation_systems ("saved_solution.gmv",
                                          equation_systems);
#endif
#ifdef LIBMESH_HAVE_EXODUS_API
      ExodusII_IO(mesh).write_equation_systems ("saved_solution.e",
                                                equation_systems);
#endif
    }
#endif // #ifndef LIBMESH_ENABLE_AMR

  return 0;
}

// Here we define the initialization routine for the
// Convection-Diffusion system.  This routine is
// responsible for applying the initial conditions to
// the system.
void init_cd (EquationSystems & es,
              const std::string & libmesh_dbg_var(system_name))
{
  // It is a good idea to make sure we are initializing
  // the proper system.
  libmesh_assert_equal_to (system_name, "Convection-Diffusion");

  // Get a reference to the Convection-Diffusion system object.
  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem>("Convection-Diffusion");

  // Project initial conditions at time 0
  es.parameters.set<Real> ("time") = system.time = 0;

  if (parsed_solution.get())
    system.project_solution(parsed_solution.get(), libmesh_nullptr);
  else
    system.project_solution(exact_value, libmesh_nullptr, es.parameters);
}



// This function defines the assembly routine which
// will be called at each time step.  It is responsible
// for computing the proper matrix entries for the
// element stiffness matrices and right-hand sides.
#ifdef LIBMESH_ENABLE_AMR
void assemble_cd (EquationSystems & es,
                  const std::string & libmesh_dbg_var(system_name))
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Convection-Diffusion");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Convection-Diffusion system object.
  TransientLinearImplicitSystem & system =
    es.get_system<TransientLinearImplicitSystem> ("Convection-Diffusion");

  // Get the Finite Element type for the first (and only)
  // variable in the system.
  FEType fe_type = system.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> fe      (FEBase::build(dim, fe_type));
  UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));

  // A Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  QGauss qrule (dim,   fe_type.default_quadrature_order());
  QGauss qface (dim-1, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule      (&qrule);
  fe_face->attach_quadrature_rule (&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.  We will start
  // with the element Jacobian * quadrature weight at each integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature
  // points.
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the DofMap
  // in future examples.
  const DofMap & dof_map = system.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Here we extract the velocity & parameters that we put in the
  // EquationSystems object.
  const RealVectorValue velocity =
    es.parameters.get<RealVectorValue> ("velocity");

  const Real diffusivity =
    es.parameters.get<Real> ("diffusivity");

  const Real dt = es.parameters.get<Real>   ("dt");

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the active_elem_iterator.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem * elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      fe->reinit (elem);

      // Zero the element matrix and right-hand side before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      Ke.resize (dof_indices.size(),
                 dof_indices.size());

      Fe.resize (dof_indices.size());

      // Now we will build the element matrix and right-hand-side.
      // Constructing the RHS requires the solution and its
      // gradient from the previous timestep.  This myst be
      // calculated at each quadrature point by summing the
      // solution degree-of-freedom values by the appropriate
      // weight functions.
      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Values to hold the old solution & its gradient.
          Number u_old = 0.;
          Gradient grad_u_old;

          // Compute the old solution & its gradient.
          for (std::size_t l=0; l<phi.size(); l++)
            {
              u_old += phi[l][qp]*system.old_solution  (dof_indices[l]);

              // This will work,
              // grad_u_old += dphi[l][qp]*system.old_solution (dof_indices[l]);
              // but we can do it without creating a temporary like this:
              grad_u_old.add_scaled (dphi[l][qp], system.old_solution (dof_indices[l]));
            }

          // Now compute the element matrix and RHS contributions.
          for (std::size_t i=0; i<phi.size(); i++)
            {
              // The RHS contribution
              Fe(i) += JxW[qp]*(
                                // Mass matrix term
                                u_old*phi[i][qp] +
                                -.5*dt*(
                                        // Convection term
                                        // (grad_u_old may be complex, so the
                                        // order here is important!)
                                        (grad_u_old*velocity)*phi[i][qp] +

                                        // Diffusion term
                                        diffusivity*(grad_u_old*dphi[i][qp]))
                                );

              for (std::size_t j=0; j<phi.size(); j++)
                {
                  // The matrix contribution
                  Ke(i,j) += JxW[qp]*(
                                      // Mass-matrix
                                      phi[i][qp]*phi[j][qp] +
                                      .5*dt*(
                                             // Convection term
                                             (velocity*dphi[j][qp])*phi[i][qp] +
                                             // Diffusion term
                                             diffusivity*(dphi[i][qp]*dphi[j][qp]))
                                      );
                }
            }
        }

      // We have now built the element matrix and RHS vector in terms
      // of the element degrees of freedom.  However, it is possible
      // that some of the element DOFs are constrained to enforce
      // solution continuity, i.e. they are not really "free".  We need
      // to constrain those DOFs in terms of non-constrained DOFs to
      // ensure a continuous solution.  The
      // DofMap::constrain_element_matrix_and_vector() method does
      // just that.
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The SparseMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.
      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);

    }
  // Finished computing the sytem matrix and right-hand side.
}
#endif // #ifdef LIBMESH_ENABLE_AMR




std::string exodus_filename(unsigned number)
{
  std::ostringstream oss;

  oss << "out_";
  oss << std::setw(3) << std::setfill('0') << number;
  oss << ".e";

  return oss.str();
}
