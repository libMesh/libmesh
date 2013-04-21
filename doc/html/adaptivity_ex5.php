<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("adaptivity_ex5",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file adaptivity_ex5.C with comments: </h1> 
<div class = "comment">
<h1>Adaptivity Example 5 - Periodic Boundary Conditions with Adaptive Mesh Refinement</h1>

<br><br>This example uses the same simple, linear transient
system as in example 10; however in this case periodic boundary
conditions are applied at the sides of the domain.

<br><br>This code also contains an example use of ParsedFunction, to
allow users to specify an exact solution on the command line.


<br><br>C++ include files that we need
</div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
        #include &lt;algorithm&gt;
        #include &lt;cstdlib&gt; // *must* precede &lt;cmath&gt; for proper std:abs() on PGI, Sun Studio CC
        #include &lt;cmath&gt;
        
</pre>
</div>
<div class = "comment">
Basic include file needed for the mesh functionality.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh.h"
        #include "libmesh/serial_mesh.h"
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
        
</pre>
</div>
<div class = "comment">
This example will solve a linear transient system,
so we need to include the \p TransientLinearImplicitSystem definition.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/transient_system.h"
        #include "libmesh/linear_implicit_system.h"
        #include "libmesh/vector_value.h"
        
</pre>
</div>
<div class = "comment">
To refine the mesh we need an \p ErrorEstimator
object to figure out which elements to refine.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/error_vector.h"
        #include "libmesh/kelly_error_estimator.h"
        
</pre>
</div>
<div class = "comment">
The definition of a geometric element
</div>

<div class ="fragment">
<pre>
        #include "libmesh/elem.h"
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
</pre>
</div>
<div class = "comment">
Function prototype.  This function will assemble the system
matrix and right-hand-side at each time step.  Note that
since the system is linear we technically do not need to
assmeble the matrix at each time step, but we will anyway.
In subsequent examples we will employ adaptive mesh refinement,
and with a changing mesh it will be necessary to rebuild the
system matrix.
</div>

<div class ="fragment">
<pre>
        void assemble_cd (EquationSystems& es,
                          const std::string& system_name);
        
</pre>
</div>
<div class = "comment">
Function prototype.  This function will initialize the system.
Initialization functions are optional for systems.  They allow
you to specify the initial values of the solution.  If an
initialization function is not provided then the default (0)
solution is provided.
</div>

<div class ="fragment">
<pre>
        void init_cd (EquationSystems& es,
                      const std::string& system_name);
        
</pre>
</div>
<div class = "comment">
Exact solution function prototype.  This gives the exact
solution as a function of space and time.  In this case the
initial condition will be taken as the exact solution at time 0,
as will the Dirichlet boundary conditions at time t.
</div>

<div class ="fragment">
<pre>
        Real exact_solution (const Real x,
                             const Real y,
                             const Real t);
        
        Number exact_value (const Point& p,
                            const Parameters& parameters,
                            const std::string&,
                            const std::string&)
        {
          return exact_solution(p(0), p(1), parameters.get&lt;Real&gt; ("time"));
        }
        
</pre>
</div>
<div class = "comment">
With --enable-fparser, the user can also optionally set their own
exact solution equations.
</div>

<div class ="fragment">
<pre>
        FunctionBase&lt;Number&gt;* parsed_solution = NULL;
        
        
</pre>
</div>
<div class = "comment">
Returns a string with 'number' formatted and placed directly
into the string in some way
</div>

<div class ="fragment">
<pre>
        std::string exodus_filename(unsigned number);
        
        
</pre>
</div>
<div class = "comment">
Begin the main program.  Note that the first
statement in the program throws an error if
you are in complex number mode, since this
example is only intended to work with real
numbers.
</div>

<div class ="fragment">
<pre>
        int main (int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Initialize libMesh.
</div>

<div class ="fragment">
<pre>
          LibMeshInit init (argc, argv);
        
        #if !defined(LIBMESH_ENABLE_AMR)
          libmesh_example_assert(false, "--enable-amr");
        #elif !defined(LIBMESH_HAVE_XDR)
</pre>
</div>
<div class = "comment">
We use XDR support in our output here
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(false, "--enable-xdr");
        #elif !defined(LIBMESH_ENABLE_PERIODIC)
          libmesh_example_assert(false, "--enable-periodic");
        #else
        
</pre>
</div>
<div class = "comment">
Our Trilinos interface does not yet support adaptive transient
problems
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(libMesh::default_solver_package() != TRILINOS_SOLVERS, "--enable-petsc");
        
</pre>
</div>
<div class = "comment">
Brief message to the user regarding the program name
and command line arguments.


<br><br>Use commandline parameter to specify if we are to
read in an initial solution or generate it ourself
</div>

<div class ="fragment">
<pre>
          std::cout &lt;&lt; "Usage:\n"
            &lt;&lt;"\t " &lt;&lt; argv[0] &lt;&lt; " -init_timestep 0\n"
            &lt;&lt; "OR\n"
            &lt;&lt;"\t " &lt;&lt; argv[0] &lt;&lt; " -read_solution -init_timestep 26\n"
            &lt;&lt; std::endl;
        
          std::cout &lt;&lt; "Running: " &lt;&lt; argv[0];
        
          for (int i=1; i&lt;argc; i++)
            std::cout &lt;&lt; " " &lt;&lt; argv[i];
        
          std::cout &lt;&lt; std::endl &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Create a GetPot object to parse the command line
</div>

<div class ="fragment">
<pre>
          GetPot command_line (argc, argv);
        
        
</pre>
</div>
<div class = "comment">
This boolean value is obtained from the command line, it is true
if the flag "-read_solution" is present, false otherwise.
It indicates whether we are going to read in
the mesh and solution files "saved_mesh.xda" and "saved_solution.xda"
or whether we are going to start from scratch by just reading
"mesh.xda"
</div>

<div class ="fragment">
<pre>
          const bool read_solution   = command_line.search("-read_solution");
        
</pre>
</div>
<div class = "comment">
This value is also obtained from the commandline and it specifies the
initial value for the t_step looping variable. We must
distinguish between the two cases here, whether we read in the
solution or we started from scratch, so that we do not overwrite the
gmv output files.
</div>

<div class ="fragment">
<pre>
          unsigned int init_timestep = 0;
        
</pre>
</div>
<div class = "comment">
Search the command line for the "init_timestep" flag and if it is
present, set init_timestep accordingly.
</div>

<div class ="fragment">
<pre>
          if(command_line.search("-init_timestep"))
            init_timestep = command_line.next(0);
          else
            {
              if (libMesh::processor_id() == 0)
                std::cerr &lt;&lt; "ERROR: Initial timestep not specified\n" &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
This handy function will print the file name, line number,
and then abort.  Currrently the library does not use C++
exception handling.
</div>

<div class ="fragment">
<pre>
              libmesh_error();
            }
        
</pre>
</div>
<div class = "comment">
This value is also obtained from the command line, and specifies
the number of time steps to take.
</div>

<div class ="fragment">
<pre>
          unsigned int n_timesteps = 0;
        
</pre>
</div>
<div class = "comment">
Again do a search on the command line for the argument
</div>

<div class ="fragment">
<pre>
          if(command_line.search("-n_timesteps"))
            n_timesteps = command_line.next(0);
          else
            {
              std::cout &lt;&lt; "ERROR: Number of timesteps not specified\n" &lt;&lt; std::endl;
              libmesh_error();
            }
        
</pre>
</div>
<div class = "comment">
The user can specify a different exact solution on the command
line, if we have an expression parser compiled in
</div>

<div class ="fragment">
<pre>
        #ifdef LIBMESH_HAVE_FPARSER
          const bool have_expression = command_line.search("-exact_solution");
        #else
          const bool have_expression = false;
        #endif
          if (have_expression)
            parsed_solution = new ParsedFunction&lt;Number&gt;(command_line.next(std::string()));
        
</pre>
</div>
<div class = "comment">
Skip this 2D example if libMesh was compiled as 1D-only.
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(2 &lt;= LIBMESH_DIM, "2D support");
        
</pre>
</div>
<div class = "comment">
Create a new mesh on the default MPI communicator.
ParallelMesh doesn't yet understand periodic BCs, plus
we still need some work on automatic parallel restarts
</div>

<div class ="fragment">
<pre>
          SerialMesh mesh(init.comm());
        
</pre>
</div>
<div class = "comment">
Create an equation systems object.
</div>

<div class ="fragment">
<pre>
          EquationSystems equation_systems (mesh);
          MeshRefinement mesh_refinement (mesh);
        
</pre>
</div>
<div class = "comment">
First we process the case where we do not read in the solution
</div>

<div class ="fragment">
<pre>
          if(!read_solution)
            {
              MeshTools::Generation::build_square(mesh, 2, 2, 0., 2., 0., 2., QUAD4);
        
</pre>
</div>
<div class = "comment">
Again do a search on the command line for an argument
</div>

<div class ="fragment">
<pre>
              unsigned int n_refinements = 5;
              if(command_line.search("-n_refinements"))
                n_refinements = command_line.next(0);
        
</pre>
</div>
<div class = "comment">
Uniformly refine the mesh 5 times
</div>

<div class ="fragment">
<pre>
              if(!read_solution)
                mesh_refinement.uniformly_refine (n_refinements);
        
</pre>
</div>
<div class = "comment">
Print information about the mesh to the screen.
</div>

<div class ="fragment">
<pre>
              mesh.print_info();
        
        
</pre>
</div>
<div class = "comment">
Declare the system and its variables.
Begin by creating a transient system
named "Convection-Diffusion".
</div>

<div class ="fragment">
<pre>
              TransientLinearImplicitSystem & system =
                equation_systems.add_system&lt;TransientLinearImplicitSystem&gt;
                ("Convection-Diffusion");
        
</pre>
</div>
<div class = "comment">
Adds the variable "u" to "Convection-Diffusion".  "u"
will be approximated using first-order approximation.
</div>

<div class ="fragment">
<pre>
              system.add_variable ("u", FIRST);
        
</pre>
</div>
<div class = "comment">
Give the system a pointer to the initialization function.
</div>

<div class ="fragment">
<pre>
              system.attach_init_function (init_cd);
            }
</pre>
</div>
<div class = "comment">
Otherwise we read in the solution and mesh
</div>

<div class ="fragment">
<pre>
          else
            {
</pre>
</div>
<div class = "comment">
Read in the mesh stored in "saved_mesh.xda"
</div>

<div class ="fragment">
<pre>
              mesh.read("saved_mesh.xdr");
        
</pre>
</div>
<div class = "comment">
Print information about the mesh to the screen.
</div>

<div class ="fragment">
<pre>
              mesh.print_info();
        
</pre>
</div>
<div class = "comment">
Read in the solution stored in "saved_solution.xda"
</div>

<div class ="fragment">
<pre>
              equation_systems.read("saved_solution.xdr", libMeshEnums::DECODE);
            }
        
</pre>
</div>
<div class = "comment">
Get a reference to the system so that we can attach things to it
</div>

<div class ="fragment">
<pre>
          TransientLinearImplicitSystem & system =
            equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;
            ("Convection-Diffusion");
        
</pre>
</div>
<div class = "comment">
Give the system a pointer to the assembly function.
</div>

<div class ="fragment">
<pre>
          system.attach_assemble_function (assemble_cd);
        
</pre>
</div>
<div class = "comment">
Creating and attaching Periodic Boundaries
</div>

<div class ="fragment">
<pre>
          DofMap & dof_map = system.get_dof_map();
        
</pre>
</div>
<div class = "comment">
Create a boundary periodic with one displaced 2.0 in the x
direction
</div>

<div class ="fragment">
<pre>
          PeriodicBoundary horz(RealVectorValue(2.0, 0., 0.));
        
</pre>
</div>
<div class = "comment">
Connect boundary ids 3 and 1 with it
</div>

<div class ="fragment">
<pre>
          horz.myboundary = 3;
          horz.pairedboundary = 1;
        
</pre>
</div>
<div class = "comment">
Add it to the PeriodicBoundaries
</div>

<div class ="fragment">
<pre>
          dof_map.add_periodic_boundary(horz);
        
        
</pre>
</div>
<div class = "comment">
Create a boundary periodic with one displaced 2.0 in the y
direction
</div>

<div class ="fragment">
<pre>
          PeriodicBoundary vert(RealVectorValue(0., 2.0, 0.));
        
</pre>
</div>
<div class = "comment">
Connect boundary ids 0 and 2 with it
</div>

<div class ="fragment">
<pre>
          vert.myboundary = 0;
          vert.pairedboundary = 2;
        
</pre>
</div>
<div class = "comment">
Add it to the PeriodicBoundaries
</div>

<div class ="fragment">
<pre>
          dof_map.add_periodic_boundary(vert);
        
</pre>
</div>
<div class = "comment">
Initialize the data structures for the equation system.
</div>

<div class ="fragment">
<pre>
          if(!read_solution)
            equation_systems.init ();
          else
            equation_systems.reinit ();
        
</pre>
</div>
<div class = "comment">
Print out the H1 norm of the initialized or saved solution, for
verification purposes:
</div>

<div class ="fragment">
<pre>
          Real H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));
        
          std::cout &lt;&lt; "Initial H1 norm = " &lt;&lt; H1norm &lt;&lt; std::endl &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Prints information about the system to the screen.
</div>

<div class ="fragment">
<pre>
          equation_systems.print_info();
        
          equation_systems.parameters.set&lt;unsigned int&gt;
            ("linear solver maximum iterations") = 250;
          equation_systems.parameters.set&lt;Real&gt;
            ("linear solver tolerance") = TOLERANCE;
        
          if(!read_solution)
          {
</pre>
</div>
<div class = "comment">
Write out the initial condition
</div>

<div class ="fragment">
<pre>
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
</pre>
</div>
<div class = "comment">
Write out the solution that was read in
</div>

<div class ="fragment">
<pre>
        #ifdef LIBMESH_HAVE_GMV
            GMVIO(mesh).write_equation_systems ("solution_read_in.gmv",
                                                equation_systems);
        #endif
        #ifdef LIBMESH_HAVE_EXODUS_API
            ExodusII_IO(mesh).write_equation_systems ("solution_read_in.e",
                                                equation_systems);
        #endif
          }
        
        
</pre>
</div>
<div class = "comment">
The Convection-Diffusion system requires that we specify
the flow velocity.  We will specify it as a RealVectorValue
data type and then use the Parameters object to pass it to
the assemble function.
</div>

<div class ="fragment">
<pre>
          equation_systems.parameters.set&lt;RealVectorValue&gt;("velocity") =
            RealVectorValue (0.8, 0.8);
        
</pre>
</div>
<div class = "comment">
The Convection-Diffusion system also requires a specified
diffusivity.  We use an isotropic (hence Real) value.
</div>

<div class ="fragment">
<pre>
          equation_systems.parameters.set&lt;Real&gt;("diffusivity") = 0.01;
        
</pre>
</div>
<div class = "comment">
Solve the system "Convection-Diffusion".  This will be done by
looping over the specified time interval and calling the
\p solve() member at each time step.  This will assemble the
system and call the linear solver.


<br><br></div>

<div class ="fragment">
<pre>
          const Real dt = 0.025;
          system.time   = init_timestep*dt;
        
        
</pre>
</div>
<div class = "comment">
Tell the MeshRefinement object about the periodic boundaries so
that it can get heuristics like level-one conformity and
unrefined island elimination right.
</div>

<div class ="fragment">
<pre>
          mesh_refinement.set_periodic_boundaries_ptr(dof_map.get_periodic_boundaries());
        
</pre>
</div>
<div class = "comment">
We do 25 timesteps both before and after writing out the
intermediate solution
</div>

<div class ="fragment">
<pre>
          for(unsigned int t_step=init_timestep;
                           t_step&lt;(init_timestep+n_timesteps);
                           t_step++)
            {
</pre>
</div>
<div class = "comment">
Increment the time counter, set the time and the
time step size as parameters in the EquationSystem.
</div>

<div class ="fragment">
<pre>
              system.time += dt;
        
              equation_systems.parameters.set&lt;Real&gt; ("time") = system.time;
              equation_systems.parameters.set&lt;Real&gt; ("dt")   = dt;
        
</pre>
</div>
<div class = "comment">
A pretty update message
</div>

<div class ="fragment">
<pre>
              std::cout &lt;&lt; " Solving time step ";
        
              {
</pre>
</div>
<div class = "comment">
Save flags to avoid polluting cout with custom precision values, etc.
</div>

<div class ="fragment">
<pre>
                std::ios_base::fmtflags os_flags = std::cout.flags();
        
                std::cout &lt;&lt; t_step
                          &lt;&lt; ", time="
                          &lt;&lt; std::setw(6)
                          &lt;&lt; std::setprecision(3)
                          &lt;&lt; std::setfill('0')
                          &lt;&lt; std::left
                          &lt;&lt; system.time
                          &lt;&lt; "..."
                          &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Restore flags
</div>

<div class ="fragment">
<pre>
                std::cout.flags(os_flags);
              }
        
</pre>
</div>
<div class = "comment">
At this point we need to update the old
solution vector.  The old solution vector
will be the current solution vector from the
previous time step.  We will do this by extracting the
system from the \p EquationSystems object and using
vector assignment.  Since only \p TransientLinearImplicitSystems
(and systems derived from them) contain old solutions
we need to specify the system type when we ask for it.
</div>

<div class ="fragment">
<pre>
              *system.old_local_solution = *system.current_local_solution;
        
</pre>
</div>
<div class = "comment">
The number of refinement steps per time step.
</div>

<div class ="fragment">
<pre>
              unsigned int max_r_steps = 1;
              if(command_line.search("-max_r_steps"))
                max_r_steps = command_line.next(0);
        
</pre>
</div>
<div class = "comment">
A refinement loop.
</div>

<div class ="fragment">
<pre>
              for (unsigned int r_step=0; r_step&lt;max_r_steps+1; r_step++)
                {
</pre>
</div>
<div class = "comment">
Assemble & solve the linear system
</div>

<div class ="fragment">
<pre>
                  system.solve();
        
</pre>
</div>
<div class = "comment">
Print out the H1 norm, for verification purposes:
</div>

<div class ="fragment">
<pre>
                  H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));
                  std::cout &lt;&lt; "H1 norm = " &lt;&lt; H1norm &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Possibly refine the mesh
</div>

<div class ="fragment">
<pre>
                  if (r_step+1 &lt;= max_r_steps)
                    {
                      std::cout &lt;&lt; "  Refining the mesh..." &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
The \p ErrorVector is a particular \p StatisticsVector
for computing error information on a finite element mesh.
</div>

<div class ="fragment">
<pre>
                      ErrorVector error;
        
</pre>
</div>
<div class = "comment">
The \p ErrorEstimator class interrogates a finite element
solution and assigns to each element a positive error value.
This value is used for deciding which elements to refine
and which to coarsen.
ErrorEstimator* error_estimator = new KellyErrorEstimator;
</div>

<div class ="fragment">
<pre>
                      KellyErrorEstimator error_estimator;
        
</pre>
</div>
<div class = "comment">
Compute the error for each active element using the provided
\p flux_jump indicator.  Note in general you will need to
provide an error estimator specifically designed for your
application.
</div>

<div class ="fragment">
<pre>
                      error_estimator.estimate_error (system,
                                                      error);
        
</pre>
</div>
<div class = "comment">
This takes the error in \p error and decides which elements
will be coarsened or refined.  Any element within 20% of the
maximum error on any element will be refined, and any
element within 7% of the minimum error on any element might
be coarsened. Note that the elements flagged for refinement
will be refined, but those flagged for coarsening _might_ be
coarsened.
</div>

<div class ="fragment">
<pre>
                      mesh_refinement.refine_fraction() = 0.80;
                      mesh_refinement.coarsen_fraction() = 0.07;
                      mesh_refinement.max_h_level() = 5;
                      mesh_refinement.flag_elements_by_error_fraction (error);
        
</pre>
</div>
<div class = "comment">
This call actually refines and coarsens the flagged
elements.
</div>

<div class ="fragment">
<pre>
                      mesh_refinement.refine_and_coarsen_elements();
        
</pre>
</div>
<div class = "comment">
This call reinitializes the \p EquationSystems object for
the newly refined mesh.  One of the steps in the
reinitialization is projecting the \p solution,
\p old_solution, etc... vectors from the old mesh to
the current one.
</div>

<div class ="fragment">
<pre>
                      equation_systems.reinit ();
                    }
                }
        
</pre>
</div>
<div class = "comment">
Again do a search on the command line for an argument
</div>

<div class ="fragment">
<pre>
              unsigned int output_freq = 10;
              if(command_line.search("-output_freq"))
                output_freq = command_line.next(0);
        
</pre>
</div>
<div class = "comment">
Output every 10 timesteps to file.
</div>

<div class ="fragment">
<pre>
              if ( (t_step+1)%output_freq == 0)
                {
</pre>
</div>
<div class = "comment">
OStringStream file_name;


<br><br></div>

<div class ="fragment">
<pre>
        #ifdef LIBMESH_HAVE_GMV
</pre>
</div>
<div class = "comment">
file_name << "out.gmv.";
OSSRealzeroright(file_name,3,0,t_step+1);

<br><br>GMVIO(mesh).write_equation_systems (file_name.str(),
equation_systems);
</div>

<div class ="fragment">
<pre>
        #endif
        #ifdef LIBMESH_HAVE_EXODUS_API
</pre>
</div>
<div class = "comment">
So... if paraview is told to open a file called out.e.{N}, it automatically tries to
open out.e.{N-1}, out.e.{N-2}, etc.  If we name the file something else, we can work
around that issue, but the right thing to do (for adaptive meshes) is to write a filename
with the adaptation step into a separate file.
</div>

<div class ="fragment">
<pre>
                  ExodusII_IO(mesh).write_equation_systems (exodus_filename(t_step+1),
                                                            equation_systems);
        #endif
                }
            }
        
          if(!read_solution)
            {
</pre>
</div>
<div class = "comment">
Print out the H1 norm of the saved solution, for verification purposes:
</div>

<div class ="fragment">
<pre>
              H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));
        
              std::cout &lt;&lt; "Final H1 norm = " &lt;&lt; H1norm &lt;&lt; std::endl &lt;&lt; std::endl;
        
              mesh.write("saved_mesh.xdr");
              equation_systems.write("saved_solution.xdr", libMeshEnums::ENCODE);
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
        
</pre>
</div>
<div class = "comment">
We might have a parser to clean up
</div>

<div class ="fragment">
<pre>
          delete parsed_solution;
        
          return 0;
        }
        
</pre>
</div>
<div class = "comment">
Here we define the initialization routine for the
Convection-Diffusion system.  This routine is
responsible for applying the initial conditions to
the system.
</div>

<div class ="fragment">
<pre>
        void init_cd (EquationSystems& es,
                      const std::string& system_name)
        {
</pre>
</div>
<div class = "comment">
It is a good idea to make sure we are initializing
the proper system.
</div>

<div class ="fragment">
<pre>
          libmesh_assert_equal_to (system_name, "Convection-Diffusion");
        
</pre>
</div>
<div class = "comment">
Get a reference to the Convection-Diffusion system object.
</div>

<div class ="fragment">
<pre>
          TransientLinearImplicitSystem & system =
            es.get_system&lt;TransientLinearImplicitSystem&gt;("Convection-Diffusion");
        
</pre>
</div>
<div class = "comment">
Project initial conditions at time 0
</div>

<div class ="fragment">
<pre>
          es.parameters.set&lt;Real&gt; ("time") = system.time = 0;
        
          if (parsed_solution)
            system.project_solution(parsed_solution, NULL);
          else
            system.project_solution(exact_value, NULL, es.parameters);
        }
        
        
        
</pre>
</div>
<div class = "comment">
This function defines the assembly routine which
will be called at each time step.  It is responsible
for computing the proper matrix entries for the
element stiffness matrices and right-hand sides.
</div>

<div class ="fragment">
<pre>
        void assemble_cd (EquationSystems& es,
                          const std::string& system_name)
        {
        #ifdef LIBMESH_ENABLE_AMR
</pre>
</div>
<div class = "comment">
It is a good idea to make sure we are assembling
the proper system.
</div>

<div class ="fragment">
<pre>
          libmesh_assert_equal_to (system_name, "Convection-Diffusion");
        
</pre>
</div>
<div class = "comment">
Get a constant reference to the mesh object.
</div>

<div class ="fragment">
<pre>
          const MeshBase& mesh = es.get_mesh();
        
</pre>
</div>
<div class = "comment">
The dimension that we are running
</div>

<div class ="fragment">
<pre>
          const unsigned int dim = mesh.mesh_dimension();
        
</pre>
</div>
<div class = "comment">
Get a reference to the Convection-Diffusion system object.
</div>

<div class ="fragment">
<pre>
          TransientLinearImplicitSystem & system =
            es.get_system&lt;TransientLinearImplicitSystem&gt; ("Convection-Diffusion");
        
</pre>
</div>
<div class = "comment">
Get the Finite Element type for the first (and only)
variable in the system.
</div>

<div class ="fragment">
<pre>
          FEType fe_type = system.variable_type(0);
        
</pre>
</div>
<div class = "comment">
Build a Finite Element object of the specified type.  Since the
\p FEBase::build() member dynamically creates memory we will
store the object as an \p AutoPtr<FEBase>.  This can be thought
of as a pointer that will clean up after itself.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;FEBase&gt; fe      (FEBase::build(dim, fe_type));
          AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
        
</pre>
</div>
<div class = "comment">
A Gauss quadrature rule for numerical integration.
Let the \p FEType object decide what order rule is appropriate.
</div>

<div class ="fragment">
<pre>
          QGauss qrule (dim,   fe_type.default_quadrature_order());
          QGauss qface (dim-1, fe_type.default_quadrature_order());
        
</pre>
</div>
<div class = "comment">
Tell the finite element object to use our quadrature rule.
</div>

<div class ="fragment">
<pre>
          fe-&gt;attach_quadrature_rule      (&qrule);
          fe_face-&gt;attach_quadrature_rule (&qface);
        
</pre>
</div>
<div class = "comment">
Here we define some references to cell-specific data that
will be used to assemble the linear system.  We will start
with the element Jacobian * quadrature weight at each integration point.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt;& JxW      = fe-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The element shape functions evaluated at the quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi = fe-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The element shape function gradients evaluated at the quadrature
points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi = fe-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
A reference to the \p DofMap object for this system.  The \p DofMap
object handles the index translation from node and element numbers
to degree of freedom numbers.  We will talk more about the \p DofMap
in future examples.
</div>

<div class ="fragment">
<pre>
          const DofMap& dof_map = system.get_dof_map();
        
</pre>
</div>
<div class = "comment">
Define data structures to contain the element matrix
and right-hand-side vector contribution.  Following
basic finite element terminology we will denote these
"Ke" and "Fe".
</div>

<div class ="fragment">
<pre>
          DenseMatrix&lt;Number&gt; Ke;
          DenseVector&lt;Number&gt; Fe;
        
</pre>
</div>
<div class = "comment">
This vector will hold the degree of freedom indices for
the element.  These define where in the global system
the element degrees of freedom get mapped.
</div>

<div class ="fragment">
<pre>
          std::vector&lt;dof_id_type&gt; dof_indices;
        
</pre>
</div>
<div class = "comment">
Here we extract the velocity & parameters that we put in the
EquationSystems object.
</div>

<div class ="fragment">
<pre>
          const RealVectorValue velocity =
            es.parameters.get&lt;RealVectorValue&gt; ("velocity");
        
          const Real diffusivity =
            es.parameters.get&lt;Real&gt; ("diffusivity");
        
          const Real dt = es.parameters.get&lt;Real&gt;   ("dt");
        
</pre>
</div>
<div class = "comment">
Now we will loop over all the elements in the mesh that
live on the local processor. We will compute the element
matrix and right-hand-side contribution.  Since the mesh
will be refined we want to only consider the ACTIVE elements,
hence we use a variant of the \p active_elem_iterator.
</div>

<div class ="fragment">
<pre>
          MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
          const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
        
          for ( ; el != end_el; ++el)
            {
</pre>
</div>
<div class = "comment">
Store a pointer to the element we are currently
working on.  This allows for nicer syntax later.
</div>

<div class ="fragment">
<pre>
              const Elem* elem = *el;
        
</pre>
</div>
<div class = "comment">
Get the degree of freedom indices for the
current element.  These define where in the global
matrix and right-hand-side this element will
contribute to.
</div>

<div class ="fragment">
<pre>
              dof_map.dof_indices (elem, dof_indices);
        
</pre>
</div>
<div class = "comment">
Compute the element-specific data for the current
element.  This involves computing the location of the
quadrature points (q_point) and the shape functions
(phi, dphi) for the current element.
</div>

<div class ="fragment">
<pre>
              fe-&gt;reinit (elem);
        
</pre>
</div>
<div class = "comment">
Zero the element matrix and right-hand side before
summing them.  We use the resize member here because
the number of degrees of freedom might have changed from
the last element.  Note that this will be the case if the
element type is different (i.e. the last element was a
triangle, now we are on a quadrilateral).
</div>

<div class ="fragment">
<pre>
              Ke.resize (dof_indices.size(),
                         dof_indices.size());
        
              Fe.resize (dof_indices.size());
        
</pre>
</div>
<div class = "comment">
Now we will build the element matrix and right-hand-side.
Constructing the RHS requires the solution and its
gradient from the previous timestep.  This myst be
calculated at each quadrature point by summing the
solution degree-of-freedom values by the appropriate
weight functions.
</div>

<div class ="fragment">
<pre>
              for (unsigned int qp=0; qp&lt;qrule.n_points(); qp++)
                {
</pre>
</div>
<div class = "comment">
Values to hold the old solution & its gradient.
</div>

<div class ="fragment">
<pre>
                  Number   u_old = 0.;
                  Gradient grad_u_old;
        
</pre>
</div>
<div class = "comment">
Compute the old solution & its gradient.
</div>

<div class ="fragment">
<pre>
                  for (unsigned int l=0; l&lt;phi.size(); l++)
                    {
                      u_old      += phi[l][qp]*system.old_solution  (dof_indices[l]);
        
</pre>
</div>
<div class = "comment">
This will work,
grad_u_old += dphi[l][qp]*system.old_solution (dof_indices[l]);
but we can do it without creating a temporary like this:
</div>

<div class ="fragment">
<pre>
                      grad_u_old.add_scaled (dphi[l][qp],system.old_solution (dof_indices[l]));
                    }
        
</pre>
</div>
<div class = "comment">
Now compute the element matrix and RHS contributions.
</div>

<div class ="fragment">
<pre>
                  for (unsigned int i=0; i&lt;phi.size(); i++)
                    {
</pre>
</div>
<div class = "comment">
The RHS contribution
</div>

<div class ="fragment">
<pre>
                      Fe(i) += JxW[qp]*(
</pre>
</div>
<div class = "comment">
Mass matrix term
</div>

<div class ="fragment">
<pre>
                                        u_old*phi[i][qp] +
                                        -.5*dt*(
</pre>
</div>
<div class = "comment">
Convection term
(grad_u_old may be complex, so the
order here is important!)
</div>

<div class ="fragment">
<pre>
                                                (grad_u_old*velocity)*phi[i][qp] +
        
</pre>
</div>
<div class = "comment">
Diffusion term
</div>

<div class ="fragment">
<pre>
                                                diffusivity*(grad_u_old*dphi[i][qp]))
                                        );
        
                      for (unsigned int j=0; j&lt;phi.size(); j++)
                        {
</pre>
</div>
<div class = "comment">
The matrix contribution
</div>

<div class ="fragment">
<pre>
                          Ke(i,j) += JxW[qp]*(
</pre>
</div>
<div class = "comment">
Mass-matrix
</div>

<div class ="fragment">
<pre>
                                              phi[i][qp]*phi[j][qp] +
                                              .5*dt*(
</pre>
</div>
<div class = "comment">
Convection term
</div>

<div class ="fragment">
<pre>
                                                     (velocity*dphi[j][qp])*phi[i][qp] +
</pre>
</div>
<div class = "comment">
Diffusion term
</div>

<div class ="fragment">
<pre>
                                                     diffusivity*(dphi[i][qp]*dphi[j][qp]))
                                              );
                        }
                    }
                }
        
</pre>
</div>
<div class = "comment">
We have now built the element matrix and RHS vector in terms
of the element degrees of freedom.  However, it is possible
that some of the element DOFs are constrained to enforce
solution continuity, i.e. they are not really "free".  We need
to constrain those DOFs in terms of non-constrained DOFs to
ensure a continuous solution.  The
\p DofMap::constrain_element_matrix_and_vector() method does
just that.
</div>

<div class ="fragment">
<pre>
              dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        
</pre>
</div>
<div class = "comment">
The element matrix and right-hand-side are now built
for this element.  Add them to the global matrix and
right-hand-side vector.  The \p SparseMatrix::add_matrix()
and \p NumericVector::add_vector() members do this for us.
</div>

<div class ="fragment">
<pre>
              system.matrix-&gt;add_matrix (Ke, dof_indices);
              system.rhs-&gt;add_vector    (Fe, dof_indices);
        
            }
</pre>
</div>
<div class = "comment">
Finished computing the sytem matrix and right-hand side.
</div>

<div class ="fragment">
<pre>
        #endif // #ifdef LIBMESH_ENABLE_AMR
        }
        
        
        
        
        std::string exodus_filename(unsigned number)
        {
          std::ostringstream oss;
        
          oss &lt;&lt; "out_";
          oss &lt;&lt; std::setw(3) &lt;&lt; std::setfill('0') &lt;&lt; number;
          oss &lt;&lt; ".e";
        
          return oss.str();
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file adaptivity_ex5.C without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;cstdlib&gt; <I><FONT COLOR="#B22222">// *must* precede &lt;cmath&gt; for proper std:abs() on PGI, Sun Studio CC
</FONT></I>  #include &lt;cmath&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/serial_mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/gmv_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/periodic_boundaries.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/periodic_boundary.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/parsed_function.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/transient_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/vector_value.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/kelly_error_estimator.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_cd (EquationSystems&amp; es,
                    <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  <B><FONT COLOR="#228B22">void</FONT></B> init_cd (EquationSystems&amp; es,
                <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  Real exact_solution (<B><FONT COLOR="#228B22">const</FONT></B> Real x,
                       <B><FONT COLOR="#228B22">const</FONT></B> Real y,
                       <B><FONT COLOR="#228B22">const</FONT></B> Real t);
  
  Number exact_value (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                      <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp; parameters,
                      <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
                      <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;)
  {
    <B><FONT COLOR="#A020F0">return</FONT></B> exact_solution(p(0), p(1), parameters.get&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>));
  }
  
  FunctionBase&lt;Number&gt;* parsed_solution = NULL;
  
  
  <B><FONT COLOR="#5F9EA0">std</FONT></B>::string exodus_filename(<B><FONT COLOR="#228B22">unsigned</FONT></B> number);
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_ENABLE_AMR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-amr&quot;</FONT></B>);
  #elif !defined(LIBMESH_HAVE_XDR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-xdr&quot;</FONT></B>);
  #elif !defined(LIBMESH_ENABLE_PERIODIC)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-periodic&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    libmesh_example_assert(libMesh::default_solver_package() != TRILINOS_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Usage:\n&quot;</FONT></B>
      &lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;\t &quot;</FONT></B> &lt;&lt; argv[0] &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; -init_timestep 0\n&quot;</FONT></B>
      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;OR\n&quot;</FONT></B>
      &lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;\t &quot;</FONT></B> &lt;&lt; argv[0] &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; -read_solution -init_timestep 26\n&quot;</FONT></B>
      &lt;&lt; std::endl;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Running: &quot;</FONT></B> &lt;&lt; argv[0];
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">int</FONT></B> i=1; i&lt;argc; i++)
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B> &lt;&lt; argv[i];
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl &lt;&lt; std::endl;
  
    GetPot command_line (argc, argv);
  
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> read_solution   = command_line.search(<B><FONT COLOR="#BC8F8F">&quot;-read_solution&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> init_timestep = 0;
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(command_line.search(<B><FONT COLOR="#BC8F8F">&quot;-init_timestep&quot;</FONT></B>))
      init_timestep = command_line.next(0);
    <B><FONT COLOR="#A020F0">else</FONT></B>
      {
        <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR: Initial timestep not specified\n&quot;</FONT></B> &lt;&lt; std::endl;
  
        libmesh_error();
      }
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_timesteps = 0;
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(command_line.search(<B><FONT COLOR="#BC8F8F">&quot;-n_timesteps&quot;</FONT></B>))
      n_timesteps = command_line.next(0);
    <B><FONT COLOR="#A020F0">else</FONT></B>
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR: Number of timesteps not specified\n&quot;</FONT></B> &lt;&lt; std::endl;
        libmesh_error();
      }
  
  #ifdef LIBMESH_HAVE_FPARSER
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> have_expression = command_line.search(<B><FONT COLOR="#BC8F8F">&quot;-exact_solution&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> have_expression = false;
  #endif
    <B><FONT COLOR="#A020F0">if</FONT></B> (have_expression)
      parsed_solution = <B><FONT COLOR="#A020F0">new</FONT></B> ParsedFunction&lt;Number&gt;(command_line.next(std::string()));
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
  
    SerialMesh mesh(init.comm());
  
    EquationSystems equation_systems (mesh);
    MeshRefinement mesh_refinement (mesh);
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(!read_solution)
      {
        <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square(mesh, 2, 2, 0., 2., 0., 2., QUAD4);
  
        <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_refinements = 5;
        <B><FONT COLOR="#A020F0">if</FONT></B>(command_line.search(<B><FONT COLOR="#BC8F8F">&quot;-n_refinements&quot;</FONT></B>))
          n_refinements = command_line.next(0);
  
        <B><FONT COLOR="#A020F0">if</FONT></B>(!read_solution)
          mesh_refinement.uniformly_refine (n_refinements);
  
        mesh.print_info();
  
  
        TransientLinearImplicitSystem &amp; system =
          equation_systems.add_system&lt;TransientLinearImplicitSystem&gt;
          (<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
        system.add_variable (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, FIRST);
  
        system.attach_init_function (init_cd);
      }
    <B><FONT COLOR="#A020F0">else</FONT></B>
      {
        mesh.read(<B><FONT COLOR="#BC8F8F">&quot;saved_mesh.xdr&quot;</FONT></B>);
  
        mesh.print_info();
  
        equation_systems.read(<B><FONT COLOR="#BC8F8F">&quot;saved_solution.xdr&quot;</FONT></B>, libMeshEnums::DECODE);
      }
  
    TransientLinearImplicitSystem &amp; system =
      equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;
      (<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    system.attach_assemble_function (assemble_cd);
  
    DofMap &amp; dof_map = system.get_dof_map();
  
    PeriodicBoundary horz(RealVectorValue(2.0, 0., 0.));
  
    horz.myboundary = 3;
    horz.pairedboundary = 1;
  
    dof_map.add_periodic_boundary(horz);
  
  
    PeriodicBoundary vert(RealVectorValue(0., 2.0, 0.));
  
    vert.myboundary = 0;
    vert.pairedboundary = 2;
  
    dof_map.add_periodic_boundary(vert);
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(!read_solution)
      equation_systems.init ();
    <B><FONT COLOR="#A020F0">else</FONT></B>
      equation_systems.reinit ();
  
    Real H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Initial H1 norm = &quot;</FONT></B> &lt;&lt; H1norm &lt;&lt; std::endl &lt;&lt; std::endl;
  
    equation_systems.print_info();
  
    equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;
      (<B><FONT COLOR="#BC8F8F">&quot;linear solver maximum iterations&quot;</FONT></B>) = 250;
    equation_systems.parameters.set&lt;Real&gt;
      (<B><FONT COLOR="#BC8F8F">&quot;linear solver tolerance&quot;</FONT></B>) = TOLERANCE;
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(!read_solution)
    {
  #ifdef LIBMESH_HAVE_GMV
      GMVIO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;out.gmv.000&quot;</FONT></B>,
                                          equation_systems);
  #endif
  #ifdef LIBMESH_HAVE_EXODUS_API
      ExodusII_IO(mesh).write_equation_systems (exodus_filename(0),
                                                equation_systems);
  #endif
    }
    <B><FONT COLOR="#A020F0">else</FONT></B>
    {
  #ifdef LIBMESH_HAVE_GMV
      GMVIO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;solution_read_in.gmv&quot;</FONT></B>,
                                          equation_systems);
  #endif
  #ifdef LIBMESH_HAVE_EXODUS_API
      ExodusII_IO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;solution_read_in.e&quot;</FONT></B>,
                                          equation_systems);
  #endif
    }
  
  
    equation_systems.parameters.set&lt;RealVectorValue&gt;(<B><FONT COLOR="#BC8F8F">&quot;velocity&quot;</FONT></B>) =
      RealVectorValue (0.8, 0.8);
  
    equation_systems.parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;diffusivity&quot;</FONT></B>) = 0.01;
  
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real dt = 0.025;
    system.time   = init_timestep*dt;
  
  
    mesh_refinement.set_periodic_boundaries_ptr(dof_map.get_periodic_boundaries());
  
    <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> t_step=init_timestep;
                     t_step&lt;(init_timestep+n_timesteps);
                     t_step++)
      {
        system.time += dt;
  
        equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>) = system.time;
        equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;dt&quot;</FONT></B>)   = dt;
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Solving time step &quot;</FONT></B>;
  
        {
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::ios_base::fmtflags os_flags = std::cout.flags();
  
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; t_step
                    &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, time=&quot;</FONT></B>
                    &lt;&lt; std::setw(6)
                    &lt;&lt; std::setprecision(3)
                    &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
                    &lt;&lt; std::left
                    &lt;&lt; system.time
                    &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;...&quot;</FONT></B>
                    &lt;&lt; std::endl;
  
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout.flags(os_flags);
        }
  
        *system.old_local_solution = *system.current_local_solution;
  
        <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_r_steps = 1;
        <B><FONT COLOR="#A020F0">if</FONT></B>(command_line.search(<B><FONT COLOR="#BC8F8F">&quot;-max_r_steps&quot;</FONT></B>))
          max_r_steps = command_line.next(0);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> r_step=0; r_step&lt;max_r_steps+1; r_step++)
          {
            system.solve();
  
            H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;H1 norm = &quot;</FONT></B> &lt;&lt; H1norm &lt;&lt; std::endl;
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (r_step+1 &lt;= max_r_steps)
              {
                <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;  Refining the mesh...&quot;</FONT></B> &lt;&lt; std::endl;
  
                ErrorVector error;
  
                KellyErrorEstimator error_estimator;
  
                error_estimator.estimate_error (system,
                                                error);
  
                mesh_refinement.refine_fraction() = 0.80;
                mesh_refinement.coarsen_fraction() = 0.07;
                mesh_refinement.max_h_level() = 5;
                mesh_refinement.flag_elements_by_error_fraction (error);
  
                mesh_refinement.refine_and_coarsen_elements();
  
                equation_systems.reinit ();
              }
          }
  
        <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> output_freq = 10;
        <B><FONT COLOR="#A020F0">if</FONT></B>(command_line.search(<B><FONT COLOR="#BC8F8F">&quot;-output_freq&quot;</FONT></B>))
          output_freq = command_line.next(0);
  
        <B><FONT COLOR="#A020F0">if</FONT></B> ( (t_step+1)%output_freq == 0)
          {
  
  #ifdef LIBMESH_HAVE_GMV
  #endif
  #ifdef LIBMESH_HAVE_EXODUS_API
            ExodusII_IO(mesh).write_equation_systems (exodus_filename(t_step+1),
                                                      equation_systems);
  #endif
          }
      }
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(!read_solution)
      {
        H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Final H1 norm = &quot;</FONT></B> &lt;&lt; H1norm &lt;&lt; std::endl &lt;&lt; std::endl;
  
        mesh.write(<B><FONT COLOR="#BC8F8F">&quot;saved_mesh.xdr&quot;</FONT></B>);
        equation_systems.write(<B><FONT COLOR="#BC8F8F">&quot;saved_solution.xdr&quot;</FONT></B>, libMeshEnums::ENCODE);
  #ifdef LIBMESH_HAVE_GMV
        GMVIO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;saved_solution.gmv&quot;</FONT></B>,
                                            equation_systems);
  #endif
  #ifdef LIBMESH_HAVE_EXODUS_API
        ExodusII_IO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;saved_solution.e&quot;</FONT></B>,
                                            equation_systems);
  #endif
      }
  #endif <I><FONT COLOR="#B22222">// #ifndef LIBMESH_ENABLE_AMR
</FONT></I>  
    <B><FONT COLOR="#A020F0">delete</FONT></B> parsed_solution;
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> init_cd (EquationSystems&amp; es,
                <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    TransientLinearImplicitSystem &amp; system =
      es.get_system&lt;TransientLinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    es.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>) = system.time = 0;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (parsed_solution)
      system.project_solution(parsed_solution, NULL);
    <B><FONT COLOR="#A020F0">else</FONT></B>
      system.project_solution(exact_value, NULL, es.parameters);
  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_cd (EquationSystems&amp; es,
                    <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
  #ifdef LIBMESH_ENABLE_AMR
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    TransientLinearImplicitSystem &amp; system =
      es.get_system&lt;TransientLinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    FEType fe_type = system.variable_type(0);
  
    AutoPtr&lt;FEBase&gt; fe      (FEBase::build(dim, fe_type));
    AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
  
    QGauss qrule (dim,   fe_type.default_quadrature_order());
    QGauss qface (dim-1, fe_type.default_quadrature_order());
  
    fe-&gt;attach_quadrature_rule      (&amp;qrule);
    fe_face-&gt;attach_quadrature_rule (&amp;qface);
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW      = fe-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
  
    DenseMatrix&lt;Number&gt; Ke;
    DenseVector&lt;Number&gt; Fe;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
  
    <B><FONT COLOR="#228B22">const</FONT></B> RealVectorValue velocity =
      es.parameters.get&lt;RealVectorValue&gt; (<B><FONT COLOR="#BC8F8F">&quot;velocity&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real diffusivity =
      es.parameters.get&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;diffusivity&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real dt = es.parameters.get&lt;Real&gt;   (<B><FONT COLOR="#BC8F8F">&quot;dt&quot;</FONT></B>);
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
  
        fe-&gt;reinit (elem);
  
        Ke.resize (dof_indices.size(),
                   dof_indices.size());
  
        Fe.resize (dof_indices.size());
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
          {
            Number   u_old = 0.;
            Gradient grad_u_old;
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> l=0; l&lt;phi.size(); l++)
              {
                u_old      += phi[l][qp]*system.old_solution  (dof_indices[l]);
  
                grad_u_old.add_scaled (dphi[l][qp],system.old_solution (dof_indices[l]));
              }
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
              {
                Fe(i) += JxW[qp]*(
                                  u_old*phi[i][qp] +
                                  -.5*dt*(
                                          (grad_u_old*velocity)*phi[i][qp] +
  
                                          diffusivity*(grad_u_old*dphi[i][qp]))
                                  );
  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
                  {
                    Ke(i,j) += JxW[qp]*(
                                        phi[i][qp]*phi[j][qp] +
                                        .5*dt*(
                                               (velocity*dphi[j][qp])*phi[i][qp] +
                                               diffusivity*(dphi[i][qp]*dphi[j][qp]))
                                        );
                  }
              }
          }
  
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
  
        system.matrix-&gt;add_matrix (Ke, dof_indices);
        system.rhs-&gt;add_vector    (Fe, dof_indices);
  
      }
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_ENABLE_AMR
</FONT></I>  }
  
  
  
  
  <B><FONT COLOR="#5F9EA0">std</FONT></B>::string exodus_filename(<B><FONT COLOR="#228B22">unsigned</FONT></B> number)
  {
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream oss;
  
    oss &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out_&quot;</FONT></B>;
    oss &lt;&lt; std::setw(3) &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>) &lt;&lt; number;
    oss &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.e&quot;</FONT></B>;
  
    <B><FONT COLOR="#A020F0">return</FONT></B> oss.str();
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex5'
***************************************************************
* Running Example adaptivity_ex5:
*  mpirun -np 4 example-devel -n_timesteps 25 -n_refinements 5 -output_freq 10 -init_timestep 0 -exact_solution '10*exp(-(pow(x-0.8*t-0.2,2)+pow(y-0.8*t-0.2,2))/(0.01*(4*t+1)))/(4*t+1)' -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
Usage:
	 /net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex5/.libs/lt-example-devel -init_timestep 0
OR
	 /net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex5/.libs/lt-example-devel -read_solution -init_timestep 26

Running: /net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex5/.libs/lt-example-devel -n_timesteps 25 -n_refinements 5 -output_freq 10 -init_timestep 0 -exact_solution 10*exp(-(pow(x-0.8*t-0.2,2) + pow(y-0.8*t-0.2,2))/(0.01*(4*t+1)))/(4*t+1) -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=4225
    n_local_nodes()=1099
  n_elem()=5460
    n_local_elem()=1398
    n_active_elem()=4096
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
  n_threads()=1
  processor_id()=0

Initial H1 norm = 17.4175

 EquationSystems
  n_systems()=1
   System #0, "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=4225
    n_local_dofs()=1099
    n_constrained_dofs()=129
    n_local_constrained_dofs()=0
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.80947
      Average Off-Processor Bandwidth <= 0.573728
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 14
    DofMap Constraints
      Number of DoF Constraints = 129
      Average DoF Constraint Length= 1
      Number of Node Constraints = 129
      Maximum Node Constraint Length= 2
      Average Node Constraint Length= 2

 Solving time step 0, time=0.0250...
H1 norm = 15.9
  Refining the mesh...
H1 norm = 15.9
 Solving time step 1, time=0.0500...
H1 norm = 14.6
  Refining the mesh...
H1 norm = 14.6
 Solving time step 2, time=0.0750...
H1 norm = 13.5
  Refining the mesh...
H1 norm = 13.5
 Solving time step 3, time=0.1000...
H1 norm = 12.6
  Refining the mesh...
H1 norm = 12.6
 Solving time step 4, time=0.1250...
H1 norm = 11.7
  Refining the mesh...
H1 norm = 11.7
 Solving time step 5, time=0.1500...
H1 norm = 11
  Refining the mesh...
H1 norm = 11
 Solving time step 6, time=0.1750...
H1 norm = 10.4
  Refining the mesh...
H1 norm = 10.4
 Solving time step 7, time=0.2000...
H1 norm = 9.82
  Refining the mesh...
H1 norm = 9.82
 Solving time step 8, time=0.2250...
H1 norm = 9.31
  Refining the mesh...
H1 norm = 9.31
 Solving time step 9, time=0.2500...
H1 norm = 8.85
  Refining the mesh...
H1 norm = 8.85
 Solving time step 10, time=0.2750...
H1 norm = 8.43
  Refining the mesh...
H1 norm = 8.43
 Solving time step 11, time=0.3000...
H1 norm = 8.05
  Refining the mesh...
H1 norm = 8.05
 Solving time step 12, time=0.3250...
H1 norm = 7.71
  Refining the mesh...
H1 norm = 7.71
 Solving time step 13, time=0.3500...
H1 norm = 7.39
  Refining the mesh...
H1 norm = 7.39
 Solving time step 14, time=0.3750...
H1 norm = 7.1
  Refining the mesh...
H1 norm = 7.1
 Solving time step 15, time=0.4000...
H1 norm = 6.83
  Refining the mesh...
H1 norm = 6.83
 Solving time step 16, time=0.4250...
H1 norm = 6.58
  Refining the mesh...
H1 norm = 6.58
 Solving time step 17, time=0.4500...
H1 norm = 6.35
  Refining the mesh...
H1 norm = 6.35
 Solving time step 18, time=0.4750...
H1 norm = 6.13
  Refining the mesh...
H1 norm = 6.13
 Solving time step 19, time=0.5000...
H1 norm = 5.93
  Refining the mesh...
H1 norm = 5.93
 Solving time step 20, time=0.5250...
H1 norm = 5.74
  Refining the mesh...
H1 norm = 5.74
 Solving time step 21, time=0.5500...
H1 norm = 5.56
  Refining the mesh...
H1 norm = 5.56
 Solving time step 22, time=0.5750...
H1 norm = 5.4
  Refining the mesh...
H1 norm = 5.4
 Solving time step 23, time=0.6000...
H1 norm = 5.24
  Refining the mesh...
H1 norm = 5.24
 Solving time step 24, time=0.6250...
H1 norm = 5.09
  Refining the mesh...
H1 norm = 5.09
Final H1 norm = 5.09


 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:46:12 2013                                                                          |
| OS:             Linux                                                                                             |
| HostName:       spark.ices.utexas.edu                                                                             |
| OS Release:     2.6.32-279.22.1.el6.x86_64                                                                        |
| OS Version:     #1 SMP Tue Feb 5 14:33:39 CST 2013                                                                |
| Machine:        x86_64                                                                                            |
| Username:       roystgnr                                                                                          |
| Configuration:  ../configure  '--enable-everything'                                                               |
|  'METHODS=devel'                                                                                                  |
|  '--prefix=/h2/roystgnr/libmesh-test'                                                                             |
|  'CXX=distcc /usr/bin/g++'                                                                                        |
|  'CC=distcc /usr/bin/gcc'                                                                                         |
|  'FC=distcc /usr/bin/gfortran'                                                                                    |
|  'F77=distcc /usr/bin/gfortran'                                                                                   |
|  'PETSC_DIR=/opt/apps/ossw/libraries/petsc/petsc-3.3-p2'                                                          |
|  'PETSC_ARCH=gcc-system-mkl-gf-10.3.12.361-mpich2-1.4.1p1-cxx-opt'                                                |
|  'SLEPC_DIR=/opt/apps/ossw/libraries/slepc/slepc-3.3-p2-petsc-3.3-p2-cxx-opt'                                     |
|  'TRILINOS_DIR=/opt/apps/ossw/libraries/trilinos/trilinos-10.12.2/sl6/gcc-system/mpich2-1.4.1p1/mkl-gf-10.3.12.361'|
|  'VTK_DIR=/opt/apps/ossw/libraries/vtk/vtk-5.10.0/sl6/gcc-system'                                                 |
|  'HDF5_DIR=/opt/apps/ossw/libraries/hdf5/hdf5-1.8.9/sl6/gcc-system'                                               |
 -------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=3.86501, Active time=3.70392                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     26        0.0220      0.000847    0.0301      0.001158    0.59     0.81     |
|   build_constraint_matrix()        7380      0.0142      0.000002    0.0142      0.000002    0.38     0.38     |
|   build_sparsity()                 26        0.0194      0.000746    0.0522      0.002009    0.52     1.41     |
|   cnstrn_elem_mat_vec()            7380      0.0185      0.000003    0.0185      0.000003    0.50     0.50     |
|   create_dof_constraints()         26        0.2381      0.009156    0.4738      0.018223    6.43     12.79    |
|   distribute_dofs()                26        0.0527      0.002028    0.1703      0.006549    1.42     4.60     |
|   dof_indices()                    63275     0.2806      0.000004    0.2806      0.000004    7.57     7.57     |
|   enforce_constraints_exactly()    76        0.0079      0.000104    0.0079      0.000104    0.21     0.21     |
|   old_dof_indices()                23508     0.1001      0.000004    0.1001      0.000004    2.70     2.70     |
|   prepare_send_list()              26        0.0002      0.000008    0.0002      0.000008    0.01     0.01     |
|   reinit()                         26        0.0855      0.003290    0.0855      0.003290    2.31     2.31     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          6         0.0087      0.001458    0.0253      0.004219    0.24     0.68     |
|   write()                          1         0.0009      0.000902    0.0028      0.002784    0.02     0.08     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               4         0.1606      0.040158    0.1606      0.040158    4.34     4.34     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        35116     0.0718      0.000002    0.0718      0.000002    1.94     1.94     |
|   init_shape_functions()           19268     0.0523      0.000003    0.0523      0.000003    1.41     1.41     |
|   inverse_map()                    73981     0.2201      0.000003    0.2201      0.000003    5.94     5.94     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             35116     0.0857      0.000002    0.0857      0.000002    2.31     2.31     |
|   compute_face_map()               9583      0.0527      0.000006    0.1211      0.000013    1.42     3.27     |
|   init_face_shape_functions()      1207      0.0035      0.000003    0.0035      0.000003    0.09     0.09     |
|   init_reference_to_physical_map() 19268     0.0777      0.000004    0.0777      0.000004    2.10     2.10     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               2         0.0398      0.019890    0.0399      0.019953    1.07     1.08     |
|                                                                                                                |
| JumpErrorEstimator                                                                                             |
|   estimate_error()                 25        0.1881      0.007526    0.7216      0.028866    5.08     19.48    |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           19140     0.0218      0.000001    0.0218      0.000001    0.59     0.59     |
|   init()                           55        0.0138      0.000251    0.0138      0.000251    0.37     0.37     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       25        0.0046      0.000184    0.0084      0.000335    0.12     0.23     |
|   find_neighbors()                 27        0.1031      0.003820    0.1182      0.004379    2.78     3.19     |
|   renumber_nodes_and_elem()        79        0.0113      0.000143    0.0113      0.000143    0.31     0.31     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0087      0.008710    0.0142      0.014193    0.24     0.38     |
|   compute_hilbert_indices()        28        0.0587      0.002095    0.0587      0.002095    1.58     1.58     |
|   find_global_indices()            28        0.0114      0.000407    0.0901      0.003216    0.31     2.43     |
|   parallel_sort()                  28        0.0025      0.000091    0.0161      0.000576    0.07     0.44     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         6         0.0002      0.000037    0.2266      0.037761    0.01     6.12     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              50        0.0078      0.000157    0.0095      0.000189    0.21     0.26     |
|   _refine_elements()               55        0.0385      0.000701    0.1047      0.001904    1.04     2.83     |
|   add_point()                      19140     0.0321      0.000002    0.0560      0.000003    0.87     1.51     |
|   make_coarsening_compatible()     111       0.1461      0.001316    0.1908      0.001719    3.94     5.15     |
|   make_flags_parallel_consistent() 75        0.0368      0.000490    0.0681      0.000908    0.99     1.84     |
|   make_refinement_compatible()     111       0.0069      0.000062    0.0094      0.000085    0.19     0.25     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0001      0.000080    0.0001      0.000080    0.00     0.00     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      27        0.1775      0.006575    0.2721      0.010077    4.79     7.35     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      143       0.0222      0.000155    0.0226      0.000158    0.60     0.61     |
|   barrier()                        1         0.0000      0.000012    0.0000      0.000012    0.00     0.00     |
|   broadcast()                      2         0.0000      0.000005    0.0000      0.000005    0.00     0.00     |
|   gather()                         14        0.0004      0.000032    0.0004      0.000032    0.01     0.01     |
|   max(bool)                        267       0.0322      0.000121    0.0322      0.000121    0.87     0.87     |
|   max(scalar)                      6520      0.0276      0.000004    0.0276      0.000004    0.75     0.75     |
|   max(vector)                      1648      0.0108      0.000007    0.0300      0.000018    0.29     0.81     |
|   min(bool)                        8252      0.0673      0.000008    0.0673      0.000008    1.82     1.82     |
|   min(scalar)                      6461      0.3005      0.000047    0.3005      0.000047    8.11     8.11     |
|   min(vector)                      1648      0.0119      0.000007    0.0354      0.000021    0.32     0.96     |
|   probe()                          1793      0.0212      0.000012    0.0212      0.000012    0.57     0.57     |
|   receive()                        1783      0.0053      0.000003    0.0265      0.000015    0.14     0.72     |
|   send()                           1744      0.0030      0.000002    0.0030      0.000002    0.08     0.08     |
|   send_receive()                   1800      0.0069      0.000004    0.0384      0.000021    0.19     1.04     |
|   sum()                            278       0.0409      0.000147    0.1594      0.000573    1.11     4.30     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           1750      0.0017      0.000001    0.0017      0.000001    0.04     0.04     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         27        0.0223      0.000826    0.0581      0.002152    0.60     1.57     |
|   set_parent_processor_ids()       27        0.0103      0.000383    0.0103      0.000383    0.28     0.28     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          50        0.1586      0.003171    0.1586      0.003171    4.28     4.28     |
|                                                                                                                |
| PointLocatorTree                                                                                               |
|   init(no master)                  50        0.0954      0.001907    0.1087      0.002174    2.57     2.93     |
|   operator()                       6165      0.0823      0.000013    0.1079      0.000017    2.22     2.91     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       75        0.0285      0.000380    0.1118      0.001491    0.77     3.02     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       50        0.1086      0.002172    0.2248      0.004495    2.93     6.07     |
|   calculate_norm()                 52        0.0507      0.000975    0.1636      0.003147    1.37     4.42     |
|   project_vector()                 76        0.1097      0.001444    0.3071      0.004041    2.96     8.29     |
|                                                                                                                |
| XdrIO                                                                                                          |
|   write()                          1         0.0024      0.002360    0.0041      0.004132    0.06     0.11     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            375016    3.7039                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example adaptivity_ex5:
*  mpirun -np 4 example-devel -n_timesteps 25 -n_refinements 5 -output_freq 10 -init_timestep 0 -exact_solution '10*exp(-(pow(x-0.8*t-0.2,2)+pow(y-0.8*t-0.2,2))/(0.01*(4*t+1)))/(4*t+1)' -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
***** Finished first 25 steps, now read in saved solution and continue *****
 
***************************************************************
* Running Example adaptivity_ex5:
*  mpirun -np 4 example-devel -read_solution -n_timesteps 25 -output_freq 10 -init_timestep 25 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
Usage:
	 /net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex5/.libs/lt-example-devel -init_timestep 0
OR
	 /net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex5/.libs/lt-example-devel -read_solution -init_timestep 26

Running: /net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex5/.libs/lt-example-devel -read_solution -n_timesteps 25 -output_freq 10 -init_timestep 25 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=737
    n_local_nodes()=203
  n_elem()=884
    n_local_elem()=234
    n_active_elem()=664
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
  n_threads()=1
  processor_id()=0

Initial H1 norm = 5.09105

 EquationSystems
  n_systems()=1
   System #0, "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=737
    n_local_dofs()=203
    n_constrained_dofs()=130
    n_local_constrained_dofs()=32
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 9.09362
      Average Off-Processor Bandwidth <= 1.1289
      Maximum  On-Processor Bandwidth <= 18
      Maximum Off-Processor Bandwidth <= 13
    DofMap Constraints
      Number of DoF Constraints = 130
      Average DoF Constraint Length= 1.92308
      Number of Node Constraints = 236
      Maximum Node Constraint Length= 5
      Average Node Constraint Length= 2.51271

 Solving time step 25, time=0.6500...
H1 norm = 4.95
  Refining the mesh...
H1 norm = 4.95
 Solving time step 26, time=0.6750...
H1 norm = 4.82
  Refining the mesh...
H1 norm = 4.82
 Solving time step 27, time=0.7000...
H1 norm = 4.69
  Refining the mesh...
H1 norm = 4.69
 Solving time step 28, time=0.7250...
H1 norm = 4.58
  Refining the mesh...
H1 norm = 4.58
 Solving time step 29, time=0.7500...
H1 norm = 4.46
  Refining the mesh...
H1 norm = 4.46
 Solving time step 30, time=0.7750...
H1 norm = 4.36
  Refining the mesh...
H1 norm = 4.36
 Solving time step 31, time=0.8000...
H1 norm = 4.25
  Refining the mesh...
H1 norm = 4.25
 Solving time step 32, time=0.8250...
H1 norm = 4.16
  Refining the mesh...
H1 norm = 4.16
 Solving time step 33, time=0.8500...
H1 norm = 4.06
  Refining the mesh...
H1 norm = 4.06
 Solving time step 34, time=0.8750...
H1 norm = 3.97
  Refining the mesh...
H1 norm = 3.97
 Solving time step 35, time=0.9000...
H1 norm = 3.89
  Refining the mesh...
H1 norm = 3.89
 Solving time step 36, time=0.9250...
H1 norm = 3.81
  Refining the mesh...
H1 norm = 3.81
 Solving time step 37, time=0.9500...
H1 norm = 3.73
  Refining the mesh...
H1 norm = 3.73
 Solving time step 38, time=0.9750...
H1 norm = 3.66
  Refining the mesh...
H1 norm = 3.66
 Solving time step 39, time=100000...
H1 norm = 3.58
  Refining the mesh...
H1 norm = 3.58
 Solving time step 40, time=1.0300...
H1 norm = 3.51
  Refining the mesh...
H1 norm = 3.51
 Solving time step 41, time=1.0500...
H1 norm = 3.45
  Refining the mesh...
H1 norm = 3.45
 Solving time step 42, time=1.0700...
H1 norm = 3.38
  Refining the mesh...
H1 norm = 3.38
 Solving time step 43, time=1.1000...
H1 norm = 3.32
  Refining the mesh...
H1 norm = 3.32
 Solving time step 44, time=1.1200...
H1 norm = 3.26
  Refining the mesh...
H1 norm = 3.26
 Solving time step 45, time=1.1500...
H1 norm = 3.21
  Refining the mesh...
H1 norm = 3.21
 Solving time step 46, time=1.1700...
H1 norm = 3.15
  Refining the mesh...
H1 norm = 3.15
 Solving time step 47, time=1.2000...
H1 norm = 3.1
  Refining the mesh...
H1 norm = 3.1
 Solving time step 48, time=1.2200...
H1 norm = 3.05
  Refining the mesh...
H1 norm = 3.05
 Solving time step 49, time=1.2500...
H1 norm = 3
  Refining the mesh...
H1 norm = 3

 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:46:16 2013                                                                          |
| OS:             Linux                                                                                             |
| HostName:       spark.ices.utexas.edu                                                                             |
| OS Release:     2.6.32-279.22.1.el6.x86_64                                                                        |
| OS Version:     #1 SMP Tue Feb 5 14:33:39 CST 2013                                                                |
| Machine:        x86_64                                                                                            |
| Username:       roystgnr                                                                                          |
| Configuration:  ../configure  '--enable-everything'                                                               |
|  'METHODS=devel'                                                                                                  |
|  '--prefix=/h2/roystgnr/libmesh-test'                                                                             |
|  'CXX=distcc /usr/bin/g++'                                                                                        |
|  'CC=distcc /usr/bin/gcc'                                                                                         |
|  'FC=distcc /usr/bin/gfortran'                                                                                    |
|  'F77=distcc /usr/bin/gfortran'                                                                                   |
|  'PETSC_DIR=/opt/apps/ossw/libraries/petsc/petsc-3.3-p2'                                                          |
|  'PETSC_ARCH=gcc-system-mkl-gf-10.3.12.361-mpich2-1.4.1p1-cxx-opt'                                                |
|  'SLEPC_DIR=/opt/apps/ossw/libraries/slepc/slepc-3.3-p2-petsc-3.3-p2-cxx-opt'                                     |
|  'TRILINOS_DIR=/opt/apps/ossw/libraries/trilinos/trilinos-10.12.2/sl6/gcc-system/mpich2-1.4.1p1/mkl-gf-10.3.12.361'|
|  'VTK_DIR=/opt/apps/ossw/libraries/vtk/vtk-5.10.0/sl6/gcc-system'                                                 |
|  'HDF5_DIR=/opt/apps/ossw/libraries/hdf5/hdf5-1.8.9/sl6/gcc-system'                                               |
 -------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=4.42663, Active time=4.2712                                                    |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     27        0.0232      0.000861    0.0326      0.001207    0.54     0.76     |
|   build_constraint_matrix()        10802     0.0170      0.000002    0.0170      0.000002    0.40     0.40     |
|   build_sparsity()                 27        0.0210      0.000776    0.0679      0.002513    0.49     1.59     |
|   cnstrn_elem_mat_vec()            10802     0.0274      0.000003    0.0274      0.000003    0.64     0.64     |
|   create_dof_constraints()         27        0.1828      0.006772    0.3157      0.011694    4.28     7.39     |
|   distribute_dofs()                27        0.0475      0.001759    0.1922      0.007117    1.11     4.50     |
|   dof_indices()                    84883     0.3012      0.000004    0.3012      0.000004    7.05     7.05     |
|   enforce_constraints_exactly()    78        0.0083      0.000107    0.0083      0.000107    0.19     0.19     |
|   old_dof_indices()                34203     0.1211      0.000004    0.1211      0.000004    2.83     2.83     |
|   prepare_send_list()              27        0.0002      0.000008    0.0002      0.000008    0.01     0.01     |
|   reinit()                         27        0.0811      0.003003    0.0811      0.003003    1.90     1.90     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          5         0.0022      0.000435    0.0102      0.002031    0.05     0.24     |
|   read()                           1         0.0033      0.003307    0.0355      0.035473    0.08     0.83     |
|   update()                         1         0.0001      0.000103    0.0001      0.000103    0.00     0.00     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               4         0.1616      0.040396    0.1616      0.040396    3.78     3.78     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        45354     0.0781      0.000002    0.0781      0.000002    1.83     1.83     |
|   init_shape_functions()           23685     0.0530      0.000002    0.0530      0.000002    1.24     1.24     |
|   inverse_map()                    78100     0.1896      0.000002    0.1896      0.000002    4.44     4.44     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             45354     0.0850      0.000002    0.0850      0.000002    1.99     1.99     |
|   compute_face_map()               11792     0.0526      0.000004    0.1237      0.000010    1.23     2.90     |
|   init_face_shape_functions()      662       0.0015      0.000002    0.0015      0.000002    0.04     0.04     |
|   init_reference_to_physical_map() 23685     0.0727      0.000003    0.0727      0.000003    1.70     1.70     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               1         0.0030      0.002986    0.0031      0.003058    0.07     0.07     |
|                                                                                                                |
| JumpErrorEstimator                                                                                             |
|   estimate_error()                 25        0.1967      0.007867    0.9685      0.038740    4.60     22.67    |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           3672      0.0036      0.000001    0.0036      0.000001    0.08     0.08     |
|   init()                           51        0.0159      0.000312    0.0159      0.000312    0.37     0.37     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       26        0.0038      0.000145    0.0077      0.000296    0.09     0.18     |
|   find_neighbors()                 26        0.0875      0.003365    0.1208      0.004647    2.05     2.83     |
|   renumber_nodes_and_elem()        78        0.0104      0.000133    0.0104      0.000133    0.24     0.24     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0087      0.008699    0.0144      0.014439    0.20     0.34     |
|   compute_hilbert_indices()        26        0.0640      0.002461    0.0640      0.002461    1.50     1.50     |
|   find_global_indices()            26        0.0111      0.000427    0.1110      0.004270    0.26     2.60     |
|   parallel_sort()                  26        0.0024      0.000094    0.0314      0.001208    0.06     0.74     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         5         0.0001      0.000029    0.1753      0.035070    0.00     4.11     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              51        0.0058      0.000113    0.0087      0.000171    0.13     0.20     |
|   _refine_elements()               51        0.0169      0.000331    0.0426      0.000835    0.40     1.00     |
|   add_point()                      3672      0.0058      0.000002    0.0098      0.000003    0.14     0.23     |
|   make_coarsening_compatible()     111       0.1593      0.001435    0.1686      0.001519    3.73     3.95     |
|   make_flags_parallel_consistent() 77        0.0429      0.000558    0.0861      0.001118    1.01     2.01     |
|   make_refinement_compatible()     111       0.0075      0.000068    0.0121      0.000109    0.18     0.28     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      26        0.1873      0.007204    0.3063      0.011779    4.39     7.17     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      140       0.0457      0.000326    0.0461      0.000329    1.07     1.08     |
|   broadcast()                      46        0.0001      0.000003    0.0001      0.000003    0.00     0.00     |
|   max(bool)                        264       0.0592      0.000224    0.0592      0.000224    1.39     1.39     |
|   max(scalar)                      6529      0.0296      0.000005    0.0296      0.000005    0.69     0.69     |
|   max(vector)                      1652      0.0104      0.000006    0.0297      0.000018    0.24     0.69     |
|   min(bool)                        8269      0.1045      0.000013    0.1045      0.000013    2.45     2.45     |
|   min(scalar)                      6472      0.6752      0.000104    0.6752      0.000104    15.81    15.81    |
|   min(vector)                      1652      0.0117      0.000007    0.0402      0.000024    0.27     0.94     |
|   probe()                          1756      0.0405      0.000023    0.0405      0.000023    0.95     0.95     |
|   receive()                        1750      0.0048      0.000003    0.0454      0.000026    0.11     1.06     |
|   send()                           1750      0.0025      0.000001    0.0025      0.000001    0.06     0.06     |
|   send_receive()                   1796      0.0062      0.000003    0.0561      0.000031    0.15     1.31     |
|   sum()                            267       0.0974      0.000365    0.4185      0.001567    2.28     9.80     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           1750      0.0015      0.000001    0.0015      0.000001    0.04     0.04     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         27        0.0198      0.000733    0.0899      0.003328    0.46     2.10     |
|   set_parent_processor_ids()       26        0.0106      0.000407    0.0106      0.000407    0.25     0.25     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          50        0.2675      0.005350    0.2675      0.005350    6.26     6.26     |
|                                                                                                                |
| PointLocatorTree                                                                                               |
|   init(no master)                  51        0.1118      0.002192    0.1288      0.002525    2.62     3.02     |
|   operator()                       3258      0.0186      0.000006    0.0293      0.000009    0.44     0.69     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       78        0.0393      0.000504    0.1639      0.002101    0.92     3.84     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       50        0.1298      0.002596    0.2719      0.005437    3.04     6.37     |
|   calculate_norm()                 51        0.0514      0.001009    0.2079      0.004076    1.20     4.87     |
|   project_vector()                 78        0.1688      0.002164    0.4148      0.005318    3.95     9.71     |
|                                                                                                                |
| XdrIO                                                                                                          |
|   read()                           1         0.0011      0.001053    0.0011      0.001120    0.02     0.03     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            415398    4.2712                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example adaptivity_ex5:
*  mpirun -np 4 example-devel -read_solution -n_timesteps 25 -output_freq 10 -init_timestep 25 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex5'
</pre>
</div>
<?php make_footer() ?>
</body>
</html>
<?php if (0) { ?>
\#Local Variables:
\#mode: html
\#End:
<?php } ?>
