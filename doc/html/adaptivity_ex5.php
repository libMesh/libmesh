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
          libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");
        
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
Create a new mesh.
ParallelMesh doesn't yet understand periodic BCs, plus
we still need some work on automatic parallel restarts
</div>

<div class ="fragment">
<pre>
          SerialMesh mesh;
        
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
              TransientLinearImplicitSystem &  system =
                equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;("Convection-Diffusion");
        
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
                  Real H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));
        
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
              TransientLinearImplicitSystem& system =
        	equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;
                  ("Convection-Diffusion");
              Real H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));
        
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

<a name="comments"></a> 
<br><br><br> <h1> The source file exact_solution.C with comments: </h1> 
<div class = "comment">
  

<br><br>This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
  

<br><br>This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
  

<br><br>You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


<br><br>

<br><br>

<br><br>C++ Includes
</div>

<div class ="fragment">
<pre>
        #include &lt;math.h&gt;
        
</pre>
</div>
<div class = "comment">
Mesh library includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh_common.h"
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        
        
        
        
        /**
         *
         */
        Real exact_solution (const Real x,
        		     const Real y,
        		     const Real t)
        {
          const Real xo = 0.2;
          const Real yo = 0.2;
          const Real u  = 0.8;
          const Real v  = 0.8;
         
          const Real num =
            pow(x - u*t - xo, 2.) +
            pow(y - v*t - yo, 2.);
        
          const Real den =
            0.01*(4.*t + 1.);
        
          return exp(-num/den)/(4.*t + 1.);
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
  
    libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  
  
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
  
    SerialMesh mesh;
  
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
        
        TransientLinearImplicitSystem &amp;  system =
          equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
        *system.old_local_solution = *system.current_local_solution;
        
        <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_r_steps = 1;
        <B><FONT COLOR="#A020F0">if</FONT></B>(command_line.search(<B><FONT COLOR="#BC8F8F">&quot;-max_r_steps&quot;</FONT></B>))
          max_r_steps = command_line.next(0);
        
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> r_step=0; r_step&lt;max_r_steps+1; r_step++)
          {
            system.solve();
  
            Real H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));
  
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
        TransientLinearImplicitSystem&amp; system =
  	equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;
            (<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
        Real H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));
  
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
<a name="nocomments"></a> 
<br><br><br> <h1> The source file exact_solution.C without comments: </h1> 
<pre> 
    
    
    
  
  
  
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_common.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  
  
  
  <I><FONT COLOR="#B22222">/**
   *
   */</FONT></I>
  Real exact_solution (<B><FONT COLOR="#228B22">const</FONT></B> Real x,
  		     <B><FONT COLOR="#228B22">const</FONT></B> Real y,
  		     <B><FONT COLOR="#228B22">const</FONT></B> Real t)
  {
    <B><FONT COLOR="#228B22">const</FONT></B> Real xo = 0.2;
    <B><FONT COLOR="#228B22">const</FONT></B> Real yo = 0.2;
    <B><FONT COLOR="#228B22">const</FONT></B> Real u  = 0.8;
    <B><FONT COLOR="#228B22">const</FONT></B> Real v  = 0.8;
   
    <B><FONT COLOR="#228B22">const</FONT></B> Real num =
      pow(x - u*t - xo, 2.) +
      pow(y - v*t - yo, 2.);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real den =
      0.01*(4.*t + 1.);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> exp(-num/den)/(4.*t + 1.);
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example adaptivity_ex5:
*  mpirun -np 12 example-devel -n_timesteps 25 -n_refinements 5 -output_freq 10 -init_timestep 0 -exact_solution '10*exp(-(pow(x-0.8*t-0.2,2)+pow(y-0.8*t-0.2,2))/(0.01*(4*t+1)))/(4*t+1)' -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Usage:
	 /workspace/libmesh/examples/adaptivity/adaptivity_ex5/.libs/lt-example-devel -init_timestep 0
OR
	 /workspace/libmesh/examples/adaptivity/adaptivity_ex5/.libs/lt-example-devel -read_solution -init_timestep 26

Running: /workspace/libmesh/examples/adaptivity/adaptivity_ex5/.libs/lt-example-devel -n_timesteps 25 -n_refinements 5 -output_freq 10 -init_timestep 0 -exact_solution 10*exp(-(pow(x-0.8*t-0.2,2) + pow(y-0.8*t-0.2,2))/(0.01*(4*t+1)))/(4*t+1) -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=4225
    n_local_nodes()=380
  n_elem()=5460
    n_local_elem()=478
    n_active_elem()=4096
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
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
    n_local_dofs()=380
    n_constrained_dofs()=129
    n_local_constrained_dofs()=20
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.66793
      Average Off-Processor Bandwidth <= 0.865325
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

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/adaptivity/adaptivity_ex5/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:00:35 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           1.171e+01      1.00016   1.171e+01
Objects:              2.571e+03      1.00000   2.571e+03
Flops:                2.855e+06      1.88692   2.179e+06  2.614e+07
Flops/sec:            2.438e+05      1.88662   1.860e+05  2.232e+06
MPI Messages:         1.019e+04      1.22473   9.007e+03  1.081e+05
MPI Message Lengths:  9.291e+05      1.18527   9.361e+01  1.012e+07
MPI Reductions:       5.931e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.1712e+01 100.0%  2.6143e+07 100.0%  1.081e+05 100.0%  9.361e+01      100.0%  5.930e+03 100.0% 

------------------------------------------------------------------------------------------------------------------------
See the 'Profiling' chapter of the users' manual for details on interpreting output.
Phase summary info:
   Count: number of times phase was executed
   Time and Flops: Max - maximum over all processors
                   Ratio - ratio of maximum to minimum over all processors
   Mess: number of messages sent
   Avg. len: average message length
   Reduct: number of global reductions
   Global: entire computation
   Stage: stages of a computation. Set stages with PetscLogStagePush() and PetscLogStagePop().
      %T - percent time in this phase         %f - percent flops in this phase
      %M - percent messages in this phase     %L - percent message lengths in this phase
      %R - percent reductions in this phase
   Total Mflop/s: 10e-6 * (sum of flops over all processors)/(max time over all processors)
------------------------------------------------------------------------------------------------------------------------
Event                Count      Time (sec)     Flops                             --- Global ---  --- Stage ---   Total
                   Max Ratio  Max     Ratio   Max  Ratio  Mess   Avg len Reduct  %T %f %M %L %R  %T %f %M %L %R Mflop/s
------------------------------------------------------------------------------------------------------------------------

--- Event Stage 0: Main Stage

VecMDot              310 1.0 3.6550e-03 1.3 1.46e+05 1.5 0.0e+00 0.0e+00 3.1e+02  0  5  0  0  5   0  5  0  0  5   390
VecNorm              410 1.0 1.3765e-02 3.9 5.37e+04 1.5 0.0e+00 0.0e+00 4.1e+02  0  2  0  0  7   0  2  0  0  7    38
VecScale             360 1.0 2.5821e-04 1.4 2.35e+04 1.5 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0   897
VecCopy              301 1.0 2.3317e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               896 1.0 6.0153e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              100 1.0 9.1240e-0363.6 1.34e+04 1.5 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0    14
VecMAXPY             360 1.0 2.4772e-04 1.4 1.87e+05 1.5 0.0e+00 0.0e+00 0.0e+00  0  7  0  0  0   0  7  0  0  0  7419
VecAssemblyBegin     804 1.0 8.4781e-02 2.3 0.00e+00 0.0 5.1e+03 7.1e+01 2.2e+03  1  0  5  4 37   1  0  5  4 37     0
VecAssemblyEnd       804 1.0 7.7271e-04 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     1039 1.0 4.6053e-03 1.1 0.00e+00 0.0 5.8e+04 1.0e+02 0.0e+00  0  0 54 58  0   0  0 54 58  0     0
VecScatterEnd       1039 1.0 9.9044e-03 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize         360 1.0 1.4027e-02 3.7 7.06e+04 1.5 0.0e+00 0.0e+00 3.6e+02  0  3  0  0  6   0  3  0  0  6    50
MatMult              360 1.0 9.9759e-03 1.7 4.26e+05 1.5 2.4e+04 5.9e+01 0.0e+00  0 16 22 14  0   0 16 22 14  0   431
MatSolve             410 1.0 1.6053e-03 1.7 1.18e+06 1.9 0.0e+00 0.0e+00 0.0e+00  0 41  0  0  0   0 41  0  0  0  6663
MatLUFactorNum        50 1.0 2.6062e-03 1.9 8.59e+05 2.6 0.0e+00 0.0e+00 0.0e+00  0 27  0  0  0   0 27  0  0  0  2681
MatILUFactorSym       50 1.0 6.2077e-03 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 1.5e+02  0  0  0  0  3   0  0  0  0  3     0
MatAssemblyBegin     100 1.0 5.9742e-02 1.8 0.00e+00 0.0 5.7e+03 2.5e+02 2.0e+02  0  0  5 14  3   0  0  5 14  3     0
MatAssemblyEnd       100 1.0 6.4476e-03 1.1 0.00e+00 0.0 3.4e+03 1.7e+01 2.1e+02  0  0  3  1  4   0  0  3  1  4     0
MatGetRowIJ           50 1.0 2.0027e-05 2.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        50 1.0 1.4875e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+02  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries       102 1.0 1.4424e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       310 1.0 4.1094e-03 1.3 2.92e+05 1.5 0.0e+00 0.0e+00 3.1e+02  0 11  0  0  5   0 11  0  0  5   698
KSPSetUp             100 1.0 1.0464e-03 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              50 1.0 5.3018e-02 1.0 2.86e+06 1.9 2.4e+04 5.9e+01 1.0e+03  0100 22 14 17   0100 22 14 17   493
PCSetUp              100 1.0 2.0391e-02 1.2 8.59e+05 2.6 0.0e+00 0.0e+00 3.0e+02  0 27  0  0  5   0 27  0  0  5   343
PCSetUpOnBlocks       50 1.0 1.6682e-02 1.3 8.59e+05 2.6 0.0e+00 0.0e+00 2.5e+02  0 27  0  0  4   0 27  0  0  4   419
PCApply              410 1.0 6.1698e-03 1.1 1.18e+06 1.9 0.0e+00 0.0e+00 0.0e+00  0 41  0  0  0   0 41  0  0  0  1734
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector  1184           1184      2994416     0
      Vector Scatter   359            359       371924     0
           Index Set   665            665       537496     0
   IS L to G Mapping   130            130        73320     0
              Matrix   128            128      1439864     0
       Krylov Solver    52             52       503360     0
      Preconditioner    52             52        46384     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.62532e-06
Average time for zero size MPI_Send(): 1.25766e-05
#PETSc Option Table entries:
-exact_solution 10*exp(-(pow(x-0.8*t-0.2,2) + pow(y-0.8*t-0.2,2))/(0.01*(4*t+1)))/(4*t+1)
-init_timestep 0
-ksp_right_pc
-log_summary
-n_refinements 5
-n_timesteps 25
-output_freq 10
-pc_type bjacobi
-sub_pc_factor_levels 4
-sub_pc_factor_zeropivot 0
-sub_pc_type ilu
#End of PETSc Option Table entries
Compiled without FORTRAN kernels
Compiled with full precision matrices (default)
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8 sizeof(PetscInt) 4
Configure run at: Thu Nov  8 11:21:02 2012
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared-libraries=1 --with-mpi-dir=/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1 --with-mumps=true --download-mumps=1 --with-metis=true --download-metis=1 --with-parmetis=true --download-parmetis=1 --with-superlu=true --download-superlu=1 --with-superludir=true --download-superlu_dist=1 --with-blacs=true --download-blacs=1 --with-scalapack=true --download-scalapack=1 --with-hypre=true --download-hypre=1 --with-blas-lib="[/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_intel_lp64.so,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_sequential.so,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_core.so]" --with-lapack-lib="[/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_lapack95_lp64.a]"
-----------------------------------------
Libraries compiled on Thu Nov  8 11:21:02 2012 on daedalus.ices.utexas.edu 
Machine characteristics: Linux-2.6.32-279.1.1.el6.x86_64-x86_64-with-redhat-6.3-Carbon
Using PETSc directory: /opt/apps/ossw/libraries/petsc/petsc-3.3-p2
Using PETSc arch: intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt
-----------------------------------------

Using C compiler: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpicxx  -wd1572 -O3   -fPIC   ${COPTFLAGS} ${CFLAGS}
Using Fortran compiler: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpif90  -fPIC -O3   ${FOPTFLAGS} ${FFLAGS} 
-----------------------------------------

Using include paths: -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/include -I/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/include
-----------------------------------------

Using C linker: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpicxx
Using Fortran linker: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpif90
Using libraries: -Wl,-rpath,/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -L/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -lpetsc -lX11 -Wl,-rpath,/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -L/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lHYPRE -lpthread -lsuperlu_dist_3.0 -lparmetis -lmetis -lscalapack -lblacs -lsuperlu_4.3 -Wl,-rpath,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64 -L/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,-rpath,/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/lib -L/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/lib -Wl,-rpath,/opt/apps/sysnet/intel/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64 -L/opt/apps/sysnet/intel/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -lmpichf90 -lifport -lifcore -lm -lm -lmpichcxx -ldl -lmpich -lopa -lmpl -lrt -lpthread -limf -lsvml -lipgo -ldecimal -lcilkrts -lstdc++ -lgcc_s -lirc -lirc_s -ldl 
-----------------------------------------


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 22:00:35 2013                                                                             |
| OS:             Linux                                                                                                |
| HostName:       hbar.ices.utexas.edu                                                                                 |
| OS Release:     2.6.32-279.1.1.el6.x86_64                                                                            |
| OS Version:     #1 SMP Tue Jul 10 11:24:23 CDT 2012                                                                  |
| Machine:        x86_64                                                                                               |
| Username:       benkirk                                                                                              |
| Configuration:  ./configure  '--enable-everything'                                                                   |
|  '--prefix=/workspace/libmesh/install'                                                                               |
|  'CXX=icpc'                                                                                                          |
|  'CC=icc'                                                                                                            |
|  'FC=ifort'                                                                                                          |
|  'F77=ifort'                                                                                                         |
|  'PETSC_DIR=/opt/apps/ossw/libraries/petsc/petsc-3.3-p2'                                                             |
|  'PETSC_ARCH=intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt'                                                |
|  'SLEPC_DIR=/opt/apps/ossw/libraries/slepc/slepc-3.3-p2-petsc-3.3-p2-cxx-opt'                                        |
|  'TRILINOS_DIR=/opt/apps/ossw/libraries/trilinos/trilinos-10.12.2/sl6/intel-12.1/mpich2-1.4.1p1/mkl-intel-10.3.12.361'|
|  'VTK_DIR=/opt/apps/ossw/libraries/vtk/vtk-5.10.0/sl6/intel-12.1'                                                    |
 ----------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=11.7923, Active time=11.3514                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     26        0.0918      0.003532    0.1747      0.006718    0.81     1.54     |
|   build_constraint_matrix()        2458      0.0223      0.000009    0.0223      0.000009    0.20     0.20     |
|   build_sparsity()                 26        0.0583      0.002242    0.1899      0.007302    0.51     1.67     |
|   cnstrn_elem_mat_vec()            2458      0.0216      0.000009    0.0216      0.000009    0.19     0.19     |
|   create_dof_constraints()         26        1.4636      0.056294    2.5433      0.097819    12.89    22.41    |
|   distribute_dofs()                26        0.3813      0.014665    1.1947      0.045950    3.36     10.52    |
|   dof_indices()                    29167     1.4005      0.000048    1.4005      0.000048    12.34    12.34    |
|   enforce_constraints_exactly()    76        0.0132      0.000173    0.0132      0.000173    0.12     0.12     |
|   old_dof_indices()                7791      0.4257      0.000055    0.4257      0.000055    3.75     3.75     |
|   prepare_send_list()              26        0.0012      0.000044    0.0012      0.000044    0.01     0.01     |
|   reinit()                         26        0.7592      0.029201    0.7592      0.029201    6.69     6.69     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          6         0.0143      0.002387    0.0658      0.010961    0.13     0.58     |
|   write()                          1         0.0073      0.007297    0.0097      0.009730    0.06     0.09     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               4         0.0489      0.012216    0.0489      0.012216    0.43     0.43     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        13238     0.1209      0.000009    0.1209      0.000009    1.06     1.06     |
|   init_shape_functions()           8030      0.0916      0.000011    0.0916      0.000011    0.81     0.81     |
|   inverse_map()                    50606     0.4091      0.000008    0.4091      0.000008    3.60     3.60     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             13238     0.1371      0.000010    0.1371      0.000010    1.21     1.21     |
|   compute_face_map()               3964      0.0846      0.000021    0.1642      0.000041    0.75     1.45     |
|   init_face_shape_functions()      1207      0.0151      0.000013    0.0151      0.000013    0.13     0.13     |
|   init_reference_to_physical_map() 8030      0.1475      0.000018    0.1475      0.000018    1.30     1.30     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               2         0.0426      0.021286    0.0428      0.021411    0.38     0.38     |
|                                                                                                                |
| JumpErrorEstimator                                                                                             |
|   estimate_error()                 25        0.2240      0.008961    1.1359      0.045435    1.97     10.01    |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           19140     0.0788      0.000004    0.0788      0.000004    0.69     0.69     |
|   init()                           55        0.0595      0.001081    0.0595      0.001081    0.52     0.52     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       25        0.0259      0.001035    0.0428      0.001714    0.23     0.38     |
|   find_neighbors()                 27        0.6884      0.025496    0.6994      0.025903    6.06     6.16     |
|   renumber_nodes_and_elem()        79        0.0399      0.000505    0.0399      0.000505    0.35     0.35     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0335      0.033486    0.0359      0.035944    0.29     0.32     |
|   compute_hilbert_indices()        28        0.1421      0.005074    0.1421      0.005074    1.25     1.25     |
|   find_global_indices()            28        0.0711      0.002541    0.2598      0.009280    0.63     2.29     |
|   parallel_sort()                  28        0.0194      0.000693    0.0269      0.000959    0.17     0.24     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         6         0.0004      0.000072    0.1586      0.026436    0.00     1.40     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              50        0.0204      0.000407    0.0213      0.000425    0.18     0.19     |
|   _refine_elements()               55        0.1624      0.002952    0.3336      0.006066    1.43     2.94     |
|   add_point()                      19140     0.0779      0.000004    0.1606      0.000008    0.69     1.41     |
|   make_coarsening_compatible()     111       0.4720      0.004253    0.6205      0.005590    4.16     5.47     |
|   make_refinement_compatible()     111       0.0229      0.000206    0.0250      0.000226    0.20     0.22     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0004      0.000354    0.0004      0.000354    0.00     0.00     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      27        1.3910      0.051517    1.6538      0.061250    12.25    14.57    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      143       0.0046      0.000032    0.0055      0.000039    0.04     0.05     |
|   barrier()                        1         0.0000      0.000022    0.0000      0.000022    0.00     0.00     |
|   broadcast()                      2         0.0000      0.000011    0.0000      0.000011    0.00     0.00     |
|   gather()                         14        0.0004      0.000025    0.0004      0.000025    0.00     0.00     |
|   max(bool)                        267       0.0550      0.000206    0.0550      0.000206    0.48     0.48     |
|   max(scalar)                      5720      0.0440      0.000008    0.0440      0.000008    0.39     0.39     |
|   max(vector)                      1448      0.0196      0.000014    0.0492      0.000034    0.17     0.43     |
|   min(bool)                        7177      0.0530      0.000007    0.0530      0.000007    0.47     0.47     |
|   min(scalar)                      5661      0.3905      0.000069    0.3905      0.000069    3.44     3.44     |
|   min(vector)                      1448      0.0206      0.000014    0.0522      0.000036    0.18     0.46     |
|   probe()                          3253      0.0348      0.000011    0.0348      0.000011    0.31     0.31     |
|   receive()                        3227      0.0215      0.000007    0.0568      0.000018    0.19     0.50     |
|   send()                           3084      0.0100      0.000003    0.0100      0.000003    0.09     0.09     |
|   send_receive()                   3140      0.0254      0.000008    0.0993      0.000032    0.22     0.87     |
|   sum()                            278       0.0147      0.000053    0.0227      0.000082    0.13     0.20     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           3106      0.0062      0.000002    0.0062      0.000002    0.05     0.05     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         27        0.0708      0.002623    0.0951      0.003523    0.62     0.84     |
|   set_parent_processor_ids()       27        0.0553      0.002050    0.0553      0.002050    0.49     0.49     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          50        0.1289      0.002578    0.1289      0.002578    1.14     1.14     |
|                                                                                                                |
| PointLocatorTree                                                                                               |
|   init(no master)                  50        0.4769      0.009538    0.4963      0.009926    4.20     4.37     |
|   operator()                       5647      0.2610      0.000046    0.3209      0.000057    2.30     2.83     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       75        0.0411      0.000548    0.4015      0.005353    0.36     3.54     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       50        0.1154      0.002307    0.3746      0.007492    1.02     3.30     |
|   calculate_norm()                 52        0.0609      0.001172    0.2996      0.005762    0.54     2.64     |
|   project_vector()                 76        0.1454      0.001913    0.8453      0.011122    1.28     7.45     |
|                                                                                                                |
| XdrIO                                                                                                          |
|   write()                          1         0.0081      0.008146    0.0139      0.013874    0.07     0.12     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            223393    11.3514                                         100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example adaptivity_ex5:
*  mpirun -np 12 example-devel -n_timesteps 25 -n_refinements 5 -output_freq 10 -init_timestep 0 -exact_solution '10*exp(-(pow(x-0.8*t-0.2,2)+pow(y-0.8*t-0.2,2))/(0.01*(4*t+1)))/(4*t+1)' -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
***** Finished first 25 steps, now read in saved solution and continue *****
 
***************************************************************
* Running Example adaptivity_ex5:
*  mpirun -np 12 example-devel -read_solution -n_timesteps 25 -output_freq 10 -init_timestep 25 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Usage:
	 /workspace/libmesh/examples/adaptivity/adaptivity_ex5/.libs/lt-example-devel -init_timestep 0
OR
	 /workspace/libmesh/examples/adaptivity/adaptivity_ex5/.libs/lt-example-devel -read_solution -init_timestep 26

Running: /workspace/libmesh/examples/adaptivity/adaptivity_ex5/.libs/lt-example-devel -read_solution -n_timesteps 25 -output_freq 10 -init_timestep 25 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=737
    n_local_nodes()=76
  n_elem()=884
    n_local_elem()=79
    n_active_elem()=664
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
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
    n_local_dofs()=76
    n_constrained_dofs()=132
    n_local_constrained_dofs()=15
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.70556
      Average Off-Processor Bandwidth <= 1.95387
      Maximum  On-Processor Bandwidth <= 20
      Maximum Off-Processor Bandwidth <= 18
    DofMap Constraints
      Number of DoF Constraints = 132
      Average DoF Constraint Length= 1.90909
      Number of Node Constraints = 237
      Maximum Node Constraint Length= 5
      Average Node Constraint Length= 2.50211

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
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/adaptivity/adaptivity_ex5/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:00:49 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           1.307e+01      1.00000   1.307e+01
Objects:              2.613e+03      1.00000   2.613e+03
Flops:                4.802e+06      1.57138   3.818e+06  4.582e+07
Flops/sec:            3.674e+05      1.57138   2.921e+05  3.505e+06
MPI Messages:         9.833e+03      1.13973   9.116e+03  1.094e+05
MPI Message Lengths:  1.118e+06      1.13050   1.154e+02  1.262e+07
MPI Reductions:       5.998e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.3071e+01 100.0%  4.5816e+07 100.0%  1.094e+05 100.0%  1.154e+02      100.0%  5.997e+03 100.0% 

------------------------------------------------------------------------------------------------------------------------
See the 'Profiling' chapter of the users' manual for details on interpreting output.
Phase summary info:
   Count: number of times phase was executed
   Time and Flops: Max - maximum over all processors
                   Ratio - ratio of maximum to minimum over all processors
   Mess: number of messages sent
   Avg. len: average message length
   Reduct: number of global reductions
   Global: entire computation
   Stage: stages of a computation. Set stages with PetscLogStagePush() and PetscLogStagePop().
      %T - percent time in this phase         %f - percent flops in this phase
      %M - percent messages in this phase     %L - percent message lengths in this phase
      %R - percent reductions in this phase
   Total Mflop/s: 10e-6 * (sum of flops over all processors)/(max time over all processors)
------------------------------------------------------------------------------------------------------------------------
Event                Count      Time (sec)     Flops                             --- Global ---  --- Stage ---   Total
                   Max Ratio  Max     Ratio   Max  Ratio  Mess   Avg len Reduct  %T %f %M %L %R  %T %f %M %L %R Mflop/s
------------------------------------------------------------------------------------------------------------------------

--- Event Stage 0: Main Stage

VecMDot              286 1.0 4.8351e-03 1.5 1.87e+05 1.4 0.0e+00 0.0e+00 2.9e+02  0  4  0  0  5   0  4  0  0  5   388
VecNorm              386 1.0 1.6516e-02 3.8 7.23e+04 1.4 0.0e+00 0.0e+00 3.9e+02  0  2  0  0  6   0  2  0  0  6    44
VecScale             336 1.0 2.7943e-04 1.6 3.15e+04 1.4 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  1129
VecCopy              308 1.0 2.5463e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               882 1.0 6.0797e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              100 1.0 8.7740e-0365.5 1.88e+04 1.4 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0    21
VecMAXPY             336 1.0 2.6584e-04 1.4 2.42e+05 1.4 0.0e+00 0.0e+00 0.0e+00  0  5  0  0  0   0  5  0  0  0  9113
VecAssemblyBegin     827 1.0 1.0912e-01 1.7 0.00e+00 0.0 5.2e+03 9.1e+01 2.2e+03  1  0  5  4 37   1  0  5  4 37     0
VecAssemblyEnd       827 1.0 8.3995e-04 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     1035 1.0 4.5924e-03 1.1 0.00e+00 0.0 5.9e+04 1.3e+02 0.0e+00  0  0 54 58  0   0  0 54 58  0     0
VecScatterEnd       1035 1.0 4.7554e-02 3.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize         336 1.0 1.6325e-02 3.8 9.44e+04 1.4 0.0e+00 0.0e+00 3.4e+02  0  2  0  0  6   0  2  0  0  6    58
MatMult              336 1.0 1.6242e-02 1.6 5.55e+05 1.3 2.3e+04 6.9e+01 0.0e+00  0 13 21 13  0   0 13 21 13  0   354
MatSolve             386 1.0 2.4102e-03 1.5 1.89e+06 1.5 0.0e+00 0.0e+00 0.0e+00  0 40  0  0  0   0 40  0  0  0  7602
MatLUFactorNum        50 1.0 4.3344e-03 1.6 1.85e+06 1.8 0.0e+00 0.0e+00 0.0e+00  0 35  0  0  0   0 35  0  0  0  3742
MatILUFactorSym       50 1.0 1.0225e-02 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 1.5e+02  0  0  0  0  3   0  0  0  0  3     0
MatAssemblyBegin     100 1.0 1.0617e-01 2.6 0.00e+00 0.0 5.7e+03 2.9e+02 2.0e+02  1  0  5 13  3   1  0  5 13  3     0
MatAssemblyEnd       100 1.0 6.7623e-03 1.1 0.00e+00 0.0 3.6e+03 1.9e+01 2.1e+02  0  0  3  1  3   0  0  3  1  3     0
MatGetRowIJ           50 1.0 1.9550e-05 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        50 1.0 1.4961e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+02  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries       104 1.0 1.5020e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       286 1.0 5.3520e-03 1.4 3.75e+05 1.4 0.0e+00 0.0e+00 2.9e+02  0  8  0  0  5   0  8  0  0  5   703
KSPSetUp             100 1.0 9.7060e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              50 1.0 6.5483e-02 1.0 4.80e+06 1.6 2.3e+04 6.9e+01 9.7e+02  0100 21 13 16   0100 21 13 16   700
PCSetUp              100 1.0 2.6102e-02 1.2 1.85e+06 1.8 0.0e+00 0.0e+00 3.0e+02  0 35  0  0  5   0 35  0  0  5   621
PCSetUpOnBlocks       50 1.0 2.2420e-02 1.3 1.85e+06 1.8 0.0e+00 0.0e+00 2.5e+02  0 35  0  0  4   0 35  0  0  4   723
PCApply              386 1.0 6.6712e-03 1.1 1.89e+06 1.5 0.0e+00 0.0e+00 0.0e+00  0 40  0  0  0   0 40  0  0  0  2747
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector  1203           1203      3417904     0
      Vector Scatter   366            366       379176     0
           Index Set   675            675       553908     0
   IS L to G Mapping   133            133        75012     0
              Matrix   131            131      2090140     0
       Krylov Solver    52             52       503360     0
      Preconditioner    52             52        46384     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 3.57628e-06
Average time for zero size MPI_Send(): 1.45833e-05
#PETSc Option Table entries:
-init_timestep 25
-ksp_right_pc
-log_summary
-n_timesteps 25
-output_freq 10
-pc_type bjacobi
-read_solution
-sub_pc_factor_levels 4
-sub_pc_factor_zeropivot 0
-sub_pc_type ilu
#End of PETSc Option Table entries
Compiled without FORTRAN kernels
Compiled with full precision matrices (default)
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8 sizeof(PetscInt) 4
Configure run at: Thu Nov  8 11:21:02 2012
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared-libraries=1 --with-mpi-dir=/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1 --with-mumps=true --download-mumps=1 --with-metis=true --download-metis=1 --with-parmetis=true --download-parmetis=1 --with-superlu=true --download-superlu=1 --with-superludir=true --download-superlu_dist=1 --with-blacs=true --download-blacs=1 --with-scalapack=true --download-scalapack=1 --with-hypre=true --download-hypre=1 --with-blas-lib="[/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_intel_lp64.so,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_sequential.so,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_core.so]" --with-lapack-lib="[/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64/libmkl_lapack95_lp64.a]"
-----------------------------------------
Libraries compiled on Thu Nov  8 11:21:02 2012 on daedalus.ices.utexas.edu 
Machine characteristics: Linux-2.6.32-279.1.1.el6.x86_64-x86_64-with-redhat-6.3-Carbon
Using PETSc directory: /opt/apps/ossw/libraries/petsc/petsc-3.3-p2
Using PETSc arch: intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt
-----------------------------------------

Using C compiler: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpicxx  -wd1572 -O3   -fPIC   ${COPTFLAGS} ${CFLAGS}
Using Fortran compiler: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpif90  -fPIC -O3   ${FOPTFLAGS} ${FFLAGS} 
-----------------------------------------

Using include paths: -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/include -I/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/include -I/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/include
-----------------------------------------

Using C linker: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpicxx
Using Fortran linker: /opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/bin/mpif90
Using libraries: -Wl,-rpath,/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -L/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -lpetsc -lX11 -Wl,-rpath,/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -L/opt/apps/ossw/libraries/petsc/petsc-3.3-p2/intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt/lib -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lHYPRE -lpthread -lsuperlu_dist_3.0 -lparmetis -lmetis -lscalapack -lblacs -lsuperlu_4.3 -Wl,-rpath,/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64 -L/opt/apps/sysnet/intel/12.1/mkl/10.3.12.361/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,-rpath,/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/lib -L/opt/apps/ossw/libraries/mpich2/mpich2-1.4.1p1/sl6/intel-12.1/lib -Wl,-rpath,/opt/apps/sysnet/intel/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64 -L/opt/apps/sysnet/intel/12.1/composer_xe_2011_sp1.7.256/compiler/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -lmpichf90 -lifport -lifcore -lm -lm -lmpichcxx -ldl -lmpich -lopa -lmpl -lrt -lpthread -limf -lsvml -lipgo -ldecimal -lcilkrts -lstdc++ -lgcc_s -lirc -lirc_s -ldl 
-----------------------------------------


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 22:00:49 2013                                                                             |
| OS:             Linux                                                                                                |
| HostName:       hbar.ices.utexas.edu                                                                                 |
| OS Release:     2.6.32-279.1.1.el6.x86_64                                                                            |
| OS Version:     #1 SMP Tue Jul 10 11:24:23 CDT 2012                                                                  |
| Machine:        x86_64                                                                                               |
| Username:       benkirk                                                                                              |
| Configuration:  ./configure  '--enable-everything'                                                                   |
|  '--prefix=/workspace/libmesh/install'                                                                               |
|  'CXX=icpc'                                                                                                          |
|  'CC=icc'                                                                                                            |
|  'FC=ifort'                                                                                                          |
|  'F77=ifort'                                                                                                         |
|  'PETSC_DIR=/opt/apps/ossw/libraries/petsc/petsc-3.3-p2'                                                             |
|  'PETSC_ARCH=intel-12.1-mkl-intel-10.3.12.361-mpich2-1.4.1p1-cxx-opt'                                                |
|  'SLEPC_DIR=/opt/apps/ossw/libraries/slepc/slepc-3.3-p2-petsc-3.3-p2-cxx-opt'                                        |
|  'TRILINOS_DIR=/opt/apps/ossw/libraries/trilinos/trilinos-10.12.2/sl6/intel-12.1/mpich2-1.4.1p1/mkl-intel-10.3.12.361'|
|  'VTK_DIR=/opt/apps/ossw/libraries/vtk/vtk-5.10.0/sl6/intel-12.1'                                                    |
 ----------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=13.0955, Active time=12.6987                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     27        0.1210      0.004483    0.2279      0.008442    0.95     1.79     |
|   build_constraint_matrix()        3608      0.0301      0.000008    0.0301      0.000008    0.24     0.24     |
|   build_sparsity()                 27        0.0766      0.002836    0.2394      0.008865    0.60     1.88     |
|   cnstrn_elem_mat_vec()            3608      0.0242      0.000007    0.0242      0.000007    0.19     0.19     |
|   create_dof_constraints()         27        1.5439      0.057180    2.3989      0.088849    12.16    18.89    |
|   distribute_dofs()                27        0.4382      0.016230    1.4267      0.052840    3.45     11.23    |
|   dof_indices()                    39169     1.8048      0.000046    1.8048      0.000046    14.21    14.21    |
|   enforce_constraints_exactly()    78        0.0159      0.000204    0.0159      0.000204    0.13     0.13     |
|   old_dof_indices()                11289     0.6097      0.000054    0.6097      0.000054    4.80     4.80     |
|   prepare_send_list()              27        0.0015      0.000055    0.0015      0.000055    0.01     0.01     |
|   reinit()                         27        0.9183      0.034010    0.9183      0.034010    7.23     7.23     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          5         0.0058      0.001162    0.0344      0.006876    0.05     0.27     |
|   read()                           1         0.0142      0.014153    0.1363      0.136318    0.11     1.07     |
|   update()                         1         0.0002      0.000226    0.0002      0.000226    0.00     0.00     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               4         0.0300      0.007503    0.0300      0.007503    0.24     0.24     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        15941     0.1536      0.000010    0.1536      0.000010    1.21     1.21     |
|   init_shape_functions()           8771      0.0997      0.000011    0.0997      0.000011    0.79     0.79     |
|   inverse_map()                    47542     0.3472      0.000007    0.3472      0.000007    2.73     2.73     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             15941     0.1489      0.000009    0.1489      0.000009    1.17     1.17     |
|   compute_face_map()               4335      0.0809      0.000019    0.1637      0.000038    0.64     1.29     |
|   init_face_shape_functions()      662       0.0080      0.000012    0.0080      0.000012    0.06     0.06     |
|   init_reference_to_physical_map() 8771      0.1420      0.000016    0.1420      0.000016    1.12     1.12     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               1         0.0063      0.006285    0.0064      0.006396    0.05     0.05     |
|                                                                                                                |
| JumpErrorEstimator                                                                                             |
|   estimate_error()                 25        0.2809      0.011234    1.3537      0.054146    2.21     10.66    |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           3672      0.0159      0.000004    0.0159      0.000004    0.12     0.12     |
|   init()                           51        0.0826      0.001619    0.0826      0.001619    0.65     0.65     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       26        0.0158      0.000608    0.0319      0.001226    0.12     0.25     |
|   find_neighbors()                 26        0.7182      0.027623    0.7623      0.029320    5.66     6.00     |
|   renumber_nodes_and_elem()        78        0.0404      0.000518    0.0404      0.000518    0.32     0.32     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0330      0.033022    0.0360      0.036044    0.26     0.28     |
|   compute_hilbert_indices()        26        0.1858      0.007145    0.1858      0.007145    1.46     1.46     |
|   find_global_indices()            26        0.0850      0.003271    0.3376      0.012983    0.67     2.66     |
|   parallel_sort()                  26        0.0210      0.000809    0.0365      0.001403    0.17     0.29     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         5         0.0003      0.000066    0.0726      0.014514    0.00     0.57     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              51        0.0180      0.000352    0.0197      0.000387    0.14     0.16     |
|   _refine_elements()               51        0.0716      0.001403    0.1138      0.002231    0.56     0.90     |
|   add_point()                      3672      0.0161      0.000004    0.0327      0.000009    0.13     0.26     |
|   make_coarsening_compatible()     111       0.5719      0.005152    0.6090      0.005486    4.50     4.80     |
|   make_refinement_compatible()     111       0.0303      0.000273    0.0349      0.000314    0.24     0.27     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      26        1.5471      0.059502    1.8922      0.072778    12.18    14.90    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      140       0.0225      0.000160    0.0231      0.000165    0.18     0.18     |
|   broadcast()                      46        0.0008      0.000017    0.0007      0.000015    0.01     0.01     |
|   max(bool)                        264       0.0548      0.000207    0.0548      0.000207    0.43     0.43     |
|   max(scalar)                      5705      0.0518      0.000009    0.0518      0.000009    0.41     0.41     |
|   max(vector)                      1446      0.0210      0.000015    0.0564      0.000039    0.17     0.44     |
|   min(bool)                        7162      0.0855      0.000012    0.0855      0.000012    0.67     0.67     |
|   min(scalar)                      5648      0.4198      0.000074    0.4198      0.000074    3.31     3.31     |
|   min(vector)                      1446      0.0224      0.000016    0.0610      0.000042    0.18     0.48     |
|   probe()                          3040      0.0457      0.000015    0.0457      0.000015    0.36     0.36     |
|   receive()                        3018      0.0201      0.000007    0.0664      0.000022    0.16     0.52     |
|   send()                           3018      0.0092      0.000003    0.0092      0.000003    0.07     0.07     |
|   send_receive()                   3048      0.0235      0.000008    0.1062      0.000035    0.18     0.84     |
|   sum()                            267       0.0377      0.000141    0.0507      0.000190    0.30     0.40     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           3018      0.0057      0.000002    0.0057      0.000002    0.04     0.04     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         27        0.0697      0.002582    0.1236      0.004578    0.55     0.97     |
|   set_parent_processor_ids()       26        0.0694      0.002669    0.0694      0.002669    0.55     0.55     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          50        0.1708      0.003417    0.1708      0.003417    1.35     1.35     |
|                                                                                                                |
| PointLocatorTree                                                                                               |
|   init(no master)                  51        0.6394      0.012537    0.6620      0.012981    5.03     5.21     |
|   operator()                       2942      0.0741      0.000025    0.1040      0.000035    0.58     0.82     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       78        0.0638      0.000818    0.6766      0.008674    0.50     5.33     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       50        0.1630      0.003259    0.5294      0.010589    1.28     4.17     |
|   calculate_norm()                 51        0.0746      0.001462    0.3912      0.007671    0.59     3.08     |
|   project_vector()                 78        0.1899      0.002435    1.2144      0.015569    1.50     9.56     |
|                                                                                                                |
| XdrIO                                                                                                          |
|   read()                           1         0.0085      0.008512    0.0088      0.008809    0.07     0.07     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            208493    12.6987                                         100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example adaptivity_ex5:
*  mpirun -np 12 example-devel -read_solution -n_timesteps 25 -output_freq 10 -init_timestep 25 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
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
