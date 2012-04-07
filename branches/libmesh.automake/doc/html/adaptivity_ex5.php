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
        #include &lt;cmath&gt;
        
</pre>
</div>
<div class = "comment">
Basic include file needed for the mesh functionality.
</div>

<div class ="fragment">
<pre>
        #include "libmesh.h"
        #include "mesh.h"
        #include "mesh_refinement.h"
        #include "gmv_io.h"
        #include "exodusII_io.h"
        #include "equation_systems.h"
        #include "fe.h"
        #include "quadrature_gauss.h"
        #include "dof_map.h"
        #include "sparse_matrix.h"
        #include "numeric_vector.h"
        #include "dense_matrix.h"
        #include "dense_vector.h"
        
        #include "periodic_boundaries.h"
        #include "mesh_generation.h"
        #include "parsed_function.h"
        
        #include "getpot.h"
        
</pre>
</div>
<div class = "comment">
Some (older) compilers do not offer full stream 
functionality, \p OStringStream works around this.
Check example 9 for details.
</div>

<div class ="fragment">
<pre>
        #include "o_string_stream.h"
        
</pre>
</div>
<div class = "comment">
This example will solve a linear transient system,
so we need to include the \p TransientLinearImplicitSystem definition.
</div>

<div class ="fragment">
<pre>
        #include "transient_system.h"
        #include "linear_implicit_system.h"
        #include "vector_value.h"
        
</pre>
</div>
<div class = "comment">
To refine the mesh we need an \p ErrorEstimator
object to figure out which elements to refine.
</div>

<div class ="fragment">
<pre>
        #include "error_vector.h"
        #include "kelly_error_estimator.h"
        
</pre>
</div>
<div class = "comment">
The definition of a geometric element
</div>

<div class ="fragment">
<pre>
        #include "elem.h"
        
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
        #elif defined(LIBMESH_ENABLE_PARMESH)
</pre>
</div>
<div class = "comment">
ParallelMesh doesn't yet understand periodic BCs, plus
we still need some work on automatic parallel restarts
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(false, "--disable-parmesh");
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
</div>

<div class ="fragment">
<pre>
          Mesh mesh;
        
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
Give the system a pointer to the matrix assembly
and initialization functions.
</div>

<div class ="fragment">
<pre>
              system.attach_assemble_function (assemble_cd);
              system.attach_init_function (init_cd);
        
</pre>
</div>
<div class = "comment">
Initialize the data structures for the equation system.
</div>

<div class ="fragment">
<pre>
              equation_systems.init ();
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
        
</pre>
</div>
<div class = "comment">
Get a reference to the system so that we can call update() on it
</div>

<div class ="fragment">
<pre>
              TransientLinearImplicitSystem & system = 
                equation_systems.get_system&lt;TransientLinearImplicitSystem&gt; 
                ("Convection-Diffusion");
        
</pre>
</div>
<div class = "comment">
We need to call update to put system in a consistent state
with the solution that was read in
</div>

<div class ="fragment">
<pre>
              system.update();
        
</pre>
</div>
<div class = "comment">
Attach the same matrix assembly function as above. Note, we do not
have to attach an init() function since we are initializing the
system by reading in "saved_solution.xda"
</div>

<div class ="fragment">
<pre>
              system.attach_assemble_function (assemble_cd);
        
</pre>
</div>
<div class = "comment">
Print out the H1 norm of the saved solution, for verification purposes:
</div>

<div class ="fragment">
<pre>
              Real H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));
        
              std::cout &lt;&lt; "Initial H1 norm = " &lt;&lt; H1norm &lt;&lt; std::endl &lt;&lt; std::endl;
            }
        
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
            ExodusII_IO(mesh).write_equation_systems ("out.e.000",
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
          TransientLinearImplicitSystem& system =
        	equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;
                  ("Convection-Diffusion");
        
          const Real dt = 0.025;
          system.time   = init_timestep*dt;
        
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
              
</pre>
</div>
<div class = "comment">
As already seen in example 9, use a work-around
for missing stream functionality (of older compilers).
Add a set of scope braces to enforce data locality.
</div>

<div class ="fragment">
<pre>
              {
                OStringStream out;
        
                OSSInt(out,2,t_step);
                out &lt;&lt; ", time=";
                OSSRealzeroleft(out,6,3,system.time);
                out &lt;&lt;  "...";
                std::cout &lt;&lt; out.str() &lt;&lt; std::endl;
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
              const unsigned int max_r_steps = 2;
              
</pre>
</div>
<div class = "comment">
A refinement loop.
</div>

<div class ="fragment">
<pre>
              for (unsigned int r_step=0; r_step&lt;max_r_steps; r_step++)
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
                  if (r_step+1 != max_r_steps)
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
                  OStringStream file_name, exodus_file_name;
        
        #ifdef LIBMESH_HAVE_GMV
                  file_name &lt;&lt; "out.gmv.";
                  OSSRealzeroright(file_name,3,0,t_step+1);
        
                  GMVIO(mesh).write_equation_systems (file_name.str(),
                                                      equation_systems);
        #endif
        #ifdef LIBMESH_HAVE_EXODUS_API
                  exodus_file_name &lt;&lt; "out.e.";
                  OSSRealzeroright(exodus_file_name,3,0,t_step+1);
                  ExodusII_IO(mesh).write_equation_systems (exodus_file_name.str(),
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
          libmesh_assert (system_name == "Convection-Diffusion");
        
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
          libmesh_assert (system_name == "Convection-Diffusion");
          
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
          const std::vector&lt;Real&gt;& JxW_face = fe_face-&gt;get_JxW();
          
</pre>
</div>
<div class = "comment">
The element shape functions evaluated at the quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi = fe-&gt;get_phi();
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& psi = fe_face-&gt;get_phi();
        
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
The XY locations of the quadrature points used for face integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Point&gt;& qface_points = fe_face-&gt;get_xyz();
            
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
          std::vector&lt;unsigned int&gt; dof_indices;
        
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
          const Real time = es.parameters.get&lt;Real&gt; ("time");
        
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
At this point the interior element integration has
been completed.  However, we have not yet addressed
boundary conditions.  For this example we will only
consider simple Dirichlet boundary conditions imposed
via the penalty method. 

<br><br>The following loops over the sides of the element.
If the element has no neighbor on a side then that
side MUST live on a boundary of the domain.
</div>

<div class ="fragment">
<pre>
              {
</pre>
</div>
<div class = "comment">
The penalty value.  
</div>

<div class ="fragment">
<pre>
                const Real penalty = 1.e10;
        
</pre>
</div>
<div class = "comment">
The following loops over the sides of the element.
If the element has no neighbor on a side then that
side MUST live on a boundary of the domain.
</div>

<div class ="fragment">
<pre>
                for (unsigned int s=0; s&lt;elem-&gt;n_sides(); s++)
                  if (elem-&gt;neighbor(s) == NULL)
                    {
                      fe_face-&gt;reinit(elem,s);
                      
                      for (unsigned int qp=0; qp&lt;qface.n_points(); qp++)
                        {
                          const Number value = 
                            parsed_solution ?
                              (*parsed_solution)(qface_points[qp], time) :
                              exact_solution (qface_points[qp](0),
                                              qface_points[qp](1), time);
                                                               
</pre>
</div>
<div class = "comment">
RHS contribution
</div>

<div class ="fragment">
<pre>
                          for (unsigned int i=0; i&lt;psi.size(); i++)
                            Fe(i) += penalty*JxW_face[qp]*value*psi[i][qp];
        
</pre>
</div>
<div class = "comment">
Matrix contribution
</div>

<div class ="fragment">
<pre>
                          for (unsigned int i=0; i&lt;psi.size(); i++)
                            for (unsigned int j=0; j&lt;psi.size(); j++)
                              Ke(i,j) += penalty*JxW_face[qp]*psi[i][qp]*psi[j][qp];
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
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
   
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;cmath&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;gmv_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;periodic_boundaries.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;parsed_function.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;o_string_stream.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;transient_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;vector_value.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;kelly_error_estimator.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  
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
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_ENABLE_AMR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-amr&quot;</FONT></B>);
  #elif !defined(LIBMESH_HAVE_XDR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-xdr&quot;</FONT></B>);
  #elif !defined(LIBMESH_ENABLE_PERIODIC)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-periodic&quot;</FONT></B>);
  #elif defined(LIBMESH_ENABLE_PARMESH)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-parmesh&quot;</FONT></B>);
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
  
    Mesh mesh;
  
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
  
        system.attach_assemble_function (assemble_cd);
        system.attach_init_function (init_cd);
  
        equation_systems.init ();
      }
    <B><FONT COLOR="#A020F0">else</FONT></B> 
      {
        mesh.read(<B><FONT COLOR="#BC8F8F">&quot;saved_mesh.xdr&quot;</FONT></B>);
  
        mesh.print_info();
  
        equation_systems.read(<B><FONT COLOR="#BC8F8F">&quot;saved_solution.xdr&quot;</FONT></B>, libMeshEnums::DECODE);
  
        TransientLinearImplicitSystem &amp; system = 
          equation_systems.get_system&lt;TransientLinearImplicitSystem&gt; 
          (<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
        system.update();
  
        system.attach_assemble_function (assemble_cd);
  
        Real H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Initial H1 norm = &quot;</FONT></B> &lt;&lt; H1norm &lt;&lt; std::endl &lt;&lt; std::endl;
      }
  
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
      ExodusII_IO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;out.e.000&quot;</FONT></B>,
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
      
  
    TransientLinearImplicitSystem&amp; system =
  	equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;
            (<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real dt = 0.025;
    system.time   = init_timestep*dt;
  
    DofMap &amp; dof_map = system.get_dof_map();
  
    
    PeriodicBoundary horz(RealVectorValue(2.0, 0., 0.));
  
    horz.myboundary = 3;
    horz.pairedboundary = 1;
  
    dof_map.add_periodic_boundary(horz);
  
    
    PeriodicBoundary vert(RealVectorValue(0., 2.0, 0.));
  
    vert.myboundary = 0;
    vert.pairedboundary = 2;
  
    dof_map.add_periodic_boundary(vert);
    
    
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
          OStringStream out;
  
          OSSInt(out,2,t_step);
          out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, time=&quot;</FONT></B>;
          OSSRealzeroleft(out,6,3,system.time);
          out &lt;&lt;  <B><FONT COLOR="#BC8F8F">&quot;...&quot;</FONT></B>;
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; out.str() &lt;&lt; std::endl;
        }
        
        TransientLinearImplicitSystem &amp;  system =
          equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
        *system.old_local_solution = *system.current_local_solution;
        
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_r_steps = 2;
        
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> r_step=0; r_step&lt;max_r_steps; r_step++)
          {
            system.solve();
  
            Real H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;H1 norm = &quot;</FONT></B> &lt;&lt; H1norm &lt;&lt; std::endl;
            
            <B><FONT COLOR="#A020F0">if</FONT></B> (r_step+1 != max_r_steps)
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
            OStringStream file_name, exodus_file_name;
  
  #ifdef LIBMESH_HAVE_GMV
            file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out.gmv.&quot;</FONT></B>;
            OSSRealzeroright(file_name,3,0,t_step+1);
  
            GMVIO(mesh).write_equation_systems (file_name.str(),
                                                equation_systems);
  #endif
  #ifdef LIBMESH_HAVE_EXODUS_API
            exodus_file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out.e.&quot;</FONT></B>;
            OSSRealzeroright(exodus_file_name,3,0,t_step+1);
            ExodusII_IO(mesh).write_equation_systems (exodus_file_name.str(),
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
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
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
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
    
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
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
    
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; psi = fe_face-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; qface_points = fe_face-&gt;get_xyz();
      
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
    
    DenseMatrix&lt;Number&gt; Ke;
    DenseVector&lt;Number&gt; Fe;
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
  
    <B><FONT COLOR="#228B22">const</FONT></B> RealVectorValue velocity =
      es.parameters.get&lt;RealVectorValue&gt; (<B><FONT COLOR="#BC8F8F">&quot;velocity&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real diffusivity =
      es.parameters.get&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;diffusivity&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real dt = es.parameters.get&lt;Real&gt;   (<B><FONT COLOR="#BC8F8F">&quot;dt&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> Real time = es.parameters.get&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>);
  
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
  
        {
          <B><FONT COLOR="#228B22">const</FONT></B> Real penalty = 1.e10;
  
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> s=0; s&lt;elem-&gt;n_sides(); s++)
            <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(s) == NULL)
              {
                fe_face-&gt;reinit(elem,s);
                
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qface.n_points(); qp++)
                  {
                    <B><FONT COLOR="#228B22">const</FONT></B> Number value = 
                      parsed_solution ?
                        (*parsed_solution)(qface_points[qp], time) :
                        exact_solution (qface_points[qp](0),
                                        qface_points[qp](1), time);
                                                         
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;psi.size(); i++)
                      Fe(i) += penalty*JxW_face[qp]*value*psi[i][qp];
  
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;psi.size(); i++)
                      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;psi.size(); j++)
                        Ke(i,j) += penalty*JxW_face[qp]*psi[i][qp]*psi[j][qp];
                  }
              } 
        } 
  
        
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        
        system.matrix-&gt;add_matrix (Ke, dof_indices);
        system.rhs-&gt;add_vector    (Fe, dof_indices);
        
      }
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_ENABLE_AMR
</FONT></I>  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Compiling C++ (in optimized mode) adaptivity_ex5.C...
Linking adaptivity_ex5-opt...
***************************************************************
* Running Example  ./adaptivity_ex5-opt [-read_solution] -n_timesteps 25 -n_refinements 5 -init_timestep [0|25]
***************************************************************
 
Usage:
	 ./adaptivity_ex5-opt -init_timestep 0
OR
	 ./adaptivity_ex5-opt -read_solution -init_timestep 26

Running: ./adaptivity_ex5-opt -n_timesteps 25 -n_refinements 5 -output_freq 10 -init_timestep 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=4225
    n_local_nodes()=4225
  n_elem()=5460
    n_local_elem()=5460
    n_active_elem()=4096
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=4225
    n_local_dofs()=4225
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.81633
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

 Solving time step  0, time=0.0250...
H1 norm = 1.58843
  Refining the mesh...
H1 norm = 1.59127
 Solving time step  1, time=0.0500...
H1 norm = 1.46336
  Refining the mesh...
H1 norm = 1.46259
 Solving time step  2, time=0.0750...
H1 norm = 1.35363
  Refining the mesh...
H1 norm = 1.35325
 Solving time step  3, time=0.1000...
H1 norm = 1.25937
  Refining the mesh...
H1 norm = 1.25874
 Solving time step  4, time=0.1250...
H1 norm = 1.17677
  Refining the mesh...
H1 norm = 1.17658
 Solving time step  5, time=0.1500...
H1 norm = 1.1046
  Refining the mesh...
H1 norm = 1.10427
 Solving time step  6, time=0.1750...
H1 norm = 1.04052
  Refining the mesh...
H1 norm = 1.04019
 Solving time step  7, time=0.2000...
H1 norm = 0.983359
  Refining the mesh...
H1 norm = 0.983071
 Solving time step  8, time=0.2250...
H1 norm = 0.932118
  Refining the mesh...
H1 norm = 0.931939
 Solving time step  9, time=0.2500...
H1 norm = 0.885955
  Refining the mesh...
H1 norm = 0.885809
 Solving time step 10, time=0.2750...
H1 norm = 0.844184
  Refining the mesh...
H1 norm = 0.844012
 Solving time step 11, time=0.3000...
H1 norm = 0.806099
  Refining the mesh...
H1 norm = 0.805971
 Solving time step 12, time=0.3250...
H1 norm = 0.771338
  Refining the mesh...
H1 norm = 0.771246
 Solving time step 13, time=0.3500...
H1 norm = 0.739486
  Refining the mesh...
H1 norm = 0.739305
 Solving time step 14, time=0.3750...
H1 norm = 0.710032
  Refining the mesh...
H1 norm = 0.709957
 Solving time step 15, time=0.4000...
H1 norm = 0.682936
  Refining the mesh...
H1 norm = 0.68284
 Solving time step 16, time=0.4250...
H1 norm = 0.657808
  Refining the mesh...
H1 norm = 0.657768
 Solving time step 17, time=0.4500...
H1 norm = 0.634539
  Refining the mesh...
H1 norm = 0.634445
 Solving time step 18, time=0.4750...
H1 norm = 0.612805
  Refining the mesh...
H1 norm = 0.612739
 Solving time step 19, time=0.5000...
H1 norm = 0.592533
  Refining the mesh...
H1 norm = 0.592455
 Solving time step 20, time=0.5250...
H1 norm = 0.573539
  Refining the mesh...
H1 norm = 0.573521
 Solving time step 21, time=0.5500...
H1 norm = 0.555804
  Refining the mesh...
H1 norm = 0.555726
 Solving time step 22, time=0.5750...
H1 norm = 0.539074
  Refining the mesh...
H1 norm = 0.539004
 Solving time step 23, time=0.6000...
H1 norm = 0.523329
  Refining the mesh...
H1 norm = 0.523287
 Solving time step 24, time=0.6250...
H1 norm = 0.508503
  Refining the mesh...
H1 norm = 0.508479
Final H1 norm = 0.508479


-------------------------------------------------------------------
| Time:           Sat Apr  7 15:55:55 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=3.0778, Active time=2.7804                                                  |
 -------------------------------------------------------------------------------------------------------------
| Event                           nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                           w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------|
|                                                                                                             |
|                                                                                                             |
| DofMap                                                                                                      |
|   add_neighbors_to_send_list()  26        0.0075      0.000289    0.0075      0.000289    0.27     0.27     |
|   build_constraint_matrix()     25348     0.0372      0.000001    0.0372      0.000001    1.34     1.34     |
|   cnstrn_elem_mat_vec()         25348     0.1031      0.000004    0.1031      0.000004    3.71     3.71     |
|   compute_sparsity()            26        0.0649      0.002495    0.0794      0.003053    2.33     2.86     |
|   create_dof_constraints()      26        0.1288      0.004952    0.2337      0.008988    4.63     8.40     |
|   distribute_dofs()             26        0.0121      0.000465    0.0481      0.001852    0.44     1.73     |
|   dof_indices()                 208010    0.1598      0.000001    0.1598      0.000001    5.75     5.75     |
|   enforce_constraints_exactly() 75        0.0047      0.000063    0.0047      0.000063    0.17     0.17     |
|   old_dof_indices()             90411     0.0548      0.000001    0.0548      0.000001    1.97     1.97     |
|   prepare_send_list()           26        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                      26        0.0360      0.001385    0.0360      0.001385    1.30     1.30     |
|                                                                                                             |
| EquationSystems                                                                                             |
|   build_solution_vector()       8         0.0155      0.001934    0.0269      0.003368    0.56     0.97     |
|   write()                       1         0.0113      0.011319    0.0117      0.011670    0.41     0.42     |
|                                                                                                             |
| ExodusII_IO                                                                                                 |
|   write_nodal_data()            4         0.0106      0.002638    0.0106      0.002638    0.38     0.38     |
|                                                                                                             |
| FE                                                                                                          |
|   compute_affine_map()          130524    0.0991      0.000001    0.0991      0.000001    3.56     3.56     |
|   compute_face_map()            36868     0.1100      0.000003    0.2204      0.000006    3.96     7.93     |
|   compute_shape_functions()     130524    0.0615      0.000000    0.0615      0.000000    2.21     2.21     |
|   init_face_shape_functions()   997       0.0016      0.000002    0.0016      0.000002    0.06     0.06     |
|   init_shape_functions()        71067     0.1776      0.000002    0.1776      0.000002    6.39     6.39     |
|   inverse_map()                 168777    0.2165      0.000001    0.2165      0.000001    7.79     7.79     |
|                                                                                                             |
| GMVIO                                                                                                       |
|   write_nodal_data()            4         0.0393      0.009821    0.0393      0.009821    1.41     1.41     |
|                                                                                                             |
| JumpErrorEstimator                                                                                          |
|   estimate_error()              25        0.4646      0.018583    1.1043      0.044171    16.71    39.72    |
|                                                                                                             |
| LocationMap                                                                                                 |
|   find()                        19140     0.0143      0.000001    0.0143      0.000001    0.52     0.52     |
|   init()                        55        0.0102      0.000185    0.0102      0.000185    0.37     0.37     |
|                                                                                                             |
| Mesh                                                                                                        |
|   contract()                    25        0.0053      0.000214    0.0092      0.000370    0.19     0.33     |
|   find_neighbors()              27        0.0756      0.002802    0.0756      0.002802    2.72     2.72     |
|   renumber_nodes_and_elem()     79        0.0089      0.000113    0.0089      0.000113    0.32     0.32     |
|                                                                                                             |
| MeshCommunication                                                                                           |
|   assign_global_indices()       1         0.0155      0.015498    0.0155      0.015511    0.56     0.56     |
|                                                                                                             |
| MeshOutput                                                                                                  |
|   write_equation_systems()      8         0.0001      0.000016    0.0769      0.009615    0.00     2.77     |
|                                                                                                             |
| MeshRefinement                                                                                              |
|   _coarsen_elements()           50        0.0041      0.000081    0.0041      0.000081    0.15     0.15     |
|   _refine_elements()            55        0.0383      0.000697    0.0832      0.001512    1.38     2.99     |
|   add_point()                   19140     0.0231      0.000001    0.0411      0.000002    0.83     1.48     |
|   make_coarsening_compatible()  61        0.0613      0.001005    0.0864      0.001416    2.20     3.11     |
|   make_refinement_compatible()  61        0.0018      0.000030    0.0019      0.000032    0.07     0.07     |
|                                                                                                             |
| MeshTools::Generation                                                                                       |
|   build_cube()                  1         0.0002      0.000169    0.0002      0.000169    0.01     0.01     |
|                                                                                                             |
| Parallel                                                                                                    |
|   allgather()                   30        0.0000      0.000000    0.0000      0.000000    0.00     0.00     |
|   receive()                     4         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   send()                        4         0.0003      0.000081    0.0003      0.000081    0.01     0.01     |
|   send_receive()                4         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|                                                                                                             |
| Partitioner                                                                                                 |
|   single_partition()            27        0.0020      0.000076    0.0020      0.000076    0.07     0.07     |
|                                                                                                             |
| PetscLinearSolver                                                                                           |
|   solve()                       50        0.0490      0.000980    0.0490      0.000980    1.76     1.76     |
|                                                                                                             |
| PointLocatorTree                                                                                            |
|   init(no master)               50        0.0689      0.001378    0.0689      0.001378    2.48     2.48     |
|   operator()                    4949      0.0367      0.000007    0.0451      0.000009    1.32     1.62     |
|                                                                                                             |
| ProjectVector                                                                                               |
|   operator()                    75        0.0847      0.001129    0.1481      0.001975    3.05     5.33     |
|                                                                                                             |
| System                                                                                                      |
|   assemble()                    50        0.2687      0.005374    0.5270      0.010541    9.66     18.96    |
|   calculate_norm()              51        0.1154      0.002262    0.1908      0.003742    4.15     6.86     |
|   project_vector()              76        0.0762      0.001002    0.2765      0.003638    2.74     9.95     |
|                                                                                                             |
| XdrIO                                                                                                       |
|   write()                       1         0.0032      0.003218    0.0032      0.003218    0.12     0.12     |
 -------------------------------------------------------------------------------------------------------------
| Totals:                         932221    2.7804                                          100.00            |
 -------------------------------------------------------------------------------------------------------------

 
***** Finished first 25 steps, now read in saved solution and continue *****
 
Usage:
	 ./adaptivity_ex5-opt -init_timestep 0
OR
	 ./adaptivity_ex5-opt -read_solution -init_timestep 26

Running: ./adaptivity_ex5-opt -read_solution -n_timesteps 25 -output_freq 10 -init_timestep 25

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=743
    n_local_nodes()=743
  n_elem()=892
    n_local_elem()=892
    n_active_elem()=670
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0

Initial H1 norm = 0.508479

 EquationSystems
  n_systems()=1
   System #0, "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=743
    n_local_dofs()=743
    n_constrained_dofs()=110
    n_local_constrained_dofs()=110
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 9.41454
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 18
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 110
      Average DoF Constraint Length= 2
      Number of Node Constraints = 217
      Maximum Node Constraint Length= 5
      Average Node Constraint Length= 2.52074

 Solving time step 25, time=0.6500...
H1 norm = 0.494524
  Refining the mesh...
H1 norm = 0.494485
 Solving time step 26, time=0.6750...
H1 norm = 0.481283
  Refining the mesh...
H1 norm = 0.481226
 Solving time step 27, time=0.7000...
H1 norm = 0.468713
  Refining the mesh...
H1 norm = 0.468697
 Solving time step 28, time=0.7250...
H1 norm = 0.456834
  Refining the mesh...
H1 norm = 0.45679
 Solving time step 29, time=0.7500...
H1 norm = 0.445515
  Refining the mesh...
H1 norm = 0.445494
 Solving time step 30, time=0.7750...
H1 norm = 0.434771
  Refining the mesh...
H1 norm = 0.434742
 Solving time step 31, time=0.8000...
H1 norm = 0.424533
  Refining the mesh...
H1 norm = 0.424515
 Solving time step 32, time=0.8250...
H1 norm = 0.414781
  Refining the mesh...
H1 norm = 0.414758
 Solving time step 33, time=0.8500...
H1 norm = 0.40547
  Refining the mesh...
H1 norm = 0.405446
 Solving time step 34, time=0.8750...
H1 norm = 0.396568
  Refining the mesh...
H1 norm = 0.396544
 Solving time step 35, time=0.9000...
H1 norm = 0.388053
  Refining the mesh...
H1 norm = 0.388038
 Solving time step 36, time=0.9250...
H1 norm = 0.37991
  Refining the mesh...
H1 norm = 0.379888
 Solving time step 37, time=0.9500...
H1 norm = 0.372095
  Refining the mesh...
H1 norm = 0.372082
 Solving time step 38, time=0.9750...
H1 norm = 0.364611
  Refining the mesh...
H1 norm = 0.364592
 Solving time step 39, time=1.0000...
H1 norm = 0.357419
  Refining the mesh...
H1 norm = 0.357396
 Solving time step 40, time=1.0300...
H1 norm = 0.350503
  Refining the mesh...
H1 norm = 0.350499
 Solving time step 41, time=1.0500...
H1 norm = 0.343875
  Refining the mesh...
H1 norm = 0.343856
 Solving time step 42, time=1.0700...
H1 norm = 0.337479
  Refining the mesh...
H1 norm = 0.337473
 Solving time step 43, time=1.1000...
H1 norm = 0.331335
  Refining the mesh...
H1 norm = 0.331322
 Solving time step 44, time=1.1200...
H1 norm = 0.32541
  Refining the mesh...
H1 norm = 0.325393
 Solving time step 45, time=1.1500...
H1 norm = 0.319688
  Refining the mesh...
H1 norm = 0.31968
 Solving time step 46, time=1.1700...
H1 norm = 0.314178
  Refining the mesh...
H1 norm = 0.314172
 Solving time step 47, time=1.2000...
H1 norm = 0.30886
  Refining the mesh...
H1 norm = 0.30885
 Solving time step 48, time=1.2200...
H1 norm = 0.303718
  Refining the mesh...
H1 norm = 0.303705
 Solving time step 49, time=1.2500...
H1 norm = 0.298746
  Refining the mesh...
H1 norm = 0.298734

-------------------------------------------------------------------
| Time:           Sat Apr  7 15:55:59 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=3.78331, Active time=3.45785                                                |
 -------------------------------------------------------------------------------------------------------------
| Event                           nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                           w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------|
|                                                                                                             |
|                                                                                                             |
| DofMap                                                                                                      |
|   add_neighbors_to_send_list()  26        0.0084      0.000325    0.0084      0.000325    0.24     0.24     |
|   build_constraint_matrix()     43136     0.0533      0.000001    0.0533      0.000001    1.54     1.54     |
|   cnstrn_elem_mat_vec()         43136     0.1281      0.000003    0.1281      0.000003    3.71     3.71     |
|   compute_sparsity()            26        0.0832      0.003200    0.1000      0.003848    2.41     2.89     |
|   create_dof_constraints()      26        0.1628      0.006263    0.2705      0.010404    4.71     7.82     |
|   distribute_dofs()             26        0.0115      0.000442    0.0510      0.001962    0.33     1.47     |
|   dof_indices()                 283516    0.2092      0.000001    0.2092      0.000001    6.05     6.05     |
|   enforce_constraints_exactly() 75        0.0058      0.000077    0.0058      0.000077    0.17     0.17     |
|   old_dof_indices()             132162    0.0823      0.000001    0.0823      0.000001    2.38     2.38     |
|   prepare_send_list()           26        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                      26        0.0395      0.001518    0.0395      0.001518    1.14     1.14     |
|                                                                                                             |
| EquationSystems                                                                                             |
|   build_solution_vector()       8         0.0090      0.001119    0.0152      0.001898    0.26     0.44     |
|   read()                        1         0.0125      0.012504    0.0365      0.036539    0.36     1.06     |
|   update()                      1         0.0001      0.000060    0.0001      0.000060    0.00     0.00     |
|                                                                                                             |
| ExodusII_IO                                                                                                 |
|   write_nodal_data()            4         0.0066      0.001646    0.0066      0.001646    0.19     0.19     |
|                                                                                                             |
| FE                                                                                                          |
|   compute_affine_map()          177698    0.1361      0.000001    0.1361      0.000001    3.94     3.94     |
|   compute_face_map()            46198     0.1398      0.000003    0.2805      0.000006    4.04     8.11     |
|   compute_shape_functions()     177698    0.0847      0.000000    0.0847      0.000000    2.45     2.45     |
|   init_face_shape_functions()   691       0.0013      0.000002    0.0013      0.000002    0.04     0.04     |
|   init_shape_functions()        90857     0.2307      0.000003    0.2307      0.000003    6.67     6.67     |
|   inverse_map()                 214824    0.2769      0.000001    0.2769      0.000001    8.01     8.01     |
|                                                                                                             |
| GMVIO                                                                                                       |
|   write_nodal_data()            4         0.0211      0.005286    0.0211      0.005286    0.61     0.61     |
|                                                                                                             |
| JumpErrorEstimator                                                                                          |
|   estimate_error()              25        0.6260      0.025039    1.4860      0.059441    18.10    42.98    |
|                                                                                                             |
| LocationMap                                                                                                 |
|   find()                        3672      0.0025      0.000001    0.0025      0.000001    0.07     0.07     |
|   init()                        50        0.0130      0.000259    0.0130      0.000259    0.37     0.37     |
|                                                                                                             |
| Mesh                                                                                                        |
|   contract()                    25        0.0023      0.000091    0.0050      0.000201    0.07     0.15     |
|   find_neighbors()              26        0.0755      0.002905    0.0755      0.002905    2.18     2.18     |
|   renumber_nodes_and_elem()     77        0.0067      0.000087    0.0067      0.000087    0.19     0.19     |
|                                                                                                             |
| MeshCommunication                                                                                           |
|   assign_global_indices()       1         0.0156      0.015609    0.0156      0.015616    0.45     0.45     |
|                                                                                                             |
| MeshOutput                                                                                                  |
|   write_equation_systems()      8         0.0001      0.000016    0.0430      0.005381    0.00     1.24     |
|                                                                                                             |
| MeshRefinement                                                                                              |
|   _coarsen_elements()           50        0.0023      0.000046    0.0023      0.000046    0.07     0.07     |
|   _refine_elements()            50        0.0134      0.000269    0.0218      0.000436    0.39     0.63     |
|   add_point()                   3672      0.0047      0.000001    0.0078      0.000002    0.14     0.22     |
|   make_coarsening_compatible()  58        0.0709      0.001222    0.0774      0.001334    2.05     2.24     |
|   make_refinement_compatible()  58        0.0022      0.000038    0.0024      0.000042    0.06     0.07     |
|                                                                                                             |
| Parallel                                                                                                    |
|   allgather()                   32        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   send_receive()                4         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                             |
| Partitioner                                                                                                 |
|   set_node_processor_ids()      1         0.0002      0.000243    0.0002      0.000243    0.01     0.01     |
|   single_partition()            26        0.0021      0.000079    0.0021      0.000079    0.06     0.06     |
|                                                                                                             |
| PetscLinearSolver                                                                                           |
|   solve()                       50        0.0605      0.001211    0.0605      0.001211    1.75     1.75     |
|                                                                                                             |
| PointLocatorTree                                                                                            |
|   init(no master)               50        0.0871      0.001742    0.0871      0.001742    2.52     2.52     |
|   operator()                    2808      0.0129      0.000005    0.0179      0.000006    0.37     0.52     |
|                                                                                                             |
| ProjectVector                                                                                               |
|   operator()                    75        0.1361      0.001815    0.2425      0.003234    3.94     7.01     |
|                                                                                                             |
| System                                                                                                      |
|   assemble()                    50        0.3736      0.007472    0.7027      0.014053    10.80    20.32    |
|   calculate_norm()              51        0.1672      0.003278    0.2775      0.005442    4.83     8.03     |
|   project_vector()              75        0.0785      0.001047    0.3850      0.005133    2.27     11.13    |
|                                                                                                             |
| XdrIO                                                                                                       |
|   read()                        1         0.0015      0.001550    0.0015      0.001550    0.04     0.04     |
 -------------------------------------------------------------------------------------------------------------
| Totals:                         1221160   3.4578                                          100.00            |
 -------------------------------------------------------------------------------------------------------------

 
Usage:
	 ./adaptivity_ex5-opt -init_timestep 0
OR
	 ./adaptivity_ex5-opt -read_solution -init_timestep 26

Running: ./adaptivity_ex5-opt -n_timesteps 25 -n_refinements 5 -output_freq 10 -init_timestep 0 -exact_solution 10*exp(-(pow(x-0.8*t-0.2,2) + pow(y-0.8*t-0.2,2))/(0.01*(4*t+1)))/(4*t+1)

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=4225
    n_local_nodes()=4225
  n_elem()=5460
    n_local_elem()=5460
    n_active_elem()=4096
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=4225
    n_local_dofs()=4225
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.81633
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

 Solving time step  0, time=0.0250...
H1 norm = 15.8843
  Refining the mesh...
H1 norm = 15.9127
 Solving time step  1, time=0.0500...
H1 norm = 14.6336
  Refining the mesh...
H1 norm = 14.6259
 Solving time step  2, time=0.0750...
H1 norm = 13.5363
  Refining the mesh...
H1 norm = 13.5325
 Solving time step  3, time=0.1000...
H1 norm = 12.5937
  Refining the mesh...
H1 norm = 12.5874
 Solving time step  4, time=0.1250...
H1 norm = 11.7677
  Refining the mesh...
H1 norm = 11.7658
 Solving time step  5, time=0.1500...
H1 norm = 11.046
  Refining the mesh...
H1 norm = 11.0427
 Solving time step  6, time=0.1750...
H1 norm = 10.4052
  Refining the mesh...
H1 norm = 10.4019
 Solving time step  7, time=0.2000...
H1 norm = 9.83359
  Refining the mesh...
H1 norm = 9.83071
 Solving time step  8, time=0.2250...
H1 norm = 9.32118
  Refining the mesh...
H1 norm = 9.31939
 Solving time step  9, time=0.2500...
H1 norm = 8.85955
  Refining the mesh...
H1 norm = 8.85809
 Solving time step 10, time=0.2750...
H1 norm = 8.44184
  Refining the mesh...
H1 norm = 8.44012
 Solving time step 11, time=0.3000...
H1 norm = 8.06099
  Refining the mesh...
H1 norm = 8.05971
 Solving time step 12, time=0.3250...
H1 norm = 7.71338
  Refining the mesh...
H1 norm = 7.71246
 Solving time step 13, time=0.3500...
H1 norm = 7.39486
  Refining the mesh...
H1 norm = 7.39305
 Solving time step 14, time=0.3750...
H1 norm = 7.10032
  Refining the mesh...
H1 norm = 7.09957
 Solving time step 15, time=0.4000...
H1 norm = 6.82936
  Refining the mesh...
H1 norm = 6.8284
 Solving time step 16, time=0.4250...
H1 norm = 6.57808
  Refining the mesh...
H1 norm = 6.57768
 Solving time step 17, time=0.4500...
H1 norm = 6.34539
  Refining the mesh...
H1 norm = 6.34445
 Solving time step 18, time=0.4750...
H1 norm = 6.12805
  Refining the mesh...
H1 norm = 6.12739
 Solving time step 19, time=0.5000...
H1 norm = 5.92533
  Refining the mesh...
H1 norm = 5.92455
 Solving time step 20, time=0.5250...
H1 norm = 5.73539
  Refining the mesh...
H1 norm = 5.73521
 Solving time step 21, time=0.5500...
H1 norm = 5.55804
  Refining the mesh...
H1 norm = 5.55726
 Solving time step 22, time=0.5750...
H1 norm = 5.39074
  Refining the mesh...
H1 norm = 5.39004
 Solving time step 23, time=0.6000...
H1 norm = 5.23329
  Refining the mesh...
H1 norm = 5.23287
 Solving time step 24, time=0.6250...
H1 norm = 5.08503
  Refining the mesh...
H1 norm = 5.08479
Final H1 norm = 5.08479


-------------------------------------------------------------------
| Time:           Sat Apr  7 15:56:02 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=3.03799, Active time=2.73411                                                |
 -------------------------------------------------------------------------------------------------------------
| Event                           nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                           w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------|
|                                                                                                             |
|                                                                                                             |
| DofMap                                                                                                      |
|   add_neighbors_to_send_list()  26        0.0074      0.000283    0.0074      0.000283    0.27     0.27     |
|   build_constraint_matrix()     25348     0.0352      0.000001    0.0352      0.000001    1.29     1.29     |
|   cnstrn_elem_mat_vec()         25348     0.0848      0.000003    0.0848      0.000003    3.10     3.10     |
|   compute_sparsity()            26        0.0574      0.002207    0.0717      0.002757    2.10     2.62     |
|   create_dof_constraints()      26        0.1221      0.004695    0.2255      0.008671    4.46     8.25     |
|   distribute_dofs()             26        0.0121      0.000465    0.0471      0.001810    0.44     1.72     |
|   dof_indices()                 208010    0.1574      0.000001    0.1574      0.000001    5.76     5.76     |
|   enforce_constraints_exactly() 75        0.0046      0.000061    0.0046      0.000061    0.17     0.17     |
|   old_dof_indices()             90411     0.0553      0.000001    0.0553      0.000001    2.02     2.02     |
|   prepare_send_list()           26        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                      26        0.0350      0.001345    0.0350      0.001345    1.28     1.28     |
|                                                                                                             |
| EquationSystems                                                                                             |
|   build_solution_vector()       8         0.0152      0.001901    0.0298      0.003731    0.56     1.09     |
|   write()                       1         0.0097      0.009687    0.0100      0.009977    0.35     0.36     |
|                                                                                                             |
| ExodusII_IO                                                                                                 |
|   write_nodal_data()            4         0.0099      0.002464    0.0099      0.002464    0.36     0.36     |
|                                                                                                             |
| FE                                                                                                          |
|   compute_affine_map()          130524    0.1004      0.000001    0.1004      0.000001    3.67     3.67     |
|   compute_face_map()            36868     0.1120      0.000003    0.2233      0.000006    4.10     8.17     |
|   compute_shape_functions()     130524    0.0612      0.000000    0.0612      0.000000    2.24     2.24     |
|   init_face_shape_functions()   997       0.0016      0.000002    0.0016      0.000002    0.06     0.06     |
|   init_shape_functions()        71067     0.1810      0.000003    0.1810      0.000003    6.62     6.62     |
|   inverse_map()                 168777    0.2228      0.000001    0.2228      0.000001    8.15     8.15     |
|                                                                                                             |
| GMVIO                                                                                                       |
|   write_nodal_data()            4         0.0383      0.009585    0.0383      0.009585    1.40     1.40     |
|                                                                                                             |
| JumpErrorEstimator                                                                                          |
|   estimate_error()              25        0.4675      0.018699    1.1180      0.044719    17.10    40.89    |
|                                                                                                             |
| LocationMap                                                                                                 |
|   find()                        19140     0.0157      0.000001    0.0157      0.000001    0.57     0.57     |
|   init()                        55        0.0100      0.000181    0.0100      0.000181    0.36     0.36     |
|                                                                                                             |
| Mesh                                                                                                        |
|   contract()                    25        0.0042      0.000167    0.0072      0.000288    0.15     0.26     |
|   find_neighbors()              27        0.0741      0.002744    0.0741      0.002744    2.71     2.71     |
|   renumber_nodes_and_elem()     79        0.0073      0.000093    0.0073      0.000093    0.27     0.27     |
|                                                                                                             |
| MeshCommunication                                                                                           |
|   assign_global_indices()       1         0.0153      0.015350    0.0154      0.015362    0.56     0.56     |
|                                                                                                             |
| MeshOutput                                                                                                  |
|   write_equation_systems()      8         0.0001      0.000014    0.0782      0.009770    0.00     2.86     |
|                                                                                                             |
| MeshRefinement                                                                                              |
|   _coarsen_elements()           50        0.0036      0.000072    0.0036      0.000072    0.13     0.13     |
|   _refine_elements()            55        0.0379      0.000689    0.0853      0.001551    1.39     3.12     |
|   add_point()                   19140     0.0236      0.000001    0.0433      0.000002    0.86     1.58     |
|   make_coarsening_compatible()  61        0.0602      0.000986    0.0850      0.001393    2.20     3.11     |
|   make_refinement_compatible()  61        0.0018      0.000029    0.0019      0.000031    0.06     0.07     |
|                                                                                                             |
| MeshTools::Generation                                                                                       |
|   build_cube()                  1         0.0001      0.000138    0.0001      0.000138    0.01     0.01     |
|                                                                                                             |
| Parallel                                                                                                    |
|   allgather()                   30        0.0000      0.000000    0.0000      0.000000    0.00     0.00     |
|   receive()                     4         0.0000      0.000005    0.0000      0.000005    0.00     0.00     |
|   send()                        4         0.0003      0.000067    0.0003      0.000067    0.01     0.01     |
|   send_receive()                4         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|                                                                                                             |
| Partitioner                                                                                                 |
|   single_partition()            27        0.0020      0.000074    0.0020      0.000074    0.07     0.07     |
|                                                                                                             |
| PetscLinearSolver                                                                                           |
|   solve()                       50        0.0473      0.000946    0.0473      0.000946    1.73     1.73     |
|                                                                                                             |
| PointLocatorTree                                                                                            |
|   init(no master)               50        0.0653      0.001306    0.0653      0.001306    2.39     2.39     |
|   operator()                    4949      0.0364      0.000007    0.0445      0.000009    1.33     1.63     |
|                                                                                                             |
| ProjectVector                                                                                               |
|   operator()                    75        0.0845      0.001126    0.1467      0.001956    3.09     5.36     |
|                                                                                                             |
| System                                                                                                      |
|   assemble()                    50        0.2608      0.005216    0.4978      0.009955    9.54     18.21    |
|   calculate_norm()              51        0.1124      0.002205    0.1879      0.003684    4.11     6.87     |
|   project_vector()              76        0.0774      0.001018    0.2744      0.003610    2.83     10.04    |
|                                                                                                             |
| XdrIO                                                                                                       |
|   write()                       1         0.0031      0.003091    0.0031      0.003091    0.11     0.11     |
 -------------------------------------------------------------------------------------------------------------
| Totals:                         932221    2.7341                                          100.00            |
 -------------------------------------------------------------------------------------------------------------

 
***** Finished first 25 steps, now read in saved solution and continue *****
 
Usage:
	 ./adaptivity_ex5-opt -init_timestep 0
OR
	 ./adaptivity_ex5-opt -read_solution -init_timestep 26

Running: ./adaptivity_ex5-opt -read_solution -n_timesteps 25 -output_freq 10 -init_timestep 25

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=743
    n_local_nodes()=743
  n_elem()=892
    n_local_elem()=892
    n_active_elem()=670
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0

Initial H1 norm = 5.08479

 EquationSystems
  n_systems()=1
   System #0, "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=743
    n_local_dofs()=743
    n_constrained_dofs()=110
    n_local_constrained_dofs()=110
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 9.41454
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 18
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 110
      Average DoF Constraint Length= 2
      Number of Node Constraints = 217
      Maximum Node Constraint Length= 5
      Average Node Constraint Length= 2.52074

 Solving time step 25, time=0.6500...
H1 norm = 4.94524
  Refining the mesh...
H1 norm = 4.94485
 Solving time step 26, time=0.6750...
H1 norm = 4.81283
  Refining the mesh...
H1 norm = 4.81226
 Solving time step 27, time=0.7000...
H1 norm = 4.68713
  Refining the mesh...
H1 norm = 4.68697
 Solving time step 28, time=0.7250...
H1 norm = 4.56834
  Refining the mesh...
H1 norm = 4.5679
 Solving time step 29, time=0.7500...
H1 norm = 4.45515
  Refining the mesh...
H1 norm = 4.45494
 Solving time step 30, time=0.7750...
H1 norm = 4.34771
  Refining the mesh...
H1 norm = 4.34742
 Solving time step 31, time=0.8000...
H1 norm = 4.24533
  Refining the mesh...
H1 norm = 4.24515
 Solving time step 32, time=0.8250...
H1 norm = 4.14781
  Refining the mesh...
H1 norm = 4.14758
 Solving time step 33, time=0.8500...
H1 norm = 4.0547
  Refining the mesh...
H1 norm = 4.05446
 Solving time step 34, time=0.8750...
H1 norm = 3.96568
  Refining the mesh...
H1 norm = 3.96544
 Solving time step 35, time=0.9000...
H1 norm = 3.88053
  Refining the mesh...
H1 norm = 3.88038
 Solving time step 36, time=0.9250...
H1 norm = 3.7991
  Refining the mesh...
H1 norm = 3.79888
 Solving time step 37, time=0.9500...
H1 norm = 3.72095
  Refining the mesh...
H1 norm = 3.72082
 Solving time step 38, time=0.9750...
H1 norm = 3.64611
  Refining the mesh...
H1 norm = 3.64592
 Solving time step 39, time=1.0000...
H1 norm = 3.57419
  Refining the mesh...
H1 norm = 3.57396
 Solving time step 40, time=1.0300...
H1 norm = 3.50503
  Refining the mesh...
H1 norm = 3.50499
 Solving time step 41, time=1.0500...
H1 norm = 3.43875
  Refining the mesh...
H1 norm = 3.43856
 Solving time step 42, time=1.0700...
H1 norm = 3.37479
  Refining the mesh...
H1 norm = 3.37473
 Solving time step 43, time=1.1000...
H1 norm = 3.31335
  Refining the mesh...
H1 norm = 3.31322
 Solving time step 44, time=1.1200...
H1 norm = 3.2541
  Refining the mesh...
H1 norm = 3.25393
 Solving time step 45, time=1.1500...
H1 norm = 3.19688
  Refining the mesh...
H1 norm = 3.1968
 Solving time step 46, time=1.1700...
H1 norm = 3.14178
  Refining the mesh...
H1 norm = 3.14172
 Solving time step 47, time=1.2000...
H1 norm = 3.0886
  Refining the mesh...
H1 norm = 3.0885
 Solving time step 48, time=1.2200...
H1 norm = 3.03718
  Refining the mesh...
H1 norm = 3.03705
 Solving time step 49, time=1.2500...
H1 norm = 2.98746
  Refining the mesh...
H1 norm = 2.98734

-------------------------------------------------------------------
| Time:           Sat Apr  7 15:56:06 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=3.80212, Active time=3.47747                                                |
 -------------------------------------------------------------------------------------------------------------
| Event                           nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                           w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------|
|                                                                                                             |
|                                                                                                             |
| DofMap                                                                                                      |
|   add_neighbors_to_send_list()  26        0.0087      0.000336    0.0087      0.000336    0.25     0.25     |
|   build_constraint_matrix()     43136     0.0535      0.000001    0.0535      0.000001    1.54     1.54     |
|   cnstrn_elem_mat_vec()         43136     0.1270      0.000003    0.1270      0.000003    3.65     3.65     |
|   compute_sparsity()            26        0.0858      0.003299    0.1032      0.003970    2.47     2.97     |
|   create_dof_constraints()      26        0.1592      0.006123    0.2664      0.010245    4.58     7.66     |
|   distribute_dofs()             26        0.0116      0.000445    0.0517      0.001987    0.33     1.49     |
|   dof_indices()                 283516    0.2030      0.000001    0.2030      0.000001    5.84     5.84     |
|   enforce_constraints_exactly() 75        0.0057      0.000076    0.0057      0.000076    0.16     0.16     |
|   old_dof_indices()             132162    0.0774      0.000001    0.0774      0.000001    2.23     2.23     |
|   prepare_send_list()           26        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                      26        0.0401      0.001541    0.0401      0.001541    1.15     1.15     |
|                                                                                                             |
| EquationSystems                                                                                             |
|   build_solution_vector()       8         0.0090      0.001131    0.0143      0.001782    0.26     0.41     |
|   read()                        1         0.0140      0.013951    0.0452      0.045173    0.40     1.30     |
|   update()                      1         0.0001      0.000059    0.0001      0.000059    0.00     0.00     |
|                                                                                                             |
| ExodusII_IO                                                                                                 |
|   write_nodal_data()            4         0.0070      0.001748    0.0070      0.001748    0.20     0.20     |
|                                                                                                             |
| FE                                                                                                          |
|   compute_affine_map()          177698    0.1368      0.000001    0.1368      0.000001    3.93     3.93     |
|   compute_face_map()            46198     0.1404      0.000003    0.2791      0.000006    4.04     8.03     |
|   compute_shape_functions()     177698    0.0867      0.000000    0.0867      0.000000    2.49     2.49     |
|   init_face_shape_functions()   691       0.0013      0.000002    0.0013      0.000002    0.04     0.04     |
|   init_shape_functions()        90857     0.2267      0.000002    0.2267      0.000002    6.52     6.52     |
|   inverse_map()                 214824    0.2721      0.000001    0.2721      0.000001    7.83     7.83     |
|                                                                                                             |
| GMVIO                                                                                                       |
|   write_nodal_data()            4         0.0213      0.005334    0.0213      0.005334    0.61     0.61     |
|                                                                                                             |
| JumpErrorEstimator                                                                                          |
|   estimate_error()              25        0.6289      0.025158    1.4796      0.059184    18.09    42.55    |
|                                                                                                             |
| LocationMap                                                                                                 |
|   find()                        3672      0.0026      0.000001    0.0026      0.000001    0.07     0.07     |
|   init()                        50        0.0130      0.000260    0.0130      0.000260    0.37     0.37     |
|                                                                                                             |
| Mesh                                                                                                        |
|   contract()                    25        0.0022      0.000089    0.0049      0.000198    0.06     0.14     |
|   find_neighbors()              26        0.0763      0.002934    0.0763      0.002934    2.19     2.19     |
|   renumber_nodes_and_elem()     77        0.0067      0.000087    0.0067      0.000087    0.19     0.19     |
|                                                                                                             |
| MeshCommunication                                                                                           |
|   assign_global_indices()       1         0.0153      0.015290    0.0153      0.015302    0.44     0.44     |
|                                                                                                             |
| MeshOutput                                                                                                  |
|   write_equation_systems()      8         0.0001      0.000015    0.0427      0.005339    0.00     1.23     |
|                                                                                                             |
| MeshRefinement                                                                                              |
|   _coarsen_elements()           50        0.0023      0.000046    0.0023      0.000046    0.07     0.07     |
|   _refine_elements()            50        0.0136      0.000271    0.0219      0.000437    0.39     0.63     |
|   add_point()                   3672      0.0045      0.000001    0.0077      0.000002    0.13     0.22     |
|   make_coarsening_compatible()  58        0.0699      0.001204    0.0764      0.001317    2.01     2.20     |
|   make_refinement_compatible()  58        0.0024      0.000041    0.0025      0.000043    0.07     0.07     |
|                                                                                                             |
| Parallel                                                                                                    |
|   allgather()                   32        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   send_receive()                4         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|                                                                                                             |
| Partitioner                                                                                                 |
|   set_node_processor_ids()      1         0.0004      0.000422    0.0004      0.000422    0.01     0.01     |
|   single_partition()            26        0.0021      0.000081    0.0021      0.000081    0.06     0.06     |
|                                                                                                             |
| PetscLinearSolver                                                                                           |
|   solve()                       50        0.0607      0.001214    0.0607      0.001214    1.75     1.75     |
|                                                                                                             |
| PointLocatorTree                                                                                            |
|   init(no master)               50        0.0857      0.001714    0.0857      0.001714    2.46     2.46     |
|   operator()                    2808      0.0130      0.000005    0.0179      0.000006    0.37     0.51     |
|                                                                                                             |
| ProjectVector                                                                                               |
|   operator()                    75        0.1381      0.001841    0.2396      0.003195    3.97     6.89     |
|                                                                                                             |
| System                                                                                                      |
|   assemble()                    50        0.3913      0.007826    0.7189      0.014378    11.25    20.67    |
|   calculate_norm()              51        0.1786      0.003503    0.2922      0.005730    5.14     8.40     |
|   project_vector()              75        0.0794      0.001059    0.3781      0.005042    2.28     10.87    |
|                                                                                                             |
| XdrIO                                                                                                       |
|   read()                        1         0.0031      0.003133    0.0031      0.003133    0.09     0.09     |
 -------------------------------------------------------------------------------------------------------------
| Totals:                         1221160   3.4775                                          100.00            |
 -------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./adaptivity_ex5-opt [-read_solution] -n_timesteps 25 -init_timestep [0|25]
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
