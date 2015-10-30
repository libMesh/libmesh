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
        #include "serial_mesh.h"
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
  #include <B><FONT COLOR="#BC8F8F">&quot;serial_mesh.h&quot;</FONT></B>
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
Linking adaptivity_ex5-opt...
***************************************************************
* Running Example  mpirun -np 6 ./adaptivity_ex5-opt [-read_solution] -n_timesteps 25 -n_refinements 5 -init_timestep [0|25] -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Usage:
	 ./adaptivity_ex5-opt -init_timestep 0
OR
	 ./adaptivity_ex5-opt -read_solution -init_timestep 26

Running: ./adaptivity_ex5-opt -n_timesteps 25 -n_refinements 5 -max_r_steps 1 -output_freq 10 -init_timestep 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=4225
    n_local_nodes()=742
  n_elem()=5460
    n_local_elem()=931
    n_active_elem()=4096
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0

Initial H1 norm = 1.74491

 EquationSystems
  n_systems()=1
   System #0, "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=4225
    n_local_dofs()=742
    n_constrained_dofs()=132
    n_local_constrained_dofs()=21
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.5283
      Average Off-Processor Bandwidth <= 0.673854
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 9
    DofMap Constraints
      Number of DoF Constraints = 129
      Average DoF Constraint Length= 1
      Number of Node Constraints = 129
      Maximum Node Constraint Length= 2
      Average Node Constraint Length= 2

 Solving time step  0, time=0.0250...
H1 norm = 1.59167
  Refining the mesh...
H1 norm = 1.59174
 Solving time step  1, time=0.0500...
H1 norm = 1.46394
  Refining the mesh...
H1 norm = 1.4632
 Solving time step  2, time=0.0750...
H1 norm = 1.35431
  Refining the mesh...
H1 norm = 1.35394
 Solving time step  3, time=0.1000...
H1 norm = 1.26011
  Refining the mesh...
H1 norm = 1.25949
 Solving time step  4, time=0.1250...
H1 norm = 1.17754
  Refining the mesh...
H1 norm = 1.17736
 Solving time step  5, time=0.1500...
H1 norm = 1.10539
  Refining the mesh...
H1 norm = 1.10506
 Solving time step  6, time=0.1750...
H1 norm = 1.0413
  Refining the mesh...
H1 norm = 1.04097
 Solving time step  7, time=0.2000...
H1 norm = 0.984134
  Refining the mesh...
H1 norm = 0.983841
 Solving time step  8, time=0.2250...
H1 norm = 0.932868
  Refining the mesh...
H1 norm = 0.932689
 Solving time step  9, time=0.2500...
H1 norm = 0.886679
  Refining the mesh...
H1 norm = 0.886531
 Solving time step 10, time=0.2750...
H1 norm = 0.844875
  Refining the mesh...
H1 norm = 0.844701
 Solving time step 11, time=0.3000...
H1 norm = 0.806755
  Refining the mesh...
H1 norm = 0.806625
 Solving time step 12, time=0.3250...
H1 norm = 0.771957
  Refining the mesh...
H1 norm = 0.771864
 Solving time step 13, time=0.3500...
H1 norm = 0.740067
  Refining the mesh...
H1 norm = 0.739885
 Solving time step 14, time=0.3750...
H1 norm = 0.710575
  Refining the mesh...
H1 norm = 0.7105
 Solving time step 15, time=0.4000...
H1 norm = 0.683442
  Refining the mesh...
H1 norm = 0.683344
 Solving time step 16, time=0.4250...
H1 norm = 0.658275
  Refining the mesh...
H1 norm = 0.658236
 Solving time step 17, time=0.4500...
H1 norm = 0.63497
  Refining the mesh...
H1 norm = 0.634876
 Solving time step 18, time=0.4750...
H1 norm = 0.613199
  Refining the mesh...
H1 norm = 0.613134
 Solving time step 19, time=0.5000...
H1 norm = 0.592893
  Refining the mesh...
H1 norm = 0.592805
 Solving time step 20, time=0.5250...
H1 norm = 0.573854
  Refining the mesh...
H1 norm = 0.573837
 Solving time step 21, time=0.5500...
H1 norm = 0.556088
  Refining the mesh...
H1 norm = 0.556018
 Solving time step 22, time=0.5750...
H1 norm = 0.539335
  Refining the mesh...
H1 norm = 0.539265
 Solving time step 23, time=0.6000...
H1 norm = 0.52356
  Refining the mesh...
H1 norm = 0.523518
 Solving time step 24, time=0.6250...
H1 norm = 0.508706
  Refining the mesh...
H1 norm = 0.50866
Final H1 norm = 0.50866

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./adaptivity_ex5-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:20:33 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           2.044e+00      1.00044   2.043e+00
Objects:              2.520e+03      1.00000   2.520e+03
Flops:                1.021e+07      1.62910   8.251e+06  4.951e+07
Flops/sec:            4.995e+06      1.62862   4.038e+06  2.423e+07
MPI Messages:         7.740e+03      1.04899   7.561e+03  4.537e+04
MPI Message Lengths:  1.288e+06      1.11049   1.596e+02  7.241e+06
MPI Reductions:       6.408e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.0432e+00 100.0%  4.9508e+07 100.0%  4.537e+04 100.0%  1.596e+02      100.0%  5.485e+03  85.6% 

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
      %T - percent time in this phase         %F - percent flops in this phase
      %M - percent messages in this phase     %L - percent message lengths in this phase
      %R - percent reductions in this phase
   Total Mflop/s: 10e-6 * (sum of flops over all processors)/(max time over all processors)
------------------------------------------------------------------------------------------------------------------------
Event                Count      Time (sec)     Flops                             --- Global ---  --- Stage ---   Total
                   Max Ratio  Max     Ratio   Max  Ratio  Mess   Avg len Reduct  %T %F %M %L %R  %T %F %M %L %R Mflop/s
------------------------------------------------------------------------------------------------------------------------

--- Event Stage 0: Main Stage

VecMDot              502 1.0 2.3740e-02 1.4 6.95e+05 1.4 0.0e+00 0.0e+00 5.0e+02  1  7  0  0  8   1  7  0  0  9   152
VecNorm              602 1.0 3.8031e-02 1.3 1.52e+05 1.4 0.0e+00 0.0e+00 6.0e+02  2  2  0  0  9   2  2  0  0 11    21
VecScale             552 1.0 3.1924e-04 1.5 6.94e+04 1.4 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  1130
VecCopy              401 1.0 2.1124e-04 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               884 1.0 4.0722e-04 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              100 1.0 6.7115e-0357.1 2.54e+04 1.4 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0    20
VecMAXPY             552 1.0 5.2595e-04 1.6 8.24e+05 1.4 0.0e+00 0.0e+00 0.0e+00  0  9  0  0  0   0  9  0  0  0  8129
VecAssemblyBegin     804 1.0 1.3380e-01 1.3 0.00e+00 0.0 2.0e+03 1.2e+02 2.2e+03  5  0  4  3 34   5  0  4  3 40     0
VecAssemblyEnd       804 1.0 6.4564e-04 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     1231 1.0 3.5672e-03 1.4 0.00e+00 0.0 2.7e+04 1.6e+02 0.0e+00  0  0 59 58  0   0  0 59 58  0     0
VecScatterEnd       1231 1.0 8.4196e-02 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  4  0  0  0  0   4  0  0  0  0     0
VecNormalize         552 1.0 3.4901e-02 1.5 2.08e+05 1.4 0.0e+00 0.0e+00 5.5e+02  1  2  0  0  9   1  2  0  0 10    31
MatMult              552 1.0 4.6558e-02 1.3 1.32e+06 1.4 1.4e+04 1.0e+02 0.0e+00  2 14 31 20  0   2 14 31 20  0   144
MatSolve             552 1.0 4.7712e-03 1.9 4.26e+06 1.7 0.0e+00 0.0e+00 0.0e+00  0 41  0  0  0   0 41  0  0  0  4263
MatLUFactorNum        50 1.0 4.7946e-03 1.8 2.87e+06 1.9 0.0e+00 0.0e+00 0.0e+00  0 27  0  0  0   0 27  0  0  0  2776
MatILUFactorSym       50 1.0 1.3419e-02 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 5.0e+01  0  0  0  0  1   0  0  0  0  1     0
MatAssemblyBegin     100 1.0 8.0574e-02 1.5 0.00e+00 0.0 2.3e+03 4.5e+02 2.0e+02  3  0  5 14  3   3  0  5 14  4     0
MatAssemblyEnd       100 1.0 1.4910e-02 1.1 0.00e+00 0.0 1.3e+03 2.9e+01 2.6e+02  1  0  3  1  4   1  0  3  1  5     0
MatGetRowIJ           50 1.0 9.6798e-04119.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        50 1.0 1.3170e-03 4.1 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+02  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries       102 1.0 1.2493e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       502 1.0 2.4331e-02 1.4 1.39e+06 1.4 0.0e+00 0.0e+00 5.0e+02  1 15  0  0  8   1 15  0  0  9   297
KSPSetup             100 1.0 4.5395e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              50 1.0 1.2992e-01 1.1 1.02e+07 1.6 1.4e+04 1.0e+02 1.3e+03  6100 31 20 20   6100 31 20 23   381
PCSetUp              100 1.0 2.2528e-02 1.8 2.87e+06 1.9 0.0e+00 0.0e+00 1.5e+02  1 27  0  0  2   1 27  0  0  3   591
PCSetUpOnBlocks       50 1.0 2.1026e-02 1.9 2.87e+06 1.9 0.0e+00 0.0e+00 1.5e+02  1 27  0  0  2   1 27  0  0  3   633
PCApply              552 1.0 8.7333e-03 1.5 4.26e+06 1.7 0.0e+00 0.0e+00 0.0e+00  0 41  0  0  0   0 41  0  0  0  2329
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec  1142           1142      3147288     0
         Vec Scatter   357            357       309876     0
           Index Set   661            661       462096     0
   IS L to G Mapping   128            128       159580     0
              Matrix   128            128      3196812     0
       Krylov Solver    52             52       490880     0
      Preconditioner    52             52        36608     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 3.9196e-05
Average time for zero size MPI_Send(): 2.33253e-05
#PETSc Option Table entries:
-init_timestep 0
-ksp_right_pc
-log_summary
-max_r_steps 1
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
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8
Configure run at: Sat May 19 03:47:23 2012
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared=1 --with-shared-libraries=1 --with-mpi-dir=/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid --with-mumps=true --download-mumps=1 --with-parmetis=true --download-parmetis=1 --with-superlu=true --download-superlu=1 --with-superludir=true --download-superlu_dist=1 --with-blacs=true --download-blacs=1 --with-scalapack=true --download-scalapack=1 --with-hypre=true --download-hypre=1 --with-blas-lib="[/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_intel_lp64.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_sequential.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_core.so]" --with-lapack-lib=/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_solver_lp64_sequential.a
-----------------------------------------
Libraries compiled on Sat May 19 03:47:23 CDT 2012 on daedalus 
Machine characteristics: Linux daedalus 2.6.32-34-generic #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011 x86_64 GNU/Linux 
Using PETSc directory: /org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5
Using PETSc arch: intel-11.1-lucid-mpich2-1.4.1-cxx-opt
-----------------------------------------
Using C compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpicxx -O3   -fPIC   
Using Fortran compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpif90 -fPIC -O3    
-----------------------------------------
Using include paths: -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/include  
------------------------------------------
Using C linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpicxx -O3 
Using Fortran linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpif90 -fPIC -O3  
Using libraries: -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -lpetsc       -lX11 -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -lHYPRE -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lblacs -lsuperlu_dist_2.4 -lparmetis -lmetis -lsuperlu_4.0 -Wl,-rpath,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t -L/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t -lmkl_solver_lp64_sequential -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -ldl -Wl,-rpath,/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/lib -L/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/lib -lmpich -lopa -lmpl -lrt -lpthread -Wl,-rpath,/opt/intel/Compiler/11.1/073/lib/intel64 -L/opt/intel/Compiler/11.1/073/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -L/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -lmpichf90 -lifport -lifcore -lm -lm -lmpichcxx -lstdc++ -lmpichcxx -lstdc++ -ldl -lmpich -lopa -lmpl -lrt -lpthread -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -ldl  
------------------------------------------

-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 6                                                |
| Time:           Fri Aug 24 15:20:33 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=2.19969, Active time=1.8436                                                    |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     26        0.0036      0.000138    0.0043      0.000166    0.19     0.23     |
|   build_constraint_matrix()        4902      0.0059      0.000001    0.0059      0.000001    0.32     0.32     |
|   build_sparsity()                 26        0.0255      0.000981    0.0325      0.001249    1.38     1.76     |
|   cnstrn_elem_mat_vec()            4902      0.0127      0.000003    0.0127      0.000003    0.69     0.69     |
|   create_dof_constraints()         26        0.1073      0.004128    0.2082      0.008008    5.82     11.29    |
|   distribute_dofs()                26        0.0094      0.000360    0.0759      0.002918    0.51     4.12     |
|   dof_indices()                    60337     0.0250      0.000000    0.0250      0.000000    1.35     1.35     |
|   enforce_constraints_exactly()    76        0.0188      0.000247    0.0188      0.000247    1.02     1.02     |
|   old_dof_indices()                15591     0.0056      0.000000    0.0056      0.000000    0.31     0.31     |
|   prepare_send_list()              26        0.0001      0.000005    0.0001      0.000005    0.01     0.01     |
|   reinit()                         26        0.0214      0.000824    0.0214      0.000824    1.16     1.16     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          8         0.0026      0.000321    0.0075      0.000939    0.14     0.41     |
|   write()                          1         0.0002      0.000232    0.0017      0.001672    0.01     0.09     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               4         0.0055      0.001376    0.0055      0.001376    0.30     0.30     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        24060     0.0183      0.000001    0.0183      0.000001    0.99     0.99     |
|   init_shape_functions()           13566     0.0147      0.000001    0.0147      0.000001    0.80     0.80     |
|   inverse_map()                    46705     0.0363      0.000001    0.0363      0.000001    1.97     1.97     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             24060     0.0167      0.000001    0.0167      0.000001    0.91     0.91     |
|   compute_face_map()               6732      0.0163      0.000002    0.0301      0.000004    0.88     1.63     |
|   init_face_shape_functions()      1199      0.0012      0.000001    0.0012      0.000001    0.07     0.07     |
|   init_reference_to_physical_map() 13566     0.0207      0.000002    0.0207      0.000002    1.12     1.12     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               4         0.0783      0.019566    0.0783      0.019566    4.25     4.25     |
|                                                                                                                |
| JumpErrorEstimator                                                                                             |
|   estimate_error()                 25        0.0646      0.002584    0.2869      0.011478    3.50     15.56    |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           19140     0.0057      0.000000    0.0057      0.000000    0.31     0.31     |
|   init()                           55        0.0050      0.000091    0.0050      0.000091    0.27     0.27     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       25        0.0020      0.000078    0.0043      0.000171    0.11     0.23     |
|   find_neighbors()                 27        0.0368      0.001362    0.0867      0.003210    1.99     4.70     |
|   renumber_nodes_and_elem()        79        0.0066      0.000084    0.0066      0.000084    0.36     0.36     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0079      0.007906    0.0109      0.010928    0.43     0.59     |
|   compute_hilbert_indices()        28        0.0303      0.001081    0.0303      0.001081    1.64     1.64     |
|   find_global_indices()            28        0.0045      0.000162    0.0925      0.003302    0.25     5.01     |
|   parallel_sort()                  28        0.0101      0.000362    0.0383      0.001368    0.55     2.08     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         8         0.0001      0.000012    0.0914      0.011422    0.01     4.96     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              50        0.0021      0.000041    0.0055      0.000110    0.11     0.30     |
|   _refine_elements()               55        0.0180      0.000327    0.0680      0.001237    0.97     3.69     |
|   add_point()                      19140     0.0131      0.000001    0.0201      0.000001    0.71     1.09     |
|   make_coarsening_compatible()     61        0.0375      0.000614    0.0513      0.000841    2.03     2.78     |
|   make_refinement_compatible()     61        0.0012      0.000019    0.0090      0.000147    0.06     0.49     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0001      0.000053    0.0001      0.000053    0.00     0.00     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      27        0.0391      0.001447    0.1443      0.005345    2.12     7.83     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      141       0.0374      0.000265    0.0374      0.000265    2.03     2.03     |
|   barrier()                        1         0.0000      0.000037    0.0000      0.000037    0.00     0.00     |
|   broadcast()                      2         0.0000      0.000010    0.0000      0.000010    0.00     0.00     |
|   gather()                         20        0.0023      0.000114    0.0023      0.000114    0.12     0.12     |
|   max(bool)                        191       0.0434      0.000227    0.0434      0.000227    2.35     2.35     |
|   max(scalar)                      80        0.1525      0.001906    0.1525      0.001906    8.27     8.27     |
|   max(vector)                      129       0.0089      0.000069    0.0089      0.000069    0.48     0.48     |
|   min(bool)                        122       0.0732      0.000600    0.0732      0.000600    3.97     3.97     |
|   min(scalar)                      25        0.0028      0.000110    0.0028      0.000110    0.15     0.15     |
|   min(vector)                      129       0.0261      0.000202    0.0261      0.000202    1.42     1.42     |
|   probe()                          1175      0.0628      0.000053    0.0628      0.000053    3.41     3.41     |
|   receive()                        1209      0.0019      0.000002    0.0649      0.000054    0.10     3.52     |
|   send()                           1144      0.0010      0.000001    0.0010      0.000001    0.05     0.05     |
|   send_receive()                   1200      0.0024      0.000002    0.0689      0.000057    0.13     3.73     |
|   sum()                            274       0.1982      0.000723    0.1982      0.000723    10.75    10.75    |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           1174      0.0005      0.000000    0.0005      0.000000    0.03     0.03     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         27        0.0044      0.000163    0.0423      0.001566    0.24     2.29     |
|   set_parent_processor_ids()       27        0.0031      0.000115    0.0031      0.000115    0.17     0.17     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          50        0.2241      0.004481    0.2241      0.004481    12.15    12.15    |
|                                                                                                                |
| PointLocatorTree                                                                                               |
|   init(no master)                  50        0.0366      0.000731    0.0569      0.001138    1.98     3.09     |
|   operator()                       5460      0.0376      0.000007    0.0437      0.000008    2.04     2.37     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       75        0.0122      0.000163    0.0181      0.000242    0.66     0.98     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       50        0.0326      0.000652    0.0665      0.001329    1.77     3.61     |
|   calculate_norm()                 52        0.0346      0.000665    0.0710      0.001366    1.88     3.85     |
|   project_vector()                 76        0.1134      0.001492    0.1573      0.002070    6.15     8.53     |
|                                                                                                                |
| XdrIO                                                                                                          |
|   write()                          1         0.0014      0.001356    0.0035      0.003540    0.07     0.19     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            267592    1.8436                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***** Finished first 25 steps, now read in saved solution and continue *****
 
Usage:
	 ./adaptivity_ex5-opt -init_timestep 0
OR
	 ./adaptivity_ex5-opt -read_solution -init_timestep 26

Running: ./adaptivity_ex5-opt -read_solution -n_timesteps 25 -max_r_steps 1 -output_freq 10 -init_timestep 25 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=737
    n_local_nodes()=138
  n_elem()=884
    n_local_elem()=162
    n_active_elem()=664
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0

Initial H1 norm = 0.50866

 EquationSystems
  n_systems()=1
   System #0, "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=737
    n_local_dofs()=138
    n_constrained_dofs()=135
    n_local_constrained_dofs()=19
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.52174
      Average Off-Processor Bandwidth <= 1.06522
      Maximum  On-Processor Bandwidth <= 15
      Maximum Off-Processor Bandwidth <= 9
    DofMap Constraints
      Number of DoF Constraints = 133
      Average DoF Constraint Length= 1.90226
      Number of Node Constraints = 236
      Maximum Node Constraint Length= 5
      Average Node Constraint Length= 2.51271

 Solving time step 25, time=0.6500...
H1 norm = 0.494674
  Refining the mesh...
H1 norm = 0.494645
 Solving time step 26, time=0.6750...
H1 norm = 0.481418
  Refining the mesh...
H1 norm = 0.481365
 Solving time step 27, time=0.7000...
H1 norm = 0.468827
  Refining the mesh...
H1 norm = 0.468812
 Solving time step 28, time=0.7250...
H1 norm = 0.456925
  Refining the mesh...
H1 norm = 0.456884
 Solving time step 29, time=0.7500...
H1 norm = 0.445586
  Refining the mesh...
H1 norm = 0.445562
 Solving time step 30, time=0.7750...
H1 norm = 0.434817
  Refining the mesh...
H1 norm = 0.434788
 Solving time step 31, time=0.8000...
H1 norm = 0.424559
  Refining the mesh...
H1 norm = 0.424541
 Solving time step 32, time=0.8250...
H1 norm = 0.414788
  Refining the mesh...
H1 norm = 0.414764
 Solving time step 33, time=0.8500...
H1 norm = 0.405457
  Refining the mesh...
H1 norm = 0.405436
 Solving time step 34, time=0.8750...
H1 norm = 0.396542
  Refining the mesh...
H1 norm = 0.396515
 Solving time step 35, time=0.9000...
H1 norm = 0.388006
  Refining the mesh...
H1 norm = 0.387991
 Solving time step 36, time=0.9250...
H1 norm = 0.379846
  Refining the mesh...
H1 norm = 0.379825
 Solving time step 37, time=0.9500...
H1 norm = 0.372017
  Refining the mesh...
H1 norm = 0.372004
 Solving time step 38, time=0.9750...
H1 norm = 0.364518
  Refining the mesh...
H1 norm = 0.364499
 Solving time step 39, time=1.0000...
H1 norm = 0.357312
  Refining the mesh...
H1 norm = 0.35729
 Solving time step 40, time=1.0300...
H1 norm = 0.350383
  Refining the mesh...
H1 norm = 0.350379
 Solving time step 41, time=1.0500...
H1 norm = 0.343743
  Refining the mesh...
H1 norm = 0.343723
 Solving time step 42, time=1.0700...
H1 norm = 0.337334
  Refining the mesh...
H1 norm = 0.337328
 Solving time step 43, time=1.1000...
H1 norm = 0.331179
  Refining the mesh...
H1 norm = 0.331165
 Solving time step 44, time=1.1200...
H1 norm = 0.325241
  Refining the mesh...
H1 norm = 0.325226
 Solving time step 45, time=1.1500...
H1 norm = 0.319511
  Refining the mesh...
H1 norm = 0.319502
 Solving time step 46, time=1.1700...
H1 norm = 0.313989
  Refining the mesh...
H1 norm = 0.313984
 Solving time step 47, time=1.2000...
H1 norm = 0.308663
  Refining the mesh...
H1 norm = 0.308651
 Solving time step 48, time=1.2200...
H1 norm = 0.30351
  Refining the mesh...
H1 norm = 0.303496
 Solving time step 49, time=1.2500...
H1 norm = 0.298527
  Refining the mesh...
H1 norm = 0.298516
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./adaptivity_ex5-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:20:36 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           2.239e+00      1.00101   2.238e+00
Objects:              2.572e+03      1.00000   2.572e+03
Flops:                1.623e+07      1.51758   1.351e+07  8.106e+07
Flops/sec:            7.248e+06      1.51747   6.037e+06  3.622e+07
MPI Messages:         8.292e+03      1.06946   8.072e+03  4.843e+04
MPI Message Lengths:  1.577e+06      1.09471   1.870e+02  9.058e+06
MPI Reductions:       6.545e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.2379e+00 100.0%  8.1061e+07 100.0%  4.843e+04 100.0%  1.870e+02      100.0%  5.598e+03  85.5% 

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
      %T - percent time in this phase         %F - percent flops in this phase
      %M - percent messages in this phase     %L - percent message lengths in this phase
      %R - percent reductions in this phase
   Total Mflop/s: 10e-6 * (sum of flops over all processors)/(max time over all processors)
------------------------------------------------------------------------------------------------------------------------
Event                Count      Time (sec)     Flops                             --- Global ---  --- Stage ---   Total
                   Max Ratio  Max     Ratio   Max  Ratio  Mess   Avg len Reduct  %T %F %M %L %R  %T %F %M %L %R Mflop/s
------------------------------------------------------------------------------------------------------------------------

--- Event Stage 0: Main Stage

VecMDot              506 1.0 2.7873e-02 1.5 1.03e+06 1.3 0.0e+00 0.0e+00 5.1e+02  1  7  0  0  8   1  7  0  0  9   193
VecNorm              606 1.0 3.8195e-02 1.4 2.18e+05 1.3 0.0e+00 0.0e+00 6.1e+02  1  1  0  0  9   1  1  0  0 11    30
VecScale             556 1.0 3.1662e-04 1.3 9.99e+04 1.3 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  1655
VecCopy              408 1.0 1.8358e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               897 1.0 3.3975e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              100 1.0 1.5306e-04 1.5 3.58e+04 1.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1227
VecMAXPY             556 1.0 9.2959e-04 1.9 1.21e+06 1.3 0.0e+00 0.0e+00 0.0e+00  0  8  0  0  0   0  8  0  0  0  6832
VecAssemblyBegin     827 1.0 1.4940e-01 1.4 0.00e+00 0.0 2.3e+03 1.5e+02 2.2e+03  6  0  5  4 34   6  0  5  4 40     0
VecAssemblyEnd       827 1.0 6.0463e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     1255 1.0 3.4311e-03 1.1 0.00e+00 0.0 2.9e+04 1.9e+02 0.0e+00  0  0 59 59  0   0  0 59 59  0     0
VecScatterEnd       1255 1.0 1.1676e-01 2.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  4  0  0  0  0   4  0  0  0  0     0
VecNormalize         556 1.0 3.4880e-02 1.5 3.00e+05 1.3 0.0e+00 0.0e+00 5.6e+02  1  2  0  0  8   1  2  0  0 10    45
MatMult              556 1.0 5.9806e-02 2.2 1.83e+06 1.4 1.5e+04 1.2e+02 0.0e+00  2 12 32 20  0   2 12 32 20  0   160
MatSolve             556 1.0 6.9680e-03 1.6 6.69e+06 1.5 0.0e+00 0.0e+00 0.0e+00  0 41  0  0  0   0 41  0  0  0  4770
MatLUFactorNum        50 1.0 7.1611e-03 1.6 5.11e+06 1.7 0.0e+00 0.0e+00 0.0e+00  0 30  0  0  0   0 30  0  0  0  3448
MatILUFactorSym       50 1.0 1.9248e-02 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 5.0e+01  1  0  0  0  1   1  0  0  0  1     0
MatAssemblyBegin     100 1.0 8.1418e-02 1.6 0.00e+00 0.0 2.6e+03 4.6e+02 2.0e+02  3  0  5 13  3   3  0  5 13  4     0
MatAssemblyEnd       100 1.0 1.5559e-02 1.1 0.00e+00 0.0 1.4e+03 3.2e+01 2.6e+02  1  0  3  0  4   1  0  3  0  5     0
MatGetRowIJ           50 1.0 1.2875e-05 2.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        50 1.0 3.7074e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+02  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries       104 1.0 1.1754e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       506 1.0 2.8652e-02 1.4 2.05e+06 1.3 0.0e+00 0.0e+00 5.1e+02  1 13  0  0  8   1 13  0  0  9   376
KSPSetup             100 1.0 4.4537e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              50 1.0 1.4495e-01 1.1 1.62e+07 1.5 1.5e+04 1.2e+02 1.3e+03  6100 32 20 19   6100 32 20 23   559
PCSetUp              100 1.0 2.9411e-02 1.4 5.11e+06 1.7 0.0e+00 0.0e+00 1.5e+02  1 30  0  0  2   1 30  0  0  3   840
PCSetUpOnBlocks       50 1.0 2.8068e-02 1.5 5.11e+06 1.7 0.0e+00 0.0e+00 1.5e+02  1 30  0  0  2   1 30  0  0  3   880
PCApply              556 1.0 1.0597e-02 1.3 6.69e+06 1.5 0.0e+00 0.0e+00 0.0e+00  0 41  0  0  0   0 41  0  0  0  3136
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec  1163           1163      3689832     0
         Vec Scatter   366            366       317688     0
           Index Set   675            675       509228     0
   IS L to G Mapping   133            133       192560     0
              Matrix   131            131      4657900     0
       Krylov Solver    52             52       490880     0
      Preconditioner    52             52        36608     0
========================================================================================================================
Average time to get PetscTime(): 1.19209e-07
Average time for MPI_Barrier(): 1.81198e-06
Average time for zero size MPI_Send(): 1.1007e-05
#PETSc Option Table entries:
-init_timestep 25
-ksp_right_pc
-log_summary
-max_r_steps 1
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
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8
Configure run at: Sat May 19 03:47:23 2012
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared=1 --with-shared-libraries=1 --with-mpi-dir=/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid --with-mumps=true --download-mumps=1 --with-parmetis=true --download-parmetis=1 --with-superlu=true --download-superlu=1 --with-superludir=true --download-superlu_dist=1 --with-blacs=true --download-blacs=1 --with-scalapack=true --download-scalapack=1 --with-hypre=true --download-hypre=1 --with-blas-lib="[/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_intel_lp64.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_sequential.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_core.so]" --with-lapack-lib=/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_solver_lp64_sequential.a
-----------------------------------------
Libraries compiled on Sat May 19 03:47:23 CDT 2012 on daedalus 
Machine characteristics: Linux daedalus 2.6.32-34-generic #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011 x86_64 GNU/Linux 
Using PETSc directory: /org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5
Using PETSc arch: intel-11.1-lucid-mpich2-1.4.1-cxx-opt
-----------------------------------------
Using C compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpicxx -O3   -fPIC   
Using Fortran compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpif90 -fPIC -O3    
-----------------------------------------
Using include paths: -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/include  
------------------------------------------
Using C linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpicxx -O3 
Using Fortran linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpif90 -fPIC -O3  
Using libraries: -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -lpetsc       -lX11 -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -lHYPRE -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lblacs -lsuperlu_dist_2.4 -lparmetis -lmetis -lsuperlu_4.0 -Wl,-rpath,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t -L/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t -lmkl_solver_lp64_sequential -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -ldl -Wl,-rpath,/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/lib -L/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/lib -lmpich -lopa -lmpl -lrt -lpthread -Wl,-rpath,/opt/intel/Compiler/11.1/073/lib/intel64 -L/opt/intel/Compiler/11.1/073/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -L/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -lmpichf90 -lifport -lifcore -lm -lm -lmpichcxx -lstdc++ -lmpichcxx -lstdc++ -ldl -lmpich -lopa -lmpl -lrt -lpthread -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -ldl  
------------------------------------------

-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 6                                                |
| Time:           Fri Aug 24 15:20:36 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=2.35644, Active time=2.04251                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     27        0.0041      0.000151    0.0049      0.000181    0.20     0.24     |
|   build_constraint_matrix()        7155      0.0054      0.000001    0.0054      0.000001    0.26     0.26     |
|   build_sparsity()                 27        0.0366      0.001357    0.0465      0.001723    1.79     2.28     |
|   cnstrn_elem_mat_vec()            7155      0.0128      0.000002    0.0128      0.000002    0.62     0.62     |
|   create_dof_constraints()         27        0.0897      0.003324    0.1514      0.005607    4.39     7.41     |
|   distribute_dofs()                27        0.0086      0.000318    0.0724      0.002683    0.42     3.55     |
|   dof_indices()                    81532     0.0277      0.000000    0.0277      0.000000    1.36     1.36     |
|   enforce_constraints_exactly()    78        0.0198      0.000254    0.0198      0.000254    0.97     0.97     |
|   old_dof_indices()                22722     0.0069      0.000000    0.0069      0.000000    0.34     0.34     |
|   prepare_send_list()              27        0.0002      0.000006    0.0002      0.000006    0.01     0.01     |
|   reinit()                         27        0.0196      0.000726    0.0196      0.000726    0.96     0.96     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          8         0.0015      0.000192    0.0103      0.001292    0.08     0.51     |
|   read()                           1         0.0093      0.009295    0.0334      0.033418    0.46     1.64     |
|   update()                         1         0.0003      0.000287    0.0003      0.000287    0.01     0.01     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               4         0.0051      0.001269    0.0051      0.001269    0.25     0.25     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        30558     0.0212      0.000001    0.0212      0.000001    1.04     1.04     |
|   init_shape_functions()           16239     0.0161      0.000001    0.0161      0.000001    0.79     0.79     |
|   inverse_map()                    54148     0.0354      0.000001    0.0354      0.000001    1.73     1.73     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             30558     0.0173      0.000001    0.0173      0.000001    0.85     0.85     |
|   compute_face_map()               8069      0.0170      0.000002    0.0318      0.000004    0.83     1.55     |
|   init_face_shape_functions()      659       0.0006      0.000001    0.0006      0.000001    0.03     0.03     |
|   init_reference_to_physical_map() 16239     0.0204      0.000001    0.0204      0.000001    1.00     1.00     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               4         0.0853      0.021326    0.0853      0.021326    4.18     4.18     |
|                                                                                                                |
| JumpErrorEstimator                                                                                             |
|   estimate_error()                 25        0.0774      0.003097    0.3618      0.014472    3.79     17.71    |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           3672      0.0011      0.000000    0.0011      0.000000    0.05     0.05     |
|   init()                           51        0.0056      0.000110    0.0056      0.000110    0.27     0.27     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       26        0.0010      0.000039    0.0026      0.000098    0.05     0.13     |
|   find_neighbors()                 26        0.0372      0.001432    0.0767      0.002950    1.82     3.75     |
|   renumber_nodes_and_elem()        78        0.0052      0.000066    0.0052      0.000066    0.25     0.25     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0067      0.006677    0.0142      0.014244    0.33     0.70     |
|   compute_hilbert_indices()        26        0.0339      0.001304    0.0339      0.001304    1.66     1.66     |
|   find_global_indices()            26        0.0044      0.000171    0.1148      0.004415    0.22     5.62     |
|   parallel_sort()                  26        0.0080      0.000308    0.0563      0.002164    0.39     2.75     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         8         0.0001      0.000008    0.1008      0.012598    0.00     4.93     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              51        0.0012      0.000024    0.0045      0.000087    0.06     0.22     |
|   _refine_elements()               51        0.0062      0.000122    0.0309      0.000606    0.30     1.51     |
|   add_point()                      3672      0.0026      0.000001    0.0039      0.000001    0.13     0.19     |
|   make_coarsening_compatible()     60        0.0428      0.000713    0.0462      0.000770    2.10     2.26     |
|   make_refinement_compatible()     60        0.0013      0.000022    0.0103      0.000172    0.06     0.51     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      26        0.0436      0.001678    0.1761      0.006774    2.14     8.62     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      141       0.0283      0.000200    0.0283      0.000200    1.38     1.38     |
|   broadcast()                      61        0.0028      0.000046    0.0028      0.000046    0.14     0.14     |
|   gather()                         1         0.0000      0.000021    0.0000      0.000021    0.00     0.00     |
|   max(bool)                        188       0.0362      0.000193    0.0362      0.000193    1.77     1.77     |
|   max(scalar)                      78        0.2005      0.002571    0.2005      0.002571    9.82     9.82     |
|   max(vector)                      129       0.0099      0.000077    0.0099      0.000077    0.48     0.48     |
|   min(bool)                        120       0.0833      0.000694    0.0833      0.000694    4.08     4.08     |
|   min(scalar)                      25        0.0029      0.000115    0.0029      0.000115    0.14     0.14     |
|   min(vector)                      129       0.0261      0.000203    0.0261      0.000203    1.28     1.28     |
|   probe()                          1090      0.0868      0.000080    0.0868      0.000080    4.25     4.25     |
|   receive()                        1090      0.0016      0.000001    0.0885      0.000081    0.08     4.33     |
|   send()                           1090      0.0007      0.000001    0.0007      0.000001    0.03     0.03     |
|   send_receive()                   1146      0.0020      0.000002    0.0918      0.000080    0.10     4.50     |
|   sum()                            263       0.2761      0.001050    0.2761      0.001050    13.52    13.52    |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           1090      0.0004      0.000000    0.0004      0.000000    0.02     0.02     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         27        0.0043      0.000161    0.0526      0.001950    0.21     2.58     |
|   set_parent_processor_ids()       26        0.0033      0.000127    0.0033      0.000127    0.16     0.16     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          50        0.2533      0.005065    0.2533      0.005065    12.40    12.40    |
|                                                                                                                |
| PointLocatorTree                                                                                               |
|   init(no master)                  51        0.0448      0.000879    0.0630      0.001235    2.19     3.08     |
|   operator()                       2839      0.0068      0.000002    0.0096      0.000003    0.33     0.47     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       78        0.0181      0.000233    0.0268      0.000344    0.89     1.31     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       50        0.0417      0.000835    0.0781      0.001563    2.04     3.83     |
|   calculate_norm()                 51        0.0420      0.000824    0.0878      0.001722    2.06     4.30     |
|   project_vector()                 78        0.1302      0.001670    0.1840      0.002359    6.38     9.01     |
|                                                                                                                |
| XdrIO                                                                                                          |
|   read()                           1         0.0010      0.000969    0.0011      0.001057    0.05     0.05     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            293069    2.0425                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
Usage:
	 ./adaptivity_ex5-opt -init_timestep 0
OR
	 ./adaptivity_ex5-opt -read_solution -init_timestep 26

Running: ./adaptivity_ex5-opt -n_timesteps 25 -n_refinements 5 -max_r_steps 1 -output_freq 10 -init_timestep 0 -exact_solution 100+10*sin(pi*x)*cos(pi*y) -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=4225
    n_local_nodes()=742
  n_elem()=5460
    n_local_elem()=931
    n_active_elem()=4096
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0

Initial H1 norm = 205.107

 EquationSystems
  n_systems()=1
   System #0, "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=4225
    n_local_dofs()=742
    n_constrained_dofs()=132
    n_local_constrained_dofs()=21
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.5283
      Average Off-Processor Bandwidth <= 0.673854
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 9
    DofMap Constraints
      Number of DoF Constraints = 129
      Average DoF Constraint Length= 1
      Number of Node Constraints = 129
      Maximum Node Constraint Length= 2
      Average Node Constraint Length= 2

 Solving time step  0, time=0.0250...
H1 norm = 205.057
  Refining the mesh...
H1 norm = 205.057
 Solving time step  1, time=0.0500...
H1 norm = 205.008
  Refining the mesh...
H1 norm = 205.008
 Solving time step  2, time=0.0750...
H1 norm = 204.96
  Refining the mesh...
H1 norm = 204.96
 Solving time step  3, time=0.1000...
H1 norm = 204.912
  Refining the mesh...
H1 norm = 204.912
 Solving time step  4, time=0.1250...
H1 norm = 204.864
  Refining the mesh...
H1 norm = 204.864
 Solving time step  5, time=0.1500...
H1 norm = 204.817
  Refining the mesh...
H1 norm = 204.817
 Solving time step  6, time=0.1750...
H1 norm = 204.77
  Refining the mesh...
H1 norm = 204.77
 Solving time step  7, time=0.2000...
H1 norm = 204.724
  Refining the mesh...
H1 norm = 204.724
 Solving time step  8, time=0.2250...
H1 norm = 204.678
  Refining the mesh...
H1 norm = 204.678
 Solving time step  9, time=0.2500...
H1 norm = 204.633
  Refining the mesh...
H1 norm = 204.632
 Solving time step 10, time=0.2750...
H1 norm = 204.588
  Refining the mesh...
H1 norm = 204.587
 Solving time step 11, time=0.3000...
H1 norm = 204.543
  Refining the mesh...
H1 norm = 204.543
 Solving time step 12, time=0.3250...
H1 norm = 204.499
  Refining the mesh...
H1 norm = 204.499
 Solving time step 13, time=0.3500...
H1 norm = 204.455
  Refining the mesh...
H1 norm = 204.455
 Solving time step 14, time=0.3750...
H1 norm = 204.412
  Refining the mesh...
H1 norm = 204.412
 Solving time step 15, time=0.4000...
H1 norm = 204.369
  Refining the mesh...
H1 norm = 204.369
 Solving time step 16, time=0.4250...
H1 norm = 204.327
  Refining the mesh...
H1 norm = 204.326
 Solving time step 17, time=0.4500...
H1 norm = 204.284
  Refining the mesh...
H1 norm = 204.284
 Solving time step 18, time=0.4750...
H1 norm = 204.243
  Refining the mesh...
H1 norm = 204.243
 Solving time step 19, time=0.5000...
H1 norm = 204.202
  Refining the mesh...
H1 norm = 204.202
 Solving time step 20, time=0.5250...
H1 norm = 204.161
  Refining the mesh...
H1 norm = 204.161
 Solving time step 21, time=0.5500...
H1 norm = 204.12
  Refining the mesh...
H1 norm = 204.12
 Solving time step 22, time=0.5750...
H1 norm = 204.08
  Refining the mesh...
H1 norm = 204.08
 Solving time step 23, time=0.6000...
H1 norm = 204.041
  Refining the mesh...
H1 norm = 204.041
 Solving time step 24, time=0.6250...
H1 norm = 204.001
  Refining the mesh...
H1 norm = 204.001
Final H1 norm = 204.001

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./adaptivity_ex5-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:20:44 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           8.342e+00      1.00072   8.337e+00
Objects:              2.510e+03      1.00000   2.510e+03
Flops:                5.241e+07      1.19009   5.005e+07  3.003e+08
Flops/sec:            6.287e+06      1.19008   6.003e+06  3.602e+07
MPI Messages:         7.808e+03      1.08002   7.496e+03  4.497e+04
MPI Message Lengths:  4.387e+06      1.06598   5.628e+02  2.531e+07
MPI Reductions:       6.166e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 8.3374e+00 100.0%  3.0032e+08 100.0%  4.497e+04 100.0%  5.628e+02      100.0%  5.243e+03  85.0% 

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
      %T - percent time in this phase         %F - percent flops in this phase
      %M - percent messages in this phase     %L - percent message lengths in this phase
      %R - percent reductions in this phase
   Total Mflop/s: 10e-6 * (sum of flops over all processors)/(max time over all processors)
------------------------------------------------------------------------------------------------------------------------
Event                Count      Time (sec)     Flops                             --- Global ---  --- Stage ---   Total
                   Max Ratio  Max     Ratio   Max  Ratio  Mess   Avg len Reduct  %T %F %M %L %R  %T %F %M %L %R Mflop/s
------------------------------------------------------------------------------------------------------------------------

--- Event Stage 0: Main Stage

VecMDot              381 1.0 2.3444e-02 1.6 2.41e+06 1.1 0.0e+00 0.0e+00 3.8e+02  0  5  0  0  6   0  5  0  0  7   588
VecNorm              481 1.0 6.8245e-02 1.4 6.71e+05 1.1 0.0e+00 0.0e+00 4.8e+02  1  1  0  0  8   1  1  0  0  9    56
VecScale             431 1.0 3.7217e-04 1.3 3.01e+05 1.1 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  4615
VecCopy              401 1.0 6.4802e-04 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               763 1.0 7.2193e-04 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              100 1.0 1.2794e-0280.3 1.40e+05 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0    62
VecMAXPY             431 1.0 1.4374e-03 1.3 2.95e+06 1.1 0.0e+00 0.0e+00 0.0e+00  0  6  0  0  0   0  6  0  0  0 11709
VecAssemblyBegin     804 1.0 3.2520e-01 1.6 0.00e+00 0.0 2.3e+03 2.9e+02 2.2e+03  3  0  5  3 35   3  0  5  3 42     0
VecAssemblyEnd       804 1.0 7.4077e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     1110 1.0 4.9620e-03 1.3 0.00e+00 0.0 2.5e+04 5.3e+02 0.0e+00  0  0 57 54  0   0  0 57 54  0     0
VecScatterEnd       1110 1.0 4.0789e-01 2.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  3  0  0  0  0   3  0  0  0  0     0
VecNormalize         431 1.0 5.8868e-02 1.6 9.02e+05 1.1 0.0e+00 0.0e+00 4.3e+02  1  2  0  0  7   1  2  0  0  8    88
MatMult              431 1.0 1.5992e-01 1.7 5.32e+06 1.1 1.2e+04 2.5e+02 0.0e+00  2 10 28 12  0   2 10 28 12  0   189
MatSolve             431 1.0 1.8283e-02 1.4 1.88e+07 1.2 0.0e+00 0.0e+00 0.0e+00  0 36  0  0  0   0 36  0  0  0  5931
MatLUFactorNum        50 1.0 3.3426e-02 1.4 2.22e+07 1.3 0.0e+00 0.0e+00 0.0e+00  0 42  0  0  0   0 42  0  0  0  3729
MatILUFactorSym       50 1.0 8.0873e-02 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 5.0e+01  1  0  0  0  1   1  0  0  0  1     0
MatAssemblyBegin     100 1.0 3.2408e-01 1.8 0.00e+00 0.0 2.8e+03 9.6e+02 2.0e+02  3  0  6 11  3   3  0  6 11  4     0
MatAssemblyEnd       100 1.0 2.3991e-02 1.1 0.00e+00 0.0 1.5e+03 6.3e+01 2.6e+02  0  0  3  0  4   0  0  3  0  5     0
MatGetRowIJ           50 1.0 1.7481e-03252.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        50 1.0 2.2762e-03 4.3 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+02  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries       102 1.0 6.3150e-0314.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       381 1.0 2.4964e-02 1.6 4.83e+06 1.1 0.0e+00 0.0e+00 3.8e+02  0  9  0  0  6   0  9  0  0  7  1105
KSPSetup             100 1.0 8.0085e-04 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              50 1.0 3.5368e-01 1.1 5.24e+07 1.2 1.2e+04 2.5e+02 1.0e+03  4100 28 12 16   4100 28 12 19   849
PCSetUp              100 1.0 1.1854e-01 1.4 2.22e+07 1.3 0.0e+00 0.0e+00 1.5e+02  1 42  0  0  2   1 42  0  0  3  1051
PCSetUpOnBlocks       50 1.0 1.1672e-01 1.4 2.22e+07 1.3 0.0e+00 0.0e+00 1.5e+02  1 42  0  0  2   1 42  0  0  3  1068
PCApply              431 1.0 2.2268e-02 1.3 1.88e+07 1.2 0.0e+00 0.0e+00 0.0e+00  0 36  0  0  0   0 36  0  0  0  4869
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec  1132           1132      9953560     0
         Vec Scatter   357            357       309876     0
           Index Set   661            661       911652     0
   IS L to G Mapping   128            128       544008     0
              Matrix   128            128     16063804     0
       Krylov Solver    52             52       490880     0
      Preconditioner    52             52        36608     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 6.38962e-06
Average time for zero size MPI_Send(): 4.06504e-05
#PETSc Option Table entries:
-exact_solution 100+10*sin(pi*x)*cos(pi*y)
-init_timestep 0
-ksp_right_pc
-log_summary
-max_r_steps 1
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
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8
Configure run at: Sat May 19 03:47:23 2012
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared=1 --with-shared-libraries=1 --with-mpi-dir=/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid --with-mumps=true --download-mumps=1 --with-parmetis=true --download-parmetis=1 --with-superlu=true --download-superlu=1 --with-superludir=true --download-superlu_dist=1 --with-blacs=true --download-blacs=1 --with-scalapack=true --download-scalapack=1 --with-hypre=true --download-hypre=1 --with-blas-lib="[/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_intel_lp64.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_sequential.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_core.so]" --with-lapack-lib=/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_solver_lp64_sequential.a
-----------------------------------------
Libraries compiled on Sat May 19 03:47:23 CDT 2012 on daedalus 
Machine characteristics: Linux daedalus 2.6.32-34-generic #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011 x86_64 GNU/Linux 
Using PETSc directory: /org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5
Using PETSc arch: intel-11.1-lucid-mpich2-1.4.1-cxx-opt
-----------------------------------------
Using C compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpicxx -O3   -fPIC   
Using Fortran compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpif90 -fPIC -O3    
-----------------------------------------
Using include paths: -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/include  
------------------------------------------
Using C linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpicxx -O3 
Using Fortran linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpif90 -fPIC -O3  
Using libraries: -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -lpetsc       -lX11 -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -lHYPRE -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lblacs -lsuperlu_dist_2.4 -lparmetis -lmetis -lsuperlu_4.0 -Wl,-rpath,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t -L/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t -lmkl_solver_lp64_sequential -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -ldl -Wl,-rpath,/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/lib -L/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/lib -lmpich -lopa -lmpl -lrt -lpthread -Wl,-rpath,/opt/intel/Compiler/11.1/073/lib/intel64 -L/opt/intel/Compiler/11.1/073/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -L/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -lmpichf90 -lifport -lifcore -lm -lm -lmpichcxx -lstdc++ -lmpichcxx -lstdc++ -ldl -lmpich -lopa -lmpl -lrt -lpthread -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -ldl  
------------------------------------------

-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 6                                                |
| Time:           Fri Aug 24 15:20:44 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=8.50028, Active time=7.57924                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     26        0.0162      0.000623    0.0175      0.000672    0.21     0.23     |
|   build_constraint_matrix()        31536     0.0162      0.000001    0.0162      0.000001    0.21     0.21     |
|   build_sparsity()                 26        0.1225      0.004711    0.1639      0.006304    1.62     2.16     |
|   cnstrn_elem_mat_vec()            31536     0.0211      0.000001    0.0211      0.000001    0.28     0.28     |
|   create_dof_constraints()         26        0.5893      0.022667    1.0752      0.041354    7.78     14.19    |
|   distribute_dofs()                26        0.0476      0.001829    0.3183      0.012244    0.63     4.20     |
|   dof_indices()                    309655    0.1403      0.000000    0.1403      0.000000    1.85     1.85     |
|   enforce_constraints_exactly()    76        0.0309      0.000407    0.0309      0.000407    0.41     0.41     |
|   old_dof_indices()                95814     0.0323      0.000000    0.0323      0.000000    0.43     0.43     |
|   prepare_send_list()              26        0.0002      0.000009    0.0002      0.000009    0.00     0.00     |
|   reinit()                         26        0.1219      0.004689    0.1219      0.004689    1.61     1.61     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          8         0.0065      0.000808    0.0417      0.005215    0.09     0.55     |
|   write()                          1         0.0008      0.000819    0.0009      0.000903    0.01     0.01     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               4         0.0112      0.002804    0.0112      0.002804    0.15     0.15     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        139208    0.1044      0.000001    0.1044      0.000001    1.38     1.38     |
|   init_shape_functions()           74910     0.0895      0.000001    0.0895      0.000001    1.18     1.18     |
|   inverse_map()                    199658    0.2122      0.000001    0.2122      0.000001    2.80     2.80     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             139208    0.1060      0.000001    0.1060      0.000001    1.40     1.40     |
|   compute_face_map()               37404     0.1173      0.000003    0.2183      0.000006    1.55     2.88     |
|   init_face_shape_functions()      6043      0.0055      0.000001    0.0055      0.000001    0.07     0.07     |
|   init_reference_to_physical_map() 74910     0.1242      0.000002    0.1242      0.000002    1.64     1.64     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               4         0.1348      0.033704    0.1348      0.033704    1.78     1.78     |
|                                                                                                                |
| JumpErrorEstimator                                                                                             |
|   estimate_error()                 25        0.4550      0.018201    1.4659      0.058638    6.00     19.34    |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           26280     0.0083      0.000000    0.0083      0.000000    0.11     0.11     |
|   init()                           55        0.0374      0.000680    0.0374      0.000680    0.49     0.49     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       25        0.0055      0.000221    0.0149      0.000594    0.07     0.20     |
|   find_neighbors()                 27        0.1799      0.006663    0.3240      0.011999    2.37     4.27     |
|   renumber_nodes_and_elem()        79        0.0322      0.000408    0.0322      0.000408    0.42     0.42     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0404      0.040403    0.0689      0.068899    0.53     0.91     |
|   compute_hilbert_indices()        28        0.1546      0.005521    0.1546      0.005521    2.04     2.04     |
|   find_global_indices()            28        0.0195      0.000696    0.4210      0.015035    0.26     5.55     |
|   parallel_sort()                  28        0.0144      0.000514    0.1721      0.006147    0.19     2.27     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         8         0.0001      0.000013    0.1879      0.023482    0.00     2.48     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              50        0.0068      0.000136    0.0158      0.000316    0.09     0.21     |
|   _refine_elements()               55        0.0392      0.000713    0.1681      0.003057    0.52     2.22     |
|   add_point()                      26280     0.0179      0.000001    0.0279      0.000001    0.24     0.37     |
|   make_coarsening_compatible()     52        0.2136      0.004108    0.2230      0.004288    2.82     2.94     |
|   make_refinement_compatible()     52        0.0044      0.000085    0.0314      0.000604    0.06     0.41     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0001      0.000071    0.0001      0.000071    0.00     0.00     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      27        0.1728      0.006400    0.6521      0.024151    2.28     8.60     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      141       0.1223      0.000867    0.1223      0.000867    1.61     1.61     |
|   barrier()                        1         0.0000      0.000019    0.0000      0.000019    0.00     0.00     |
|   broadcast()                      2         0.0000      0.000010    0.0000      0.000010    0.00     0.00     |
|   gather()                         20        0.0289      0.001445    0.0289      0.001445    0.38     0.38     |
|   max(bool)                        182       0.1484      0.000816    0.1484      0.000816    1.96     1.96     |
|   max(scalar)                      80        0.5447      0.006808    0.5447      0.006808    7.19     7.19     |
|   max(vector)                      129       0.0075      0.000058    0.0075      0.000058    0.10     0.10     |
|   min(bool)                        104       0.2671      0.002568    0.2671      0.002568    3.52     3.52     |
|   min(scalar)                      25        0.0030      0.000122    0.0030      0.000122    0.04     0.04     |
|   min(vector)                      129       0.1053      0.000817    0.1053      0.000817    1.39     1.39     |
|   probe()                          1175      0.2961      0.000252    0.2961      0.000252    3.91     3.91     |
|   receive()                        1209      0.0030      0.000002    0.2992      0.000247    0.04     3.95     |
|   send()                           1144      0.0017      0.000001    0.0017      0.000001    0.02     0.02     |
|   send_receive()                   1200      0.0025      0.000002    0.3039      0.000253    0.03     4.01     |
|   sum()                            274       0.7495      0.002735    0.7495      0.002735    9.89     9.89     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           1174      0.0005      0.000000    0.0005      0.000000    0.01     0.01     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         27        0.0224      0.000828    0.1960      0.007258    0.29     2.59     |
|   set_parent_processor_ids()       27        0.0163      0.000605    0.0163      0.000605    0.22     0.22     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          50        0.6938      0.013876    0.6938      0.013876    9.15     9.15     |
|                                                                                                                |
| PointLocatorTree                                                                                               |
|   init(no master)                  50        0.2165      0.004330    0.2633      0.005266    2.86     3.47     |
|   operator()                       14033     0.1402      0.000010    0.1635      0.000012    1.85     2.16     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       75        0.0779      0.001038    0.1153      0.001537    1.03     1.52     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       50        0.1819      0.003637    0.3028      0.006056    2.40     4.00     |
|   calculate_norm()                 52        0.1784      0.003431    0.3619      0.006959    2.35     4.77     |
|   project_vector()                 76        0.3124      0.004111    0.4824      0.006348    4.12     6.37     |
|                                                                                                                |
| XdrIO                                                                                                          |
|   write()                          1         0.0059      0.005873    0.0434      0.043400    0.08     0.57     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1214662   7.5792                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***** Finished first 25 steps, now read in saved solution and continue *****
 
Usage:
	 ./adaptivity_ex5-opt -init_timestep 0
OR
	 ./adaptivity_ex5-opt -read_solution -init_timestep 26

Running: ./adaptivity_ex5-opt -read_solution -n_timesteps 25 -max_r_steps 1 -output_freq 10 -init_timestep 25 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=4045
    n_local_nodes()=704
  n_elem()=5172
    n_local_elem()=882
    n_active_elem()=3880
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0

Initial H1 norm = 248.368

 EquationSystems
  n_systems()=1
   System #0, "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=4045
    n_local_dofs()=704
    n_constrained_dofs()=214
    n_local_constrained_dofs()=60
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.61932
      Average Off-Processor Bandwidth <= 0.632102
      Maximum  On-Processor Bandwidth <= 14
      Maximum Off-Processor Bandwidth <= 6
    DofMap Constraints
      Number of DoF Constraints = 214
      Average DoF Constraint Length= 1.43458
      Number of Node Constraints = 308
      Maximum Node Constraint Length= 4
      Average Node Constraint Length= 2.31169

 Solving time step 25, time=0.6500...
H1 norm = 265.97
  Refining the mesh...
H1 norm = 268.126
 Solving time step 26, time=0.6750...
H1 norm = 269.68
  Refining the mesh...
H1 norm = 269.597
 Solving time step 27, time=0.7000...
H1 norm = 271.154
  Refining the mesh...
H1 norm = 270.895
 Solving time step 28, time=0.7250...
H1 norm = 271.987
  Refining the mesh...
H1 norm = 271.857
 Solving time step 29, time=0.7500...
H1 norm = 272.766
  Refining the mesh...
H1 norm = 272.766
 Solving time step 30, time=0.7750...
H1 norm = 273.472
  Refining the mesh...
H1 norm = 273.472
 Solving time step 31, time=0.8000...
H1 norm = 274.022
  Refining the mesh...
H1 norm = 274.022
 Solving time step 32, time=0.8250...
H1 norm = 274.522
  Refining the mesh...
H1 norm = 274.506
 Solving time step 33, time=0.8500...
H1 norm = 274.951
  Refining the mesh...
H1 norm = 274.951
 Solving time step 34, time=0.8750...
H1 norm = 275.309
  Refining the mesh...
H1 norm = 275.309
 Solving time step 35, time=0.9000...
H1 norm = 275.567
  Refining the mesh...
H1 norm = 275.567
 Solving time step 36, time=0.9250...
H1 norm = 275.709
  Refining the mesh...
H1 norm = 275.709
 Solving time step 37, time=0.9500...
H1 norm = 275.723
  Refining the mesh...
H1 norm = 275.723
 Solving time step 38, time=0.9750...
H1 norm = 275.621
  Refining the mesh...
H1 norm = 275.621
 Solving time step 39, time=1.0000...
H1 norm = 275.425
  Refining the mesh...
H1 norm = 275.425
 Solving time step 40, time=1.0300...
H1 norm = 275.161
  Refining the mesh...
H1 norm = 275.147
 Solving time step 41, time=1.0500...
H1 norm = 274.831
  Refining the mesh...
H1 norm = 274.831
 Solving time step 42, time=1.0700...
H1 norm = 274.483
  Refining the mesh...
H1 norm = 274.483
 Solving time step 43, time=1.1000...
H1 norm = 274.104
  Refining the mesh...
H1 norm = 274.104
 Solving time step 44, time=1.1200...
H1 norm = 273.701
  Refining the mesh...
H1 norm = 273.701
 Solving time step 45, time=1.1500...
H1 norm = 273.28
  Refining the mesh...
H1 norm = 273.28
 Solving time step 46, time=1.1700...
H1 norm = 272.85
  Refining the mesh...
H1 norm = 272.85
 Solving time step 47, time=1.2000...
H1 norm = 272.416
  Refining the mesh...
H1 norm = 272.416
 Solving time step 48, time=1.2200...
H1 norm = 271.985
  Refining the mesh...
H1 norm = 271.795
 Solving time step 49, time=1.2500...
H1 norm = 271.37
  Refining the mesh...
H1 norm = 271.357
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./adaptivity_ex5-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:20:46 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.588e+00      1.00036   1.588e+00
Objects:              2.566e+03      1.00156   2.563e+03
Flops:                3.100e+06      1.23857   2.804e+06  1.683e+07
Flops/sec:            1.952e+06      1.23834   1.766e+06  1.060e+07
MPI Messages:         6.440e+03      1.20733   5.812e+03  3.488e+04
MPI Message Lengths:  7.874e+05      1.15781   1.246e+02  4.344e+06
MPI Reductions:       6.177e+03      1.00065

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.5876e+00 100.0%  1.6827e+07 100.0%  3.488e+04 100.0%  1.246e+02      100.0%  5.227e+03  84.7% 

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
      %T - percent time in this phase         %F - percent flops in this phase
      %M - percent messages in this phase     %L - percent message lengths in this phase
      %R - percent reductions in this phase
   Total Mflop/s: 10e-6 * (sum of flops over all processors)/(max time over all processors)
------------------------------------------------------------------------------------------------------------------------
Event                Count      Time (sec)     Flops                             --- Global ---  --- Stage ---   Total
                   Max Ratio  Max     Ratio   Max  Ratio  Mess   Avg len Reduct  %T %F %M %L %R  %T %F %M %L %R Mflop/s
------------------------------------------------------------------------------------------------------------------------

--- Event Stage 0: Main Stage

VecMDot              320 1.0 2.0000e-02 1.8 2.60e+05 1.3 0.0e+00 0.0e+00 3.2e+02  1  8  0  0  5   1  8  0  0  6    71
VecNorm              420 1.0 2.4503e-02 1.3 6.16e+04 1.3 0.0e+00 0.0e+00 4.2e+02  1  2  0  0  7   1  2  0  0  8    14
VecScale             370 1.0 2.0051e-04 1.5 2.75e+04 1.3 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0   754
VecCopy              391 1.0 1.5330e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               677 1.0 2.0266e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               83 1.0 1.3304e-04 1.7 1.16e+04 1.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   476
VecMAXPY             353 1.0 2.2149e-04 1.3 3.10e+05 1.3 0.0e+00 0.0e+00 0.0e+00  0 10  0  0  0   0 10  0  0  0  7698
VecAssemblyBegin     827 1.0 1.2848e-01 1.3 0.00e+00 0.0 1.9e+03 8.0e+01 2.2e+03  7  0  5  3 36   7  0  5  3 43     0
VecAssemblyEnd       827 1.0 5.8746e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     1069 1.0 2.3715e-03 1.3 0.00e+00 0.0 2.0e+04 1.2e+02 0.0e+00  0  0 56 54  0   0  0 56 54  0     0
VecScatterEnd       1069 1.0 6.6670e-02 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  3  0  0  0  0   3  0  0  0  0     0
VecNormalize         370 1.0 2.3312e-02 1.4 8.26e+04 1.3 0.0e+00 0.0e+00 3.7e+02  1  3  0  0  6   1  3  0  0  7    19
MatMult              370 1.0 2.8821e-02 2.1 5.37e+05 1.3 9.9e+03 8.7e+01 0.0e+00  1 17 28 20  0   1 17 28 20  0   102
MatSolve             353 1.0 1.5254e-03 1.5 1.20e+06 1.2 0.0e+00 0.0e+00 0.0e+00  0 38  0  0  0   0 38  0  0  0  4238
MatLUFactorNum        50 1.0 1.6768e-03 1.3 7.23e+05 1.4 0.0e+00 0.0e+00 0.0e+00  0 22  0  0  0   0 22  0  0  0  2237
MatILUFactorSym       50 1.0 5.6727e-03 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 5.0e+01  0  0  0  0  1   0  0  0  0  1     0
MatAssemblyBegin     100 1.0 4.0241e-02 1.9 0.00e+00 0.0 2.4e+03 3.3e+02 2.0e+02  2  0  7 18  3   2  0  7 18  4     0
MatAssemblyEnd       100 1.0 1.8474e-02 1.2 0.00e+00 0.0 1.4e+03 2.4e+01 2.6e+02  1  0  4  1  4   1  0  4  1  5     0
MatGetRowIJ           50 1.0 4.1008e-05 3.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        50 1.0 3.5882e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+02  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries       104 1.0 1.0633e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       320 1.0 2.0382e-02 1.8 5.22e+05 1.3 0.0e+00 0.0e+00 3.2e+02  1 17  0  0  5   1 17  0  0  6   141
KSPSetup             100 1.0 3.7551e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              50 1.0 7.4887e-02 1.1 3.10e+06 1.2 9.9e+03 8.7e+01 8.9e+02  5100 28 20 14   5100 28 20 17   225
PCSetUp              100 1.0 9.5568e-03 1.3 7.23e+05 1.4 0.0e+00 0.0e+00 1.5e+02  1 22  0  0  2   1 22  0  0  3   392
PCSetUpOnBlocks       50 1.0 8.3933e-03 1.4 7.23e+05 1.4 0.0e+00 0.0e+00 1.5e+02  0 22  0  0  2   0 22  0  0  3   447
PCApply              353 1.0 3.8767e-03 1.3 1.20e+06 1.2 0.0e+00 0.0e+00 0.0e+00  0 38  0  0  0   0 38  0  0  0  1668
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec  1153           1153      2573856     0
         Vec Scatter   366            366       317688     0
           Index Set   675            675       410624     0
   IS L to G Mapping   133            133       123776     0
              Matrix   131            131      1561720     0
       Krylov Solver    52             52       490880     0
      Preconditioner    52             52        36608     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 2.19345e-06
Average time for zero size MPI_Send(): 2.60274e-05
#PETSc Option Table entries:
-init_timestep 25
-ksp_right_pc
-log_summary
-max_r_steps 1
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
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8
Configure run at: Sat May 19 03:47:23 2012
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared=1 --with-shared-libraries=1 --with-mpi-dir=/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid --with-mumps=true --download-mumps=1 --with-parmetis=true --download-parmetis=1 --with-superlu=true --download-superlu=1 --with-superludir=true --download-superlu_dist=1 --with-blacs=true --download-blacs=1 --with-scalapack=true --download-scalapack=1 --with-hypre=true --download-hypre=1 --with-blas-lib="[/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_intel_lp64.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_sequential.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_core.so]" --with-lapack-lib=/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_solver_lp64_sequential.a
-----------------------------------------
Libraries compiled on Sat May 19 03:47:23 CDT 2012 on daedalus 
Machine characteristics: Linux daedalus 2.6.32-34-generic #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011 x86_64 GNU/Linux 
Using PETSc directory: /org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5
Using PETSc arch: intel-11.1-lucid-mpich2-1.4.1-cxx-opt
-----------------------------------------
Using C compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpicxx -O3   -fPIC   
Using Fortran compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpif90 -fPIC -O3    
-----------------------------------------
Using include paths: -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/include  
------------------------------------------
Using C linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpicxx -O3 
Using Fortran linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpif90 -fPIC -O3  
Using libraries: -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -lpetsc       -lX11 -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -lHYPRE -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lblacs -lsuperlu_dist_2.4 -lparmetis -lmetis -lsuperlu_4.0 -Wl,-rpath,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t -L/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t -lmkl_solver_lp64_sequential -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -ldl -Wl,-rpath,/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/lib -L/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/lib -lmpich -lopa -lmpl -lrt -lpthread -Wl,-rpath,/opt/intel/Compiler/11.1/073/lib/intel64 -L/opt/intel/Compiler/11.1/073/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -L/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -lmpichf90 -lifport -lifcore -lm -lm -lmpichcxx -lstdc++ -lmpichcxx -lstdc++ -ldl -lmpich -lopa -lmpl -lrt -lpthread -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -ldl  
------------------------------------------

-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 6                                                |
| Time:           Fri Aug 24 15:20:46 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=1.68061, Active time=1.43536                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     27        0.0040      0.000146    0.0042      0.000156    0.28     0.29     |
|   build_constraint_matrix()        2375      0.0029      0.000001    0.0029      0.000001    0.20     0.20     |
|   build_sparsity()                 27        0.0191      0.000708    0.0237      0.000877    1.33     1.65     |
|   cnstrn_elem_mat_vec()            2375      0.0092      0.000004    0.0092      0.000004    0.64     0.64     |
|   create_dof_constraints()         27        0.1098      0.004065    0.2190      0.008112    7.65     15.26    |
|   distribute_dofs()                27        0.0117      0.000432    0.0575      0.002131    0.81     4.01     |
|   dof_indices()                    41708     0.0166      0.000000    0.0166      0.000000    1.16     1.16     |
|   enforce_constraints_exactly()    78        0.0200      0.000257    0.0200      0.000257    1.40     1.40     |
|   old_dof_indices()                10851     0.0037      0.000000    0.0037      0.000000    0.26     0.26     |
|   prepare_send_list()              27        0.0001      0.000002    0.0001      0.000002    0.00     0.00     |
|   reinit()                         27        0.0208      0.000769    0.0208      0.000769    1.45     1.45     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          8         0.0020      0.000252    0.0138      0.001721    0.14     0.96     |
|   read()                           1         0.0160      0.015992    0.1477      0.147738    1.11     10.29    |
|   update()                         1         0.0006      0.000580    0.0006      0.000580    0.04     0.04     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               4         0.0043      0.001087    0.0043      0.001087    0.30     0.30     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        15904     0.0111      0.000001    0.0111      0.000001    0.78     0.78     |
|   init_shape_functions()           10609     0.0098      0.000001    0.0098      0.000001    0.68     0.68     |
|   inverse_map()                    41729     0.0307      0.000001    0.0307      0.000001    2.14     2.14     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             15904     0.0115      0.000001    0.0115      0.000001    0.80     0.80     |
|   compute_face_map()               5254      0.0133      0.000003    0.0235      0.000004    0.92     1.64     |
|   init_face_shape_functions()      2233      0.0022      0.000001    0.0022      0.000001    0.15     0.15     |
|   init_reference_to_physical_map() 10609     0.0277      0.000003    0.0277      0.000003    1.93     1.93     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               4         0.0703      0.017573    0.0703      0.017573    4.90     4.90     |
|                                                                                                                |
| JumpErrorEstimator                                                                                             |
|   estimate_error()                 25        0.0615      0.002458    0.1403      0.005610    4.28     9.77     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           48        0.0000      0.000000    0.0000      0.000000    0.00     0.00     |
|   init()                           51        0.0033      0.000066    0.0033      0.000066    0.23     0.23     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       26        0.0017      0.000066    0.0036      0.000137    0.12     0.25     |
|   find_neighbors()                 9         0.0159      0.001769    0.0242      0.002693    1.11     1.69     |
|   renumber_nodes_and_elem()        44        0.0035      0.000080    0.0035      0.000080    0.25     0.25     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0350      0.035013    0.0740      0.073959    2.44     5.15     |
|   compute_hilbert_indices()        9         0.0106      0.001177    0.0106      0.001177    0.74     0.74     |
|   find_global_indices()            9         0.0020      0.000225    0.0379      0.004215    0.14     2.64     |
|   parallel_sort()                  9         0.0137      0.001522    0.0190      0.002106    0.95     1.32     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         8         0.0001      0.000018    0.0886      0.011069    0.01     6.17     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              51        0.0029      0.000058    0.0051      0.000101    0.20     0.36     |
|   _refine_elements()               51        0.0020      0.000040    0.0109      0.000214    0.14     0.76     |
|   add_point()                      48        0.0000      0.000001    0.0001      0.000001    0.00     0.00     |
|   make_coarsening_compatible()     51        0.0432      0.000847    0.0613      0.001203    3.01     4.27     |
|   make_refinement_compatible()     51        0.0004      0.000007    0.0075      0.000146    0.03     0.52     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      9         0.0140      0.001554    0.0579      0.006429    0.97     4.03     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      107       0.0572      0.000535    0.0572      0.000535    3.99     3.99     |
|   broadcast()                      61        0.0089      0.000147    0.0089      0.000146    0.62     0.62     |
|   gather()                         1         0.0005      0.000507    0.0005      0.000507    0.04     0.04     |
|   max(bool)                        179       0.0166      0.000093    0.0166      0.000093    1.16     1.16     |
|   max(scalar)                      61        0.1488      0.002440    0.1488      0.002440    10.37    10.37    |
|   max(vector)                      78        0.0065      0.000083    0.0065      0.000083    0.45     0.45     |
|   min(bool)                        102       0.0308      0.000301    0.0308      0.000301    2.14     2.14     |
|   min(scalar)                      25        0.0010      0.000039    0.0010      0.000039    0.07     0.07     |
|   min(vector)                      78        0.0121      0.000155    0.0121      0.000155    0.84     0.84     |
|   probe()                          750       0.0839      0.000112    0.0839      0.000112    5.85     5.85     |
|   receive()                        750       0.0013      0.000002    0.0853      0.000114    0.09     5.94     |
|   send()                           750       0.0006      0.000001    0.0006      0.000001    0.04     0.04     |
|   send_receive()                   772       0.0016      0.000002    0.0879      0.000114    0.11     6.13     |
|   sum()                            212       0.0583      0.000275    0.0583      0.000275    4.06     4.06     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           750       0.0003      0.000000    0.0003      0.000000    0.02     0.02     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         10        0.0027      0.000268    0.0632      0.006321    0.19     4.40     |
|   set_parent_processor_ids()       9         0.0014      0.000154    0.0014      0.000154    0.10     0.10     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          50        0.1389      0.002778    0.1389      0.002778    9.68     9.68     |
|                                                                                                                |
| PointLocatorTree                                                                                               |
|   init(no master)                  34        0.0195      0.000572    0.0308      0.000905    1.36     2.14     |
|   operator()                       7296      0.0499      0.000007    0.0578      0.000008    3.48     4.02     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       78        0.0083      0.000106    0.0118      0.000151    0.58     0.82     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       50        0.0159      0.000318    0.0351      0.000702    1.11     2.44     |
|   calculate_norm()                 51        0.0302      0.000592    0.0472      0.000925    2.10     3.29     |
|   project_vector()                 78        0.1091      0.001399    0.1465      0.001878    7.60     10.21    |
|                                                                                                                |
| XdrIO                                                                                                          |
|   read()                           1         0.0038      0.003820    0.0126      0.012564    0.27     0.88     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            172642    1.4354                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 6 ./adaptivity_ex5-opt [-read_solution] -n_timesteps 25 -init_timestep [0|25] -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
