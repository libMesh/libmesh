<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex24",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 24 - Periodic Boundary Conditions with Adaptive Mesh Refinement</h1>

<br><br>This example uses the same simple, linear transient
system as in example 10; however in this case periodic boundary
conditions are applied at the sides of the domain.
 

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
              mesh.read("saved_mesh.xda");
        
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
              equation_systems.read("saved_solution.xda", libMeshEnums::READ);
        
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
            GMVIO(mesh).write_equation_systems ("out.gmv.000",
                                                equation_systems);
            ExodusII_IO(mesh).write_equation_systems ("out.e.000",
                                                equation_systems);
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
            GMVIO(mesh).write_equation_systems ("solution_read_in.gmv",
                                                equation_systems);
            ExodusII_IO(mesh).write_equation_systems ("solution_read_in.e",
                                                equation_systems);
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
</div>

<div class ="fragment">
<pre>
          const Real dt = 0.025;
          Real time     = init_timestep*dt;
        
</pre>
</div>
<div class = "comment">
Creating and attaching Periodic Boundaries
</div>

<div class ="fragment">
<pre>
          TransientLinearImplicitSystem& system =
        	equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;
                  ("Convection-Diffusion");
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
              time += dt;
        
              equation_systems.parameters.set&lt;Real&gt; ("time") = time;
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
                OSSRealzeroleft(out,6,3,time);
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
        
                  file_name &lt;&lt; "out.gmv.";
                  OSSRealzeroright(file_name,3,0,t_step+1);
        
                  GMVIO(mesh).write_equation_systems (file_name.str(),
                                                      equation_systems);
                  exodus_file_name &lt;&lt; "out.e.";
                  OSSRealzeroright(exodus_file_name,3,0,t_step+1);
                  ExodusII_IO(mesh).write_equation_systems (exodus_file_name.str(),
                                                      equation_systems);
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
        
              mesh.write("saved_mesh.xda");
              equation_systems.write("saved_solution.xda", libMeshEnums::WRITE);
              GMVIO(mesh).write_equation_systems ("saved_solution.gmv",
                                                  equation_systems);
              ExodusII_IO(mesh).write_equation_systems ("saved_solution.e",
                                                  equation_systems);
            }
        #endif // #ifndef LIBMESH_ENABLE_AMR
          
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
          es.parameters.set&lt;Real&gt; ("time") = 0;
          
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
                          const Number value = exact_solution (qface_points[qp](0),
                                                               qface_points[qp](1),
                                                               time);
                                                               
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
right-hand-side vector.  The \p PetscMatrix::add_matrix()
and \p PetscVector::add_vector() members do this for us.
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
  
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_ENABLE_AMR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-amr&quot;</FONT></B>);
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
        mesh.read(<B><FONT COLOR="#BC8F8F">&quot;saved_mesh.xda&quot;</FONT></B>);
  
        mesh.print_info();
  
        equation_systems.read(<B><FONT COLOR="#BC8F8F">&quot;saved_solution.xda&quot;</FONT></B>, libMeshEnums::READ);
  
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
      GMVIO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;out.gmv.000&quot;</FONT></B>,
                                          equation_systems);
      ExodusII_IO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;out.e.000&quot;</FONT></B>,
                                          equation_systems);
    }
    <B><FONT COLOR="#A020F0">else</FONT></B>
    {
      GMVIO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;solution_read_in.gmv&quot;</FONT></B>,
                                          equation_systems);
      ExodusII_IO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;solution_read_in.e&quot;</FONT></B>,
                                          equation_systems);
    }
    
  
    equation_systems.parameters.set&lt;RealVectorValue&gt;(<B><FONT COLOR="#BC8F8F">&quot;velocity&quot;</FONT></B>) = 
      RealVectorValue (0.8, 0.8);
  
    equation_systems.parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;diffusivity&quot;</FONT></B>) = 0.01;
      
    <B><FONT COLOR="#228B22">const</FONT></B> Real dt = 0.025;
    Real time     = init_timestep*dt;
  
    TransientLinearImplicitSystem&amp; system =
  	equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;
            (<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
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
        time += dt;
  
        equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>) = time;
        equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;dt&quot;</FONT></B>)   = dt;
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Solving time step &quot;</FONT></B>;
        
        {
          OStringStream out;
  
          OSSInt(out,2,t_step);
          out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, time=&quot;</FONT></B>;
          OSSRealzeroleft(out,6,3,time);
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
  
            file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out.gmv.&quot;</FONT></B>;
            OSSRealzeroright(file_name,3,0,t_step+1);
  
            GMVIO(mesh).write_equation_systems (file_name.str(),
                                                equation_systems);
            exodus_file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out.e.&quot;</FONT></B>;
            OSSRealzeroright(exodus_file_name,3,0,t_step+1);
            ExodusII_IO(mesh).write_equation_systems (exodus_file_name.str(),
                                                equation_systems);
          }
      }
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(!read_solution)
      {
        TransientLinearImplicitSystem&amp; system =
  	equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;
            (<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
        Real H1norm = system.calculate_norm(*system.solution, SystemNorm(H1));
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Final H1 norm = &quot;</FONT></B> &lt;&lt; H1norm &lt;&lt; std::endl &lt;&lt; std::endl;
  
        mesh.write(<B><FONT COLOR="#BC8F8F">&quot;saved_mesh.xda&quot;</FONT></B>);
        equation_systems.write(<B><FONT COLOR="#BC8F8F">&quot;saved_solution.xda&quot;</FONT></B>, libMeshEnums::WRITE);
        GMVIO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;saved_solution.gmv&quot;</FONT></B>,
                                            equation_systems);
        ExodusII_IO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;saved_solution.e&quot;</FONT></B>,
                                            equation_systems);
      }
  #endif <I><FONT COLOR="#B22222">// #ifndef LIBMESH_ENABLE_AMR
</FONT></I>    
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> init_cd (EquationSystems&amp; es,
                <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    TransientLinearImplicitSystem &amp; system =
      es.get_system&lt;TransientLinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    es.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>) = 0;
    
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
                    <B><FONT COLOR="#228B22">const</FONT></B> Number value = exact_solution (qface_points[qp](0),
                                                         qface_points[qp](1),
                                                         time);
                                                         
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
***************************************************************
* Running Example  mpirun -np 2 ./ex24-opt [-read_solution] -n_timesteps 25 -n_refinements 5 -init_timestep [0|25] -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Usage:
	 ./ex24-opt -init_timestep 0
OR
	 ./ex24-opt -read_solution -init_timestep 26

Running: ./ex24-opt -n_timesteps 25 -n_refinements 5 -output_freq 10 -init_timestep 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=4225
    n_local_nodes()=2154
  n_elem()=5460
    n_local_elem()=2763
    n_active_elem()=4096
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=4225
    n_local_dofs()=2154
    n_constrained_dofs()=0
    n_vectors()=3

 Solving time step  0, time=0.0250...
H1 norm = 1.58843
  Refining the mesh...
H1 norm = 1.59126
 Solving time step  1, time=0.0500...
H1 norm = 1.46335
  Refining the mesh...
H1 norm = 1.46259
 Solving time step  2, time=0.0750...
H1 norm = 1.35362
  Refining the mesh...
H1 norm = 1.3532
 Solving time step  3, time=0.1000...
H1 norm = 1.26706
  Refining the mesh...
H1 norm = 1.25834
 Solving time step  4, time=0.1250...
H1 norm = 1.19156
  Refining the mesh...
H1 norm = 1.1753
 Solving time step  5, time=0.1500...
H1 norm = 1.11401
  Refining the mesh...
H1 norm = 1.10286
 Solving time step  6, time=0.1750...
H1 norm = 1.0402
  Refining the mesh...
H1 norm = 1.03894
 Solving time step  7, time=0.2000...
H1 norm = 0.980851
  Refining the mesh...
H1 norm = 0.981935
 Solving time step  8, time=0.2250...
H1 norm = 0.930851
  Refining the mesh...
H1 norm = 0.93087
 Solving time step  9, time=0.2500...
H1 norm = 0.884844
  Refining the mesh...
H1 norm = 0.884818
 Solving time step 10, time=0.2750...
H1 norm = 0.843725
  Refining the mesh...
H1 norm = 0.843011
 Solving time step 11, time=0.3000...
H1 norm = 0.805383
  Refining the mesh...
H1 norm = 0.805019
 Solving time step 12, time=0.3250...
H1 norm = 0.770451
  Refining the mesh...
H1 norm = 0.770366
 Solving time step 13, time=0.3500...
H1 norm = 0.739221
  Refining the mesh...
H1 norm = 0.739206
 Solving time step 14, time=0.3750...
H1 norm = 0.709996
  Refining the mesh...
H1 norm = 0.70976
 Solving time step 15, time=0.4000...
H1 norm = 0.682995
  Refining the mesh...
H1 norm = 0.682983
 Solving time step 16, time=0.4250...
H1 norm = 0.657861
  Refining the mesh...
H1 norm = 0.657899
 Solving time step 17, time=0.4500...
H1 norm = 0.633984
  Refining the mesh...
H1 norm = 0.634595
 Solving time step 18, time=0.4750...
H1 norm = 0.616379
  Refining the mesh...
H1 norm = 0.616357
 Solving time step 19, time=0.5000...
H1 norm = 0.594936
  Refining the mesh...
H1 norm = 0.594251
 Solving time step 20, time=0.5250...
H1 norm = 0.575096
  Refining the mesh...
H1 norm = 0.575097
 Solving time step 21, time=0.5500...
H1 norm = 0.557235
  Refining the mesh...
H1 norm = 0.557228
 Solving time step 22, time=0.5750...
H1 norm = 0.540147
  Refining the mesh...
H1 norm = 0.540141
 Solving time step 23, time=0.6000...
H1 norm = 0.523817
  Refining the mesh...
H1 norm = 0.52381
 Solving time step 24, time=0.6250...
H1 norm = 0.509187
  Refining the mesh...
H1 norm = 0.509185
Final H1 norm = 0.509185

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex24-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Thu Feb  3 12:13:03 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.521e+00      1.00039   1.521e+00
Objects:              2.586e+03      1.00000   2.586e+03
Flops:                2.113e+07      1.20970   1.930e+07  3.860e+07
Flops/sec:            1.390e+07      1.20923   1.269e+07  2.539e+07
MPI Messages:         2.284e+03      1.00000   2.284e+03  4.569e+03
MPI Message Lengths:  1.787e+06      1.04285   7.662e+02  3.501e+06
MPI Reductions:       5.694e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.5205e+00 100.0%  3.8602e+07 100.0%  4.569e+03 100.0%  7.662e+02      100.0%  4.773e+03  83.8% 

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

VecMDot               76 1.0 6.0773e-04 1.5 9.25e+04 1.1 0.0e+00 0.0e+00 7.6e+01  0  0  0  0  1   0  0  0  0  2   288
VecNorm              176 1.0 9.5630e-04 1.1 1.22e+05 1.1 0.0e+00 0.0e+00 1.8e+02  0  1  0  0  3   0  1  0  0  4   242
VecScale             126 1.0 7.5102e-05 1.2 4.38e+04 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1105
VecCopy              240 1.0 1.1444e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               438 1.0 1.4901e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               90 1.0 3.7289e-03 1.0 6.25e+04 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0    32
VecMAXPY             116 1.0 6.3181e-05 1.0 1.46e+05 1.1 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  4364
VecAssemblyBegin     702 1.0 6.2032e-03 1.3 0.00e+00 0.0 2.4e+02 3.3e+02 1.9e+03  0  0  5  2 33   0  0  5  2 39     0
VecAssemblyEnd       702 1.0 3.1137e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      879 1.0 1.1129e-03 1.1 0.00e+00 0.0 1.2e+03 7.7e+02 0.0e+00  0  0 27 28  0   0  0 27 28  0     0
VecScatterEnd        879 1.0 2.6488e-02 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
VecNormalize         126 1.0 9.4104e-04 1.1 1.31e+05 1.1 0.0e+00 0.0e+00 1.3e+02  0  1  0  0  2   0  1  0  0  3   265
MatMult              126 1.0 2.6956e-02 1.8 8.09e+05 1.1 2.5e+02 5.3e+02 0.0e+00  1  4  6  4  0   1  4  6  4  0    57
MatSolve             116 1.0 2.8610e-03 1.3 3.47e+06 1.2 0.0e+00 0.0e+00 0.0e+00  0 17  0  0  0   0 17  0  0  0  2237
MatLUFactorNum        50 1.0 1.9330e-02 1.2 1.64e+07 1.2 0.0e+00 0.0e+00 0.0e+00  1 77  0  0  0   1 77  0  0  0  1541
MatILUFactorSym       50 1.0 5.3654e-02 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 5.0e+01  3  0  0  0  1   3  0  0  0  1     0
MatAssemblyBegin     100 1.0 1.1640e-02 1.4 0.00e+00 0.0 2.8e+02 1.6e+03 2.0e+02  1  0  6 13  4   1  0  6 13  4     0
MatAssemblyEnd       100 1.0 3.6428e-03 1.1 0.00e+00 0.0 1.0e+02 1.3e+02 2.6e+02  0  0  2  0  4   0  0  2  0  5     0
MatGetRowIJ           50 1.0 8.8215e-06 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        50 1.0 3.2640e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+02  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries       102 1.0 1.7047e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        76 1.0 7.0858e-04 1.4 1.85e+05 1.1 0.0e+00 0.0e+00 7.6e+01  0  1  0  0  1   0  1  0  0  2   495
KSPSetup             100 1.0 3.1948e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              50 1.0 1.0038e-01 1.0 2.11e+07 1.2 2.5e+02 5.3e+02 4.0e+02  7100  6  4  7   7100  6  4  8   385
PCSetUp              100 1.0 7.5412e-02 1.2 1.64e+07 1.2 0.0e+00 0.0e+00 1.5e+02  5 77  0  0  3   5 77  0  0  3   395
PCSetUpOnBlocks       50 1.0 7.4395e-02 1.2 1.64e+07 1.2 0.0e+00 0.0e+00 1.5e+02  5 77  0  0  3   5 77  0  0  3   400
PCApply              116 1.0 3.5706e-03 1.2 3.47e+06 1.2 0.0e+00 0.0e+00 0.0e+00  0 17  0  0  0   0 17  0  0  0  1792
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec   910            910      3524176     0
         Vec Scatter   506            506       439208     0
           Index Set   810            810       680180     0
   IS L to G Mapping   128            128       292700     0
              Matrix   128            128     10858296     0
       Krylov Solver    52             52       490880     0
      Preconditioner    52             52        36608     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 6.19888e-07
Average time for zero size MPI_Send(): 6.07967e-06
#PETSc Option Table entries:
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
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8
Configure run at: Fri Oct 15 13:01:23 2010
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared=1 --with-mpi-dir=/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid --with-mumps=true --download-mumps=ifneeded --with-parmetis=true --download-parmetis=ifneeded --with-superlu=true --download-superlu=ifneeded --with-superludir=true --download-superlu_dist=ifneeded --with-blacs=true --download-blacs=ifneeded --with-scalapack=true --download-scalapack=ifneeded --with-hypre=true --download-hypre=ifneeded --with-blas-lib="[/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_intel_lp64.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_sequential.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_core.so]" --with-lapack-lib=/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_solver_lp64_sequential.a
-----------------------------------------
Libraries compiled on Fri Oct 15 13:01:23 CDT 2010 on atreides 
Machine characteristics: Linux atreides 2.6.32-25-generic #44-Ubuntu SMP Fri Sep 17 20:05:27 UTC 2010 x86_64 GNU/Linux 
Using PETSc directory: /org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5
Using PETSc arch: gcc-4.5-lucid-mpich2-1.2.1-cxx-opt
-----------------------------------------
Using C compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpicxx -Wall -Wwrite-strings -Wno-strict-aliasing -O3   -fPIC   
Using Fortran compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpif90 -fPIC -Wall -Wno-unused-variable -O3    
-----------------------------------------
Using include paths: -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/include  
------------------------------------------
Using C linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpicxx -Wall -Wwrite-strings -Wno-strict-aliasing -O3 
Using Fortran linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpif90 -fPIC -Wall -Wno-unused-variable -O3  
Using libraries: -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -lpetsc       -lX11 -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -lHYPRE -lsuperlu_dist_2.4 -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -lscalapack -lblacs -lsuperlu_4.0 -Wl,-rpath,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t -L/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t -lmkl_solver_lp64_sequential -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -Wl,-rpath,/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/lib -L/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/lib -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib/gcc/x86_64-unknown-linux-gnu/4.5.1 -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib/gcc/x86_64-unknown-linux-gnu/4.5.1 -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib64 -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib64 -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib -ldl -lmpich -lopa -lpthread -lrt -lgcc_s -lmpichf90 -lgfortran -lm -lm -lmpichcxx -lstdc++ -ldl -lmpich -lopa -lpthread -lrt -lgcc_s -ldl  
------------------------------------------

-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 2                                                |
| Time:           Thu Feb  3 12:13:03 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=1.52856, Active time=1.42636                                                |
 -------------------------------------------------------------------------------------------------------------
| Event                           nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                           w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------|
|                                                                                                             |
|                                                                                                             |
| DofMap                                                                                                      |
|   add_neighbors_to_send_list()  26        0.0053      0.000203    0.0059      0.000227    0.37     0.41     |
|   build_constraint_matrix()     12628     0.0101      0.000001    0.0101      0.000001    0.71     0.71     |
|   cnstrn_elem_mat_vec()         12628     0.0095      0.000001    0.0095      0.000001    0.67     0.67     |
|   compute_sparsity()            26        0.0254      0.000976    0.0316      0.001217    1.78     2.22     |
|   create_dof_constraints()      26        0.1235      0.004749    0.1657      0.006373    8.66     11.62    |
|   distribute_dofs()             26        0.0154      0.000594    0.0463      0.001782    1.08     3.25     |
|   dof_indices()                 119298    0.0412      0.000000    0.0412      0.000000    2.89     2.89     |
|   enforce_constraints_exactly() 75        0.0028      0.000038    0.0028      0.000038    0.20     0.20     |
|   old_dof_indices()             43074     0.0126      0.000000    0.0126      0.000000    0.88     0.88     |
|   prepare_send_list()           26        0.0001      0.000005    0.0001      0.000005    0.01     0.01     |
|   reinit()                      26        0.0299      0.001150    0.0299      0.001150    2.10     2.10     |
|                                                                                                             |
| EquationSystems                                                                                             |
|   write()                       1         0.0038      0.003843    0.0039      0.003863    0.27     0.27     |
|                                                                                                             |
| FE                                                                                                          |
|   compute_affine_map()          66018     0.0328      0.000000    0.0328      0.000000    2.30     2.30     |
|   compute_face_map()            18837     0.0387      0.000002    0.0757      0.000004    2.71     5.31     |
|   compute_shape_functions()     66018     0.0156      0.000000    0.0156      0.000000    1.10     1.10     |
|   init_face_shape_functions()   17636     0.0083      0.000000    0.0083      0.000000    0.58     0.58     |
|   init_shape_functions()        36435     0.0715      0.000002    0.0715      0.000002    5.01     5.01     |
|   inverse_map()                 196686    0.1668      0.000001    0.1668      0.000001    11.69    11.69    |
|                                                                                                             |
| GMVIO                                                                                                       |
|   write_nodal_data()            4         0.0170      0.004256    0.0170      0.004256    1.19     1.19     |
|                                                                                                             |
| JumpErrorEstimator                                                                                          |
|   estimate_error()              25        0.1676      0.006704    0.4171      0.016684    11.75    29.24    |
|                                                                                                             |
| LocationMap                                                                                                 |
|   find()                        19116     0.0064      0.000000    0.0064      0.000000    0.45     0.45     |
|   init()                        55        0.0065      0.000119    0.0065      0.000119    0.46     0.46     |
|                                                                                                             |
| Mesh                                                                                                        |
|   contract()                    25        0.0036      0.000144    0.0053      0.000210    0.25     0.37     |
|   find_neighbors()              27        0.0464      0.001718    0.0465      0.001723    3.25     3.26     |
|   renumber_nodes_and_elem()     79        0.0037      0.000047    0.0037      0.000047    0.26     0.26     |
|                                                                                                             |
| MeshCommunication                                                                                           |
|   assign_global_indices()       1         0.0068      0.006775    0.0068      0.006843    0.47     0.48     |
|   compute_hilbert_indices()     28        0.0275      0.000984    0.0275      0.000984    1.93     1.93     |
|   find_global_indices()         28        0.0053      0.000188    0.0352      0.001257    0.37     2.47     |
|   parallel_sort()               28        0.0016      0.000058    0.0019      0.000067    0.11     0.13     |
|                                                                                                             |
| MeshRefinement                                                                                              |
|   _coarsen_elements()           50        0.0031      0.000063    0.0032      0.000065    0.22     0.23     |
|   _refine_elements()            55        0.0198      0.000359    0.0411      0.000748    1.39     2.88     |
|   add_point()                   19116     0.0121      0.000001    0.0198      0.000001    0.85     1.39     |
|   make_coarsening_compatible()  63        0.1351      0.002145    0.1893      0.003004    9.47     13.27    |
|   make_refinement_compatible()  63        0.0021      0.000034    0.0025      0.000040    0.15     0.18     |
|                                                                                                             |
| MeshTools::Generation                                                                                       |
|   build_cube()                  1         0.0000      0.000047    0.0000      0.000047    0.00     0.00     |
|                                                                                                             |
| MetisPartitioner                                                                                            |
|   partition()                   27        0.0236      0.000873    0.0585      0.002168    1.65     4.10     |
|                                                                                                             |
| Parallel                                                                                                    |
|   allgather()                   112       0.0004      0.000003    0.0004      0.000003    0.03     0.03     |
|   barrier()                     1         0.0000      0.000008    0.0000      0.000008    0.00     0.00     |
|   broadcast()                   2         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   gather()                      20        0.0001      0.000003    0.0001      0.000003    0.00     0.00     |
|   max(bool)                     193       0.0006      0.000003    0.0006      0.000003    0.04     0.04     |
|   max(scalar)                   156       0.0007      0.000005    0.0007      0.000005    0.05     0.05     |
|   max(vector)                   129       0.0002      0.000002    0.0002      0.000002    0.02     0.02     |
|   min(bool)                     126       0.0029      0.000023    0.0029      0.000023    0.20     0.20     |
|   min(scalar)                   25        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   min(vector)                   129       0.0006      0.000005    0.0006      0.000005    0.04     0.04     |
|   probe()                       235       0.0004      0.000002    0.0004      0.000002    0.03     0.03     |
|   receive()                     245       0.0004      0.000002    0.0008      0.000003    0.03     0.06     |
|   send()                        232       0.0002      0.000001    0.0002      0.000001    0.01     0.01     |
|   send_receive()                288       0.0005      0.000002    0.0016      0.000006    0.04     0.12     |
|   sum()                         196       0.0079      0.000040    0.0079      0.000040    0.55     0.55     |
|   wait()                        228       0.0001      0.000000    0.0001      0.000000    0.01     0.01     |
|                                                                                                             |
| Partitioner                                                                                                 |
|   set_node_processor_ids()      27        0.0080      0.000296    0.0084      0.000312    0.56     0.59     |
|   set_parent_processor_ids()    27        0.0021      0.000079    0.0021      0.000079    0.15     0.15     |
|                                                                                                             |
| PetscLinearSolver                                                                                           |
|   solve()                       50        0.1145      0.002291    0.1145      0.002291    8.03     8.03     |
|                                                                                                             |
| ProjectVector                                                                                               |
|   operator()                    75        0.0282      0.000376    0.0431      0.000575    1.98     3.02     |
|                                                                                                             |
| System                                                                                                      |
|   assemble()                    50        0.0770      0.001540    0.1321      0.002643    5.40     9.26     |
|   calculate_norm()              51        0.0359      0.000705    0.0561      0.001099    2.52     3.93     |
|   project_vector()              76        0.0408      0.000537    0.0963      0.001267    2.86     6.75     |
|                                                                                                             |
| XdrIO                                                                                                       |
|   write()                       1         0.0030      0.002963    0.0031      0.003127    0.21     0.22     |
 -------------------------------------------------------------------------------------------------------------
| Totals:                         630955    1.4264                                          100.00            |
 -------------------------------------------------------------------------------------------------------------

 
***** Finished first 25 steps, now read in saved solution and continue *****
 
Usage:
	 ./ex24-opt -init_timestep 0
OR
	 ./ex24-opt -read_solution -init_timestep 26

Running: ./ex24-opt -read_solution -n_timesteps 25 -output_freq 10 -init_timestep 25 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=737
    n_local_nodes()=384
  n_elem()=884
    n_local_elem()=456
    n_active_elem()=664
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

Initial H1 norm = 0.509185

 EquationSystems
  n_systems()=1
   System "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=737
    n_local_dofs()=384
    n_constrained_dofs()=110
    n_vectors()=3

 Solving time step 25, time=0.6500...
H1 norm = 0.495054
  Refining the mesh...
H1 norm = 0.495
 Solving time step 26, time=0.6750...
H1 norm = 0.481405
  Refining the mesh...
H1 norm = 0.481395
 Solving time step 27, time=0.7000...
H1 norm = 0.468861
  Refining the mesh...
H1 norm = 0.468833
 Solving time step 28, time=0.7250...
H1 norm = 0.456975
  Refining the mesh...
H1 norm = 0.45693
 Solving time step 29, time=0.7500...
H1 norm = 0.44567
  Refining the mesh...
H1 norm = 0.44567
 Solving time step 30, time=0.7750...
H1 norm = 0.434922
  Refining the mesh...
H1 norm = 0.434915
 Solving time step 31, time=0.8000...
H1 norm = 0.424677
  Refining the mesh...
H1 norm = 0.424677
 Solving time step 32, time=0.8250...
H1 norm = 0.414928
  Refining the mesh...
H1 norm = 0.414923
 Solving time step 33, time=0.8500...
H1 norm = 0.405638
  Refining the mesh...
H1 norm = 0.40561
 Solving time step 34, time=0.8750...
H1 norm = 0.396739
  Refining the mesh...
H1 norm = 0.396711
 Solving time step 35, time=0.9000...
H1 norm = 0.388224
  Refining the mesh...
H1 norm = 0.388223
 Solving time step 36, time=0.9250...
H1 norm = 0.38009
  Refining the mesh...
H1 norm = 0.380065
 Solving time step 37, time=0.9500...
H1 norm = 0.372274
  Refining the mesh...
H1 norm = 0.372272
 Solving time step 38, time=0.9750...
H1 norm = 0.364797
  Refining the mesh...
H1 norm = 0.364779
 Solving time step 39, time=1.0000...
H1 norm = 0.357609
  Refining the mesh...
H1 norm = 0.357605
 Solving time step 40, time=1.0300...
H1 norm = 0.3507
  Refining the mesh...
H1 norm = 0.350699
 Solving time step 41, time=1.0500...
H1 norm = 0.344075
  Refining the mesh...
H1 norm = 0.344072
 Solving time step 42, time=1.0700...
H1 norm = 0.337687
  Refining the mesh...
H1 norm = 0.337686
 Solving time step 43, time=1.1000...
H1 norm = 0.331543
  Refining the mesh...
H1 norm = 0.331524
 Solving time step 44, time=1.1200...
H1 norm = 0.325566
  Refining the mesh...
H1 norm = 0.325598
 Solving time step 45, time=1.1500...
H1 norm = 0.319891
  Refining the mesh...
H1 norm = 0.319889
 Solving time step 46, time=1.1700...
H1 norm = 0.31438
  Refining the mesh...
H1 norm = 0.31438
 Solving time step 47, time=1.2000...
H1 norm = 0.309061
  Refining the mesh...
H1 norm = 0.309059
 Solving time step 48, time=1.2200...
H1 norm = 0.303919
  Refining the mesh...
H1 norm = 0.303918
 Solving time step 49, time=1.2500...
H1 norm = 0.298949
  Refining the mesh...
H1 norm = 0.298936
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex24-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Thu Feb  3 12:13:05 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.781e+00      1.00049   1.780e+00
Objects:              2.798e+03      1.00000   2.798e+03
Flops:                4.180e+07      1.20899   3.819e+07  7.638e+07
Flops/sec:            2.348e+07      1.20840   2.145e+07  4.291e+07
MPI Messages:         2.398e+03      1.00000   2.398e+03  4.796e+03
MPI Message Lengths:  2.321e+06      1.04365   9.475e+02  4.544e+06
MPI Reductions:       5.882e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.7801e+00 100.0%  7.6381e+07 100.0%  4.796e+03 100.0%  9.475e+02      100.0%  4.961e+03  84.3% 

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

VecMDot              167 1.0 1.4305e-03 1.7 5.28e+05 1.1 0.0e+00 0.0e+00 1.7e+02  0  1  0  0  3   0  1  0  0  3   703
VecNorm              267 1.0 1.6289e-03 1.2 2.67e+05 1.1 0.0e+00 0.0e+00 2.7e+02  0  1  0  0  5   0  1  0  0  5   311
VecScale             217 1.0 1.2207e-04 1.2 1.08e+05 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1691
VecCopy              235 1.0 1.3423e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               519 1.0 1.8454e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               85 1.0 8.7976e-05 1.0 8.39e+04 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1814
VecMAXPY             202 1.0 2.1577e-04 1.1 6.96e+05 1.1 0.0e+00 0.0e+00 0.0e+00  0  2  0  0  0   0  2  0  0  0  6139
VecAssemblyBegin     703 1.0 5.9037e-03 1.1 0.00e+00 0.0 2.4e+02 3.8e+02 1.9e+03  0  0  5  2 32   0  0  5  2 38     0
VecAssemblyEnd       703 1.0 3.2592e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      972 1.0 1.3740e-03 1.1 0.00e+00 0.0 1.4e+03 9.5e+02 0.0e+00  0  0 30 30  0   0  0 30 30  0     0
VecScatterEnd        972 1.0 2.2917e-02 5.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
VecNormalize         217 1.0 1.6897e-03 1.2 3.25e+05 1.1 0.0e+00 0.0e+00 2.2e+02  0  1  0  0  4   0  1  0  0  4   366
MatMult              217 1.0 2.4532e-02 4.0 1.98e+06 1.1 4.3e+02 6.0e+02 0.0e+00  1  5  9  6  0   1  5  9  6  0   153
MatSolve             202 1.0 6.1259e-03 1.1 9.41e+06 1.1 0.0e+00 0.0e+00 0.0e+00  0 23  0  0  0   0 23  0  0  0  2886
MatLUFactorNum        50 1.0 3.2010e-02 1.2 2.87e+07 1.2 0.0e+00 0.0e+00 0.0e+00  2 68  0  0  0   2 68  0  0  0  1616
MatILUFactorSym       50 1.0 8.4386e-02 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 5.0e+01  4  0  0  0  1   4  0  0  0  1     0
MatAssemblyBegin     100 1.0 4.6849e-03 1.4 0.00e+00 0.0 2.8e+02 1.7e+03 2.0e+02  0  0  6 11  3   0  0  6 11  4     0
MatAssemblyEnd       100 1.0 4.2849e-03 1.1 0.00e+00 0.0 1.0e+02 1.5e+02 2.6e+02  0  0  2  0  4   0  0  2  0  5     0
MatGetRowIJ           50 1.0 6.9141e-06 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        50 1.0 3.6716e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+02  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries       102 1.0 1.9073e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       167 1.0 1.6911e-03 1.5 1.06e+06 1.1 0.0e+00 0.0e+00 1.7e+02  0  3  0  0  3   0  3  0  0  3  1189
KSPSetup             100 1.0 3.2210e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              50 1.0 1.3661e-01 1.0 4.18e+07 1.2 4.3e+02 6.0e+02 5.8e+02  8100  9  6 10   8100  9  6 12   559
PCSetUp              100 1.0 1.1884e-01 1.2 2.87e+07 1.2 0.0e+00 0.0e+00 1.5e+02  6 68  0  0  3   6 68  0  0  3   435
PCSetUpOnBlocks       50 1.0 1.1785e-01 1.2 2.87e+07 1.2 0.0e+00 0.0e+00 1.5e+02  6 68  0  0  3   6 68  0  0  3   439
PCApply              202 1.0 7.3452e-03 1.1 9.41e+06 1.1 0.0e+00 0.0e+00 0.0e+00  0 23  0  0  0   0 23  0  0  0  2407
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec  1120           1120      5412528     0
         Vec Scatter   507            507       440076     0
           Index Set   811            811       786016     0
   IS L to G Mapping   128            128       364140     0
              Matrix   128            128     16219232     0
       Krylov Solver    52             52       490880     0
      Preconditioner    52             52        36608     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 6.19888e-07
Average time for zero size MPI_Send(): 5.96046e-06
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
sizeof(short) 2 sizeof(int) 4 sizeof(long) 8 sizeof(void*) 8 sizeof(PetscScalar) 8
Configure run at: Fri Oct 15 13:01:23 2010
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared=1 --with-mpi-dir=/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid --with-mumps=true --download-mumps=ifneeded --with-parmetis=true --download-parmetis=ifneeded --with-superlu=true --download-superlu=ifneeded --with-superludir=true --download-superlu_dist=ifneeded --with-blacs=true --download-blacs=ifneeded --with-scalapack=true --download-scalapack=ifneeded --with-hypre=true --download-hypre=ifneeded --with-blas-lib="[/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_intel_lp64.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_sequential.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_core.so]" --with-lapack-lib=/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_solver_lp64_sequential.a
-----------------------------------------
Libraries compiled on Fri Oct 15 13:01:23 CDT 2010 on atreides 
Machine characteristics: Linux atreides 2.6.32-25-generic #44-Ubuntu SMP Fri Sep 17 20:05:27 UTC 2010 x86_64 GNU/Linux 
Using PETSc directory: /org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5
Using PETSc arch: gcc-4.5-lucid-mpich2-1.2.1-cxx-opt
-----------------------------------------
Using C compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpicxx -Wall -Wwrite-strings -Wno-strict-aliasing -O3   -fPIC   
Using Fortran compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpif90 -fPIC -Wall -Wno-unused-variable -O3    
-----------------------------------------
Using include paths: -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/include  
------------------------------------------
Using C linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpicxx -Wall -Wwrite-strings -Wno-strict-aliasing -O3 
Using Fortran linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpif90 -fPIC -Wall -Wno-unused-variable -O3  
Using libraries: -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -lpetsc       -lX11 -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -lHYPRE -lsuperlu_dist_2.4 -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -lscalapack -lblacs -lsuperlu_4.0 -Wl,-rpath,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t -L/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t -lmkl_solver_lp64_sequential -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -Wl,-rpath,/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/lib -L/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/lib -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib/gcc/x86_64-unknown-linux-gnu/4.5.1 -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib/gcc/x86_64-unknown-linux-gnu/4.5.1 -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib64 -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib64 -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib -ldl -lmpich -lopa -lpthread -lrt -lgcc_s -lmpichf90 -lgfortran -lm -lm -lmpichcxx -lstdc++ -ldl -lmpich -lopa -lpthread -lrt -lgcc_s -ldl  
------------------------------------------

-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 2                                                |
| Time:           Thu Feb  3 12:13:05 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=1.78902, Active time=1.68098                                                |
 -------------------------------------------------------------------------------------------------------------
| Event                           nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                           w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------|
|                                                                                                             |
|                                                                                                             |
| DofMap                                                                                                      |
|   add_neighbors_to_send_list()  26        0.0066      0.000253    0.0075      0.000287    0.39     0.44     |
|   build_constraint_matrix()     21576     0.0142      0.000001    0.0142      0.000001    0.85     0.85     |
|   cnstrn_elem_mat_vec()         21576     0.0150      0.000001    0.0150      0.000001    0.89     0.89     |
|   compute_sparsity()            26        0.0345      0.001328    0.0425      0.001635    2.05     2.53     |
|   create_dof_constraints()      26        0.1578      0.006068    0.1873      0.007204    9.38     11.14    |
|   distribute_dofs()             26        0.0173      0.000666    0.0539      0.002075    1.03     3.21     |
|   dof_indices()                 162543    0.0541      0.000000    0.0541      0.000000    3.22     3.22     |
|   enforce_constraints_exactly() 75        0.0031      0.000042    0.0031      0.000042    0.19     0.19     |
|   old_dof_indices()             65970     0.0186      0.000000    0.0186      0.000000    1.11     1.11     |
|   prepare_send_list()           26        0.0002      0.000006    0.0002      0.000006    0.01     0.01     |
|   reinit()                      26        0.0355      0.001367    0.0355      0.001367    2.11     2.11     |
|                                                                                                             |
| EquationSystems                                                                                             |
|   read()                        1         0.0053      0.005329    0.0163      0.016303    0.32     0.97     |
|   update()                      1         0.0002      0.000206    0.0002      0.000206    0.01     0.01     |
|                                                                                                             |
| FE                                                                                                          |
|   compute_affine_map()          89420     0.0443      0.000000    0.0443      0.000000    2.64     2.64     |
|   compute_face_map()            23369     0.0479      0.000002    0.0940      0.000004    2.85     5.59     |
|   compute_shape_functions()     89420     0.0219      0.000000    0.0219      0.000000    1.30     1.30     |
|   init_face_shape_functions()   22436     0.0107      0.000000    0.0107      0.000000    0.64     0.64     |
|   init_shape_functions()        46037     0.0897      0.000002    0.0897      0.000002    5.34     5.34     |
|   inverse_map()                 166620    0.1417      0.000001    0.1417      0.000001    8.43     8.43     |
|                                                                                                             |
| GMVIO                                                                                                       |
|   write_nodal_data()            4         0.0100      0.002499    0.0100      0.002499    0.59     0.59     |
|                                                                                                             |
| JumpErrorEstimator                                                                                          |
|   estimate_error()              25        0.2172      0.008688    0.5527      0.022106    12.92    32.88    |
|                                                                                                             |
| LocationMap                                                                                                 |
|   find()                        3672      0.0013      0.000000    0.0013      0.000000    0.07     0.07     |
|   init()                        50        0.0086      0.000171    0.0086      0.000171    0.51     0.51     |
|                                                                                                             |
| Mesh                                                                                                        |
|   contract()                    25        0.0024      0.000096    0.0038      0.000153    0.14     0.23     |
|   find_neighbors()              26        0.0499      0.001918    0.0504      0.001938    2.97     3.00     |
|   renumber_nodes_and_elem()     77        0.0036      0.000047    0.0036      0.000047    0.21     0.21     |
|                                                                                                             |
| MeshCommunication                                                                                           |
|   assign_global_indices()       1         0.0070      0.006968    0.0071      0.007057    0.41     0.42     |
|   compute_hilbert_indices()     26        0.0358      0.001379    0.0358      0.001379    2.13     2.13     |
|   find_global_indices()         26        0.0063      0.000241    0.0446      0.001714    0.37     2.65     |
|   parallel_sort()               26        0.0017      0.000066    0.0019      0.000075    0.10     0.12     |
|                                                                                                             |
| MeshRefinement                                                                                              |
|   _coarsen_elements()           50        0.0029      0.000058    0.0030      0.000060    0.17     0.18     |
|   _refine_elements()            50        0.0096      0.000192    0.0140      0.000280    0.57     0.83     |
|   add_point()                   3672      0.0024      0.000001    0.0039      0.000001    0.14     0.23     |
|   make_coarsening_compatible()  58        0.1436      0.002476    0.1505      0.002595    8.54     8.96     |
|   make_refinement_compatible()  58        0.0025      0.000043    0.0029      0.000050    0.15     0.17     |
|                                                                                                             |
| MetisPartitioner                                                                                            |
|   partition()                   26        0.0285      0.001098    0.0733      0.002819    1.70     4.36     |
|                                                                                                             |
| Parallel                                                                                                    |
|   allgather()                   111       0.0005      0.000005    0.0005      0.000005    0.03     0.03     |
|   broadcast()                   50        0.0001      0.000002    0.0001      0.000001    0.01     0.00     |
|   gather()                      1         0.0000      0.000005    0.0000      0.000005    0.00     0.00     |
|   max(bool)                     183       0.0006      0.000003    0.0006      0.000003    0.04     0.04     |
|   max(scalar)                   152       0.0044      0.000029    0.0044      0.000029    0.26     0.26     |
|   max(vector)                   127       0.0003      0.000002    0.0003      0.000002    0.02     0.02     |
|   min(bool)                     116       0.0015      0.000013    0.0015      0.000013    0.09     0.09     |
|   min(scalar)                   25        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   min(vector)                   127       0.0005      0.000004    0.0005      0.000004    0.03     0.03     |
|   probe()                       214       0.0003      0.000002    0.0003      0.000002    0.02     0.02     |
|   receive()                     214       0.0004      0.000002    0.0008      0.000004    0.03     0.05     |
|   send()                        214       0.0002      0.000001    0.0002      0.000001    0.01     0.01     |
|   send_receive()                270       0.0005      0.000002    0.0016      0.000006    0.03     0.09     |
|   sum()                         184       0.0113      0.000062    0.0113      0.000062    0.67     0.67     |
|   wait()                        214       0.0001      0.000000    0.0001      0.000000    0.00     0.00     |
|                                                                                                             |
| Partitioner                                                                                                 |
|   set_node_processor_ids()      27        0.0091      0.000336    0.0095      0.000353    0.54     0.57     |
|   set_parent_processor_ids()    26        0.0027      0.000104    0.0027      0.000104    0.16     0.16     |
|                                                                                                             |
| PetscLinearSolver                                                                                           |
|   solve()                       50        0.1468      0.002936    0.1468      0.002936    8.73     8.73     |
|                                                                                                             |
| ProjectVector                                                                                               |
|   operator()                    75        0.0471      0.000629    0.0712      0.000950    2.80     4.24     |
|                                                                                                             |
| System                                                                                                      |
|   assemble()                    50        0.1091      0.002183    0.1801      0.003601    6.49     10.71    |
|   calculate_norm()              51        0.0505      0.000990    0.0786      0.001541    3.00     4.68     |
|   project_vector()              75        0.0404      0.000539    0.1267      0.001689    2.40     7.54     |
|                                                                                                             |
| XdrIO                                                                                                       |
|   read()                        1         0.0024      0.002407    0.0024      0.002440    0.14     0.15     |
 -------------------------------------------------------------------------------------------------------------
| Totals:                         719629    1.6810                                          100.00            |
 -------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 2 ./ex24-opt [-read_solution] -n_timesteps 25 -init_timestep [0|25] -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
