<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("adaptivity_ex2",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file adaptivity_ex2.C with comments: </h1> 
<div class = "comment">
<h1>Adaptivity Example 2 - Solving a Transient System with Adaptive Mesh Refinement</h1>

<br><br>This example shows how a simple, linear transient
system can be solved in parallel.  The system is simple
scalar convection-diffusion with a specified external
velocity.  The initial condition is given, and the
solution is advanced in time with a standard Crank-Nicolson
time-stepping strategy.

<br><br>Also, we use this example to demonstrate writing out and reading
in of solutions. We do 25 time steps, then save the solution
and do another 25 time steps starting from the saved solution.
 

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
        #include "libmesh/equation_systems.h"
        #include "libmesh/fe.h"
        #include "libmesh/quadrature_gauss.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/dense_vector.h"
        
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
        
        #ifndef LIBMESH_ENABLE_AMR
          libmesh_example_assert(false, "--enable-amr");
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
We still need some work on automatic parallel restarts with
ParallelMesh
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
</pre>
</div>
<div class = "comment">
Read the mesh from file.
</div>

<div class ="fragment">
<pre>
              mesh.read ("mesh.xda");
        
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
</pre>
</div>
<div class = "comment">
Write out the initial condition
</div>

<div class ="fragment">
<pre>
            GMVIO(mesh).write_equation_systems ("out.gmv.000",
                                                equation_systems);
          else
</pre>
</div>
<div class = "comment">
Write out the solution that was read in
</div>

<div class ="fragment">
<pre>
            GMVIO(mesh).write_equation_systems ("solution_read_in.gmv",
                                                equation_systems);
        
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


<br><br>Since only \p TransientLinearImplicitSystems (and systems
derived from them) contain old solutions, to use the
old_local_solution later we now need to specify the system
type when we ask for it.
</div>

<div class ="fragment">
<pre>
          TransientLinearImplicitSystem &  system =
            equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;("Convection-Diffusion");
        
          const Real dt = 0.025;
          system.time   = init_timestep*dt;
          
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
Add a set of scope braces to enforce data locality.
</div>

<div class ="fragment">
<pre>
              {
                std::ostringstream out;
        
                out &lt;&lt; std::setw(2)
                    &lt;&lt; std::right
                    &lt;&lt; t_step
                    &lt;&lt; ", time="
                    &lt;&lt; std::fixed
                    &lt;&lt; std::setw(6)
                    &lt;&lt; std::setprecision(3)
                    &lt;&lt; std::setfill('0')
                    &lt;&lt; std::left
                    &lt;&lt; system.time
                    &lt;&lt;  "...";
        
                std::cout &lt;&lt; out.str() &lt;&lt; std::endl;
              }
              
</pre>
</div>
<div class = "comment">
At this point we need to update the old
solution vector.  The old solution vector
will be the current solution vector from the
previous time step.


<br><br></div>

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
                  std::ostringstream file_name;
        
                  file_name &lt;&lt; "out.gmv."
                            &lt;&lt; std::setw(3)
                            &lt;&lt; std::setfill('0')
                            &lt;&lt; std::right
                            &lt;&lt; t_step+1;
        
                  GMVIO(mesh).write_equation_systems (file_name.str(),
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
                                                               system.time);
                                                               
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
<br><br><br> <h1> The source file adaptivity_ex2.C without comments: </h1> 
<pre> 
   
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;cstdlib&gt; <I><FONT COLOR="#B22222">// *must* precede &lt;cmath&gt; for proper std:abs() on PGI, Sun Studio CC
</FONT></I>  #include &lt;cmath&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/serial_mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/gmv_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  
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
  
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #ifndef LIBMESH_ENABLE_AMR
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-amr&quot;</FONT></B>);
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
  
    SerialMesh mesh;
  
    EquationSystems equation_systems (mesh);
    MeshRefinement mesh_refinement (mesh);
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(!read_solution)
      {
        mesh.read (<B><FONT COLOR="#BC8F8F">&quot;mesh.xda&quot;</FONT></B>);
  
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
      GMVIO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;out.gmv.000&quot;</FONT></B>,
                                          equation_systems);
    <B><FONT COLOR="#A020F0">else</FONT></B>
      GMVIO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;solution_read_in.gmv&quot;</FONT></B>,
                                          equation_systems);
  
    equation_systems.parameters.set&lt;RealVectorValue&gt;(<B><FONT COLOR="#BC8F8F">&quot;velocity&quot;</FONT></B>) = 
      RealVectorValue (0.8, 0.8);
  
    equation_systems.parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;diffusivity&quot;</FONT></B>) = 0.01;
      
  
    TransientLinearImplicitSystem &amp;  system =
      equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real dt = 0.025;
    system.time   = init_timestep*dt;
    
    <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> t_step=init_timestep; 
                     t_step&lt;(init_timestep+n_timesteps); 
                     t_step++)
      {
        system.time += dt;
  
        equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>) = system.time;
        equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;dt&quot;</FONT></B>)   = dt;
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Solving time step &quot;</FONT></B>;
        
        {
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream out;
  
          out &lt;&lt; std::setw(2)
              &lt;&lt; std::right
              &lt;&lt; t_step
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, time=&quot;</FONT></B>
              &lt;&lt; std::fixed
              &lt;&lt; std::setw(6)
              &lt;&lt; std::setprecision(3)
              &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
              &lt;&lt; std::left
              &lt;&lt; system.time
              &lt;&lt;  <B><FONT COLOR="#BC8F8F">&quot;...&quot;</FONT></B>;
  
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; out.str() &lt;&lt; std::endl;
        }
        
  
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
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream file_name;
  
            file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out.gmv.&quot;</FONT></B>
                      &lt;&lt; std::setw(3)
                      &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
                      &lt;&lt; std::right
                      &lt;&lt; t_step+1;
  
            GMVIO(mesh).write_equation_systems (file_name.str(),
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
      }
  #endif <I><FONT COLOR="#B22222">// #ifndef LIBMESH_ENABLE_AMR
</FONT></I>    
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> init_cd (EquationSystems&amp; es,
                <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    TransientLinearImplicitSystem &amp; system =
      es.get_system&lt;TransientLinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    es.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>) = system.time = 0;
    
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
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
    
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; psi = fe_face-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; qface_points = fe_face-&gt;get_xyz();
      
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
                                                         system.time);
                                                         
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
* Running Example adaptivity_ex2:
*  mpirun -np 12 example-devel -n_timesteps 25 -n_refinements 5 -output_freq 10 -init_timestep 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Usage:
	 /workspace/libmesh/examples/adaptivity/adaptivity_ex2/.libs/lt-example-devel -init_timestep 0
OR
	 /workspace/libmesh/examples/adaptivity/adaptivity_ex2/.libs/lt-example-devel -read_solution -init_timestep 26

Running: /workspace/libmesh/examples/adaptivity/adaptivity_ex2/.libs/lt-example-devel -n_timesteps 25 -n_refinements 5 -output_freq 10 -init_timestep 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=6273
    n_local_nodes()=491
  n_elem()=13650
    n_local_elem()=1229
    n_active_elem()=10240
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
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
    n_dofs()=6273
    n_local_dofs()=491
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.45943
      Average Off-Processor Bandwidth <= 0.323609
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 7
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

 Solving time step  0, time=0.0250...
H1 norm = 1.58843
  Refining the mesh...
H1 norm = 1.58839
 Solving time step  1, time=0.0500...
H1 norm = 1.46061
  Refining the mesh...
H1 norm = 1.45992
 Solving time step  2, time=0.0750...
H1 norm = 1.35107
  Refining the mesh...
H1 norm = 1.35069
 Solving time step  3, time=0.1000...
H1 norm = 1.25698
  Refining the mesh...
H1 norm = 1.25635
 Solving time step  4, time=0.1250...
H1 norm = 1.17458
  Refining the mesh...
H1 norm = 1.1744
 Solving time step  5, time=0.1500...
H1 norm = 1.10264
  Refining the mesh...
H1 norm = 1.10224
 Solving time step  6, time=0.1750...
H1 norm = 1.03868
  Refining the mesh...
H1 norm = 1.03853
 Solving time step  7, time=0.2000...
H1 norm = 0.981978
  Refining the mesh...
H1 norm = 0.981584
 Solving time step  8, time=0.2250...
H1 norm = 0.930848
  Refining the mesh...
H1 norm = 0.930668
 Solving time step  9, time=0.2500...
H1 norm = 0.884889
  Refining the mesh...
H1 norm = 0.884744
 Solving time step 10, time=0.2750...
H1 norm = 0.84331
  Refining the mesh...
H1 norm = 0.843132
 Solving time step 11, time=0.3000...
H1 norm = 0.805394
  Refining the mesh...
H1 norm = 0.805264
 Solving time step 12, time=0.3250...
H1 norm = 0.770793
  Refining the mesh...
H1 norm = 0.7707
 Solving time step 13, time=0.3500...
H1 norm = 0.739086
  Refining the mesh...
H1 norm = 0.738904
 Solving time step 14, time=0.3750...
H1 norm = 0.709764
  Refining the mesh...
H1 norm = 0.709687
 Solving time step 15, time=0.4000...
H1 norm = 0.682786
  Refining the mesh...
H1 norm = 0.68269
 Solving time step 16, time=0.4250...
H1 norm = 0.657767
  Refining the mesh...
H1 norm = 0.657729
 Solving time step 17, time=0.4500...
H1 norm = 0.634598
  Refining the mesh...
H1 norm = 0.634502
 Solving time step 18, time=0.4750...
H1 norm = 0.61295
  Refining the mesh...
H1 norm = 0.612906
 Solving time step 19, time=0.5000...
H1 norm = 0.592784
  Refining the mesh...
H1 norm = 0.592699
 Solving time step 20, time=0.5250...
H1 norm = 0.573854
  Refining the mesh...
H1 norm = 0.573819
 Solving time step 21, time=0.5500...
H1 norm = 0.556154
  Refining the mesh...
H1 norm = 0.55611
 Solving time step 22, time=0.5750...
H1 norm = 0.539527
  Refining the mesh...
H1 norm = 0.539462
 Solving time step 23, time=0.6000...
H1 norm = 0.523854
  Refining the mesh...
H1 norm = 0.523811
 Solving time step 24, time=0.6250...
H1 norm = 0.509112
  Refining the mesh...
H1 norm = 0.509057
Final H1 norm = 0.509057

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/adaptivity/adaptivity_ex2/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 21:58:27 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           1.139e+01      1.00002   1.139e+01
Objects:              2.555e+03      1.00000   2.555e+03
Flops:                2.616e+06      1.83413   2.020e+06  2.424e+07
Flops/sec:            2.296e+05      1.83413   1.773e+05  2.128e+06
MPI Messages:         8.698e+03      1.26656   7.864e+03  9.437e+04
MPI Message Lengths:  7.931e+05      1.19158   9.377e+01  8.849e+06
MPI Reductions:       5.898e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.1392e+01 100.0%  2.4239e+07 100.0%  9.437e+04 100.0%  9.377e+01      100.0%  5.897e+03 100.0% 

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

VecMDot              301 1.0 3.5844e-03 1.5 1.29e+05 1.6 0.0e+00 0.0e+00 3.0e+02  0  5  0  0  5   0  5  0  0  5   341
VecNorm              401 1.0 1.4075e-02 2.8 5.01e+04 1.6 0.0e+00 0.0e+00 4.0e+02  0  2  0  0  7   0  2  0  0  7    35
VecScale             351 1.0 2.5249e-04 1.7 2.16e+04 1.6 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0   830
VecCopy              301 1.0 2.5249e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               885 1.0 5.7530e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              100 1.0 1.6033e-02 2.3 1.37e+04 1.5 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0     9
VecMAXPY             351 1.0 1.9455e-04 1.3 1.66e+05 1.6 0.0e+00 0.0e+00 0.0e+00  0  7  0  0  0   0  7  0  0  0  8159
VecAssemblyBegin     802 1.0 1.1708e-01 3.1 0.00e+00 0.0 4.1e+03 7.9e+01 2.2e+03  1  0  4  4 37   1  0  4  4 37     0
VecAssemblyEnd       802 1.0 7.5579e-04 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     1028 1.0 4.1199e-03 1.1 0.00e+00 0.0 5.0e+04 1.0e+02 0.0e+00  0  0 53 57  0   0  0 53 57  0     0
VecScatterEnd       1028 1.0 1.0186e-02 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize         351 1.0 1.3615e-02 3.0 6.49e+04 1.6 0.0e+00 0.0e+00 3.5e+02  0  3  0  0  6   0  3  0  0  6    46
MatMult              351 1.0 7.2086e-03 1.4 3.66e+05 1.6 1.8e+04 5.7e+01 0.0e+00  0 15 20 12  0   0 15 20 12  0   492
MatSolve             401 1.0 1.6134e-03 1.7 1.03e+06 1.9 0.0e+00 0.0e+00 0.0e+00  0 39  0  0  0   0 39  0  0  0  5894
MatLUFactorNum        50 1.0 2.6040e-03 1.8 8.36e+05 2.0 0.0e+00 0.0e+00 0.0e+00  0 31  0  0  0   0 31  0  0  0  2895
MatILUFactorSym       50 1.0 5.7592e-03 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 1.5e+02  0  0  0  0  3   0  0  0  0  3     0
MatAssemblyBegin     100 1.0 1.1373e-01 1.8 0.00e+00 0.0 4.3e+03 2.8e+02 2.0e+02  1  0  5 13  3   1  0  5 13  3     0
MatAssemblyEnd       100 1.0 5.6927e-03 1.2 0.00e+00 0.0 2.7e+03 1.7e+01 2.1e+02  0  0  3  1  4   0  0  3  1  4     0
MatGetRowIJ           50 1.0 3.2663e-05 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        50 1.0 1.4513e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+02  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries       102 1.0 1.4257e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       301 1.0 3.9859e-03 1.4 2.59e+05 1.6 0.0e+00 0.0e+00 3.0e+02  0 10  0  0  5   0 10  0  0  5   618
KSPSetUp             100 1.0 1.0448e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              50 1.0 5.9328e-02 1.0 2.62e+06 1.8 1.8e+04 5.7e+01 1.0e+03  1100 20 12 17   1100 20 12 17   409
PCSetUp              100 1.0 1.9998e-02 1.2 8.36e+05 2.0 0.0e+00 0.0e+00 3.0e+02  0 31  0  0  5   0 31  0  0  5   377
PCSetUpOnBlocks       50 1.0 1.6297e-02 1.2 8.36e+05 2.0 0.0e+00 0.0e+00 2.5e+02  0 31  0  0  4   0 31  0  0  4   463
PCApply              401 1.0 5.9686e-03 1.1 1.03e+06 1.9 0.0e+00 0.0e+00 0.0e+00  0 39  0  0  0   0 39  0  0  0  1593
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector  1172           1172      2904992     0
      Vector Scatter   357            357       369852     0
           Index Set   663            663       532112     0
   IS L to G Mapping   130            130        73320     0
              Matrix   128            128      1492604     0
       Krylov Solver    52             52       503360     0
      Preconditioner    52             52        46384     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 5.62191e-05
Average time for zero size MPI_Send(): 2.87493e-05
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
| Time:           Thu Jan 31 21:58:27 2013                                                                             |
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
| libMesh Performance: Alive time=11.6414, Active time=10.9636                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     26        0.1305      0.005021    0.2120      0.008152    1.19     1.93     |
|   build_constraint_matrix()        2442      0.0228      0.000009    0.0228      0.000009    0.21     0.21     |
|   build_sparsity()                 26        0.0698      0.002686    0.2173      0.008357    0.64     1.98     |
|   cnstrn_elem_mat_vec()            2442      0.0120      0.000005    0.0120      0.000005    0.11     0.11     |
|   create_dof_constraints()         26        0.3193      0.012282    0.6211      0.023889    2.91     5.67     |
|   distribute_dofs()                26        0.5046      0.019409    1.5935      0.061290    4.60     14.53    |
|   dof_indices()                    31729     1.3934      0.000044    1.3934      0.000044    12.71    12.71    |
|   enforce_constraints_exactly()    75        0.0118      0.000157    0.0118      0.000157    0.11     0.11     |
|   old_dof_indices()                10194     0.4962      0.000049    0.4962      0.000049    4.53     4.53     |
|   prepare_send_list()              26        0.0011      0.000042    0.0011      0.000042    0.01     0.01     |
|   reinit()                         26        1.0172      0.039125    1.0172      0.039125    9.28     9.28     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          4         0.0159      0.003980    0.0645      0.016115    0.15     0.59     |
|   write()                          1         0.0061      0.006062    0.0069      0.006854    0.06     0.06     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        14056     0.1245      0.000009    0.1245      0.000009    1.14     1.14     |
|   init_shape_functions()           7429      0.0775      0.000010    0.0775      0.000010    0.71     0.71     |
|   inverse_map()                    25243     0.1654      0.000007    0.1654      0.000007    1.51     1.51     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             14056     0.1229      0.000009    0.1229      0.000009    1.12     1.12     |
|   compute_face_map()               3723      0.0621      0.000017    0.1229      0.000033    0.57     1.12     |
|   init_face_shape_functions()      63        0.0007      0.000012    0.0007      0.000012    0.01     0.01     |
|   init_reference_to_physical_map() 7429      0.0775      0.000010    0.0775      0.000010    0.71     0.71     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               4         0.0795      0.019875    0.0800      0.019994    0.73     0.73     |
|                                                                                                                |
| JumpErrorEstimator                                                                                             |
|   estimate_error()                 25        0.2806      0.011224    1.2793      0.051171    2.56     11.67    |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           35664     0.1346      0.000004    0.1346      0.000004    1.23     1.23     |
|   init()                           55        0.0623      0.001133    0.0623      0.001133    0.57     0.57     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       25        0.0567      0.002267    0.1007      0.004029    0.52     0.92     |
|   find_neighbors()                 27        1.0252      0.037970    1.0394      0.038496    9.35     9.48     |
|   renumber_nodes_and_elem()        25        0.0440      0.001761    0.0440      0.001761    0.40     0.40     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0361      0.036082    0.0386      0.038566    0.33     0.35     |
|   broadcast()                      1         0.0011      0.001097    0.0017      0.001697    0.01     0.02     |
|   compute_hilbert_indices()        28        0.2093      0.007474    0.2093      0.007474    1.91     1.91     |
|   find_global_indices()            28        0.1068      0.003815    0.3652      0.013044    0.97     3.33     |
|   parallel_sort()                  28        0.0194      0.000693    0.0269      0.000962    0.18     0.25     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         4         0.0003      0.000067    0.1452      0.036293    0.00     1.32     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              50        0.0322      0.000644    0.0334      0.000667    0.29     0.30     |
|   _refine_elements()               55        0.3009      0.005471    0.5977      0.010867    2.74     5.45     |
|   add_point()                      35664     0.1350      0.000004    0.2774      0.000008    1.23     2.53     |
|   make_coarsening_compatible()     113       0.6162      0.005453    0.6162      0.005453    5.62     5.62     |
|   make_refinement_compatible()     113       0.0322      0.000285    0.0331      0.000293    0.29     0.30     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      27        1.6401      0.060745    2.0081      0.074373    14.96    18.32    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      143       0.0105      0.000073    0.0117      0.000082    0.10     0.11     |
|   barrier()                        1         0.0000      0.000026    0.0000      0.000026    0.00     0.00     |
|   broadcast()                      11        0.0004      0.000038    0.0003      0.000031    0.00     0.00     |
|   gather()                         13        0.0004      0.000028    0.0004      0.000028    0.00     0.00     |
|   max(bool)                        269       0.0190      0.000071    0.0190      0.000071    0.17     0.17     |
|   max(scalar)                      5129      0.0352      0.000007    0.0352      0.000007    0.32     0.32     |
|   max(vector)                      1251      0.0167      0.000013    0.0415      0.000033    0.15     0.38     |
|   min(bool)                        6492      0.0438      0.000007    0.0438      0.000007    0.40     0.40     |
|   min(scalar)                      5069      0.3531      0.000070    0.3531      0.000070    3.22     3.22     |
|   min(vector)                      1251      0.0177      0.000014    0.0440      0.000035    0.16     0.40     |
|   probe()                          3242      0.0522      0.000016    0.0522      0.000016    0.48     0.48     |
|   receive()                        3216      0.0211      0.000007    0.0738      0.000023    0.19     0.67     |
|   send()                           3084      0.0097      0.000003    0.0097      0.000003    0.09     0.09     |
|   send_receive()                   3140      0.0241      0.000008    0.1141      0.000036    0.22     1.04     |
|   sum()                            271       0.0712      0.000263    0.0798      0.000294    0.65     0.73     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           3106      0.0056      0.000002    0.0056      0.000002    0.05     0.05     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         27        0.0967      0.003581    0.1406      0.005207    0.88     1.28     |
|   set_parent_processor_ids()       27        0.0877      0.003248    0.0877      0.003248    0.80     0.80     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          50        0.1810      0.003619    0.1810      0.003619    1.65     1.65     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       75        0.0463      0.000618    0.4333      0.005777    0.42     3.95     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       50        0.1388      0.002775    0.4441      0.008882    1.27     4.05     |
|   calculate_norm()                 51        0.0668      0.001309    0.3801      0.007454    0.61     3.47     |
|   project_vector()                 76        0.2103      0.002767    1.0174      0.013387    1.92     9.28     |
|                                                                                                                |
| XdrIO                                                                                                          |
|   write()                          1         0.0118      0.011750    0.0157      0.015711    0.11     0.14     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            227024    10.9636                                         100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example adaptivity_ex2:
*  mpirun -np 12 example-devel -n_timesteps 25 -n_refinements 5 -output_freq 10 -init_timestep 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
***** Finished first 25 steps, now read in saved solution and continue *****
 
***************************************************************
* Running Example adaptivity_ex2:
*  mpirun -np 12 example-devel -read_solution -n_timesteps 25 -output_freq 10 -init_timestep 25 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Usage:
	 /workspace/libmesh/examples/adaptivity/adaptivity_ex2/.libs/lt-example-devel -init_timestep 0
OR
	 /workspace/libmesh/examples/adaptivity/adaptivity_ex2/.libs/lt-example-devel -read_solution -init_timestep 26

Running: /workspace/libmesh/examples/adaptivity/adaptivity_ex2/.libs/lt-example-devel -read_solution -n_timesteps 25 -output_freq 10 -init_timestep 25 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=713
    n_local_nodes()=62
  n_elem()=1018
    n_local_elem()=112
    n_active_elem()=766
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

Initial H1 norm = 0.509057

 EquationSystems
  n_systems()=1
   System #0, "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=713
    n_local_dofs()=62
    n_constrained_dofs()=122
    n_local_constrained_dofs()=19
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.35344
      Average Off-Processor Bandwidth <= 1.53436
      Maximum  On-Processor Bandwidth <= 19
      Maximum Off-Processor Bandwidth <= 12
    DofMap Constraints
      Number of DoF Constraints = 122
      Average DoF Constraint Length= 2
      Number of Node Constraints = 239
      Maximum Node Constraint Length= 5
      Average Node Constraint Length= 2.53138

 Solving time step 25, time=0.6500...
H1 norm = 0.495174
  Refining the mesh...
H1 norm = 0.495112
 Solving time step 26, time=0.6750...
H1 norm = 0.481983
  Refining the mesh...
H1 norm = 0.481907
 Solving time step 27, time=0.7000...
H1 norm = 0.46947
  Refining the mesh...
H1 norm = 0.469405
 Solving time step 28, time=0.7250...
H1 norm = 0.457612
  Refining the mesh...
H1 norm = 0.45754
 Solving time step 29, time=0.7500...
H1 norm = 0.446345
  Refining the mesh...
H1 norm = 0.446298
 Solving time step 30, time=0.7750...
H1 norm = 0.435654
  Refining the mesh...
H1 norm = 0.435614
 Solving time step 31, time=0.8000...
H1 norm = 0.425499
  Refining the mesh...
H1 norm = 0.425442
 Solving time step 32, time=0.8250...
H1 norm = 0.415779
  Refining the mesh...
H1 norm = 0.415738
 Solving time step 33, time=0.8500...
H1 norm = 0.406507
  Refining the mesh...
H1 norm = 0.406469
 Solving time step 34, time=0.8750...
H1 norm = 0.397619
  Refining the mesh...
H1 norm = 0.397581
 Solving time step 35, time=0.9000...
H1 norm = 0.389087
  Refining the mesh...
H1 norm = 0.389069
 Solving time step 36, time=0.9250...
H1 norm = 0.380921
  Refining the mesh...
H1 norm = 0.380885
 Solving time step 37, time=0.9500...
H1 norm = 0.373064
  Refining the mesh...
H1 norm = 0.373044
 Solving time step 38, time=0.9750...
H1 norm = 0.365554
  Refining the mesh...
H1 norm = 0.36553
 Solving time step 39, time=1.0000...
H1 norm = 0.358358
  Refining the mesh...
H1 norm = 0.358315
 Solving time step 40, time=1.0250...
H1 norm = 0.351439
  Refining the mesh...
H1 norm = 0.351413
 Solving time step 41, time=1.0500...
H1 norm = 0.34482
  Refining the mesh...
H1 norm = 0.344802
 Solving time step 42, time=1.0750...
H1 norm = 0.338466
  Refining the mesh...
H1 norm = 0.338451
 Solving time step 43, time=1.1000...
H1 norm = 0.332347
  Refining the mesh...
H1 norm = 0.332331
 Solving time step 44, time=1.1250...
H1 norm = 0.32644
  Refining the mesh...
H1 norm = 0.32642
 Solving time step 45, time=1.1500...
H1 norm = 0.320721
  Refining the mesh...
H1 norm = 0.320709
 Solving time step 46, time=1.1750...
H1 norm = 0.3152
  Refining the mesh...
H1 norm = 0.315188
 Solving time step 47, time=1.2000...
H1 norm = 0.309858
  Refining the mesh...
H1 norm = 0.309843
 Solving time step 48, time=1.2250...
H1 norm = 0.304687
  Refining the mesh...
H1 norm = 0.30468
 Solving time step 49, time=1.2500...
H1 norm = 0.2997
  Refining the mesh...
H1 norm = 0.299685
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/adaptivity/adaptivity_ex2/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 21:58:44 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           1.631e+01      1.00000   1.631e+01
Objects:              2.555e+03      1.00000   2.555e+03
Flops:                6.334e+06      1.44462   5.227e+06  6.272e+07
Flops/sec:            3.883e+05      1.44462   3.205e+05  3.846e+06
MPI Messages:         9.004e+03      1.20931   8.086e+03  9.704e+04
MPI Message Lengths:  1.229e+06      1.11967   1.448e+02  1.405e+07
MPI Reductions:       5.859e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.6311e+01 100.0%  6.2724e+07 100.0%  9.704e+04 100.0%  1.448e+02      100.0%  5.858e+03 100.0% 

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

VecMDot              287 1.0 3.7153e-03 1.1 2.18e+05 1.3 0.0e+00 0.0e+00 2.9e+02  0  4  0  0  5   0  4  0  0  5   633
VecNorm              387 1.0 1.7657e-02 2.1 8.45e+04 1.2 0.0e+00 0.0e+00 3.9e+02  0  1  0  0  7   0  1  0  0  7    52
VecScale             337 1.0 2.6035e-04 1.3 3.68e+04 1.2 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  1522
VecCopy              302 1.0 2.7132e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               871 1.0 6.8665e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              100 1.0 8.6045e-0359.9 2.20e+04 1.2 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0    27
VecMAXPY             337 1.0 3.0303e-04 1.2 2.82e+05 1.3 0.0e+00 0.0e+00 0.0e+00  0  5  0  0  0   0  5  0  0  0 10024
VecAssemblyBegin     803 1.0 9.9451e-02 1.7 0.00e+00 0.0 4.4e+03 1.3e+02 2.2e+03  1  0  5  4 37   1  0  5  4 37     0
VecAssemblyEnd       803 1.0 8.2541e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     1015 1.0 4.4456e-03 1.1 0.00e+00 0.0 5.1e+04 1.6e+02 0.0e+00  0  0 52 57  0   0  0 52 57  0     0
VecScatterEnd       1015 1.0 1.8093e-02 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize         337 1.0 1.7727e-02 2.1 1.10e+05 1.2 0.0e+00 0.0e+00 3.4e+02  0  2  0  0  6   0  2  0  0  6    67
MatMult              337 1.0 1.8552e-02 1.5 5.88e+05 1.3 1.8e+04 7.7e+01 0.0e+00  0 10 19 10  0   0 10 19 10  0   337
MatSolve             387 1.0 2.7087e-03 1.3 2.42e+06 1.4 0.0e+00 0.0e+00 0.0e+00  0 39  0  0  0   0 39  0  0  0  9041
MatLUFactorNum        50 1.0 5.7805e-03 1.5 2.70e+06 1.7 0.0e+00 0.0e+00 0.0e+00  0 40  0  0  0   0 40  0  0  0  4333
MatILUFactorSym       50 1.0 1.3196e-02 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 1.5e+02  0  0  0  0  3   0  0  0  0  3     0
MatAssemblyBegin     100 1.0 1.5637e-01 1.9 0.00e+00 0.0 4.6e+03 4.0e+02 2.0e+02  1  0  5 13  3   1  0  5 13  3     0
MatAssemblyEnd       100 1.0 6.5057e-03 1.1 0.00e+00 0.0 2.8e+03 2.1e+01 2.1e+02  0  0  3  0  4   0  0  3  0  4     0
MatGetRowIJ           50 1.0 1.3828e-05 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        50 1.0 1.4324e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+02  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries       102 1.0 1.7452e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       287 1.0 4.2262e-03 1.1 4.37e+05 1.3 0.0e+00 0.0e+00 2.9e+02  0  8  0  0  5   0  8  0  0  5  1116
KSPSetUp             100 1.0 1.0200e-03 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              50 1.0 7.5028e-02 1.0 6.33e+06 1.4 1.8e+04 7.7e+01 9.8e+02  0100 19 10 17   0100 19 10 17   836
PCSetUp              100 1.0 3.0409e-02 1.2 2.70e+06 1.7 0.0e+00 0.0e+00 3.0e+02  0 40  0  0  5   0 40  0  0  5   824
PCSetUpOnBlocks       50 1.0 2.6854e-02 1.3 2.70e+06 1.7 0.0e+00 0.0e+00 2.5e+02  0 40  0  0  4   0 40  0  0  4   933
PCApply              387 1.0 6.9437e-03 1.1 2.42e+06 1.4 0.0e+00 0.0e+00 0.0e+00  0 39  0  0  0   0 39  0  0  0  3527
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector  1180           1180      3696120     0
      Vector Scatter   355            355       367780     0
           Index Set   659            659       547400     0
   IS L to G Mapping   128            128        72192     0
              Matrix   128            128      2505612     0
       Krylov Solver    52             52       503360     0
      Preconditioner    52             52        46384     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 5.19753e-06
Average time for zero size MPI_Send(): 1.33316e-05
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
| Time:           Thu Jan 31 21:58:44 2013                                                                             |
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
| libMesh Performance: Alive time=16.3322, Active time=15.7699                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     26        0.1913      0.007356    0.3016      0.011599    1.21     1.91     |
|   build_constraint_matrix()        7286      0.0687      0.000009    0.0687      0.000009    0.44     0.44     |
|   build_sparsity()                 26        0.1239      0.004767    0.3457      0.013295    0.79     2.19     |
|   cnstrn_elem_mat_vec()            7286      0.0443      0.000006    0.0443      0.000006    0.28     0.28     |
|   create_dof_constraints()         26        0.6010      0.023114    1.2310      0.047345    3.81     7.81     |
|   distribute_dofs()                26        0.6539      0.025149    2.1902      0.084240    4.15     13.89    |
|   dof_indices()                    61639     2.4718      0.000040    2.4718      0.000040    15.67    15.67    |
|   enforce_constraints_exactly()    75        0.0153      0.000203    0.0153      0.000203    0.10     0.10     |
|   old_dof_indices()                22275     0.9588      0.000043    0.9588      0.000043    6.08     6.08     |
|   prepare_send_list()              26        0.0015      0.000056    0.0015      0.000056    0.01     0.01     |
|   reinit()                         26        1.4192      0.054584    1.4192      0.054584    9.00     9.00     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          4         0.0084      0.002099    0.0397      0.009916    0.05     0.25     |
|   read()                           1         0.0177      0.017733    0.1473      0.147273    0.11     0.93     |
|   update()                         1         0.0003      0.000266    0.0003      0.000266    0.00     0.00     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        26276     0.2165      0.000008    0.2165      0.000008    1.37     1.37     |
|   init_shape_functions()           12101     0.1144      0.000009    0.1144      0.000009    0.73     0.73     |
|   inverse_map()                    45463     0.2605      0.000006    0.2605      0.000006    1.65     1.65     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             26276     0.2313      0.000009    0.2313      0.000009    1.47     1.47     |
|   compute_face_map()               5848      0.0938      0.000016    0.1767      0.000030    0.60     1.12     |
|   init_face_shape_functions()      141       0.0011      0.000008    0.0011      0.000008    0.01     0.01     |
|   init_reference_to_physical_map() 12101     0.1054      0.000009    0.1054      0.000009    0.67     0.67     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               4         0.0406      0.010148    0.0416      0.010392    0.26     0.26     |
|                                                                                                                |
| JumpErrorEstimator                                                                                             |
|   estimate_error()                 25        0.4410      0.017639    1.9260      0.077038    2.80     12.21    |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           7080      0.0275      0.000004    0.0275      0.000004    0.17     0.17     |
|   init()                           50        0.0975      0.001950    0.0975      0.001950    0.62     0.62     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       25        0.0283      0.001131    0.0535      0.002139    0.18     0.34     |
|   find_neighbors()                 26        1.2147      0.046719    1.2472      0.047968    7.70     7.91     |
|   renumber_nodes_and_elem()        77        0.0669      0.000868    0.0669      0.000868    0.42     0.42     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0355      0.035487    0.0382      0.038198    0.23     0.24     |
|   compute_hilbert_indices()        26        0.3553      0.013666    0.3553      0.013666    2.25     2.25     |
|   find_global_indices()            26        0.1598      0.006144    0.5884      0.022629    1.01     3.73     |
|   parallel_sort()                  26        0.0245      0.000941    0.0367      0.001412    0.16     0.23     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         4         0.0003      0.000065    0.0824      0.020610    0.00     0.52     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              50        0.0337      0.000674    0.0374      0.000748    0.21     0.24     |
|   _refine_elements()               50        0.1386      0.002772    0.2061      0.004123    0.88     1.31     |
|   add_point()                      7080      0.0289      0.000004    0.0578      0.000008    0.18     0.37     |
|   make_coarsening_compatible()     118       1.1387      0.009650    1.1387      0.009650    7.22     7.22     |
|   make_refinement_compatible()     118       0.0693      0.000587    0.0732      0.000621    0.44     0.46     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      26        2.0200      0.077691    2.6153      0.100589    12.81    16.58    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      137       0.0490      0.000358    0.0496      0.000362    0.31     0.31     |
|   broadcast()                      44        0.0008      0.000018    0.0007      0.000016    0.01     0.00     |
|   max(bool)                        269       0.0320      0.000119    0.0320      0.000119    0.20     0.20     |
|   max(scalar)                      5101      0.0371      0.000007    0.0371      0.000007    0.24     0.24     |
|   max(vector)                      1245      0.0163      0.000013    0.0403      0.000032    0.10     0.26     |
|   min(bool)                        6474      0.0811      0.000013    0.0811      0.000013    0.51     0.51     |
|   min(scalar)                      5045      0.5935      0.000118    0.5935      0.000118    3.76     3.76     |
|   min(vector)                      1245      0.0173      0.000014    0.0426      0.000034    0.11     0.27     |
|   probe()                          2974      0.0690      0.000023    0.0690      0.000023    0.44     0.44     |
|   receive()                        2952      0.0205      0.000007    0.0900      0.000031    0.13     0.57     |
|   send()                           2952      0.0090      0.000003    0.0090      0.000003    0.06     0.06     |
|   send_receive()                   2982      0.0235      0.000008    0.1297      0.000044    0.15     0.82     |
|   sum()                            261       0.0453      0.000174    0.0555      0.000213    0.29     0.35     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           2952      0.0057      0.000002    0.0057      0.000002    0.04     0.04     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         27        0.1105      0.004093    0.1723      0.006381    0.70     1.09     |
|   set_parent_processor_ids()       26        0.1327      0.005106    0.1327      0.005106    0.84     0.84     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          50        0.2400      0.004800    0.2400      0.004800    1.52     1.52     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       75        0.1137      0.001516    1.0801      0.014402    0.72     6.85     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       50        0.2815      0.005630    0.9203      0.018405    1.78     5.84     |
|   calculate_norm()                 51        0.1275      0.002499    0.6404      0.012556    0.81     4.06     |
|   project_vector()                 75        0.2336      0.003114    1.8291      0.024387    1.48     11.60    |
|                                                                                                                |
| XdrIO                                                                                                          |
|   read()                           1         0.0106      0.010635    0.0109      0.010929    0.07     0.07     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            276729    15.7699                                         100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example adaptivity_ex2:
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
