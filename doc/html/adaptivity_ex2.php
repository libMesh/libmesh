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
Create a new mesh on the default MPI communicator.
We still need some work on automatic parallel restarts with
ParallelMesh
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
  
    SerialMesh mesh(init.comm());
  
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
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex2'
***************************************************************
* Running Example adaptivity_ex2:
*  mpirun -np 4 example-devel -n_timesteps 25 -n_refinements 5 -output_freq 10 -init_timestep 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
Usage:
	 /net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex2/.libs/lt-example-devel -init_timestep 0
OR
	 /net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex2/.libs/lt-example-devel -read_solution -init_timestep 26

Running: /net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex2/.libs/lt-example-devel -n_timesteps 25 -n_refinements 5 -output_freq 10 -init_timestep 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=6273
    n_local_nodes()=1789
  n_elem()=13650
    n_local_elem()=3277
    n_active_elem()=10240
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
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
    n_local_dofs()=1789
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.53164
      Average Off-Processor Bandwidth <= 0.127212
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 5
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
H1 norm = 0.843309
  Refining the mesh...
H1 norm = 0.843131
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
H1 norm = 0.592783
  Refining the mesh...
H1 norm = 0.592698
 Solving time step 20, time=0.5250...
H1 norm = 0.573854
  Refining the mesh...
H1 norm = 0.573819
 Solving time step 21, time=0.5500...
H1 norm = 0.556154
  Refining the mesh...
H1 norm = 0.556109
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


 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:45:24 2013                                                                          |
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
| libMesh Performance: Alive time=4.02198, Active time=3.79428                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     26        0.0259      0.000996    0.0327      0.001258    0.68     0.86     |
|   build_constraint_matrix()        7406      0.0095      0.000001    0.0095      0.000001    0.25     0.25     |
|   build_sparsity()                 26        0.0199      0.000764    0.0645      0.002482    0.52     1.70     |
|   cnstrn_elem_mat_vec()            7406      0.0118      0.000002    0.0118      0.000002    0.31     0.31     |
|   create_dof_constraints()         26        0.0519      0.001998    0.0800      0.003077    1.37     2.11     |
|   distribute_dofs()                26        0.0606      0.002331    0.2250      0.008655    1.60     5.93     |
|   dof_indices()                    74816     0.2320      0.000003    0.2320      0.000003    6.12     6.12     |
|   enforce_constraints_exactly()    75        0.0066      0.000088    0.0066      0.000088    0.17     0.17     |
|   old_dof_indices()                31512     0.0872      0.000003    0.0872      0.000003    2.30     2.30     |
|   prepare_send_list()              26        0.0002      0.000007    0.0002      0.000007    0.00     0.00     |
|   reinit()                         26        0.0947      0.003644    0.0947      0.003644    2.50     2.50     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          4         0.0080      0.002009    0.0236      0.005907    0.21     0.62     |
|   write()                          1         0.0015      0.001458    0.0035      0.003533    0.04     0.09     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        42008     0.0586      0.000001    0.0586      0.000001    1.54     1.54     |
|   init_shape_functions()           22651     0.0410      0.000002    0.0410      0.000002    1.08     1.08     |
|   inverse_map()                    55483     0.1020      0.000002    0.1020      0.000002    2.69     2.69     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             42008     0.0647      0.000002    0.0647      0.000002    1.71     1.71     |
|   compute_face_map()               11364     0.0413      0.000004    0.0907      0.000008    1.09     2.39     |
|   init_face_shape_functions()      203       0.0004      0.000002    0.0004      0.000002    0.01     0.01     |
|   init_reference_to_physical_map() 22651     0.0495      0.000002    0.0495      0.000002    1.31     1.31     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               4         0.0633      0.015824    0.0636      0.015894    1.67     1.68     |
|                                                                                                                |
| JumpErrorEstimator                                                                                             |
|   estimate_error()                 25        0.1660      0.006638    0.8787      0.035150    4.37     23.16    |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           35664     0.0417      0.000001    0.0417      0.000001    1.10     1.10     |
|   init()                           55        0.0114      0.000208    0.0114      0.000208    0.30     0.30     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       25        0.0070      0.000279    0.0141      0.000565    0.18     0.37     |
|   find_neighbors()                 27        0.1236      0.004579    0.1649      0.006109    3.26     4.35     |
|   renumber_nodes_and_elem()        25        0.0072      0.000286    0.0072      0.000286    0.19     0.19     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0094      0.009403    0.0146      0.014614    0.25     0.39     |
|   broadcast()                      1         0.0002      0.000218    0.0005      0.000501    0.01     0.01     |
|   compute_hilbert_indices()        28        0.0816      0.002914    0.0816      0.002914    2.15     2.15     |
|   find_global_indices()            28        0.0156      0.000557    0.1325      0.004732    0.41     3.49     |
|   parallel_sort()                  28        0.0028      0.000100    0.0303      0.001082    0.07     0.80     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         4         0.0001      0.000031    0.0876      0.021907    0.00     2.31     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              50        0.0089      0.000177    0.0151      0.000301    0.23     0.40     |
|   _refine_elements()               55        0.0718      0.001305    0.2000      0.003636    1.89     5.27     |
|   add_point()                      35664     0.0607      0.000002    0.1062      0.000003    1.60     2.80     |
|   make_coarsening_compatible()     113       0.1443      0.001277    0.1443      0.001277    3.80     3.80     |
|   make_flags_parallel_consistent() 75        0.0387      0.000516    0.0885      0.001180    1.02     2.33     |
|   make_refinement_compatible()     113       0.0068      0.000060    0.0122      0.000108    0.18     0.32     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      27        0.2863      0.010602    0.4282      0.015858    7.54     11.28    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      143       0.0498      0.000349    0.0503      0.000352    1.31     1.33     |
|   barrier()                        1         0.0000      0.000014    0.0000      0.000014    0.00     0.00     |
|   broadcast()                      29        0.0001      0.000005    0.0001      0.000003    0.00     0.00     |
|   gather()                         13        0.0006      0.000046    0.0006      0.000046    0.02     0.02     |
|   max(bool)                        269       0.0399      0.000148    0.0399      0.000148    1.05     1.05     |
|   max(scalar)                      5929      0.0318      0.000005    0.0318      0.000005    0.84     0.84     |
|   max(vector)                      1451      0.0092      0.000006    0.0278      0.000019    0.24     0.73     |
|   min(bool)                        7567      0.1282      0.000017    0.1282      0.000017    3.38     3.38     |
|   min(scalar)                      5869      0.6131      0.000104    0.6131      0.000104    16.16    16.16    |
|   min(vector)                      1451      0.0102      0.000007    0.0311      0.000021    0.27     0.82     |
|   probe()                          1790      0.0497      0.000028    0.0497      0.000028    1.31     1.31     |
|   receive()                        1780      0.0044      0.000002    0.0540      0.000030    0.11     1.42     |
|   send()                           1744      0.0024      0.000001    0.0024      0.000001    0.06     0.06     |
|   send_receive()                   1800      0.0057      0.000003    0.0637      0.000035    0.15     1.68     |
|   sum()                            271       0.0979      0.000361    0.4648      0.001715    2.58     12.25    |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           1750      0.0013      0.000001    0.0013      0.000001    0.03     0.03     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         27        0.0264      0.000978    0.0984      0.003645    0.70     2.59     |
|   set_parent_processor_ids()       27        0.0142      0.000526    0.0142      0.000526    0.37     0.37     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          50        0.2388      0.004775    0.2388      0.004775    6.29     6.29     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       75        0.0242      0.000322    0.0938      0.001251    0.64     2.47     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       50        0.0948      0.001895    0.1954      0.003908    2.50     5.15     |
|   calculate_norm()                 51        0.0412      0.000807    0.1787      0.003504    1.09     4.71     |
|   project_vector()                 76        0.1913      0.002517    0.3705      0.004876    5.04     9.77     |
|                                                                                                                |
| XdrIO                                                                                                          |
|   write()                          1         0.0044      0.004397    0.0064      0.006354    0.12     0.17     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            421970    3.7943                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example adaptivity_ex2:
*  mpirun -np 4 example-devel -n_timesteps 25 -n_refinements 5 -output_freq 10 -init_timestep 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
***** Finished first 25 steps, now read in saved solution and continue *****
 
***************************************************************
* Running Example adaptivity_ex2:
*  mpirun -np 4 example-devel -read_solution -n_timesteps 25 -output_freq 10 -init_timestep 25 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
Usage:
	 /net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex2/.libs/lt-example-devel -init_timestep 0
OR
	 /net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex2/.libs/lt-example-devel -read_solution -init_timestep 26

Running: /net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex2/.libs/lt-example-devel -read_solution -n_timesteps 25 -output_freq 10 -init_timestep 25 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=713
    n_local_nodes()=212
  n_elem()=1018
    n_local_elem()=242
    n_active_elem()=766
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
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
    n_local_dofs()=212
    n_constrained_dofs()=122
    n_local_constrained_dofs()=25
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.75035
      Average Off-Processor Bandwidth <= 0.692847
      Maximum  On-Processor Bandwidth <= 19
      Maximum Off-Processor Bandwidth <= 8
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
H1 norm = 0.425441
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
H1 norm = 0.373043
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

 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:45:30 2013                                                                          |
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
| libMesh Performance: Alive time=6.04325, Active time=5.83294                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     26        0.0338      0.001302    0.0441      0.001698    0.58     0.76     |
|   build_constraint_matrix()        20657     0.0267      0.000001    0.0267      0.000001    0.46     0.46     |
|   build_sparsity()                 26        0.0302      0.001161    0.1057      0.004065    0.52     1.81     |
|   cnstrn_elem_mat_vec()            20657     0.0386      0.000002    0.0386      0.000002    0.66     0.66     |
|   create_dof_constraints()         26        0.1008      0.003878    0.1688      0.006491    1.73     2.89     |
|   distribute_dofs()                26        0.0664      0.002552    0.2837      0.010910    1.14     4.86     |
|   dof_indices()                    141269    0.4231      0.000003    0.4231      0.000003    7.25     7.25     |
|   enforce_constraints_exactly()    75        0.0089      0.000119    0.0089      0.000119    0.15     0.15     |
|   old_dof_indices()                64149     0.1861      0.000003    0.1861      0.000003    3.19     3.19     |
|   prepare_send_list()              26        0.0002      0.000009    0.0002      0.000009    0.00     0.00     |
|   reinit()                         26        0.1196      0.004600    0.1196      0.004600    2.05     2.05     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          4         0.0035      0.000863    0.0132      0.003289    0.06     0.23     |
|   read()                           1         0.0049      0.004923    0.0386      0.038614    0.08     0.66     |
|   update()                         1         0.0001      0.000097    0.0001      0.000097    0.00     0.00     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        77027     0.1115      0.000001    0.1115      0.000001    1.91     1.91     |
|   init_shape_functions()           37014     0.0705      0.000002    0.0705      0.000002    1.21     1.21     |
|   inverse_map()                    93589     0.1739      0.000002    0.1739      0.000002    2.98     2.98     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             77027     0.1273      0.000002    0.1273      0.000002    2.18     2.18     |
|   compute_face_map()               17929     0.0687      0.000004    0.1485      0.000008    1.18     2.55     |
|   init_face_shape_functions()      419       0.0007      0.000002    0.0007      0.000002    0.01     0.01     |
|   init_reference_to_physical_map() 37014     0.0809      0.000002    0.0809      0.000002    1.39     1.39     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               4         0.0371      0.009274    0.0374      0.009341    0.64     0.64     |
|                                                                                                                |
| JumpErrorEstimator                                                                                             |
|   estimate_error()                 25        0.2873      0.011494    1.3839      0.055355    4.93     23.73    |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           7080      0.0065      0.000001    0.0065      0.000001    0.11     0.11     |
|   init()                           50        0.0180      0.000359    0.0180      0.000359    0.31     0.31     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       25        0.0059      0.000237    0.0115      0.000460    0.10     0.20     |
|   find_neighbors()                 26        0.1376      0.005291    0.1976      0.007601    2.36     3.39     |
|   renumber_nodes_and_elem()        77        0.0155      0.000202    0.0155      0.000202    0.27     0.27     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0169      0.016942    0.0174      0.017363    0.29     0.30     |
|   compute_hilbert_indices()        26        0.1158      0.004454    0.1158      0.004454    1.99     1.99     |
|   find_global_indices()            26        0.0201      0.000774    0.2101      0.008082    0.34     3.60     |
|   parallel_sort()                  26        0.0033      0.000125    0.0683      0.002626    0.06     1.17     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         4         0.0001      0.000033    0.0510      0.012739    0.00     0.87     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              50        0.0095      0.000191    0.0157      0.000313    0.16     0.27     |
|   _refine_elements()               50        0.0297      0.000593    0.0775      0.001550    0.51     1.33     |
|   add_point()                      7080      0.0108      0.000002    0.0180      0.000003    0.18     0.31     |
|   make_coarsening_compatible()     118       0.2960      0.002508    0.2960      0.002508    5.07     5.07     |
|   make_flags_parallel_consistent() 75        0.0717      0.000957    0.1367      0.001823    1.23     2.34     |
|   make_refinement_compatible()     118       0.0155      0.000132    0.0255      0.000216    0.27     0.44     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      26        0.3213      0.012356    0.5441      0.020926    5.51     9.33     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      137       0.0662      0.000483    0.0666      0.000486    1.14     1.14     |
|   broadcast()                      44        0.0002      0.000004    0.0002      0.000004    0.00     0.00     |
|   max(bool)                        269       0.0602      0.000224    0.0602      0.000224    1.03     1.03     |
|   max(scalar)                      5901      0.0327      0.000006    0.0327      0.000006    0.56     0.56     |
|   max(vector)                      1445      0.0094      0.000006    0.0271      0.000019    0.16     0.46     |
|   min(bool)                        7549      0.1830      0.000024    0.1830      0.000024    3.14     3.14     |
|   min(scalar)                      5845      0.9475      0.000162    0.9475      0.000162    16.24    16.24    |
|   min(vector)                      1445      0.0105      0.000007    0.0308      0.000021    0.18     0.53     |
|   probe()                          1714      0.0709      0.000041    0.0709      0.000041    1.22     1.22     |
|   receive()                        1708      0.0052      0.000003    0.0763      0.000045    0.09     1.31     |
|   send()                           1708      0.0028      0.000002    0.0028      0.000002    0.05     0.05     |
|   send_receive()                   1754      0.0061      0.000003    0.0872      0.000050    0.10     1.49     |
|   sum()                            261       0.2031      0.000778    0.7246      0.002776    3.48     12.42    |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           1708      0.0015      0.000001    0.0015      0.000001    0.03     0.03     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         27        0.0311      0.001152    0.1801      0.006671    0.53     3.09     |
|   set_parent_processor_ids()       26        0.0183      0.000706    0.0183      0.000706    0.31     0.31     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          50        0.4405      0.008809    0.4405      0.008809    7.55     7.55     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       75        0.0672      0.000896    0.2588      0.003450    1.15     4.44     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       50        0.2072      0.004145    0.4390      0.008781    3.55     7.53     |
|   calculate_norm()                 51        0.0844      0.001655    0.3657      0.007170    1.45     6.27     |
|   project_vector()                 75        0.2858      0.003811    0.6619      0.008825    4.90     11.35    |
|                                                                                                                |
| XdrIO                                                                                                          |
|   read()                           1         0.0038      0.003773    0.0039      0.003862    0.06     0.07     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            633718    5.8329                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example adaptivity_ex2:
*  mpirun -np 4 example-devel -read_solution -n_timesteps 25 -output_freq 10 -init_timestep 25 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex2'
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
