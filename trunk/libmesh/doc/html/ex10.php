<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex10",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 10 - Solving a Transient System with Adaptive Mesh Refinement</h1>

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
        #include "equation_systems.h"
        #include "fe.h"
        #include "quadrature_gauss.h"
        #include "dof_map.h"
        #include "sparse_matrix.h"
        #include "numeric_vector.h"
        #include "dense_matrix.h"
        #include "dense_vector.h"
        
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
        
        #ifdef LIBMESH_ENABLE_PARMESH
</pre>
</div>
<div class = "comment">
We still need some work on automatic parallel restarts
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(false, "--disable-parmesh");
        #endif
        
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
</div>

<div class ="fragment">
<pre>
          const Real dt = 0.025;
          Real time     = init_timestep*dt;
          
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
                  OStringStream file_name;
        
                  file_name &lt;&lt; "out.gmv.";
                  OSSRealzeroright(file_name,3,0,t_step+1);
        
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
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example  ./ex10-dbg [-read_solution] -n_timesteps 25 -n_refinements 5 -init_timestep [0|25]
***************************************************************
 
Usage:
	 ./ex10-dbg -init_timestep 0
OR
	 ./ex10-dbg -read_solution -init_timestep 26

Running: ./ex10-dbg -n_timesteps 25 -n_refinements 5 -output_freq 10 -init_timestep 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=6273
    n_local_nodes()=6273
  n_elem()=13650
    n_local_elem()=13650
    n_active_elem()=10240
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
    n_dofs()=6273
    n_local_dofs()=6273
    n_constrained_dofs()=0
    n_vectors()=3
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.57038
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0

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

WARNING! There are options you set that were not used!
WARNING! could be spelling mistake, etc!
Option left: name:-init_timestep value: 0
Option left: name:-n_refinements value: 5
Option left: name:-n_timesteps value: 25
Option left: name:-output_freq value: 10

 ---------------------------------------------------------------------------- 
| Reference count information                                                |
 ---------------------------------------------------------------------------- 
| N7libMesh10Parameters5ValueE reference count information:
|  Creations:    6
|  Destructions: 6
| N7libMesh12LinearSolverIdEE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh12SparseMatrixIdEE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh13NumericVectorIdEE reference count information:
|  Creations:    258
|  Destructions: 258
| N7libMesh15EquationSystemsE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh4ElemE reference count information:
|  Creations:    232245
|  Destructions: 232245
| N7libMesh4NodeE reference count information:
|  Creations:    6989
|  Destructions: 6989
| N7libMesh5QBaseE reference count information:
|  Creations:    940
|  Destructions: 940
| N7libMesh6DofMapE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh6FEBaseE reference count information:
|  Creations:    378
|  Destructions: 378
| N7libMesh6SystemE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh9DofObjectE reference count information:
|  Creations:    291916
|  Destructions: 291916
 ---------------------------------------------------------------------------- 

-------------------------------------------------------------------
| Time:           Thu Dec  1 21:21:11 2011                         |
| OS:             Linux                                            |
| HostName:       dknez-laptop                                     |
| OS Release:     3.0.0-13-generic                                 |
| OS Version:     #22-Ubuntu SMP Wed Nov 2 13:27:26 UTC 2011       |
| Machine:        x86_64                                           |
| Username:       dknez                                            |
| Configuration:  ./configure run on Tue Nov 22 16:02:41 EST 2011  |
-------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=25.7582, Active time=24.9602                                                |
 -------------------------------------------------------------------------------------------------------------
| Event                           nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                           w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------|
|                                                                                                             |
|                                                                                                             |
| DofMap                                                                                                      |
|   add_neighbors_to_send_list()  26        0.2094      0.008054    0.2094      0.008054    0.84     0.84     |
|   build_constraint_matrix()     29356     0.2893      0.000010    0.2893      0.000010    1.16     1.16     |
|   cnstrn_elem_mat_vec()         29356     0.1919      0.000007    0.1919      0.000007    0.77     0.77     |
|   compute_sparsity()            26        2.7213      0.104666    2.9897      0.114990    10.90    11.98    |
|   create_dof_constraints()      26        0.2691      0.010350    0.3992      0.015353    1.08     1.60     |
|   distribute_dofs()             26        0.7061      0.027158    1.6521      0.063542    2.83     6.62     |
|   dof_indices()                 267936    2.7678      0.000010    2.7678      0.000010    11.09    11.09    |
|   enforce_constraints_exactly() 75        0.0379      0.000505    0.0379      0.000505    0.15     0.15     |
|   old_dof_indices()             121164    1.2267      0.000010    1.2267      0.000010    4.91     4.91     |
|   prepare_send_list()           26        0.0002      0.000006    0.0002      0.000006    0.00     0.00     |
|   reinit()                      26        0.9459      0.036380    0.9459      0.036380    3.79     3.79     |
|                                                                                                             |
| EquationSystems                                                                                             |
|   build_solution_vector()       4         0.0723      0.018083    0.1867      0.046678    0.29     0.75     |
|   write()                       1         0.1980      0.198018    0.1981      0.198069    0.79     0.79     |
|                                                                                                             |
| FE                                                                                                          |
|   compute_affine_map()          168850    1.1632      0.000007    1.1632      0.000007    4.66     4.66     |
|   compute_face_map()            45464     0.6349      0.000014    1.2362      0.000027    2.54     4.95     |
|   compute_shape_functions()     168850    1.1522      0.000007    1.1522      0.000007    4.62     4.62     |
|   init_face_shape_functions()   434       0.0029      0.000007    0.0029      0.000007    0.01     0.01     |
|   init_shape_functions()        89895     1.4674      0.000016    1.4674      0.000016    5.88     5.88     |
|   inverse_map()                 196675    1.1708      0.000006    1.1708      0.000006    4.69     4.69     |
|                                                                                                             |
| GMVIO                                                                                                       |
|   write_nodal_data()            4         0.0949      0.023727    0.0949      0.023727    0.38     0.38     |
|                                                                                                             |
| JumpErrorEstimator                                                                                          |
|   estimate_error()              25        2.6843      0.107374    7.9070      0.316279    10.75    31.68    |
|                                                                                                             |
| LocationMap                                                                                                 |
|   find()                        35664     0.1043      0.000003    0.1043      0.000003    0.42     0.42     |
|   init()                        55        0.0494      0.000899    0.0494      0.000899    0.20     0.20     |
|                                                                                                             |
| Mesh                                                                                                        |
|   contract()                    25        0.0895      0.003581    0.1224      0.004894    0.36     0.49     |
|   find_neighbors()              27        1.5200      0.056296    1.5200      0.056296    6.09     6.09     |
|   renumber_nodes_and_elem()     77        0.0792      0.001029    0.0792      0.001029    0.32     0.32     |
|                                                                                                             |
| MeshCommunication                                                                                           |
|   assign_global_indices()       1         0.2284      0.228435    0.2285      0.228506    0.92     0.92     |
|                                                                                                             |
| MeshOutput                                                                                                  |
|   write_equation_systems()      4         0.0001      0.000028    0.2817      0.070435    0.00     1.13     |
|                                                                                                             |
| MeshRefinement                                                                                              |
|   _coarsen_elements()           50        0.0315      0.000631    0.0315      0.000631    0.13     0.13     |
|   _refine_elements()            55        0.2234      0.004062    0.5180      0.009418    0.90     2.08     |
|   add_point()                   35664     0.1683      0.000005    0.2841      0.000008    0.67     1.14     |
|   make_coarsening_compatible()  113       0.6096      0.005395    0.6096      0.005395    2.44     2.44     |
|   make_refinement_compatible()  113       0.0373      0.000330    0.0373      0.000330    0.15     0.15     |
|                                                                                                             |
| Parallel                                                                                                    |
|   allgather()                   30        0.0001      0.000003    0.0001      0.000003    0.00     0.00     |
|   receive()                     4         0.0000      0.000005    0.0000      0.000005    0.00     0.00     |
|   send()                        4         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   send_receive()                4         0.0001      0.000013    0.0001      0.000013    0.00     0.00     |
|                                                                                                             |
| Partitioner                                                                                                 |
|   single_partition()            27        0.0252      0.000933    0.0252      0.000933    0.10     0.10     |
|                                                                                                             |
| PetscLinearSolver                                                                                           |
|   solve()                       50        0.1164      0.002329    0.1164      0.002329    0.47     0.47     |
|                                                                                                             |
| ProjectVector                                                                                               |
|   operator()                    75        0.5758      0.007677    1.5628      0.020837    2.31     6.26     |
|                                                                                                             |
| System                                                                                                      |
|   assemble()                    50        1.8549      0.037098    3.6022      0.072044    7.43     14.43    |
|   calculate_norm()              51        0.5975      0.011716    1.6949      0.033234    2.39     6.79     |
|   project_vector()              76        0.6218      0.008182    3.1010      0.040802    2.49     12.42    |
|                                                                                                             |
| XdrIO                                                                                                       |
|   write()                       1         0.0205      0.020480    0.0205      0.020480    0.08     0.08     |
 -------------------------------------------------------------------------------------------------------------
| Totals:                         1190465   24.9602                                         100.00            |
 -------------------------------------------------------------------------------------------------------------

 
***** Finished first 25 steps, now read in saved solution and continue *****
 
Usage:
	 ./ex10-dbg -init_timestep 0
OR
	 ./ex10-dbg -read_solution -init_timestep 26

Running: ./ex10-dbg -read_solution -n_timesteps 25 -output_freq 10 -init_timestep 25

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=713
    n_local_nodes()=713
  n_elem()=1018
    n_local_elem()=1018
    n_active_elem()=766
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
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
    n_local_dofs()=713
    n_constrained_dofs()=122
    n_vectors()=3
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 9.10659
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 17
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 122    DofMap Constraints
      Number of Constraints = 122
      Average DoF Constraint Length= 2
      Number of Node Constraints = 239
      Maximum Node Constraint Length= 4
      Average Node Constraint Length= 2.04184

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
 Solving time step 40, time=1.0300...
H1 norm = 0.351439
  Refining the mesh...
H1 norm = 0.351413
 Solving time step 41, time=1.0500...
H1 norm = 0.34482
  Refining the mesh...
make[1]: *** [run] Interrupt
