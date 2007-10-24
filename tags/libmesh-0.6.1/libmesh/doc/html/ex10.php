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
solution is advanced in time with a standard Crank-Nicholson
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
          libMesh::init (argc, argv);
        
        #ifndef ENABLE_AMR
          std::cerr &lt;&lt; "ERROR: This example requires libMesh to be\n"
                    &lt;&lt; "compiled with AMR support!"
                    &lt;&lt; std::endl;
          return 0;
        #else
        
          {    
        
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
              error();
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
              std::cerr &lt;&lt; "ERROR: Number of timesteps not specified\n" &lt;&lt; std::endl;
              error();
            }
        
        
</pre>
</div>
<div class = "comment">
Create a two-dimensional mesh.
</div>

<div class ="fragment">
<pre>
            Mesh mesh (2);
        
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
Uniformly refine the mesh 5 times
</div>

<div class ="fragment">
<pre>
              if(!read_solution)
                mesh_refinement.uniformly_refine (5);
        
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
Output evey 10 timesteps to file.
</div>

<div class ="fragment">
<pre>
                if ( (t_step+1)%10 == 0)
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
                mesh.write("saved_mesh.xda");
                equation_systems.write("saved_solution.xda", libMeshEnums::WRITE);
              }
          }
        #endif // #ifndef ENABLE_AMR
          
</pre>
</div>
<div class = "comment">
All done.  
</div>

<div class ="fragment">
<pre>
          return libMesh::close ();
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
          assert (system_name == "Convection-Diffusion");
        
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
        #ifdef ENABLE_AMR
</pre>
</div>
<div class = "comment">
It is a good idea to make sure we are assembling
the proper system.
</div>

<div class ="fragment">
<pre>
          assert (system_name == "Convection-Diffusion");
          
</pre>
</div>
<div class = "comment">
Get a constant reference to the mesh object.
</div>

<div class ="fragment">
<pre>
          const Mesh& mesh = es.get_mesh();
          
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
                                                0.01*(grad_u_old*dphi[i][qp]))     
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
                                                     0.01*(dphi[i][qp]*dphi[j][qp]))      
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
        #endif // #ifdef ENABLE_AMR
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
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;o_string_stream.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;transient_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;vector_value.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;kelly_error_estimator.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  
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
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::init (argc, argv);
  
  #ifndef ENABLE_AMR
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR: This example requires libMesh to be\n&quot;</FONT></B>
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;compiled with AMR support!&quot;</FONT></B>
              &lt;&lt; std::endl;
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    {    
  
  
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
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR: Initial timestep not specified\n&quot;</FONT></B> &lt;&lt; std::endl;
  
        error();
      }
  
  
  
      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_timesteps = 0;
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(command_line.search(<B><FONT COLOR="#BC8F8F">&quot;-n_timesteps&quot;</FONT></B>))
        n_timesteps = command_line.next(0);
      <B><FONT COLOR="#A020F0">else</FONT></B>
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR: Number of timesteps not specified\n&quot;</FONT></B> &lt;&lt; std::endl;
        error();
      }
  
  
      Mesh mesh (2);
  
      EquationSystems equation_systems (mesh);
      MeshRefinement mesh_refinement (mesh);
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(!read_solution)
      {
        mesh.read (<B><FONT COLOR="#BC8F8F">&quot;mesh.xda&quot;</FONT></B>);
  
        <B><FONT COLOR="#A020F0">if</FONT></B>(!read_solution)
          mesh_refinement.uniformly_refine (5);
  
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
      
      <B><FONT COLOR="#228B22">const</FONT></B> Real dt = 0.025;
      Real time     = init_timestep*dt;
      
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
  	
  	<B><FONT COLOR="#A020F0">if</FONT></B> ( (t_step+1)%10 == 0)
  	  {
  	    OStringStream file_name;
  
  	    file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out.gmv.&quot;</FONT></B>;
  	    OSSRealzeroright(file_name,3,0,t_step+1);
  
  	    GMVIO(mesh).write_equation_systems (file_name.str(),
  						equation_systems);
  	  }
        }
  
        <B><FONT COLOR="#A020F0">if</FONT></B>(!read_solution)
        {
          mesh.write(<B><FONT COLOR="#BC8F8F">&quot;saved_mesh.xda&quot;</FONT></B>);
          equation_systems.write(<B><FONT COLOR="#BC8F8F">&quot;saved_solution.xda&quot;</FONT></B>, libMeshEnums::WRITE);
        }
    }
  #endif <I><FONT COLOR="#B22222">// #ifndef ENABLE_AMR
</FONT></I>    
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close ();
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> init_cd (EquationSystems&amp; es,
  	      <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    TransientLinearImplicitSystem &amp; system =
      es.get_system&lt;TransientLinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    es.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>) = 0;
    
    system.project_solution(exact_value, NULL, es.parameters);
  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_cd (EquationSystems&amp; es,
  		  <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
  #ifdef ENABLE_AMR
    assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> Mesh&amp; mesh = es.get_mesh();
    
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
  					
  					0.01*(grad_u_old*dphi[i][qp]))     
  				);
  	      
  	      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
  		{
  		  Ke(i,j) += JxW[qp]*(
  				      phi[i][qp]*phi[j][qp] + 
  				      .5*dt*(
  					     (velocity*dphi[j][qp])*phi[i][qp] +
  					     0.01*(dphi[i][qp]*dphi[j][qp]))      
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
  #endif <I><FONT COLOR="#B22222">// #ifdef ENABLE_AMR
</FONT></I>  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example  ./ex10-devel
***************************************************************
 
Usage:
	 ./ex10-devel -init_timestep 0
OR
	 ./ex10-devel -read_solution -init_timestep 26

Running: ./ex10-devel -n_timesteps 25 -init_timestep 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=6273
  n_elem()=13650
   n_local_elem()=13650
   n_active_elem()=10240
  n_subdomains()=1
  n_processors()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE" 
    Approximation Orders="FIRST" 
    n_dofs()=6273
    n_local_dofs()=6273
    n_constrained_dofs()=0
    n_vectors()=3

 Solving time step  0, time=0.0250...
  Refining the mesh...
 Solving time step  1, time=0.0500...
  Refining the mesh...
 Solving time step  2, time=0.0750...
  Refining the mesh...
 Solving time step  3, time=0.1000...
  Refining the mesh...
 Solving time step  4, time=0.1250...
  Refining the mesh...
 Solving time step  5, time=0.1500...
  Refining the mesh...
 Solving time step  6, time=0.1750...
  Refining the mesh...
 Solving time step  7, time=0.2000...
  Refining the mesh...
 Solving time step  8, time=0.2250...
  Refining the mesh...
 Solving time step  9, time=0.2500...
  Refining the mesh...
 Solving time step 10, time=0.2750...
  Refining the mesh...
 Solving time step 11, time=0.3000...
  Refining the mesh...
 Solving time step 12, time=0.3250...
  Refining the mesh...
 Solving time step 13, time=0.3500...
  Refining the mesh...
 Solving time step 14, time=0.3750...
  Refining the mesh...
 Solving time step 15, time=0.4000...
  Refining the mesh...
 Solving time step 16, time=0.4250...
  Refining the mesh...
 Solving time step 17, time=0.4500...
  Refining the mesh...
 Solving time step 18, time=0.4750...
  Refining the mesh...
 Solving time step 19, time=0.5000...
  Refining the mesh...
 Solving time step 20, time=0.5250...
  Refining the mesh...
 Solving time step 21, time=0.5500...
  Refining the mesh...
 Solving time step 22, time=0.5750...
  Refining the mesh...
 Solving time step 23, time=0.6000...
  Refining the mesh...
 Solving time step 24, time=0.6250...
  Refining the mesh...
 
***** Finished first 25 steps, now read in saved solution and continue *****
 
Usage:
	 ./ex10-devel -init_timestep 0
OR
	 ./ex10-devel -read_solution -init_timestep 26

Running: ./ex10-devel -read_solution -n_timesteps 25 -init_timestep 25

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=713
  n_elem()=1018
   n_local_elem()=1018
   n_active_elem()=766
  n_subdomains()=1
  n_processors()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE" 
    Approximation Orders="FIRST" 
    n_dofs()=713
    n_local_dofs()=713
    n_constrained_dofs()=122
    n_vectors()=3

 Solving time step 25, time=0.6500...
  Refining the mesh...
 Solving time step 26, time=0.6750...
  Refining the mesh...
 Solving time step 27, time=0.7000...
  Refining the mesh...
 Solving time step 28, time=0.7250...
  Refining the mesh...
 Solving time step 29, time=0.7500...
  Refining the mesh...
 Solving time step 30, time=0.7750...
  Refining the mesh...
 Solving time step 31, time=0.8000...
  Refining the mesh...
 Solving time step 32, time=0.8250...
  Refining the mesh...
 Solving time step 33, time=0.8500...
  Refining the mesh...
 Solving time step 34, time=0.8750...
  Refining the mesh...
 Solving time step 35, time=0.9000...
  Refining the mesh...
 Solving time step 36, time=0.9250...
  Refining the mesh...
 Solving time step 37, time=0.9500...
  Refining the mesh...
 Solving time step 38, time=0.9750...
  Refining the mesh...
 Solving time step 39, time=1.0000...
  Refining the mesh...
 Solving time step 40, time=1.0300...
  Refining the mesh...
 Solving time step 41, time=1.0500...
  Refining the mesh...
 Solving time step 42, time=1.0700...
  Refining the mesh...
 Solving time step 43, time=1.1000...
  Refining the mesh...
 Solving time step 44, time=1.1200...
  Refining the mesh...
 Solving time step 45, time=1.1500...
  Refining the mesh...
 Solving time step 46, time=1.1700...
  Refining the mesh...
 Solving time step 47, time=1.2000...
  Refining the mesh...
 Solving time step 48, time=1.2200...
  Refining the mesh...
 Solving time step 49, time=1.2500...
  Refining the mesh...
 
***************************************************************
* Done Running Example  ./ex10-devel
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
