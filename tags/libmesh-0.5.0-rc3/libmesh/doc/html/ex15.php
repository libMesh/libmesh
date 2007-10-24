<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex15",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 15 - Biharmonic Equation</h1>

<br><br>This example solves the Biharmonic equation on a square domain
using a Galerkin formulation with C1 elements approximating the
H^2_0 function space.
The initial mesh contains two TRI6 elements.
The mesh is provided in the standard libMesh ASCII format file
named "domain.xda".  In addition, an input file named "ex15.in"
is provided which allows the user to set several parameters for
the solution so that the problem can be re-run without a
re-compile.  The solution technique employed is to have a
refinement loop with a linear solve inside followed by a
refinement of the grid and projection of the solution to the new grid
In the final loop iteration, there is no additional
refinement after the solve.  In the input file "ex15.in", the variable
"max_r_steps" controls the number of refinement steps, and
"max_r_level" controls the maximum element refinement level.


<br><br>LibMesh include files.
</div>

<div class ="fragment">
<pre>
        #include "mesh.h"
        #include "equation_systems.h"
        #include "linear_implicit_system.h"
        #include "gmv_io.h"
        #include "tecplot_io.h"
        #include "fe.h"
        #include "quadrature_clough.h"
        #include "dense_matrix.h"
        #include "dense_vector.h"
        #include "sparse_matrix.h"
        #include "mesh_refinement.h"
        #include "error_vector.h"
        #include "kelly_error_estimator.h"
        #include "getpot.h"
        #include "exact_solution.h"
        
</pre>
</div>
<div class = "comment">
Function prototype.  This is the function that will assemble
the linear system for our Biharmonic problem.  Note that the
function will take the \p EquationSystems object and the
name of the system we are assembling as input.  From the
\p EquationSystems object we have acess to the \p Mesh and
other objects we might need.
</div>

<div class ="fragment">
<pre>
        void assemble_biharmonic(EquationSystems& es,
                              const std::string& system_name);
        
        
</pre>
</div>
<div class = "comment">
Prototype for calculation of the exact solution.  Useful
for setting boundary conditions.
</div>

<div class ="fragment">
<pre>
        Number exact_solution(const Point& p,
        		      const Real,          // time, not needed
        		      const std::string&,  // sys_name, not needed
        		      const std::string&); // unk_name, not needed);
        
</pre>
</div>
<div class = "comment">
Prototype for calculation of the gradient of the exact solution.  
Necessary for setting boundary conditions in H^2_0 and testing
H^1 convergence of the solution
</div>

<div class ="fragment">
<pre>
        Gradient exact_derivative(const Point& p,
        			  const Real,          // time, not needed
        			  const std::string&,  // sys_name, not needed
        			  const std::string&); // unk_name, not needed);
        
        Number forcing_function(const Point& p);
        
        
        
        int main(int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Initialize libMesh.
</div>

<div class ="fragment">
<pre>
          libMesh::init (argc, argv);
        
        #ifndef ENABLE_SECOND_DERIVATIVES
        
          std::cerr &lt;&lt; "ERROR: This example requires the library to be "
        	    &lt;&lt; "compiled with second derivatives support!"
        	    &lt;&lt; std::endl;
          here();
        
          return 0;
        
        #else
        
          {
</pre>
</div>
<div class = "comment">
Set the dimensionality of the mesh = 2
</div>

<div class ="fragment">
<pre>
            const unsigned int dim = 2;
        
</pre>
</div>
<div class = "comment">
Create a two-dimensional mesh.
</div>

<div class ="fragment">
<pre>
            Mesh mesh (dim);
            
</pre>
</div>
<div class = "comment">
Parse the input file
</div>

<div class ="fragment">
<pre>
            GetPot input_file("ex15.in");
        
</pre>
</div>
<div class = "comment">
Read in parameters from the input file
</div>

<div class ="fragment">
<pre>
            const unsigned int max_r_steps = input_file("max_r_steps", 3);
        
            std::cerr.setf(std::ios::scientific);
            std::cerr.precision(3);
            
</pre>
</div>
<div class = "comment">
Output file for plotting the error 
</div>

<div class ="fragment">
<pre>
            std::string output_file = "clough_uniform.m";
            
            std::ofstream out (output_file.c_str());
            out &lt;&lt; "% dofs     L2-error     H1-error" &lt;&lt; std::endl;
            out &lt;&lt; "e = [" &lt;&lt; std::endl;
            
</pre>
</div>
<div class = "comment">
Read in the mesh
</div>

<div class ="fragment">
<pre>
            mesh.read("domain.xda");
        
</pre>
</div>
<div class = "comment">
Mesh Refinement object
</div>

<div class ="fragment">
<pre>
            MeshRefinement mesh_refinement(mesh);
            
</pre>
</div>
<div class = "comment">
Create an equation systems object.
</div>

<div class ="fragment">
<pre>
            EquationSystems equation_systems (mesh);
        
</pre>
</div>
<div class = "comment">
Declare the system and its variables.
</div>

<div class ="fragment">
<pre>
            {
</pre>
</div>
<div class = "comment">
Creates a system named "Biharmonic"
</div>

<div class ="fragment">
<pre>
              equation_systems.add_system&lt;LinearImplicitSystem&gt; ("Biharmonic");
        
</pre>
</div>
<div class = "comment">
Adds the variable "u" to "Biharmonic".  "u"
will be approximated using Clough-Tocher cubic C1 triangles
</div>

<div class ="fragment">
<pre>
              equation_systems("Biharmonic").add_variable("u", THIRD, CLOUGH);
              
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly
function.
</div>

<div class ="fragment">
<pre>
              equation_systems("Biharmonic").attach_assemble_function
        		      (assemble_biharmonic);
              
</pre>
</div>
<div class = "comment">
Initialize the data structures for the equation system.
</div>

<div class ="fragment">
<pre>
              equation_systems.init();
        
</pre>
</div>
<div class = "comment">
Set linear solver max iterations
</div>

<div class ="fragment">
<pre>
              equation_systems.parameters.set&lt;unsigned int&gt;
        		      ("linear solver maximum iterations") = 1000;
        
</pre>
</div>
<div class = "comment">
Linear solver tolerance.
</div>

<div class ="fragment">
<pre>
              equation_systems.parameters.set&lt;Real&gt;
        		      ("linear solver tolerance") = 1.e-13;
              
</pre>
</div>
<div class = "comment">
Prints information about the system to the screen.
</div>

<div class ="fragment">
<pre>
              equation_systems.print_info();
            }
        
</pre>
</div>
<div class = "comment">
Construct ExactSolution object and attach function to compute exact solution
</div>

<div class ="fragment">
<pre>
            ExactSolution exact_sol(equation_systems);
            exact_sol.attach_exact_value(exact_solution);
            exact_sol.attach_exact_deriv(exact_derivative);
        
</pre>
</div>
<div class = "comment">
Convenient reference to the system
</div>

<div class ="fragment">
<pre>
            LinearImplicitSystem& system = 
              equation_systems.get_system&lt;LinearImplicitSystem&gt;("Biharmonic");
        
</pre>
</div>
<div class = "comment">
A refinement loop.
</div>

<div class ="fragment">
<pre>
            for (unsigned int r_step=0; r_step&lt;max_r_steps; r_step++)
              {
        	std::cout &lt;&lt; "Beginning Solve " &lt;&lt; r_step &lt;&lt; std::endl;
        	
</pre>
</div>
<div class = "comment">
Solve the system "Biharmonic", just like example 2.
</div>

<div class ="fragment">
<pre>
                equation_systems("Biharmonic").solve();
        
        	std::cout &lt;&lt; "System has: " &lt;&lt; equation_systems.n_active_dofs()
        		  &lt;&lt; " degrees of freedom."
        		  &lt;&lt; std::endl;
        	
        	std::cout &lt;&lt; "Linear solver converged at step: "
        		  &lt;&lt; system.n_linear_iterations()
        		  &lt;&lt; ", final residual: "
        		  &lt;&lt; system.final_linear_residual()
        		  &lt;&lt; std::endl;
        	
</pre>
</div>
<div class = "comment">
Compute the error.
</div>

<div class ="fragment">
<pre>
                exact_sol.compute_error("Biharmonic", "u");
        
</pre>
</div>
<div class = "comment">
Print out the error values
</div>

<div class ="fragment">
<pre>
                std::cout &lt;&lt; "L2-Error is: "
        		  &lt;&lt; exact_sol.l2_error("Biharmonic", "u")
        		  &lt;&lt; std::endl;
        	std::cout &lt;&lt; "H1-Error is: "
        		  &lt;&lt; exact_sol.h1_error("Biharmonic", "u")
        		  &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Print to output file
</div>

<div class ="fragment">
<pre>
                out &lt;&lt; equation_systems.n_active_dofs() &lt;&lt; " "
        	    &lt;&lt; exact_sol.l2_error("Biharmonic", "u") &lt;&lt; " "
        	    &lt;&lt; exact_sol.h1_error("Biharmonic", "u") &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Possibly refine the mesh - Clough Tocher elements currently
do not support hanging nodes, so we use uniform refinement
</div>

<div class ="fragment">
<pre>
                if (r_step+1 != max_r_steps)
        	  {
        	    std::cout &lt;&lt; "  Refining the mesh..." &lt;&lt; std::endl;
        
                    mesh_refinement.uniformly_refine(1);
        	    
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
Write out the solution
After solving the system write the solution
to a GMV-formatted plot file.
</div>

<div class ="fragment">
<pre>
            GMVIO (mesh).write_equation_systems ("domain.gmv",
            					 equation_systems);
        
</pre>
</div>
<div class = "comment">
Close up the output file.
</div>

<div class ="fragment">
<pre>
            out &lt;&lt; "];" &lt;&lt; std::endl;
            out &lt;&lt; "hold on" &lt;&lt; std::endl;
            out &lt;&lt; "plot(e(:,1), e(:,2), 'bo-');" &lt;&lt; std::endl;
            out &lt;&lt; "plot(e(:,1), e(:,3), 'ro-');" &lt;&lt; std::endl;
</pre>
</div>
<div class = "comment">
out << "set(gca,'XScale', 'Log');" << std::endl;
out << "set(gca,'YScale', 'Log');" << std::endl;
</div>

<div class ="fragment">
<pre>
            out &lt;&lt; "xlabel('dofs');" &lt;&lt; std::endl;
            out  &lt;&lt; "title('Clough-Tocher elements');" &lt;&lt; std::endl;
            out &lt;&lt; "legend('L2-error', 'H1-error');" &lt;&lt; std::endl;
</pre>
</div>
<div class = "comment">
out << "disp('L2-error linear fit');" << std::endl;
out << "polyfit(log10(e(:,1)), log10(e(:,2)), 1)" << std::endl;
out << "disp('H1-error linear fit');" << std::endl;
out << "polyfit(log10(e(:,1)), log10(e(:,3)), 1)" << std::endl;
</div>

<div class ="fragment">
<pre>
          }
        
          
</pre>
</div>
<div class = "comment">
All done.  
</div>

<div class ="fragment">
<pre>
          return libMesh::close ();
        #endif
        }
        
        
        
        
</pre>
</div>
<div class = "comment">
We now define the exact solution
</div>

<div class ="fragment">
<pre>
        Number exact_solution(const Point& p,
        		    const Real,         // time, not needed
        		    const std::string&, // sys_name, not needed
        		    const std::string&) // unk_name, not needed
        {
          const Real x = p(0);
          const Real y = p(1);
          
          return sin(x*y);
        }
        
        
        Number forcing_function(const Point& p)
        {
          const Real x = p(0);
          const Real y = p(1);
        
          return (x*x + y*y) * (x*x + y*y) * sin(x*y) -
        		  8*x*y*cos(x*y) - 4*sin(x*y);
        }
        
        
        
        
</pre>
</div>
<div class = "comment">
We now define the gradient of the exact solution
</div>

<div class ="fragment">
<pre>
        Gradient exact_derivative(const Point& p,
        			  const Real,         // time, not needed
        			  const std::string&, // sys_name, not needed
        			  const std::string&) // unk_name, not needed
        {
</pre>
</div>
<div class = "comment">
Gradient value to be returned.
</div>

<div class ="fragment">
<pre>
          Gradient gradu;
          
</pre>
</div>
<div class = "comment">
x and y coordinates in space
</div>

<div class ="fragment">
<pre>
          const Real x = p(0);
          const Real y = p(1);
        
          gradu(0) = y * cos(x*y);
          gradu(1) = x * cos(x*y);
        
          return gradu;
        }
        
        
        
        
        
        
</pre>
</div>
<div class = "comment">
We now define the matrix assembly function for the
Biharmonic system.  We need to first compute element
matrices and right-hand sides, and then take into
account the boundary conditions, which will be handled
via a penalty method.
</div>

<div class ="fragment">
<pre>
        void assemble_biharmonic(EquationSystems& es,
                              const std::string& system_name)
        {
</pre>
</div>
<div class = "comment">
It is a good idea to make sure we are assembling
the proper system.
</div>

<div class ="fragment">
<pre>
          assert (system_name == "Biharmonic");
        
        #ifdef ENABLE_SECOND_DERIVATIVES
        
</pre>
</div>
<div class = "comment">
Declare a performance log.  Give it a descriptive
string to identify what part of the code we are
logging, since there may be many PerfLogs in an
application.
</div>

<div class ="fragment">
<pre>
          PerfLog perf_log ("Matrix Assembly",false);
          
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
Get a reference to the LinearImplicitSystem we are solving
</div>

<div class ="fragment">
<pre>
          LinearImplicitSystem& system = es.get_system&lt;LinearImplicitSystem&gt;("Biharmonic");
          
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
Get a constant reference to the Finite Element type
for the first (and only) variable in the system.
</div>

<div class ="fragment">
<pre>
          FEType fe_type = dof_map.variable_type(0);
        
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
          AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
          
</pre>
</div>
<div class = "comment">
A 7th order Clough quadrature rule for numerical integration.
With 2D triangles, the Clough quadrature rule puts a Gaussian
quadrature rule on each of the 3 subelements
</div>

<div class ="fragment">
<pre>
          QClough qrule (dim, SEVENTH);
        
</pre>
</div>
<div class = "comment">
Tell the finite element object to use our quadrature rule.
</div>

<div class ="fragment">
<pre>
          fe-&gt;attach_quadrature_rule (&qrule);
        
</pre>
</div>
<div class = "comment">
Declare a special finite element object for
boundary integration.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
        	      
</pre>
</div>
<div class = "comment">
Boundary integration requires another quadraure rule,
with dimensionality one less than the dimensionality
of the element.
In 1D, the Clough and Gauss quadrature rules are identical.
</div>

<div class ="fragment">
<pre>
          QClough qface(dim-1, SIXTH);
          
</pre>
</div>
<div class = "comment">
Tell the finte element object to use our
quadrature rule.
</div>

<div class ="fragment">
<pre>
          fe_face-&gt;attach_quadrature_rule (&qface);
        
</pre>
</div>
<div class = "comment">
Here we define some references to cell-specific data that
will be used to assemble the linear system.
We begin with the element Jacobian * quadrature weight at each
integration point.   
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt;& JxW = fe-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The physical XY locations of the quadrature points on the element.
These might be useful for evaluating spatially varying material
properties at the quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Point&gt;& q_point = fe-&gt;get_xyz();
        
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
The element shape function second derivatives evaluated at the
quadrature points.  Note that for the simple biharmonic, shape
function first derivatives are unnecessary.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealTensor&gt; &gt;& d2phi = fe-&gt;get_d2phi();
        
</pre>
</div>
<div class = "comment">
Define data structures to contain the element matrix
and right-hand-side vector contribution.  Following
basic finite element terminology we will denote these
"Ke" and "Fe". More detail is in example 3.
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
Now we will loop over all the elements in the mesh.  We will
compute the element matrix and right-hand-side contribution.  See
example 3 for a discussion of the element iterators.


<br><br></div>

<div class ="fragment">
<pre>
          MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
          const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
          
          for ( ; el != end_el; ++el)
            {
</pre>
</div>
<div class = "comment">
Start logging the shape function initialization.
This is done through a simple function call with
the name of the event to log.
</div>

<div class ="fragment">
<pre>
              perf_log.start_event("elem init");      
        
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
summing them.
</div>

<div class ="fragment">
<pre>
              Ke.resize (dof_indices.size(),
        		 dof_indices.size());
        
              Fe.resize (dof_indices.size());
        
</pre>
</div>
<div class = "comment">
Stop logging the shape function initialization.
If you forget to stop logging an event the PerfLog
object will probably catch the error and abort.
</div>

<div class ="fragment">
<pre>
              perf_log.stop_event("elem init");      
        
</pre>
</div>
<div class = "comment">
Now we will build the element matrix.  This involves
a double loop to integrate laplacians of the test funcions
(i) against laplacians of the trial functions (j).

<br><br>This step is why we need the Clough-Tocher elements -
these C1 differentiable elements have square-integrable
second derivatives.

<br><br>Now start logging the element matrix computation
</div>

<div class ="fragment">
<pre>
              perf_log.start_event ("Ke");
        
              for (unsigned int qp=0; qp&lt;qrule.n_points(); qp++)
        	for (unsigned int i=0; i&lt;phi.size(); i++)
        	  for (unsigned int j=0; j&lt;phi.size(); j++)
        	    Ke(i,j) += JxW[qp]*(d2phi[i][qp](0,0)+d2phi[i][qp](1,1))
        			    *(d2phi[j][qp](0,0)+d2phi[j][qp](1,1));
        	    
        
</pre>
</div>
<div class = "comment">
Stop logging the matrix computation
</div>

<div class ="fragment">
<pre>
              perf_log.stop_event ("Ke");
        
        
</pre>
</div>
<div class = "comment">
At this point the interior element integration has
been completed.  However, we have not yet addressed
boundary conditions.  For this example we will only
consider simple Dirichlet boundary conditions imposed
via the penalty method.  Note that this is a fourth-order
problem: Dirichlet boundary conditions include *both*
boundary values and boundary normal fluxes.
</div>

<div class ="fragment">
<pre>
              {
</pre>
</div>
<div class = "comment">
Start logging the boundary condition computation
</div>

<div class ="fragment">
<pre>
                perf_log.start_event ("BCs");
        
</pre>
</div>
<div class = "comment">
The penalty value.  
</div>

<div class ="fragment">
<pre>
                const Real penalty = 1e10;
        	const Real penalty2 = 1e10;
        
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
</pre>
</div>
<div class = "comment">
The value of the shape functions at the quadrature
points.
</div>

<div class ="fragment">
<pre>
                      const std::vector&lt;std::vector&lt;Real&gt; &gt;&  phi_face =
        			      fe_face-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The value of the shape function derivatives at the
quadrature points.
</div>

<div class ="fragment">
<pre>
                      const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi_face =
        			      fe_face-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
The Jacobian * Quadrature Weight at the quadrature
points on the face.
</div>

<div class ="fragment">
<pre>
                      const std::vector&lt;Real&gt;& JxW_face = fe_face-&gt;get_JxW();
                                                                                       
</pre>
</div>
<div class = "comment">
The XYZ locations (in physical space) of the
quadrature points on the face.  This is where
we will interpolate the boundary value function.
</div>

<div class ="fragment">
<pre>
                      const std::vector&lt;Point &gt;& qface_point = fe_face-&gt;get_xyz();
        
        	      const std::vector&lt;Point&gt;& face_normals =
        			      fe_face-&gt;get_normals();
        
</pre>
</div>
<div class = "comment">
Compute the shape function values on the element
face.
</div>

<div class ="fragment">
<pre>
                      fe_face-&gt;reinit(elem, s);
                                                                                        
</pre>
</div>
<div class = "comment">
Loop over the face quagrature points for integration.
</div>

<div class ="fragment">
<pre>
                      for (unsigned int qp=0; qp&lt;qface.n_points(); qp++)
                        {
</pre>
</div>
<div class = "comment">
The boundary value.
</div>

<div class ="fragment">
<pre>
                          const Number value = exact_solution(qface_point[qp],
        						      0, "null",
        						      "void");
        		  const Gradient flux =
        				  exact_derivative(qface_point[qp], 0,
        						   "null", "void");
        
</pre>
</div>
<div class = "comment">
Matrix contribution of the L2 projection.
Note that the basis function values are
integrated against test function values while
basis fluxes are integrated against test function
fluxes.
</div>

<div class ="fragment">
<pre>
                          for (unsigned int i=0; i&lt;phi_face.size(); i++)
                            for (unsigned int j=0; j&lt;phi_face.size(); j++)
        		      Ke(i,j) += JxW_face[qp] *
        				 (penalty * phi_face[i][qp] *
        				  phi_face[j][qp] + penalty2
        				  * (dphi_face[i][qp] *
        				  face_normals[qp]) *
        				  (dphi_face[j][qp] *
        				   face_normals[qp]));
        
</pre>
</div>
<div class = "comment">
Right-hand-side contribution of the L2
projection.
</div>

<div class ="fragment">
<pre>
                          for (unsigned int i=0; i&lt;phi_face.size(); i++)
                            Fe(i) += JxW_face[qp] *
        				    (penalty * value * phi_face[i][qp]
        				     + penalty2 * 
        				     (flux * face_normals[qp])
        				    * (dphi_face[i][qp]
        				       * face_normals[qp]));
                        }
        	    } 
        	
</pre>
</div>
<div class = "comment">
Stop logging the boundary condition computation
</div>

<div class ="fragment">
<pre>
                perf_log.stop_event ("BCs");
              } 
        
              for (unsigned int qp=0; qp&lt;qrule.n_points(); qp++)
        	for (unsigned int i=0; i&lt;phi.size(); i++)
        	  Fe(i) += JxW[qp]*phi[i][qp]*forcing_function(q_point[qp]);
        	    
              
        
</pre>
</div>
<div class = "comment">
The element matrix and right-hand-side are now built
for this element.  Add them to the global matrix and
right-hand-side vector.  The \p PetscMatrix::add_matrix()
and \p PetscVector::add_vector() members do this for us.
Start logging the insertion of the local (element)
matrix and vector into the global matrix and vector
</div>

<div class ="fragment">
<pre>
              perf_log.start_event ("matrix insertion");
        
              dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
              system.matrix-&gt;add_matrix (Ke, dof_indices);
              system.rhs-&gt;add_vector    (Fe, dof_indices);
        
</pre>
</div>
<div class = "comment">
Stop logging the insertion of the local (element)
matrix and vector into the global matrix and vector
</div>

<div class ="fragment">
<pre>
              perf_log.stop_event ("matrix insertion");
            }
        
</pre>
</div>
<div class = "comment">
That's it.  We don't need to do anything else to the
PerfLog.  When it goes out of scope (at this function return)
it will print its log to the screen. Pretty easy, huh?


<br><br></div>

<div class ="fragment">
<pre>
          #else
        
        #endif
        
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;mesh.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;equation_systems.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;linear_implicit_system.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;gmv_io.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;tecplot_io.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;fe.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;quadrature_clough.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_matrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_vector.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;sparse_matrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;mesh_refinement.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;error_vector.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;kelly_error_estimator.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;getpot.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;exact_solution.h&quot;</FONT></B>
  
  <FONT COLOR="#228B22"><B>void</FONT></B> assemble_biharmonic(EquationSystems&amp; es,
                        <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name);
  
  
  Number exact_solution(<FONT COLOR="#228B22"><B>const</FONT></B> Point&amp; p,
  		      <FONT COLOR="#228B22"><B>const</FONT></B> Real,          <I><FONT COLOR="#B22222">// time, not needed
</FONT></I>  		      <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp;,  <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>  		      <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp;); <I><FONT COLOR="#B22222">// unk_name, not needed);
</FONT></I>  
  Gradient exact_derivative(<FONT COLOR="#228B22"><B>const</FONT></B> Point&amp; p,
  			  <FONT COLOR="#228B22"><B>const</FONT></B> Real,          <I><FONT COLOR="#B22222">// time, not needed
</FONT></I>  			  <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp;,  <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>  			  <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp;); <I><FONT COLOR="#B22222">// unk_name, not needed);
</FONT></I>  
  Number forcing_function(<FONT COLOR="#228B22"><B>const</FONT></B> Point&amp; p);
  
  
  
  <FONT COLOR="#228B22"><B>int</FONT></B> main(<FONT COLOR="#228B22"><B>int</FONT></B> argc, <FONT COLOR="#228B22"><B>char</FONT></B>** argv)
  {
    libMesh::init (argc, argv);
  
  #ifndef ENABLE_SECOND_DERIVATIVES
  
    std::cerr &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;ERROR: This example requires the library to be &quot;</FONT></B>
  	    &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;compiled with second derivatives support!&quot;</FONT></B>
  	    &lt;&lt; std::endl;
    here();
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    {
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = 2;
  
      Mesh mesh (dim);
      
      GetPot input_file(<FONT COLOR="#BC8F8F"><B>&quot;ex15.in&quot;</FONT></B>);
  
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> max_r_steps = input_file(<FONT COLOR="#BC8F8F"><B>&quot;max_r_steps&quot;</FONT></B>, 3);
  
      std::cerr.setf(std::ios::scientific);
      std::cerr.precision(3);
      
      std::string output_file = <FONT COLOR="#BC8F8F"><B>&quot;clough_uniform.m&quot;</FONT></B>;
      
      std::ofstream out (output_file.c_str());
      out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;% dofs     L2-error     H1-error&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;e = [&quot;</FONT></B> &lt;&lt; std::endl;
      
      mesh.read(<FONT COLOR="#BC8F8F"><B>&quot;domain.xda&quot;</FONT></B>);
  
      MeshRefinement mesh_refinement(mesh);
      
      EquationSystems equation_systems (mesh);
  
      {
        equation_systems.add_system&lt;LinearImplicitSystem&gt; (<FONT COLOR="#BC8F8F"><B>&quot;Biharmonic&quot;</FONT></B>);
  
        equation_systems(<FONT COLOR="#BC8F8F"><B>&quot;Biharmonic&quot;</FONT></B>).add_variable(<FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>, THIRD, CLOUGH);
        
        equation_systems(<FONT COLOR="#BC8F8F"><B>&quot;Biharmonic&quot;</FONT></B>).attach_assemble_function
  		      (assemble_biharmonic);
        
        equation_systems.init();
  
        equation_systems.parameters.set&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt;
  		      (<FONT COLOR="#BC8F8F"><B>&quot;linear solver maximum iterations&quot;</FONT></B>) = 1000;
  
        equation_systems.parameters.set&lt;Real&gt;
  		      (<FONT COLOR="#BC8F8F"><B>&quot;linear solver tolerance&quot;</FONT></B>) = 1.e-13;
        
        equation_systems.print_info();
      }
  
      ExactSolution exact_sol(equation_systems);
      exact_sol.attach_exact_value(exact_solution);
      exact_sol.attach_exact_deriv(exact_derivative);
  
      LinearImplicitSystem&amp; system = 
        equation_systems.get_system&lt;LinearImplicitSystem&gt;(<FONT COLOR="#BC8F8F"><B>&quot;Biharmonic&quot;</FONT></B>);
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> r_step=0; r_step&lt;max_r_steps; r_step++)
        {
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Beginning Solve &quot;</FONT></B> &lt;&lt; r_step &lt;&lt; std::endl;
  	
  	equation_systems(<FONT COLOR="#BC8F8F"><B>&quot;Biharmonic&quot;</FONT></B>).solve();
  
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;System has: &quot;</FONT></B> &lt;&lt; equation_systems.n_active_dofs()
  		  &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; degrees of freedom.&quot;</FONT></B>
  		  &lt;&lt; std::endl;
  	
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Linear solver converged at step: &quot;</FONT></B>
  		  &lt;&lt; system.n_linear_iterations()
  		  &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;, final residual: &quot;</FONT></B>
  		  &lt;&lt; system.final_linear_residual()
  		  &lt;&lt; std::endl;
  	
  	exact_sol.compute_error(<FONT COLOR="#BC8F8F"><B>&quot;Biharmonic&quot;</FONT></B>, <FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>);
  
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;L2-Error is: &quot;</FONT></B>
  		  &lt;&lt; exact_sol.l2_error(<FONT COLOR="#BC8F8F"><B>&quot;Biharmonic&quot;</FONT></B>, <FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>)
  		  &lt;&lt; std::endl;
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;H1-Error is: &quot;</FONT></B>
  		  &lt;&lt; exact_sol.h1_error(<FONT COLOR="#BC8F8F"><B>&quot;Biharmonic&quot;</FONT></B>, <FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>)
  		  &lt;&lt; std::endl;
  
  	out &lt;&lt; equation_systems.n_active_dofs() &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; &quot;</FONT></B>
  	    &lt;&lt; exact_sol.l2_error(<FONT COLOR="#BC8F8F"><B>&quot;Biharmonic&quot;</FONT></B>, <FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>) &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; &quot;</FONT></B>
  	    &lt;&lt; exact_sol.h1_error(<FONT COLOR="#BC8F8F"><B>&quot;Biharmonic&quot;</FONT></B>, <FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>) &lt;&lt; std::endl;
  
  	<B><FONT COLOR="#A020F0">if</FONT></B> (r_step+1 != max_r_steps)
  	  {
  	    std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;  Refining the mesh...&quot;</FONT></B> &lt;&lt; std::endl;
  
              mesh_refinement.uniformly_refine(1);
  	    
  	    equation_systems.reinit ();
  	  }
        }	    
      
      
  
      
      GMVIO (mesh).write_equation_systems (<FONT COLOR="#BC8F8F"><B>&quot;domain.gmv&quot;</FONT></B>,
      					 equation_systems);
  
      out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;];&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;hold on&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;plot(e(:,1), e(:,2), 'bo-');&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;plot(e(:,1), e(:,3), 'ro-');&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;xlabel('dofs');&quot;</FONT></B> &lt;&lt; std::endl;
      out  &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;title('Clough-Tocher elements');&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;legend('L2-error', 'H1-error');&quot;</FONT></B> &lt;&lt; std::endl;
    }
  
    
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close ();
  #endif
  }
  
  
  
  
  Number exact_solution(<FONT COLOR="#228B22"><B>const</FONT></B> Point&amp; p,
  		    <FONT COLOR="#228B22"><B>const</FONT></B> Real,         <I><FONT COLOR="#B22222">// time, not needed
</FONT></I>  		    <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp;, <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>  		    <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp;) <I><FONT COLOR="#B22222">// unk_name, not needed
</FONT></I>  {
    <FONT COLOR="#228B22"><B>const</FONT></B> Real x = p(0);
    <FONT COLOR="#228B22"><B>const</FONT></B> Real y = p(1);
    
    <B><FONT COLOR="#A020F0">return</FONT></B> sin(x*y);
  }
  
  
  Number forcing_function(<FONT COLOR="#228B22"><B>const</FONT></B> Point&amp; p)
  {
    <FONT COLOR="#228B22"><B>const</FONT></B> Real x = p(0);
    <FONT COLOR="#228B22"><B>const</FONT></B> Real y = p(1);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> (x*x + y*y) * (x*x + y*y) * sin(x*y) -
  		  8*x*y*cos(x*y) - 4*sin(x*y);
  }
  
  
  
  
  Gradient exact_derivative(<FONT COLOR="#228B22"><B>const</FONT></B> Point&amp; p,
  			  <FONT COLOR="#228B22"><B>const</FONT></B> Real,         <I><FONT COLOR="#B22222">// time, not needed
</FONT></I>  			  <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp;, <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>  			  <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp;) <I><FONT COLOR="#B22222">// unk_name, not needed
</FONT></I>  {
    Gradient gradu;
    
    <FONT COLOR="#228B22"><B>const</FONT></B> Real x = p(0);
    <FONT COLOR="#228B22"><B>const</FONT></B> Real y = p(1);
  
    gradu(0) = y * cos(x*y);
    gradu(1) = x * cos(x*y);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> gradu;
  }
  
  
  
  
  
  
  <FONT COLOR="#228B22"><B>void</FONT></B> assemble_biharmonic(EquationSystems&amp; es,
                        <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name)
  {
    assert (system_name == <FONT COLOR="#BC8F8F"><B>&quot;Biharmonic&quot;</FONT></B>);
  
  #ifdef ENABLE_SECOND_DERIVATIVES
  
    PerfLog perf_log (<FONT COLOR="#BC8F8F"><B>&quot;Matrix Assembly&quot;</FONT></B>,false);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> Mesh&amp; mesh = es.get_mesh();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = mesh.mesh_dimension();
  
    LinearImplicitSystem&amp; system = es.get_system&lt;LinearImplicitSystem&gt;(<FONT COLOR="#BC8F8F"><B>&quot;Biharmonic&quot;</FONT></B>);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
  
    FEType fe_type = dof_map.variable_type(0);
  
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    QClough qrule (dim, SEVENTH);
  
    fe-&gt;attach_quadrature_rule (&amp;qrule);
  
    AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
  	      
    QClough qface(dim-1, SIXTH);
    
    fe_face-&gt;attach_quadrature_rule (&amp;qface);
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Point&gt;&amp; q_point = fe-&gt;get_xyz();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;RealTensor&gt; &gt;&amp; d2phi = fe-&gt;get_d2phi();
  
    DenseMatrix&lt;Number&gt; Ke;
    DenseVector&lt;Number&gt; Fe;
  
    std::vector&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt; dof_indices;
  
  
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    <FONT COLOR="#228B22"><B>const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
    
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        perf_log.start_event(<FONT COLOR="#BC8F8F"><B>&quot;elem init&quot;</FONT></B>);      
  
        <FONT COLOR="#228B22"><B>const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
  
        fe-&gt;reinit (elem);
  
        Ke.resize (dof_indices.size(),
  		 dof_indices.size());
  
        Fe.resize (dof_indices.size());
  
        perf_log.stop_event(<FONT COLOR="#BC8F8F"><B>&quot;elem init&quot;</FONT></B>);      
  
        perf_log.start_event (<FONT COLOR="#BC8F8F"><B>&quot;Ke&quot;</FONT></B>);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi.size(); i++)
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j&lt;phi.size(); j++)
  	    Ke(i,j) += JxW[qp]*(d2phi[i][qp](0,0)+d2phi[i][qp](1,1))
  			    *(d2phi[j][qp](0,0)+d2phi[j][qp](1,1));
  	    
  
        perf_log.stop_event (<FONT COLOR="#BC8F8F"><B>&quot;Ke&quot;</FONT></B>);
  
  
        {
  	perf_log.start_event (<FONT COLOR="#BC8F8F"><B>&quot;BCs&quot;</FONT></B>);
  
  	<FONT COLOR="#228B22"><B>const</FONT></B> Real penalty = 1e10;
  	<FONT COLOR="#228B22"><B>const</FONT></B> Real penalty2 = 1e10;
  
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> s=0; s&lt;elem-&gt;n_sides(); s++)
  	  <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(s) == NULL)
  	    {
  	      <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp;  phi_face =
  			      fe_face-&gt;get_phi();
  
                <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi_face =
  			      fe_face-&gt;get_dphi();
  
                <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
                                                                                 
                <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Point &gt;&amp; qface_point = fe_face-&gt;get_xyz();
  
  	      <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Point&gt;&amp; face_normals =
  			      fe_face-&gt;get_normals();
  
                fe_face-&gt;reinit(elem, s);
                                                                                  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp&lt;qface.n_points(); qp++)
                  {
  		  <FONT COLOR="#228B22"><B>const</FONT></B> Number value = exact_solution(qface_point[qp],
  						      0, <FONT COLOR="#BC8F8F"><B>&quot;null&quot;</FONT></B>,
  						      <FONT COLOR="#BC8F8F"><B>&quot;void&quot;</FONT></B>);
  		  <FONT COLOR="#228B22"><B>const</FONT></B> Gradient flux =
  				  exact_derivative(qface_point[qp], 0,
  						   <FONT COLOR="#BC8F8F"><B>&quot;null&quot;</FONT></B>, <FONT COLOR="#BC8F8F"><B>&quot;void&quot;</FONT></B>);
  
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi_face.size(); i++)
                      <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j&lt;phi_face.size(); j++)
  		      Ke(i,j) += JxW_face[qp] *
  				 (penalty * phi_face[i][qp] *
  				  phi_face[j][qp] + penalty2
  				  * (dphi_face[i][qp] *
  				  face_normals[qp]) *
  				  (dphi_face[j][qp] *
  				   face_normals[qp]));
  
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi_face.size(); i++)
                      Fe(i) += JxW_face[qp] *
  				    (penalty * value * phi_face[i][qp]
  				     + penalty2 * 
  				     (flux * face_normals[qp])
  				    * (dphi_face[i][qp]
  				       * face_normals[qp]));
                  }
  	    } 
  	
  	perf_log.stop_event (<FONT COLOR="#BC8F8F"><B>&quot;BCs&quot;</FONT></B>);
        } 
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi.size(); i++)
  	  Fe(i) += JxW[qp]*phi[i][qp]*forcing_function(q_point[qp]);
  	    
        
  
        perf_log.start_event (<FONT COLOR="#BC8F8F"><B>&quot;matrix insertion&quot;</FONT></B>);
  
        dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
        system.matrix-&gt;add_matrix (Ke, dof_indices);
        system.rhs-&gt;add_vector    (Fe, dof_indices);
  
        perf_log.stop_event (<FONT COLOR="#BC8F8F"><B>&quot;matrix insertion&quot;</FONT></B>);
      }
  
  
    #<B><FONT COLOR="#A020F0">else</FONT></B>
  
  #endif
  
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example  ./ex15
***************************************************************
 
 EquationSystems
  n_systems()=1
   System "Biharmonic"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="21", "12" 
    Infinite Element Mapping="0" 
    Approximation Orders="3", "3" 
    n_dofs()=12
    n_local_dofs()=12
    n_constrained_dofs()=0
    n_vectors()=1

Beginning Solve 0
System has: 12 degrees of freedom.
Linear solver converged at step: 1, final residual: 1.03663e-15
L2-Error is: 0.000246843
H1-Error is: 0.00268853
  Refining the mesh...
Beginning Solve 1
System has: 27 degrees of freedom.
Linear solver converged at step: 9, final residual: 2.04817e-14
L2-Error is: 7.83194e-05
H1-Error is: 0.000799048
  Refining the mesh...
Beginning Solve 2
System has: 75 degrees of freedom.
Linear solver converged at step: 15, final residual: 1.21094e-13
L2-Error is: 9.00289e-06
H1-Error is: 0.000151356

 ---------------------------------------------------------------------------- 
| Reference count information                                                |
 ---------------------------------------------------------------------------- 
| 12LinearSolverIdE reference count information:
|  Creations:    1
|  Destructions: 1
| 12SparseMatrixIdE reference count information:
|  Creations:    1
|  Destructions: 1
| 13NumericVectorIdE reference count information:
|  Creations:    5
|  Destructions: 5
| 4Elem reference count information:
|  Creations:    90
|  Destructions: 90
| 4Node reference count information:
|  Creations:    45
|  Destructions: 45
| 5QBase reference count information:
|  Creations:    17
|  Destructions: 17
| 6DofMap reference count information:
|  Creations:    1
|  Destructions: 1
| 6FEBase reference count information:
|  Creations:    11
|  Destructions: 11
| 6System reference count information:
|  Creations:    1
|  Destructions: 1
| N10Parameters5ValueE reference count information:
|  Creations:    2
|  Destructions: 2
 ---------------------------------------------------------------------------- 
 
***************************************************************
* Done Running Example  ./ex15
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
