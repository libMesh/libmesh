<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex14",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 14 - Laplace Equation in the L-Shaped Domain</h1>

<br><br>This example solves the Laplace equation on the classic "L-shaped"
domain with adaptive mesh refinement.  In this case, the exact
solution is u(r,\theta) = r^{2/3} * \sin ( (2/3) * \theta), but
the standard Kelly error indicator is used to estimate the error.
The initial mesh contains three QUAD9 elements which represent the
standard quadrants I, II, and III of the domain [-1,1]x[-1,1],
i.e.
Element 0: [-1,0]x[ 0,1]
Element 1: [ 0,1]x[ 0,1]
Element 2: [-1,0]x[-1,0]
The mesh is provided in the standard libMesh ASCII format file
named "lshaped.xda".  In addition, an input file named "ex14.in"
is provided which allows the user to set several parameters for
the solution so that the problem can be re-run without a
re-compile.  The solution technique employed is to have a
refinement loop with a linear solve inside followed by a
refinement of the grid and projection of the solution to the new grid
In the final loop iteration, there is no additional
refinement after the solve.  In the input file "ex14.in", the variable
"max_r_steps" controls the number of refinement steps,
"max_r_level" controls the maximum element refinement level, and
"refine_percentage" / "coarsen_percentage" determine the number of
elements which will be refined / coarsened at each step.


<br><br>LibMesh include files.
</div>

<div class ="fragment">
<pre>
        #include "mesh.h"
        #include "equation_systems.h"
        #include "implicit_system.h"
        #include "gmv_io.h"
        #include "tecplot_io.h"
        #include "fe.h"
        #include "quadrature_gauss.h"
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
the linear system for our Laplace problem.  Note that the
function will take the \p EquationSystems object and the
name of the system we are assembling as input.  From the
\p EquationSystems object we have acess to the \p Mesh and
other objects we might need.
</div>

<div class ="fragment">
<pre>
        void assemble_laplace(EquationSystems& es,
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
</div>

<div class ="fragment">
<pre>
        Gradient exact_derivative(const Point& p,
        			  const Real,          // time, not needed
        			  const std::string&,  // sys_name, not needed
        			  const std::string&); // unk_name, not needed);
        
        
        
        
        
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
            GetPot input_file("ex14.in");
        
</pre>
</div>
<div class = "comment">
Read in parameters from the input file
</div>

<div class ="fragment">
<pre>
            const unsigned int max_r_steps = input_file("max_r_steps", 3);
            const unsigned int max_r_level = input_file("max_r_level", 3);
            const Real refine_percentage   = input_file("refine_percentage", 0.5);
            const Real coarsen_percentage  = input_file("coarsen_percentage", 0.5);
            const unsigned int uniform_refine = input_file("uniform_refine",0);
            const std::string approx_order    = input_file("approx_order", "FIRST");
            
</pre>
</div>
<div class = "comment">
Output file for plotting the error as a function of
the number of degrees of freedom.
</div>

<div class ="fragment">
<pre>
            std::string output_file = "bi";
            if (approx_order == "FIRST")
              output_file += "linear_";
            else
              output_file += "quadratic_";
        
            if (uniform_refine == 0)
              output_file += "adaptive.m";
            else
              output_file += "uniform.m";
            
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
            mesh.read("lshaped.xda");
        
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
Creates a system named "Laplace"
</div>

<div class ="fragment">
<pre>
              equation_systems.add_system&lt;ImplicitSystem&gt; ("Laplace");
        
</pre>
</div>
<div class = "comment">
Adds the variable "u" to "Laplace".  "u"
will be approximated using second-order approximation.
</div>

<div class ="fragment">
<pre>
              if (approx_order == "FIRST")
        	equation_systems("Laplace").add_variable("u", FIRST);
              else
        	equation_systems("Laplace").add_variable("u", SECOND);
              
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly
function.
</div>

<div class ="fragment">
<pre>
              equation_systems("Laplace").attach_assemble_function (assemble_laplace);
              
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
              equation_systems.set_parameter("linear solver maximum iterations") = 100;
        
</pre>
</div>
<div class = "comment">
Linear solver tolerance.
</div>

<div class ="fragment">
<pre>
              equation_systems.set_parameter("linear solver tolerance") = 1.e-12;
              
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
            ImplicitSystem& system = equation_systems.get_system&lt;ImplicitSystem&gt;("Laplace");
        
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
Solve the system "Laplace", just like example 2.
</div>

<div class ="fragment">
<pre>
                equation_systems("Laplace").solve();
        
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
                exact_sol.compute_error("Laplace", "u");
        
</pre>
</div>
<div class = "comment">
Print out the error values
</div>

<div class ="fragment">
<pre>
                std::cout &lt;&lt; "L2-Error is: "
        		  &lt;&lt; exact_sol.l2_error("Laplace", "u")
        		  &lt;&lt; std::endl;
        	std::cout &lt;&lt; "H1-Error is: "
        		  &lt;&lt; exact_sol.h1_error("Laplace", "u")
        		  &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Print to output file
</div>

<div class ="fragment">
<pre>
                out &lt;&lt; equation_systems.n_active_dofs() &lt;&lt; " "
        	    &lt;&lt; exact_sol.l2_error("Laplace", "u") &lt;&lt; " "
        	    &lt;&lt; exact_sol.h1_error("Laplace", "u") &lt;&lt; std::endl;
        
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
        
        	    if (uniform_refine == 0)
                      {
        
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
element within 10% of the minimum error on any element might
be coarsened. Note that the elements flagged for refinement
will be refined, but those flagged for coarsening _might_ be
coarsened.
</div>

<div class ="fragment">
<pre>
                        mesh_refinement.flag_elements_by_error_fraction (error,
        								 refine_percentage,
        								 coarsen_percentage,
        								 max_r_level);
        		
</pre>
</div>
<div class = "comment">
This call actually refines and coarsens the flagged
elements.
</div>

<div class ="fragment">
<pre>
                        mesh_refinement.refine_and_coarsen_elements();
        	      }
        
        	    else if (uniform_refine == 1)
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
            GMVIO (mesh).write_equation_systems ("lshaped.gmv",
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
            if (approx_order == "FIRST")
              out &lt;&lt; "title('Bilinear elements');" &lt;&lt; std::endl;
            else
              out  &lt;&lt; "title('Biquadratic elements');" &lt;&lt; std::endl;
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
        }
        
        
        
        
</pre>
</div>
<div class = "comment">
We now define the exact solution, being careful
to obtain an angle from atan2 in the correct
quadrant.
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
          
</pre>
</div>
<div class = "comment">
The boundary value, given by the exact solution,
u_exact = r^(2/3)*sin(2*theta/3).
</div>

<div class ="fragment">
<pre>
          Real theta = atan2(y,x);
        
</pre>
</div>
<div class = "comment">
Make sure 0 <= theta <= 2*pi
</div>

<div class ="fragment">
<pre>
          if (theta &lt; 0)
            theta += 2. * libMesh::pi;
        		  
          return pow(x*x + y*y, 1./3.)*sin(2./3.*theta);
        }
        
        
        
        
        
</pre>
</div>
<div class = "comment">
We now define the gradient of the exact solution, again being careful
to obtain an angle from atan2 in the correct
quadrant.
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
        
</pre>
</div>
<div class = "comment">
We can't compute the gradient at x=0, it is not defined.
</div>

<div class ="fragment">
<pre>
          assert (x != 0.);
        
</pre>
</div>
<div class = "comment">
For convenience...
</div>

<div class ="fragment">
<pre>
          const Real tt = 2./3.;
          const Real ot = 1./3.;
          
</pre>
</div>
<div class = "comment">
The value of the radius, squared
</div>

<div class ="fragment">
<pre>
          const Real r2 = x*x + y*y;
        
</pre>
</div>
<div class = "comment">
The boundary value, given by the exact solution,
u_exact = r^(2/3)*sin(2*theta/3).
</div>

<div class ="fragment">
<pre>
          Real theta = atan2(y,x);
          
</pre>
</div>
<div class = "comment">
Make sure 0 <= theta <= 2*pi
</div>

<div class ="fragment">
<pre>
          if (theta &lt; 0)
            theta += 2. * libMesh::pi;
        
</pre>
</div>
<div class = "comment">
du/dx
</div>

<div class ="fragment">
<pre>
          gradu(0) = tt*x*pow(r2,-tt)*sin(tt*theta) - pow(r2,ot)*cos(tt*theta)*tt/(1.+y*y/x/x)*y/x/x;
        
</pre>
</div>
<div class = "comment">
du/dy
</div>

<div class ="fragment">
<pre>
          gradu(1) = tt*y*pow(r2,-tt)*sin(tt*theta) + pow(r2,ot)*cos(tt*theta)*tt/(1.+y*y/x/x)*1./x;
            
          return gradu;
        }
        
        
        
        
        
        
</pre>
</div>
<div class = "comment">
We now define the matrix assembly function for the
Laplace system.  We need to first compute element
matrices and right-hand sides, and then take into
account the boundary conditions, which will be handled
via a penalty method.
</div>

<div class ="fragment">
<pre>
        void assemble_laplace(EquationSystems& es,
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
          assert (system_name == "Laplace");
        
        
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
Get a reference to the ImplicitSystem we are solving
</div>

<div class ="fragment">
<pre>
          ImplicitSystem& system = es.get_system&lt;ImplicitSystem&gt;("Laplace");
          
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
A 5th order Gauss quadrature rule for numerical integration.
</div>

<div class ="fragment">
<pre>
          QGauss qrule (dim, FIFTH);
        
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
Boundary integration requires one quadraure rule,
with dimensionality one less than the dimensionality
of the element.
</div>

<div class ="fragment">
<pre>
          QGauss qface(dim-1, FIFTH);
          
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
const std::vector<Point>& q_point = fe->get_xyz();


<br><br>The element shape functions evaluated at the quadrature points.
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
example 3 for a discussion of the element iterators.  Here we use
the \p const_active_local_elem_iterator to indicate we only want
to loop over elements that are assigned to the local processor
which are "active" in the sense of AMR.  This allows each
processor to compute its components of the global matrix for
active elements while ignoring parent elements which have been
refined.
const_active_local_elem_iterator           el (mesh.elements_begin());
const const_active_local_elem_iterator end_el (mesh.elements_end());


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
a double loop to integrate the test funcions (i) against
the trial functions (j).

<br><br>Now start logging the element matrix computation
</div>

<div class ="fragment">
<pre>
              perf_log.start_event ("Ke");
        
              for (unsigned int qp=0; qp&lt;qrule.n_points(); qp++)
        	for (unsigned int i=0; i&lt;phi.size(); i++)
        	  for (unsigned int j=0; j&lt;phi.size(); j++)
        	    Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
        	    
        
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
via the penalty method. This is discussed at length in
example 3.
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
</pre>
</div>
<div class = "comment">
Get a pointer to the side.
</div>

<div class ="fragment">
<pre>
                      AutoPtr&lt;Elem&gt; side (elem-&gt;build_side(s));
        	      
</pre>
</div>
<div class = "comment">
Loop over the nodes on the side with NULL neighbor.
</div>

<div class ="fragment">
<pre>
                      for (unsigned int ns=0; ns&lt;side-&gt;n_nodes(); ns++)
        		{
</pre>
</div>
<div class = "comment">
const Real x = side->point(ns)(0);
const Real y = side->point(ns)(1);
</div>

<div class ="fragment">
<pre>
                          const Number value = exact_solution(side-&gt;point(ns),
        						      0.,
        						      "null",
        						      "void");
        		  
</pre>
</div>
<div class = "comment">
std::cout << "(x,y,bc)=("
<< x << ","
<< y << ","
<< value << ")"
<< std::endl;
		  

<br><br>Find the node on the element matching this node on
the side.  That defines where in the element matrix
the boundary condition will be applied.
</div>

<div class ="fragment">
<pre>
                          for (unsigned int n=0; n&lt;dof_indices.size(); n++)
        		    if (elem-&gt;node(n) == side-&gt;node(ns))
        		      {
</pre>
</div>
<div class = "comment">
Matrix contribution.
</div>

<div class ="fragment">
<pre>
                                Ke(n,n) += penalty;
        			
</pre>
</div>
<div class = "comment">
Right-hand-side contribution.
</div>

<div class ="fragment">
<pre>
                                Fe(n) += penalty*value;
        		      }
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
Start logging the insertion of the local (element)
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
</div>

<div class ="fragment">
<pre>
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;mesh.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;equation_systems.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;implicit_system.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;gmv_io.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;tecplot_io.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;fe.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_matrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_vector.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;sparse_matrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;mesh_refinement.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;error_vector.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;kelly_error_estimator.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;getpot.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;exact_solution.h&quot;</FONT></B>
  
  <FONT COLOR="#228B22"><B>void</FONT></B> assemble_laplace(EquationSystems&amp; es,
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
  
  
  
  
  <FONT COLOR="#228B22"><B>int</FONT></B> main(<FONT COLOR="#228B22"><B>int</FONT></B> argc, <FONT COLOR="#228B22"><B>char</FONT></B>** argv)
  {
    libMesh::init (argc, argv);
  
    {
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = 2;
  
      Mesh mesh (dim);
      
      GetPot input_file(<FONT COLOR="#BC8F8F"><B>&quot;ex14.in&quot;</FONT></B>);
  
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> max_r_steps = input_file(<FONT COLOR="#BC8F8F"><B>&quot;max_r_steps&quot;</FONT></B>, 3);
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> max_r_level = input_file(<FONT COLOR="#BC8F8F"><B>&quot;max_r_level&quot;</FONT></B>, 3);
      <FONT COLOR="#228B22"><B>const</FONT></B> Real refine_percentage   = input_file(<FONT COLOR="#BC8F8F"><B>&quot;refine_percentage&quot;</FONT></B>, 0.5);
      <FONT COLOR="#228B22"><B>const</FONT></B> Real coarsen_percentage  = input_file(<FONT COLOR="#BC8F8F"><B>&quot;coarsen_percentage&quot;</FONT></B>, 0.5);
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> uniform_refine = input_file(<FONT COLOR="#BC8F8F"><B>&quot;uniform_refine&quot;</FONT></B>,0);
      <FONT COLOR="#228B22"><B>const</FONT></B> std::string approx_order    = input_file(<FONT COLOR="#BC8F8F"><B>&quot;approx_order&quot;</FONT></B>, <FONT COLOR="#BC8F8F"><B>&quot;FIRST&quot;</FONT></B>);
      
      std::string output_file = <FONT COLOR="#BC8F8F"><B>&quot;bi&quot;</FONT></B>;
      <B><FONT COLOR="#A020F0">if</FONT></B> (approx_order == <FONT COLOR="#BC8F8F"><B>&quot;FIRST&quot;</FONT></B>)
        output_file += <FONT COLOR="#BC8F8F"><B>&quot;linear_&quot;</FONT></B>;
      <B><FONT COLOR="#A020F0">else</FONT></B>
        output_file += <FONT COLOR="#BC8F8F"><B>&quot;quadratic_&quot;</FONT></B>;
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (uniform_refine == 0)
        output_file += <FONT COLOR="#BC8F8F"><B>&quot;adaptive.m&quot;</FONT></B>;
      <B><FONT COLOR="#A020F0">else</FONT></B>
        output_file += <FONT COLOR="#BC8F8F"><B>&quot;uniform.m&quot;</FONT></B>;
      
      std::ofstream out (output_file.c_str());
      out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;% dofs     L2-error     H1-error&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;e = [&quot;</FONT></B> &lt;&lt; std::endl;
      
      mesh.read(<FONT COLOR="#BC8F8F"><B>&quot;lshaped.xda&quot;</FONT></B>);
  
      MeshRefinement mesh_refinement(mesh);
      
      EquationSystems equation_systems (mesh);
  
      {
        equation_systems.add_system&lt;ImplicitSystem&gt; (<FONT COLOR="#BC8F8F"><B>&quot;Laplace&quot;</FONT></B>);
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (approx_order == <FONT COLOR="#BC8F8F"><B>&quot;FIRST&quot;</FONT></B>)
  	equation_systems(<FONT COLOR="#BC8F8F"><B>&quot;Laplace&quot;</FONT></B>).add_variable(<FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>, FIRST);
        <B><FONT COLOR="#A020F0">else</FONT></B>
  	equation_systems(<FONT COLOR="#BC8F8F"><B>&quot;Laplace&quot;</FONT></B>).add_variable(<FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>, SECOND);
        
        equation_systems(<FONT COLOR="#BC8F8F"><B>&quot;Laplace&quot;</FONT></B>).attach_assemble_function (assemble_laplace);
        
        equation_systems.init();
  
        equation_systems.set_parameter(<FONT COLOR="#BC8F8F"><B>&quot;linear solver maximum iterations&quot;</FONT></B>) = 100;
  
        equation_systems.set_parameter(<FONT COLOR="#BC8F8F"><B>&quot;linear solver tolerance&quot;</FONT></B>) = 1.e-12;
        
        equation_systems.print_info();
      }
  
      ExactSolution exact_sol(equation_systems);
      exact_sol.attach_exact_value(exact_solution);
      exact_sol.attach_exact_deriv(exact_derivative);
  
      ImplicitSystem&amp; system = equation_systems.get_system&lt;ImplicitSystem&gt;(<FONT COLOR="#BC8F8F"><B>&quot;Laplace&quot;</FONT></B>);
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> r_step=0; r_step&lt;max_r_steps; r_step++)
        {
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Beginning Solve &quot;</FONT></B> &lt;&lt; r_step &lt;&lt; std::endl;
  	
  	equation_systems(<FONT COLOR="#BC8F8F"><B>&quot;Laplace&quot;</FONT></B>).solve();
  
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;System has: &quot;</FONT></B> &lt;&lt; equation_systems.n_active_dofs()
  		  &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; degrees of freedom.&quot;</FONT></B>
  		  &lt;&lt; std::endl;
  	
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Linear solver converged at step: &quot;</FONT></B>
  		  &lt;&lt; system.n_linear_iterations()
  		  &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;, final residual: &quot;</FONT></B>
  		  &lt;&lt; system.final_linear_residual()
  		  &lt;&lt; std::endl;
  	
  	exact_sol.compute_error(<FONT COLOR="#BC8F8F"><B>&quot;Laplace&quot;</FONT></B>, <FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>);
  
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;L2-Error is: &quot;</FONT></B>
  		  &lt;&lt; exact_sol.l2_error(<FONT COLOR="#BC8F8F"><B>&quot;Laplace&quot;</FONT></B>, <FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>)
  		  &lt;&lt; std::endl;
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;H1-Error is: &quot;</FONT></B>
  		  &lt;&lt; exact_sol.h1_error(<FONT COLOR="#BC8F8F"><B>&quot;Laplace&quot;</FONT></B>, <FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>)
  		  &lt;&lt; std::endl;
  
  	out &lt;&lt; equation_systems.n_active_dofs() &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; &quot;</FONT></B>
  	    &lt;&lt; exact_sol.l2_error(<FONT COLOR="#BC8F8F"><B>&quot;Laplace&quot;</FONT></B>, <FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>) &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; &quot;</FONT></B>
  	    &lt;&lt; exact_sol.h1_error(<FONT COLOR="#BC8F8F"><B>&quot;Laplace&quot;</FONT></B>, <FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>) &lt;&lt; std::endl;
  
  	<B><FONT COLOR="#A020F0">if</FONT></B> (r_step+1 != max_r_steps)
  	  {
  	    std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;  Refining the mesh...&quot;</FONT></B> &lt;&lt; std::endl;
  
  	    <B><FONT COLOR="#A020F0">if</FONT></B> (uniform_refine == 0)
                {
  
  		ErrorVector error;
  		
  		KellyErrorEstimator error_estimator;
  		
  		error_estimator.estimate_error (system,
  						error);
  		
  		mesh_refinement.flag_elements_by_error_fraction (error,
  								 refine_percentage,
  								 coarsen_percentage,
  								 max_r_level);
  		
  		mesh_refinement.refine_and_coarsen_elements();
  	      }
  
  	    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (uniform_refine == 1)
                mesh_refinement.uniformly_refine(1);
  	    
  	    equation_systems.reinit ();
  	  }
        }	    
      
      
  
      
      GMVIO (mesh).write_equation_systems (<FONT COLOR="#BC8F8F"><B>&quot;lshaped.gmv&quot;</FONT></B>,
      					 equation_systems);
  
      out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;];&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;hold on&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;plot(e(:,1), e(:,2), 'bo-');&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;plot(e(:,1), e(:,3), 'ro-');&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;xlabel('dofs');&quot;</FONT></B> &lt;&lt; std::endl;
      <B><FONT COLOR="#A020F0">if</FONT></B> (approx_order == <FONT COLOR="#BC8F8F"><B>&quot;FIRST&quot;</FONT></B>)
        out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;title('Bilinear elements');&quot;</FONT></B> &lt;&lt; std::endl;
      <B><FONT COLOR="#A020F0">else</FONT></B>
        out  &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;title('Biquadratic elements');&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;legend('L2-error', 'H1-error');&quot;</FONT></B> &lt;&lt; std::endl;
    }
  
    
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close ();
  }
  
  
  
  
  Number exact_solution(<FONT COLOR="#228B22"><B>const</FONT></B> Point&amp; p,
  		    <FONT COLOR="#228B22"><B>const</FONT></B> Real,         <I><FONT COLOR="#B22222">// time, not needed
</FONT></I>  		    <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp;, <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>  		    <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp;) <I><FONT COLOR="#B22222">// unk_name, not needed
</FONT></I>  {
    <FONT COLOR="#228B22"><B>const</FONT></B> Real x = p(0);
    <FONT COLOR="#228B22"><B>const</FONT></B> Real y = p(1);
    
    Real theta = atan2(y,x);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (theta &lt; 0)
      theta += 2. * libMesh::pi;
  		  
    <B><FONT COLOR="#A020F0">return</FONT></B> pow(x*x + y*y, 1./3.)*sin(2./3.*theta);
  }
  
  
  
  
  
  Gradient exact_derivative(<FONT COLOR="#228B22"><B>const</FONT></B> Point&amp; p,
  			  <FONT COLOR="#228B22"><B>const</FONT></B> Real,         <I><FONT COLOR="#B22222">// time, not needed
</FONT></I>  			  <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp;, <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>  			  <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp;) <I><FONT COLOR="#B22222">// unk_name, not needed
</FONT></I>  {
    Gradient gradu;
    
    <FONT COLOR="#228B22"><B>const</FONT></B> Real x = p(0);
    <FONT COLOR="#228B22"><B>const</FONT></B> Real y = p(1);
  
    assert (x != 0.);
  
    <FONT COLOR="#228B22"><B>const</FONT></B> Real tt = 2./3.;
    <FONT COLOR="#228B22"><B>const</FONT></B> Real ot = 1./3.;
    
    <FONT COLOR="#228B22"><B>const</FONT></B> Real r2 = x*x + y*y;
  
    Real theta = atan2(y,x);
    
    <B><FONT COLOR="#A020F0">if</FONT></B> (theta &lt; 0)
      theta += 2. * libMesh::pi;
  
    gradu(0) = tt*x*pow(r2,-tt)*sin(tt*theta) - pow(r2,ot)*cos(tt*theta)*tt/(1.+y*y/x/x)*y/x/x;
  
    gradu(1) = tt*y*pow(r2,-tt)*sin(tt*theta) + pow(r2,ot)*cos(tt*theta)*tt/(1.+y*y/x/x)*1./x;
      
    <B><FONT COLOR="#A020F0">return</FONT></B> gradu;
  }
  
  
  
  
  
  
  <FONT COLOR="#228B22"><B>void</FONT></B> assemble_laplace(EquationSystems&amp; es,
                        <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name)
  {
    assert (system_name == <FONT COLOR="#BC8F8F"><B>&quot;Laplace&quot;</FONT></B>);
  
  
    PerfLog perf_log (<FONT COLOR="#BC8F8F"><B>&quot;Matrix Assembly&quot;</FONT></B>,false);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> Mesh&amp; mesh = es.get_mesh();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = mesh.mesh_dimension();
  
    ImplicitSystem&amp; system = es.get_system&lt;ImplicitSystem&gt;(<FONT COLOR="#BC8F8F"><B>&quot;Laplace&quot;</FONT></B>);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
  
    FEType fe_type = dof_map.variable_type(0);
  
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    QGauss qrule (dim, FIFTH);
  
    fe-&gt;attach_quadrature_rule (&amp;qrule);
  
    AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
  	      
    QGauss qface(dim-1, FIFTH);
    
    fe_face-&gt;attach_quadrature_rule (&amp;qface);
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
  
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
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
  	    Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
  	    
  
        perf_log.stop_event (<FONT COLOR="#BC8F8F"><B>&quot;Ke&quot;</FONT></B>);
  
  
        {
  	perf_log.start_event (<FONT COLOR="#BC8F8F"><B>&quot;BCs&quot;</FONT></B>);
  
  	<FONT COLOR="#228B22"><B>const</FONT></B> Real penalty = 1.e10;
  
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> s=0; s&lt;elem-&gt;n_sides(); s++)
  	  <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(s) == NULL)
  	    {
  	      AutoPtr&lt;Elem&gt; side (elem-&gt;build_side(s));
  	      
  	      <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> ns=0; ns&lt;side-&gt;n_nodes(); ns++)
  		{
  		  <FONT COLOR="#228B22"><B>const</FONT></B> Number value = exact_solution(side-&gt;point(ns),
  						      0.,
  						      <FONT COLOR="#BC8F8F"><B>&quot;null&quot;</FONT></B>,
  						      <FONT COLOR="#BC8F8F"><B>&quot;void&quot;</FONT></B>);
  		  
  		  
  		  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n=0; n&lt;dof_indices.size(); n++)
  		    <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;node(n) == side-&gt;node(ns))
  		      {
  			Ke(n,n) += penalty;
  			
  			Fe(n) += penalty*value;
  		      }
  		} 
  	    } 
  	
  	perf_log.stop_event (<FONT COLOR="#BC8F8F"><B>&quot;BCs&quot;</FONT></B>);
        } 
        
  
        perf_log.start_event (<FONT COLOR="#BC8F8F"><B>&quot;matrix insertion&quot;</FONT></B>);
  
        dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
        system.matrix-&gt;add_matrix (Ke, dof_indices);
        system.rhs-&gt;add_vector    (Fe, dof_indices);
  
        perf_log.stop_event (<FONT COLOR="#BC8F8F"><B>&quot;matrix insertion&quot;</FONT></B>);
      }
  
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Compiling C++ (in debug mode) ex14.C...
Linking ex14...
/home/peterson/code/libmesh/contrib/tecplot/lib/i686-pc-linux-gnu/tecio.a(tecxxx.o)(.text+0x1a7): In function `tecini':
: the use of `mktemp' is dangerous, better use `mkstemp'
/usr/bin/ld: warning: libstdc++.so.5, needed by /usr/local/petsc/petsc-2.1.5/lib/libg_complex/linux/libpetscsles.so, may conflict with libstdc++.so.6
***************************************************************
* Running Example  ./ex14
***************************************************************
 
 EquationSystems
  n_systems()=1
   System "Laplace"
    Type "Implicit"
    Variables="u" 
    Finite Element Types="0", "12" 
    Infinite Element Mapping="0" 
    Approximation Orders="2", "3" 
    n_dofs()=21
    n_local_dofs()=21
    n_constrained_dofs()=0
    n_vectors()=1
  n_parameters()=2
   Parameters:
    "linear solver maximum iterations"=100
    "linear solver tolerance"=1e-12

Beginning Solve 0
System has: 21 degrees of freedom.
Linear solver converged at step: 5, final residual: 3.10189e-16
L2-Error is: (0.0151553,0)
H1-Error is: (0.125274,0)
  Refining the mesh...
Beginning Solve 1
System has: 65 degrees of freedom.
Linear solver converged at step: 13, final residual: 3.32844e-14
L2-Error is: (0.00515182,0)
H1-Error is: (0.077772,0)
  Refining the mesh...
Beginning Solve 2
System has: 97 degrees of freedom.
Linear solver converged at step: 19, final residual: 4.249e-15
L2-Error is: (0.00198908,0)
H1-Error is: (0.049429,0)
  Refining the mesh...
Beginning Solve 3
System has: 129 degrees of freedom.
Linear solver converged at step: 22, final residual: 2.36346e-14
L2-Error is: (0.000892486,0)
H1-Error is: (0.0320038,0)
  Refining the mesh...
Beginning Solve 4
System has: 161 degrees of freedom.
Linear solver converged at step: 25, final residual: 1.227e-14
L2-Error is: (0.000566126,0)
H1-Error is: (0.0215045,0)
  Refining the mesh...
Beginning Solve 5
System has: 193 degrees of freedom.
Linear solver converged at step: 28, final residual: 6.04725e-15
L2-Error is: (0.000492378,0)
H1-Error is: (0.0154804,0)
  Refining the mesh...
Beginning Solve 6
System has: 225 degrees of freedom.
Linear solver converged at step: 30, final residual: 4.60327e-15
L2-Error is: (0.000477307,0)
H1-Error is: (0.0122975,0)
  Refining the mesh...
Beginning Solve 7
System has: 257 degrees of freedom.
Linear solver converged at step: 32, final residual: 3.9266e-15
L2-Error is: (0.000473809,0)
H1-Error is: (0.0107766,0)
  Refining the mesh...
Beginning Solve 8
System has: 313 degrees of freedom.
Linear solver converged at step: 42, final residual: 3.45761e-15
L2-Error is: (0.000391032,0)
H1-Error is: (0.00901463,0)
  Refining the mesh...
Beginning Solve 9
System has: 417 degrees of freedom.
Linear solver converged at step: 50, final residual: 3.90649e-15
L2-Error is: (0.000288213,0)
H1-Error is: (0.00692399,0)
  Refining the mesh...
Beginning Solve 10
System has: 601 degrees of freedom.
Linear solver converged at step: 55, final residual: 2.97739e-15
L2-Error is: (0.000103614,0)
H1-Error is: (0.00417121,0)
  Refining the mesh...
Beginning Solve 11
System has: 701 degrees of freedom.
Linear solver converged at step: 63, final residual: 1.05751e-15
L2-Error is: (0.000100765,0)
H1-Error is: (0.00353862,0)
  Refining the mesh...
Beginning Solve 12
System has: 801 degrees of freedom.
Linear solver converged at step: 72, final residual: 7.28591e-16
L2-Error is: (0.000100503,0)
H1-Error is: (0.00325363,0)
  Refining the mesh...
Beginning Solve 13
System has: 1089 degrees of freedom.
Linear solver converged at step: 80, final residual: 1.34075e-15
L2-Error is: (4.26027e-05,0)
H1-Error is: (0.00226665,0)
  Refining the mesh...
Beginning Solve 14
System has: 1365 degrees of freedom.
Linear solver converged at step: 95, final residual: 8.71036e-16
L2-Error is: (3.47025e-05,0)
H1-Error is: (0.00178852,0)

 ---------------------------------------------------------------------------- 
| Reference count information                                                |
 ---------------------------------------------------------------------------- 
| 12SparseMatrixISt7complexIdEE reference count information:
|  Creations:    1
|  Destructions: 1
| 13NumericVectorISt7complexIdEE reference count information:
|  Creations:    17
|  Destructions: 17
| 21LinearSolverInterfaceISt7complexIdEE reference count information:
|  Creations:    1
|  Destructions: 1
| 4Elem reference count information:
|  Creations:    15331
|  Destructions: 15331
| 4Node reference count information:
|  Creations:    1609
|  Destructions: 1609
| 5QBase reference count information:
|  Creations:    89
|  Destructions: 89
| 6DofMap reference count information:
|  Creations:    1
|  Destructions: 1
| 6FEBase reference count information:
|  Creations:    73
|  Destructions: 73
| 6System reference count information:
|  Creations:    1
|  Destructions: 1
 ---------------------------------------------------------------------------- 
 
***************************************************************
* Done Running Example  ./ex14
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
