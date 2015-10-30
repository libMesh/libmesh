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
        #include "linear_implicit_system.h"
        #include "gmv_io.h"
        #include "tecplot_io.h"
        #include "fe.h"
        #include "quadrature_gauss.h"
        #include "dense_matrix.h"
        #include "dense_vector.h"
        #include "sparse_matrix.h"
        #include "mesh_refinement.h"
        #include "error_vector.h"
        #include "exact_error_estimator.h"
        #include "kelly_error_estimator.h"
        #include "patch_recovery_error_estimator.h"
        #include "uniform_refinement_estimator.h"
        #include "hp_coarsentest.h"
        #include "hp_singular.h"
        #include "mesh_generation.h"
        #include "mesh_modification.h"
        #include "getpot.h"
        #include "exact_solution.h"
        #include "dof_map.h"
        #include "numeric_vector.h"
        #include "elem.h"
        #include "string_to_enum.h"
        
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
        		      const Parameters&,   // EquationSystem parameters, not needed
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
        			  const Parameters&,   // EquationSystems parameters, not needed
        			  const std::string&,  // sys_name, not needed
        			  const std::string&); // unk_name, not needed);
        
        
</pre>
</div>
<div class = "comment">
These are non-const because the input file may change it,
It is global because our exact_* functions use it.


<br><br>Set the dimensionality of the mesh
</div>

<div class ="fragment">
<pre>
        unsigned int dim = 2;
        
</pre>
</div>
<div class = "comment">
Choose whether or not to use the singular solution
</div>

<div class ="fragment">
<pre>
        bool singularity = true;
        
        
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
            const unsigned int max_r_steps    = input_file("max_r_steps", 3);
            const unsigned int max_r_level    = input_file("max_r_level", 3);
            const Real refine_percentage      = input_file("refine_percentage", 0.5);
            const Real coarsen_percentage     = input_file("coarsen_percentage", 0.5);
            const unsigned int uniform_refine = input_file("uniform_refine",0);
            const std::string refine_type     = input_file("refinement_type", "h");
            const std::string approx_type     = input_file("approx_type", "LAGRANGE");
            const unsigned int approx_order   = input_file("approx_order", 1);
            const std::string element_type    = input_file("element_type", "tensor");
            const int extra_error_quadrature  = input_file("extra_error_quadrature", 0);
            const int max_linear_iterations   = input_file("max_linear_iterations", 5000);
            dim = input_file("dimension", 2);
            const std::string indicator_type = input_file("indicator_type", "kelly");
            singularity = input_file("singularity", true);
            
</pre>
</div>
<div class = "comment">
Output file for plotting the error as a function of
the number of degrees of freedom.
</div>

<div class ="fragment">
<pre>
            std::string approx_name = "";
            if (element_type == "tensor")
              approx_name += "bi";
            if (approx_order == 1)
              approx_name += "linear";
            else if (approx_order == 2)
              approx_name += "quadratic";
            else if (approx_order == 3)
              approx_name += "cubic";
            else if (approx_order == 4)
              approx_name += "quartic";
        
            std::string output_file = approx_name;
            output_file += "_";
            output_file += refine_type;
            if (uniform_refine == 0)
              output_file += "_adaptive.m";
            else
              output_file += "_uniform.m";
            
            std::ofstream out (output_file.c_str());
            out &lt;&lt; "% dofs     L2-error     H1-error" &lt;&lt; std::endl;
            out &lt;&lt; "e = [" &lt;&lt; std::endl;
            
</pre>
</div>
<div class = "comment">
Create an n-dimensional mesh.
</div>

<div class ="fragment">
<pre>
            Mesh mesh (dim);
            
</pre>
</div>
<div class = "comment">
Read in the mesh
</div>

<div class ="fragment">
<pre>
            if (dim == 1)
              MeshTools::Generation::build_line(mesh,1,-1.,0.);
            else if (dim == 2)
              mesh.read("lshaped.xda");
            else
              mesh.read("lshaped3D.xda");
        
</pre>
</div>
<div class = "comment">
Use triangles if the config file says so
</div>

<div class ="fragment">
<pre>
            if (element_type == "simplex")
              MeshTools::Modification::all_tri(mesh);
        
</pre>
</div>
<div class = "comment">
We used first order elements to describe the geometry,
but we may need second order elements to hold the degrees
of freedom
</div>

<div class ="fragment">
<pre>
            if (approx_order &gt; 1 || refine_type != "h")
              mesh.all_second_order();
        
</pre>
</div>
<div class = "comment">
Mesh Refinement object
</div>

<div class ="fragment">
<pre>
            MeshRefinement mesh_refinement(mesh);
            mesh_refinement.refine_fraction() = refine_percentage;
            mesh_refinement.coarsen_fraction() = coarsen_percentage;
            mesh_refinement.max_h_level() = max_r_level;
        
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
              LinearImplicitSystem& system =
        	equation_systems.add_system&lt;LinearImplicitSystem&gt; ("Laplace");
              
</pre>
</div>
<div class = "comment">
Adds the variable "u" to "Laplace", using 
the finite element type and order specified
in the config file
</div>

<div class ="fragment">
<pre>
              system.add_variable("u", static_cast&lt;Order&gt;(approx_order),
                                  Utility::string_to_enum&lt;FEFamily&gt;(approx_type));
        
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly
function.
</div>

<div class ="fragment">
<pre>
              system.attach_assemble_function (assemble_laplace);
        
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
              equation_systems.parameters.set&lt;unsigned int&gt;("linear solver maximum iterations")
                = max_linear_iterations;
        
</pre>
</div>
<div class = "comment">
Linear solver tolerance.
</div>

<div class ="fragment">
<pre>
              equation_systems.parameters.set&lt;Real&gt;("linear solver tolerance") =
                TOLERANCE * TOLERANCE * TOLERANCE;
              
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
Construct ExactSolution object and attach solution functions
</div>

<div class ="fragment">
<pre>
            ExactSolution exact_sol(equation_systems);
            exact_sol.attach_exact_value(exact_solution);
            exact_sol.attach_exact_deriv(exact_derivative);
        
</pre>
</div>
<div class = "comment">
Use higher quadrature order for more accurate error results
</div>

<div class ="fragment">
<pre>
            exact_sol.extra_quadrature_order(extra_error_quadrature);
        
</pre>
</div>
<div class = "comment">
Convenient reference to the system
</div>

<div class ="fragment">
<pre>
            LinearImplicitSystem& system =
              equation_systems.get_system&lt;LinearImplicitSystem&gt;("Laplace");
        
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
                system.solve();
        
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
        		
                        if (indicator_type == "exact")
                          {
</pre>
</div>
<div class = "comment">
The \p ErrorEstimator class interrogates a
finite element solution and assigns to each
element a positive error value.
This value is used for deciding which elements to
refine and which to coarsen.
For these simple test problems, we can use
numerical quadrature of the exact error between
the approximate and analytic solutions.
However, for real problems, we would need an error
indicator which only relies on the approximate
solution.
</div>

<div class ="fragment">
<pre>
                            ExactErrorEstimator error_estimator;
        
                            error_estimator.attach_exact_value(exact_solution);
                            error_estimator.attach_exact_deriv(exact_derivative);
        
</pre>
</div>
<div class = "comment">
We optimize in H1 norm
</div>

<div class ="fragment">
<pre>
                            error_estimator.sobolev_order() = 1;
        
</pre>
</div>
<div class = "comment">
Compute the error for each active element using
the provided indicator.  Note in general you
will need to provide an error estimator
specifically designed for your application.
</div>

<div class ="fragment">
<pre>
                            error_estimator.estimate_error (system, error);
                          }
                        else if (indicator_type == "patch")
                          {
</pre>
</div>
<div class = "comment">
The patch recovery estimator should give a
good estimate of the solution interpolation
error.
</div>

<div class ="fragment">
<pre>
                            PatchRecoveryErrorEstimator error_estimator;
        
        		    error_estimator.estimate_error (system, error);
                          }
                        else if (indicator_type == "uniform")
                          {
</pre>
</div>
<div class = "comment">
Error indication based on uniform refinement
is reliable, but very expensive.
</div>

<div class ="fragment">
<pre>
                            UniformRefinementEstimator error_estimator;
        
        		    error_estimator.estimate_error (system, error);
                          }
                        else
                          {
                            assert (indicator_type == "kelly");
        
</pre>
</div>
<div class = "comment">
The Kelly error estimator is based on 
an error bound for the Poisson problem
on linear elements, but is useful for
driving adaptive refinement in many problems
</div>

<div class ="fragment">
<pre>
                            KellyErrorEstimator error_estimator;
        
        		    error_estimator.estimate_error (system, error);
                          }
        		
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
                        mesh_refinement.flag_elements_by_error_fraction (error);
        
</pre>
</div>
<div class = "comment">
If we are doing adaptive p refinement, we want
elements flagged for that instead.
</div>

<div class ="fragment">
<pre>
                        if (refine_type == "p")
                          mesh_refinement.switch_h_to_p_refinement();
</pre>
</div>
<div class = "comment">
If we are doing "matched hp" refinement, we
flag elements for both h and p
</div>

<div class ="fragment">
<pre>
                        if (refine_type == "matchedhp")
                          mesh_refinement.add_p_to_h_refinement();
</pre>
</div>
<div class = "comment">
If we are doing hp refinement, we 
try switching some elements from h to p
</div>

<div class ="fragment">
<pre>
                        if (refine_type == "hp")
        	          {
        		    HPCoarsenTest hpselector;
                            hpselector.select_refinement(system);
        	          }
</pre>
</div>
<div class = "comment">
If we are doing "singular hp" refinement, we 
try switching most elements from h to p
</div>

<div class ="fragment">
<pre>
                        if (refine_type == "singularhp")
        	          {
</pre>
</div>
<div class = "comment">
This only differs from p refinement for
the singular problem
</div>

<div class ="fragment">
<pre>
                            assert (singularity);
        		    HPSingularity hpselector;
</pre>
</div>
<div class = "comment">
Our only singular point is at the origin
</div>

<div class ="fragment">
<pre>
                            hpselector.singular_points.push_back(Point());
                            hpselector.select_refinement(system);
        	          }
        		
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
                      {
                        if (refine_type == "h" || refine_type == "hp" ||
                            refine_type == "matchedhp")
                          mesh_refinement.uniformly_refine(1);
                        if (refine_type == "p" || refine_type == "hp" ||
                            refine_type == "matchedhp")
                          mesh_refinement.uniformly_p_refine(1);
                      }
        	    
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
            out &lt;&lt; "title('" &lt;&lt; approx_name &lt;&lt; " elements');" &lt;&lt; std::endl;
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
We now define the exact solution, being careful
to obtain an angle from atan2 in the correct
quadrant.
</div>

<div class ="fragment">
<pre>
        Number exact_solution(const Point& p,
        		      const Parameters&,  // parameters, not needed
        		      const std::string&, // sys_name, not needed
        		      const std::string&) // unk_name, not needed
        {
          const Real x = p(0);
          const Real y = (dim &gt; 1) ? p(1) : 0.;
          
          if (singularity)
            {
</pre>
</div>
<div class = "comment">
The exact solution to the singular problem,
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
Make the 3D solution similar
</div>

<div class ="fragment">
<pre>
              const Real z = (dim &gt; 2) ? p(2) : 0;
        		  
              return pow(x*x + y*y, 1./3.)*sin(2./3.*theta) + z;
            }
          else
            {
</pre>
</div>
<div class = "comment">
The exact solution to a nonsingular problem,
good for testing ideal convergence rates
</div>

<div class ="fragment">
<pre>
              const Real z = (dim &gt; 2) ? p(2) : 0;
        
              return cos(x) * exp(y) * (1. - z);
            }
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
        			  const Parameters&,  // parameters, not needed
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
          const Real y = dim &gt; 1 ? p(1) : 0.;
        
          if (singularity)
            {
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
              if (dim &gt; 1)
                gradu(1) = tt*y*pow(r2,-tt)*sin(tt*theta) + pow(r2,ot)*cos(tt*theta)*tt/(1.+y*y/x/x)*1./x;
        
              if (dim &gt; 2)
                gradu(2) = 1.;
            }
          else
            {
              const Real z = (dim &gt; 2) ? p(2) : 0;
        
              gradu(0) = -sin(x) * exp(y) * (1. - z);
              if (dim &gt; 1)
                gradu(1) = cos(x) * exp(y) * (1. - z);
              if (dim &gt; 2)
                gradu(2) = -cos(x) * exp(y);
            }
        
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
        #ifdef ENABLE_AMR
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
Get a reference to the LinearImplicitSystem we are solving
</div>

<div class ="fragment">
<pre>
          LinearImplicitSystem& system = es.get_system&lt;LinearImplicitSystem&gt;("Laplace");
          
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
          AutoPtr&lt;FEBase&gt; fe      (FEBase::build(dim, fe_type));
          AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
          
</pre>
</div>
<div class = "comment">
Quadrature rules for numerical integration.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;QBase&gt; qrule(fe_type.default_quadrature_rule(dim));
          AutoPtr&lt;QBase&gt; qface(fe_type.default_quadrature_rule(dim-1));
        
</pre>
</div>
<div class = "comment">
Tell the finite element object to use our quadrature rule.
</div>

<div class ="fragment">
<pre>
          fe-&gt;attach_quadrature_rule      (qrule.get());
          fe_face-&gt;attach_quadrature_rule (qface.get());
        
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
          const std::vector&lt;Real&gt;& JxW      = fe-&gt;get_JxW();
          const std::vector&lt;Real&gt;& JxW_face = fe_face-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The physical XY locations of the quadrature points on the element.
These might be useful for evaluating spatially varying material
properties or forcing functions at the quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Point&gt;& q_point = fe-&gt;get_xyz();
        
</pre>
</div>
<div class = "comment">
The element shape functions evaluated at the quadrature points.
For this simple problem we usually only need them on element
boundaries.
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
        
              for (unsigned int qp=0; qp&lt;qrule-&gt;n_points(); qp++)
        	for (unsigned int i=0; i&lt;dphi.size(); i++)
        	  for (unsigned int j=0; j&lt;dphi.size(); j++)
        	    Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
        
</pre>
</div>
<div class = "comment">
We need a forcing function to make the 1D case interesting
</div>

<div class ="fragment">
<pre>
              if (dim == 1)
                for (unsigned int qp=0; qp&lt;qrule-&gt;n_points(); qp++)
                  {
                    Real x = q_point[qp](0);
                    Real f = singularity ? sqrt(3.)/9.*pow(-x, -4./3.) :
                                           cos(x);
                    for (unsigned int i=0; i&lt;dphi.size(); ++i)
                      Fe(i) += JxW[qp]*phi[i][qp]*f;
                  }
        
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
via the penalty method.

<br><br>This approach adds the L2 projection of the boundary
data in penalty form to the weak statement.  This is
a more generic approach for applying Dirichlet BCs
which is applicable to non-Lagrange finite element
discretizations.
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
        	      fe_face-&gt;reinit(elem,s);
        	      
        	      for (unsigned int qp=0; qp&lt;qface-&gt;n_points(); qp++)
        		{
        		  const Number value = exact_solution (qface_points[qp],
        						       es.parameters,
        						       "null",
        						       "void");
        
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
        #endif // #ifdef ENABLE_AMR
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;gmv_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;tecplot_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exact_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;kelly_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;patch_recovery_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;uniform_refinement_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;hp_coarsentest.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;hp_singular.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_modification.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exact_solution.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;string_to_enum.h&quot;</FONT></B>
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_laplace(EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  
  Number exact_solution(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  		      <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,   <I><FONT COLOR="#B22222">// EquationSystem parameters, not needed
</FONT></I>  		      <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,  <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>  		      <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;); <I><FONT COLOR="#B22222">// unk_name, not needed);
</FONT></I>  
  Gradient exact_derivative(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  			  <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,   <I><FONT COLOR="#B22222">// EquationSystems parameters, not needed
</FONT></I>  			  <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,  <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>  			  <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;); <I><FONT COLOR="#B22222">// unk_name, not needed);
</FONT></I>  
  
  
  <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = 2;
  
  <B><FONT COLOR="#228B22">bool</FONT></B> singularity = true;
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main(<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::init (argc, argv);
  
  #ifndef ENABLE_AMR
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR: This example requires libMesh to be\n&quot;</FONT></B>
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;compiled with AMR support!&quot;</FONT></B>
              &lt;&lt; std::endl;
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    {
      GetPot input_file(<B><FONT COLOR="#BC8F8F">&quot;ex14.in&quot;</FONT></B>);
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_r_steps    = input_file(<B><FONT COLOR="#BC8F8F">&quot;max_r_steps&quot;</FONT></B>, 3);
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_r_level    = input_file(<B><FONT COLOR="#BC8F8F">&quot;max_r_level&quot;</FONT></B>, 3);
      <B><FONT COLOR="#228B22">const</FONT></B> Real refine_percentage      = input_file(<B><FONT COLOR="#BC8F8F">&quot;refine_percentage&quot;</FONT></B>, 0.5);
      <B><FONT COLOR="#228B22">const</FONT></B> Real coarsen_percentage     = input_file(<B><FONT COLOR="#BC8F8F">&quot;coarsen_percentage&quot;</FONT></B>, 0.5);
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> uniform_refine = input_file(<B><FONT COLOR="#BC8F8F">&quot;uniform_refine&quot;</FONT></B>,0);
      <B><FONT COLOR="#228B22">const</FONT></B> std::string refine_type     = input_file(<B><FONT COLOR="#BC8F8F">&quot;refinement_type&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;h&quot;</FONT></B>);
      <B><FONT COLOR="#228B22">const</FONT></B> std::string approx_type     = input_file(<B><FONT COLOR="#BC8F8F">&quot;approx_type&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;LAGRANGE&quot;</FONT></B>);
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> approx_order   = input_file(<B><FONT COLOR="#BC8F8F">&quot;approx_order&quot;</FONT></B>, 1);
      <B><FONT COLOR="#228B22">const</FONT></B> std::string element_type    = input_file(<B><FONT COLOR="#BC8F8F">&quot;element_type&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;tensor&quot;</FONT></B>);
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> extra_error_quadrature  = input_file(<B><FONT COLOR="#BC8F8F">&quot;extra_error_quadrature&quot;</FONT></B>, 0);
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_linear_iterations   = input_file(<B><FONT COLOR="#BC8F8F">&quot;max_linear_iterations&quot;</FONT></B>, 5000);
      dim = input_file(<B><FONT COLOR="#BC8F8F">&quot;dimension&quot;</FONT></B>, 2);
      <B><FONT COLOR="#228B22">const</FONT></B> std::string indicator_type = input_file(<B><FONT COLOR="#BC8F8F">&quot;indicator_type&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;kelly&quot;</FONT></B>);
      singularity = input_file(<B><FONT COLOR="#BC8F8F">&quot;singularity&quot;</FONT></B>, true);
      
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string approx_name = <B><FONT COLOR="#BC8F8F">&quot;&quot;</FONT></B>;
      <B><FONT COLOR="#A020F0">if</FONT></B> (element_type == <B><FONT COLOR="#BC8F8F">&quot;tensor&quot;</FONT></B>)
        approx_name += <B><FONT COLOR="#BC8F8F">&quot;bi&quot;</FONT></B>;
      <B><FONT COLOR="#A020F0">if</FONT></B> (approx_order == 1)
        approx_name += <B><FONT COLOR="#BC8F8F">&quot;linear&quot;</FONT></B>;
      <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (approx_order == 2)
        approx_name += <B><FONT COLOR="#BC8F8F">&quot;quadratic&quot;</FONT></B>;
      <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (approx_order == 3)
        approx_name += <B><FONT COLOR="#BC8F8F">&quot;cubic&quot;</FONT></B>;
      <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (approx_order == 4)
        approx_name += <B><FONT COLOR="#BC8F8F">&quot;quartic&quot;</FONT></B>;
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string output_file = approx_name;
      output_file += <B><FONT COLOR="#BC8F8F">&quot;_&quot;</FONT></B>;
      output_file += refine_type;
      <B><FONT COLOR="#A020F0">if</FONT></B> (uniform_refine == 0)
        output_file += <B><FONT COLOR="#BC8F8F">&quot;_adaptive.m&quot;</FONT></B>;
      <B><FONT COLOR="#A020F0">else</FONT></B>
        output_file += <B><FONT COLOR="#BC8F8F">&quot;_uniform.m&quot;</FONT></B>;
      
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream out (output_file.c_str());
      out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;% dofs     L2-error     H1-error&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;e = [&quot;</FONT></B> &lt;&lt; std::endl;
      
      Mesh mesh (dim);
      
      <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 1)
        <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_line(mesh,1,-1.,0.);
      <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 2)
        mesh.read(<B><FONT COLOR="#BC8F8F">&quot;lshaped.xda&quot;</FONT></B>);
      <B><FONT COLOR="#A020F0">else</FONT></B>
        mesh.read(<B><FONT COLOR="#BC8F8F">&quot;lshaped3D.xda&quot;</FONT></B>);
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (element_type == <B><FONT COLOR="#BC8F8F">&quot;simplex&quot;</FONT></B>)
        <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Modification::all_tri(mesh);
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (approx_order &gt; 1 || refine_type != <B><FONT COLOR="#BC8F8F">&quot;h&quot;</FONT></B>)
        mesh.all_second_order();
  
      MeshRefinement mesh_refinement(mesh);
      mesh_refinement.refine_fraction() = refine_percentage;
      mesh_refinement.coarsen_fraction() = coarsen_percentage;
      mesh_refinement.max_h_level() = max_r_level;
  
      EquationSystems equation_systems (mesh);
  
      {
        LinearImplicitSystem&amp; system =
  	equation_systems.add_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>);
        
        system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, static_cast&lt;Order&gt;(approx_order),
                            <B><FONT COLOR="#5F9EA0">Utility</FONT></B>::string_to_enum&lt;FEFamily&gt;(approx_type));
  
        system.attach_assemble_function (assemble_laplace);
  
        equation_systems.init();
  
        equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;linear solver maximum iterations&quot;</FONT></B>)
          = max_linear_iterations;
  
        equation_systems.parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;linear solver tolerance&quot;</FONT></B>) =
          TOLERANCE * TOLERANCE * TOLERANCE;
        
        equation_systems.print_info();
      }
  
      ExactSolution exact_sol(equation_systems);
      exact_sol.attach_exact_value(exact_solution);
      exact_sol.attach_exact_deriv(exact_derivative);
  
      exact_sol.extra_quadrature_order(extra_error_quadrature);
  
      LinearImplicitSystem&amp; system =
        equation_systems.get_system&lt;LinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>);
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> r_step=0; r_step&lt;max_r_steps; r_step++)
        {
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Beginning Solve &quot;</FONT></B> &lt;&lt; r_step &lt;&lt; std::endl;
  	
  	system.solve();
  
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;System has: &quot;</FONT></B> &lt;&lt; equation_systems.n_active_dofs()
  		  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; degrees of freedom.&quot;</FONT></B>
  		  &lt;&lt; std::endl;
  
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Linear solver converged at step: &quot;</FONT></B>
  		  &lt;&lt; system.n_linear_iterations()
  		  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, final residual: &quot;</FONT></B>
  		  &lt;&lt; system.final_linear_residual()
  		  &lt;&lt; std::endl;
  	
  	exact_sol.compute_error(<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>);
  
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;L2-Error is: &quot;</FONT></B>
  		  &lt;&lt; exact_sol.l2_error(<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)
  		  &lt;&lt; std::endl;
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;H1-Error is: &quot;</FONT></B>
  		  &lt;&lt; exact_sol.h1_error(<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)
  		  &lt;&lt; std::endl;
  
  	out &lt;&lt; equation_systems.n_active_dofs() &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B>
  	    &lt;&lt; exact_sol.l2_error(<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>) &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B>
  	    &lt;&lt; exact_sol.h1_error(<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>) &lt;&lt; std::endl;
  
  	<B><FONT COLOR="#A020F0">if</FONT></B> (r_step+1 != max_r_steps)
  	  {
  	    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;  Refining the mesh...&quot;</FONT></B> &lt;&lt; std::endl;
  
  	    <B><FONT COLOR="#A020F0">if</FONT></B> (uniform_refine == 0)
                {
  
  		ErrorVector error;
  		
                  <B><FONT COLOR="#A020F0">if</FONT></B> (indicator_type == <B><FONT COLOR="#BC8F8F">&quot;exact&quot;</FONT></B>)
                    {
                      ExactErrorEstimator error_estimator;
  
                      error_estimator.attach_exact_value(exact_solution);
                      error_estimator.attach_exact_deriv(exact_derivative);
  
                      error_estimator.sobolev_order() = 1;
  
  		    error_estimator.estimate_error (system, error);
                    }
                  <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (indicator_type == <B><FONT COLOR="#BC8F8F">&quot;patch&quot;</FONT></B>)
                    {
  		    PatchRecoveryErrorEstimator error_estimator;
  
  		    error_estimator.estimate_error (system, error);
                    }
                  <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (indicator_type == <B><FONT COLOR="#BC8F8F">&quot;uniform&quot;</FONT></B>)
                    {
                      UniformRefinementEstimator error_estimator;
  
  		    error_estimator.estimate_error (system, error);
                    }
                  <B><FONT COLOR="#A020F0">else</FONT></B>
                    {
                      assert (indicator_type == <B><FONT COLOR="#BC8F8F">&quot;kelly&quot;</FONT></B>);
  
  		    KellyErrorEstimator error_estimator;
  
  		    error_estimator.estimate_error (system, error);
                    }
  		
  		mesh_refinement.flag_elements_by_error_fraction (error);
  
                  <B><FONT COLOR="#A020F0">if</FONT></B> (refine_type == <B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>)
                    mesh_refinement.switch_h_to_p_refinement();
                  <B><FONT COLOR="#A020F0">if</FONT></B> (refine_type == <B><FONT COLOR="#BC8F8F">&quot;matchedhp&quot;</FONT></B>)
                    mesh_refinement.add_p_to_h_refinement();
                  <B><FONT COLOR="#A020F0">if</FONT></B> (refine_type == <B><FONT COLOR="#BC8F8F">&quot;hp&quot;</FONT></B>)
  	          {
  		    HPCoarsenTest hpselector;
                      hpselector.select_refinement(system);
  	          }
                  <B><FONT COLOR="#A020F0">if</FONT></B> (refine_type == <B><FONT COLOR="#BC8F8F">&quot;singularhp&quot;</FONT></B>)
  	          {
                      assert (singularity);
  		    HPSingularity hpselector;
                      hpselector.singular_points.push_back(Point());
                      hpselector.select_refinement(system);
  	          }
  		
  		mesh_refinement.refine_and_coarsen_elements();
  	      }
  
  	    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (uniform_refine == 1)
                {
                  <B><FONT COLOR="#A020F0">if</FONT></B> (refine_type == <B><FONT COLOR="#BC8F8F">&quot;h&quot;</FONT></B> || refine_type == <B><FONT COLOR="#BC8F8F">&quot;hp&quot;</FONT></B> ||
                      refine_type == <B><FONT COLOR="#BC8F8F">&quot;matchedhp&quot;</FONT></B>)
                    mesh_refinement.uniformly_refine(1);
                  <B><FONT COLOR="#A020F0">if</FONT></B> (refine_type == <B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B> || refine_type == <B><FONT COLOR="#BC8F8F">&quot;hp&quot;</FONT></B> ||
                      refine_type == <B><FONT COLOR="#BC8F8F">&quot;matchedhp&quot;</FONT></B>)
                    mesh_refinement.uniformly_p_refine(1);
                }
  	    
  	    equation_systems.reinit ();
  	  }
        }	    
      
      
  
      
      GMVIO (mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;lshaped.gmv&quot;</FONT></B>,
      					 equation_systems);
  
      out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;];&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;hold on&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;plot(e(:,1), e(:,2), 'bo-');&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;plot(e(:,1), e(:,3), 'ro-');&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;xlabel('dofs');&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;title('&quot;</FONT></B> &lt;&lt; approx_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; elements');&quot;</FONT></B> &lt;&lt; std::endl;
      out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;legend('L2-error', 'H1-error');&quot;</FONT></B> &lt;&lt; std::endl;
    }
  #endif <I><FONT COLOR="#B22222">// #ifndef ENABLE_AMR
</FONT></I>    
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close ();
  }
  
  
  
  
  Number exact_solution(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  		      <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,  <I><FONT COLOR="#B22222">// parameters, not needed
</FONT></I>  		      <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>  		      <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;) <I><FONT COLOR="#B22222">// unk_name, not needed
</FONT></I>  {
    <B><FONT COLOR="#228B22">const</FONT></B> Real x = p(0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real y = (dim &gt; 1) ? p(1) : 0.;
    
    <B><FONT COLOR="#A020F0">if</FONT></B> (singularity)
      {
        Real theta = atan2(y,x);
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (theta &lt; 0)
          theta += 2. * libMesh::pi;
  
        <B><FONT COLOR="#228B22">const</FONT></B> Real z = (dim &gt; 2) ? p(2) : 0;
  		  
        <B><FONT COLOR="#A020F0">return</FONT></B> pow(x*x + y*y, 1./3.)*sin(2./3.*theta) + z;
      }
    <B><FONT COLOR="#A020F0">else</FONT></B>
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Real z = (dim &gt; 2) ? p(2) : 0;
  
        <B><FONT COLOR="#A020F0">return</FONT></B> cos(x) * exp(y) * (1. - z);
      }
  }
  
  
  
  
  
  Gradient exact_derivative(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
  			  <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,  <I><FONT COLOR="#B22222">// parameters, not needed
</FONT></I>  			  <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>  			  <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;) <I><FONT COLOR="#B22222">// unk_name, not needed
</FONT></I>  {
    Gradient gradu;
    
    <B><FONT COLOR="#228B22">const</FONT></B> Real x = p(0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real y = dim &gt; 1 ? p(1) : 0.;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (singularity)
      {
        assert (x != 0.);
  
        <B><FONT COLOR="#228B22">const</FONT></B> Real tt = 2./3.;
        <B><FONT COLOR="#228B22">const</FONT></B> Real ot = 1./3.;
    
        <B><FONT COLOR="#228B22">const</FONT></B> Real r2 = x*x + y*y;
  
        Real theta = atan2(y,x);
    
        <B><FONT COLOR="#A020F0">if</FONT></B> (theta &lt; 0)
          theta += 2. * libMesh::pi;
  
        gradu(0) = tt*x*pow(r2,-tt)*sin(tt*theta) - pow(r2,ot)*cos(tt*theta)*tt/(1.+y*y/x/x)*y/x/x;
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (dim &gt; 1)
          gradu(1) = tt*y*pow(r2,-tt)*sin(tt*theta) + pow(r2,ot)*cos(tt*theta)*tt/(1.+y*y/x/x)*1./x;
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (dim &gt; 2)
          gradu(2) = 1.;
      }
    <B><FONT COLOR="#A020F0">else</FONT></B>
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Real z = (dim &gt; 2) ? p(2) : 0;
  
        gradu(0) = -sin(x) * exp(y) * (1. - z);
        <B><FONT COLOR="#A020F0">if</FONT></B> (dim &gt; 1)
          gradu(1) = cos(x) * exp(y) * (1. - z);
        <B><FONT COLOR="#A020F0">if</FONT></B> (dim &gt; 2)
          gradu(2) = -cos(x) * exp(y);
      }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> gradu;
  }
  
  
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_laplace(EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
  #ifdef ENABLE_AMR
    assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>);
  
  
    PerfLog perf_log (<B><FONT COLOR="#BC8F8F">&quot;Matrix Assembly&quot;</FONT></B>,false);
    
    <B><FONT COLOR="#228B22">const</FONT></B> Mesh&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    LinearImplicitSystem&amp; system = es.get_system&lt;LinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
  
    FEType fe_type = dof_map.variable_type(0);
  
    AutoPtr&lt;FEBase&gt; fe      (FEBase::build(dim, fe_type));
    AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
    
    AutoPtr&lt;QBase&gt; qrule(fe_type.default_quadrature_rule(dim));
    AutoPtr&lt;QBase&gt; qface(fe_type.default_quadrature_rule(dim-1));
  
    fe-&gt;attach_quadrature_rule      (qrule.get());
    fe_face-&gt;attach_quadrature_rule (qface.get());
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW      = fe-&gt;get_JxW();
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; q_point = fe-&gt;get_xyz();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; psi = fe_face-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; qface_points = fe_face-&gt;get_xyz();
  
    DenseMatrix&lt;Number&gt; Ke;
    DenseVector&lt;Number&gt; Fe;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
    
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        perf_log.start_event(<B><FONT COLOR="#BC8F8F">&quot;elem init&quot;</FONT></B>);      
  
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
  
        fe-&gt;reinit (elem);
  
        Ke.resize (dof_indices.size(),
  		 dof_indices.size());
  
        Fe.resize (dof_indices.size());
  
        perf_log.stop_event(<B><FONT COLOR="#BC8F8F">&quot;elem init&quot;</FONT></B>);      
  
        perf_log.start_event (<B><FONT COLOR="#BC8F8F">&quot;Ke&quot;</FONT></B>);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule-&gt;n_points(); qp++)
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;dphi.size(); i++)
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;dphi.size(); j++)
  	    Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 1)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule-&gt;n_points(); qp++)
            {
              Real x = q_point[qp](0);
              Real f = singularity ? sqrt(3.)/9.*pow(-x, -4./3.) :
                                     cos(x);
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;dphi.size(); ++i)
                Fe(i) += JxW[qp]*phi[i][qp]*f;
            }
  
        perf_log.stop_event (<B><FONT COLOR="#BC8F8F">&quot;Ke&quot;</FONT></B>);
  
  
        {
  	perf_log.start_event (<B><FONT COLOR="#BC8F8F">&quot;BCs&quot;</FONT></B>);
  
  	<B><FONT COLOR="#228B22">const</FONT></B> Real penalty = 1.e10;
  
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> s=0; s&lt;elem-&gt;n_sides(); s++)
  	  <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(s) == NULL)
  	    {
  	      fe_face-&gt;reinit(elem,s);
  	      
  	      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qface-&gt;n_points(); qp++)
  		{
  		  <B><FONT COLOR="#228B22">const</FONT></B> Number value = exact_solution (qface_points[qp],
  						       es.parameters,
  						       <B><FONT COLOR="#BC8F8F">&quot;null&quot;</FONT></B>,
  						       <B><FONT COLOR="#BC8F8F">&quot;void&quot;</FONT></B>);
  
  		  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;psi.size(); i++)
  		    Fe(i) += penalty*JxW_face[qp]*value*psi[i][qp];
  
  		  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;psi.size(); i++)
  		    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;psi.size(); j++)
  		      Ke(i,j) += penalty*JxW_face[qp]*psi[i][qp]*psi[j][qp];
  		}
  	    } 
  	
  	perf_log.stop_event (<B><FONT COLOR="#BC8F8F">&quot;BCs&quot;</FONT></B>);
        } 
        
  
        perf_log.start_event (<B><FONT COLOR="#BC8F8F">&quot;matrix insertion&quot;</FONT></B>);
  
        dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
        system.matrix-&gt;add_matrix (Ke, dof_indices);
        system.rhs-&gt;add_vector    (Fe, dof_indices);
  
        perf_log.stop_event (<B><FONT COLOR="#BC8F8F">&quot;matrix insertion&quot;</FONT></B>);
      }
  
  #endif <I><FONT COLOR="#B22222">// #ifdef ENABLE_AMR
</FONT></I>  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example  ./ex14-devel
***************************************************************
 
 EquationSystems
  n_systems()=1
   System "Laplace"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE" 
    Approximation Orders="SECOND" 
    n_dofs()=21
    n_local_dofs()=21
    n_constrained_dofs()=0
    n_vectors()=1

Beginning Solve 0
System has: 21 degrees of freedom.
Linear solver converged at step: 32, final residual: 1.19434e-18
L2-Error is: 0.0150899
H1-Error is: 0.125332
  Refining the mesh...
Beginning Solve 1
System has: 65 degrees of freedom.
Linear solver converged at step: 20, final residual: 4.1598e-18
L2-Error is: 0.00515425
H1-Error is: 0.0777803
  Refining the mesh...
Beginning Solve 2
System has: 97 degrees of freedom.
Linear solver converged at step: 22, final residual: 3.12684e-18
L2-Error is: 0.00199025
H1-Error is: 0.0494318
  Refining the mesh...
Beginning Solve 3
System has: 129 degrees of freedom.
Linear solver converged at step: 29, final residual: 3.7861e-18
L2-Error is: 0.000890215
H1-Error is: 0.0320056
  Refining the mesh...
Beginning Solve 4
System has: 161 degrees of freedom.
Linear solver converged at step: 30, final residual: 3.91646e-18
L2-Error is: 0.000559523
H1-Error is: 0.0215069
  Refining the mesh...
Beginning Solve 5
System has: 193 degrees of freedom.
Linear solver converged at step: 37, final residual: 2.81216e-18
L2-Error is: 0.000483386
H1-Error is: 0.0154837
  Refining the mesh...
Beginning Solve 6
System has: 225 degrees of freedom.
Linear solver converged at step: 42, final residual: 2.80178e-18
L2-Error is: 0.000467457
H1-Error is: 0.0123018
  Refining the mesh...
Beginning Solve 7
System has: 257 degrees of freedom.
Linear solver converged at step: 44, final residual: 3.93336e-18
L2-Error is: 0.000463656
H1-Error is: 0.0107815
  Refining the mesh...
Beginning Solve 8
System has: 341 degrees of freedom.
Linear solver converged at step: 52, final residual: 4.58993e-18
L2-Error is: 0.000301328
H1-Error is: 0.00820153
  Refining the mesh...
Beginning Solve 9
System has: 445 degrees of freedom.
Linear solver converged at step: 60, final residual: 3.05606e-18
L2-Error is: 0.000182515
H1-Error is: 0.00601168
 
***************************************************************
* Done Running Example  ./ex14-devel
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
