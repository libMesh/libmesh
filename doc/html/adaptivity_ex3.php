<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("adaptivity_ex3",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file adaptivity_ex3.C with comments: </h1> 
<div class = "comment">
<h1>Adaptivity Example 3 - Laplace Equation in the L-Shaped Domain</h1>

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
named "lshaped.xda".  In addition, an input file named "adaptivity_ex3.in"
is provided which allows the user to set several parameters for
the solution so that the problem can be re-run without a
re-compile.  The solution technique employed is to have a
refinement loop with a linear solve inside followed by a
refinement of the grid and projection of the solution to the new grid
In the final loop iteration, there is no additional
refinement after the solve.  In the input file "adaptivity_ex3.in",
the variable "max_r_steps" controls the number of refinement steps,
"max_r_level" controls the maximum element refinement level, and
"refine_percentage" / "coarsen_percentage" determine the number of
elements which will be refined / coarsened at each step.


<br><br>LibMesh include files.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/mesh.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/linear_implicit_system.h"
        #include "libmesh/exodusII_io.h"
        #include "libmesh/tecplot_io.h"
        #include "libmesh/fe.h"
        #include "libmesh/quadrature_gauss.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/dense_vector.h"
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/mesh_refinement.h"
        #include "libmesh/error_vector.h"
        #include "libmesh/exact_error_estimator.h"
        #include "libmesh/kelly_error_estimator.h"
        #include "libmesh/patch_recovery_error_estimator.h"
        #include "libmesh/uniform_refinement_estimator.h"
        #include "libmesh/hp_coarsentest.h"
        #include "libmesh/hp_singular.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/mesh_modification.h"
        #include "libmesh/perf_log.h"
        #include "libmesh/getpot.h"
        #include "libmesh/exact_solution.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/elem.h"
        #include "libmesh/string_to_enum.h"
        
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
          LibMeshInit init (argc, argv);
        
</pre>
</div>
<div class = "comment">
Skip adaptive examples on a non-adaptive libMesh build
</div>

<div class ="fragment">
<pre>
        #ifndef LIBMESH_ENABLE_AMR
          libmesh_example_assert(false, "--enable-amr");
        #else
        
</pre>
</div>
<div class = "comment">
Parse the input file
</div>

<div class ="fragment">
<pre>
          GetPot input_file("adaptivity_ex3.in");
        
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
          const bool output_intermediate    = input_file("output_intermediate", false);
          dim = input_file("dimension", 2);
          const std::string indicator_type = input_file("indicator_type", "kelly");
          singularity = input_file("singularity", true);
        
</pre>
</div>
<div class = "comment">
Skip higher-dimensional examples on a lower-dimensional libMesh build
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(dim &lt;= LIBMESH_DIM, "2D/3D support");
        
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
Create a mesh, with dimension to be overridden later, on the
default MPI communicator.
</div>

<div class ="fragment">
<pre>
          Mesh mesh(init.comm());
        
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
            std::pow(TOLERANCE, 2.5);
        
</pre>
</div>
<div class = "comment">
Prints information about the system to the screen.
</div>

<div class ="fragment">
<pre>
          equation_systems.print_info();
        
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
        
        #ifdef LIBMESH_HAVE_EXODUS_API
</pre>
</div>
<div class = "comment">
After solving the system write the solution
to a ExodusII-formatted plot file.
</div>

<div class ="fragment">
<pre>
              if (output_intermediate)
                {
                  std::ostringstream outfile;
                  outfile &lt;&lt; "lshaped_" &lt;&lt; r_step &lt;&lt; ".e";
                  ExodusII_IO (mesh).write_equation_systems (outfile.str(),
                                                       equation_systems);
                }
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
        
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
We optimize in H1 norm, the default
error_estimator.error_norm = H1;


<br><br>Compute the error for each active element using
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
                          libmesh_assert_equal_to (indicator_type, "kelly");
        
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
Write out the error distribution
</div>

<div class ="fragment">
<pre>
                      std::ostringstream ss;
        	      ss &lt;&lt; r_step;
        #ifdef LIBMESH_HAVE_EXODUS_API
        	      std::string error_output = "error_"+ss.str()+".e";
        #else
        	      std::string error_output = "error_"+ss.str()+".gmv";
        #endif
                      error.plot_error( error_output, mesh );
        
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
                          libmesh_assert (singularity);
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
        
        #ifdef LIBMESH_HAVE_EXODUS_API
</pre>
</div>
<div class = "comment">
Write out the solution
After solving the system write the solution
to a ExodusII-formatted plot file.
</div>

<div class ="fragment">
<pre>
          ExodusII_IO (mesh).write_equation_systems ("lshaped.e",
                                               equation_systems);
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
        
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
        #endif // #ifndef LIBMESH_ENABLE_AMR
        
</pre>
</div>
<div class = "comment">
All done.
</div>

<div class ="fragment">
<pre>
          return 0;
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
              libmesh_assert_not_equal_to (x, 0.);
        
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
        #ifdef LIBMESH_ENABLE_AMR
</pre>
</div>
<div class = "comment">
It is a good idea to make sure we are assembling
the proper system.
</div>

<div class ="fragment">
<pre>
          libmesh_assert_equal_to (system_name, "Laplace");
        
        
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
          const MeshBase& mesh = es.get_mesh();
        
</pre>
</div>
<div class = "comment">
The dimension that we are running
</div>

<div class ="fragment">
<pre>
          const unsigned int mesh_dim = mesh.mesh_dimension();
        
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
          AutoPtr&lt;FEBase&gt; fe      (FEBase::build(mesh_dim, fe_type));
          AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(mesh_dim, fe_type));
        
</pre>
</div>
<div class = "comment">
Quadrature rules for numerical integration.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;QBase&gt; qrule(fe_type.default_quadrature_rule(mesh_dim));
          AutoPtr&lt;QBase&gt; qface(fe_type.default_quadrature_rule(mesh_dim-1));
        
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
          std::vector&lt;dof_id_type&gt; dof_indices;
        
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
              perf_log.push("elem init");
        
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
              perf_log.pop("elem init");
        
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
              perf_log.push ("Ke");
        
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
              if (mesh_dim == 1)
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
              perf_log.pop ("Ke");
        
        
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
                perf_log.push ("BCs");
        
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
                perf_log.pop ("BCs");
              }
        
        
</pre>
</div>
<div class = "comment">
The element matrix and right-hand-side are now built
for this element.  Add them to the global matrix and
right-hand-side vector.  The \p SparseMatrix::add_matrix()
and \p NumericVector::add_vector() members do this for us.
Start logging the insertion of the local (element)
matrix and vector into the global matrix and vector
</div>

<div class ="fragment">
<pre>
              perf_log.push ("matrix insertion");
        
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
              perf_log.pop ("matrix insertion");
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
        #endif // #ifdef LIBMESH_ENABLE_AMR
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file adaptivity_ex3.C without comments: </h1> 
<pre> 
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/tecplot_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exact_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/kelly_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/patch_recovery_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/uniform_refinement_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/hp_coarsentest.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/hp_singular.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_modification.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/perf_log.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exact_solution.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/string_to_enum.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_laplace(EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  
  Number exact_solution(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                        <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,   <I><FONT COLOR="#B22222">// EquationSystem parameters, not needed
</FONT></I>                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,  <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;); <I><FONT COLOR="#B22222">// unk_name, not needed);
</FONT></I>  
  Gradient exact_derivative(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                            <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,   <I><FONT COLOR="#B22222">// EquationSystems parameters, not needed
</FONT></I>                            <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,  <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>                            <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;); <I><FONT COLOR="#B22222">// unk_name, not needed);
</FONT></I>  
  
  
  <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = 2;
  
  <B><FONT COLOR="#228B22">bool</FONT></B> singularity = true;
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main(<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #ifndef LIBMESH_ENABLE_AMR
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-amr&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    GetPot input_file(<B><FONT COLOR="#BC8F8F">&quot;adaptivity_ex3.in&quot;</FONT></B>);
  
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
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> output_intermediate    = input_file(<B><FONT COLOR="#BC8F8F">&quot;output_intermediate&quot;</FONT></B>, false);
    dim = input_file(<B><FONT COLOR="#BC8F8F">&quot;dimension&quot;</FONT></B>, 2);
    <B><FONT COLOR="#228B22">const</FONT></B> std::string indicator_type = input_file(<B><FONT COLOR="#BC8F8F">&quot;indicator_type&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;kelly&quot;</FONT></B>);
    singularity = input_file(<B><FONT COLOR="#BC8F8F">&quot;singularity&quot;</FONT></B>, true);
  
    libmesh_example_assert(dim &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D/3D support&quot;</FONT></B>);
  
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
  
    Mesh mesh(init.comm());
  
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
  
    LinearImplicitSystem&amp; system =
      equation_systems.add_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>);
  
    system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, static_cast&lt;Order&gt;(approx_order),
                        <B><FONT COLOR="#5F9EA0">Utility</FONT></B>::string_to_enum&lt;FEFamily&gt;(approx_type));
  
    system.attach_assemble_function (assemble_laplace);
  
    equation_systems.init();
  
    equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;linear solver maximum iterations&quot;</FONT></B>)
      = max_linear_iterations;
  
    equation_systems.parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;linear solver tolerance&quot;</FONT></B>) =
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::pow(TOLERANCE, 2.5);
  
    equation_systems.print_info();
  
    ExactSolution exact_sol(equation_systems);
    exact_sol.attach_exact_value(exact_solution);
    exact_sol.attach_exact_deriv(exact_derivative);
  
    exact_sol.extra_quadrature_order(extra_error_quadrature);
  
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
  
  #ifdef LIBMESH_HAVE_EXODUS_API
        <B><FONT COLOR="#A020F0">if</FONT></B> (output_intermediate)
          {
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream outfile;
            outfile &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;lshaped_&quot;</FONT></B> &lt;&lt; r_step &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.e&quot;</FONT></B>;
            ExodusII_IO (mesh).write_equation_systems (outfile.str(),
                                                 equation_systems);
          }
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>  
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
                    libmesh_assert_equal_to (indicator_type, <B><FONT COLOR="#BC8F8F">&quot;kelly&quot;</FONT></B>);
  
                    KellyErrorEstimator error_estimator;
  
                    error_estimator.estimate_error (system, error);
                  }
  
                <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream ss;
  	      ss &lt;&lt; r_step;
  #ifdef LIBMESH_HAVE_EXODUS_API
  	      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string error_output = <B><FONT COLOR="#BC8F8F">&quot;error_&quot;</FONT></B>+ss.str()+<B><FONT COLOR="#BC8F8F">&quot;.e&quot;</FONT></B>;
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  	      <B><FONT COLOR="#5F9EA0">std</FONT></B>::string error_output = <B><FONT COLOR="#BC8F8F">&quot;error_&quot;</FONT></B>+ss.str()+<B><FONT COLOR="#BC8F8F">&quot;.gmv&quot;</FONT></B>;
  #endif
                error.plot_error( error_output, mesh );
  
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
                    libmesh_assert (singularity);
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
  
  #ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO (mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;lshaped.e&quot;</FONT></B>,
                                         equation_systems);
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>  
    out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;];&quot;</FONT></B> &lt;&lt; std::endl;
    out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;hold on&quot;</FONT></B> &lt;&lt; std::endl;
    out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;plot(e(:,1), e(:,2), 'bo-');&quot;</FONT></B> &lt;&lt; std::endl;
    out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;plot(e(:,1), e(:,3), 'ro-');&quot;</FONT></B> &lt;&lt; std::endl;
    out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;xlabel('dofs');&quot;</FONT></B> &lt;&lt; std::endl;
    out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;title('&quot;</FONT></B> &lt;&lt; approx_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; elements');&quot;</FONT></B> &lt;&lt; std::endl;
    out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;legend('L2-error', 'H1-error');&quot;</FONT></B> &lt;&lt; std::endl;
  #endif <I><FONT COLOR="#B22222">// #ifndef LIBMESH_ENABLE_AMR
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  
  
  
  Number exact_solution(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                        <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,  <I><FONT COLOR="#B22222">// parameters, not needed
</FONT></I>                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;) <I><FONT COLOR="#B22222">// unk_name, not needed
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
</FONT></I>                            <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>                            <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;) <I><FONT COLOR="#B22222">// unk_name, not needed
</FONT></I>  {
    Gradient gradu;
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real x = p(0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real y = dim &gt; 1 ? p(1) : 0.;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (singularity)
      {
        libmesh_assert_not_equal_to (x, 0.);
  
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
  #ifdef LIBMESH_ENABLE_AMR
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>);
  
  
    PerfLog perf_log (<B><FONT COLOR="#BC8F8F">&quot;Matrix Assembly&quot;</FONT></B>,false);
  
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> mesh_dim = mesh.mesh_dimension();
  
    LinearImplicitSystem&amp; system = es.get_system&lt;LinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Laplace&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
  
    FEType fe_type = dof_map.variable_type(0);
  
    AutoPtr&lt;FEBase&gt; fe      (FEBase::build(mesh_dim, fe_type));
    AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(mesh_dim, fe_type));
  
    AutoPtr&lt;QBase&gt; qrule(fe_type.default_quadrature_rule(mesh_dim));
    AutoPtr&lt;QBase&gt; qface(fe_type.default_quadrature_rule(mesh_dim-1));
  
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
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        perf_log.push(<B><FONT COLOR="#BC8F8F">&quot;elem init&quot;</FONT></B>);
  
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
  
        fe-&gt;reinit (elem);
  
        Ke.resize (dof_indices.size(),
                   dof_indices.size());
  
        Fe.resize (dof_indices.size());
  
        perf_log.pop(<B><FONT COLOR="#BC8F8F">&quot;elem init&quot;</FONT></B>);
  
        perf_log.push (<B><FONT COLOR="#BC8F8F">&quot;Ke&quot;</FONT></B>);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule-&gt;n_points(); qp++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;dphi.size(); i++)
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;dphi.size(); j++)
              Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (mesh_dim == 1)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule-&gt;n_points(); qp++)
            {
              Real x = q_point[qp](0);
              Real f = singularity ? sqrt(3.)/9.*pow(-x, -4./3.) :
                                     cos(x);
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;dphi.size(); ++i)
                Fe(i) += JxW[qp]*phi[i][qp]*f;
            }
  
        perf_log.pop (<B><FONT COLOR="#BC8F8F">&quot;Ke&quot;</FONT></B>);
  
  
        {
          perf_log.push (<B><FONT COLOR="#BC8F8F">&quot;BCs&quot;</FONT></B>);
  
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
  
          perf_log.pop (<B><FONT COLOR="#BC8F8F">&quot;BCs&quot;</FONT></B>);
        }
  
  
        perf_log.push (<B><FONT COLOR="#BC8F8F">&quot;matrix insertion&quot;</FONT></B>);
  
        dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
        system.matrix-&gt;add_matrix (Ke, dof_indices);
        system.rhs-&gt;add_vector    (Fe, dof_indices);
  
        perf_log.pop (<B><FONT COLOR="#BC8F8F">&quot;matrix insertion&quot;</FONT></B>);
      }
  
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_ENABLE_AMR
</FONT></I>  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex3'
***************************************************************
* Running Example adaptivity_ex3:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
 EquationSystems
  n_systems()=1
   System #0, "Laplace"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" 
    n_dofs()=21
    n_local_dofs()=9
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.66667
      Average Off-Processor Bandwidth <= 3.14286
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 12
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

Beginning Solve 0
System has: 21 degrees of freedom.
Linear solver converged at step: 10, final residual: 6.95538e-16
L2-Error is: 0.0150899
H1-Error is: 0.125332
  Refining the mesh...
Beginning Solve 1
System has: 65 degrees of freedom.
Linear solver converged at step: 21, final residual: 1.29417e-15
L2-Error is: 0.00515425
H1-Error is: 0.0777803
  Refining the mesh...
Beginning Solve 2
System has: 97 degrees of freedom.
Linear solver converged at step: 26, final residual: 1.91529e-15
L2-Error is: 0.00199025
H1-Error is: 0.0494318
  Refining the mesh...
Beginning Solve 3
System has: 129 degrees of freedom.
Linear solver converged at step: 29, final residual: 1.15387e-15
L2-Error is: 0.000890215
H1-Error is: 0.0320056
  Refining the mesh...
Beginning Solve 4
System has: 161 degrees of freedom.
Linear solver converged at step: 29, final residual: 1.71187e-15
L2-Error is: 0.000559523
H1-Error is: 0.0215069
  Refining the mesh...
Beginning Solve 5
System has: 193 degrees of freedom.
Linear solver converged at step: 30, final residual: 1.35959e-15
L2-Error is: 0.000483386
H1-Error is: 0.0154837
  Refining the mesh...
Beginning Solve 6
System has: 225 degrees of freedom.
Linear solver converged at step: 37, final residual: 4.81002e-15
L2-Error is: 0.000467457
H1-Error is: 0.0123018
  Refining the mesh...
Beginning Solve 7
System has: 257 degrees of freedom.
Linear solver converged at step: 38, final residual: 2.18805e-15
L2-Error is: 0.000463656
H1-Error is: 0.0107815
  Refining the mesh...
Beginning Solve 8
System has: 341 degrees of freedom.
Linear solver converged at step: 36, final residual: 5.80633e-15
L2-Error is: 0.000301328
H1-Error is: 0.00820153
  Refining the mesh...
Beginning Solve 9
System has: 445 degrees of freedom.
Linear solver converged at step: 43, final residual: 6.00646e-15
L2-Error is: 0.000182515
H1-Error is: 0.00601168

 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:45:42 2013                                                                          |
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
| libMesh Performance: Alive time=0.722761, Active time=0.296766                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     19        0.0023      0.000122    0.0046      0.000243    0.78     1.56     |
|   build_constraint_matrix()        115       0.0006      0.000005    0.0006      0.000005    0.21     0.21     |
|   build_sparsity()                 10        0.0019      0.000186    0.0049      0.000493    0.63     1.66     |
|   cnstrn_elem_mat_vec()            115       0.0017      0.000014    0.0017      0.000014    0.56     0.56     |
|   create_dof_constraints()         19        0.0101      0.000529    0.0188      0.000987    3.39     6.32     |
|   distribute_dofs()                19        0.0071      0.000375    0.0207      0.001088    2.40     6.97     |
|   dof_indices()                    2305      0.0160      0.000007    0.0160      0.000007    5.39     5.39     |
|   enforce_constraints_exactly()    8         0.0012      0.000152    0.0012      0.000152    0.41     0.41     |
|   old_dof_indices()                236       0.0024      0.000010    0.0024      0.000010    0.79     0.79     |
|   prepare_send_list()              19        0.0001      0.000005    0.0001      0.000005    0.03     0.03     |
|   reinit()                         19        0.0076      0.000402    0.0076      0.000402    2.58     2.58     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0002      0.000181    0.0007      0.000695    0.06     0.23     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0463      0.046302    0.0463      0.046302    15.60    15.60    |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        670       0.0055      0.000008    0.0055      0.000008    1.86     1.86     |
|   init_shape_functions()           452       0.0034      0.000008    0.0034      0.000008    1.15     1.15     |
|   inverse_map()                    3947      0.0137      0.000003    0.0137      0.000003    4.63     4.63     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             670       0.0036      0.000005    0.0036      0.000005    1.22     1.22     |
|   compute_face_map()               242       0.0028      0.000011    0.0073      0.000030    0.94     2.46     |
|   init_face_shape_functions()      18        0.0001      0.000005    0.0001      0.000005    0.03     0.03     |
|   init_reference_to_physical_map() 452       0.0059      0.000013    0.0059      0.000013    1.99     1.99     |
|                                                                                                                |
| JumpErrorEstimator                                                                                             |
|   estimate_error()                 9         0.0069      0.000769    0.0306      0.003404    2.33     10.32    |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           740       0.0013      0.000002    0.0013      0.000002    0.44     0.44     |
|   init()                           18        0.0020      0.000111    0.0020      0.000111    0.67     0.67     |
|                                                                                                                |
| Mesh                                                                                                           |
|   all_first_order()                9         0.0009      0.000097    0.0009      0.000097    0.29     0.29     |
|   all_second_order()               1         0.0001      0.000055    0.0001      0.000055    0.02     0.02     |
|   contract()                       9         0.0003      0.000028    0.0008      0.000090    0.08     0.27     |
|   find_neighbors()                 29        0.0068      0.000234    0.0089      0.000308    2.29     3.01     |
|   renumber_nodes_and_elem()        9         0.0006      0.000062    0.0006      0.000062    0.19     0.19     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   broadcast()                      1         0.0002      0.000169    0.0005      0.000488    0.06     0.16     |
|   compute_hilbert_indices()        21        0.0035      0.000168    0.0035      0.000168    1.19     1.19     |
|   find_global_indices()            21        0.0012      0.000056    0.0096      0.000457    0.39     3.23     |
|   parallel_sort()                  21        0.0016      0.000077    0.0024      0.000116    0.54     0.82     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0001      0.000053    0.0471      0.047113    0.02     15.88    |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              18        0.0002      0.000013    0.0003      0.000017    0.08     0.11     |
|   _refine_elements()               18        0.0022      0.000122    0.0054      0.000301    0.74     1.82     |
|   add_point()                      740       0.0016      0.000002    0.0030      0.000004    0.55     1.01     |
|   make_coarsening_compatible()     36        0.0046      0.000127    0.0046      0.000127    1.54     1.54     |
|   make_flags_parallel_consistent() 27        0.0026      0.000098    0.0104      0.000387    0.89     3.52     |
|   make_refinement_compatible()     36        0.0003      0.000009    0.0005      0.000013    0.10     0.16     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      20        0.0145      0.000727    0.0255      0.001274    4.90     8.59     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      99        0.0008      0.000008    0.0011      0.000011    0.28     0.38     |
|   broadcast()                      27        0.0001      0.000003    0.0001      0.000002    0.03     0.02     |
|   max(bool)                        92        0.0005      0.000006    0.0005      0.000006    0.18     0.18     |
|   max(scalar)                      2783      0.0110      0.000004    0.0110      0.000004    3.72     3.72     |
|   max(vector)                      665       0.0046      0.000007    0.0119      0.000018    1.56     4.01     |
|   min(bool)                        3445      0.0120      0.000003    0.0120      0.000003    4.04     4.04     |
|   min(scalar)                      2711      0.0128      0.000005    0.0128      0.000005    4.32     4.32     |
|   min(vector)                      665       0.0050      0.000008    0.0124      0.000019    1.69     4.18     |
|   probe()                          912       0.0024      0.000003    0.0024      0.000003    0.82     0.82     |
|   receive()                        912       0.0028      0.000003    0.0054      0.000006    0.95     1.81     |
|   send()                           912       0.0017      0.000002    0.0017      0.000002    0.57     0.57     |
|   send_receive()                   954       0.0038      0.000004    0.0121      0.000013    1.28     4.08     |
|   sum()                            143       0.0011      0.000008    0.0023      0.000016    0.38     0.78     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           912       0.0009      0.000001    0.0009      0.000001    0.32     0.32     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         29        0.0036      0.000123    0.0079      0.000274    1.21     2.68     |
|   set_parent_processor_ids()       20        0.0005      0.000027    0.0005      0.000027    0.18     0.18     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          10        0.0383      0.003828    0.0383      0.003828    12.90    12.90    |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       9         0.0013      0.000143    0.0050      0.000558    0.43     1.69     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       10        0.0065      0.000652    0.0159      0.001586    2.20     5.34     |
|   project_vector()                 9         0.0029      0.000322    0.0115      0.001277    0.98     3.87     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            26440     0.2968                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example adaptivity_ex3:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex3'
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
