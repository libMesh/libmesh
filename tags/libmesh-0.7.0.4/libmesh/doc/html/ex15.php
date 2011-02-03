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

<br><br>This example solves the Biharmonic equation on a square or cube,
using a Galerkin formulation with C1 elements approximating the
H^2_0 function space.
The initial mesh contains two TRI6, one QUAD9 or one HEX27
An input file named "ex15.in"
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
        #include "exodusII_io.h"
        #include "fe.h"
        #include "quadrature.h"
        #include "dense_matrix.h"
        #include "dense_vector.h"
        #include "sparse_matrix.h"
        #include "mesh_generation.h"
        #include "mesh_modification.h"
        #include "mesh_refinement.h"
        #include "error_vector.h"
        #include "fourth_error_estimators.h"
        #include "getpot.h"
        #include "exact_solution.h"
        #include "dof_map.h"
        #include "numeric_vector.h"
        #include "elem.h"
        #include "tensor_value.h"
        #include "perf_log.h"
        
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
Prototypes for calculation of the exact solution.  Necessary
for setting boundary conditions.
</div>

<div class ="fragment">
<pre>
        Number exact_2D_solution(const Point& p,
                                 const Parameters&,   // parameters, not needed
                                 const std::string&,  // sys_name, not needed
                                 const std::string&); // unk_name, not needed);
        
        Number exact_3D_solution(const Point& p,
          const Parameters&, const std::string&, const std::string&);
        
</pre>
</div>
<div class = "comment">
Prototypes for calculation of the gradient of the exact solution.  
Necessary for setting boundary conditions in H^2_0 and testing
H^1 convergence of the solution
</div>

<div class ="fragment">
<pre>
        Gradient exact_2D_derivative(const Point& p,
          const Parameters&, const std::string&, const std::string&);
        
        Gradient exact_3D_derivative(const Point& p,
          const Parameters&, const std::string&, const std::string&);
        
        Tensor exact_2D_hessian(const Point& p,
          const Parameters&, const std::string&, const std::string&);
        
        Tensor exact_3D_hessian(const Point& p,
          const Parameters&, const std::string&, const std::string&);
        
        Number forcing_function_2D(const Point& p);
        
        Number forcing_function_3D(const Point& p);
        
</pre>
</div>
<div class = "comment">
Pointers to dimension-independent functions
</div>

<div class ="fragment">
<pre>
        Number (*exact_solution)(const Point& p,
          const Parameters&, const std::string&, const std::string&);
        Gradient (*exact_derivative)(const Point& p,
          const Parameters&, const std::string&, const std::string&);
        Tensor (*exact_hessian)(const Point& p,
          const Parameters&, const std::string&, const std::string&);
        Number (*forcing_function)(const Point& p);
        
        
        
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
This example requires Adaptive Mesh Refinement support
</div>

<div class ="fragment">
<pre>
        #ifndef LIBMESH_ENABLE_AMR
          libmesh_example_assert(false, "--enable-amr");
        #else
        
</pre>
</div>
<div class = "comment">
This example requires second derivative calculation support
</div>

<div class ="fragment">
<pre>
        #ifndef LIBMESH_ENABLE_SECOND_DERIVATIVES
          libmesh_example_assert(false, "--enable-second");
        #else
        
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
          const unsigned int max_r_level = input_file("max_r_level", 10);
          const unsigned int max_r_steps = input_file("max_r_steps", 4);
          const std::string approx_type  = input_file("approx_type",
                                                      "HERMITE");
          const unsigned int uniform_refine =
                          input_file("uniform_refine", 0);
          const Real refine_percentage =
                          input_file("refine_percentage", 0.5);
          const Real coarsen_percentage =
                          input_file("coarsen_percentage", 0.5);
          const unsigned int dim =
                          input_file("dimension", 2);
          const unsigned int max_linear_iterations =
                          input_file("max_linear_iterations", 10000);
        
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
We have only defined 2 and 3 dimensional problems
</div>

<div class ="fragment">
<pre>
          libmesh_assert (dim == 2 || dim == 3);
        
</pre>
</div>
<div class = "comment">
Currently only the Hermite cubics give a 3D C^1 basis
</div>

<div class ="fragment">
<pre>
          libmesh_assert (dim == 2 || approx_type == "HERMITE");
        
</pre>
</div>
<div class = "comment">
Create a mesh.
</div>

<div class ="fragment">
<pre>
          Mesh mesh;
          
</pre>
</div>
<div class = "comment">
Output file for plotting the error 
</div>

<div class ="fragment">
<pre>
          std::string output_file = "";
        
          if (dim == 2)
            output_file += "2D_";
          else if (dim == 3)
            output_file += "3D_";
        
          if (approx_type == "HERMITE")
            output_file += "hermite_";
          else if (approx_type == "SECOND")
            output_file += "reducedclough_";
          else
            output_file += "clough_";
        
          if (uniform_refine == 0)
            output_file += "adaptive";
          else
            output_file += "uniform";
        
          std::string exd_file = output_file;
          exd_file += ".exd";
          output_file += ".m";
        
          std::ofstream out (output_file.c_str());
          out &lt;&lt; "% dofs     L2-error     H1-error      H2-error\n"
              &lt;&lt; "e = [\n";
          
</pre>
</div>
<div class = "comment">
Set up the dimension-dependent coarse mesh and solution
We build more than one cell so as to avoid bugs on fewer than 
4 processors in 2D or 8 in 3D.
</div>

<div class ="fragment">
<pre>
          if (dim == 2)
            {
              MeshTools::Generation::build_square(mesh, 2, 2);
              exact_solution = &exact_2D_solution;
              exact_derivative = &exact_2D_derivative;
              exact_hessian = &exact_2D_hessian;
              forcing_function = &forcing_function_2D;
            }
          else if (dim == 3)
            {
              MeshTools::Generation::build_cube(mesh, 2, 2, 2);
              exact_solution = &exact_3D_solution;
              exact_derivative = &exact_3D_derivative;
              exact_hessian = &exact_3D_hessian;
              forcing_function = &forcing_function_3D;
            }
        
</pre>
</div>
<div class = "comment">
Convert the mesh to second order: necessary for computing with
Clough-Tocher elements, useful for getting slightly less 
broken visualization output with Hermite elements
</div>

<div class ="fragment">
<pre>
          mesh.all_second_order();
        
</pre>
</div>
<div class = "comment">
Convert it to triangles if necessary
</div>

<div class ="fragment">
<pre>
          if (approx_type != "HERMITE")
            MeshTools::Modification::all_tri(mesh);
        
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
Create a system named "Biharmonic"
</div>

<div class ="fragment">
<pre>
          LinearImplicitSystem& system =
            equation_systems.add_system&lt;LinearImplicitSystem&gt; ("Biharmonic");
        
</pre>
</div>
<div class = "comment">
Adds the variable "u" to "Biharmonic".  "u"
will be approximated using Hermite tensor product squares
or (possibly reduced) cubic Clough-Tocher triangles
</div>

<div class ="fragment">
<pre>
          if (approx_type == "HERMITE")
            system.add_variable("u", THIRD, HERMITE);
          else if (approx_type == "SECOND")
            system.add_variable("u", SECOND, CLOUGH);
          else if (approx_type == "CLOUGH")
            system.add_variable("u", THIRD, CLOUGH);
          else
            libmesh_error();
        
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly
function.
</div>

<div class ="fragment">
<pre>
          system.attach_assemble_function (assemble_biharmonic);
              
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
            ("linear solver maximum iterations") = max_linear_iterations;
        
</pre>
</div>
<div class = "comment">
Linear solver tolerance.
</div>

<div class ="fragment">
<pre>
          equation_systems.parameters.set&lt;Real&gt;
            ("linear solver tolerance") = TOLERANCE * TOLERANCE;
              
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
Construct ExactSolution object and attach function to compute exact solution
</div>

<div class ="fragment">
<pre>
          ExactSolution exact_sol(equation_systems);
          exact_sol.attach_exact_value(exact_solution);
          exact_sol.attach_exact_deriv(exact_derivative);
          exact_sol.attach_exact_hessian(exact_hessian);
        
</pre>
</div>
<div class = "comment">
Construct zero solution object, useful for computing solution norms
Attaching "zero_solution" functions is unnecessary
</div>

<div class ="fragment">
<pre>
          ExactSolution zero_sol(equation_systems);
        
</pre>
</div>
<div class = "comment">
A refinement loop.
</div>

<div class ="fragment">
<pre>
          for (unsigned int r_step=0; r_step&lt;max_r_steps; r_step++)
            {
              mesh.print_info();
              equation_systems.print_info();
        
              std::cout &lt;&lt; "Beginning Solve " &lt;&lt; r_step &lt;&lt; std::endl;
              
</pre>
</div>
<div class = "comment">
Solve the system "Biharmonic", just like example 2.
</div>

<div class ="fragment">
<pre>
              system.solve();
        
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
Compute the norm.
</div>

<div class ="fragment">
<pre>
              zero_sol.compute_error("Biharmonic", "u");
        
</pre>
</div>
<div class = "comment">
Print out the error values
</div>

<div class ="fragment">
<pre>
              std::cout &lt;&lt; "L2-Norm is: "
                        &lt;&lt; zero_sol.l2_error("Biharmonic", "u")
                        &lt;&lt; std::endl;
              std::cout &lt;&lt; "H1-Norm is: "
                        &lt;&lt; zero_sol.h1_error("Biharmonic", "u")
                        &lt;&lt; std::endl;
              std::cout &lt;&lt; "H2-Norm is: "
                        &lt;&lt; zero_sol.h2_error("Biharmonic", "u")
                        &lt;&lt; std::endl
                        &lt;&lt; std::endl;
              std::cout &lt;&lt; "L2-Error is: "
                        &lt;&lt; exact_sol.l2_error("Biharmonic", "u")
                        &lt;&lt; std::endl;
              std::cout &lt;&lt; "H1-Error is: "
                        &lt;&lt; exact_sol.h1_error("Biharmonic", "u")
                        &lt;&lt; std::endl;
              std::cout &lt;&lt; "H2-Error is: "
                        &lt;&lt; exact_sol.h2_error("Biharmonic", "u")
                        &lt;&lt; std::endl
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
                  &lt;&lt; exact_sol.h1_error("Biharmonic", "u") &lt;&lt; " "
                  &lt;&lt; exact_sol.h2_error("Biharmonic", "u") &lt;&lt; std::endl;
        
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
                      ErrorVector error;
                      LaplacianErrorEstimator error_estimator;
        
                      error_estimator.estimate_error(system, error);
                      mesh_refinement.flag_elements_by_elem_fraction (error);
        
                      std::cout &lt;&lt; "Mean Error: " &lt;&lt; error.mean() &lt;&lt;
                                      std::endl;
                      std::cout &lt;&lt; "Error Variance: " &lt;&lt; error.variance() &lt;&lt;
                                      std::endl;
        
                      mesh_refinement.refine_and_coarsen_elements();
                    }
                  else
                    {
                      mesh_refinement.uniformly_refine(1);
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
          ExodusII_IO (mesh).write_equation_systems (exd_file,
                                               equation_systems);
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
        
</pre>
</div>
<div class = "comment">
Close up the output file.
</div>

<div class ="fragment">
<pre>
          out &lt;&lt; "];\n"
              &lt;&lt; "hold on\n"
              &lt;&lt; "loglog(e(:,1), e(:,2), 'bo-');\n"
              &lt;&lt; "loglog(e(:,1), e(:,3), 'ro-');\n"
              &lt;&lt; "loglog(e(:,1), e(:,4), 'go-');\n"
              &lt;&lt; "xlabel('log(dofs)');\n"
              &lt;&lt; "ylabel('log(error)');\n"
              &lt;&lt; "title('C1 " &lt;&lt; approx_type &lt;&lt; " elements');\n"
              &lt;&lt; "legend('L2-error', 'H1-error', 'H2-error');\n";
          
</pre>
</div>
<div class = "comment">
All done.  
</div>

<div class ="fragment">
<pre>
          return 0;
        #endif // #ifndef LIBMESH_ENABLE_SECOND_DERIVATIVES
        #endif // #ifndef LIBMESH_ENABLE_AMR
        }
        
        
        
        Number exact_2D_solution(const Point& p,
                                 const Parameters&,  // parameters, not needed
                                 const std::string&, // sys_name, not needed
                                 const std::string&) // unk_name, not needed
        {
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
analytic solution value
</div>

<div class ="fragment">
<pre>
          return 256.*(x-x*x)*(x-x*x)*(y-y*y)*(y-y*y);
        }
        
        
</pre>
</div>
<div class = "comment">
We now define the gradient of the exact solution
</div>

<div class ="fragment">
<pre>
        Gradient exact_2D_derivative(const Point& p,
                                     const Parameters&,  // parameters, not needed
                                     const std::string&, // sys_name, not needed
                                     const std::string&) // unk_name, not needed
        {
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
First derivatives to be returned.
</div>

<div class ="fragment">
<pre>
          Gradient gradu;
        
          gradu(0) = 256.*2.*(x-x*x)*(1-2*x)*(y-y*y)*(y-y*y);
          gradu(1) = 256.*2.*(x-x*x)*(x-x*x)*(y-y*y)*(1-2*y);
        
          return gradu;
        }
        
        
</pre>
</div>
<div class = "comment">
We now define the hessian of the exact solution
</div>

<div class ="fragment">
<pre>
        Tensor exact_2D_hessian(const Point& p,
                                const Parameters&,  // parameters, not needed
                                const std::string&, // sys_name, not needed
                                const std::string&) // unk_name, not needed
        {
</pre>
</div>
<div class = "comment">
Second derivatives to be returned.
</div>

<div class ="fragment">
<pre>
          Tensor hessu;
          
</pre>
</div>
<div class = "comment">
x and y coordinates in space
</div>

<div class ="fragment">
<pre>
          const Real x = p(0);
          const Real y = p(1);
        
          hessu(0,0) = 256.*2.*(1-6.*x+6.*x*x)*(y-y*y)*(y-y*y);
          hessu(0,1) = 256.*4.*(x-x*x)*(1.-2.*x)*(y-y*y)*(1.-2.*y);
          hessu(1,1) = 256.*2.*(x-x*x)*(x-x*x)*(1.-6.*y+6.*y*y);
        
</pre>
</div>
<div class = "comment">
Hessians are always symmetric
</div>

<div class ="fragment">
<pre>
          hessu(1,0) = hessu(0,1);
          return hessu;
        }
        
        
        
        Number forcing_function_2D(const Point& p)
        {
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
Equals laplacian(laplacian(u))
</div>

<div class ="fragment">
<pre>
          return 256. * 8. * (3.*((y-y*y)*(y-y*y)+(x-x*x)*(x-x*x))
                 + (1.-6.*x+6.*x*x)*(1.-6.*y+6.*y*y));
        }
        
        
        
        Number exact_3D_solution(const Point& p,
                                 const Parameters&,  // parameters, not needed
                                 const std::string&, // sys_name, not needed
                                 const std::string&) // unk_name, not needed
        {
</pre>
</div>
<div class = "comment">
xyz coordinates in space
</div>

<div class ="fragment">
<pre>
          const Real x = p(0);
          const Real y = p(1);
          const Real z = p(2);
          
</pre>
</div>
<div class = "comment">
analytic solution value
</div>

<div class ="fragment">
<pre>
          return 4096.*(x-x*x)*(x-x*x)*(y-y*y)*(y-y*y)*(z-z*z)*(z-z*z);
        }
        
        
        Gradient exact_3D_derivative(const Point& p,
                                     const Parameters&,  // parameters, not needed
                                     const std::string&, // sys_name, not needed
                                     const std::string&) // unk_name, not needed
        {
</pre>
</div>
<div class = "comment">
First derivatives to be returned.
</div>

<div class ="fragment">
<pre>
          Gradient gradu;
          
</pre>
</div>
<div class = "comment">
xyz coordinates in space
</div>

<div class ="fragment">
<pre>
          const Real x = p(0);
          const Real y = p(1);
          const Real z = p(2);
        
          gradu(0) = 4096.*2.*(x-x*x)*(1.-2.*x)*(y-y*y)*(y-y*y)*(z-z*z)*(z-z*z);
          gradu(1) = 4096.*2.*(x-x*x)*(x-x*x)*(y-y*y)*(1.-2.*y)*(z-z*z)*(z-z*z);
          gradu(2) = 4096.*2.*(x-x*x)*(x-x*x)*(y-y*y)*(y-y*y)*(z-z*z)*(1.-2.*z);
        
          return gradu;
        }
        
        
</pre>
</div>
<div class = "comment">
We now define the hessian of the exact solution
</div>

<div class ="fragment">
<pre>
        Tensor exact_3D_hessian(const Point& p,
                                const Parameters&,  // parameters, not needed
                                const std::string&, // sys_name, not needed
                                const std::string&) // unk_name, not needed
        {
</pre>
</div>
<div class = "comment">
Second derivatives to be returned.
</div>

<div class ="fragment">
<pre>
          Tensor hessu;
          
</pre>
</div>
<div class = "comment">
xyz coordinates in space
</div>

<div class ="fragment">
<pre>
          const Real x = p(0);
          const Real y = p(1);
          const Real z = p(2);
        
          hessu(0,0) = 4096.*(2.-12.*x+12.*x*x)*(y-y*y)*(y-y*y)*(z-z*z)*(z-z*z);
          hessu(0,1) = 4096.*4.*(x-x*x)*(1.-2.*x)*(y-y*y)*(1.-2.*y)*(z-z*z)*(z-z*z);
          hessu(0,2) = 4096.*4.*(x-x*x)*(1.-2.*x)*(y-y*y)*(y-y*y)*(z-z*z)*(1.-2.*z);
          hessu(1,1) = 4096.*(x-x*x)*(x-x*x)*(2.-12.*y+12.*y*y)*(z-z*z)*(z-z*z);
          hessu(1,2) = 4096.*4.*(x-x*x)*(x-x*x)*(y-y*y)*(1.-2.*y)*(z-z*z)*(1.-2.*z);
          hessu(2,2) = 4096.*(x-x*x)*(x-x*x)*(y-y*y)*(y-y*y)*(2.-12.*z+12.*z*z);
        
</pre>
</div>
<div class = "comment">
Hessians are always symmetric
</div>

<div class ="fragment">
<pre>
          hessu(1,0) = hessu(0,1);
          hessu(2,0) = hessu(0,2);
          hessu(2,1) = hessu(1,2);
        
          return hessu;
        }
        
        
        
        Number forcing_function_3D(const Point& p)
        {
</pre>
</div>
<div class = "comment">
xyz coordinates in space
</div>

<div class ="fragment">
<pre>
          const Real x = p(0);
          const Real y = p(1);
          const Real z = p(2);
        
</pre>
</div>
<div class = "comment">
Equals laplacian(laplacian(u))
</div>

<div class ="fragment">
<pre>
          return 4096. * 8. * (3.*((y-y*y)*(y-y*y)*(x-x*x)*(x-x*x) +
                                   (z-z*z)*(z-z*z)*(x-x*x)*(x-x*x) +
                                   (z-z*z)*(z-z*z)*(y-y*y)*(y-y*y)) +
                 (1.-6.*x+6.*x*x)*(1.-6.*y+6.*y*y)*(z-z*z)*(z-z*z) +
                 (1.-6.*x+6.*x*x)*(1.-6.*z+6.*z*z)*(y-y*y)*(y-y*y) +
                 (1.-6.*y+6.*y*y)*(1.-6.*z+6.*z*z)*(x-x*x)*(x-x*x));
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
        #ifdef LIBMESH_ENABLE_AMR
        #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        
</pre>
</div>
<div class = "comment">
It is a good idea to make sure we are assembling
the proper system.
</div>

<div class ="fragment">
<pre>
          libmesh_assert (system_name == "Biharmonic");
        
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
Quadrature rule for numerical integration.
With 2D triangles, the Clough quadrature rule puts a Gaussian
quadrature rule on each of the 3 subelements
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;QBase&gt; qrule(fe_type.default_quadrature_rule(dim));
        
</pre>
</div>
<div class = "comment">
Tell the finite element object to use our quadrature rule.
</div>

<div class ="fragment">
<pre>
          fe-&gt;attach_quadrature_rule (qrule.get());
        
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
          AutoPtr&lt;QBase&gt; qface(fe_type.default_quadrature_rule(dim-1));
        
</pre>
</div>
<div class = "comment">
Tell the finte element object to use our
quadrature rule.
</div>

<div class ="fragment">
<pre>
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
For efficiency we will compute shape function laplacians n times,
not n^2
</div>

<div class ="fragment">
<pre>
          std::vector&lt;Real&gt; shape_laplacian;
        
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
Make sure there is enough room in this cache
</div>

<div class ="fragment">
<pre>
              shape_laplacian.resize(dof_indices.size());
        
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
a double loop to integrate laplacians of the test funcions
(i) against laplacians of the trial functions (j).

<br><br>This step is why we need the Clough-Tocher elements -
these C1 differentiable elements have square-integrable
second derivatives.

<br><br>Now start logging the element matrix computation
</div>

<div class ="fragment">
<pre>
              perf_log.push ("Ke");
        
              for (unsigned int qp=0; qp&lt;qrule-&gt;n_points(); qp++)
                {
                  for (unsigned int i=0; i&lt;phi.size(); i++)
                    {
                      shape_laplacian[i] = d2phi[i][qp](0,0)+d2phi[i][qp](1,1);
                      if (dim == 3)
                         shape_laplacian[i] += d2phi[i][qp](2,2);
                    }
                  for (unsigned int i=0; i&lt;phi.size(); i++)
                    for (unsigned int j=0; j&lt;phi.size(); j++)
                      Ke(i,j) += JxW[qp]*
                                 shape_laplacian[i]*shape_laplacian[j];
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
                perf_log.push ("BCs");
        
</pre>
</div>
<div class = "comment">
The penalty values, for solution boundary trace and flux.  
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
                      const std::vector&lt;Point&gt;& qface_point = fe_face-&gt;get_xyz();
        
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
                      for (unsigned int qp=0; qp&lt;qface-&gt;n_points(); qp++)
                        {
</pre>
</div>
<div class = "comment">
The boundary value.
</div>

<div class ="fragment">
<pre>
                          Number value = exact_solution(qface_point[qp],
                                                        es.parameters, "null",
                                                        "void");
                          Gradient flux = exact_2D_derivative(qface_point[qp],
                                                              es.parameters,
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
                perf_log.pop ("BCs");
              } 
        
              for (unsigned int qp=0; qp&lt;qrule-&gt;n_points(); qp++)
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
              perf_log.push ("matrix insertion");
        
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
              perf_log.pop ("matrix insertion");
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
        
        #endif // #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
        #endif // #ifdef LIBMESH_ENABLE_AMR
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_modification.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fourth_error_estimators.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exact_solution.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;tensor_value.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;perf_log.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_biharmonic(EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  
  Number exact_2D_solution(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                           <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,   <I><FONT COLOR="#B22222">// parameters, not needed
</FONT></I>                           <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,  <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>                           <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;); <I><FONT COLOR="#B22222">// unk_name, not needed);
</FONT></I>  
  Number exact_3D_solution(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
    <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;);
  
  Gradient exact_2D_derivative(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
    <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;);
  
  Gradient exact_3D_derivative(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
    <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;);
  
  Tensor exact_2D_hessian(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
    <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;);
  
  Tensor exact_3D_hessian(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
    <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;);
  
  Number forcing_function_2D(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p);
  
  Number forcing_function_3D(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p);
  
  Number (*exact_solution)(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
    <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;);
  Gradient (*exact_derivative)(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
    <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;);
  Tensor (*exact_hessian)(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
    <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;);
  Number (*forcing_function)(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p);
  
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main(<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #ifndef LIBMESH_ENABLE_AMR
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-amr&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
  #ifndef LIBMESH_ENABLE_SECOND_DERIVATIVES
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-second&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    GetPot input_file(<B><FONT COLOR="#BC8F8F">&quot;ex15.in&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_r_level = input_file(<B><FONT COLOR="#BC8F8F">&quot;max_r_level&quot;</FONT></B>, 10);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_r_steps = input_file(<B><FONT COLOR="#BC8F8F">&quot;max_r_steps&quot;</FONT></B>, 4);
    <B><FONT COLOR="#228B22">const</FONT></B> std::string approx_type  = input_file(<B><FONT COLOR="#BC8F8F">&quot;approx_type&quot;</FONT></B>,
                                                <B><FONT COLOR="#BC8F8F">&quot;HERMITE&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> uniform_refine =
                    input_file(<B><FONT COLOR="#BC8F8F">&quot;uniform_refine&quot;</FONT></B>, 0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real refine_percentage =
                    input_file(<B><FONT COLOR="#BC8F8F">&quot;refine_percentage&quot;</FONT></B>, 0.5);
    <B><FONT COLOR="#228B22">const</FONT></B> Real coarsen_percentage =
                    input_file(<B><FONT COLOR="#BC8F8F">&quot;coarsen_percentage&quot;</FONT></B>, 0.5);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim =
                    input_file(<B><FONT COLOR="#BC8F8F">&quot;dimension&quot;</FONT></B>, 2);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_linear_iterations =
                    input_file(<B><FONT COLOR="#BC8F8F">&quot;max_linear_iterations&quot;</FONT></B>, 10000);
  
    libmesh_example_assert(dim &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D/3D support&quot;</FONT></B>);
    
    libmesh_assert (dim == 2 || dim == 3);
  
    libmesh_assert (dim == 2 || approx_type == <B><FONT COLOR="#BC8F8F">&quot;HERMITE&quot;</FONT></B>);
  
    Mesh mesh;
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string output_file = <B><FONT COLOR="#BC8F8F">&quot;&quot;</FONT></B>;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 2)
      output_file += <B><FONT COLOR="#BC8F8F">&quot;2D_&quot;</FONT></B>;
    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
      output_file += <B><FONT COLOR="#BC8F8F">&quot;3D_&quot;</FONT></B>;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (approx_type == <B><FONT COLOR="#BC8F8F">&quot;HERMITE&quot;</FONT></B>)
      output_file += <B><FONT COLOR="#BC8F8F">&quot;hermite_&quot;</FONT></B>;
    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (approx_type == <B><FONT COLOR="#BC8F8F">&quot;SECOND&quot;</FONT></B>)
      output_file += <B><FONT COLOR="#BC8F8F">&quot;reducedclough_&quot;</FONT></B>;
    <B><FONT COLOR="#A020F0">else</FONT></B>
      output_file += <B><FONT COLOR="#BC8F8F">&quot;clough_&quot;</FONT></B>;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (uniform_refine == 0)
      output_file += <B><FONT COLOR="#BC8F8F">&quot;adaptive&quot;</FONT></B>;
    <B><FONT COLOR="#A020F0">else</FONT></B>
      output_file += <B><FONT COLOR="#BC8F8F">&quot;uniform&quot;</FONT></B>;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string exd_file = output_file;
    exd_file += <B><FONT COLOR="#BC8F8F">&quot;.exd&quot;</FONT></B>;
    output_file += <B><FONT COLOR="#BC8F8F">&quot;.m&quot;</FONT></B>;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream out (output_file.c_str());
    out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;% dofs     L2-error     H1-error      H2-error\n&quot;</FONT></B>
        &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;e = [\n&quot;</FONT></B>;
    
    <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 2)
      {
        <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square(mesh, 2, 2);
        exact_solution = &amp;exact_2D_solution;
        exact_derivative = &amp;exact_2D_derivative;
        exact_hessian = &amp;exact_2D_hessian;
        forcing_function = &amp;forcing_function_2D;
      }
    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
      {
        <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube(mesh, 2, 2, 2);
        exact_solution = &amp;exact_3D_solution;
        exact_derivative = &amp;exact_3D_derivative;
        exact_hessian = &amp;exact_3D_hessian;
        forcing_function = &amp;forcing_function_3D;
      }
  
    mesh.all_second_order();
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (approx_type != <B><FONT COLOR="#BC8F8F">&quot;HERMITE&quot;</FONT></B>)
      <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Modification::all_tri(mesh);
  
    MeshRefinement mesh_refinement(mesh);
    mesh_refinement.refine_fraction() = refine_percentage;
    mesh_refinement.coarsen_fraction() = coarsen_percentage;
    mesh_refinement.max_h_level() = max_r_level;
  
    EquationSystems equation_systems (mesh);
  
    LinearImplicitSystem&amp; system =
      equation_systems.add_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Biharmonic&quot;</FONT></B>);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (approx_type == <B><FONT COLOR="#BC8F8F">&quot;HERMITE&quot;</FONT></B>)
      system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, THIRD, HERMITE);
    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (approx_type == <B><FONT COLOR="#BC8F8F">&quot;SECOND&quot;</FONT></B>)
      system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, SECOND, CLOUGH);
    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (approx_type == <B><FONT COLOR="#BC8F8F">&quot;CLOUGH&quot;</FONT></B>)
      system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, THIRD, CLOUGH);
    <B><FONT COLOR="#A020F0">else</FONT></B>
      libmesh_error();
  
    system.attach_assemble_function (assemble_biharmonic);
        
    equation_systems.init();
  
    equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;
      (<B><FONT COLOR="#BC8F8F">&quot;linear solver maximum iterations&quot;</FONT></B>) = max_linear_iterations;
  
    equation_systems.parameters.set&lt;Real&gt;
      (<B><FONT COLOR="#BC8F8F">&quot;linear solver tolerance&quot;</FONT></B>) = TOLERANCE * TOLERANCE;
        
    equation_systems.print_info();
  
    ExactSolution exact_sol(equation_systems);
    exact_sol.attach_exact_value(exact_solution);
    exact_sol.attach_exact_deriv(exact_derivative);
    exact_sol.attach_exact_hessian(exact_hessian);
  
    ExactSolution zero_sol(equation_systems);
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> r_step=0; r_step&lt;max_r_steps; r_step++)
      {
        mesh.print_info();
        equation_systems.print_info();
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Beginning Solve &quot;</FONT></B> &lt;&lt; r_step &lt;&lt; std::endl;
        
        system.solve();
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Linear solver converged at step: &quot;</FONT></B>
                  &lt;&lt; system.n_linear_iterations()
                  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, final residual: &quot;</FONT></B>
                  &lt;&lt; system.final_linear_residual()
                  &lt;&lt; std::endl;
  
        exact_sol.compute_error(<B><FONT COLOR="#BC8F8F">&quot;Biharmonic&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>);
        zero_sol.compute_error(<B><FONT COLOR="#BC8F8F">&quot;Biharmonic&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>);
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;L2-Norm is: &quot;</FONT></B>
                  &lt;&lt; zero_sol.l2_error(<B><FONT COLOR="#BC8F8F">&quot;Biharmonic&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)
                  &lt;&lt; std::endl;
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;H1-Norm is: &quot;</FONT></B>
                  &lt;&lt; zero_sol.h1_error(<B><FONT COLOR="#BC8F8F">&quot;Biharmonic&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)
                  &lt;&lt; std::endl;
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;H2-Norm is: &quot;</FONT></B>
                  &lt;&lt; zero_sol.h2_error(<B><FONT COLOR="#BC8F8F">&quot;Biharmonic&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)
                  &lt;&lt; std::endl
                  &lt;&lt; std::endl;
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;L2-Error is: &quot;</FONT></B>
                  &lt;&lt; exact_sol.l2_error(<B><FONT COLOR="#BC8F8F">&quot;Biharmonic&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)
                  &lt;&lt; std::endl;
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;H1-Error is: &quot;</FONT></B>
                  &lt;&lt; exact_sol.h1_error(<B><FONT COLOR="#BC8F8F">&quot;Biharmonic&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)
                  &lt;&lt; std::endl;
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;H2-Error is: &quot;</FONT></B>
                  &lt;&lt; exact_sol.h2_error(<B><FONT COLOR="#BC8F8F">&quot;Biharmonic&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)
                  &lt;&lt; std::endl
                  &lt;&lt; std::endl;
  
        out &lt;&lt; equation_systems.n_active_dofs() &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B>
            &lt;&lt; exact_sol.l2_error(<B><FONT COLOR="#BC8F8F">&quot;Biharmonic&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>) &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B>
            &lt;&lt; exact_sol.h1_error(<B><FONT COLOR="#BC8F8F">&quot;Biharmonic&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>) &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B>
            &lt;&lt; exact_sol.h2_error(<B><FONT COLOR="#BC8F8F">&quot;Biharmonic&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>) &lt;&lt; std::endl;
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (r_step+1 != max_r_steps)
          {
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;  Refining the mesh...&quot;</FONT></B> &lt;&lt; std::endl;
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (uniform_refine == 0)
              {
                ErrorVector error;
                LaplacianErrorEstimator error_estimator;
  
                error_estimator.estimate_error(system, error);
                mesh_refinement.flag_elements_by_elem_fraction (error);
  
                <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Mean Error: &quot;</FONT></B> &lt;&lt; error.mean() &lt;&lt;
                                <B><FONT COLOR="#5F9EA0">std</FONT></B>::endl;
                <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error Variance: &quot;</FONT></B> &lt;&lt; error.variance() &lt;&lt;
                                <B><FONT COLOR="#5F9EA0">std</FONT></B>::endl;
  
                mesh_refinement.refine_and_coarsen_elements();
              }
            <B><FONT COLOR="#A020F0">else</FONT></B>
              {
                mesh_refinement.uniformly_refine(1);
              }
                
            equation_systems.reinit ();
          }
      }            
    
  #ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO (mesh).write_equation_systems (exd_file,
                                         equation_systems);
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>  
    out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;];\n&quot;</FONT></B>
        &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;hold on\n&quot;</FONT></B>
        &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;loglog(e(:,1), e(:,2), 'bo-');\n&quot;</FONT></B>
        &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;loglog(e(:,1), e(:,3), 'ro-');\n&quot;</FONT></B>
        &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;loglog(e(:,1), e(:,4), 'go-');\n&quot;</FONT></B>
        &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;xlabel('log(dofs)');\n&quot;</FONT></B>
        &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ylabel('log(error)');\n&quot;</FONT></B>
        &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;title('C1 &quot;</FONT></B> &lt;&lt; approx_type &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; elements');\n&quot;</FONT></B>
        &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;legend('L2-error', 'H1-error', 'H2-error');\n&quot;</FONT></B>;
    
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  #endif <I><FONT COLOR="#B22222">// #ifndef LIBMESH_ENABLE_SECOND_DERIVATIVES
</FONT></I>  #endif <I><FONT COLOR="#B22222">// #ifndef LIBMESH_ENABLE_AMR
</FONT></I>  }
  
  
  
  Number exact_2D_solution(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                           <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,  <I><FONT COLOR="#B22222">// parameters, not needed
</FONT></I>                           <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>                           <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;) <I><FONT COLOR="#B22222">// unk_name, not needed
</FONT></I>  {
    <B><FONT COLOR="#228B22">const</FONT></B> Real x = p(0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real y = p(1);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 256.*(x-x*x)*(x-x*x)*(y-y*y)*(y-y*y);
  }
  
  
  Gradient exact_2D_derivative(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                               <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,  <I><FONT COLOR="#B22222">// parameters, not needed
</FONT></I>                               <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>                               <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;) <I><FONT COLOR="#B22222">// unk_name, not needed
</FONT></I>  {
    <B><FONT COLOR="#228B22">const</FONT></B> Real x = p(0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real y = p(1);
  
    Gradient gradu;
  
    gradu(0) = 256.*2.*(x-x*x)*(1-2*x)*(y-y*y)*(y-y*y);
    gradu(1) = 256.*2.*(x-x*x)*(x-x*x)*(y-y*y)*(1-2*y);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> gradu;
  }
  
  
  Tensor exact_2D_hessian(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                          <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,  <I><FONT COLOR="#B22222">// parameters, not needed
</FONT></I>                          <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>                          <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;) <I><FONT COLOR="#B22222">// unk_name, not needed
</FONT></I>  {
    Tensor hessu;
    
    <B><FONT COLOR="#228B22">const</FONT></B> Real x = p(0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real y = p(1);
  
    hessu(0,0) = 256.*2.*(1-6.*x+6.*x*x)*(y-y*y)*(y-y*y);
    hessu(0,1) = 256.*4.*(x-x*x)*(1.-2.*x)*(y-y*y)*(1.-2.*y);
    hessu(1,1) = 256.*2.*(x-x*x)*(x-x*x)*(1.-6.*y+6.*y*y);
  
    hessu(1,0) = hessu(0,1);
    <B><FONT COLOR="#A020F0">return</FONT></B> hessu;
  }
  
  
  
  Number forcing_function_2D(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p)
  {
    <B><FONT COLOR="#228B22">const</FONT></B> Real x = p(0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real y = p(1);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 256. * 8. * (3.*((y-y*y)*(y-y*y)+(x-x*x)*(x-x*x))
           + (1.-6.*x+6.*x*x)*(1.-6.*y+6.*y*y));
  }
  
  
  
  Number exact_3D_solution(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                           <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,  <I><FONT COLOR="#B22222">// parameters, not needed
</FONT></I>                           <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>                           <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;) <I><FONT COLOR="#B22222">// unk_name, not needed
</FONT></I>  {
    <B><FONT COLOR="#228B22">const</FONT></B> Real x = p(0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real y = p(1);
    <B><FONT COLOR="#228B22">const</FONT></B> Real z = p(2);
    
    <B><FONT COLOR="#A020F0">return</FONT></B> 4096.*(x-x*x)*(x-x*x)*(y-y*y)*(y-y*y)*(z-z*z)*(z-z*z);
  }
  
  
  Gradient exact_3D_derivative(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                               <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,  <I><FONT COLOR="#B22222">// parameters, not needed
</FONT></I>                               <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>                               <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;) <I><FONT COLOR="#B22222">// unk_name, not needed
</FONT></I>  {
    Gradient gradu;
    
    <B><FONT COLOR="#228B22">const</FONT></B> Real x = p(0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real y = p(1);
    <B><FONT COLOR="#228B22">const</FONT></B> Real z = p(2);
  
    gradu(0) = 4096.*2.*(x-x*x)*(1.-2.*x)*(y-y*y)*(y-y*y)*(z-z*z)*(z-z*z);
    gradu(1) = 4096.*2.*(x-x*x)*(x-x*x)*(y-y*y)*(1.-2.*y)*(z-z*z)*(z-z*z);
    gradu(2) = 4096.*2.*(x-x*x)*(x-x*x)*(y-y*y)*(y-y*y)*(z-z*z)*(1.-2.*z);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> gradu;
  }
  
  
  Tensor exact_3D_hessian(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                          <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp;,  <I><FONT COLOR="#B22222">// parameters, not needed
</FONT></I>                          <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>                          <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;) <I><FONT COLOR="#B22222">// unk_name, not needed
</FONT></I>  {
    Tensor hessu;
    
    <B><FONT COLOR="#228B22">const</FONT></B> Real x = p(0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real y = p(1);
    <B><FONT COLOR="#228B22">const</FONT></B> Real z = p(2);
  
    hessu(0,0) = 4096.*(2.-12.*x+12.*x*x)*(y-y*y)*(y-y*y)*(z-z*z)*(z-z*z);
    hessu(0,1) = 4096.*4.*(x-x*x)*(1.-2.*x)*(y-y*y)*(1.-2.*y)*(z-z*z)*(z-z*z);
    hessu(0,2) = 4096.*4.*(x-x*x)*(1.-2.*x)*(y-y*y)*(y-y*y)*(z-z*z)*(1.-2.*z);
    hessu(1,1) = 4096.*(x-x*x)*(x-x*x)*(2.-12.*y+12.*y*y)*(z-z*z)*(z-z*z);
    hessu(1,2) = 4096.*4.*(x-x*x)*(x-x*x)*(y-y*y)*(1.-2.*y)*(z-z*z)*(1.-2.*z);
    hessu(2,2) = 4096.*(x-x*x)*(x-x*x)*(y-y*y)*(y-y*y)*(2.-12.*z+12.*z*z);
  
    hessu(1,0) = hessu(0,1);
    hessu(2,0) = hessu(0,2);
    hessu(2,1) = hessu(1,2);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> hessu;
  }
  
  
  
  Number forcing_function_3D(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p)
  {
    <B><FONT COLOR="#228B22">const</FONT></B> Real x = p(0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real y = p(1);
    <B><FONT COLOR="#228B22">const</FONT></B> Real z = p(2);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 4096. * 8. * (3.*((y-y*y)*(y-y*y)*(x-x*x)*(x-x*x) +
                             (z-z*z)*(z-z*z)*(x-x*x)*(x-x*x) +
                             (z-z*z)*(z-z*z)*(y-y*y)*(y-y*y)) +
           (1.-6.*x+6.*x*x)*(1.-6.*y+6.*y*y)*(z-z*z)*(z-z*z) +
           (1.-6.*x+6.*x*x)*(1.-6.*z+6.*z*z)*(y-y*y)*(y-y*y) +
           (1.-6.*y+6.*y*y)*(1.-6.*z+6.*z*z)*(x-x*x)*(x-x*x));
  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_biharmonic(EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
  #ifdef LIBMESH_ENABLE_AMR
  #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Biharmonic&quot;</FONT></B>);
  
    PerfLog perf_log (<B><FONT COLOR="#BC8F8F">&quot;Matrix Assembly&quot;</FONT></B>,false);
    
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    LinearImplicitSystem&amp; system = es.get_system&lt;LinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Biharmonic&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
  
    FEType fe_type = dof_map.variable_type(0);
  
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    AutoPtr&lt;QBase&gt; qrule(fe_type.default_quadrature_rule(dim));
  
    fe-&gt;attach_quadrature_rule (qrule.get());
  
    AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
                
    AutoPtr&lt;QBase&gt; qface(fe_type.default_quadrature_rule(dim-1));
  
    fe_face-&gt;attach_quadrature_rule (qface.get());
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; q_point = fe-&gt;get_xyz();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealTensor&gt; &gt;&amp; d2phi = fe-&gt;get_d2phi();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Real&gt; shape_laplacian;
  
    DenseMatrix&lt;Number&gt; Ke;
    DenseVector&lt;Number&gt; Fe;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
  
  
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
  
        shape_laplacian.resize(dof_indices.size());
  
        perf_log.pop(<B><FONT COLOR="#BC8F8F">&quot;elem init&quot;</FONT></B>);      
  
        perf_log.push (<B><FONT COLOR="#BC8F8F">&quot;Ke&quot;</FONT></B>);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule-&gt;n_points(); qp++)
          {
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
              {
                shape_laplacian[i] = d2phi[i][qp](0,0)+d2phi[i][qp](1,1);
                <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 3)
                   shape_laplacian[i] += d2phi[i][qp](2,2);
              }
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
                Ke(i,j) += JxW[qp]*
                           shape_laplacian[i]*shape_laplacian[j];
          }
  
        perf_log.pop (<B><FONT COLOR="#BC8F8F">&quot;Ke&quot;</FONT></B>);
  
  
        {
          perf_log.push (<B><FONT COLOR="#BC8F8F">&quot;BCs&quot;</FONT></B>);
  
          <B><FONT COLOR="#228B22">const</FONT></B> Real penalty = 1e10;
          <B><FONT COLOR="#228B22">const</FONT></B> Real penalty2 = 1e10;
  
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> s=0; s&lt;elem-&gt;n_sides(); s++)
            <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(s) == NULL)
              {
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp;  phi_face =
                                fe_face-&gt;get_phi();
  
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi_face =
                                fe_face-&gt;get_dphi();
  
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
                                                                                 
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; qface_point = fe_face-&gt;get_xyz();
  
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; face_normals =
                                fe_face-&gt;get_normals();
  
                fe_face-&gt;reinit(elem, s);
                                                                                  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qface-&gt;n_points(); qp++)
                  {
                    Number value = exact_solution(qface_point[qp],
                                                  es.parameters, <B><FONT COLOR="#BC8F8F">&quot;null&quot;</FONT></B>,
                                                  <B><FONT COLOR="#BC8F8F">&quot;void&quot;</FONT></B>);
                    Gradient flux = exact_2D_derivative(qface_point[qp],
                                                        es.parameters,
                                                        <B><FONT COLOR="#BC8F8F">&quot;null&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;void&quot;</FONT></B>);
  
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi_face.size(); i++)
                      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi_face.size(); j++)
                        Ke(i,j) += JxW_face[qp] *
                                   (penalty * phi_face[i][qp] *
                                    phi_face[j][qp] + penalty2
                                    * (dphi_face[i][qp] *
                                    face_normals[qp]) *
                                    (dphi_face[j][qp] *
                                     face_normals[qp]));
  
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi_face.size(); i++)
                      Fe(i) += JxW_face[qp] *
                                      (penalty * value * phi_face[i][qp]
                                       + penalty2 * 
                                       (flux * face_normals[qp])
                                      * (dphi_face[i][qp]
                                         * face_normals[qp]));
  
                  }
              } 
          
          perf_log.pop (<B><FONT COLOR="#BC8F8F">&quot;BCs&quot;</FONT></B>);
        } 
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule-&gt;n_points(); qp++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
            Fe(i) += JxW[qp]*phi[i][qp]*forcing_function(q_point[qp]);
  
        perf_log.push (<B><FONT COLOR="#BC8F8F">&quot;matrix insertion&quot;</FONT></B>);
  
        dof_map.constrain_element_matrix_and_vector(Ke, Fe, dof_indices);
        system.matrix-&gt;add_matrix (Ke, dof_indices);
        system.rhs-&gt;add_vector    (Fe, dof_indices);
  
        perf_log.pop (<B><FONT COLOR="#BC8F8F">&quot;matrix insertion&quot;</FONT></B>);
      }
  
  
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
</FONT></I>  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_ENABLE_AMR
</FONT></I>  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example  mpirun -np 2 ./ex15-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
 EquationSystems
  n_systems()=1
   System "Biharmonic"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="HERMITE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="THIRD", "THIRD" 
    n_dofs()=36
    n_local_dofs()=24
    n_constrained_dofs()=0
    n_vectors()=1

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=25
    n_local_nodes()=15
  n_elem()=4
    n_local_elem()=2
    n_active_elem()=4
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Biharmonic"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="HERMITE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="THIRD", "THIRD" 
    n_dofs()=36
    n_local_dofs()=24
    n_constrained_dofs()=0
    n_vectors()=1

Beginning Solve 0
Linear solver converged at step: 9, final residual: 2.42815e-11
L2-Norm is: 0.384025
H1-Norm is: 1.98976
H2-Norm is: 14.3417

L2-Error is: 0.0335358
H1-Error is: 0.267039
H2-Error is: 3.51162

  Refining the mesh...
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=81
    n_local_nodes()=45
  n_elem()=20
    n_local_elem()=10
    n_active_elem()=16
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Biharmonic"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="HERMITE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="THIRD", "THIRD" 
    n_dofs()=100
    n_local_dofs()=60
    n_constrained_dofs()=0
    n_vectors()=1

Beginning Solve 1
Linear solver converged at step: 15, final residual: 2.80063e-11
L2-Norm is: 0.404988
H1-Norm is: 2.02995
H2-Norm is: 14.7459

L2-Error is: 0.0020746
H1-Error is: 0.0316727
H2-Error is: 0.822125

  Refining the mesh...
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=289
    n_local_nodes()=153
  n_elem()=84
    n_local_elem()=42
    n_active_elem()=64
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Biharmonic"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="HERMITE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="THIRD", "THIRD" 
    n_dofs()=324
    n_local_dofs()=180
    n_constrained_dofs()=0
    n_vectors()=1

Beginning Solve 2
Linear solver converged at step: 21, final residual: 2.21056e-11
L2-Norm is: 0.406264
H1-Norm is: 2.03164
H2-Norm is: 14.7676

L2-Error is: 0.000129445
H1-Error is: 0.00390589
H2-Error is: 0.202531

  Refining the mesh...
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1089
    n_local_nodes()=561
  n_elem()=340
    n_local_elem()=170
    n_active_elem()=256
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Biharmonic"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="HERMITE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="THIRD", "THIRD" 
    n_dofs()=1156
    n_local_dofs()=612
    n_constrained_dofs()=0
    n_vectors()=1

Beginning Solve 3
Linear solver converged at step: 89, final residual: 1.96136e-11
L2-Norm is: 0.406344
H1-Norm is: 2.03174
H2-Norm is: 14.7689

L2-Error is: 8.07721e-06
H1-Error is: 0.000486566
H2-Error is: 0.050454

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex15-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Thu Feb  3 12:11:08 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           4.247e-01      1.00254   4.242e-01
Objects:              2.490e+02      1.00000   2.490e+02
Flops:                1.899e+07      1.18962   1.748e+07  3.496e+07
Flops/sec:            4.472e+07      1.18661   4.120e+07  8.240e+07
MPI Messages:         2.250e+02      1.00000   2.250e+02  4.500e+02
MPI Message Lengths:  1.187e+05      1.02261   5.215e+02  2.347e+05
MPI Reductions:       5.610e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 4.2414e-01 100.0%  3.4956e+07 100.0%  4.500e+02 100.0%  5.215e+02      100.0%  4.880e+02  87.0% 

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

VecMDot              134 1.0 1.7250e-03 1.9 1.77e+06 1.1 0.0e+00 0.0e+00 1.3e+02  0 10  0  0 24   0 10  0  0 27  1930
VecNorm              144 1.0 7.3934e-04 1.1 1.25e+05 1.1 0.0e+00 0.0e+00 1.4e+02  0  1  0  0 26   0  1  0  0 30   317
VecScale             140 1.0 6.6757e-05 1.0 6.15e+04 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1729
VecCopy               21 1.0 8.5831e-06 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               168 1.0 5.2214e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               12 1.0 3.4089e-03 1.0 8.40e+03 1.1 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     5
VecMAXPY             140 1.0 4.8304e-04 1.1 1.89e+06 1.1 0.0e+00 0.0e+00 0.0e+00  0 10  0  0  0   0 10  0  0  0  7361
VecAssemblyBegin      30 1.0 3.9012e-0310.0 0.00e+00 0.0 1.4e+01 3.1e+02 7.2e+01  1  0  3  2 13   1  0  3  2 15     0
VecAssemblyEnd        30 1.0 2.4080e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      163 1.0 2.4390e-04 1.1 0.00e+00 0.0 3.0e+02 4.2e+02 0.0e+00  0  0 68 55  0   0  0 68 55  0     0
VecScatterEnd        163 1.0 1.6778e-03 8.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize         140 1.0 8.7976e-04 1.1 1.84e+05 1.1 0.0e+00 0.0e+00 1.4e+02  0  1  0  0 25   0  1  0  0 29   394
MatMult              140 1.0 3.4900e-03 1.6 4.00e+06 1.1 2.8e+02 4.3e+02 0.0e+00  1 21 62 51  0   1 21 62 51  0  2148
MatSolve             140 1.0 3.8195e-03 1.2 8.99e+06 1.2 0.0e+00 0.0e+00 0.0e+00  1 47  0  0  0   1 47  0  0  0  4305
MatLUFactorNum         4 1.0 1.6959e-03 1.3 2.15e+06 1.3 0.0e+00 0.0e+00 0.0e+00  0 11  0  0  0   0 11  0  0  0  2223
MatILUFactorSym        4 1.0 4.9901e-03 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 4.0e+00  1  0  0  0  1   1  0  0  0  1     0
MatAssemblyBegin       8 1.0 4.2415e-04 1.8 0.00e+00 0.0 1.2e+01 5.1e+03 1.6e+01  0  0  3 26  3   0  0  3 26  3     0
MatAssemblyEnd         8 1.0 7.1096e-04 1.1 0.00e+00 0.0 1.6e+01 7.0e+01 3.2e+01  0  0  4  0  6   0  0  4  0  7     0
MatGetRowIJ            4 1.0 9.5367e-07 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         4 1.0 7.4863e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.6e+01  0  0  0  0  3   0  0  0  0  3     0
MatZeroEntries        12 1.0 8.1301e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       134 1.0 2.2163e-03 1.6 3.54e+06 1.1 0.0e+00 0.0e+00 1.3e+02  0 19  0  0 24   0 19  0  0 27  3006
KSPSetup               8 1.0 6.8188e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               4 1.0 2.0028e-02 1.0 1.90e+07 1.2 2.8e+02 4.3e+02 3.0e+02  5100 62 51 53   5100 62 51 61  1745
PCSetUp                8 1.0 7.1108e-03 1.3 2.15e+06 1.3 0.0e+00 0.0e+00 2.0e+01  1 11  0  0  4   1 11  0  0  4   530
PCSetUpOnBlocks        4 1.0 6.8619e-03 1.3 2.15e+06 1.3 0.0e+00 0.0e+00 2.0e+01  1 11  0  0  4   1 11  0  0  4   549
PCApply              140 1.0 4.5435e-03 1.2 8.99e+06 1.2 0.0e+00 0.0e+00 0.0e+00  1 47  0  0  0   1 47  0  0  0  3619
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec   143            143       456576     0
         Vec Scatter    21             21        18228     0
           Index Set    49             49        38800     0
   IS L to G Mapping     4              4         5664     0
              Matrix    16             16      1098576     0
       Krylov Solver     8              8        75520     0
      Preconditioner     8              8         5632     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 1.23978e-06
Average time for zero size MPI_Send(): 5.48363e-06
#PETSc Option Table entries:
-ksp_right_pc
-log_summary
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
| Time:           Thu Feb  3 12:11:08 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.43281, Active time=0.391137                                              |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 4         0.0002      0.000040    0.0002      0.000048    0.04     0.05     |
|   compute_sparsity()           4         0.0015      0.000376    0.0017      0.000424    0.38     0.43     |
|   create_dof_constraints()     4         0.0004      0.000103    0.0004      0.000103    0.11     0.11     |
|   distribute_dofs()            4         0.0006      0.000151    0.0015      0.000371    0.15     0.38     |
|   dof_indices()                1206      0.0007      0.000001    0.0007      0.000001    0.19     0.19     |
|   old_dof_indices()            336       0.0002      0.000001    0.0002      0.000001    0.05     0.05     |
|   prepare_send_list()          4         0.0000      0.000004    0.0000      0.000004    0.00     0.00     |
|   reinit()                     4         0.0008      0.000188    0.0008      0.000188    0.19     0.19     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         906       0.0014      0.000002    0.0014      0.000002    0.36     0.36     |
|   compute_face_map()           60        0.0003      0.000005    0.0007      0.000012    0.07     0.18     |
|   compute_shape_functions()    408       0.0007      0.000002    0.0007      0.000002    0.18     0.18     |
|   init_face_shape_functions()  60        0.0001      0.000001    0.0001      0.000001    0.01     0.01     |
|   init_shape_functions()       408       0.3343      0.000819    0.3343      0.000819    85.48    85.48    |
|   inverse_map()                4320      0.0066      0.000002    0.0066      0.000002    1.69     1.69     |
|                                                                                                            |
| LocationMap                                                                                                |
|   find()                       1680      0.0009      0.000001    0.0009      0.000001    0.23     0.23     |
|   init()                       6         0.0003      0.000057    0.0003      0.000057    0.09     0.09     |
|                                                                                                            |
| Mesh                                                                                                       |
|   all_second_order()           1         0.0001      0.000060    0.0001      0.000060    0.02     0.02     |
|   contract()                   3         0.0000      0.000009    0.0001      0.000021    0.01     0.02     |
|   find_neighbors()             5         0.0007      0.000149    0.0008      0.000160    0.19     0.20     |
|   renumber_nodes_and_elem()    13        0.0001      0.000008    0.0001      0.000008    0.03     0.03     |
|                                                                                                            |
| MeshCommunication                                                                                          |
|   compute_hilbert_indices()    6         0.0006      0.000107    0.0006      0.000107    0.16     0.16     |
|   find_global_indices()        6         0.0002      0.000033    0.0014      0.000226    0.05     0.35     |
|   parallel_sort()              6         0.0003      0.000050    0.0004      0.000061    0.08     0.09     |
|                                                                                                            |
| MeshRefinement                                                                                             |
|   _coarsen_elements()          3         0.0000      0.000008    0.0000      0.000009    0.01     0.01     |
|   _refine_elements()           6         0.0018      0.000294    0.0041      0.000684    0.45     1.05     |
|   add_point()                  1680      0.0012      0.000001    0.0022      0.000001    0.31     0.57     |
|   make_coarsening_compatible() 3         0.0001      0.000049    0.0001      0.000049    0.04     0.04     |
|   make_refinement_compatible() 3         0.0000      0.000001    0.0000      0.000003    0.00     0.00     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0001      0.000060    0.0001      0.000060    0.02     0.02     |
|                                                                                                            |
| MetisPartitioner                                                                                           |
|   partition()                  5         0.0008      0.000158    0.0018      0.000358    0.20     0.46     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  22        0.0001      0.000003    0.0001      0.000003    0.02     0.02     |
|   broadcast()                  4         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   gather()                     4         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   max(bool)                    16        0.0000      0.000003    0.0000      0.000003    0.01     0.01     |
|   max(scalar)                  20        0.0008      0.000038    0.0008      0.000038    0.19     0.19     |
|   max(vector)                  6         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   min(bool)                    6         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   min(vector)                  6         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   probe()                      38        0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|   receive()                    38        0.0001      0.000002    0.0001      0.000003    0.02     0.03     |
|   send()                       38        0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|   send_receive()               50        0.0001      0.000002    0.0003      0.000006    0.02     0.07     |
|   sum()                        29        0.0006      0.000022    0.0006      0.000022    0.16     0.16     |
|   wait()                       38        0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|                                                                                                            |
| Partitioner                                                                                                |
|   set_node_processor_ids()     5         0.0003      0.000063    0.0004      0.000076    0.08     0.10     |
|   set_parent_processor_ids()   5         0.0000      0.000008    0.0000      0.000008    0.01     0.01     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      4         0.0217      0.005429    0.0217      0.005429    5.55     5.55     |
|                                                                                                            |
| ProjectVector                                                                                              |
|   operator()                   3         0.0043      0.001449    0.3111      0.103693    1.11     79.53    |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   4         0.0036      0.000912    0.0244      0.006091    0.93     6.23     |
|   project_vector()             3         0.0041      0.001381    0.3153      0.105106    1.06     80.62    |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        11494     0.3911                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./ex15-opt
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
