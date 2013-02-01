<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("systems_of_equations_ex3",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file systems_of_equations_ex3.C with comments: </h1> 
<div class = "comment">
<h1>Systems Example 3 - Navier-Stokes with SCALAR Lagrange Multiplier</h1>

<br><br>This example shows how the transient Navier-Stokes problem from
example 13 can be solved using a scalar Lagrange multiplier formulation to
constrain the integral of the pressure variable, rather than pinning the
pressure at a single point.


<br><br>C++ include files that we need
</div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
        #include &lt;algorithm&gt;
        #include &lt;sstream&gt;
        #include &lt;math.h&gt;
        
</pre>
</div>
<div class = "comment">
Basic include file needed for the mesh functionality.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh.h"
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/exodusII_io.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/fe.h"
        #include "libmesh/quadrature_gauss.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/dense_vector.h"
        #include "libmesh/linear_implicit_system.h"
        #include "libmesh/transient_system.h"
        #include "libmesh/perf_log.h"
        #include "libmesh/boundary_info.h"
        #include "libmesh/utility.h"
        
</pre>
</div>
<div class = "comment">
For systems of equations the \p DenseSubMatrix
and \p DenseSubVector provide convenient ways for
assembling the element matrix and vector on a
component-by-component basis.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/dense_submatrix.h"
        #include "libmesh/dense_subvector.h"
        
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
matrix and right-hand-side.
</div>

<div class ="fragment">
<pre>
        void assemble_stokes (EquationSystems& es,
                              const std::string& system_name);
        
</pre>
</div>
<div class = "comment">
The main program.
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
Trilinos solver NaNs by default on the zero pressure block.
We'll skip this example for now.
</div>

<div class ="fragment">
<pre>
          if (libMesh::default_solver_package() == TRILINOS_SOLVERS)
            {
              std::cout &lt;&lt; "We skip example 13 when using the Trilinos solvers.\n"
                        &lt;&lt; std::endl;
              return 0;
            }
        
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
Use the MeshTools::Generation mesh generator to create a uniform
2D grid on the square [-1,1]^2.  We instruct the mesh generator
to build a mesh of 8x8 \p Quad9 elements in 2D.  Building these
higher-order elements allows us to use higher-order
approximation, as in example 3.
</div>

<div class ="fragment">
<pre>
          MeshTools::Generation::build_square (mesh,
                                               20, 20,
                                               0., 1.,
                                               0., 1.,
                                               QUAD9);
          
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
Create an equation systems object.
</div>

<div class ="fragment">
<pre>
          EquationSystems equation_systems (mesh);
          
</pre>
</div>
<div class = "comment">
Declare the system and its variables.
Creates a transient system named "Navier-Stokes"
</div>

<div class ="fragment">
<pre>
          TransientLinearImplicitSystem & system = 
            equation_systems.add_system&lt;TransientLinearImplicitSystem&gt; ("Navier-Stokes");
          
</pre>
</div>
<div class = "comment">
Add the variables "u" & "v" to "Navier-Stokes".  They
will be approximated using second-order approximation.
</div>

<div class ="fragment">
<pre>
          system.add_variable ("u", SECOND);
          system.add_variable ("v", SECOND);
        
</pre>
</div>
<div class = "comment">
Add the variable "p" to "Navier-Stokes". This will
be approximated with a first-order basis,
providing an LBB-stable pressure-velocity pair.
</div>

<div class ="fragment">
<pre>
          system.add_variable ("p", FIRST);
        
</pre>
</div>
<div class = "comment">
Add a scalar Lagrange multiplier to constrain the
pressure to have zero mean.
</div>

<div class ="fragment">
<pre>
          system.add_variable ("alpha", FIRST, SCALAR);
        
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly
function.
</div>

<div class ="fragment">
<pre>
          system.attach_assemble_function (assemble_stokes);
          
</pre>
</div>
<div class = "comment">
Initialize the data structures for the equation system.
</div>

<div class ="fragment">
<pre>
          equation_systems.init ();
        
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
Create a performance-logging object for this example
</div>

<div class ="fragment">
<pre>
          PerfLog perf_log("Systems Example 3");
          
</pre>
</div>
<div class = "comment">
Get a reference to the Stokes system to use later.
</div>

<div class ="fragment">
<pre>
          TransientLinearImplicitSystem&  navier_stokes_system =
                equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;("Navier-Stokes");
        
</pre>
</div>
<div class = "comment">
Now we begin the timestep loop to compute the time-accurate
solution of the equations.
</div>

<div class ="fragment">
<pre>
          const Real dt = 0.01;
          navier_stokes_system.time     = 0.0;
          const unsigned int n_timesteps = 15;
        
</pre>
</div>
<div class = "comment">
The number of steps and the stopping criterion are also required
for the nonlinear iterations.
</div>

<div class ="fragment">
<pre>
          const unsigned int n_nonlinear_steps = 15;
          const Real nonlinear_tolerance       = 1.e-3;
        
</pre>
</div>
<div class = "comment">
We also set a standard linear solver flag in the EquationSystems object
which controls the maxiumum number of linear solver iterations allowed.
</div>

<div class ="fragment">
<pre>
          equation_systems.parameters.set&lt;unsigned int&gt;("linear solver maximum iterations") = 250;
          
</pre>
</div>
<div class = "comment">
Tell the system of equations what the timestep is by using
the set_parameter function.  The matrix assembly routine can
then reference this parameter.
</div>

<div class ="fragment">
<pre>
          equation_systems.parameters.set&lt;Real&gt; ("dt")   = dt;
        
</pre>
</div>
<div class = "comment">
The first thing to do is to get a copy of the solution at
the current nonlinear iteration.  This value will be used to
determine if we can exit the nonlinear loop.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;NumericVector&lt;Number&gt; &gt;
            last_nonlinear_soln (navier_stokes_system.solution-&gt;clone());
        
          for (unsigned int t_step=0; t_step&lt;n_timesteps; ++t_step)
            {
</pre>
</div>
<div class = "comment">
Incremenet the time counter, set the time and the
time step size as parameters in the EquationSystem.
</div>

<div class ="fragment">
<pre>
              navier_stokes_system.time += dt;
        
</pre>
</div>
<div class = "comment">
A pretty update message
</div>

<div class ="fragment">
<pre>
              std::cout &lt;&lt; "\n\n*** Solving time step " &lt;&lt; t_step &lt;&lt; 
                           ", time = " &lt;&lt; navier_stokes_system.time &lt;&lt;
                           " ***" &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Now we need to update the solution vector from the
previous time step.  This is done directly through
the reference to the Stokes system.
</div>

<div class ="fragment">
<pre>
              *navier_stokes_system.old_local_solution = *navier_stokes_system.current_local_solution;
        
</pre>
</div>
<div class = "comment">
At the beginning of each solve, reset the linear solver tolerance
to a "reasonable" starting value.
</div>

<div class ="fragment">
<pre>
              const Real initial_linear_solver_tol = 1.e-6;
              equation_systems.parameters.set&lt;Real&gt; ("linear solver tolerance") = initial_linear_solver_tol;
        
</pre>
</div>
<div class = "comment">
Now we begin the nonlinear loop
</div>

<div class ="fragment">
<pre>
              for (unsigned int l=0; l&lt;n_nonlinear_steps; ++l)
                {
</pre>
</div>
<div class = "comment">
Update the nonlinear solution.
</div>

<div class ="fragment">
<pre>
                  last_nonlinear_soln-&gt;zero();
                  last_nonlinear_soln-&gt;add(*navier_stokes_system.solution);
                  
</pre>
</div>
<div class = "comment">
Assemble & solve the linear system.
</div>

<div class ="fragment">
<pre>
                  perf_log.push("linear solve");
                  equation_systems.get_system("Navier-Stokes").solve();
                  perf_log.pop("linear solve");
        
</pre>
</div>
<div class = "comment">
Compute the difference between this solution and the last
nonlinear iterate.
</div>

<div class ="fragment">
<pre>
                  last_nonlinear_soln-&gt;add (-1., *navier_stokes_system.solution);
        
</pre>
</div>
<div class = "comment">
Close the vector before computing its norm
</div>

<div class ="fragment">
<pre>
                  last_nonlinear_soln-&gt;close();
        
</pre>
</div>
<div class = "comment">
Compute the l2 norm of the difference
</div>

<div class ="fragment">
<pre>
                  const Real norm_delta = last_nonlinear_soln-&gt;l2_norm();
        
</pre>
</div>
<div class = "comment">
How many iterations were required to solve the linear system?
</div>

<div class ="fragment">
<pre>
                  const unsigned int n_linear_iterations = navier_stokes_system.n_linear_iterations();
                  
</pre>
</div>
<div class = "comment">
What was the final residual of the linear system?
</div>

<div class ="fragment">
<pre>
                  const Real final_linear_residual = navier_stokes_system.final_linear_residual();
                  
</pre>
</div>
<div class = "comment">
Print out convergence information for the linear and
nonlinear iterations.
</div>

<div class ="fragment">
<pre>
                  std::cout &lt;&lt; "Linear solver converged at step: "
                            &lt;&lt; n_linear_iterations
                            &lt;&lt; ", final residual: "
                            &lt;&lt; final_linear_residual
                            &lt;&lt; "  Nonlinear convergence: ||u - u_old|| = "
                            &lt;&lt; norm_delta
                            &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Terminate the solution iteration if the difference between
this nonlinear iterate and the last is sufficiently small, AND
if the most recent linear system was solved to a sufficient tolerance.
</div>

<div class ="fragment">
<pre>
                  if ((norm_delta &lt; nonlinear_tolerance) &&
                      (navier_stokes_system.final_linear_residual() &lt; nonlinear_tolerance))
                    {
                      std::cout &lt;&lt; " Nonlinear solver converged at step "
                                &lt;&lt; l
                                &lt;&lt; std::endl;
                      break;
                    }
                  
</pre>
</div>
<div class = "comment">
Otherwise, decrease the linear system tolerance.  For the inexact Newton
method, the linear solver tolerance needs to decrease as we get closer to
the solution to ensure quadratic convergence.  The new linear solver tolerance
is chosen (heuristically) as the square of the previous linear system residual norm.
Real flr2 = final_linear_residual*final_linear_residual;
</div>

<div class ="fragment">
<pre>
                  equation_systems.parameters.set&lt;Real&gt; ("linear solver tolerance") =
                    std::min(Utility::pow&lt;2&gt;(final_linear_residual), initial_linear_solver_tol);
        
                } // end nonlinear loop
              
</pre>
</div>
<div class = "comment">
Write out every nth timestep to file.
</div>

<div class ="fragment">
<pre>
              const unsigned int write_interval = 1;
              
        #ifdef LIBMESH_HAVE_EXODUS_API
              if ((t_step+1)%write_interval == 0)
                {
                  std::ostringstream file_name;
        
</pre>
</div>
<div class = "comment">
We write the file in the ExodusII format.
</div>

<div class ="fragment">
<pre>
                  file_name &lt;&lt; "out_"
                            &lt;&lt; std::setw(3)
                            &lt;&lt; std::setfill('0')
                            &lt;&lt; std::right
                            &lt;&lt; t_step + 1
                            &lt;&lt; ".e";
                  
                  ExodusII_IO(mesh).write_equation_systems (file_name.str(),
                                                      equation_systems);
                }
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
            } // end timestep loop.
          
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
The matrix assembly function to be called at each time step to
prepare for the linear solve.
</div>

<div class ="fragment">
<pre>
        void assemble_stokes (EquationSystems& es,
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
          libmesh_assert_equal_to (system_name, "Navier-Stokes");
          
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
Get a reference to the Stokes system object.
</div>

<div class ="fragment">
<pre>
          TransientLinearImplicitSystem & navier_stokes_system =
            es.get_system&lt;TransientLinearImplicitSystem&gt; ("Navier-Stokes");
        
</pre>
</div>
<div class = "comment">
Numeric ids corresponding to each variable in the system
</div>

<div class ="fragment">
<pre>
          const unsigned int u_var = navier_stokes_system.variable_number ("u");
          const unsigned int v_var = navier_stokes_system.variable_number ("v");
          const unsigned int p_var = navier_stokes_system.variable_number ("p");
          const unsigned int alpha_var = navier_stokes_system.variable_number ("alpha");
          
</pre>
</div>
<div class = "comment">
Get the Finite Element type for "u".  Note this will be
the same as the type for "v".
</div>

<div class ="fragment">
<pre>
          FEType fe_vel_type = navier_stokes_system.variable_type(u_var);
          
</pre>
</div>
<div class = "comment">
Get the Finite Element type for "p".
</div>

<div class ="fragment">
<pre>
          FEType fe_pres_type = navier_stokes_system.variable_type(p_var);
        
</pre>
</div>
<div class = "comment">
Build a Finite Element object of the specified type for
the velocity variables.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;FEBase&gt; fe_vel  (FEBase::build(dim, fe_vel_type));
            
</pre>
</div>
<div class = "comment">
Build a Finite Element object of the specified type for
the pressure variables.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;FEBase&gt; fe_pres (FEBase::build(dim, fe_pres_type));
          
</pre>
</div>
<div class = "comment">
A Gauss quadrature rule for numerical integration.
Let the \p FEType object decide what order rule is appropriate.
</div>

<div class ="fragment">
<pre>
          QGauss qrule (dim, fe_vel_type.default_quadrature_order());
        
</pre>
</div>
<div class = "comment">
Tell the finite element objects to use our quadrature rule.
</div>

<div class ="fragment">
<pre>
          fe_vel-&gt;attach_quadrature_rule (&qrule);
          fe_pres-&gt;attach_quadrature_rule (&qrule);
          
</pre>
</div>
<div class = "comment">
Here we define some references to cell-specific data that
will be used to assemble the linear system.

<br><br>The element Jacobian * quadrature weight at each integration point.   
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt;& JxW = fe_vel-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The element shape functions evaluated at the quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi = fe_vel-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The element shape function gradients for the velocity
variables evaluated at the quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi = fe_vel-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
The element shape functions for the pressure variable
evaluated at the quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& psi = fe_pres-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The value of the linear shape function gradients at the quadrature points
const std::vector<std::vector<RealGradient> >& dpsi = fe_pres->get_dphi();
  

<br><br>A reference to the \p DofMap object for this system.  The \p DofMap
object handles the index translation from node and element numbers
to degree of freedom numbers.  We will talk more about the \p DofMap
in future examples.
</div>

<div class ="fragment">
<pre>
          const DofMap & dof_map = navier_stokes_system.get_dof_map();
        
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
        
          DenseSubMatrix&lt;Number&gt;
            Kuu(Ke), Kuv(Ke), Kup(Ke),
            Kvu(Ke), Kvv(Ke), Kvp(Ke),
            Kpu(Ke), Kpv(Ke), Kpp(Ke);
          DenseSubMatrix&lt;Number&gt; Kalpha_p(Ke), Kp_alpha(Ke);
        
          DenseSubVector&lt;Number&gt;
            Fu(Fe),
            Fv(Fe),
            Fp(Fe);
        
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
          std::vector&lt;dof_id_type&gt; dof_indices_u;
          std::vector&lt;dof_id_type&gt; dof_indices_v;
          std::vector&lt;dof_id_type&gt; dof_indices_p;
          std::vector&lt;dof_id_type&gt; dof_indices_alpha;
        
</pre>
</div>
<div class = "comment">
Find out what the timestep size parameter is from the system, and
the value of theta for the theta method.  We use implicit Euler (theta=1)
for this simulation even though it is only first-order accurate in time.
The reason for this decision is that the second-order Crank-Nicolson
method is notoriously oscillatory for problems with discontinuous
initial data such as the lid-driven cavity.  Therefore,
we sacrifice accuracy in time for stability, but since the solution
reaches steady state relatively quickly we can afford to take small
timesteps.  If you monitor the initial nonlinear residual for this
simulation, you should see that it is monotonically decreasing in time.
</div>

<div class ="fragment">
<pre>
          const Real dt    = es.parameters.get&lt;Real&gt;("dt");
</pre>
</div>
<div class = "comment">
const Real time  = es.parameters.get<Real>("time");
</div>

<div class ="fragment">
<pre>
          const Real theta = 1.;
            
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
              dof_map.dof_indices (elem, dof_indices_u, u_var);
              dof_map.dof_indices (elem, dof_indices_v, v_var);
              dof_map.dof_indices (elem, dof_indices_p, p_var);
              dof_map.dof_indices (elem, dof_indices_alpha, alpha_var);
        
              const unsigned int n_dofs   = dof_indices.size();
              const unsigned int n_u_dofs = dof_indices_u.size(); 
              const unsigned int n_v_dofs = dof_indices_v.size(); 
              const unsigned int n_p_dofs = dof_indices_p.size();
              
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
              fe_vel-&gt;reinit  (elem);
              fe_pres-&gt;reinit (elem);
        
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
              Ke.resize (n_dofs, n_dofs);
              Fe.resize (n_dofs);
        
</pre>
</div>
<div class = "comment">
Reposition the submatrices...  The idea is this:

<br><br>-           -          -  -
| Kuu Kuv Kup |        | Fu |
Ke = | Kvu Kvv Kvp |;  Fe = | Fv |
| Kpu Kpv Kpp |        | Fp |
-           -          -  -

<br><br>The \p DenseSubMatrix.repostition () member takes the
(row_offset, column_offset, row_size, column_size).

<br><br>Similarly, the \p DenseSubVector.reposition () member
takes the (row_offset, row_size)
</div>

<div class ="fragment">
<pre>
              Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
              Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
              Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);
              
              Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
              Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
              Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);
              
              Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
              Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
              Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);
        
</pre>
</div>
<div class = "comment">
Also, add a row and a column to constrain the pressure
</div>

<div class ="fragment">
<pre>
              Kp_alpha.reposition (p_var*n_u_dofs, p_var*n_u_dofs+n_p_dofs, n_p_dofs, 1);
              Kalpha_p.reposition (p_var*n_u_dofs+n_p_dofs, p_var*n_u_dofs, 1, n_p_dofs);
        
        
              Fu.reposition (u_var*n_u_dofs, n_u_dofs);
              Fv.reposition (v_var*n_u_dofs, n_v_dofs);
              Fp.reposition (p_var*n_u_dofs, n_p_dofs);
        
</pre>
</div>
<div class = "comment">
Now we will build the element matrix and right-hand-side.
Constructing the RHS requires the solution and its
gradient from the previous timestep.  This must be
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
Values to hold the solution & its gradient at the previous timestep.
</div>

<div class ="fragment">
<pre>
                  Number   u = 0., u_old = 0.;
                  Number   v = 0., v_old = 0.;
                  Number   p_old = 0.;
                  Gradient grad_u, grad_u_old;
                  Gradient grad_v, grad_v_old;
                  
</pre>
</div>
<div class = "comment">
Compute the velocity & its gradient from the previous timestep
and the old Newton iterate.
</div>

<div class ="fragment">
<pre>
                  for (unsigned int l=0; l&lt;n_u_dofs; l++)
                    {
</pre>
</div>
<div class = "comment">
From the old timestep:
</div>

<div class ="fragment">
<pre>
                      u_old += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_u[l]);
                      v_old += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_v[l]);
                      grad_u_old.add_scaled (dphi[l][qp],navier_stokes_system.old_solution (dof_indices_u[l]));
                      grad_v_old.add_scaled (dphi[l][qp],navier_stokes_system.old_solution (dof_indices_v[l]));
        
</pre>
</div>
<div class = "comment">
From the previous Newton iterate:
</div>

<div class ="fragment">
<pre>
                      u += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_u[l]); 
                      v += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_v[l]);
                      grad_u.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_u[l]));
                      grad_v.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_v[l]));
                    }
        
</pre>
</div>
<div class = "comment">
Compute the old pressure value at this quadrature point.
</div>

<div class ="fragment">
<pre>
                  for (unsigned int l=0; l&lt;n_p_dofs; l++)
                    {
                      p_old += psi[l][qp]*navier_stokes_system.old_solution (dof_indices_p[l]);
                    }
        
</pre>
</div>
<div class = "comment">
Definitions for convenience.  It is sometimes simpler to do a
dot product if you have the full vector at your disposal.
</div>

<div class ="fragment">
<pre>
                  const NumberVectorValue U_old (u_old, v_old);
                  const NumberVectorValue U     (u,     v);
                  const Number  u_x = grad_u(0);
                  const Number  u_y = grad_u(1);
                  const Number  v_x = grad_v(0);
                  const Number  v_y = grad_v(1);
                  
</pre>
</div>
<div class = "comment">
First, an i-loop over the velocity degrees of freedom.
We know that n_u_dofs == n_v_dofs so we can compute contributions
for both at the same time.
</div>

<div class ="fragment">
<pre>
                  for (unsigned int i=0; i&lt;n_u_dofs; i++)
                    {
                      Fu(i) += JxW[qp]*(u_old*phi[i][qp] -                            // mass-matrix term 
                                        (1.-theta)*dt*(U_old*grad_u_old)*phi[i][qp] + // convection term
                                        (1.-theta)*dt*p_old*dphi[i][qp](0)  -         // pressure term on rhs
                                        (1.-theta)*dt*(grad_u_old*dphi[i][qp]) +      // diffusion term on rhs
                                        theta*dt*(U*grad_u)*phi[i][qp]);              // Newton term
        
                        
                      Fv(i) += JxW[qp]*(v_old*phi[i][qp] -                             // mass-matrix term
                                        (1.-theta)*dt*(U_old*grad_v_old)*phi[i][qp] +  // convection term
                                        (1.-theta)*dt*p_old*dphi[i][qp](1) -           // pressure term on rhs
                                        (1.-theta)*dt*(grad_v_old*dphi[i][qp]) +       // diffusion term on rhs
                                        theta*dt*(U*grad_v)*phi[i][qp]);               // Newton term
                    
        
</pre>
</div>
<div class = "comment">
Note that the Fp block is identically zero unless we are using
some kind of artificial compressibility scheme...


<br><br>Matrix contributions for the uu and vv couplings.
</div>

<div class ="fragment">
<pre>
                      for (unsigned int j=0; j&lt;n_u_dofs; j++)
                        {
                          Kuu(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp] +                // mass matrix term
                                               theta*dt*(dphi[i][qp]*dphi[j][qp]) +   // diffusion term
                                               theta*dt*(U*dphi[j][qp])*phi[i][qp] +  // convection term
                                               theta*dt*u_x*phi[i][qp]*phi[j][qp]);   // Newton term
        
                          Kuv(i,j) += JxW[qp]*theta*dt*u_y*phi[i][qp]*phi[j][qp];     // Newton term
                          
                          Kvv(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp] +                // mass matrix term
                                               theta*dt*(dphi[i][qp]*dphi[j][qp]) +   // diffusion term
                                               theta*dt*(U*dphi[j][qp])*phi[i][qp] +  // convection term
                                               theta*dt*v_y*phi[i][qp]*phi[j][qp]);   // Newton term
        
                          Kvu(i,j) += JxW[qp]*theta*dt*v_x*phi[i][qp]*phi[j][qp];     // Newton term
                        }
        
</pre>
</div>
<div class = "comment">
Matrix contributions for the up and vp couplings.
</div>

<div class ="fragment">
<pre>
                      for (unsigned int j=0; j&lt;n_p_dofs; j++)
                        {
                          Kup(i,j) += JxW[qp]*(-theta*dt*psi[j][qp]*dphi[i][qp](0));
                          Kvp(i,j) += JxW[qp]*(-theta*dt*psi[j][qp]*dphi[i][qp](1));
                        }
                    }
        
</pre>
</div>
<div class = "comment">
Now an i-loop over the pressure degrees of freedom.  This code computes
the matrix entries due to the continuity equation.  Note: To maintain a
symmetric matrix, we may (or may not) multiply the continuity equation by
negative one.  Here we do not.
</div>

<div class ="fragment">
<pre>
                  for (unsigned int i=0; i&lt;n_p_dofs; i++)
                    {
                      Kp_alpha(i,0) += JxW[qp]*psi[i][qp];
                      Kalpha_p(0,i) += JxW[qp]*psi[i][qp];
                      for (unsigned int j=0; j&lt;n_u_dofs; j++)
                        {
                          Kpu(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](0);
                          Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1);
                        }
                    }
                } // end of the quadrature point qp-loop
        
              
</pre>
</div>
<div class = "comment">
At this point the interior element integration has
been completed.  However, we have not yet addressed
boundary conditions.  For this example we will only
consider simple Dirichlet boundary conditions imposed
via the penalty method. The penalty method used here
is equivalent (for Lagrange basis functions) to lumping
the matrix resulting from the L2 projection penalty
approach introduced in example 3.
</div>

<div class ="fragment">
<pre>
              {
</pre>
</div>
<div class = "comment">
The penalty value.  \f$ \frac{1}{\epsilon} \f$
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
                      AutoPtr&lt;Elem&gt; side (elem-&gt;build_side(s));
                                    
</pre>
</div>
<div class = "comment">
Loop over the nodes on the side.
</div>

<div class ="fragment">
<pre>
                      for (unsigned int ns=0; ns&lt;side-&gt;n_nodes(); ns++)
                        {
</pre>
</div>
<div class = "comment">
Boundary ids are set internally by
build_square().
0=bottom
1=right
2=top
3=left
                   

<br><br>Set u = 1 on the top boundary, 0 everywhere else
</div>

<div class ="fragment">
<pre>
                          const Real u_value = (mesh.boundary_info-&gt;has_boundary_id(elem,s,2)) ? 1. : 0.;
                          
</pre>
</div>
<div class = "comment">
Set v = 0 everywhere
</div>

<div class ="fragment">
<pre>
                          const Real v_value = 0.;
                          
</pre>
</div>
<div class = "comment">
Find the node on the element matching this node on
the side.  That defined where in the element matrix
the boundary condition will be applied.
</div>

<div class ="fragment">
<pre>
                          for (unsigned int n=0; n&lt;elem-&gt;n_nodes(); n++)
                            if (elem-&gt;node(n) == side-&gt;node(ns))
                              {
</pre>
</div>
<div class = "comment">
Matrix contribution.
</div>

<div class ="fragment">
<pre>
                                Kuu(n,n) += penalty;
                                Kvv(n,n) += penalty;
                                            
</pre>
</div>
<div class = "comment">
Right-hand-side contribution.
</div>

<div class ="fragment">
<pre>
                                Fu(n) += penalty*u_value;
                                Fv(n) += penalty*v_value;
                              }
                        } // end face node loop          
                    } // end if (elem-&gt;neighbor(side) == NULL)
              } // end boundary condition section          
        
</pre>
</div>
<div class = "comment">
If this assembly program were to be used on an adaptive mesh,
we would have to apply any hanging node constraint equations
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
              navier_stokes_system.matrix-&gt;add_matrix (Ke, dof_indices);
              navier_stokes_system.rhs-&gt;add_vector    (Fe, dof_indices);
            } // end of element loop
        
</pre>
</div>
<div class = "comment">
We can set the mean of the pressure by setting Falpha
</div>

<div class ="fragment">
<pre>
          navier_stokes_system.rhs-&gt;add(navier_stokes_system.rhs-&gt;size()-1,10.);
        
</pre>
</div>
<div class = "comment">
That's it.
</div>

<div class ="fragment">
<pre>
          return;
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file systems_of_equations_ex3.C without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;sstream&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/transient_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/perf_log.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/boundary_info.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/utility.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_submatrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_subvector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_stokes (EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
    
    <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::default_solver_package() == TRILINOS_SOLVERS)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;We skip example 13 when using the Trilinos solvers.\n&quot;</FONT></B>
                  &lt;&lt; std::endl;
        <B><FONT COLOR="#A020F0">return</FONT></B> 0;
      }
  
    Mesh mesh;
    
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square (mesh,
                                         20, 20,
                                         0., 1.,
                                         0., 1.,
                                         QUAD9);
    
    mesh.print_info();
    
    EquationSystems equation_systems (mesh);
    
    TransientLinearImplicitSystem &amp; system = 
      equation_systems.add_system&lt;TransientLinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Navier-Stokes&quot;</FONT></B>);
    
    system.add_variable (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, SECOND);
    system.add_variable (<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>, SECOND);
  
    system.add_variable (<B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>, FIRST);
  
    system.add_variable (<B><FONT COLOR="#BC8F8F">&quot;alpha&quot;</FONT></B>, FIRST, SCALAR);
  
    system.attach_assemble_function (assemble_stokes);
    
    equation_systems.init ();
  
    equation_systems.print_info();
  
    PerfLog perf_log(<B><FONT COLOR="#BC8F8F">&quot;Systems Example 3&quot;</FONT></B>);
    
    TransientLinearImplicitSystem&amp;  navier_stokes_system =
          equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Navier-Stokes&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real dt = 0.01;
    navier_stokes_system.time     = 0.0;
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_timesteps = 15;
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_nonlinear_steps = 15;
    <B><FONT COLOR="#228B22">const</FONT></B> Real nonlinear_tolerance       = 1.e-3;
  
    equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;linear solver maximum iterations&quot;</FONT></B>) = 250;
    
    equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;dt&quot;</FONT></B>)   = dt;
  
    AutoPtr&lt;NumericVector&lt;Number&gt; &gt;
      last_nonlinear_soln (navier_stokes_system.solution-&gt;clone());
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> t_step=0; t_step&lt;n_timesteps; ++t_step)
      {
        navier_stokes_system.time += dt;
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n\n*** Solving time step &quot;</FONT></B> &lt;&lt; t_step &lt;&lt; 
                     <B><FONT COLOR="#BC8F8F">&quot;, time = &quot;</FONT></B> &lt;&lt; navier_stokes_system.time &lt;&lt;
                     <B><FONT COLOR="#BC8F8F">&quot; ***&quot;</FONT></B> &lt;&lt; std::endl;
  
        *navier_stokes_system.old_local_solution = *navier_stokes_system.current_local_solution;
  
        <B><FONT COLOR="#228B22">const</FONT></B> Real initial_linear_solver_tol = 1.e-6;
        equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;linear solver tolerance&quot;</FONT></B>) = initial_linear_solver_tol;
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> l=0; l&lt;n_nonlinear_steps; ++l)
          {
            last_nonlinear_soln-&gt;zero();
            last_nonlinear_soln-&gt;add(*navier_stokes_system.solution);
            
            perf_log.push(<B><FONT COLOR="#BC8F8F">&quot;linear solve&quot;</FONT></B>);
            equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Navier-Stokes&quot;</FONT></B>).solve();
            perf_log.pop(<B><FONT COLOR="#BC8F8F">&quot;linear solve&quot;</FONT></B>);
  
            last_nonlinear_soln-&gt;add (-1., *navier_stokes_system.solution);
  
            last_nonlinear_soln-&gt;close();
  
            <B><FONT COLOR="#228B22">const</FONT></B> Real norm_delta = last_nonlinear_soln-&gt;l2_norm();
  
            <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_linear_iterations = navier_stokes_system.n_linear_iterations();
            
            <B><FONT COLOR="#228B22">const</FONT></B> Real final_linear_residual = navier_stokes_system.final_linear_residual();
            
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Linear solver converged at step: &quot;</FONT></B>
                      &lt;&lt; n_linear_iterations
                      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, final residual: &quot;</FONT></B>
                      &lt;&lt; final_linear_residual
                      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;  Nonlinear convergence: ||u - u_old|| = &quot;</FONT></B>
                      &lt;&lt; norm_delta
                      &lt;&lt; std::endl;
  
            <B><FONT COLOR="#A020F0">if</FONT></B> ((norm_delta &lt; nonlinear_tolerance) &amp;&amp;
                (navier_stokes_system.final_linear_residual() &lt; nonlinear_tolerance))
              {
                <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Nonlinear solver converged at step &quot;</FONT></B>
                          &lt;&lt; l
                          &lt;&lt; std::endl;
                <B><FONT COLOR="#A020F0">break</FONT></B>;
              }
            
            equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;linear solver tolerance&quot;</FONT></B>) =
              <B><FONT COLOR="#5F9EA0">std</FONT></B>::min(Utility::pow&lt;2&gt;(final_linear_residual), initial_linear_solver_tol);
  
          } <I><FONT COLOR="#B22222">// end nonlinear loop
</FONT></I>        
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> write_interval = 1;
        
  #ifdef LIBMESH_HAVE_EXODUS_API
        <B><FONT COLOR="#A020F0">if</FONT></B> ((t_step+1)%write_interval == 0)
          {
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream file_name;
  
            file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out_&quot;</FONT></B>
                      &lt;&lt; std::setw(3)
                      &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
                      &lt;&lt; std::right
                      &lt;&lt; t_step + 1
                      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.e&quot;</FONT></B>;
            
            ExodusII_IO(mesh).write_equation_systems (file_name.str(),
                                                equation_systems);
          }
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>      } <I><FONT COLOR="#B22222">// end timestep loop.
</FONT></I>    
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_stokes (EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Navier-Stokes&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
    
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
    
    TransientLinearImplicitSystem &amp; navier_stokes_system =
      es.get_system&lt;TransientLinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Navier-Stokes&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = navier_stokes_system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> v_var = navier_stokes_system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> p_var = navier_stokes_system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> alpha_var = navier_stokes_system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;alpha&quot;</FONT></B>);
    
    FEType fe_vel_type = navier_stokes_system.variable_type(u_var);
    
    FEType fe_pres_type = navier_stokes_system.variable_type(p_var);
  
    AutoPtr&lt;FEBase&gt; fe_vel  (FEBase::build(dim, fe_vel_type));
      
    AutoPtr&lt;FEBase&gt; fe_pres (FEBase::build(dim, fe_pres_type));
    
    QGauss qrule (dim, fe_vel_type.default_quadrature_order());
  
    fe_vel-&gt;attach_quadrature_rule (&amp;qrule);
    fe_pres-&gt;attach_quadrature_rule (&amp;qrule);
    
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe_vel-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe_vel-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe_vel-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; psi = fe_pres-&gt;get_phi();
  
    
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap &amp; dof_map = navier_stokes_system.get_dof_map();
  
    DenseMatrix&lt;Number&gt; Ke;
    DenseVector&lt;Number&gt; Fe;
  
    DenseSubMatrix&lt;Number&gt;
      Kuu(Ke), Kuv(Ke), Kup(Ke),
      Kvu(Ke), Kvv(Ke), Kvp(Ke),
      Kpu(Ke), Kpv(Ke), Kpp(Ke);
    DenseSubMatrix&lt;Number&gt; Kalpha_p(Ke), Kp_alpha(Ke);
  
    DenseSubVector&lt;Number&gt;
      Fu(Fe),
      Fv(Fe),
      Fp(Fe);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices_u;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices_v;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices_p;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices_alpha;
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real dt    = es.parameters.get&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;dt&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> Real theta = 1.;
      
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
    
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {    
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        dof_map.dof_indices (elem, dof_indices_u, u_var);
        dof_map.dof_indices (elem, dof_indices_v, v_var);
        dof_map.dof_indices (elem, dof_indices_p, p_var);
        dof_map.dof_indices (elem, dof_indices_alpha, alpha_var);
  
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_dofs   = dof_indices.size();
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = dof_indices_u.size(); 
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_v_dofs = dof_indices_v.size(); 
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_p_dofs = dof_indices_p.size();
        
        fe_vel-&gt;reinit  (elem);
        fe_pres-&gt;reinit (elem);
  
        Ke.resize (n_dofs, n_dofs);
        Fe.resize (n_dofs);
  
        Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
        Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);
        
        Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
        Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
        Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);
        
        Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
        Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
        Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);
  
        Kp_alpha.reposition (p_var*n_u_dofs, p_var*n_u_dofs+n_p_dofs, n_p_dofs, 1);
        Kalpha_p.reposition (p_var*n_u_dofs+n_p_dofs, p_var*n_u_dofs, 1, n_p_dofs);
  
  
        Fu.reposition (u_var*n_u_dofs, n_u_dofs);
        Fv.reposition (v_var*n_u_dofs, n_v_dofs);
        Fp.reposition (p_var*n_u_dofs, n_p_dofs);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
          {
            Number   u = 0., u_old = 0.;
            Number   v = 0., v_old = 0.;
            Number   p_old = 0.;
            Gradient grad_u, grad_u_old;
            Gradient grad_v, grad_v_old;
            
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> l=0; l&lt;n_u_dofs; l++)
              {
                u_old += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_u[l]);
                v_old += phi[l][qp]*navier_stokes_system.old_solution (dof_indices_v[l]);
                grad_u_old.add_scaled (dphi[l][qp],navier_stokes_system.old_solution (dof_indices_u[l]));
                grad_v_old.add_scaled (dphi[l][qp],navier_stokes_system.old_solution (dof_indices_v[l]));
  
                u += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_u[l]); 
                v += phi[l][qp]*navier_stokes_system.current_solution (dof_indices_v[l]);
                grad_u.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_u[l]));
                grad_v.add_scaled (dphi[l][qp],navier_stokes_system.current_solution (dof_indices_v[l]));
              }
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> l=0; l&lt;n_p_dofs; l++)
              {
                p_old += psi[l][qp]*navier_stokes_system.old_solution (dof_indices_p[l]);
              }
  
            <B><FONT COLOR="#228B22">const</FONT></B> NumberVectorValue U_old (u_old, v_old);
            <B><FONT COLOR="#228B22">const</FONT></B> NumberVectorValue U     (u,     v);
            <B><FONT COLOR="#228B22">const</FONT></B> Number  u_x = grad_u(0);
            <B><FONT COLOR="#228B22">const</FONT></B> Number  u_y = grad_u(1);
            <B><FONT COLOR="#228B22">const</FONT></B> Number  v_x = grad_v(0);
            <B><FONT COLOR="#228B22">const</FONT></B> Number  v_y = grad_v(1);
            
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_u_dofs; i++)
              {
                Fu(i) += JxW[qp]*(u_old*phi[i][qp] -                            <I><FONT COLOR="#B22222">// mass-matrix term 
</FONT></I>                                  (1.-theta)*dt*(U_old*grad_u_old)*phi[i][qp] + <I><FONT COLOR="#B22222">// convection term
</FONT></I>                                  (1.-theta)*dt*p_old*dphi[i][qp](0)  -         <I><FONT COLOR="#B22222">// pressure term on rhs
</FONT></I>                                  (1.-theta)*dt*(grad_u_old*dphi[i][qp]) +      <I><FONT COLOR="#B22222">// diffusion term on rhs
</FONT></I>                                  theta*dt*(U*grad_u)*phi[i][qp]);              <I><FONT COLOR="#B22222">// Newton term
</FONT></I>  
                  
                Fv(i) += JxW[qp]*(v_old*phi[i][qp] -                             <I><FONT COLOR="#B22222">// mass-matrix term
</FONT></I>                                  (1.-theta)*dt*(U_old*grad_v_old)*phi[i][qp] +  <I><FONT COLOR="#B22222">// convection term
</FONT></I>                                  (1.-theta)*dt*p_old*dphi[i][qp](1) -           <I><FONT COLOR="#B22222">// pressure term on rhs
</FONT></I>                                  (1.-theta)*dt*(grad_v_old*dphi[i][qp]) +       <I><FONT COLOR="#B22222">// diffusion term on rhs
</FONT></I>                                  theta*dt*(U*grad_v)*phi[i][qp]);               <I><FONT COLOR="#B22222">// Newton term
</FONT></I>              
  
  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_u_dofs; j++)
                  {
                    Kuu(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp] +                <I><FONT COLOR="#B22222">// mass matrix term
</FONT></I>                                         theta*dt*(dphi[i][qp]*dphi[j][qp]) +   <I><FONT COLOR="#B22222">// diffusion term
</FONT></I>                                         theta*dt*(U*dphi[j][qp])*phi[i][qp] +  <I><FONT COLOR="#B22222">// convection term
</FONT></I>                                         theta*dt*u_x*phi[i][qp]*phi[j][qp]);   <I><FONT COLOR="#B22222">// Newton term
</FONT></I>  
                    Kuv(i,j) += JxW[qp]*theta*dt*u_y*phi[i][qp]*phi[j][qp];     <I><FONT COLOR="#B22222">// Newton term
</FONT></I>                    
                    Kvv(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp] +                <I><FONT COLOR="#B22222">// mass matrix term
</FONT></I>                                         theta*dt*(dphi[i][qp]*dphi[j][qp]) +   <I><FONT COLOR="#B22222">// diffusion term
</FONT></I>                                         theta*dt*(U*dphi[j][qp])*phi[i][qp] +  <I><FONT COLOR="#B22222">// convection term
</FONT></I>                                         theta*dt*v_y*phi[i][qp]*phi[j][qp]);   <I><FONT COLOR="#B22222">// Newton term
</FONT></I>  
                    Kvu(i,j) += JxW[qp]*theta*dt*v_x*phi[i][qp]*phi[j][qp];     <I><FONT COLOR="#B22222">// Newton term
</FONT></I>                  }
  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_p_dofs; j++)
                  {
                    Kup(i,j) += JxW[qp]*(-theta*dt*psi[j][qp]*dphi[i][qp](0));
                    Kvp(i,j) += JxW[qp]*(-theta*dt*psi[j][qp]*dphi[i][qp](1));
                  }
              }
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_p_dofs; i++)
              {
                Kp_alpha(i,0) += JxW[qp]*psi[i][qp];
                Kalpha_p(0,i) += JxW[qp]*psi[i][qp];
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_u_dofs; j++)
                  {
                    Kpu(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](0);
                    Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1);
                  }
              }
          } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>  
        
        {
          <B><FONT COLOR="#228B22">const</FONT></B> Real penalty = 1.e10;
                    
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> s=0; s&lt;elem-&gt;n_sides(); s++)
            <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(s) == NULL)
              {
                AutoPtr&lt;Elem&gt; side (elem-&gt;build_side(s));
                              
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> ns=0; ns&lt;side-&gt;n_nodes(); ns++)
                  {
                     
                    <B><FONT COLOR="#228B22">const</FONT></B> Real u_value = (mesh.boundary_info-&gt;has_boundary_id(elem,s,2)) ? 1. : 0.;
                    
                    <B><FONT COLOR="#228B22">const</FONT></B> Real v_value = 0.;
                    
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n=0; n&lt;elem-&gt;n_nodes(); n++)
                      <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;node(n) == side-&gt;node(ns))
                        {
                          Kuu(n,n) += penalty;
                          Kvv(n,n) += penalty;
                                      
                          Fu(n) += penalty*u_value;
                          Fv(n) += penalty*v_value;
                        }
                  } <I><FONT COLOR="#B22222">// end face node loop          
</FONT></I>              } <I><FONT COLOR="#B22222">// end if (elem-&gt;neighbor(side) == NULL)
</FONT></I>        } <I><FONT COLOR="#B22222">// end boundary condition section          
</FONT></I>  
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
  
        navier_stokes_system.matrix-&gt;add_matrix (Ke, dof_indices);
        navier_stokes_system.rhs-&gt;add_vector    (Fe, dof_indices);
      } <I><FONT COLOR="#B22222">// end of element loop
</FONT></I>  
    navier_stokes_system.rhs-&gt;add(navier_stokes_system.rhs-&gt;size()-1,10.);
  
    <B><FONT COLOR="#A020F0">return</FONT></B>;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example systems_of_equations_ex3:
*  mpirun -np 12 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1681
    n_local_nodes()=165
  n_elem()=400
    n_local_elem()=34
    n_active_elem()=400
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Navier-Stokes"
    Type "TransientLinearImplicit"
    Variables={ "u" "v" } "p" "alpha" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" "SCALAR", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" "CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" "FIRST", "THIRD" "FIRST", "THIRD" 
    n_dofs()=3804
    n_local_dofs()=379
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 35.1872
      Average Off-Processor Bandwidth <= 7.65405
      Maximum  On-Processor Bandwidth <= 253
      Maximum Off-Processor Bandwidth <= 3551
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0



*** Solving time step 0, time = 0.01 ***
Linear solver converged at step: 51, final residual: 0.00830395  Nonlinear convergence: ||u - u_old|| = 2535.68
Linear solver converged at step: 17, final residual: 0.00913569  Nonlinear convergence: ||u - u_old|| = 0.683391
Linear solver converged at step: 0, final residual: 0.009136  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.009136  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.009136  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.009136  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.009136  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.009136  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.009136  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.009136  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.009136  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.009136  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.009136  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.009136  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.009136  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 1, time = 0.02 ***
Linear solver converged at step: 45, final residual: 0.0103502  Nonlinear convergence: ||u - u_old|| = 28.9276
Linear solver converged at step: 10, final residual: 0.0107284  Nonlinear convergence: ||u - u_old|| = 0.032819
Linear solver converged at step: 0, final residual: 0.0107284  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107284  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107284  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107284  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107284  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107284  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107284  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107284  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107284  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107284  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107284  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107284  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107284  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 2, time = 0.03 ***
Linear solver converged at step: 40, final residual: 0.00980235  Nonlinear convergence: ||u - u_old|| = 5.48348
Linear solver converged at step: 8, final residual: 0.0107636  Nonlinear convergence: ||u - u_old|| = 0.0121513
Linear solver converged at step: 0, final residual: 0.0107636  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107636  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107636  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107636  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107636  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107636  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107636  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107636  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107636  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107636  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107636  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107636  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107636  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 3, time = 0.04 ***
Linear solver converged at step: 28, final residual: 0.0106529  Nonlinear convergence: ||u - u_old|| = 2.02368
Linear solver converged at step: 6, final residual: 0.0107949  Nonlinear convergence: ||u - u_old|| = 0.00156818
Linear solver converged at step: 0, final residual: 0.0107949  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107949  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107949  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107949  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107949  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107949  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107949  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107949  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107949  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107949  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107949  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107949  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107949  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 4, time = 0.05 ***
Linear solver converged at step: 27, final residual: 0.0080356  Nonlinear convergence: ||u - u_old|| = 0.908583
Linear solver converged at step: 0, final residual: 0.00805066  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00805066  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00805066  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00805066  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00805066  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00805066  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00805066  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00805066  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00805066  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00805066  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00805066  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00805066  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00805066  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00805066  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 5, time = 0.06 ***
Linear solver converged at step: 25, final residual: 0.00901576  Nonlinear convergence: ||u - u_old|| = 0.436349
Linear solver converged at step: 0, final residual: 0.00901236  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00901236  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00901236  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00901236  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00901236  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00901236  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00901236  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00901236  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00901236  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00901236  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00901236  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00901236  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00901236  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00901236  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 6, time = 0.07 ***
Linear solver converged at step: 23, final residual: 0.0100393  Nonlinear convergence: ||u - u_old|| = 0.239729
Linear solver converged at step: 0, final residual: 0.0100309  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0100309  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0100309  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0100309  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0100309  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0100309  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0100309  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0100309  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0100309  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0100309  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0100309  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0100309  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0100309  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0100309  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 7, time = 0.08 ***
Linear solver converged at step: 21, final residual: 0.0093876  Nonlinear convergence: ||u - u_old|| = 0.134052
Linear solver converged at step: 0, final residual: 0.00938764  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00938764  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00938764  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00938764  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00938764  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00938764  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00938764  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00938764  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00938764  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00938764  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00938764  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00938764  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00938764  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.00938764  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 8, time = 0.09 ***
Linear solver converged at step: 18, final residual: 0.0102795  Nonlinear convergence: ||u - u_old|| = 0.0665637
Linear solver converged at step: 0, final residual: 0.0102792  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102792  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102792  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102792  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102792  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102792  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102792  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102792  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102792  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102792  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102792  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102792  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102792  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102792  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 9, time = 0.1 ***
Linear solver converged at step: 15, final residual: 0.0107515  Nonlinear convergence: ||u - u_old|| = 0.0355621
Linear solver converged at step: 0, final residual: 0.0107514  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107514  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107514  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107514  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107514  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107514  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107514  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107514  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107514  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107514  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107514  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107514  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107514  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0107514  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 10, time = 0.11 ***
Linear solver converged at step: 15, final residual: 0.0104437  Nonlinear convergence: ||u - u_old|| = 0.0295371
Linear solver converged at step: 0, final residual: 0.0104436  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0104436  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0104436  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0104436  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0104436  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0104436  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0104436  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0104436  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0104436  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0104436  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0104436  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0104436  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0104436  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0104436  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 11, time = 0.12 ***
Linear solver converged at step: 10, final residual: 0.0108009  Nonlinear convergence: ||u - u_old|| = 0.0196576
Linear solver converged at step: 0, final residual: 0.0108008  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0108008  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0108008  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0108008  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0108008  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0108008  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0108008  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0108008  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0108008  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0108008  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0108008  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0108008  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0108008  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0108008  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 12, time = 0.13 ***
Linear solver converged at step: 12, final residual: 0.0102339  Nonlinear convergence: ||u - u_old|| = 0.0171197
Linear solver converged at step: 0, final residual: 0.0102339  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102339  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102339  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102339  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102339  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102339  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102339  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102339  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102339  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102339  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102339  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102339  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102339  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102339  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 13, time = 0.14 ***
Linear solver converged at step: 7, final residual: 0.0102967  Nonlinear convergence: ||u - u_old|| = 0.00818877
Linear solver converged at step: 0, final residual: 0.0102967  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102967  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102967  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102967  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102967  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102967  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102967  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102967  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102967  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102967  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102967  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102967  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102967  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0102967  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 14, time = 0.15 ***
Linear solver converged at step: 5, final residual: 0.010699  Nonlinear convergence: ||u - u_old|| = 0.00520447
Linear solver converged at step: 0, final residual: 0.010699  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.010699  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.010699  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.010699  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.010699  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.010699  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.010699  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.010699  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.010699  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.010699  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.010699  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.010699  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.010699  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.010699  Nonlinear convergence: ||u - u_old|| = 0

 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 22:14:20 2013                                                                             |
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
 -----------------------------------------------------------------------------------------------------------
| Systems Example 3 Performance: Alive time=18.6212, Active time=18.2194                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| linear solve                  225       18.2194     0.080975    18.2194     0.080975    100.00   100.00   |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       225       18.2194                                         100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/systems_of_equations/systems_of_equations_ex3/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:14:20 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           1.890e+01      1.00020   1.890e+01
Objects:              1.657e+03      1.00000   1.657e+03
Flops:                9.786e+08      2.32944   7.057e+08  8.468e+09
Flops/sec:            5.179e+07      2.32944   3.735e+07  4.481e+08
MPI Messages:         1.706e+04      4.99648   7.398e+03  8.878e+04
MPI Message Lengths:  3.276e+07      4.51605   1.592e+03  1.414e+08
MPI Reductions:       6.726e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.8897e+01 100.0%  8.4685e+09 100.0%  8.878e+04 100.0%  1.592e+03      100.0%  6.725e+03 100.0% 

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

VecMDot              383 1.0 1.5422e-02 3.6 3.23e+06 1.5 0.0e+00 0.0e+00 3.8e+02  0  0  0  0  6   0  0  0  0  6  2103
VecNorm             1061 1.0 6.9877e-0175.5 8.04e+05 1.5 0.0e+00 0.0e+00 1.1e+03  1  0  0  0 16   1  0  0  0 16    12
VecScale             611 1.0 8.1897e-04 1.7 2.32e+05 1.5 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2838
VecCopy              469 1.0 9.0313e-04 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet              1316 1.0 1.0006e-03 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              700 1.0 2.0205e-02 2.5 5.31e+05 1.5 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   264
VecMAXPY             405 1.0 2.1396e-03 1.5 3.53e+06 1.5 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0 16542
VecAssemblyBegin     916 1.0 9.5586e-02 4.3 0.00e+00 0.0 1.5e+04 2.5e+02 2.7e+03  0  0 17  3 41   0  0 17  3 41     0
VecAssemblyEnd       916 1.0 3.4006e-03 3.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      851 1.0 7.8137e-03 1.8 0.00e+00 0.0 5.0e+04 6.3e+02 0.0e+00  0  0 56 22  0   0  0 56 22  0     0
VecScatterEnd        851 1.0 1.9810e+00652.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  4  0  0  0  0   4  0  0  0  0     0
VecNormalize         611 1.0 6.9132e-01100.9 6.95e+05 1.5 0.0e+00 0.0e+00 6.1e+02  1  0  0  0  9   1  0  0  0  9    10
MatMult              611 1.0 1.9986e+0098.8 1.81e+07 1.5 3.7e+04 7.0e+02 0.0e+00  4  2 41 18  0   4  2 41 18  0    94
MatSolve             836 1.0 5.7770e-02 1.7 8.02e+07 1.7 0.0e+00 0.0e+00 0.0e+00  0  9  0  0  0   0  9  0  0  0 13518
MatLUFactorNum       225 1.0 8.6310e-01 2.4 8.75e+08 2.5 0.0e+00 0.0e+00 0.0e+00  3 88  0  0  0   3 88  0  0  0  8592
MatILUFactorSym      225 1.0 2.6720e+00 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 6.8e+02 10  0  0  0 10  10  0  0  0 10     0
MatAssemblyBegin     450 1.0 1.3946e+00 5.9 0.00e+00 0.0 2.3e+04 4.6e+03 9.0e+02  5  0 26 75 13   5  0 26 75 13     0
MatAssemblyEnd       450 1.0 2.0228e-0115.3 0.00e+00 0.0 1.2e+02 1.8e+02 8.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ          225 1.0 8.0824e-05 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering       225 1.0 1.6046e-02 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 9.0e+02  0  0  0  0 13   0  0  0  0 13     0
MatZeroEntries       227 1.0 3.8443e-03 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       383 1.0 1.7610e-02 2.7 6.47e+06 1.5 0.0e+00 0.0e+00 3.8e+02  0  1  0  0  6   0  1  0  0  6  3686
KSPSetUp             450 1.0 2.0919e-03 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve             225 1.0 3.7292e+00 1.0 9.78e+08 2.3 3.7e+04 7.0e+02 2.8e+03 20100 41 18 42  20100 41 18 42  2269
PCSetUp              450 1.0 3.6121e+00 2.2 8.75e+08 2.5 0.0e+00 0.0e+00 1.6e+03 14 88  0  0 23  14 88  0  0 23  2053
PCSetUpOnBlocks      225 1.0 3.6105e+00 2.2 8.75e+08 2.5 0.0e+00 0.0e+00 1.6e+03 14 88  0  0 23  14 88  0  0 23  2054
PCApply              836 1.0 7.0836e-02 1.6 8.02e+07 1.7 0.0e+00 0.0e+00 0.0e+00  0  9  0  0  0   0  9  0  0  0 11025
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector   276            276      1228936     0
      Vector Scatter     6              6         6216     0
           Index Set  1137           1137      1709784     0
   IS L to G Mapping     5              5         2820     0
              Matrix   228            228    121949428     0
       Krylov Solver     2              2        19360     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 3.19481e-06
Average time for zero size MPI_Send(): 1.32521e-05
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

 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=19.0152, Active time=18.7757                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   SCALAR_dof_indices()             15873     0.1117      0.000007    0.1117      0.000007    0.59     0.59     |
|   add_neighbors_to_send_list()     1         0.0149      0.014948    0.0276      0.027619    0.08     0.15     |
|   build_sparsity()                 1         0.0141      0.014146    0.0324      0.032424    0.08     0.17     |
|   create_dof_constraints()         1         0.0015      0.001519    0.0015      0.001519    0.01     0.01     |
|   distribute_dofs()                1         0.0363      0.036294    0.0896      0.089579    0.19     0.48     |
|   dof_indices()                    40353     6.9097      0.000171    7.0258      0.000174    36.80    37.42    |
|   prepare_send_list()              1         0.0002      0.000180    0.0002      0.000180    0.00     0.00     |
|   reinit()                         1         0.0515      0.051497    0.0515      0.051497    0.27     0.27     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          15        0.0283      0.001890    0.2624      0.017491    0.15     1.40     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               15        0.0813      0.005418    0.0813      0.005418    0.43     0.43     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        15300     0.5158      0.000034    0.5158      0.000034    2.75     2.75     |
|   init_shape_functions()           450       0.0260      0.000058    0.0260      0.000058    0.14     0.14     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             15300     0.2764      0.000018    0.2764      0.000018    1.47     1.47     |
|   init_reference_to_physical_map() 450       0.0577      0.000128    0.0577      0.000128    0.31     0.31     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0169      0.016851    0.0179      0.017855    0.09     0.10     |
|   renumber_nodes_and_elem()        2         0.0009      0.000440    0.0009      0.000440    0.00     0.00     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0068      0.003417    0.0068      0.003417    0.04     0.04     |
|   find_global_indices()            2         0.0030      0.001486    0.0146      0.007313    0.02     0.08     |
|   parallel_sort()                  2         0.0030      0.001476    0.0035      0.001762    0.02     0.02     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         15        0.0012      0.000082    0.3469      0.023130    0.01     1.85     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0046      0.004608    0.0046      0.004608    0.02     0.02     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0430      0.042960    0.0494      0.049399    0.23     0.26     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      9         0.0004      0.000041    0.0004      0.000044    0.00     0.00     |
|   max(bool)                        1         0.0000      0.000007    0.0000      0.000007    0.00     0.00     |
|   max(scalar)                      315       0.0030      0.000010    0.0030      0.000010    0.02     0.02     |
|   max(vector)                      66        0.0010      0.000016    0.0029      0.000043    0.01     0.02     |
|   min(bool)                        373       0.0031      0.000008    0.0031      0.000008    0.02     0.02     |
|   min(scalar)                      309       0.0191      0.000062    0.0191      0.000062    0.10     0.10     |
|   min(vector)                      66        0.0011      0.000017    0.0029      0.000045    0.01     0.02     |
|   probe()                          132       0.0017      0.000013    0.0017      0.000013    0.01     0.01     |
|   receive()                        132       0.0009      0.000007    0.0027      0.000020    0.00     0.01     |
|   send()                           132       0.0005      0.000004    0.0005      0.000004    0.00     0.00     |
|   send_receive()                   136       0.0013      0.000009    0.0048      0.000035    0.01     0.03     |
|   sum()                            62        0.0058      0.000093    0.0076      0.000122    0.03     0.04     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           132       0.0003      0.000002    0.0003      0.000002    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0018      0.001784    0.0026      0.002551    0.01     0.01     |
|   set_parent_processor_ids()       1         0.0008      0.000762    0.0008      0.000762    0.00     0.00     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          225       4.6870      0.020831    4.6870      0.020831    24.96    24.96    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       225       5.8432      0.025970    13.5121     0.060054    31.12    71.97    |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            90105     18.7757                                         100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example systems_of_equations_ex3:
*  mpirun -np 12 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
