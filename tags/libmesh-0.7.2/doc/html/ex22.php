<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex22",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 22 - Unsteady Navier-Stokes Equations - SCALAR variables</h1>

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
        #include &lt;math.h&gt;
        
</pre>
</div>
<div class = "comment">
Basic include file needed for the mesh functionality.
</div>

<div class ="fragment">
<pre>
        #include "libmesh.h"
        #include "mesh.h"
        #include "mesh_generation.h"
        #include "exodusII_io.h"
        #include "equation_systems.h"
        #include "fe.h"
        #include "quadrature_gauss.h"
        #include "dof_map.h"
        #include "sparse_matrix.h"
        #include "numeric_vector.h"
        #include "dense_matrix.h"
        #include "dense_vector.h"
        #include "linear_implicit_system.h"
        #include "transient_system.h"
        #include "perf_log.h"
        #include "boundary_info.h"
        #include "utility.h"
        
</pre>
</div>
<div class = "comment">
Some (older) compilers do not offer full stream 
functionality, OStringStream works around this.
</div>

<div class ="fragment">
<pre>
        #include "o_string_stream.h"
        
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
        #include "dense_submatrix.h"
        #include "dense_subvector.h"
        
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
          PerfLog perf_log("Example 13");
          
</pre>
</div>
<div class = "comment">
Now we begin the timestep loop to compute the time-accurate
solution of the equations.
</div>

<div class ="fragment">
<pre>
          const Real dt = 0.01;
          Real time     = 0.0;
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
Get a reference to the Stokes system to use later.
</div>

<div class ="fragment">
<pre>
          TransientLinearImplicitSystem&  navier_stokes_system =
                equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;("Navier-Stokes");
        
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
              time += dt;
        
</pre>
</div>
<div class = "comment">
Let the system of equations know the current time.
This might be necessary for a time-dependent forcing
function for example.
</div>

<div class ="fragment">
<pre>
              equation_systems.parameters.set&lt;Real&gt; ("time") = time;
        
</pre>
</div>
<div class = "comment">
A pretty update message
</div>

<div class ="fragment">
<pre>
              std::cout &lt;&lt; "\n\n*** Solving time step " &lt;&lt; t_step &lt;&lt; ", time = " &lt;&lt; time &lt;&lt; " ***" &lt;&lt; std::endl;
        
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
                  OStringStream file_name;
        
</pre>
</div>
<div class = "comment">
We write the file in the ExodusII format.
</div>

<div class ="fragment">
<pre>
                  file_name &lt;&lt; "out_";
                  OSSRealzeroright(file_name,3,0, t_step + 1);
                  file_name &lt;&lt; ".exd";
                  
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
          libmesh_assert (system_name == "Navier-Stokes");
          
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
          std::vector&lt;unsigned int&gt; dof_indices;
          std::vector&lt;unsigned int&gt; dof_indices_u;
          std::vector&lt;unsigned int&gt; dof_indices_v;
          std::vector&lt;unsigned int&gt; dof_indices_p;
          std::vector&lt;unsigned int&gt; dof_indices_alpha;
        
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
</pre>
</div>
<div class = "comment">
Get the boundary ID for side 's'.
These are set internally by build_square().
0=bottom
1=right
2=top
3=left
</div>

<div class ="fragment">
<pre>
                      short int bc_id = mesh.boundary_info-&gt;boundary_id (elem,s);
                      if (bc_id==BoundaryInfo::invalid_id)
                          libmesh_error();
        
                      
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
Get the boundary values.
                   

<br><br>Set u = 1 on the top boundary, 0 everywhere else
</div>

<div class ="fragment">
<pre>
                          const Real u_value = (bc_id==2) ? 1. : 0.;
                          
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
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;transient_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;perf_log.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;boundary_info.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;utility.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;o_string_stream.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_submatrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_subvector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  
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
  
    PerfLog perf_log(<B><FONT COLOR="#BC8F8F">&quot;Example 13&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> Real dt = 0.01;
    Real time     = 0.0;
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_timesteps = 15;
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_nonlinear_steps = 15;
    <B><FONT COLOR="#228B22">const</FONT></B> Real nonlinear_tolerance       = 1.e-3;
  
    equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;linear solver maximum iterations&quot;</FONT></B>) = 250;
    
    equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;dt&quot;</FONT></B>)   = dt;
  
    TransientLinearImplicitSystem&amp;  navier_stokes_system =
          equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Navier-Stokes&quot;</FONT></B>);
  
    AutoPtr&lt;NumericVector&lt;Number&gt; &gt;
      last_nonlinear_soln (navier_stokes_system.solution-&gt;clone());
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> t_step=0; t_step&lt;n_timesteps; ++t_step)
      {
        time += dt;
  
        equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>) = time;
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n\n*** Solving time step &quot;</FONT></B> &lt;&lt; t_step &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, time = &quot;</FONT></B> &lt;&lt; time &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; ***&quot;</FONT></B> &lt;&lt; std::endl;
  
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
            OStringStream file_name;
  
            file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out_&quot;</FONT></B>;
            OSSRealzeroright(file_name,3,0, t_step + 1);
            file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.exd&quot;</FONT></B>;
            
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
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Navier-Stokes&quot;</FONT></B>);
    
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
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices_u;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices_v;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices_p;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices_alpha;
  
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
                <B><FONT COLOR="#228B22">short</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> bc_id = mesh.boundary_info-&gt;boundary_id (elem,s);
                <B><FONT COLOR="#A020F0">if</FONT></B> (bc_id==BoundaryInfo::invalid_id)
                    libmesh_error();
  
                
                AutoPtr&lt;Elem&gt; side (elem-&gt;build_side(s));
                              
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> ns=0; ns&lt;side-&gt;n_nodes(); ns++)
                  {
                     
                    <B><FONT COLOR="#228B22">const</FONT></B> Real u_value = (bc_id==2) ? 1. : 0.;
                    
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
Compiling C++ (in optimized mode) ex22.C...
Linking ex22-opt...
***************************************************************
* Running Example  mpirun -np 2 ./ex22-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1681
    n_local_nodes()=863
  n_elem()=400
    n_local_elem()=200
    n_active_elem()=400
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Navier-Stokes"
    Type "TransientLinearImplicit"
    Variables="u" "v" "p" "alpha" 
    Finite Element Types="LAGRANGE" "LAGRANGE" "LAGRANGE" "SCALAR" 
    Approximation Orders="SECOND" "SECOND" "FIRST" "FIRST" 
    n_dofs()=3804
    n_local_dofs()=1958
    n_constrained_dofs()=0
    n_vectors()=3



*** Solving time step 0, time = 0.01 ***
Linear solver converged at step: 1, final residual: 0.176607  Nonlinear convergence: ||u - u_old|| = 692.321
Linear solver converged at step: 0, final residual: 0.176619  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176619  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176619  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176619  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176619  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176619  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176619  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176619  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176619  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176619  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176619  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176619  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176619  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176619  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 1, time = 0.02 ***
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 2, time = 0.03 ***
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 3, time = 0.04 ***
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 4, time = 0.05 ***
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 5, time = 0.06 ***
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 6, time = 0.07 ***
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 7, time = 0.08 ***
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 8, time = 0.09 ***
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 9, time = 0.1 ***
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 10, time = 0.11 ***
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 11, time = 0.12 ***
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 12, time = 0.13 ***
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 13, time = 0.14 ***
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 14, time = 0.15 ***
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.176894  Nonlinear convergence: ||u - u_old|| = 0

-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 2                                                |
| Time:           Tue Feb 22 12:23:41 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Thu Feb 17 16:40:46 CST 2011  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| Example 13 Performance: Alive time=94.0339, Active time=93.9979                                           |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| linear solve                  225       93.9979     0.417768    93.9979     0.417768    100.00   100.00   |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       225       93.9979                                         100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex22-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Tue Feb 22 12:23:41 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           9.406e+01      1.00002   9.406e+01
Objects:              1.392e+03      1.00000   1.392e+03
Flops:                4.905e+10      1.14177   4.600e+10  9.201e+10
Flops/sec:            5.215e+08      1.14175   4.891e+08  9.782e+08
MPI Messages:         2.484e+03      1.00000   2.484e+03  4.968e+03
MPI Message Lengths:  2.289e+07      1.00888   9.175e+03  4.558e+07
MPI Reductions:       6.586e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 9.4059e+01 100.0%  9.2010e+10 100.0%  4.968e+03 100.0%  9.175e+03      100.0%  6.116e+03  92.9% 

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

VecMDot                1 1.0 3.1948e-05 2.0 3.92e+03 1.1 0.0e+00 0.0e+00 1.0e+00  0  0  0  0  0   0  0  0  0  0   238
VecNorm              676 1.0 3.0718e-02 3.6 2.65e+06 1.1 0.0e+00 0.0e+00 6.8e+02  0  0  0  0 10   0  0  0  0 11   167
VecScale             226 1.0 4.3917e-04 1.1 4.43e+05 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1958
VecCopy              467 1.0 1.2431e-03 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               460 1.0 5.4622e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              676 1.0 5.2910e-03 1.1 2.65e+06 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   972
VecMAXPY               2 1.0 5.9605e-06 1.2 7.83e+03 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2553
VecAssemblyBegin     900 1.0 7.0212e-03 1.1 0.00e+00 0.0 9.0e+02 1.2e+03 2.7e+03  0  0 18  2 41   0  0 18  2 44     0
VecAssemblyEnd       900 1.0 8.5139e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      676 1.0 4.0455e-03 1.1 0.00e+00 0.0 9.0e+02 5.4e+03 0.0e+00  0  0 18 11  0   0  0 18 11  0     0
VecScatterEnd        676 1.0 9.5120e+007155.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  5  0  0  0  0   5  0  0  0  0     0
VecNormalize         226 1.0 2.7648e-02 5.0 1.33e+06 1.1 0.0e+00 0.0e+00 2.3e+02  0  0  0  0  3   0  0  0  0  4    93
MatMult              226 1.0 9.5411e+00297.0 3.52e+07 1.0 4.5e+02 8.7e+03 0.0e+00  5  0  9  9  0   5  0  9  9  0     7
MatSolve               2 1.0 2.6710e-03 1.1 3.33e+06 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2378
MatLUFactorNum       225 1.0 3.7326e+01 1.1 4.90e+10 1.1 0.0e+00 0.0e+00 0.0e+00 37100  0  0  0  37100  0  0  0  2463
MatILUFactorSym      225 1.0 5.4315e+01 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 2.2e+02 55  0  0  0  3  55  0  0  0  4     0
MatAssemblyBegin     450 1.0 9.7404e-02 3.6 0.00e+00 0.0 1.4e+03 2.4e+04 9.0e+02  0  0 27 70 14   0  0 27 70 15     0
MatAssemblyEnd       450 1.0 1.0440e-01 1.0 0.00e+00 0.0 4.0e+00 2.2e+03 4.6e+02  0  0  0  0  7   0  0  0  0  7     0
MatGetRowIJ          225 1.0 5.9366e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering       225 1.0 3.0253e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 4.5e+02  0  0  0  0  7   0  0  0  0  7     0
MatZeroEntries       227 1.0 1.2238e-02 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog         1 1.0 5.8174e-05 1.6 7.83e+03 1.1 0.0e+00 0.0e+00 1.0e+00  0  0  0  0  0   0  0  0  0  0   262
KSPSetup             450 1.0 1.1706e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve             225 1.0 9.1850e+01 1.0 4.90e+10 1.1 4.5e+02 8.7e+03 1.1e+03 98100  9  9 17  98100  9  9 18  1002
PCSetUp              450 1.0 9.1803e+01 1.1 4.90e+10 1.1 0.0e+00 0.0e+00 6.8e+02 93100  0  0 10  93100  0  0 11  1001
PCSetUpOnBlocks      225 1.0 9.1802e+01 1.1 4.90e+10 1.1 0.0e+00 0.0e+00 6.8e+02 93100  0  0 10  93100  0  0 11  1001
PCApply                2 1.0 2.7690e-03 1.1 3.33e+06 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2294
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    20             20       236944     0
         Vec Scatter   229            229       198772     0
           Index Set   908            908      5762144     0
   IS L to G Mapping     3              3        27264     0
              Matrix   228            228   2249465464     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 1.19209e-07
Average time for MPI_Barrier(): 8.10623e-07
Average time for zero size MPI_Send(): 5.00679e-06
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
 
***************************************************************
* Done Running Example  mpirun -np 2 ./ex22-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
