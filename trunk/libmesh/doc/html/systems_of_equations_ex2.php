<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("systems_of_equations_ex2",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Systems Example 2 - Unsteady Nonlinear Navier-Stokes</h1>

<br><br>This example shows how a simple, unsteady, nonlinear system of equations
can be solved in parallel.  The system of equations are the familiar
Navier-Stokes equations for low-speed incompressible fluid flow.  This
example introduces the concept of the inner nonlinear loop for each
timestep, and requires a good deal of linear algebra number-crunching
at each step.  If you have a ExodusII viewer such as ParaView installed,
the script movie.sh in this directory will also take appropriate screen
shots of each of the solution files in the time sequence.  These rgb files
can then be animated with the "animate" utility of ImageMagick if it is
installed on your system.  On a PIII 1GHz machine in debug mode, this
example takes a little over a minute to run.  If you would like to see
a more detailed time history, or compute more timesteps, that is certainly
possible by changing the n_timesteps and dt variables below.


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
                                               QUAD4);
        
          mesh.all_second_order();
          
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
          PerfLog perf_log("Systems Example 2");
          
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
Incremenet the time counter, set the time step size as
a parameter in the EquationSystem.
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
                  file_name &lt;&lt; ".e";
                  
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
                      boundary_id_type bc_id = mesh.boundary_info-&gt;boundary_id (elem,s);
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
                
</pre>
</div>
<div class = "comment">
Pin the pressure to zero at global node number "pressure_node".
This effectively removes the non-trivial null space of constant
pressure solutions.
</div>

<div class ="fragment">
<pre>
                const bool pin_pressure = true;
                if (pin_pressure)
                  {
                    const unsigned int pressure_node = 0;
                    const Real p_value               = 0.0;
                    for (unsigned int c=0; c&lt;elem-&gt;n_nodes(); c++)
                      if (elem-&gt;node(c) == pressure_node)
                        {
                          Kpp(c,c) += penalty;
                          Fp(c)    += penalty*p_value;
                        }
                  }
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
                                         QUAD4);
  
    mesh.all_second_order();
    
    mesh.print_info();
    
    EquationSystems equation_systems (mesh);
    
    TransientLinearImplicitSystem &amp; system = 
      equation_systems.add_system&lt;TransientLinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Navier-Stokes&quot;</FONT></B>);
    
    system.add_variable (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, SECOND);
    system.add_variable (<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>, SECOND);
  
    system.add_variable (<B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>, FIRST);
  
    system.attach_assemble_function (assemble_stokes);
    
    equation_systems.init ();
  
    equation_systems.print_info();
  
    PerfLog perf_log(<B><FONT COLOR="#BC8F8F">&quot;Systems Example 2&quot;</FONT></B>);
    
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
            OStringStream file_name;
  
            file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out_&quot;</FONT></B>;
            OSSRealzeroright(file_name,3,0, t_step + 1);
            file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.e&quot;</FONT></B>;
            
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
  
    DenseSubVector&lt;Number&gt;
      Fu(Fe),
      Fv(Fe),
      Fp(Fe);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices_u;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices_v;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices_p;
  
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
                boundary_id_type bc_id = mesh.boundary_info-&gt;boundary_id (elem,s);
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
</FONT></I>          
          <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> pin_pressure = true;
          <B><FONT COLOR="#A020F0">if</FONT></B> (pin_pressure)
            {
              <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> pressure_node = 0;
              <B><FONT COLOR="#228B22">const</FONT></B> Real p_value               = 0.0;
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> c=0; c&lt;elem-&gt;n_nodes(); c++)
                <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;node(c) == pressure_node)
                  {
                    Kpp(c,c) += penalty;
                    Fp(c)    += penalty*p_value;
                  }
            }
        } <I><FONT COLOR="#B22222">// end boundary condition section          
</FONT></I>        
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
  
        navier_stokes_system.matrix-&gt;add_matrix (Ke, dof_indices);
        navier_stokes_system.rhs-&gt;add_vector    (Fe, dof_indices);
      } <I><FONT COLOR="#B22222">// end of element loop
</FONT></I>    
    <B><FONT COLOR="#A020F0">return</FONT></B>;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Linking systems_of_equations_ex2-opt...
***************************************************************
* Running Example  mpirun -np 6 ./systems_of_equations_ex2-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1681
    n_local_nodes()=305
  n_elem()=400
    n_local_elem()=67
    n_active_elem()=400
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Navier-Stokes"
    Type "TransientLinearImplicit"
    Variables="u" "v" "p" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" "CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" "SECOND", "THIRD" "FIRST", "THIRD" 
    n_dofs()=3803
    n_local_dofs()=696
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 36.7155
      Average Off-Processor Bandwidth <= 2.72557
      Maximum  On-Processor Bandwidth <= 59
      Maximum Off-Processor Bandwidth <= 37
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0



*** Solving time step 0, time = 0.01 ***
Linear solver converged at step: 1, final residual: 0.0507375  Nonlinear convergence: ||u - u_old|| = 312.175
Linear solver converged at step: 0, final residual: 0.0507458  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0507458  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0507458  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0507458  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0507458  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0507458  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0507458  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0507458  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0507458  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0507458  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0507458  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0507458  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0507458  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.0507458  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 1, time = 0.02 ***
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 2, time = 0.03 ***
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 3, time = 0.04 ***
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 4, time = 0.05 ***
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 5, time = 0.06 ***
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 6, time = 0.07 ***
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 7, time = 0.08 ***
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 8, time = 0.09 ***
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 9, time = 0.1 ***
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 10, time = 0.11 ***
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 11, time = 0.12 ***
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 12, time = 0.13 ***
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 13, time = 0.14 ***
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0


*** Solving time step 14, time = 0.15 ***
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0
Linear solver converged at step: 0, final residual: 0.051332  Nonlinear convergence: ||u - u_old|| = 0

-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 6                                                |
| Time:           Fri Aug 24 15:27:00 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| Systems Example 2 Performance: Alive time=36.3396, Active time=35.6202                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| linear solve                  225       35.6202     0.158312    35.6202     0.158312    100.00   100.00   |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       225       35.6202                                         100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./systems_of_equations_ex2-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:27:00 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           3.638e+01      1.00195   3.632e+01
Objects:              9.420e+02      1.00000   9.420e+02
Flops:                8.594e+09      3.23269   4.749e+09  2.849e+10
Flops/sec:            2.367e+08      3.23267   1.308e+08  7.846e+08
MPI Messages:         4.269e+03      1.94399   3.192e+03  1.915e+04
MPI Message Lengths:  1.116e+07      1.73569   2.586e+03  4.951e+07
MPI Reductions:       5.959e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 3.6321e+01 100.0%  2.8495e+10 100.0%  1.915e+04 100.0%  2.586e+03      100.0%  5.489e+03  92.1% 

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

VecMDot                1 1.0 8.6808e-0445.5 1.39e+03 1.2 0.0e+00 0.0e+00 1.0e+00  0  0  0  0  0   0  0  0  0  0     9
VecNorm              676 1.0 6.3057e+00 3.1 9.41e+05 1.2 0.0e+00 0.0e+00 6.8e+02 10  0  0  0 11  10  0  0  0 12     1
VecScale             226 1.0 3.7956e-04 1.7 1.57e+05 1.2 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2264
VecCopy              692 1.0 1.2591e-03 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               460 1.0 4.0174e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              676 1.0 1.4250e-02 1.8 9.41e+05 1.2 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   361
VecMAXPY               2 1.0 2.1458e-06 2.2 2.78e+03 1.2 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  7089
VecAssemblyBegin     916 1.0 1.6580e+00 6.4 0.00e+00 0.0 4.0e+03 4.3e+02 2.7e+03  3  0 21  3 46   3  0 21  3 50     0
VecAssemblyEnd       916 1.0 1.4341e-03 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      466 1.0 3.9847e-03 1.5 0.00e+00 0.0 8.9e+03 5.8e+02 0.0e+00  0  0 46 10  0   0  0 46 10  0     0
VecScatterEnd        466 1.0 1.2920e+01 2.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 23  0  0  0  0  23  0  0  0  0     0
VecNormalize         226 1.0 5.4440e+00 5.1 4.72e+05 1.2 0.0e+00 0.0e+00 2.3e+02  8  0  0  0  4   8  0  0  0  4     0
MatMult              226 1.0 1.2897e+01 2.4 1.23e+07 1.2 4.1e+03 4.5e+02 0.0e+00 23  0 21  4  0  23  0 21  4  0     5
MatSolve               2 1.0 1.2450e-03 1.4 8.97e+05 1.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  3728
MatLUFactorNum       225 1.0 8.5741e+00 2.9 8.58e+09 3.2 0.0e+00 0.0e+00 0.0e+00 13100  0  0  0  13100  0  0  0  3314
MatILUFactorSym      225 1.0 1.8219e+01 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 2.2e+02 44  0  0  0  4  44  0  0  0  4     0
MatAssemblyBegin     450 1.0 9.4196e-01 1.4 0.00e+00 0.0 6.1e+03 7.0e+03 9.0e+02  2  0 32 86 15   2  0 32 86 16     0
MatAssemblyEnd       450 1.0 1.0122e-01 1.1 0.00e+00 0.0 3.6e+01 1.2e+02 4.6e+02  0  0  0  0  8   0  0  0  0  8     0
MatGetRowIJ          225 1.0 6.8426e-05 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering       225 1.0 2.7695e-03 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 4.5e+02  0  0  0  0  8   0  0  0  0  8     0
MatZeroEntries       227 1.0 7.5724e-03 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog         1 1.0 8.9788e-0420.0 2.78e+03 1.2 0.0e+00 0.0e+00 1.0e+00  0  0  0  0  0   0  0  0  0  0    17
KSPSetup             450 1.0 1.0324e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve             225 1.0 3.3310e+01 1.0 8.59e+09 3.2 4.1e+03 4.5e+02 1.1e+03 90100 21  4 19  90100 21  4 21   855
PCSetUp              450 1.0 2.4817e+01 1.4 8.58e+09 3.2 0.0e+00 0.0e+00 6.8e+02 57100  0  0 11  57100  0  0 12  1145
PCSetUpOnBlocks      225 1.0 2.4816e+01 1.4 8.58e+09 3.2 0.0e+00 0.0e+00 6.8e+02 57100  0  0 11  57100  0  0 12  1145
PCApply                2 1.0 1.3850e-03 1.3 8.97e+05 1.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  3351
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    20             20       105248     0
         Vec Scatter     4              4         3472     0
           Index Set   683            683      2237520     0
   IS L to G Mapping     3              3        11952     0
              Matrix   228            228    607385620     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 4.03881e-05
Average time for zero size MPI_Send(): 3.9657e-05
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
Configure run at: Sat May 19 03:47:23 2012
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared=1 --with-shared-libraries=1 --with-mpi-dir=/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid --with-mumps=true --download-mumps=1 --with-parmetis=true --download-parmetis=1 --with-superlu=true --download-superlu=1 --with-superludir=true --download-superlu_dist=1 --with-blacs=true --download-blacs=1 --with-scalapack=true --download-scalapack=1 --with-hypre=true --download-hypre=1 --with-blas-lib="[/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_intel_lp64.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_sequential.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_core.so]" --with-lapack-lib=/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_solver_lp64_sequential.a
-----------------------------------------
Libraries compiled on Sat May 19 03:47:23 CDT 2012 on daedalus 
Machine characteristics: Linux daedalus 2.6.32-34-generic #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011 x86_64 GNU/Linux 
Using PETSc directory: /org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5
Using PETSc arch: intel-11.1-lucid-mpich2-1.4.1-cxx-opt
-----------------------------------------
Using C compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpicxx -O3   -fPIC   
Using Fortran compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpif90 -fPIC -O3    
-----------------------------------------
Using include paths: -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/include  
------------------------------------------
Using C linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpicxx -O3 
Using Fortran linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpif90 -fPIC -O3  
Using libraries: -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -lpetsc       -lX11 -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -lHYPRE -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lblacs -lsuperlu_dist_2.4 -lparmetis -lmetis -lsuperlu_4.0 -Wl,-rpath,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t -L/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t -lmkl_solver_lp64_sequential -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -ldl -Wl,-rpath,/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/lib -L/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/lib -lmpich -lopa -lmpl -lrt -lpthread -Wl,-rpath,/opt/intel/Compiler/11.1/073/lib/intel64 -L/opt/intel/Compiler/11.1/073/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -L/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -lmpichf90 -lifport -lifcore -lm -lm -lmpichcxx -lstdc++ -lmpichcxx -lstdc++ -ldl -lmpich -lopa -lmpl -lrt -lpthread -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -ldl  
------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=36.5358, Active time=35.4971                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0001      0.000092    0.0001      0.000129    0.00     0.00     |
|   build_sparsity()                 1         0.0020      0.002001    0.0025      0.002451    0.01     0.01     |
|   create_dof_constraints()         1         0.0001      0.000099    0.0001      0.000099    0.00     0.00     |
|   distribute_dofs()                1         0.0006      0.000581    0.0044      0.004434    0.00     0.01     |
|   dof_indices()                    63756     0.0397      0.000001    0.0397      0.000001    0.11     0.11     |
|   prepare_send_list()              1         0.0000      0.000019    0.0000      0.000019    0.00     0.00     |
|   reinit()                         1         0.0010      0.000979    0.0010      0.000979    0.00     0.00     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          15        0.0044      0.000290    0.0166      0.001109    0.01     0.05     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               15        0.0151      0.001007    0.0151      0.001007    0.04     0.04     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        30150     0.0839      0.000003    0.0839      0.000003    0.24     0.24     |
|   init_shape_functions()           450       0.0027      0.000006    0.0027      0.000006    0.01     0.01     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             30150     0.0373      0.000001    0.0373      0.000001    0.11     0.11     |
|   init_reference_to_physical_map() 450       0.0060      0.000013    0.0060      0.000013    0.02     0.02     |
|                                                                                                                |
| Mesh                                                                                                           |
|   all_second_order()               1         0.0026      0.002610    0.0026      0.002610    0.01     0.01     |
|   find_neighbors()                 2         0.0006      0.000307    0.0008      0.000420    0.00     0.00     |
|   renumber_nodes_and_elem()        4         0.0001      0.000037    0.0001      0.000037    0.00     0.00     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        3         0.0050      0.001652    0.0050      0.001652    0.01     0.01     |
|   find_global_indices()            3         0.0004      0.000143    0.0100      0.003334    0.00     0.03     |
|   parallel_sort()                  3         0.0014      0.000452    0.0043      0.001436    0.00     0.01     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         15        0.0002      0.000011    0.0319      0.002128    0.00     0.09     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0002      0.000151    0.0002      0.000151    0.00     0.00     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      2         0.0028      0.001416    0.0095      0.004742    0.01     0.03     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      10        0.0026      0.000259    0.0026      0.000259    0.01     0.01     |
|   broadcast()                      1         0.0000      0.000012    0.0000      0.000012    0.00     0.00     |
|   gather()                         1         0.0000      0.000005    0.0000      0.000005    0.00     0.00     |
|   max(bool)                        1         0.0000      0.000014    0.0000      0.000014    0.00     0.00     |
|   max(scalar)                      3         0.0003      0.000092    0.0003      0.000092    0.00     0.00     |
|   max(vector)                      3         0.0001      0.000028    0.0001      0.000028    0.00     0.00     |
|   min(vector)                      3         0.0001      0.000049    0.0001      0.000049    0.00     0.00     |
|   probe()                          70        0.0041      0.000059    0.0041      0.000059    0.01     0.01     |
|   receive()                        70        0.0001      0.000002    0.0042      0.000060    0.00     0.01     |
|   send()                           70        0.0001      0.000001    0.0001      0.000001    0.00     0.00     |
|   send_receive()                   76        0.0001      0.000002    0.0045      0.000059    0.00     0.01     |
|   sum()                            55        0.0139      0.000252    0.0139      0.000252    0.04     0.04     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           70        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         2         0.0002      0.000106    0.0037      0.001857    0.00     0.01     |
|   set_parent_processor_ids()       2         0.0001      0.000025    0.0001      0.000025    0.00     0.00     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          225       34.4071     0.152921    34.4071     0.152921    96.93    96.93    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       225       0.8621      0.003832    1.0387      0.004616    2.43     2.93     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            125913    35.4971                                         100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 6 ./systems_of_equations_ex2-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
