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
          PerfLog perf_log("Example 13");
          
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
  
    PerfLog perf_log(<B><FONT COLOR="#BC8F8F">&quot;Example 13&quot;</FONT></B>);
    
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
Compiling C++ (in optimized mode) systems_of_equations_ex2.C...
Linking systems_of_equations_ex2-opt...
***************************************************************
* Running Example  ./systems_of_equations_ex2-opt
***************************************************************
 
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1681
    n_local_nodes()=1681
  n_elem()=400
    n_local_elem()=400
    n_active_elem()=400
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
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
    n_local_dofs()=3803
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 38.9716
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 59
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0



*** Solving time step 0, time = 0.01 ***
Linear solver converged at step: 236, final residual: 0.000214866  Nonlinear convergence: ||u - u_old|| = 282.453
Linear solver converged at step: 90, final residual: 9.99752e-06  Nonlinear convergence: ||u - u_old|| = 0.743588
Linear solver converged at step: 250, final residual: 1.04092e-07  Nonlinear convergence: ||u - u_old|| = 0.0151698
Linear solver converged at step: 250, final residual: 6.12561e-10  Nonlinear convergence: ||u - u_old|| = 0.000178243
 Nonlinear solver converged at step 3


*** Solving time step 1, time = 0.02 ***
Linear solver converged at step: 148, final residual: 0.000196444  Nonlinear convergence: ||u - u_old|| = 30.6968
Linear solver converged at step: 74, final residual: 8.13684e-06  Nonlinear convergence: ||u - u_old|| = 0.194999
Linear solver converged at step: 250, final residual: 3.57979e-08  Nonlinear convergence: ||u - u_old|| = 0.00887051
Linear solver converged at step: 250, final residual: 1.73324e-10  Nonlinear convergence: ||u - u_old|| = 5.87181e-05
 Nonlinear solver converged at step 3


*** Solving time step 2, time = 0.03 ***
Linear solver converged at step: 77, final residual: 0.000223335  Nonlinear convergence: ||u - u_old|| = 6.00665
Linear solver converged at step: 118, final residual: 9.3062e-06  Nonlinear convergence: ||u - u_old|| = 0.268525
Linear solver converged at step: 250, final residual: 3.68148e-08  Nonlinear convergence: ||u - u_old|| = 0.00315838
Linear solver converged at step: 250, final residual: 1.04561e-09  Nonlinear convergence: ||u - u_old|| = 6.4103e-05
 Nonlinear solver converged at step 3


*** Solving time step 3, time = 0.04 ***
Linear solver converged at step: 72, final residual: 0.000210493  Nonlinear convergence: ||u - u_old|| = 2.18249
Linear solver converged at step: 159, final residual: 9.54006e-06  Nonlinear convergence: ||u - u_old|| = 0.260477
Linear solver converged at step: 250, final residual: 1.75501e-07  Nonlinear convergence: ||u - u_old|| = 0.0104015
Linear solver converged at step: 250, final residual: 3.10351e-09  Nonlinear convergence: ||u - u_old|| = 0.000301456
 Nonlinear solver converged at step 3


*** Solving time step 4, time = 0.05 ***
Linear solver converged at step: 58, final residual: 0.000192963  Nonlinear convergence: ||u - u_old|| = 1.02052
Linear solver converged at step: 119, final residual: 7.33059e-06  Nonlinear convergence: ||u - u_old|| = 0.103956
Linear solver converged at step: 225, final residual: 1.16642e-08  Nonlinear convergence: ||u - u_old|| = 0.00534664
Linear solver converged at step: 250, final residual: 1.57996e-12  Nonlinear convergence: ||u - u_old|| = 3.60887e-06
 Nonlinear solver converged at step 3


*** Solving time step 5, time = 0.06 ***
Linear solver converged at step: 52, final residual: 0.000214167  Nonlinear convergence: ||u - u_old|| = 0.476324
Linear solver converged at step: 172, final residual: 1.00119e-05  Nonlinear convergence: ||u - u_old|| = 0.202403
Linear solver converged at step: 158, final residual: 2.16558e-08  Nonlinear convergence: ||u - u_old|| = 0.0173926
Linear solver converged at step: 250, final residual: 2.42714e-11  Nonlinear convergence: ||u - u_old|| = 2.11938e-05
 Nonlinear solver converged at step 3


*** Solving time step 6, time = 0.07 ***
Linear solver converged at step: 21, final residual: 0.000225026  Nonlinear convergence: ||u - u_old|| = 0.256854
Linear solver converged at step: 117, final residual: 1.12946e-05  Nonlinear convergence: ||u - u_old|| = 0.147309
Linear solver converged at step: 250, final residual: 1.54083e-06  Nonlinear convergence: ||u - u_old|| = 0.0128016
Linear solver converged at step: 250, final residual: 2.37644e-09  Nonlinear convergence: ||u - u_old|| = 0.00239445
Linear solver converged at step: 250, final residual: 2.6763e-11  Nonlinear convergence: ||u - u_old|| = 3.10462e-06
 Nonlinear solver converged at step 4


*** Solving time step 7, time = 0.08 ***
Linear solver converged at step: 19, final residual: 0.000192228  Nonlinear convergence: ||u - u_old|| = 0.152223
Linear solver converged at step: 67, final residual: 7.02714e-06  Nonlinear convergence: ||u - u_old|| = 0.0850414
Linear solver converged at step: 194, final residual: 1.08542e-08  Nonlinear convergence: ||u - u_old|| = 0.00357158
Linear solver converged at step: 250, final residual: 3.18135e-12  Nonlinear convergence: ||u - u_old|| = 1.67286e-05
 Nonlinear solver converged at step 3


*** Solving time step 8, time = 0.09 ***
Linear solver converged at step: 15, final residual: 0.000204876  Nonlinear convergence: ||u - u_old|| = 0.086678
Linear solver converged at step: 97, final residual: 8.76746e-06  Nonlinear convergence: ||u - u_old|| = 0.0454467
Linear solver converged at step: 137, final residual: 1.7204e-08  Nonlinear convergence: ||u - u_old|| = 0.00314359
Linear solver converged at step: 250, final residual: 3.57292e-12  Nonlinear convergence: ||u - u_old|| = 1.96891e-05
 Nonlinear solver converged at step 3


*** Solving time step 9, time = 0.1 ***
Linear solver converged at step: 11, final residual: 0.00022001  Nonlinear convergence: ||u - u_old|| = 0.0505868
Linear solver converged at step: 69, final residual: 1.08054e-05  Nonlinear convergence: ||u - u_old|| = 0.0182682
Linear solver converged at step: 250, final residual: 3.49668e-08  Nonlinear convergence: ||u - u_old|| = 0.0127024
Linear solver converged at step: 250, final residual: 1.35967e-10  Nonlinear convergence: ||u - u_old|| = 5.5716e-05
 Nonlinear solver converged at step 3


*** Solving time step 10, time = 0.11 ***
Linear solver converged at step: 11, final residual: 0.000139573  Nonlinear convergence: ||u - u_old|| = 0.0322373
Linear solver converged at step: 104, final residual: 4.3327e-06  Nonlinear convergence: ||u - u_old|| = 0.0126849
Linear solver converged at step: 250, final residual: 3.63694e-08  Nonlinear convergence: ||u - u_old|| = 0.00565483
Linear solver converged at step: 250, final residual: 3.17858e-11  Nonlinear convergence: ||u - u_old|| = 5.24272e-05
 Nonlinear solver converged at step 3


*** Solving time step 11, time = 0.12 ***
Linear solver converged at step: 10, final residual: 0.000145884  Nonlinear convergence: ||u - u_old|| = 0.0205125
Linear solver converged at step: 48, final residual: 4.75038e-06  Nonlinear convergence: ||u - u_old|| = 0.0078103
Linear solver converged at step: 250, final residual: 1.24897e-08  Nonlinear convergence: ||u - u_old|| = 0.00689682
Linear solver converged at step: 250, final residual: 3.2272e-11  Nonlinear convergence: ||u - u_old|| = 2.14012e-05
 Nonlinear solver converged at step 3


*** Solving time step 12, time = 0.13 ***
Linear solver converged at step: 10, final residual: 9.34715e-05  Nonlinear convergence: ||u - u_old|| = 0.0136016
Linear solver converged at step: 68, final residual: 1.95662e-06  Nonlinear convergence: ||u - u_old|| = 0.00578812
Linear solver converged at step: 250, final residual: 4.80563e-09  Nonlinear convergence: ||u - u_old|| = 0.00253896
Linear solver converged at step: 250, final residual: 6.17056e-11  Nonlinear convergence: ||u - u_old|| = 9.37577e-06
 Nonlinear solver converged at step 3


*** Solving time step 13, time = 0.14 ***
Linear solver converged at step: 9, final residual: 0.000149994  Nonlinear convergence: ||u - u_old|| = 0.00910172
Linear solver converged at step: 19, final residual: 3.88794e-06  Nonlinear convergence: ||u - u_old|| = 0.00478718
Linear solver converged at step: 250, final residual: 1.08721e-08  Nonlinear convergence: ||u - u_old|| = 0.00216698
Linear solver converged at step: 250, final residual: 1.2336e-09  Nonlinear convergence: ||u - u_old|| = 1.52877e-05
 Nonlinear solver converged at step 3


*** Solving time step 14, time = 0.15 ***
Linear solver converged at step: 8, final residual: 0.00015955  Nonlinear convergence: ||u - u_old|| = 0.00594443
Linear solver converged at step: 17, final residual: 5.39679e-06  Nonlinear convergence: ||u - u_old|| = 0.00369955
Linear solver converged at step: 250, final residual: 4.63355e-08  Nonlinear convergence: ||u - u_old|| = 0.00112717
Linear solver converged at step: 250, final residual: 2.63014e-11  Nonlinear convergence: ||u - u_old|| = 7.54831e-05
 Nonlinear solver converged at step 3

-------------------------------------------------------------------
| Time:           Sat Apr  7 16:04:00 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| Example 13 Performance: Alive time=20.1578, Active time=20.074                                            |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| linear solve                  61        20.0740     0.329082    20.0740     0.329082    100.00   100.00   |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       61        20.0740                                         100.00            |
 -----------------------------------------------------------------------------------------------------------

 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=20.2554, Active time=20.1144                                               |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 1         0.0003      0.000263    0.0003      0.000263    0.00     0.00     |
|   compute_sparsity()           1         0.0098      0.009779    0.0121      0.012069    0.05     0.06     |
|   create_dof_constraints()     1         0.0003      0.000257    0.0003      0.000257    0.00     0.00     |
|   distribute_dofs()            1         0.0012      0.001156    0.0029      0.002933    0.01     0.01     |
|   dof_indices()                116000    0.1205      0.000001    0.1205      0.000001    0.60     0.60     |
|   prepare_send_list()          1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                     1         0.0018      0.001776    0.0018      0.001776    0.01     0.01     |
|                                                                                                            |
| EquationSystems                                                                                            |
|   build_solution_vector()      15        0.0227      0.001512    0.0465      0.003100    0.11     0.23     |
|                                                                                                            |
| ExodusII_IO                                                                                                |
|   write_nodal_data()           15        0.0276      0.001837    0.0276      0.001837    0.14     0.14     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         48800     0.0803      0.000002    0.0803      0.000002    0.40     0.40     |
|   compute_shape_functions()    48800     0.0337      0.000001    0.0337      0.000001    0.17     0.17     |
|   init_shape_functions()       122       0.0047      0.000038    0.0047      0.000038    0.02     0.02     |
|                                                                                                            |
| Mesh                                                                                                       |
|   all_second_order()           1         0.0038      0.003832    0.0038      0.003832    0.02     0.02     |
|   find_neighbors()             2         0.0029      0.001464    0.0029      0.001464    0.01     0.01     |
|   renumber_nodes_and_elem()    4         0.0003      0.000081    0.0003      0.000081    0.00     0.00     |
|                                                                                                            |
| MeshOutput                                                                                                 |
|   write_equation_systems()     15        0.0003      0.000022    0.0744      0.004960    0.00     0.37     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0008      0.000822    0.0008      0.000822    0.00     0.00     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  1         0.0000      0.000000    0.0000      0.000000    0.00     0.00     |
|                                                                                                            |
| Partitioner                                                                                                |
|   single_partition()           2         0.0002      0.000108    0.0002      0.000108    0.00     0.00     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      61        17.9257     0.293864    17.9257     0.293864    89.12    89.12    |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   61        1.8776      0.030781    2.1255      0.034845    9.33     10.57    |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        213906    20.1144                                         100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./systems_of_equations_ex2-opt
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
