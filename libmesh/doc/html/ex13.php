<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex13",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 13 - Unsteady Navier-Stokes Equations - Unsteady Nonlinear Systems of Equations</h1>

<br><br>This example shows how a simple, unsteady, nonlinear system of equations
can be solved in parallel.  The system of equations are the familiar
Navier-Stokes equations for low-speed incompressible fluid flow.  This
example introduces the concept of the inner nonlinear loop for each
timestep, and requires a good deal of linear algebra number-crunching
at each step.  If you have the General Mesh Viewer (GMV) installed,
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
        #include "gmv_io.h"
        #include "equation_systems.h"
        #include "fe.h"
        #include "quadrature_gauss.h"
        #include "dof_map.h"
        #include "sparse_matrix.h"
        #include "numeric_vector.h"
        #include "dense_matrix.h"
        #include "dense_vector.h"
        #include "implicit_system.h"
        #include "error_vector.h"
        #include "error_estimator.h"
        #include "transient_system.h"
        #include "perf_log.h"
        
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
Use the internal mesh generator to create a uniform
grid on the square [-1,1]^D.  We instruct the mesh generator
to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
elements in 3D.  Building these higher-order elements allows
us to use higher-order approximation, as in example 3.
</div>

<div class ="fragment">
<pre>
            mesh.build_square (20, 20,
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
</div>

<div class ="fragment">
<pre>
            {
</pre>
</div>
<div class = "comment">
Creates a transient system named "Navier-Stokes"
</div>

<div class ="fragment">
<pre>
              TransientImplicitSystem & system = 
        	equation_systems.add_system&lt;TransientImplicitSystem&gt; ("Navier-Stokes");
              
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
        
              equation_systems.set_parameter("linear solver maximum iterations") = 250;
              equation_systems.set_parameter("linear solver tolerance") = 1.e-3;
              
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
            const Real dt = 0.005;
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
            const Real nonlinear_tolerance       = 1.e-2;
            
</pre>
</div>
<div class = "comment">
Tell the system of equations what the timestep is by using
the set_parameter function.  The matrix assembly routine can
then reference this parameter.
</div>

<div class ="fragment">
<pre>
            equation_systems.set_parameter ("dt")   = dt;
        
</pre>
</div>
<div class = "comment">
Get a reference to the Stokes system to use later.
</div>

<div class ="fragment">
<pre>
            TransientImplicitSystem&  stokes_system =
        	  equation_systems.get_system&lt;TransientImplicitSystem&gt;("Navier-Stokes");
        
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
              last_nonlinear_soln (stokes_system.solution-&gt;clone());
        
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
                equation_systems.set_parameter ("time") = time;
        
</pre>
</div>
<div class = "comment">
A pretty update message
</div>

<div class ="fragment">
<pre>
                std::cout &lt;&lt; " Solving time step " &lt;&lt; t_step &lt;&lt; ", time = " &lt;&lt; time &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Now we need to update the solution vector from the
previous time step.  This is done directly through
the reference to the Stokes system.
</div>

<div class ="fragment">
<pre>
                *stokes_system.old_local_solution = *stokes_system.solution;
        
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
        	    last_nonlinear_soln-&gt;add(*stokes_system.solution);
        	    
</pre>
</div>
<div class = "comment">
Assemble & solve the linear system.
</div>

<div class ="fragment">
<pre>
                    perf_log.start_event("linear solve");
        	    equation_systems("Navier-Stokes").solve();
        	    perf_log.stop_event("linear solve");
        
</pre>
</div>
<div class = "comment">
Compute the difference between this solution and the last
nonlinear iterate.
</div>

<div class ="fragment">
<pre>
                    last_nonlinear_soln-&gt;add (-1., *stokes_system.solution);
        
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
Print out convergence information for the linear and
nonlinear iterations.
</div>

<div class ="fragment">
<pre>
                    std::cout &lt;&lt; "Linear solver converged at step: "
        		      &lt;&lt; stokes_system.n_linear_iterations()
        		      &lt;&lt; ", final residual: "
        		      &lt;&lt; stokes_system.final_linear_residual()
        		      &lt;&lt; "  Nonlinear convergence: ||u - u_old|| = "
        		      &lt;&lt; norm_delta
        		      &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Terminate the solution iteration if the difference between
this iteration and the last is sufficiently small.
</div>

<div class ="fragment">
<pre>
                    if (norm_delta &lt; nonlinear_tolerance)
        	      {
        		std::cout &lt;&lt; " Nonlinear solver converged at step "
        			  &lt;&lt; l
        			  &lt;&lt; std::endl;
        		break;
        	      }
        	    
        	  } // end nonlinear loop
        	
</pre>
</div>
<div class = "comment">
Write out every nth timestep to file.
</div>

<div class ="fragment">
<pre>
                const unsigned int write_interval = 1;
        	
        	if ((t_step+1)%write_interval == 0)
        	  {
        	    OStringStream file_name;
        
</pre>
</div>
<div class = "comment">
We write the file name in the gmv auto-read format.
</div>

<div class ="fragment">
<pre>
                    file_name &lt;&lt; "out.gmv.";
        	    OSSRealzeroright(file_name,3,0, t_step + 1);
        	    
        	    GMVIO(mesh).write_equation_systems (file_name.str(),
        						equation_systems);
        	  }
              } // end timestep loop.
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
          assert (system_name == "Navier-Stokes");
          
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
Get a reference to the Stokes system object.
</div>

<div class ="fragment">
<pre>
          TransientImplicitSystem & stokes_system =
            es.get_system&lt;TransientImplicitSystem&gt; ("Navier-Stokes");
        
</pre>
</div>
<div class = "comment">
Numeric ids corresponding to each variable in the system
</div>

<div class ="fragment">
<pre>
          const unsigned int u_var = stokes_system.variable_number ("u");
          const unsigned int v_var = stokes_system.variable_number ("v");
          const unsigned int p_var = stokes_system.variable_number ("p");
          
</pre>
</div>
<div class = "comment">
Get the Finite Element type for "u".  Note this will be
the same as the type for "v".
</div>

<div class ="fragment">
<pre>
          FEType fe_vel_type = stokes_system.variable_type(u_var);
          
</pre>
</div>
<div class = "comment">
Get the Finite Element type for "p".
</div>

<div class ="fragment">
<pre>
          FEType fe_pres_type = stokes_system.variable_type(p_var);
        
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
          const DofMap & dof_map = stokes_system.get_dof_map();
        
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
          const Real dt    = es.parameter("dt");
</pre>
</div>
<div class = "comment">
const Real time  = es.parameter("time");
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
          const_active_local_elem_iterator           el (mesh.elements_begin());
          const const_active_local_elem_iterator end_el (mesh.elements_end());
          
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
| Kpu Kpv Kpp |        | Fv |
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
                      u_old += phi[l][qp]*stokes_system.old_solution (dof_indices_u[l]);
        	      v_old += phi[l][qp]*stokes_system.old_solution (dof_indices_v[l]);
        	      grad_u_old.add_scaled (dphi[l][qp],stokes_system.old_solution (dof_indices_u[l]));
        	      grad_v_old.add_scaled (dphi[l][qp],stokes_system.old_solution (dof_indices_v[l]));
        
</pre>
</div>
<div class = "comment">
From the previous Newton iterate:
</div>

<div class ="fragment">
<pre>
                      u += phi[l][qp]*stokes_system.current_solution (dof_indices_u[l]); 
        	      v += phi[l][qp]*stokes_system.current_solution (dof_indices_v[l]);
        	      grad_u.add_scaled (dphi[l][qp],stokes_system.current_solution (dof_indices_u[l]));
        	      grad_v.add_scaled (dphi[l][qp],stokes_system.current_solution (dof_indices_v[l]));
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
        	      p_old += psi[l][qp]*stokes_system.old_solution (dof_indices_p[l]);
        	    }
        
</pre>
</div>
<div class = "comment">
Definitions for convenience.  It is sometimes simpler to do a
dot product if you have the full vector at your disposal.
</div>

<div class ="fragment">
<pre>
                  const Point U_old (u_old, v_old);
        	  const Point U     (u,     v);
        	  const Real  u_x = grad_u(0);
        	  const Real  u_y = grad_u(1);
        	  const Real  v_x = grad_v(0);
        	  const Real  v_y = grad_v(1);
        	  
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
The location on the boundary of the current
node.
		   

<br><br>const Real xf = side->point(ns)(0);
</div>

<div class ="fragment">
<pre>
                          const Real yf = side-&gt;point(ns)(1);
        		  
</pre>
</div>
<div class = "comment">
The penalty value.  \f$ \frac{1}{\epsilon \f$
</div>

<div class ="fragment">
<pre>
                          const Real penalty = 1.e10;
        		  
</pre>
</div>
<div class = "comment">
The boundary values.
		   

<br><br>Set u = 1 on the top boundary, 0 everywhere else
</div>

<div class ="fragment">
<pre>
                          const Real u_value = (yf &gt; .99) ? 1. : 0.;
        		  
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
              stokes_system.matrix-&gt;add_matrix (Ke, dof_indices);
              stokes_system.rhs-&gt;add_vector    (Fe, dof_indices);
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
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;libmesh.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;mesh.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;gmv_io.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;equation_systems.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;fe.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dof_map.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;sparse_matrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;numeric_vector.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_matrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_vector.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;implicit_system.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;error_vector.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;error_estimator.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;transient_system.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;perf_log.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;o_string_stream.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_submatrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_subvector.h&quot;</FONT></B>
  
  <FONT COLOR="#228B22"><B>void</FONT></B> assemble_stokes (EquationSystems&amp; es,
  		      <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name);
  
  <FONT COLOR="#228B22"><B>int</FONT></B> main (<FONT COLOR="#228B22"><B>int</FONT></B> argc, <FONT COLOR="#228B22"><B>char</FONT></B>** argv)
  {
    libMesh::init (argc, argv);
    {    
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = 2;     
      
      Mesh mesh (dim);
      
      mesh.build_square (20, 20,
  		       0., 1.,
  		       0., 1.,
  		       QUAD9);
      
      mesh.print_info();
      
      EquationSystems equation_systems (mesh);
      
      {
        TransientImplicitSystem &amp; system = 
  	equation_systems.add_system&lt;TransientImplicitSystem&gt; (<FONT COLOR="#BC8F8F"><B>&quot;Navier-Stokes&quot;</FONT></B>);
        
        system.add_variable (<FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>, SECOND);
        system.add_variable (<FONT COLOR="#BC8F8F"><B>&quot;v&quot;</FONT></B>, SECOND);
  
        system.add_variable (<FONT COLOR="#BC8F8F"><B>&quot;p&quot;</FONT></B>, FIRST);
  
        system.attach_assemble_function (assemble_stokes);
        
        equation_systems.init ();
  
        equation_systems.set_parameter(<FONT COLOR="#BC8F8F"><B>&quot;linear solver maximum iterations&quot;</FONT></B>) = 250;
        equation_systems.set_parameter(<FONT COLOR="#BC8F8F"><B>&quot;linear solver tolerance&quot;</FONT></B>) = 1.e-3;
        
        equation_systems.print_info();
      }
  
      PerfLog perf_log(<FONT COLOR="#BC8F8F"><B>&quot;Example 13&quot;</FONT></B>);
      
      <FONT COLOR="#228B22"><B>const</FONT></B> Real dt = 0.005;
      Real time     = 0.0;
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_timesteps = 15;
  
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_nonlinear_steps = 15;
      <FONT COLOR="#228B22"><B>const</FONT></B> Real nonlinear_tolerance       = 1.e-2;
      
      equation_systems.set_parameter (<FONT COLOR="#BC8F8F"><B>&quot;dt&quot;</FONT></B>)   = dt;
  
      TransientImplicitSystem&amp;  stokes_system =
  	  equation_systems.get_system&lt;TransientImplicitSystem&gt;(<FONT COLOR="#BC8F8F"><B>&quot;Navier-Stokes&quot;</FONT></B>);
  
      AutoPtr&lt;NumericVector&lt;Number&gt; &gt;
        last_nonlinear_soln (stokes_system.solution-&gt;clone());
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> t_step=0; t_step&lt;n_timesteps; ++t_step)
        {
  	time += dt;
  
  	equation_systems.set_parameter (<FONT COLOR="#BC8F8F"><B>&quot;time&quot;</FONT></B>) = time;
  
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; Solving time step &quot;</FONT></B> &lt;&lt; t_step &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;, time = &quot;</FONT></B> &lt;&lt; time &lt;&lt; std::endl;
  
  	*stokes_system.old_local_solution = *stokes_system.solution;
  
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> l=0; l&lt;n_nonlinear_steps; ++l)
  	  {
  	    last_nonlinear_soln-&gt;zero();
  	    last_nonlinear_soln-&gt;add(*stokes_system.solution);
  	    
  	    perf_log.start_event(<FONT COLOR="#BC8F8F"><B>&quot;linear solve&quot;</FONT></B>);
  	    equation_systems(<FONT COLOR="#BC8F8F"><B>&quot;Navier-Stokes&quot;</FONT></B>).solve();
  	    perf_log.stop_event(<FONT COLOR="#BC8F8F"><B>&quot;linear solve&quot;</FONT></B>);
  
  	    last_nonlinear_soln-&gt;add (-1., *stokes_system.solution);
  
  	    last_nonlinear_soln-&gt;close();
  
  	    <FONT COLOR="#228B22"><B>const</FONT></B> Real norm_delta = last_nonlinear_soln-&gt;l2_norm();
  
  	    std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Linear solver converged at step: &quot;</FONT></B>
  		      &lt;&lt; stokes_system.n_linear_iterations()
  		      &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;, final residual: &quot;</FONT></B>
  		      &lt;&lt; stokes_system.final_linear_residual()
  		      &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;  Nonlinear convergence: ||u - u_old|| = &quot;</FONT></B>
  		      &lt;&lt; norm_delta
  		      &lt;&lt; std::endl;
  
  	    <B><FONT COLOR="#A020F0">if</FONT></B> (norm_delta &lt; nonlinear_tolerance)
  	      {
  		std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; Nonlinear solver converged at step &quot;</FONT></B>
  			  &lt;&lt; l
  			  &lt;&lt; std::endl;
  		<B><FONT COLOR="#A020F0">break</FONT></B>;
  	      }
  	    
  	  } <I><FONT COLOR="#B22222">// end nonlinear loop
</FONT></I>  	
  	<FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> write_interval = 1;
  	
  	<B><FONT COLOR="#A020F0">if</FONT></B> ((t_step+1)%write_interval == 0)
  	  {
  	    OStringStream file_name;
  
  	    file_name &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;out.gmv.&quot;</FONT></B>;
  	    OSSRealzeroright(file_name,3,0, t_step + 1);
  	    
  	    GMVIO(mesh).write_equation_systems (file_name.str(),
  						equation_systems);
  	  }
        } <I><FONT COLOR="#B22222">// end timestep loop.
</FONT></I>    }
    
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close ();
  }
  
  
  
  
  
  
  <FONT COLOR="#228B22"><B>void</FONT></B> assemble_stokes (EquationSystems&amp; es,
  		      <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name)
  {
    assert (system_name == <FONT COLOR="#BC8F8F"><B>&quot;Navier-Stokes&quot;</FONT></B>);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> Mesh&amp; mesh = es.get_mesh();
    
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = mesh.mesh_dimension();
    
    TransientImplicitSystem &amp; stokes_system =
      es.get_system&lt;TransientImplicitSystem&gt; (<FONT COLOR="#BC8F8F"><B>&quot;Navier-Stokes&quot;</FONT></B>);
  
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> u_var = stokes_system.variable_number (<FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>);
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> v_var = stokes_system.variable_number (<FONT COLOR="#BC8F8F"><B>&quot;v&quot;</FONT></B>);
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> p_var = stokes_system.variable_number (<FONT COLOR="#BC8F8F"><B>&quot;p&quot;</FONT></B>);
    
    FEType fe_vel_type = stokes_system.variable_type(u_var);
    
    FEType fe_pres_type = stokes_system.variable_type(p_var);
  
    AutoPtr&lt;FEBase&gt; fe_vel  (FEBase::build(dim, fe_vel_type));
      
    AutoPtr&lt;FEBase&gt; fe_pres (FEBase::build(dim, fe_pres_type));
    
    QGauss qrule (dim, fe_vel_type.default_quadrature_order());
  
    fe_vel-&gt;attach_quadrature_rule (&amp;qrule);
    fe_pres-&gt;attach_quadrature_rule (&amp;qrule);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe_vel-&gt;get_JxW();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe_vel-&gt;get_phi();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe_vel-&gt;get_dphi();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; psi = fe_pres-&gt;get_phi();
  
    
    <FONT COLOR="#228B22"><B>const</FONT></B> DofMap &amp; dof_map = stokes_system.get_dof_map();
  
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
  
    std::vector&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt; dof_indices;
    std::vector&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt; dof_indices_u;
    std::vector&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt; dof_indices_v;
    std::vector&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt; dof_indices_p;
  
    <FONT COLOR="#228B22"><B>const</FONT></B> Real dt    = es.parameter(<FONT COLOR="#BC8F8F"><B>&quot;dt&quot;</FONT></B>);
    <FONT COLOR="#228B22"><B>const</FONT></B> Real theta = 1.;
      
    const_active_local_elem_iterator           el (mesh.elements_begin());
    <FONT COLOR="#228B22"><B>const</FONT></B> const_active_local_elem_iterator end_el (mesh.elements_end());
    
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {    
        <FONT COLOR="#228B22"><B>const</FONT></B> Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        dof_map.dof_indices (elem, dof_indices_u, u_var);
        dof_map.dof_indices (elem, dof_indices_v, v_var);
        dof_map.dof_indices (elem, dof_indices_p, p_var);
  
        <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_dofs   = dof_indices.size();
        <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_u_dofs = dof_indices_u.size(); 
        <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_v_dofs = dof_indices_v.size(); 
        <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_p_dofs = dof_indices_p.size();
        
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
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
  	{
  	  Number   u = 0., u_old = 0.;
  	  Number   v = 0., v_old = 0.;
  	  Number   p_old = 0.;
  	  Gradient grad_u, grad_u_old;
  	  Gradient grad_v, grad_v_old;
  	  
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> l=0; l&lt;n_u_dofs; l++)
  	    {
  	      u_old += phi[l][qp]*stokes_system.old_solution (dof_indices_u[l]);
  	      v_old += phi[l][qp]*stokes_system.old_solution (dof_indices_v[l]);
  	      grad_u_old.add_scaled (dphi[l][qp],stokes_system.old_solution (dof_indices_u[l]));
  	      grad_v_old.add_scaled (dphi[l][qp],stokes_system.old_solution (dof_indices_v[l]));
  
  	      u += phi[l][qp]*stokes_system.current_solution (dof_indices_u[l]); 
  	      v += phi[l][qp]*stokes_system.current_solution (dof_indices_v[l]);
  	      grad_u.add_scaled (dphi[l][qp],stokes_system.current_solution (dof_indices_u[l]));
  	      grad_v.add_scaled (dphi[l][qp],stokes_system.current_solution (dof_indices_v[l]));
  	    }
  
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> l=0; l&lt;n_p_dofs; l++)
  	    {
  	      p_old += psi[l][qp]*stokes_system.old_solution (dof_indices_p[l]);
  	    }
  
  	  <FONT COLOR="#228B22"><B>const</FONT></B> Point U_old (u_old, v_old);
  	  <FONT COLOR="#228B22"><B>const</FONT></B> Point U     (u,     v);
  	  <FONT COLOR="#228B22"><B>const</FONT></B> Real  u_x = grad_u(0);
  	  <FONT COLOR="#228B22"><B>const</FONT></B> Real  u_y = grad_u(1);
  	  <FONT COLOR="#228B22"><B>const</FONT></B> Real  v_x = grad_v(0);
  	  <FONT COLOR="#228B22"><B>const</FONT></B> Real  v_y = grad_v(1);
  	  
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;n_u_dofs; i++)
  	    {
  	      Fu(i) += JxW[qp]*(u_old*phi[i][qp] -                            <I><FONT COLOR="#B22222">// mass-matrix term 
</FONT></I>  				(1.-theta)*dt*(U_old*grad_u_old)*phi[i][qp] + <I><FONT COLOR="#B22222">// convection term
</FONT></I>  				(1.-theta)*dt*p_old*dphi[i][qp](0)  -         <I><FONT COLOR="#B22222">// pressure term on rhs
</FONT></I>  				(1.-theta)*dt*(grad_u_old*dphi[i][qp]) +      <I><FONT COLOR="#B22222">// diffusion term on rhs
</FONT></I>  				theta*dt*(U*grad_u)*phi[i][qp]);              <I><FONT COLOR="#B22222">// Newton term
</FONT></I>  
  		
  	      Fv(i) += JxW[qp]*(v_old*phi[i][qp] -                             <I><FONT COLOR="#B22222">// mass-matrix term
</FONT></I>  				(1.-theta)*dt*(U_old*grad_v_old)*phi[i][qp] +  <I><FONT COLOR="#B22222">// convection term
</FONT></I>  				(1.-theta)*dt*p_old*dphi[i][qp](1) -           <I><FONT COLOR="#B22222">// pressure term on rhs
</FONT></I>  				(1.-theta)*dt*(grad_v_old*dphi[i][qp]) +       <I><FONT COLOR="#B22222">// diffusion term on rhs
</FONT></I>  				theta*dt*(U*grad_v)*phi[i][qp]);               <I><FONT COLOR="#B22222">// Newton term
</FONT></I>  	    
  
  
  	      <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j&lt;n_u_dofs; j++)
  		{
  		  Kuu(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp] +                <I><FONT COLOR="#B22222">// mass matrix term
</FONT></I>  				       theta*dt*(dphi[i][qp]*dphi[j][qp]) +   <I><FONT COLOR="#B22222">// diffusion term
</FONT></I>  				       theta*dt*(U*dphi[j][qp])*phi[i][qp] +  <I><FONT COLOR="#B22222">// convection term
</FONT></I>  				       theta*dt*u_x*phi[i][qp]*phi[j][qp]);   <I><FONT COLOR="#B22222">// Newton term
</FONT></I>  
  		  Kuv(i,j) += JxW[qp]*theta*dt*u_y*phi[i][qp]*phi[j][qp];     <I><FONT COLOR="#B22222">// Newton term
</FONT></I>  		  
  		  Kvv(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp] +                <I><FONT COLOR="#B22222">// mass matrix term
</FONT></I>  				       theta*dt*(dphi[i][qp]*dphi[j][qp]) +   <I><FONT COLOR="#B22222">// diffusion term
</FONT></I>  				       theta*dt*(U*dphi[j][qp])*phi[i][qp] +  <I><FONT COLOR="#B22222">// convection term
</FONT></I>  				       theta*dt*v_y*phi[i][qp]*phi[j][qp]);   <I><FONT COLOR="#B22222">// Newton term
</FONT></I>  
  		  Kvu(i,j) += JxW[qp]*theta*dt*v_x*phi[i][qp]*phi[j][qp];     <I><FONT COLOR="#B22222">// Newton term
</FONT></I>  		}
  
  	      <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j&lt;n_p_dofs; j++)
  		{
  		  Kup(i,j) += JxW[qp]*(-theta*dt*psi[j][qp]*dphi[i][qp](0));
  		  Kvp(i,j) += JxW[qp]*(-theta*dt*psi[j][qp]*dphi[i][qp](1));
  		}
  	    }
  
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;n_p_dofs; i++)
  	    {
  	      <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j&lt;n_u_dofs; j++)
  		{
  		  Kpu(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](0);
  		  Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1);
  		}
  	    }
  	} <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>  
        
        {
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> s=0; s&lt;elem-&gt;n_sides(); s++)
  	  <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(s) == NULL)
  	    {
  	      AutoPtr&lt;Elem&gt; side (elem-&gt;build_side(s));
  	      	      
  	      <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> ns=0; ns&lt;side-&gt;n_nodes(); ns++)
  		{
  		   
  		  <FONT COLOR="#228B22"><B>const</FONT></B> Real yf = side-&gt;point(ns)(1);
  		  
  		  <FONT COLOR="#228B22"><B>const</FONT></B> Real penalty = 1.e10;
  		  
  		   
  		  <FONT COLOR="#228B22"><B>const</FONT></B> Real u_value = (yf &gt; .99) ? 1. : 0.;
  		  
  		  <FONT COLOR="#228B22"><B>const</FONT></B> Real v_value = 0.;
  		  
  		  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n=0; n&lt;elem-&gt;n_nodes(); n++)
  		    <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;node(n) == side-&gt;node(ns))
  		      {
  			Kuu(n,n) += penalty;
  			Kvv(n,n) += penalty;
  		  		  
  			Fu(n) += penalty*u_value;
  			Fv(n) += penalty*v_value;
  		      }
  		} <I><FONT COLOR="#B22222">// end face node loop	  
</FONT></I>  	    } <I><FONT COLOR="#B22222">// end if (elem-&gt;neighbor(side) == NULL)
</FONT></I>        } <I><FONT COLOR="#B22222">// end boundary condition section	  
</FONT></I>        
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        
        stokes_system.matrix-&gt;add_matrix (Ke, dof_indices);
        stokes_system.rhs-&gt;add_vector    (Fe, dof_indices);
      } <I><FONT COLOR="#B22222">// end of element loop
</FONT></I>    
    <B><FONT COLOR="#A020F0">return</FONT></B>;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example  ./ex13
***************************************************************
 
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1681
  n_elem()=400
   n_local_elem()=400
   n_active_elem()=400
  n_subdomains()=1
  n_processors()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Navier-Stokes"
    Type "TransientImplicit"
    Variables="u" "v" "p" 
    Finite Element Types="0" "0" "0" 
    Approximation Orders="2" "2" "1" 
    n_dofs()=3803
    n_local_dofs()=3803
    n_constrained_dofs()=0
    n_vectors()=1
  n_parameters()=2
   Parameters:
    "linear solver maximum iterations"=250
    "linear solver tolerance"=0.001

 Solving time step 0, time = 0.005
Linear solver converged at step: 16, final residual: 0.121329  Nonlinear convergence: ||u - u_old|| = 223.471
Linear solver converged at step: 24, final residual: 0.000317615  Nonlinear convergence: ||u - u_old|| = 6.5317
Linear solver converged at step: 28, final residual: 2.58051e-07  Nonlinear convergence: ||u - u_old|| = 0.0072933
 Nonlinear solver converged at step 2
 Solving time step 1, time = 0.01
Linear solver converged at step: 23, final residual: 0.0025053  Nonlinear convergence: ||u - u_old|| = 35.4701
Linear solver converged at step: 24, final residual: 1.47282e-05  Nonlinear convergence: ||u - u_old|| = 0.0466908
Linear solver converged at step: 24, final residual: 1.27369e-08  Nonlinear convergence: ||u - u_old|| = 0.00027294
 Nonlinear solver converged at step 2
 Solving time step 2, time = 0.015
Linear solver converged at step: 24, final residual: 0.000790085  Nonlinear convergence: ||u - u_old|| = 9.79414
Linear solver converged at step: 24, final residual: 3.41055e-06  Nonlinear convergence: ||u - u_old|| = 0.0208939
Linear solver converged at step: 28, final residual: 2.96023e-09  Nonlinear convergence: ||u - u_old|| = 0.000135564
 Nonlinear solver converged at step 2
 Solving time step 3, time = 0.02
Linear solver converged at step: 24, final residual: 0.000382917  Nonlinear convergence: ||u - u_old|| = 4.07167
Linear solver converged at step: 25, final residual: 1.42401e-06  Nonlinear convergence: ||u - u_old|| = 0.00672638
 Nonlinear solver converged at step 1
 Solving time step 4, time = 0.025
Linear solver converged at step: 23, final residual: 0.000247016  Nonlinear convergence: ||u - u_old|| = 2.13266
Linear solver converged at step: 25, final residual: 8.20034e-07  Nonlinear convergence: ||u - u_old|| = 0.00772513
 Nonlinear solver converged at step 1
 Solving time step 5, time = 0.03
Linear solver converged at step: 23, final residual: 0.000126994  Nonlinear convergence: ||u - u_old|| = 1.26206
Linear solver converged at step: 27, final residual: 4.97111e-07  Nonlinear convergence: ||u - u_old|| = 0.0057915
 Nonlinear solver converged at step 1
 Solving time step 6, time = 0.035
Linear solver converged at step: 24, final residual: 8.03056e-05  Nonlinear convergence: ||u - u_old|| = 0.805059
Linear solver converged at step: 24, final residual: 2.93058e-07  Nonlinear convergence: ||u - u_old|| = 0.00147397
 Nonlinear solver converged at step 1
 Solving time step 7, time = 0.04
Linear solver converged at step: 24, final residual: 7.94782e-05  Nonlinear convergence: ||u - u_old|| = 0.53619
Linear solver converged at step: 23, final residual: 1.75705e-07  Nonlinear convergence: ||u - u_old|| = 0.000919154
 Nonlinear solver converged at step 1
 Solving time step 8, time = 0.045
Linear solver converged at step: 23, final residual: 6.92094e-05  Nonlinear convergence: ||u - u_old|| = 0.368994
Linear solver converged at step: 27, final residual: 9.13022e-08  Nonlinear convergence: ||u - u_old|| = 0.00104694
 Nonlinear solver converged at step 1
 Solving time step 9, time = 0.05
Linear solver converged at step: 23, final residual: 4.73186e-05  Nonlinear convergence: ||u - u_old|| = 0.260159
Linear solver converged at step: 27, final residual: 6.18452e-08  Nonlinear convergence: ||u - u_old|| = 0.000762091
 Nonlinear solver converged at step 1
 Solving time step 10, time = 0.055
Linear solver converged at step: 22, final residual: 4.01872e-05  Nonlinear convergence: ||u - u_old|| = 0.187318
Linear solver converged at step: 27, final residual: 3.72042e-08  Nonlinear convergence: ||u - u_old|| = 0.000711503
 Nonlinear solver converged at step 1
 Solving time step 11, time = 0.06
Linear solver converged at step: 22, final residual: 2.7971e-05  Nonlinear convergence: ||u - u_old|| = 0.137235
Linear solver converged at step: 28, final residual: 2.40094e-08  Nonlinear convergence: ||u - u_old|| = 0.000617341
 Nonlinear solver converged at step 1
 Solving time step 12, time = 0.065
Linear solver converged at step: 21, final residual: 2.51357e-05  Nonlinear convergence: ||u - u_old|| = 0.102129
Linear solver converged at step: 28, final residual: 2.11058e-08  Nonlinear convergence: ||u - u_old|| = 0.00086522
 Nonlinear solver converged at step 1
 Solving time step 13, time = 0.07
Linear solver converged at step: 21, final residual: 1.78319e-05  Nonlinear convergence: ||u - u_old|| = 0.07677
Linear solver converged at step: 28, final residual: 1.60542e-08  Nonlinear convergence: ||u - u_old|| = 0.000734319
 Nonlinear solver converged at step 1
 Solving time step 14, time = 0.075
Linear solver converged at step: 21, final residual: 1.27677e-05  Nonlinear convergence: ||u - u_old|| = 0.0582848
Linear solver converged at step: 28, final residual: 1.17971e-08  Nonlinear convergence: ||u - u_old|| = 0.000574086
 Nonlinear solver converged at step 1

 ----------------------------------------------------------------------------
| Time:           Mon Apr 19 12:06:05 2004
| OS:             Linux
| HostName:       arthur
| OS Release      2.4.20-19.9smp
| OS Version:     #1 SMP Tue Jul 15 17:04:18 EDT 2003
| Machine:        i686
| Username:       peterson
 ----------------------------------------------------------------------------
 ----------------------------------------------------------------------------
| Example 13 Performance: Alive time=53.0525, Active time=50.2831
 ----------------------------------------------------------------------------
| Event                         nCalls  Total       Avg         Percent of   |
|                                       Time        Time        Active Time  |
|----------------------------------------------------------------------------|
|                                                                            |
| linear solve                  33      50.2831     1.523732    100.00       |
 ----------------------------------------------------------------------------
| Totals:                       33      50.2831                 100.00       |
 ----------------------------------------------------------------------------


 ---------------------------------------------------------------------------- 
| Reference count information                                                |
 ---------------------------------------------------------------------------- 
| 12SparseMatrixIdE reference count information:
| Creations:    1
| Destructions: 1
| 13NumericVectorIdE reference count information:
| Creations:    6
| Destructions: 6
| 21LinearSolverInterfaceIdE reference count information:
| Creations:    1
| Destructions: 1
| 4Elem reference count information:
| Creations:    4640
| Destructions: 4640
| 4Node reference count information:
| Creations:    1681
| Destructions: 1681
| 5QBase reference count information:
| Creations:    66
| Destructions: 66
| 6DofMap reference count information:
| Creations:    1
| Destructions: 1
| 6FEBase reference count information:
| Creations:    66
| Destructions: 66
| 6System reference count information:
| Creations:    1
| Destructions: 1
 ---------------------------------------------------------------------------- 
 
***************************************************************
* Done Running Example  ./ex13
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
