<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex18",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 18 - Unsteady Navier-Stokes Equations with DiffSystem</h1>

<br><br>This example shows how the transient nonlinear problem from
example 13 can be solved using the DiffSystem class framework to
simplify the user-implemented equations.


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
Basic include files
</div>

<div class ="fragment">
<pre>
        #include "equation_systems.h"
        #include "fe_base.h"
        #include "gmv_io.h"
        #include "libmesh.h"
        #include "mesh.h"
        #include "mesh_generation.h"
        #include "quadrature.h"
        
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
DiffSystem framework files
</div>

<div class ="fragment">
<pre>
        #include "diff_solver.h"
        #include "fem_system.h"
        #include "euler_solver.h"
        
</pre>
</div>
<div class = "comment">
The Navier-Stokes system class.
FEMSystem, TimeSolver and  NewtonSolver will handle most tasks,
but we must specify element residuals


<br><br>The interval between our timesteps
</div>

<div class ="fragment">
<pre>
        const Real deltat = 0.005;
        
</pre>
</div>
<div class = "comment">
And the number of timesteps to take
</div>

<div class ="fragment">
<pre>
        const unsigned int n_timesteps = 15;
        
        class NavierSystem : public FEMSystem
        {
        public:
</pre>
</div>
<div class = "comment">
Constructor
</div>

<div class ="fragment">
<pre>
          NavierSystem(EquationSystems& es,
                       const std::string& name,
                       const unsigned int number)
          : FEMSystem(es, name, number), Reynolds(1.) {}
        
</pre>
</div>
<div class = "comment">
System initialization
</div>

<div class ="fragment">
<pre>
          virtual void init_data ();
        
</pre>
</div>
<div class = "comment">
Element residual and jacobian calculations
Time dependent parts
</div>

<div class ="fragment">
<pre>
          virtual bool element_time_derivative (bool request_jacobian);
        
</pre>
</div>
<div class = "comment">
Constraint parts
</div>

<div class ="fragment">
<pre>
          virtual bool element_constraint (bool request_jacobian);
          virtual bool side_constraint (bool request_jacobian);
        
</pre>
</div>
<div class = "comment">
Finite elements for the velocity and pressure on element interiors
</div>

<div class ="fragment">
<pre>
          FEBase *fe_velocity, *fe_pressure;
        
</pre>
</div>
<div class = "comment">
Finite element for the velocity on element sides
</div>

<div class ="fragment">
<pre>
          FEBase *fe_side_vel;
        
</pre>
</div>
<div class = "comment">
The Reynolds number to solve for
</div>

<div class ="fragment">
<pre>
          Real Reynolds;
        };
        
        
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
Create a two-dimensional mesh.
</div>

<div class ="fragment">
<pre>
            Mesh mesh (2);
            
</pre>
</div>
<div class = "comment">
Use the MeshTools::Generation mesh generator to create a uniform
grid on the square [-1,1]^D.  We instruct the mesh generator
to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
elements in 3D.  Building these higher-order elements allows
us to use higher-order approximation, as in example 3.
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
Declare the system "Navier-Stokes" and its variables.
</div>

<div class ="fragment">
<pre>
            NavierSystem & stokes_system = 
              equation_systems.add_system&lt;NavierSystem&gt; ("Navier-Stokes");
        
</pre>
</div>
<div class = "comment">
Solve this as a time-dependent system
</div>

<div class ="fragment">
<pre>
            stokes_system.time_solver =
               AutoPtr&lt;TimeSolver&gt;(new EulerSolver(stokes_system));
        
</pre>
</div>
<div class = "comment">
Initialize the system
</div>

<div class ="fragment">
<pre>
            equation_systems.init ();
        
</pre>
</div>
<div class = "comment">
Set the time stepping options
</div>

<div class ="fragment">
<pre>
            stokes_system.deltat = deltat;
        
</pre>
</div>
<div class = "comment">
And the nonlinear solver options
</div>

<div class ="fragment">
<pre>
            DiffSolver &solver = *stokes_system.time_solver-&gt;diff_solver;
            solver.max_nonlinear_iterations = 15;
            solver.relative_step_tolerance = 0.02;
        
</pre>
</div>
<div class = "comment">
And the linear solver options
</div>

<div class ="fragment">
<pre>
            solver.max_linear_iterations = 250;
            solver.initial_linear_tolerance = std::sqrt(TOLERANCE);
        
</pre>
</div>
<div class = "comment">
Print information about the system to the screen.
</div>

<div class ="fragment">
<pre>
            equation_systems.print_info();
        
</pre>
</div>
<div class = "comment">
Now we begin the timestep loop to compute the time-accurate
solution of the equations.
</div>

<div class ="fragment">
<pre>
            for (unsigned int t_step=0; t_step != n_timesteps; ++t_step)
              {
</pre>
</div>
<div class = "comment">
A pretty update message
</div>

<div class ="fragment">
<pre>
                std::cout &lt;&lt; " Solving time step " &lt;&lt; t_step &lt;&lt; ", time = "
                          &lt;&lt; stokes_system.time &lt;&lt; std::endl;
        
                stokes_system.solve();
        
                stokes_system.time_solver-&gt;advance_timestep();
        
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
              }
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
        
        
        
        void NavierSystem::init_data ()
        {
</pre>
</div>
<div class = "comment">
Add the velocity components "u" & "v".  They
will be approximated using second-order approximation.
</div>

<div class ="fragment">
<pre>
          this-&gt;add_variable ("u", SECOND);
          this-&gt;add_variable ("v", SECOND);
        
</pre>
</div>
<div class = "comment">
Add the pressure variable "p". This will
be approximated with a first-order basis,
providing an LBB-stable pressure-velocity pair.
</div>

<div class ="fragment">
<pre>
          this-&gt;add_variable ("p", FIRST);
        
</pre>
</div>
<div class = "comment">
Do the parent's initialization after variables are defined
</div>

<div class ="fragment">
<pre>
          FEMSystem::init_data();
        
</pre>
</div>
<div class = "comment">
Tell the system to march u and v forward in time, but 
leave p as a constraint only
</div>

<div class ="fragment">
<pre>
          this-&gt;time_evolving(0);
          this-&gt;time_evolving(1);
        
</pre>
</div>
<div class = "comment">
Get references to the finite elements we need
</div>

<div class ="fragment">
<pre>
          fe_velocity = element_fe[this-&gt;variable_type(0)];
          fe_pressure = element_fe[this-&gt;variable_type(2)];
          fe_side_vel = side_fe[this-&gt;variable_type(0)];
        
</pre>
</div>
<div class = "comment">
Now make sure we have requested all the data
we need to build the linear system.
</div>

<div class ="fragment">
<pre>
          fe_velocity-&gt;get_JxW();
          fe_velocity-&gt;get_phi();
          fe_velocity-&gt;get_dphi();
        
          fe_pressure-&gt;get_phi();
        
          fe_side_vel-&gt;get_JxW();
          fe_side_vel-&gt;get_phi();
          fe_side_vel-&gt;get_xyz();
        }
        
        
        bool NavierSystem::element_time_derivative (bool request_jacobian)
        {
</pre>
</div>
<div class = "comment">
First we get some references to cell-specific data that
will be used to assemble the linear system.


<br><br>Element Jacobian * quadrature weights for interior integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW = fe_velocity-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The velocity shape functions at interior quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt; &phi = fe_velocity-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The velocity shape function gradients at interior
quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt; &dphi =
            fe_velocity-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
The pressure shape functions at interior
quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt; &psi = fe_pressure-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = dof_indices_var[0].size(); 
          assert (n_u_dofs == dof_indices_var[1].size()); 
          const unsigned int n_p_dofs = dof_indices_var[2].size();
        
</pre>
</div>
<div class = "comment">
The subvectors and submatrices we need to fill:
</div>

<div class ="fragment">
<pre>
          DenseSubMatrix&lt;Number&gt; &Kuu = *elem_subjacobians[0][0];
          DenseSubMatrix&lt;Number&gt; &Kvv = *elem_subjacobians[1][1];
          DenseSubMatrix&lt;Number&gt; &Kuv = *elem_subjacobians[0][1];
          DenseSubMatrix&lt;Number&gt; &Kvu = *elem_subjacobians[1][0];
          DenseSubMatrix&lt;Number&gt; &Kup = *elem_subjacobians[0][2];
          DenseSubMatrix&lt;Number&gt; &Kvp = *elem_subjacobians[1][2];
          DenseSubVector&lt;Number&gt; &Fu = *elem_subresiduals[0];
          DenseSubVector&lt;Number&gt; &Fv = *elem_subresiduals[1];
              
</pre>
</div>
<div class = "comment">
Now we will build the element Jacobian and residual.
Constructing the residual requires the solution and its
gradient from the previous timestep.  This must be
calculated at each quadrature point by summing the
solution degree-of-freedom values by the appropriate
weight functions.
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = element_qrule-&gt;n_points();
        
          for (unsigned int qp=0; qp != n_qpoints; qp++)
            {
</pre>
</div>
<div class = "comment">
Compute the solution & its gradient at the old Newton iterate
</div>

<div class ="fragment">
<pre>
              Number u = interior_value(0, qp),
                     v = interior_value(1, qp),
                     p = interior_value(2, qp);
              Gradient grad_u = interior_gradient(0, qp),
                       grad_v = interior_gradient(1, qp);
        
</pre>
</div>
<div class = "comment">
Definitions for convenience.  It is sometimes simpler to do a
dot product if you have the full vector at your disposal.
</div>

<div class ="fragment">
<pre>
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
              for (unsigned int i=0; i != n_u_dofs; i++)
                {
                  Fu(i) += JxW[qp] *
                           (-Reynolds*(U*grad_u)*phi[i][qp] + // convection term
                            p*dphi[i][qp](0) -                // pressure term
                            (grad_u*dphi[i][qp]));            // diffusion term
                    
                  Fv(i) += JxW[qp] *
                           (-Reynolds*(U*grad_v)*phi[i][qp] + // convection term
                            p*dphi[i][qp](1) -                // pressure term
                            (grad_v*dphi[i][qp]));            // diffusion term
        
</pre>
</div>
<div class = "comment">
Note that the Fp block is identically zero unless we are using
some kind of artificial compressibility scheme...


<br><br>Matrix contributions for the uu and vv couplings.
</div>

<div class ="fragment">
<pre>
                  if (request_jacobian)
                    for (unsigned int j=0; j != n_u_dofs; j++)
                      {
                        Kuu(i,j) += JxW[qp] *
         /* convection term */      (-Reynolds*(U*dphi[j][qp])*phi[i][qp] -
         /* diffusion term  */       (dphi[i][qp]*dphi[j][qp]) -
         /* Newton term     */       Reynolds*u_x*phi[i][qp]*phi[j][qp]);
        
                        Kuv(i,j) += JxW[qp] *
         /* Newton term     */      -Reynolds*u_y*phi[i][qp]*phi[j][qp];
        
                        Kvv(i,j) += JxW[qp] *
         /* convection term */      (-Reynolds*(U*dphi[j][qp])*phi[i][qp] -
         /* diffusion term  */       (dphi[i][qp]*dphi[j][qp]) -
         /* Newton term     */       Reynolds*v_y*phi[i][qp]*phi[j][qp]);
        
                        Kvu(i,j) += JxW[qp] * 
         /* Newton term     */      -Reynolds*v_x*phi[i][qp]*phi[j][qp];
                    }
        
</pre>
</div>
<div class = "comment">
Matrix contributions for the up and vp couplings.
</div>

<div class ="fragment">
<pre>
                  if (request_jacobian)
                    for (unsigned int j=0; j != n_p_dofs; j++)
                      {
                        Kup(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](0);
                        Kvp(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](1);
                      }
                }
            } // end of the quadrature point qp-loop
          
          return request_jacobian;
        }
        
        
        
        bool NavierSystem::element_constraint (bool request_jacobian)
        {
</pre>
</div>
<div class = "comment">
Here we define some references to cell-specific data that
will be used to assemble the linear system.


<br><br>Element Jacobian * quadrature weight for interior integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW = fe_velocity-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The velocity shape function gradients at interior
quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt; &dphi =
            fe_velocity-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
The pressure shape functions at interior
quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt; &psi = fe_pressure-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in each variable
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = dof_indices_var[0].size();
          const unsigned int n_p_dofs = dof_indices_var[2].size();
        
</pre>
</div>
<div class = "comment">
The subvectors and submatrices we need to fill:
</div>

<div class ="fragment">
<pre>
          DenseSubMatrix&lt;Number&gt; &Kpu = *elem_subjacobians[2][0];
          DenseSubMatrix&lt;Number&gt; &Kpv = *elem_subjacobians[2][1];
          DenseSubVector&lt;Number&gt; &Fp = *elem_subresiduals[2];
        
</pre>
</div>
<div class = "comment">
Add the constraint given by the continuity equation
</div>

<div class ="fragment">
<pre>
          unsigned int n_qpoints = element_qrule-&gt;n_points();
          for (unsigned int qp=0; qp != n_qpoints; qp++)
            {
</pre>
</div>
<div class = "comment">
Compute the velocity gradient at the old Newton iterate
</div>

<div class ="fragment">
<pre>
              Gradient grad_u = interior_gradient(0, qp),
                       grad_v = interior_gradient(1, qp);
        
</pre>
</div>
<div class = "comment">
Now a loop over the pressure degrees of freedom.  This
computes the contributions of the continuity equation.
</div>

<div class ="fragment">
<pre>
              for (unsigned int i=0; i != n_p_dofs; i++)
                {
                  Fp(i) += JxW[qp] * psi[i][qp] *
                           (grad_u(0) + grad_v(1));
        
                  if (request_jacobian)
                    for (unsigned int j=0; j != n_u_dofs; j++)
                      {
                        Kpu(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](0);
                        Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1);
                      }
                }
            } // end of the quadrature point qp-loop
        
          return request_jacobian;
        }
        
              
        
        bool NavierSystem::side_constraint (bool request_jacobian)
        {
</pre>
</div>
<div class = "comment">
Here we define some references to cell-specific data that
will be used to assemble the linear system.


<br><br>Element Jacobian * quadrature weight for side integration
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt; &JxW_side = fe_side_vel-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The velocity shape functions at side quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt; &phi_side =
            fe_side_vel-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The XYZ locations (in physical space) of the
quadrature points on the side.  This is where
we will interpolate the boundary value function.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Point &gt; &qside_point = fe_side_vel-&gt;get_xyz();
        
</pre>
</div>
<div class = "comment">
The number of local degrees of freedom in u and v
</div>

<div class ="fragment">
<pre>
          const unsigned int n_u_dofs = dof_indices_var[0].size(); 
        
</pre>
</div>
<div class = "comment">
The subvectors and submatrices we need to fill:
</div>

<div class ="fragment">
<pre>
          DenseSubMatrix&lt;Number&gt; &Kuu = *elem_subjacobians[0][0];
          DenseSubMatrix&lt;Number&gt; &Kvv = *elem_subjacobians[1][1];
        
          DenseSubVector&lt;Number&gt; &Fu = *elem_subresiduals[0];
          DenseSubVector&lt;Number&gt; &Fv = *elem_subresiduals[1];
        
</pre>
</div>
<div class = "comment">
At this point the interior element integration has
been completed.  However, we have not yet addressed
boundary conditions.  For this example we will only
consider simple Dirichlet boundary conditions imposed
at each timestep via the penalty method.


<br><br>The penalty value.  \f$ \frac{1}{\epsilon} \f$
</div>

<div class ="fragment">
<pre>
          const Real penalty = 1.e10;
        
          unsigned int n_sidepoints = side_qrule-&gt;n_points();
          for (unsigned int qp=0; qp != n_sidepoints; qp++)
            {
</pre>
</div>
<div class = "comment">
Compute the solution at the old Newton iterate
</div>

<div class ="fragment">
<pre>
              Number u = side_value(0, qp),
                     v = side_value(1, qp);
        
</pre>
</div>
<div class = "comment">
The location of the current boundary quadrature point
const Real xf = qside_point[qp](0);
</div>

<div class ="fragment">
<pre>
              const Real yf = qside_point[qp](1);
        
</pre>
</div>
<div class = "comment">
Set u = 1 on the top boundary, 0 everywhere else
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
        
              for (unsigned int i=0; i != n_u_dofs; i++)
                {
                  Fu(i) += JxW_side[qp] * penalty *
                           (u - u_value) * phi_side[i][qp];
                  Fv(i) += JxW_side[qp] * penalty *
                           (v - v_value) * phi_side[i][qp];
        
                  if (request_jacobian)
                    for (unsigned int j=0; j != n_u_dofs; j++)
                      {
                        Kuu(i,j) += JxW_side[qp] * penalty *
                                    phi_side[i][qp] * phi_side[j][qp];
                        Kvv(i,j) += JxW_side[qp] * penalty *
                                    phi_side[i][qp] * phi_side[j][qp];
                      }
                }
            }
        
          return request_jacobian;
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;equation_systems.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;fe_base.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;gmv_io.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;libmesh.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;mesh.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;mesh_generation.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;quadrature.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;o_string_stream.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;diff_solver.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;fem_system.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;euler_solver.h&quot;</FONT></B>
  
  
  <FONT COLOR="#228B22"><B>const</FONT></B> Real deltat = 0.005;
  
  <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_timesteps = 15;
  
  <FONT COLOR="#228B22"><B>class</FONT></B> NavierSystem : <FONT COLOR="#228B22"><B>public</FONT></B> FEMSystem
  {
  <FONT COLOR="#228B22"><B>public</FONT></B>:
    NavierSystem(EquationSystems&amp; es,
                 <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; name,
                 <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> number)
    : FEMSystem(es, name, number), Reynolds(1.) {}
  
    <FONT COLOR="#228B22"><B>virtual</FONT></B> <FONT COLOR="#228B22"><B>void</FONT></B> init_data ();
  
    <FONT COLOR="#228B22"><B>virtual</FONT></B> <FONT COLOR="#228B22"><B>bool</FONT></B> element_time_derivative (<FONT COLOR="#228B22"><B>bool</FONT></B> request_jacobian);
  
    <FONT COLOR="#228B22"><B>virtual</FONT></B> <FONT COLOR="#228B22"><B>bool</FONT></B> element_constraint (<FONT COLOR="#228B22"><B>bool</FONT></B> request_jacobian);
    <FONT COLOR="#228B22"><B>virtual</FONT></B> <FONT COLOR="#228B22"><B>bool</FONT></B> side_constraint (<FONT COLOR="#228B22"><B>bool</FONT></B> request_jacobian);
  
    FEBase *fe_velocity, *fe_pressure;
  
    FEBase *fe_side_vel;
  
    Real Reynolds;
  };
  
  
  <FONT COLOR="#228B22"><B>int</FONT></B> main (<FONT COLOR="#228B22"><B>int</FONT></B> argc, <FONT COLOR="#228B22"><B>char</FONT></B>** argv)
  {
    libMesh::init (argc, argv);
    {    
      Mesh mesh (2);
      
      MeshTools::Generation::build_square (mesh,
                                           20, 20,
                                           0., 1.,
                                           0., 1.,
                                           QUAD9);
      
      mesh.print_info();
      
      EquationSystems equation_systems (mesh);
      
      NavierSystem &amp; stokes_system = 
        equation_systems.add_system&lt;NavierSystem&gt; (<FONT COLOR="#BC8F8F"><B>&quot;Navier-Stokes&quot;</FONT></B>);
  
      stokes_system.time_solver =
         AutoPtr&lt;TimeSolver&gt;(<B><FONT COLOR="#A020F0">new</FONT></B> EulerSolver(stokes_system));
  
      equation_systems.init ();
  
      stokes_system.deltat = deltat;
  
      DiffSolver &amp;solver = *stokes_system.time_solver-&gt;diff_solver;
      solver.max_nonlinear_iterations = 15;
      solver.relative_step_tolerance = 0.02;
  
      solver.max_linear_iterations = 250;
      solver.initial_linear_tolerance = std::sqrt(TOLERANCE);
  
      equation_systems.print_info();
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> t_step=0; t_step != n_timesteps; ++t_step)
        {
          std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; Solving time step &quot;</FONT></B> &lt;&lt; t_step &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;, time = &quot;</FONT></B>
                    &lt;&lt; stokes_system.time &lt;&lt; std::endl;
  
          stokes_system.solve();
  
          stokes_system.time_solver-&gt;advance_timestep();
  
          <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> write_interval = 1;
  
          <B><FONT COLOR="#A020F0">if</FONT></B> ((t_step+1)%write_interval == 0)
            {
              OStringStream file_name;
  
              file_name &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;out.gmv.&quot;</FONT></B>;
              OSSRealzeroright(file_name,3,0, t_step + 1);
  
              GMVIO(mesh).write_equation_systems (file_name.str(),
                                                  equation_systems);
            }
        }
    }
    
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close ();
  }
  
  
  
  <FONT COLOR="#228B22"><B>void</FONT></B> NavierSystem::init_data ()
  {
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>, SECOND);
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<FONT COLOR="#BC8F8F"><B>&quot;v&quot;</FONT></B>, SECOND);
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;add_variable (<FONT COLOR="#BC8F8F"><B>&quot;p&quot;</FONT></B>, FIRST);
  
    FEMSystem::init_data();
  
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;time_evolving(0);
    <B><FONT COLOR="#A020F0">this</FONT></B>-&gt;time_evolving(1);
  
    fe_velocity = element_fe[<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;variable_type(0)];
    fe_pressure = element_fe[<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;variable_type(2)];
    fe_side_vel = side_fe[<B><FONT COLOR="#A020F0">this</FONT></B>-&gt;variable_type(0)];
  
    fe_velocity-&gt;get_JxW();
    fe_velocity-&gt;get_phi();
    fe_velocity-&gt;get_dphi();
  
    fe_pressure-&gt;get_phi();
  
    fe_side_vel-&gt;get_JxW();
    fe_side_vel-&gt;get_phi();
    fe_side_vel-&gt;get_xyz();
  }
  
  
  <FONT COLOR="#228B22"><B>bool</FONT></B> NavierSystem::element_time_derivative (<FONT COLOR="#228B22"><B>bool</FONT></B> request_jacobian)
  {
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = fe_velocity-&gt;get_JxW();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt; &amp;phi = fe_velocity-&gt;get_phi();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt; &amp;dphi =
      fe_velocity-&gt;get_dphi();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt; &amp;psi = fe_pressure-&gt;get_phi();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_u_dofs = dof_indices_var[0].size(); 
    assert (n_u_dofs == dof_indices_var[1].size()); 
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_p_dofs = dof_indices_var[2].size();
  
    DenseSubMatrix&lt;Number&gt; &amp;Kuu = *elem_subjacobians[0][0];
    DenseSubMatrix&lt;Number&gt; &amp;Kvv = *elem_subjacobians[1][1];
    DenseSubMatrix&lt;Number&gt; &amp;Kuv = *elem_subjacobians[0][1];
    DenseSubMatrix&lt;Number&gt; &amp;Kvu = *elem_subjacobians[1][0];
    DenseSubMatrix&lt;Number&gt; &amp;Kup = *elem_subjacobians[0][2];
    DenseSubMatrix&lt;Number&gt; &amp;Kvp = *elem_subjacobians[1][2];
    DenseSubVector&lt;Number&gt; &amp;Fu = *elem_subresiduals[0];
    DenseSubVector&lt;Number&gt; &amp;Fv = *elem_subresiduals[1];
        
    <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_qpoints = element_qrule-&gt;n_points();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp != n_qpoints; qp++)
      {
        Number u = interior_value(0, qp),
               v = interior_value(1, qp),
               p = interior_value(2, qp);
        Gradient grad_u = interior_gradient(0, qp),
                 grad_v = interior_gradient(1, qp);
  
        <FONT COLOR="#228B22"><B>const</FONT></B> NumberVectorValue U     (u,     v);
        <FONT COLOR="#228B22"><B>const</FONT></B> Number  u_x = grad_u(0);
        <FONT COLOR="#228B22"><B>const</FONT></B> Number  u_y = grad_u(1);
        <FONT COLOR="#228B22"><B>const</FONT></B> Number  v_x = grad_v(0);
        <FONT COLOR="#228B22"><B>const</FONT></B> Number  v_y = grad_v(1);
            
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i != n_u_dofs; i++)
          {
            Fu(i) += JxW[qp] *
                     (-Reynolds*(U*grad_u)*phi[i][qp] + <I><FONT COLOR="#B22222">// convection term
</FONT></I>                      p*dphi[i][qp](0) -                <I><FONT COLOR="#B22222">// pressure term
</FONT></I>                      (grad_u*dphi[i][qp]));            <I><FONT COLOR="#B22222">// diffusion term
</FONT></I>              
            Fv(i) += JxW[qp] *
                     (-Reynolds*(U*grad_v)*phi[i][qp] + <I><FONT COLOR="#B22222">// convection term
</FONT></I>                      p*dphi[i][qp](1) -                <I><FONT COLOR="#B22222">// pressure term
</FONT></I>                      (grad_v*dphi[i][qp]));            <I><FONT COLOR="#B22222">// diffusion term
</FONT></I>  
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (request_jacobian)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j != n_u_dofs; j++)
                {
                  Kuu(i,j) += JxW[qp] *
   <I><FONT COLOR="#B22222">/* convection term */</FONT></I>      (-Reynolds*(U*dphi[j][qp])*phi[i][qp] -
   <I><FONT COLOR="#B22222">/* diffusion term  */</FONT></I>       (dphi[i][qp]*dphi[j][qp]) -
   <I><FONT COLOR="#B22222">/* Newton term     */</FONT></I>       Reynolds*u_x*phi[i][qp]*phi[j][qp]);
  
                  Kuv(i,j) += JxW[qp] *
   <I><FONT COLOR="#B22222">/* Newton term     */</FONT></I>      -Reynolds*u_y*phi[i][qp]*phi[j][qp];
  
                  Kvv(i,j) += JxW[qp] *
   <I><FONT COLOR="#B22222">/* convection term */</FONT></I>      (-Reynolds*(U*dphi[j][qp])*phi[i][qp] -
   <I><FONT COLOR="#B22222">/* diffusion term  */</FONT></I>       (dphi[i][qp]*dphi[j][qp]) -
   <I><FONT COLOR="#B22222">/* Newton term     */</FONT></I>       Reynolds*v_y*phi[i][qp]*phi[j][qp]);
  
                  Kvu(i,j) += JxW[qp] * 
   <I><FONT COLOR="#B22222">/* Newton term     */</FONT></I>      -Reynolds*v_x*phi[i][qp]*phi[j][qp];
              }
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (request_jacobian)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j != n_p_dofs; j++)
                {
                  Kup(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](0);
                  Kvp(i,j) += JxW[qp]*psi[j][qp]*dphi[i][qp](1);
                }
          }
      } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>    
    <B><FONT COLOR="#A020F0">return</FONT></B> request_jacobian;
  }
  
  
  
  <FONT COLOR="#228B22"><B>bool</FONT></B> NavierSystem::element_constraint (<FONT COLOR="#228B22"><B>bool</FONT></B> request_jacobian)
  {
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt; &amp;JxW = fe_velocity-&gt;get_JxW();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt; &amp;dphi =
      fe_velocity-&gt;get_dphi();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt; &amp;psi = fe_pressure-&gt;get_phi();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_u_dofs = dof_indices_var[0].size();
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_p_dofs = dof_indices_var[2].size();
  
    DenseSubMatrix&lt;Number&gt; &amp;Kpu = *elem_subjacobians[2][0];
    DenseSubMatrix&lt;Number&gt; &amp;Kpv = *elem_subjacobians[2][1];
    DenseSubVector&lt;Number&gt; &amp;Fp = *elem_subresiduals[2];
  
    <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_qpoints = element_qrule-&gt;n_points();
    <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp != n_qpoints; qp++)
      {
        Gradient grad_u = interior_gradient(0, qp),
                 grad_v = interior_gradient(1, qp);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i != n_p_dofs; i++)
          {
            Fp(i) += JxW[qp] * psi[i][qp] *
                     (grad_u(0) + grad_v(1));
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (request_jacobian)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j != n_u_dofs; j++)
                {
                  Kpu(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](0);
                  Kpv(i,j) += JxW[qp]*psi[i][qp]*dphi[j][qp](1);
                }
          }
      } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> request_jacobian;
  }
  
        
  
  <FONT COLOR="#228B22"><B>bool</FONT></B> NavierSystem::side_constraint (<FONT COLOR="#228B22"><B>bool</FONT></B> request_jacobian)
  {
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt; &amp;JxW_side = fe_side_vel-&gt;get_JxW();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt; &amp;phi_side =
      fe_side_vel-&gt;get_phi();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Point &gt; &amp;qside_point = fe_side_vel-&gt;get_xyz();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_u_dofs = dof_indices_var[0].size(); 
  
    DenseSubMatrix&lt;Number&gt; &amp;Kuu = *elem_subjacobians[0][0];
    DenseSubMatrix&lt;Number&gt; &amp;Kvv = *elem_subjacobians[1][1];
  
    DenseSubVector&lt;Number&gt; &amp;Fu = *elem_subresiduals[0];
    DenseSubVector&lt;Number&gt; &amp;Fv = *elem_subresiduals[1];
  
  
    <FONT COLOR="#228B22"><B>const</FONT></B> Real penalty = 1.e10;
  
    <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_sidepoints = side_qrule-&gt;n_points();
    <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp != n_sidepoints; qp++)
      {
        Number u = side_value(0, qp),
               v = side_value(1, qp);
  
        <FONT COLOR="#228B22"><B>const</FONT></B> Real yf = qside_point[qp](1);
  
        <FONT COLOR="#228B22"><B>const</FONT></B> Real u_value = (yf &gt; .99) ? 1. : 0.;
                
        <FONT COLOR="#228B22"><B>const</FONT></B> Real v_value = 0.;
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i != n_u_dofs; i++)
          {
            Fu(i) += JxW_side[qp] * penalty *
                     (u - u_value) * phi_side[i][qp];
            Fv(i) += JxW_side[qp] * penalty *
                     (v - v_value) * phi_side[i][qp];
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (request_jacobian)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j != n_u_dofs; j++)
                {
                  Kuu(i,j) += JxW_side[qp] * penalty *
                              phi_side[i][qp] * phi_side[j][qp];
                  Kvv(i,j) += JxW_side[qp] * penalty *
                              phi_side[i][qp] * phi_side[j][qp];
                }
          }
      }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> request_jacobian;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example  ./ex18
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
    Type "PDE"
    Variables="u" "v" "p" 
    Finite Element Types="LAGRANGE" "LAGRANGE" "LAGRANGE" 
    Approximation Orders="SECOND" "SECOND" "FIRST" 
    n_dofs()=3803
    n_local_dofs()=3803
    n_constrained_dofs()=0
    n_vectors()=3

 Solving time step 0, time = 0
 Solving time step 1, time = 0.005
 Solving time step 2, time = 0.01
 Solving time step 3, time = 0.015
 Solving time step 4, time = 0.02
 Solving time step 5, time = 0.025
 Solving time step 6, time = 0.03
 Solving time step 7, time = 0.035
 Solving time step 8, time = 0.04
 Solving time step 9, time = 0.045
 Solving time step 10, time = 0.05
 Solving time step 11, time = 0.055
 Solving time step 12, time = 0.06
 Solving time step 13, time = 0.065
 Solving time step 14, time = 0.07
 
***************************************************************
* Done Running Example  ./ex18
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
