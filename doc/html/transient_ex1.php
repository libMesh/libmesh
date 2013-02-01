<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("transient_ex1",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file exact_solution.C with comments: </h1> 
<div class = "comment">
  

<br><br>This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
  

<br><br>This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
  

<br><br>You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


<br><br>

<br><br>

<br><br>C++ Includes
</div>

<div class ="fragment">
<pre>
        #include &lt;math.h&gt;
        
</pre>
</div>
<div class = "comment">
Mesh library includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh_common.h"
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        
        
        
        
        /**
         *
         */
        Real exact_solution (const Real x,
        		     const Real y,
        		     const Real t)
        {
          const Real xo = 0.2;
          const Real yo = 0.2;
          const Real u  = 0.8;
          const Real v  = 0.8;
         
          const Real num =
            pow(x - u*t - xo, 2.) +
            pow(y - v*t - yo, 2.);
        
          const Real den =
            0.01*(4.*t + 1.);
        
          return exp(-num/den)/(4.*t + 1.);
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file transient_ex1.C with comments: </h1> 
<div class = "comment">
<h1>Transient Example 1 - Solving a Transient Linear System in Parallel</h1>

<br><br>This example shows how a simple, linear transient
system can be solved in parallel.  The system is simple
scalar convection-diffusion with a specified external
velocity.  The initial condition is given, and the
solution is advanced in time with a standard Crank-Nicolson
time-stepping strategy.


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
        #include "libmesh/mesh_refinement.h"
        #include "libmesh/gmv_io.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/fe.h"
        #include "libmesh/quadrature_gauss.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/dense_vector.h"
        
</pre>
</div>
<div class = "comment">
This example will solve a linear transient system,
so we need to include the \p TransientLinearImplicitSystem definition.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/linear_implicit_system.h"
        #include "libmesh/transient_system.h"
        #include "libmesh/vector_value.h"
        
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
matrix and right-hand-side at each time step.  Note that
since the system is linear we technically do not need to
assmeble the matrix at each time step, but we will anyway.
In subsequent examples we will employ adaptive mesh refinement,
and with a changing mesh it will be necessary to rebuild the
system matrix.
</div>

<div class ="fragment">
<pre>
        void assemble_cd (EquationSystems& es,
                          const std::string& system_name);
        
</pre>
</div>
<div class = "comment">
Function prototype.  This function will initialize the system.
Initialization functions are optional for systems.  They allow
you to specify the initial values of the solution.  If an
initialization function is not provided then the default (0)
solution is provided.
</div>

<div class ="fragment">
<pre>
        void init_cd (EquationSystems& es,
                      const std::string& system_name);
        
</pre>
</div>
<div class = "comment">
Exact solution function prototype.  This gives the exact
solution as a function of space and time.  In this case the
initial condition will be taken as the exact solution at time 0,
as will the Dirichlet boundary conditions at time t.
</div>

<div class ="fragment">
<pre>
        Real exact_solution (const Real x,
                             const Real y,
                             const Real t);
        
        Number exact_value (const Point& p,
                            const Parameters& parameters,
                            const std::string&,
                            const std::string&)
        {
          return exact_solution(p(0), p(1), parameters.get&lt;Real&gt; ("time"));
        }
        
        
        
</pre>
</div>
<div class = "comment">
We can now begin the main program.  Note that this
example will fail if you are using complex numbers
since it was designed to be run only with real numbers.
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
This example requires Adaptive Mesh Refinement support - although
it only refines uniformly, the refinement code used is the same
underneath
</div>

<div class ="fragment">
<pre>
        #ifndef LIBMESH_ENABLE_AMR
          libmesh_example_assert(false, "--enable-amr");
        #else
        
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
Read the mesh from file.  This is the coarse mesh that will be used
in example 10 to demonstrate adaptive mesh refinement.  Here we will
simply read it in and uniformly refine it 5 times before we compute
with it.
</div>

<div class ="fragment">
<pre>
          Mesh mesh;
              
          mesh.read ("mesh.xda");
          
</pre>
</div>
<div class = "comment">
Create a MeshRefinement object to handle refinement of our mesh.
This class handles all the details of mesh refinement and coarsening.
</div>

<div class ="fragment">
<pre>
          MeshRefinement mesh_refinement (mesh);
          
</pre>
</div>
<div class = "comment">
Uniformly refine the mesh 5 times.  This is the
first time we use the mesh refinement capabilities
of the library.
</div>

<div class ="fragment">
<pre>
          mesh_refinement.uniformly_refine (5);
          
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
Add a transient system to the EquationSystems
object named "Convection-Diffusion".
</div>

<div class ="fragment">
<pre>
          TransientLinearImplicitSystem & system = 
            equation_systems.add_system&lt;TransientLinearImplicitSystem&gt; ("Convection-Diffusion");
            
</pre>
</div>
<div class = "comment">
Adds the variable "u" to "Convection-Diffusion".  "u"
will be approximated using first-order approximation.
</div>

<div class ="fragment">
<pre>
          system.add_variable ("u", FIRST);
            
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly
and initialization functions.
</div>

<div class ="fragment">
<pre>
          system.attach_assemble_function (assemble_cd);
          system.attach_init_function (init_cd);
            
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
Write out the initial conditions.
</div>

<div class ="fragment">
<pre>
          GMVIO(mesh).write_equation_systems ("out_000.gmv",
                                              equation_systems);
          
</pre>
</div>
<div class = "comment">
The Convection-Diffusion system requires that we specify
the flow velocity.  We will specify it as a RealVectorValue
data type and then use the Parameters object to pass it to
the assemble function.
</div>

<div class ="fragment">
<pre>
          equation_systems.parameters.set&lt;RealVectorValue&gt;("velocity") = 
            RealVectorValue (0.8, 0.8);
          
</pre>
</div>
<div class = "comment">
Solve the system "Convection-Diffusion".  This will be done by
looping over the specified time interval and calling the
solve() member at each time step.  This will assemble the
system and call the linear solver.
</div>

<div class ="fragment">
<pre>
          const Real dt = 0.025;
          system.time   = 0.;
          
          for (unsigned int t_step = 0; t_step &lt; 50; t_step++)
            {
</pre>
</div>
<div class = "comment">
Incremenet the time counter, set the time and the
time step size as parameters in the EquationSystem.
</div>

<div class ="fragment">
<pre>
              system.time += dt;
        
              equation_systems.parameters.set&lt;Real&gt; ("time") = system.time;
              equation_systems.parameters.set&lt;Real&gt; ("dt")   = dt;
        
</pre>
</div>
<div class = "comment">
A pretty update message
</div>

<div class ="fragment">
<pre>
              std::cout &lt;&lt; " Solving time step ";
              
</pre>
</div>
<div class = "comment">
Since some compilers fail to offer full stream
functionality, libMesh offers a string stream
to work around this.  Note that for other compilers,
this is just a set of preprocessor macros and therefore
should cost nothing (compared to a hand-coded string stream).
We use additional curly braces here simply to enforce data
locality.
</div>

<div class ="fragment">
<pre>
              {
                std::ostringstream out;
        
                out &lt;&lt; std::setw(2)
                    &lt;&lt; std::right
                    &lt;&lt; t_step
                    &lt;&lt; ", time="
                    &lt;&lt; std::fixed
                    &lt;&lt; std::setw(6)
                    &lt;&lt; std::setprecision(3)
                    &lt;&lt; std::setfill('0')
                    &lt;&lt; std::left
                    &lt;&lt; system.time
                    &lt;&lt;  "...";
        
                std::cout &lt;&lt; out.str() &lt;&lt; std::endl;
              }
              
</pre>
</div>
<div class = "comment">
At this point we need to update the old
solution vector.  The old solution vector
will be the current solution vector from the
previous time step.  We will do this by extracting the
system from the \p EquationSystems object and using
vector assignment.  Since only \p TransientSystems
(and systems derived from them) contain old solutions
we need to specify the system type when we ask for it.
</div>

<div class ="fragment">
<pre>
              TransientLinearImplicitSystem&  system =
                equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;("Convection-Diffusion");
        
              *system.old_local_solution = *system.current_local_solution;
              
</pre>
</div>
<div class = "comment">
Assemble & solve the linear system
</div>

<div class ="fragment">
<pre>
              equation_systems.get_system("Convection-Diffusion").solve();
              
</pre>
</div>
<div class = "comment">
Output evey 10 timesteps to file.
</div>

<div class ="fragment">
<pre>
              if ( (t_step+1)%10 == 0)
                {
                  std::ostringstream file_name;
        
                  file_name &lt;&lt; "out_"
                            &lt;&lt; std::setw(3)
                            &lt;&lt; std::setfill('0')
                            &lt;&lt; std::right
                            &lt;&lt; t_step+1
                            &lt;&lt; ".gmv";
        
                  GMVIO(mesh).write_equation_systems (file_name.str(),
                                                      equation_systems);
                }
            }
        #endif // #ifdef LIBMESH_ENABLE_AMR
        
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
We now define the function which provides the
initialization routines for the "Convection-Diffusion"
system.  This handles things like setting initial
conditions and boundary conditions.
</div>

<div class ="fragment">
<pre>
        void init_cd (EquationSystems& es,
                      const std::string& system_name)
        {
</pre>
</div>
<div class = "comment">
It is a good idea to make sure we are initializing
the proper system.
</div>

<div class ="fragment">
<pre>
          libmesh_assert_equal_to (system_name, "Convection-Diffusion");
        
</pre>
</div>
<div class = "comment">
Get a reference to the Convection-Diffusion system object.
</div>

<div class ="fragment">
<pre>
          TransientLinearImplicitSystem & system =
            es.get_system&lt;TransientLinearImplicitSystem&gt;("Convection-Diffusion");
        
</pre>
</div>
<div class = "comment">
Project initial conditions at time 0
</div>

<div class ="fragment">
<pre>
          es.parameters.set&lt;Real&gt; ("time") = system.time = 0;
          
          system.project_solution(exact_value, NULL, es.parameters);
        }
        
        
        
</pre>
</div>
<div class = "comment">
Now we define the assemble function which will be used
by the EquationSystems object at each timestep to assemble
the linear system for solution.
</div>

<div class ="fragment">
<pre>
        void assemble_cd (EquationSystems& es,
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
          libmesh_assert_equal_to (system_name, "Convection-Diffusion");
          
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
Get a reference to the Convection-Diffusion system object.
</div>

<div class ="fragment">
<pre>
          TransientLinearImplicitSystem & system =
            es.get_system&lt;TransientLinearImplicitSystem&gt; ("Convection-Diffusion");
          
</pre>
</div>
<div class = "comment">
Get a constant reference to the Finite Element type
for the first (and only) variable in the system.
</div>

<div class ="fragment">
<pre>
          FEType fe_type = system.variable_type(0);
          
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
A Gauss quadrature rule for numerical integration.
Let the \p FEType object decide what order rule is appropriate.
</div>

<div class ="fragment">
<pre>
          QGauss qrule (dim,   fe_type.default_quadrature_order());
          QGauss qface (dim-1, fe_type.default_quadrature_order());
        
</pre>
</div>
<div class = "comment">
Tell the finite element object to use our quadrature rule.
</div>

<div class ="fragment">
<pre>
          fe-&gt;attach_quadrature_rule      (&qrule);
          fe_face-&gt;attach_quadrature_rule (&qface);
        
</pre>
</div>
<div class = "comment">
Here we define some references to cell-specific data that
will be used to assemble the linear system.  We will start
with the element Jacobian * quadrature weight at each integration point.   
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt;& JxW      = fe-&gt;get_JxW();
          const std::vector&lt;Real&gt;& JxW_face = fe_face-&gt;get_JxW();
          
</pre>
</div>
<div class = "comment">
The element shape functions evaluated at the quadrature points.
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
Define data structures to contain the element matrix
and right-hand-side vector contribution.  Following
basic finite element terminology we will denote these
"Ke" and "Fe".
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
Here we extract the velocity & parameters that we put in the
EquationSystems object.
</div>

<div class ="fragment">
<pre>
          const RealVectorValue velocity =
            es.parameters.get&lt;RealVectorValue&gt; ("velocity");
        
          const Real dt = es.parameters.get&lt;Real&gt;   ("dt");
        
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
Now we will build the element matrix and right-hand-side.
Constructing the RHS requires the solution and its
gradient from the previous timestep.  This myst be
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
Values to hold the old solution & its gradient.
</div>

<div class ="fragment">
<pre>
                  Number   u_old = 0.;
                  Gradient grad_u_old;
                  
</pre>
</div>
<div class = "comment">
Compute the old solution & its gradient.
</div>

<div class ="fragment">
<pre>
                  for (unsigned int l=0; l&lt;phi.size(); l++)
                    {
                      u_old      += phi[l][qp]*system.old_solution  (dof_indices[l]);
                      
</pre>
</div>
<div class = "comment">
This will work,
grad_u_old += dphi[l][qp]*system.old_solution (dof_indices[l]);
but we can do it without creating a temporary like this:
</div>

<div class ="fragment">
<pre>
                      grad_u_old.add_scaled (dphi[l][qp],system.old_solution (dof_indices[l]));
                    }
                  
</pre>
</div>
<div class = "comment">
Now compute the element matrix and RHS contributions.
</div>

<div class ="fragment">
<pre>
                  for (unsigned int i=0; i&lt;phi.size(); i++)
                    {
</pre>
</div>
<div class = "comment">
The RHS contribution
</div>

<div class ="fragment">
<pre>
                      Fe(i) += JxW[qp]*(
</pre>
</div>
<div class = "comment">
Mass matrix term
</div>

<div class ="fragment">
<pre>
                                        u_old*phi[i][qp] +
                                        -.5*dt*(
</pre>
</div>
<div class = "comment">
Convection term
(grad_u_old may be complex, so the
order here is important!)
</div>

<div class ="fragment">
<pre>
                                                (grad_u_old*velocity)*phi[i][qp] +
                                                
</pre>
</div>
<div class = "comment">
Diffusion term
</div>

<div class ="fragment">
<pre>
                                                0.01*(grad_u_old*dphi[i][qp]))     
                                        );
                      
                      for (unsigned int j=0; j&lt;phi.size(); j++)
                        {
</pre>
</div>
<div class = "comment">
The matrix contribution
</div>

<div class ="fragment">
<pre>
                          Ke(i,j) += JxW[qp]*(
</pre>
</div>
<div class = "comment">
Mass-matrix
</div>

<div class ="fragment">
<pre>
                                              phi[i][qp]*phi[j][qp] + 
        
                                              .5*dt*(
</pre>
</div>
<div class = "comment">
Convection term
</div>

<div class ="fragment">
<pre>
                                                     (velocity*dphi[j][qp])*phi[i][qp] +
                                                     
</pre>
</div>
<div class = "comment">
Diffusion term
</div>

<div class ="fragment">
<pre>
                                                     0.01*(dphi[i][qp]*dphi[j][qp]))      
                                              );
                        } 
                    } 
                } 
        
</pre>
</div>
<div class = "comment">
At this point the interior element integration has
been completed.  However, we have not yet addressed
boundary conditions.  For this example we will only
consider simple Dirichlet boundary conditions imposed
via the penalty method. 

<br><br>The following loops over the sides of the element.
If the element has no neighbor on a side then that
side MUST live on a boundary of the domain.
</div>

<div class ="fragment">
<pre>
              {
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
                      
                      for (unsigned int qp=0; qp&lt;qface.n_points(); qp++)
                        {
                          const Number value = exact_solution (qface_points[qp](0),
                                                               qface_points[qp](1),
                                                               system.time);
                                                               
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
              } 
        
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
              system.matrix-&gt;add_matrix (Ke, dof_indices);
              system.rhs-&gt;add_vector    (Fe, dof_indices);
            }
          
</pre>
</div>
<div class = "comment">
That concludes the system matrix assembly routine.
</div>

<div class ="fragment">
<pre>
        #endif // #ifdef LIBMESH_ENABLE_AMR
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file exact_solution.C without comments: </h1> 
<pre> 
    
    
    
  
  
  
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_common.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  
  
  
  <I><FONT COLOR="#B22222">/**
   *
   */</FONT></I>
  Real exact_solution (<B><FONT COLOR="#228B22">const</FONT></B> Real x,
  		     <B><FONT COLOR="#228B22">const</FONT></B> Real y,
  		     <B><FONT COLOR="#228B22">const</FONT></B> Real t)
  {
    <B><FONT COLOR="#228B22">const</FONT></B> Real xo = 0.2;
    <B><FONT COLOR="#228B22">const</FONT></B> Real yo = 0.2;
    <B><FONT COLOR="#228B22">const</FONT></B> Real u  = 0.8;
    <B><FONT COLOR="#228B22">const</FONT></B> Real v  = 0.8;
   
    <B><FONT COLOR="#228B22">const</FONT></B> Real num =
      pow(x - u*t - xo, 2.) +
      pow(y - v*t - yo, 2.);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real den =
      0.01*(4.*t + 1.);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> exp(-num/den)/(4.*t + 1.);
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file transient_ex1.C without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;sstream&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/gmv_io.h&quot;</FONT></B>
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
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/vector_value.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_cd (EquationSystems&amp; es,
                    <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  <B><FONT COLOR="#228B22">void</FONT></B> init_cd (EquationSystems&amp; es,
                <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  Real exact_solution (<B><FONT COLOR="#228B22">const</FONT></B> Real x,
                       <B><FONT COLOR="#228B22">const</FONT></B> Real y,
                       <B><FONT COLOR="#228B22">const</FONT></B> Real t);
  
  Number exact_value (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                      <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp; parameters,
                      <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,
                      <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;)
  {
    <B><FONT COLOR="#A020F0">return</FONT></B> exact_solution(p(0), p(1), parameters.get&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>));
  }
  
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #ifndef LIBMESH_ENABLE_AMR
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-amr&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
  
    Mesh mesh;
        
    mesh.read (<B><FONT COLOR="#BC8F8F">&quot;mesh.xda&quot;</FONT></B>);
    
    MeshRefinement mesh_refinement (mesh);
    
    mesh_refinement.uniformly_refine (5);
    
    mesh.print_info();
    
    EquationSystems equation_systems (mesh);
    
    TransientLinearImplicitSystem &amp; system = 
      equation_systems.add_system&lt;TransientLinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
      
    system.add_variable (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, FIRST);
      
    system.attach_assemble_function (assemble_cd);
    system.attach_init_function (init_cd);
      
    equation_systems.init ();
      
    equation_systems.print_info();
      
    GMVIO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;out_000.gmv&quot;</FONT></B>,
                                        equation_systems);
    
    equation_systems.parameters.set&lt;RealVectorValue&gt;(<B><FONT COLOR="#BC8F8F">&quot;velocity&quot;</FONT></B>) = 
      RealVectorValue (0.8, 0.8);
    
    <B><FONT COLOR="#228B22">const</FONT></B> Real dt = 0.025;
    system.time   = 0.;
    
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> t_step = 0; t_step &lt; 50; t_step++)
      {
        system.time += dt;
  
        equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>) = system.time;
        equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;dt&quot;</FONT></B>)   = dt;
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Solving time step &quot;</FONT></B>;
        
        {
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream out;
  
          out &lt;&lt; std::setw(2)
              &lt;&lt; std::right
              &lt;&lt; t_step
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, time=&quot;</FONT></B>
              &lt;&lt; std::fixed
              &lt;&lt; std::setw(6)
              &lt;&lt; std::setprecision(3)
              &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
              &lt;&lt; std::left
              &lt;&lt; system.time
              &lt;&lt;  <B><FONT COLOR="#BC8F8F">&quot;...&quot;</FONT></B>;
  
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; out.str() &lt;&lt; std::endl;
        }
        
        TransientLinearImplicitSystem&amp;  system =
          equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
        *system.old_local_solution = *system.current_local_solution;
        
        equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>).solve();
        
        <B><FONT COLOR="#A020F0">if</FONT></B> ( (t_step+1)%10 == 0)
          {
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream file_name;
  
            file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out_&quot;</FONT></B>
                      &lt;&lt; std::setw(3)
                      &lt;&lt; std::setfill(<B><FONT COLOR="#BC8F8F">'0'</FONT></B>)
                      &lt;&lt; std::right
                      &lt;&lt; t_step+1
                      &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.gmv&quot;</FONT></B>;
  
            GMVIO(mesh).write_equation_systems (file_name.str(),
                                                equation_systems);
          }
      }
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_ENABLE_AMR
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> init_cd (EquationSystems&amp; es,
                <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    TransientLinearImplicitSystem &amp; system =
      es.get_system&lt;TransientLinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    es.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>) = system.time = 0;
    
    system.project_solution(exact_value, NULL, es.parameters);
  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_cd (EquationSystems&amp; es,
                    <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
  #ifdef LIBMESH_ENABLE_AMR
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
    
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
    
    TransientLinearImplicitSystem &amp; system =
      es.get_system&lt;TransientLinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
    
    FEType fe_type = system.variable_type(0);
    
    AutoPtr&lt;FEBase&gt; fe      (FEBase::build(dim, fe_type));
    AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
    
    QGauss qrule (dim,   fe_type.default_quadrature_order());
    QGauss qface (dim-1, fe_type.default_quadrature_order());
  
    fe-&gt;attach_quadrature_rule      (&amp;qrule);
    fe_face-&gt;attach_quadrature_rule (&amp;qface);
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW      = fe-&gt;get_JxW();
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
    
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; psi = fe_face-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; qface_points = fe_face-&gt;get_xyz();
    
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
  
    DenseMatrix&lt;Number&gt; Ke;
    DenseVector&lt;Number&gt; Fe;
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
  
    <B><FONT COLOR="#228B22">const</FONT></B> RealVectorValue velocity =
      es.parameters.get&lt;RealVectorValue&gt; (<B><FONT COLOR="#BC8F8F">&quot;velocity&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real dt = es.parameters.get&lt;Real&gt;   (<B><FONT COLOR="#BC8F8F">&quot;dt&quot;</FONT></B>);
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {    
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        fe-&gt;reinit (elem);
        
        Ke.resize (dof_indices.size(),
                   dof_indices.size());
  
        Fe.resize (dof_indices.size());
        
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
          {
            Number   u_old = 0.;
            Gradient grad_u_old;
            
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> l=0; l&lt;phi.size(); l++)
              {
                u_old      += phi[l][qp]*system.old_solution  (dof_indices[l]);
                
                grad_u_old.add_scaled (dphi[l][qp],system.old_solution (dof_indices[l]));
              }
            
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
              {
                Fe(i) += JxW[qp]*(
                                  u_old*phi[i][qp] +
                                  -.5*dt*(
                                          (grad_u_old*velocity)*phi[i][qp] +
                                          
                                          0.01*(grad_u_old*dphi[i][qp]))     
                                  );
                
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
                  {
                    Ke(i,j) += JxW[qp]*(
                                        phi[i][qp]*phi[j][qp] + 
  
                                        .5*dt*(
                                               (velocity*dphi[j][qp])*phi[i][qp] +
                                               
                                               0.01*(dphi[i][qp]*dphi[j][qp]))      
                                        );
                  } 
              } 
          } 
  
        {
          <B><FONT COLOR="#228B22">const</FONT></B> Real penalty = 1.e10;
  
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> s=0; s&lt;elem-&gt;n_sides(); s++)
            <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(s) == NULL)
              {
                fe_face-&gt;reinit(elem,s);
                
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qface.n_points(); qp++)
                  {
                    <B><FONT COLOR="#228B22">const</FONT></B> Number value = exact_solution (qface_points[qp](0),
                                                         qface_points[qp](1),
                                                         system.time);
                                                         
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;psi.size(); i++)
                      Fe(i) += penalty*JxW_face[qp]*value*psi[i][qp];
  
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;psi.size(); i++)
                      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;psi.size(); j++)
                        Ke(i,j) += penalty*JxW_face[qp]*psi[i][qp]*psi[j][qp];
                  }
              } 
        } 
  
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
  
        system.matrix-&gt;add_matrix (Ke, dof_indices);
        system.rhs-&gt;add_vector    (Fe, dof_indices);
      }
    
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_ENABLE_AMR
</FONT></I>  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example transient_ex1:
*  mpirun -np 12 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=6273
    n_local_nodes()=491
  n_elem()=13650
    n_local_elem()=1229
    n_active_elem()=10240
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=6273
    n_local_dofs()=491
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=3
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.45943
      Average Off-Processor Bandwidth <= 0.323609
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 7
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

 Solving time step  0, time=0.0250...
 Solving time step  1, time=0.0500...
 Solving time step  2, time=0.0750...
 Solving time step  3, time=0.1000...
 Solving time step  4, time=0.1250...
 Solving time step  5, time=0.1500...
 Solving time step  6, time=0.1750...
 Solving time step  7, time=0.2000...
 Solving time step  8, time=0.2250...
 Solving time step  9, time=0.2500...
 Solving time step 10, time=0.2750...
 Solving time step 11, time=0.3000...
 Solving time step 12, time=0.3250...
 Solving time step 13, time=0.3500...
 Solving time step 14, time=0.3750...
 Solving time step 15, time=0.4000...
 Solving time step 16, time=0.4250...
 Solving time step 17, time=0.4500...
 Solving time step 18, time=0.4750...
 Solving time step 19, time=0.5000...
 Solving time step 20, time=0.5250...
 Solving time step 21, time=0.5500...
 Solving time step 22, time=0.5750...
 Solving time step 23, time=0.6000...
 Solving time step 24, time=0.6250...
 Solving time step 25, time=0.6500...
 Solving time step 26, time=0.6750...
 Solving time step 27, time=0.7000...
 Solving time step 28, time=0.7250...
 Solving time step 29, time=0.7500...
 Solving time step 30, time=0.7750...
 Solving time step 31, time=0.8000...
 Solving time step 32, time=0.8250...
 Solving time step 33, time=0.8500...
 Solving time step 34, time=0.8750...
 Solving time step 35, time=0.9000...
 Solving time step 36, time=0.9250...
 Solving time step 37, time=0.9500...
 Solving time step 38, time=0.9750...
 Solving time step 39, time=1.0000...
 Solving time step 40, time=1.0250...
 Solving time step 41, time=1.0500...
 Solving time step 42, time=1.0750...
 Solving time step 43, time=1.1000...
 Solving time step 44, time=1.1250...
 Solving time step 45, time=1.1500...
 Solving time step 46, time=1.1750...
 Solving time step 47, time=1.2000...
 Solving time step 48, time=1.2250...
 Solving time step 49, time=1.2500...
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/transient/transient_ex1/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:21:10 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           7.606e+00      1.00003   7.605e+00
Objects:              3.220e+02      1.00000   3.220e+02
Flops:                6.051e+07      1.84220   3.988e+07  4.786e+08
Flops/sec:            7.957e+06      1.84225   5.244e+06  6.293e+07
MPI Messages:         5.568e+03      2.00000   3.867e+03  4.640e+04
MPI Message Lengths:  1.267e+06      2.15735   2.177e+02  1.010e+07
MPI Reductions:       2.495e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 7.6055e+00 100.0%  4.7860e+08 100.0%  4.640e+04 100.0%  2.177e+02      100.0%  2.494e+03 100.0% 

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

VecMDot              640 1.0 1.8373e-02 2.8 6.46e+06 1.7 0.0e+00 0.0e+00 6.4e+02  0 12  0  0 26   0 12  0  0 26  3061
VecNorm              740 1.0 4.8631e-0211.0 1.07e+06 1.7 0.0e+00 0.0e+00 7.4e+02  0  2  0  0 30   0  2  0  0 30   191
VecScale             690 1.0 5.8055e-04 1.2 4.97e+05 1.7 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  7456
VecCopy              151 1.0 2.4414e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               848 1.0 7.6938e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              100 1.0 1.1456e-02 1.1 1.44e+05 1.7 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   110
VecMAXPY             690 1.0 4.6268e-03 1.6 7.38e+06 1.7 0.0e+00 0.0e+00 0.0e+00  0 13  0  0  0   0 13  0  0  0 13902
VecAssemblyBegin     202 1.0 2.9962e-0137.7 0.00e+00 0.0 2.5e+03 3.3e+02 6.1e+02  4  0  5  8 24   4  0  5  8 24     0
VecAssemblyEnd       202 1.0 4.1866e-04 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      791 1.0 3.8319e-03 1.5 0.00e+00 0.0 4.0e+04 1.6e+02 0.0e+00  0  0 85 61  0   0  0 85 61  0     0
VecScatterEnd        791 1.0 4.7305e-0227.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize         690 1.0 4.8599e-02 9.3 1.49e+06 1.7 0.0e+00 0.0e+00 6.9e+02  0  3  0  0 28   0  3  0  0 28   267
MatMult              690 1.0 5.6543e-02 4.1 8.20e+06 2.1 3.4e+04 1.5e+02 0.0e+00  0 13 74 50  0   0 13 74 50  0  1082
MatSolve             740 1.0 2.4829e-02 1.6 2.56e+07 1.9 0.0e+00 0.0e+00 0.0e+00  0 42  0  0  0   0 42  0  0  0  8038
MatLUFactorNum        50 1.0 2.9424e-02 1.9 1.12e+07 2.2 0.0e+00 0.0e+00 0.0e+00  0 17  0  0  0   0 17  0  0  0  2800
MatILUFactorSym       50 1.0 5.9617e-02 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 1.5e+02  0  0  0  0  6   0  0  0  0  6     0
MatAssemblyBegin     100 1.0 3.6495e-0126.2 0.00e+00 0.0 3.8e+03 8.2e+02 2.0e+02  3  0  8 30  8   3  0  8 30  8     0
MatAssemblyEnd       100 1.0 5.5621e-03 2.4 0.00e+00 0.0 1.0e+02 3.9e+01 8.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ           50 1.0 2.1935e-05 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        50 1.0 1.7724e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+02  0  0  0  0  4   0  0  0  0  4     0
MatZeroEntries        52 1.0 2.7490e-04 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       640 1.0 2.1851e-02 2.0 1.29e+07 1.7 0.0e+00 0.0e+00 6.4e+02  0 24  0  0 26   0 24  0  0 26  5150
KSPSetUp             100 1.0 6.2919e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              50 1.0 1.8000e-01 1.0 6.05e+07 1.8 3.4e+04 1.5e+02 1.6e+03  2100 74 50 65   2100 74 50 65  2659
PCSetUp              100 1.0 1.0066e-01 1.8 1.12e+07 2.2 0.0e+00 0.0e+00 2.5e+02  1 17  0  0 10   1 17  0  0 10   819
PCSetUpOnBlocks       50 1.0 9.9889e-02 1.8 1.12e+07 2.2 0.0e+00 0.0e+00 2.5e+02  1 17  0  0 10   1 17  0  0 10   825
PCApply              740 1.0 3.2255e-02 1.4 2.56e+07 1.9 0.0e+00 0.0e+00 0.0e+00  0 42  0  0  0   0 42  0  0  0  6187
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    91             91       461320     0
      Vector Scatter     6              6         6216     0
           Index Set   162            162       220480     0
   IS L to G Mapping     5              5         2820     0
              Matrix    53             53      6271260     0
       Krylov Solver     2              2        19360     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 1.19209e-07
Average time for MPI_Barrier(): 4.81606e-06
Average time for zero size MPI_Send(): 1.43449e-05
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


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 22:21:10 2013                                                                             |
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
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=7.7609, Active time=7.43228                                                    |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0420      0.042038    0.0506      0.050569    0.57     0.68     |
|   build_sparsity()                 1         0.0210      0.021007    0.0648      0.064809    0.28     0.87     |
|   create_dof_constraints()         1         0.0461      0.046069    0.0461      0.046069    0.62     0.62     |
|   distribute_dofs()                1         0.1445      0.144493    0.4767      0.476730    1.94     6.41     |
|   dof_indices()                    52053     2.2193      0.000043    2.2193      0.000043    29.86    29.86    |
|   prepare_send_list()              1         0.0001      0.000113    0.0001      0.000113    0.00     0.00     |
|   reinit()                         1         0.3051      0.305103    0.3051      0.305103    4.11     4.11     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          6         0.0793      0.013211    0.3317      0.055282    1.07     4.46     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        45550     0.4217      0.000009    0.4217      0.000009    5.67     5.67     |
|   init_shape_functions()           900       0.0057      0.000006    0.0057      0.000006    0.08     0.08     |
|   inverse_map()                    1700      0.0107      0.000006    0.0107      0.000006    0.14     0.14     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             45550     0.4011      0.000009    0.4011      0.000009    5.40     5.40     |
|   compute_face_map()               850       0.0141      0.000017    0.0252      0.000030    0.19     0.34     |
|   init_face_shape_functions()      50        0.0008      0.000017    0.0008      0.000017    0.01     0.01     |
|   init_reference_to_physical_map() 900       0.0084      0.000009    0.0084      0.000009    0.11     0.11     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               6         0.3875      0.064580    0.3882      0.064707    5.21     5.22     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           32736     0.1201      0.000004    0.1201      0.000004    1.62     1.62     |
|   init()                           5         0.0036      0.000726    0.0036      0.000726    0.05     0.05     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 2         0.2709      0.135437    0.2880      0.143984    3.64     3.87     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   broadcast()                      1         0.0011      0.001127    0.0017      0.001722    0.02     0.02     |
|   compute_hilbert_indices()        3         0.0814      0.027149    0.0814      0.027149    1.10     1.10     |
|   find_global_indices()            3         0.0367      0.012233    0.1297      0.043227    0.49     1.74     |
|   parallel_sort()                  3         0.0037      0.001227    0.0064      0.002144    0.05     0.09     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         6         0.0005      0.000083    0.7216      0.120274    0.01     9.71     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _refine_elements()               5         0.2304      0.046089    0.4978      0.099555    3.10     6.70     |
|   add_point()                      32736     0.1216      0.000004    0.2487      0.000008    1.64     3.35     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      2         0.4168      0.208407    0.5442      0.272112    5.61     7.32     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      11        0.0169      0.001535    0.0173      0.001569    0.23     0.23     |
|   broadcast()                      9         0.0004      0.000044    0.0003      0.000037    0.01     0.00     |
|   max(bool)                        6         0.0123      0.002043    0.0123      0.002043    0.16     0.16     |
|   max(scalar)                      272       0.0035      0.000013    0.0035      0.000013    0.05     0.05     |
|   max(vector)                      61        0.0010      0.000016    0.0027      0.000044    0.01     0.04     |
|   min(bool)                        322       0.0029      0.000009    0.0029      0.000009    0.04     0.04     |
|   min(scalar)                      264       0.0743      0.000281    0.0743      0.000281    1.00     1.00     |
|   min(vector)                      61        0.0011      0.000018    0.0033      0.000055    0.01     0.04     |
|   probe()                          176       0.0120      0.000068    0.0120      0.000068    0.16     0.16     |
|   receive()                        176       0.0014      0.000008    0.0134      0.000076    0.02     0.18     |
|   send()                           176       0.0007      0.000004    0.0007      0.000004    0.01     0.01     |
|   send_receive()                   182       0.0016      0.000009    0.0162      0.000089    0.02     0.22     |
|   sum()                            39        0.0054      0.000138    0.0293      0.000751    0.07     0.39     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           176       0.0004      0.000002    0.0004      0.000002    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         2         0.0263      0.013135    0.0357      0.017872    0.35     0.48     |
|   set_parent_processor_ids()       2         0.0311      0.015535    0.0311      0.015535    0.42     0.42     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          50        0.5567      0.011134    0.5567      0.011134    7.49     7.49     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       50        1.2697      0.025395    4.0762      0.081524    17.08    54.84    |
|   project_vector()                 1         0.0204      0.020400    0.0598      0.059826    0.27     0.80     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            215109    7.4323                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example transient_ex1:
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
