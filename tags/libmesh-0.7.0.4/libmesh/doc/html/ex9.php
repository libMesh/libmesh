<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex9",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 9 - Solving a Transient Linear System in Parallel</h1>

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
        #include "mesh_refinement.h"
        #include "gmv_io.h"
        #include "equation_systems.h"
        #include "fe.h"
        #include "quadrature_gauss.h"
        #include "dof_map.h"
        #include "sparse_matrix.h"
        #include "numeric_vector.h"
        #include "dense_matrix.h"
        #include "dense_vector.h"
        
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
This example will solve a linear transient system,
so we need to include the \p TransientLinearImplicitSystem definition.
</div>

<div class ="fragment">
<pre>
        #include "linear_implicit_system.h"
        #include "transient_system.h"
        #include "vector_value.h"
        
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
              
          mesh.read ("../ex10/mesh.xda");
          
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
          Real time     = 0.;
          
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
              time += dt;
        
              equation_systems.parameters.set&lt;Real&gt; ("time") = time;
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
                OStringStream out;
        
                OSSInt(out,2,t_step);
                out &lt;&lt; ", time=";
                OSSRealzeroleft(out,6,3,time);
                out &lt;&lt;  "...";
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
                  OStringStream file_name;
        
                  file_name &lt;&lt; "out_";
                  OSSRealzeroright(file_name,3,0,t_step+1);
                  file_name &lt;&lt; ".gmv";
        
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
          libmesh_assert (system_name == "Convection-Diffusion");
        
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
          es.parameters.set&lt;Real&gt; ("time") = 0;
          
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
          libmesh_assert (system_name == "Convection-Diffusion");
          
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
          std::vector&lt;unsigned int&gt; dof_indices;
        
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
          const Real time = es.parameters.get&lt;Real&gt; ("time");
        
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
                                                               time);
                                                               
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
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;gmv_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;o_string_stream.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;transient_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;vector_value.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  
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
        
    mesh.read (<B><FONT COLOR="#BC8F8F">&quot;../ex10/mesh.xda&quot;</FONT></B>);
    
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
    Real time     = 0.;
    
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> t_step = 0; t_step &lt; 50; t_step++)
      {
        time += dt;
  
        equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>) = time;
        equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;dt&quot;</FONT></B>)   = dt;
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; Solving time step &quot;</FONT></B>;
        
        {
          OStringStream out;
  
          OSSInt(out,2,t_step);
          out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, time=&quot;</FONT></B>;
          OSSRealzeroleft(out,6,3,time);
          out &lt;&lt;  <B><FONT COLOR="#BC8F8F">&quot;...&quot;</FONT></B>;
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; out.str() &lt;&lt; std::endl;
        }
        
        TransientLinearImplicitSystem&amp;  system =
          equation_systems.get_system&lt;TransientLinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
        *system.old_local_solution = *system.current_local_solution;
        
        equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>).solve();
        
        <B><FONT COLOR="#A020F0">if</FONT></B> ( (t_step+1)%10 == 0)
          {
            OStringStream file_name;
  
            file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out_&quot;</FONT></B>;
            OSSRealzeroright(file_name,3,0,t_step+1);
            file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.gmv&quot;</FONT></B>;
  
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
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    TransientLinearImplicitSystem &amp; system =
      es.get_system&lt;TransientLinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
  
    es.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>) = 0;
    
    system.project_solution(exact_value, NULL, es.parameters);
  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_cd (EquationSystems&amp; es,
                    <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
  #ifdef LIBMESH_ENABLE_AMR
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Convection-Diffusion&quot;</FONT></B>);
    
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
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
  
    <B><FONT COLOR="#228B22">const</FONT></B> RealVectorValue velocity =
      es.parameters.get&lt;RealVectorValue&gt; (<B><FONT COLOR="#BC8F8F">&quot;velocity&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real dt = es.parameters.get&lt;Real&gt;   (<B><FONT COLOR="#BC8F8F">&quot;dt&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> Real time = es.parameters.get&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>);
  
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
                                                         time);
                                                         
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
* Running Example  mpirun -np 2 ./ex9-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=6273
    n_local_nodes()=3171
  n_elem()=13650
    n_local_elem()=6839
    n_active_elem()=10240
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Convection-Diffusion"
    Type "TransientLinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=6273
    n_local_dofs()=3171
    n_constrained_dofs()=0
    n_vectors()=3

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
 Solving time step 40, time=1.0300...
 Solving time step 41, time=1.0500...
 Solving time step 42, time=1.0800...
 Solving time step 43, time=1.1000...
 Solving time step 44, time=1.1200...
 Solving time step 45, time=1.1500...
 Solving time step 46, time=1.1700...
 Solving time step 47, time=1.2000...
 Solving time step 48, time=1.2200...
 Solving time step 49, time=1.2500...
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex9-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Thu Feb  3 12:14:14 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.824e+00      1.01035   1.815e+00
Objects:              3.530e+02      1.00000   3.530e+02
Flops:                1.390e+08      1.00538   1.386e+08  2.772e+08
Flops/sec:            7.656e+07      1.00495   7.637e+07  1.527e+08
MPI Messages:         7.460e+02      1.00000   7.460e+02  1.492e+03
MPI Message Lengths:  1.884e+06      1.01517   2.507e+03  3.740e+06
MPI Reductions:       1.925e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.8151e+00 100.0%  2.7722e+08 100.0%  1.492e+03 100.0%  2.507e+03      100.0%  1.805e+03  93.8% 

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

VecMDot              308 1.0 3.5172e-03 1.1 8.16e+06 1.0 0.0e+00 0.0e+00 3.1e+02  0  6  0  0 16   0  6  0  0 17  4590
VecNorm              408 1.0 9.7795e-03 1.0 2.59e+06 1.0 0.0e+00 0.0e+00 4.1e+02  1  2  0  0 21   1  2  0  0 23   523
VecScale             358 1.0 4.3201e-04 1.0 1.14e+06 1.0 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  5198
VecCopy              200 1.0 4.4799e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               465 1.0 5.5814e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY              100 1.0 3.9990e-03 1.0 6.34e+05 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   314
VecMAXPY             358 1.0 2.6245e-03 1.0 1.01e+07 1.0 0.0e+00 0.0e+00 0.0e+00  0  7  0  0  0   0  7  0  0  0  7625
VecAssemblyBegin     152 1.0 1.4167e-03 1.0 0.00e+00 0.0 1.0e+02 1.0e+03 4.6e+02  0  0  7  3 24   0  0  7  3 25     0
VecAssemblyEnd       152 1.0 1.2326e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      460 1.0 1.1163e-03 1.0 0.00e+00 0.0 8.2e+02 5.9e+02 0.0e+00  0  0 55 13  0   0  0 55 13  0     0
VecScatterEnd        460 1.0 3.2108e-03 4.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize         358 1.0 9.3317e-03 1.0 3.41e+06 1.0 0.0e+00 0.0e+00 3.6e+02  1  2  0  0 19   1  2  0  0 20   722
MatMult              358 1.0 1.5127e-02 1.2 1.61e+07 1.0 7.2e+02 5.5e+02 0.0e+00  1 11 48 11  0   1 11 48 11  0  2099
MatSolve             358 1.0 3.9714e-02 1.0 5.37e+07 1.0 0.0e+00 0.0e+00 0.0e+00  2 39  0  0  0   2 39  0  0  0  2689
MatLUFactorNum        50 1.0 7.9542e-02 1.0 4.73e+07 1.0 0.0e+00 0.0e+00 0.0e+00  4 34  0  0  0   4 34  0  0  0  1181
MatILUFactorSym       50 1.0 1.9598e-01 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 5.0e+01 11  0  0  0  3  11  0  0  0  3     0
MatAssemblyBegin     100 1.0 9.6035e-0225.3 0.00e+00 0.0 1.5e+02 3.1e+03 2.0e+02  3  0 10 12 10   3  0 10 12 11     0
MatAssemblyEnd       100 1.0 7.1003e-03 1.1 0.00e+00 0.0 4.0e+00 1.4e+02 1.1e+02  0  0  0  0  6   0  0  0  0  6     0
MatGetRowIJ           50 1.0 1.6451e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering        50 1.0 9.4891e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+02  0  0  0  0  5   0  0  0  0  6     0
MatZeroEntries        52 1.0 6.4778e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       308 1.0 5.7874e-03 1.1 1.63e+07 1.0 0.0e+00 0.0e+00 3.1e+02  0 12  0  0 16   0 12  0  0 17  5580
KSPSetup             100 1.0 8.8215e-05 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve              50 1.0 3.5689e-01 1.0 1.39e+08 1.0 7.2e+02 5.5e+02 8.7e+02 20100 48 11 45  20100 48 11 48   777
PCSetUp              100 1.0 2.8055e-01 1.0 4.73e+07 1.0 0.0e+00 0.0e+00 1.5e+02 15 34  0  0  8  15 34  0  0  8   335
PCSetUpOnBlocks       50 1.0 2.8022e-01 1.0 4.73e+07 1.0 0.0e+00 0.0e+00 1.5e+02 15 34  0  0  8  15 34  0  0  8   335
PCApply              358 1.0 4.2087e-02 1.0 5.37e+07 1.0 0.0e+00 0.0e+00 0.0e+00  2 39  0  0  0   2 39  0  0  0  2537
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    29             29       598928     0
         Vec Scatter    55             55        47740     0
           Index Set   209            209      2012376     0
   IS L to G Mapping     3              3        40104     0
              Matrix    53             53     46427072     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 1.38283e-06
Average time for zero size MPI_Send(): 6.4373e-06
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
| Time:           Thu Feb  3 12:14:14 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=1.83248, Active time=1.71482                                               |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 1         0.0019      0.001936    0.0020      0.001996    0.11     0.12     |
|   compute_sparsity()           1         0.0096      0.009636    0.0131      0.013082    0.56     0.76     |
|   create_dof_constraints()     1         0.0030      0.003026    0.0030      0.003026    0.18     0.18     |
|   distribute_dofs()            1         0.0059      0.005887    0.0145      0.014489    0.34     0.84     |
|   dof_indices()                302309    0.1133      0.000000    0.1133      0.000000    6.61     6.61     |
|   prepare_send_list()          1         0.0000      0.000010    0.0000      0.000010    0.00     0.00     |
|   reinit()                     1         0.0085      0.008455    0.0085      0.008455    0.49     0.49     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         262450    0.1180      0.000000    0.1180      0.000000    6.88     6.88     |
|   compute_face_map()           6400      0.0121      0.000002    0.0225      0.000004    0.71     1.31     |
|   compute_shape_functions()    262450    0.0559      0.000000    0.0559      0.000000    3.26     3.26     |
|   init_face_shape_functions()  650       0.0004      0.000001    0.0004      0.000001    0.02     0.02     |
|   init_shape_functions()       6500      0.0081      0.000001    0.0081      0.000001    0.47     0.47     |
|   inverse_map()                25600     0.0179      0.000001    0.0179      0.000001    1.05     1.05     |
|                                                                                                            |
| GMVIO                                                                                                      |
|   write_nodal_data()           6         0.1127      0.018784    0.1127      0.018784    6.57     6.57     |
|                                                                                                            |
| LocationMap                                                                                                |
|   find()                       32736     0.0102      0.000000    0.0102      0.000000    0.59     0.59     |
|   init()                       5         0.0003      0.000069    0.0003      0.000069    0.02     0.02     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             2         0.0198      0.009887    0.0198      0.009909    1.15     1.16     |
|   renumber_nodes_and_elem()    2         0.0007      0.000359    0.0007      0.000359    0.04     0.04     |
|                                                                                                            |
| MeshCommunication                                                                                          |
|   broadcast_bcs()              1         0.0000      0.000005    0.0000      0.000006    0.00     0.00     |
|   broadcast_mesh()             1         0.0001      0.000068    0.0001      0.000086    0.00     0.01     |
|   compute_hilbert_indices()    3         0.0171      0.005713    0.0171      0.005713    1.00     1.00     |
|   find_global_indices()        3         0.0023      0.000764    0.0203      0.006771    0.13     1.18     |
|   parallel_sort()              3         0.0006      0.000202    0.0007      0.000225    0.04     0.04     |
|                                                                                                            |
| MeshRefinement                                                                                             |
|   _refine_elements()           5         0.0253      0.005062    0.0608      0.012162    1.48     3.55     |
|   add_point()                  32736     0.0198      0.000001    0.0321      0.000001    1.15     1.87     |
|                                                                                                            |
| MetisPartitioner                                                                                           |
|   partition()                  2         0.0116      0.005792    0.0315      0.015752    0.68     1.84     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  7         0.0001      0.000010    0.0001      0.000010    0.00     0.00     |
|   broadcast()                  9         0.0000      0.000003    0.0000      0.000002    0.00     0.00     |
|   gather()                     1         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   max(bool)                    5         0.0011      0.000212    0.0011      0.000212    0.06     0.06     |
|   max(scalar)                  5         0.0001      0.000018    0.0001      0.000018    0.01     0.01     |
|   max(vector)                  3         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   min(vector)                  3         0.0000      0.000005    0.0000      0.000005    0.00     0.00     |
|   probe()                      14        0.0001      0.000008    0.0001      0.000008    0.01     0.01     |
|   receive()                    14        0.0001      0.000007    0.0002      0.000015    0.01     0.01     |
|   send()                       14        0.0001      0.000004    0.0001      0.000004    0.00     0.00     |
|   send_receive()               20        0.0000      0.000002    0.0003      0.000017    0.00     0.02     |
|   sum()                        27        0.0011      0.000042    0.0011      0.000042    0.07     0.07     |
|   wait()                       14        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                            |
| Partitioner                                                                                                |
|   set_node_processor_ids()     2         0.0030      0.001486    0.0030      0.001505    0.17     0.18     |
|   set_parent_processor_ids()   2         0.0012      0.000593    0.0012      0.000593    0.07     0.07     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      50        0.3697      0.007393    0.3697      0.007393    21.56    21.56    |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   50        0.7534      0.015067    1.1208      0.022417    43.93    65.36    |
|   project_vector()             1         0.0097      0.009705    0.0115      0.011464    0.57     0.67     |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        932111    1.7148                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 2 ./ex9-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
