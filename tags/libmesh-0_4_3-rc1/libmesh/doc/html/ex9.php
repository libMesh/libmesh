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
solution is advanced in time with a standard Crank-Nicholson
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
so we need to include the \p TransientImplicitSystem definition.
</div>

<div class ="fragment">
<pre>
        #include "transient_system.h"
        #include "vector_value.h"
        
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
Read the mesh from file.  This is the coarse mesh that will be used
in example 10 to demonstrate adaptive mesh refinement.  Here we will
simply read it in and uniformly refine it 5 times before we compute
with it.
</div>

<div class ="fragment">
<pre>
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
            TransientImplicitSystem & system = 
              equation_systems.add_system&lt;TransientImplicitSystem&gt; ("Convection-Diffusion");
              
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
data type and then use the DataMap object to pass it to
the assemble function.  The DataMap is a convenient way
to encapsulate various data types.
</div>

<div class ="fragment">
<pre>
            RealVectorValue velocity (0.8, 0.8);
        
            equation_systems.data_map.add_data ("velocity", velocity);
            
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
        
        	equation_systems.set_parameter ("time") = time;
        	equation_systems.set_parameter ("dt")   = dt;
        
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
        	  out &lt;&lt;  "..." &lt;&lt; std::endl;
        	  std::cout &lt;&lt; out.str();
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
                TransientImplicitSystem&  system =
        	  equation_systems.get_system&lt;TransientImplicitSystem&gt;("Convection-Diffusion");
        
        	*system.old_local_solution = *system.current_local_solution;
        	
</pre>
</div>
<div class = "comment">
Assemble & solve the linear system
</div>

<div class ="fragment">
<pre>
                equation_systems("Convection-Diffusion").solve();
        	
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
          assert (system_name == "Convection-Diffusion");
            
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
Get a reference to the Convection-Diffusion system object.
</div>

<div class ="fragment">
<pre>
          TransientImplicitSystem & system =
            es.get_system&lt;TransientImplicitSystem&gt; ("Convection-Diffusion");
          
</pre>
</div>
<div class = "comment">
Get a reference to the \p DofMap for this system.
</div>

<div class ="fragment">
<pre>
          const DofMap& dof_map = system.get_dof_map();
          
</pre>
</div>
<div class = "comment">
Get a reference to the solution vector.
</div>

<div class ="fragment">
<pre>
          NumericVector&lt;Number&gt;& solution = *system.solution;
          
</pre>
</div>
<div class = "comment">
A vector to hold the global DOF indices for this element.
</div>

<div class ="fragment">
<pre>
          std::vector&lt;unsigned int&gt; dof_indices;
        
</pre>
</div>
<div class = "comment">
Loop over the active local elements and compute the initial value
of the solution at the element degrees of freedom.  Assign
these initial values to the solution vector.  There is a small
catch, however...  We only want to assign the components that
live on the local processor, hence there will be an if-test
in the loop.
</div>

<div class ="fragment">
<pre>
          const_active_local_elem_iterator       elem_it (mesh.elements_begin());
          const const_active_local_elem_iterator elem_end(mesh.elements_end());
        
          for ( ; elem_it != elem_end; ++elem_it)
            {
              const Elem* elem = *elem_it;
        
              dof_map.dof_indices (elem, dof_indices);
              
</pre>
</div>
<div class = "comment">
For these Lagrange-elements the number
of degrees of freedom should be <= the number
of nodes.
</div>

<div class ="fragment">
<pre>
              assert (dof_indices.size() &lt;= elem-&gt;n_nodes());
        
</pre>
</div>
<div class = "comment">
Loop over the element DOFs, compute the initial
value if the DOF is local to the processor.
</div>

<div class ="fragment">
<pre>
              for (unsigned int i=0; i&lt;dof_indices.size(); i++)
        	if ((dof_indices[i] &gt;= solution.first_local_index()) &&
        	    (dof_indices[i] &lt;  solution.last_local_index()))
        	  {
        	    const Point&  p = elem-&gt;point (i);
        	    const Real    x = p(0);
        	    const Real    y = p(1);
        	    const Real time = 0.;
        	    
        	    solution.set (dof_indices[i],
        			  exact_solution (x,y,time));	    
        	  }	 
            }
          
</pre>
</div>
<div class = "comment">
The initial solution has now been set for the local
solution components.  However, the local matrix assembly
will likely depend on solution components that live on
other processors.  We need to get those components, and
the TransientSystem::update() member will do that
for us.
</div>

<div class ="fragment">
<pre>
          system.update ();
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
</pre>
</div>
<div class = "comment">
It is a good idea to make sure we are assembling
the proper system.
</div>

<div class ="fragment">
<pre>
          assert (system_name == "Convection-Diffusion");
          
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
Get a reference to the Convection-Diffusion system object.
</div>

<div class ="fragment">
<pre>
          TransientImplicitSystem & system =
            es.get_system&lt;TransientImplicitSystem&gt; ("Convection-Diffusion");
          
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
          AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
          
</pre>
</div>
<div class = "comment">
A Gauss quadrature rule for numerical integration.
Let the \p FEType object decide what order rule is appropriate.
</div>

<div class ="fragment">
<pre>
          QGauss qrule (dim, fe_type.default_quadrature_order());
        
</pre>
</div>
<div class = "comment">
Tell the finite element object to use our quadrature rule.
</div>

<div class ="fragment">
<pre>
          fe-&gt;attach_quadrature_rule (&qrule);
          
</pre>
</div>
<div class = "comment">
Here we define some references to cell-specific data that
will be used to assemble the linear system.  We will start
with the element Jacobian * quadrature weight at each integration point.   
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt;& JxW = fe-&gt;get_JxW();
          
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
The element shape function gradients evaluated at the quadrature
points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi = fe-&gt;get_dphi();
          
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
            es.data_map.get_data&lt;RealVectorValue&gt; ("velocity");
        
          const Real dt = es.parameter   ("dt");
          const Real time = es.parameter ("time");
        
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
via the penalty method. The penalty method used here
is equivalent (for Lagrange basis functions) to lumping
the matrix resulting from the L2 projection penalty
approach introduced in example 3.

<br><br>The following loops over the sides of the element.
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
</div>

<div class ="fragment">
<pre>
                        const Real xf = side-&gt;point(ns)(0);
        		const Real yf = side-&gt;point(ns)(1);
        		  
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
The boundary value.
</div>

<div class ="fragment">
<pre>
                        const Real value = exact_solution(xf, yf, time);
        
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
                              Ke(n,n) += penalty;
        		  		  
</pre>
</div>
<div class = "comment">
Right-hand-side contribution.
</div>

<div class ="fragment">
<pre>
                              Fe(n) += penalty*value;
        		    }
        	      } 
        	  }
        
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
  #include <FONT COLOR="#BC8F8F"><B>&quot;mesh_refinement.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;gmv_io.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;equation_systems.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;fe.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dof_map.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;sparse_matrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;numeric_vector.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_matrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_vector.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;o_string_stream.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;transient_system.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;vector_value.h&quot;</FONT></B>
  
  <FONT COLOR="#228B22"><B>void</FONT></B> assemble_cd (EquationSystems&amp; es,
  		  <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name);
  
  <FONT COLOR="#228B22"><B>void</FONT></B> init_cd (EquationSystems&amp; es,
  	      <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name);
  
  Real exact_solution (<FONT COLOR="#228B22"><B>const</FONT></B> Real x,
  		     <FONT COLOR="#228B22"><B>const</FONT></B> Real y,
  		     <FONT COLOR="#228B22"><B>const</FONT></B> Real t);
  
  <FONT COLOR="#228B22"><B>int</FONT></B> main (<FONT COLOR="#228B22"><B>int</FONT></B> argc, <FONT COLOR="#228B22"><B>char</FONT></B>** argv)
  {
    libMesh::init (argc, argv);
  
    {    
      Mesh mesh (2);
          
      mesh.read (<FONT COLOR="#BC8F8F"><B>&quot;../ex10/mesh.xda&quot;</FONT></B>);
      
      MeshRefinement mesh_refinement (mesh);
      
      mesh_refinement.uniformly_refine (5);
      
      mesh.print_info();
      
      EquationSystems equation_systems (mesh);
      
      TransientImplicitSystem &amp; system = 
        equation_systems.add_system&lt;TransientImplicitSystem&gt; (<FONT COLOR="#BC8F8F"><B>&quot;Convection-Diffusion&quot;</FONT></B>);
        
      system.add_variable (<FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>, FIRST);
        
      system.attach_assemble_function (assemble_cd);
      system.attach_init_function (init_cd);
        
      equation_systems.init ();
        
      equation_systems.print_info();
        
      GMVIO(mesh).write_equation_systems (<FONT COLOR="#BC8F8F"><B>&quot;out_000.gmv&quot;</FONT></B>,
  					equation_systems);
      
      RealVectorValue velocity (0.8, 0.8);
  
      equation_systems.data_map.add_data (<FONT COLOR="#BC8F8F"><B>&quot;velocity&quot;</FONT></B>, velocity);
      
      <FONT COLOR="#228B22"><B>const</FONT></B> Real dt = 0.025;
      Real time     = 0.;
      
      <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> t_step = 0; t_step &lt; 50; t_step++)
        {
  	time += dt;
  
  	equation_systems.set_parameter (<FONT COLOR="#BC8F8F"><B>&quot;time&quot;</FONT></B>) = time;
  	equation_systems.set_parameter (<FONT COLOR="#BC8F8F"><B>&quot;dt&quot;</FONT></B>)   = dt;
  
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; Solving time step &quot;</FONT></B>;
  	
  	{
  	  OStringStream out;
  
  	  OSSInt(out,2,t_step);
  	  out &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;, time=&quot;</FONT></B>;
  	  OSSRealzeroleft(out,6,3,time);
  	  out &lt;&lt;  <FONT COLOR="#BC8F8F"><B>&quot;...&quot;</FONT></B> &lt;&lt; std::endl;
  	  std::cout &lt;&lt; out.str();
  	}
  	
  	TransientImplicitSystem&amp;  system =
  	  equation_systems.get_system&lt;TransientImplicitSystem&gt;(<FONT COLOR="#BC8F8F"><B>&quot;Convection-Diffusion&quot;</FONT></B>);
  
  	*system.old_local_solution = *system.current_local_solution;
  	
  	equation_systems(<FONT COLOR="#BC8F8F"><B>&quot;Convection-Diffusion&quot;</FONT></B>).solve();
  	
  	<B><FONT COLOR="#A020F0">if</FONT></B> ( (t_step+1)%10 == 0)
  	  {
  	    OStringStream file_name;
  
  	    file_name &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;out_&quot;</FONT></B>;
  	    OSSRealzeroright(file_name,3,0,t_step+1);
  	    file_name &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;.gmv&quot;</FONT></B>;
  
  	    GMVIO(mesh).write_equation_systems (file_name.str(),
  						equation_systems);
  	  }
        }
    }
  
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close ();
  }
  
  <FONT COLOR="#228B22"><B>void</FONT></B> init_cd (EquationSystems&amp; es,
  	      <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name)
  {
    assert (system_name == <FONT COLOR="#BC8F8F"><B>&quot;Convection-Diffusion&quot;</FONT></B>);
      
    <FONT COLOR="#228B22"><B>const</FONT></B> Mesh&amp; mesh = es.get_mesh();
    
    TransientImplicitSystem &amp; system =
      es.get_system&lt;TransientImplicitSystem&gt; (<FONT COLOR="#BC8F8F"><B>&quot;Convection-Diffusion&quot;</FONT></B>);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
    
    NumericVector&lt;Number&gt;&amp; solution = *system.solution;
    
    std::vector&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt; dof_indices;
  
    const_active_local_elem_iterator       elem_it (mesh.elements_begin());
    <FONT COLOR="#228B22"><B>const</FONT></B> const_active_local_elem_iterator elem_end(mesh.elements_end());
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; elem_it != elem_end; ++elem_it)
      {
        <FONT COLOR="#228B22"><B>const</FONT></B> Elem* elem = *elem_it;
  
        dof_map.dof_indices (elem, dof_indices);
        
        assert (dof_indices.size() &lt;= elem-&gt;n_nodes());
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;dof_indices.size(); i++)
  	<B><FONT COLOR="#A020F0">if</FONT></B> ((dof_indices[i] &gt;= solution.first_local_index()) &amp;&amp;
  	    (dof_indices[i] &lt;  solution.last_local_index()))
  	  {
  	    <FONT COLOR="#228B22"><B>const</FONT></B> Point&amp;  p = elem-&gt;point (i);
  	    <FONT COLOR="#228B22"><B>const</FONT></B> Real    x = p(0);
  	    <FONT COLOR="#228B22"><B>const</FONT></B> Real    y = p(1);
  	    <FONT COLOR="#228B22"><B>const</FONT></B> Real time = 0.;
  	    
  	    solution.set (dof_indices[i],
  			  exact_solution (x,y,time));	    
  	  }	 
      }
    
    system.update ();
  }
  
  
  
  <FONT COLOR="#228B22"><B>void</FONT></B> assemble_cd (EquationSystems&amp; es,
  		  <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name)
  {
    assert (system_name == <FONT COLOR="#BC8F8F"><B>&quot;Convection-Diffusion&quot;</FONT></B>);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> Mesh&amp; mesh = es.get_mesh();
    
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = mesh.mesh_dimension();
    
    TransientImplicitSystem &amp; system =
      es.get_system&lt;TransientImplicitSystem&gt; (<FONT COLOR="#BC8F8F"><B>&quot;Convection-Diffusion&quot;</FONT></B>);
    
    FEType fe_type = system.variable_type(0);
    
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    QGauss qrule (dim, fe_type.default_quadrature_order());
  
    fe-&gt;attach_quadrature_rule (&amp;qrule);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
    
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
    
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
    
    <FONT COLOR="#228B22"><B>const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
  
    DenseMatrix&lt;Number&gt; Ke;
    DenseVector&lt;Number&gt; Fe;
    
    std::vector&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt; dof_indices;
  
    <FONT COLOR="#228B22"><B>const</FONT></B> RealVectorValue velocity =
      es.data_map.get_data&lt;RealVectorValue&gt; (<FONT COLOR="#BC8F8F"><B>&quot;velocity&quot;</FONT></B>);
  
    <FONT COLOR="#228B22"><B>const</FONT></B> Real dt = es.parameter   (<FONT COLOR="#BC8F8F"><B>&quot;dt&quot;</FONT></B>);
    <FONT COLOR="#228B22"><B>const</FONT></B> Real time = es.parameter (<FONT COLOR="#BC8F8F"><B>&quot;time&quot;</FONT></B>);
  
    const_active_local_elem_iterator           el (mesh.elements_begin());
    <FONT COLOR="#228B22"><B>const</FONT></B> const_active_local_elem_iterator end_el (mesh.elements_end());
    
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {    
        <FONT COLOR="#228B22"><B>const</FONT></B> Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        fe-&gt;reinit (elem);
        
        Ke.resize (dof_indices.size(),
  		 dof_indices.size());
  
        Fe.resize (dof_indices.size());
        
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
  	{
  	  Number   u_old = 0.;
  	  Gradient grad_u_old;
  	  
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> l=0; l&lt;phi.size(); l++)
  	    {
  	      u_old      += phi[l][qp]*system.old_solution  (dof_indices[l]);
  	      
  	      grad_u_old.add_scaled (dphi[l][qp],system.old_solution (dof_indices[l]));
  	    }
  	  
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi.size(); i++)
  	    {
  	      Fe(i) += JxW[qp]*(
  				u_old*phi[i][qp] +
  				-.5*dt*(
  					(grad_u_old*velocity)*phi[i][qp] +
  					
  					0.01*(grad_u_old*dphi[i][qp]))     
  				);
  	      
  	      <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j&lt;phi.size(); j++)
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
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> s=0; s&lt;elem-&gt;n_sides(); s++)
  	<B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(s) == NULL)
  	  {
  	    AutoPtr&lt;Elem&gt; side (elem-&gt;build_side(s));
  	    
  	    
  	    <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> ns=0; ns&lt;side-&gt;n_nodes(); ns++)
  	      {
  		<FONT COLOR="#228B22"><B>const</FONT></B> Real xf = side-&gt;point(ns)(0);
  		<FONT COLOR="#228B22"><B>const</FONT></B> Real yf = side-&gt;point(ns)(1);
  		  
  		<FONT COLOR="#228B22"><B>const</FONT></B> Real penalty = 1.e10;
  		  
  		<FONT COLOR="#228B22"><B>const</FONT></B> Real value = exact_solution(xf, yf, time);
  
  		<B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n=0; n&lt;elem-&gt;n_nodes(); n++)
  		  <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;node(n) == side-&gt;node(ns))
  		    {
  		      Ke(n,n) += penalty;
  		  		  
  		      Fe(n) += penalty*value;
  		    }
  	      } 
  	  }
  
        system.matrix-&gt;add_matrix (Ke, dof_indices);
        system.rhs-&gt;add_vector    (Fe, dof_indices);
      }
    
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Compiling C++ (in debug mode) exact_solution.C...
Linking ex9...
/home/peterson/code/libmesh/contrib/tecplot/lib/i686-pc-linux-gnu/tecio.a(tecxxx.o)(.text+0x1a7): In function `tecini':
: the use of `mktemp' is dangerous, better use `mkstemp'
***************************************************************
* Running Example  ./ex9
***************************************************************
 
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=6273
  n_elem()=13650
   n_local_elem()=13650
   n_active_elem()=10240
  n_subdomains()=1
  n_processors()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Convection-Diffusion"
    Type "TransientImplicit"
    Variables="u" 
    Finite Element Types="0" 
    Approximation Orders="1" 
    n_dofs()=6273
    n_local_dofs()=6273
    n_constrained_dofs()=0
    n_vectors()=1
  n_parameters()=2
   Parameters:
    "linear solver maximum iterations"=5000
    "linear solver tolerance"=1e-12

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

 ---------------------------------------------------------------------------- 
| Reference count information                                                |
 ---------------------------------------------------------------------------- 
| 12SparseMatrixIdE reference count information:
| Creations:    1
| Destructions: 1
| 13NumericVectorIdE reference count information:
| Creations:    5
| Destructions: 5
| 21LinearSolverInterfaceIdE reference count information:
| Creations:    1
| Destructions: 1
| 4Elem reference count information:
| Creations:    70112
| Destructions: 70112
| 4Node reference count information:
| Creations:    6273
| Destructions: 6273
| 5QBase reference count information:
| Creations:    100
| Destructions: 100
| 6DofMap reference count information:
| Creations:    1
| Destructions: 1
| 6FEBase reference count information:
| Creations:    50
| Destructions: 50
| 6System reference count information:
| Creations:    1
| Destructions: 1
 ---------------------------------------------------------------------------- 
 
***************************************************************
* Done Running Example  ./ex9
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
