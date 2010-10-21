<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex19",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 19 - Solving the 2D Young Laplace Problem using nonlinear solvers</h1>

<br><br>This example shows how to use the NonlinearImplicitSystem class
to efficiently solve nonlinear problems in parallel.

<br><br>In nonlinear systems, we aim at finding x that satisfy R(x) = 0. 
In nonlinear finite element analysis, the residual is typically 
of the form R(x) = K(x)*x - f, with K(x) the system matrix and f 
the "right-hand-side". The NonlinearImplicitSystem class expects  
two callback functions to compute the residual R and its Jacobian 
for the Newton iterations. Here, we just approximate 
the true Jacobian by K(x).

<br><br>You can turn on preconditining of the matrix free system using the
jacobian by passing "-pre" on the command line.  Currently this only
work with Petsc so this isn't used by using "make run"

<br><br>This example also runs with the experimental Trilinos NOX solvers by specifying 
the --use-trilinos command line argument.
 

<br><br>

<br><br>C++ include files that we need
</div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
        #include &lt;algorithm&gt;
        #include &lt;cmath&gt;
        
</pre>
</div>
<div class = "comment">
Various include files needed for the mesh & solver functionality.
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
        #include "elem.h"
        #include "string_to_enum.h"
        #include "getpot.h"
        
</pre>
</div>
<div class = "comment">
The nonlinear solver and system we will be using
</div>

<div class ="fragment">
<pre>
        #include "nonlinear_solver.h"
        #include "nonlinear_implicit_system.h"
        
</pre>
</div>
<div class = "comment">
Necessary for programmatically setting petsc options
</div>

<div class ="fragment">
<pre>
        #ifdef LIBMESH_HAVE_PETSC
        #include &lt;petsc.h&gt;
        #endif
        
</pre>
</div>
<div class = "comment">
A reference to our equation system
</div>

<div class ="fragment">
<pre>
        EquationSystems *_equation_system = NULL;
        
</pre>
</div>
<div class = "comment">
Let-s define the physical parameters of the equation
</div>

<div class ="fragment">
<pre>
        const Real kappa = 1.;
        const Real sigma = 0.2;
        
        
</pre>
</div>
<div class = "comment">
This function computes the Jacobian K(x)
</div>

<div class ="fragment">
<pre>
        void compute_jacobian (const NumericVector&lt;Number&gt;& soln,
        		       SparseMatrix&lt;Number&gt;&  jacobian)
        {
</pre>
</div>
<div class = "comment">
Get a reference to the equation system.
</div>

<div class ="fragment">
<pre>
          EquationSystems &es = *_equation_system;
        
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
Get a reference to the NonlinearImplicitSystem we are solving
</div>

<div class ="fragment">
<pre>
          NonlinearImplicitSystem& system = 
            es.get_system&lt;NonlinearImplicitSystem&gt;("Laplace-Young");
          
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
A 5th order Gauss quadrature rule for numerical integration.
</div>

<div class ="fragment">
<pre>
          QGauss qrule (dim, FIFTH);
        
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
Define data structures to contain the Jacobian element matrix.
Following basic finite element terminology we will denote these
"Ke". More detail is in example 3.
</div>

<div class ="fragment">
<pre>
          DenseMatrix&lt;Number&gt; Ke;
        
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
Now we will loop over all the elements in the mesh.
We will compute the element Jacobian contribution.
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
Zero the element Jacobian before
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
                   
</pre>
</div>
<div class = "comment">
Now we will build the element Jacobian.  This involves
a double loop to integrate the test funcions (i) against
the trial functions (j). Note that the Jacobian depends
on the current solution x, which we access using the soln
vector.

<br><br></div>

<div class ="fragment">
<pre>
              for (unsigned int qp=0; qp&lt;qrule.n_points(); qp++)
        	{
        	  RealGradient grad_u = 0;
            
        	  for (unsigned int i=0; i&lt;phi.size(); i++)
        	    grad_u += dphi[i][qp]*soln(dof_indices[i]);
        	  
        	  const Real K = 1./std::sqrt(1. + grad_u*grad_u);
        	  
        	  for (unsigned int i=0; i&lt;phi.size(); i++)
        	    for (unsigned int j=0; j&lt;phi.size(); j++)
        	      Ke(i,j) += JxW[qp]*(
        				  K*(dphi[i][qp]*dphi[j][qp]) +
        				  kappa*phi[i][qp]*phi[j][qp]
        				  );
        	}
              
              dof_map.constrain_element_matrix (Ke, dof_indices);
              
</pre>
</div>
<div class = "comment">
Add the element matrix to the system Jacobian.
</div>

<div class ="fragment">
<pre>
              jacobian.add_matrix (Ke, dof_indices);
            }
        
</pre>
</div>
<div class = "comment">
That's it.
</div>

<div class ="fragment">
<pre>
        }
        
        
</pre>
</div>
<div class = "comment">
Here we compute the residual R(x) = K(x)*x - f. The current solution
x is passed in the soln vector
</div>

<div class ="fragment">
<pre>
        void compute_residual (const NumericVector&lt;Number&gt;& soln,
        		       NumericVector&lt;Number&gt;& residual)
        {
          EquationSystems &es = *_equation_system;
        
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
          libmesh_assert (dim == 2);
        
</pre>
</div>
<div class = "comment">
Get a reference to the NonlinearImplicitSystem we are solving
</div>

<div class ="fragment">
<pre>
          NonlinearImplicitSystem& system = 
            es.get_system&lt;NonlinearImplicitSystem&gt;("Laplace-Young");
          
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
A 5th order Gauss quadrature rule for numerical integration.
</div>

<div class ="fragment">
<pre>
          QGauss qrule (dim, FIFTH);
        
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
Declare a special finite element object for
boundary integration.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
        	      
</pre>
</div>
<div class = "comment">
Boundary integration requires one quadraure rule,
with dimensionality one less than the dimensionality
of the element.
</div>

<div class ="fragment">
<pre>
          QGauss qface(dim-1, FIFTH);
          
</pre>
</div>
<div class = "comment">
Tell the finte element object to use our
quadrature rule.
</div>

<div class ="fragment">
<pre>
          fe_face-&gt;attach_quadrature_rule (&qface);
        
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
Define data structures to contain the resdual contributions
</div>

<div class ="fragment">
<pre>
          DenseVector&lt;Number&gt; Re;
        
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
Now we will loop over all the elements in the mesh.
We will compute the element residual.
</div>

<div class ="fragment">
<pre>
          residual.zero();
        
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
We use the resize member here because
the number of degrees of freedom might have changed from
the last element.  Note that this will be the case if the
element type is different (i.e. the last element was a
triangle, now we are on a quadrilateral).
</div>

<div class ="fragment">
<pre>
              Re.resize (dof_indices.size());
              
</pre>
</div>
<div class = "comment">
Now we will build the residual. This involves
the construction of the matrix K and multiplication of it
with the current solution x. We rearrange this into two loops: 
In the first, we calculate only the contribution of  
K_ij*x_j which is independent of the row i. In the second loops,
we multiply with the row-dependent part and add it to the element
residual.


<br><br></div>

<div class ="fragment">
<pre>
              for (unsigned int qp=0; qp&lt;qrule.n_points(); qp++)
        	{
        	  Real u = 0;
        	  RealGradient grad_u = 0;
        	  
        	  for (unsigned int j=0; j&lt;phi.size(); j++)
        	    {
        	      u      += phi[j][qp]*soln(dof_indices[j]);
        	      grad_u += dphi[j][qp]*soln(dof_indices[j]);
        	    }
        	  
        	  const Real K = 1./std::sqrt(1. + grad_u*grad_u);
        	  
        	  for (unsigned int i=0; i&lt;phi.size(); i++)
        	    Re(i) += JxW[qp]*(
        			      K*(dphi[i][qp]*grad_u) +
        			      kappa*phi[i][qp]*u
        			      );
        	}
        
</pre>
</div>
<div class = "comment">
At this point the interior element integration has
been completed.  However, we have not yet addressed
boundary conditions.
      

<br><br>The following loops over the sides of the element.
If the element has no neighbor on a side then that
side MUST live on a boundary of the domain.
</div>

<div class ="fragment">
<pre>
              for (unsigned int side=0; side&lt;elem-&gt;n_sides(); side++)
        	if (elem-&gt;neighbor(side) == NULL)
        	  {
</pre>
</div>
<div class = "comment">
The value of the shape functions at the quadrature
points.
</div>

<div class ="fragment">
<pre>
                    const std::vector&lt;std::vector&lt;Real&gt; &gt;&  phi_face = fe_face-&gt;get_phi();
        
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
Compute the shape function values on the element face.
</div>

<div class ="fragment">
<pre>
                    fe_face-&gt;reinit(elem, side);
        
</pre>
</div>
<div class = "comment">
Loop over the face quadrature points for integration.
</div>

<div class ="fragment">
<pre>
                    for (unsigned int qp=0; qp&lt;qface.n_points(); qp++)
        	      {
</pre>
</div>
<div class = "comment">
This is the right-hand-side contribution (f),
which has to be subtracted from the current residual
</div>

<div class ="fragment">
<pre>
                        for (unsigned int i=0; i&lt;phi_face.size(); i++)
        		  Re(i) -= JxW_face[qp]*sigma*phi_face[i][qp];
        	      } 
        	  }
              
              dof_map.constrain_element_vector (Re, dof_indices);
              residual.add_vector (Re, dof_indices);
            }
        
</pre>
</div>
<div class = "comment">
That's it.  
</div>

<div class ="fragment">
<pre>
        }
        
        
        
</pre>
</div>
<div class = "comment">
Begin the main program.
</div>

<div class ="fragment">
<pre>
        int main (int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Initialize libMesh and any dependent libaries, like in example 2.
</div>

<div class ="fragment">
<pre>
          LibMeshInit init (argc, argv);
         
</pre>
</div>
<div class = "comment">
Braces are used to force object scope, like in example 2
</div>

<div class ="fragment">
<pre>
          {
</pre>
</div>
<div class = "comment">
Create a GetPot object to parse the command line
</div>

<div class ="fragment">
<pre>
            GetPot command_line (argc, argv);
            
</pre>
</div>
<div class = "comment">
Check for proper calling arguments.
</div>

<div class ="fragment">
<pre>
            if (argc &lt; 3)
              {
                if (libMesh::processor_id() == 0)
        	  std::cerr &lt;&lt; "Usage:\n"
        		    &lt;&lt;"\t " &lt;&lt; argv[0] &lt;&lt; " -r 2"
        		    &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
This handy function will print the file name, line number,
and then abort.  Currrently the library does not use C++
exception handling.
</div>

<div class ="fragment">
<pre>
                error();
              }
            
</pre>
</div>
<div class = "comment">
Brief message to the user regarding the program name
and command line arguments.
</div>

<div class ="fragment">
<pre>
            else 
              {
        	std::cout &lt;&lt; "Running " &lt;&lt; argv[0];
        	
        	for (int i=1; i&lt;argc; i++)
        	  std::cout &lt;&lt; " " &lt;&lt; argv[i];
        	
        	std::cout &lt;&lt; std::endl &lt;&lt; std::endl;
              }
            
        
</pre>
</div>
<div class = "comment">
The dimension of our problem 
</div>

<div class ="fragment">
<pre>
            const unsigned int dim = 2;
            
</pre>
</div>
<div class = "comment">
Read number of refinements 
</div>

<div class ="fragment">
<pre>
            int nr = 2;
            if ( command_line.search(1, "-r") )
              nr = command_line.next(nr);
            
</pre>
</div>
<div class = "comment">
Read FE order from command line
</div>

<div class ="fragment">
<pre>
            std::string order = "FIRST"; 
            if ( command_line.search(2, "-Order", "-o") )
              order = command_line.next(order);
        
</pre>
</div>
<div class = "comment">
Read FE Family from command line
</div>

<div class ="fragment">
<pre>
            std::string family = "LAGRANGE"; 
            if ( command_line.search(2, "-FEFamily", "-f") )
              family = command_line.next(family);
            
</pre>
</div>
<div class = "comment">
Cannot use dicontinuous basis.
</div>

<div class ="fragment">
<pre>
            if ((family == "MONOMIAL") || (family == "XYZ"))
              {
        	std::cout &lt;&lt; "ex19 currently requires a C^0 (or higher) FE basis." &lt;&lt; std::endl;
        	error();
              }
        
            if ( command_line.search(1, "-pre") )
              {
        #ifdef LIBMESH_HAVE_PETSC
</pre>
</div>
<div class = "comment">
Use the jacobian for preconditioning.
</div>

<div class ="fragment">
<pre>
                PetscOptionsSetValue("-snes_mf_operator",PETSC_NULL);
        #else
                std::cerr&lt;&lt;"Must be using PetsC to use jacobian based preconditioning"&lt;&lt;std::endl;
        
</pre>
</div>
<div class = "comment">
returning zero so that "make run" won't fail if we ever enable this capability there.
</div>

<div class ="fragment">
<pre>
                return 0;
        #endif //LIBMESH_HAVE_PETSC
              }  
              
</pre>
</div>
<div class = "comment">
Create a mesh with user-defined dimension.
</div>

<div class ="fragment">
<pre>
            Mesh mesh (dim);    
            mesh.read ("lshaped.xda");
        
            if (order != "FIRST")
              mesh.all_second_order();
        
            MeshRefinement(mesh).uniformly_refine(nr);
        
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
            _equation_system = &equation_systems;
            
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
Creates a system named "Laplace-Young"
</div>

<div class ="fragment">
<pre>
              NonlinearImplicitSystem& system =
        	equation_systems.add_system&lt;NonlinearImplicitSystem&gt; ("Laplace-Young");
        
        
</pre>
</div>
<div class = "comment">
Here we specify the tolerance for the nonlinear solver and 
the maximum of nonlinear iterations. 
</div>

<div class ="fragment">
<pre>
              equation_systems.parameters.set&lt;Real&gt;         ("nonlinear solver tolerance")          = 1.e-12;
              equation_systems.parameters.set&lt;unsigned int&gt; ("nonlinear solver maximum iterations") = 50;
        
              
</pre>
</div>
<div class = "comment">
Adds the variable "u" to "Laplace-Young".  "u"
will be approximated using second-order approximation.
</div>

<div class ="fragment">
<pre>
              system.add_variable("u",
        			  Utility::string_to_enum&lt;Order&gt;   (order),
        			  Utility::string_to_enum&lt;FEFamily&gt;(family));
        
</pre>
</div>
<div class = "comment">
Give the system a pointer to the functions that update 
the residual and Jacobian.
      

<br><br></div>

<div class ="fragment">
<pre>
              system.nonlinear_solver-&gt;residual = compute_residual;
              system.nonlinear_solver-&gt;jacobian = compute_jacobian;
        
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
Prints information about the system to the screen.
</div>

<div class ="fragment">
<pre>
              equation_systems.print_info();
            }
            
</pre>
</div>
<div class = "comment">
Solve the system "Laplace-Young"
</div>

<div class ="fragment">
<pre>
            equation_systems.get_system("Laplace-Young").solve();
        
</pre>
</div>
<div class = "comment">
After solving the system write the solution
</div>

<div class ="fragment">
<pre>
            GMVIO (mesh).write_equation_systems ("out.gmv", 
        					 equation_systems);
          }
          
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

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
   
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;cmath&gt;
  
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
  #include <FONT COLOR="#BC8F8F"><B>&quot;elem.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;string_to_enum.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;getpot.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;nonlinear_solver.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;nonlinear_implicit_system.h&quot;</FONT></B>
  
  #ifdef LIBMESH_HAVE_PETSC
  #include &lt;petsc.h&gt;
  #endif
  
  EquationSystems *_equation_system = NULL;
  
  <FONT COLOR="#228B22"><B>const</FONT></B> Real kappa = 1.;
  <FONT COLOR="#228B22"><B>const</FONT></B> Real sigma = 0.2;
  
  
  <FONT COLOR="#228B22"><B>void</FONT></B> compute_jacobian (<FONT COLOR="#228B22"><B>const</FONT></B> NumericVector&lt;Number&gt;&amp; soln,
  		       SparseMatrix&lt;Number&gt;&amp;  jacobian)
  {
    EquationSystems &amp;es = *_equation_system;
  
    <FONT COLOR="#228B22"><B>const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = mesh.mesh_dimension();
  
    NonlinearImplicitSystem&amp; system = 
      es.get_system&lt;NonlinearImplicitSystem&gt;(<FONT COLOR="#BC8F8F"><B>&quot;Laplace-Young&quot;</FONT></B>);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
  
    FEType fe_type = dof_map.variable_type(0);
  
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    QGauss qrule (dim, FIFTH);
  
    fe-&gt;attach_quadrature_rule (&amp;qrule);
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
    
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    DenseMatrix&lt;Number&gt; Ke;
  
    std::vector&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt; dof_indices;
  
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    <FONT COLOR="#228B22"><B>const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        <FONT COLOR="#228B22"><B>const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
  
        fe-&gt;reinit (elem);
  
        Ke.resize (dof_indices.size(),
  		 dof_indices.size());
             
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
  	{
  	  RealGradient grad_u = 0;
      
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi.size(); i++)
  	    grad_u += dphi[i][qp]*soln(dof_indices[i]);
  	  
  	  <FONT COLOR="#228B22"><B>const</FONT></B> Real K = 1./std::sqrt(1. + grad_u*grad_u);
  	  
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi.size(); i++)
  	    <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j&lt;phi.size(); j++)
  	      Ke(i,j) += JxW[qp]*(
  				  K*(dphi[i][qp]*dphi[j][qp]) +
  				  kappa*phi[i][qp]*phi[j][qp]
  				  );
  	}
        
        dof_map.constrain_element_matrix (Ke, dof_indices);
        
        jacobian.add_matrix (Ke, dof_indices);
      }
  
  }
  
  
  <FONT COLOR="#228B22"><B>void</FONT></B> compute_residual (<FONT COLOR="#228B22"><B>const</FONT></B> NumericVector&lt;Number&gt;&amp; soln,
  		       NumericVector&lt;Number&gt;&amp; residual)
  {
    EquationSystems &amp;es = *_equation_system;
  
    <FONT COLOR="#228B22"><B>const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = mesh.mesh_dimension();
    libmesh_assert (dim == 2);
  
    NonlinearImplicitSystem&amp; system = 
      es.get_system&lt;NonlinearImplicitSystem&gt;(<FONT COLOR="#BC8F8F"><B>&quot;Laplace-Young&quot;</FONT></B>);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
  
    FEType fe_type = dof_map.variable_type(0);
  
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    QGauss qrule (dim, FIFTH);
  
    fe-&gt;attach_quadrature_rule (&amp;qrule);
  
    AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
  	      
    QGauss qface(dim-1, FIFTH);
    
    fe_face-&gt;attach_quadrature_rule (&amp;qface);
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
    
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    DenseVector&lt;Number&gt; Re;
  
    std::vector&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt; dof_indices;
  
    residual.zero();
  
    MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
    <FONT COLOR="#228B22"><B>const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        <FONT COLOR="#228B22"><B>const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
  
        fe-&gt;reinit (elem);
  
        Re.resize (dof_indices.size());
        
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
  	{
  	  Real u = 0;
  	  RealGradient grad_u = 0;
  	  
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j&lt;phi.size(); j++)
  	    {
  	      u      += phi[j][qp]*soln(dof_indices[j]);
  	      grad_u += dphi[j][qp]*soln(dof_indices[j]);
  	    }
  	  
  	  <FONT COLOR="#228B22"><B>const</FONT></B> Real K = 1./std::sqrt(1. + grad_u*grad_u);
  	  
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi.size(); i++)
  	    Re(i) += JxW[qp]*(
  			      K*(dphi[i][qp]*grad_u) +
  			      kappa*phi[i][qp]*u
  			      );
  	}
  
        
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> side=0; side&lt;elem-&gt;n_sides(); side++)
  	<B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(side) == NULL)
  	  {
  	    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp;  phi_face = fe_face-&gt;get_phi();
  
  	    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
  
  	    fe_face-&gt;reinit(elem, side);
  
  	    <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp&lt;qface.n_points(); qp++)
  	      {
  		<B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi_face.size(); i++)
  		  Re(i) -= JxW_face[qp]*sigma*phi_face[i][qp];
  	      } 
  	  }
        
        dof_map.constrain_element_vector (Re, dof_indices);
        residual.add_vector (Re, dof_indices);
      }
  
  }
  
  
  
  <FONT COLOR="#228B22"><B>int</FONT></B> main (<FONT COLOR="#228B22"><B>int</FONT></B> argc, <FONT COLOR="#228B22"><B>char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
   
    {
      GetPot command_line (argc, argv);
      
      <B><FONT COLOR="#A020F0">if</FONT></B> (argc &lt; 3)
        {
          <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
  	  std::cerr &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Usage:\n&quot;</FONT></B>
  		    &lt;&lt;<FONT COLOR="#BC8F8F"><B>&quot;\t &quot;</FONT></B> &lt;&lt; argv[0] &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; -r 2&quot;</FONT></B>
  		    &lt;&lt; std::endl;
  
  	error();
        }
      
      <B><FONT COLOR="#A020F0">else</FONT></B> 
        {
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Running &quot;</FONT></B> &lt;&lt; argv[0];
  	
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>int</FONT></B> i=1; i&lt;argc; i++)
  	  std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; &quot;</FONT></B> &lt;&lt; argv[i];
  	
  	std::cout &lt;&lt; std::endl &lt;&lt; std::endl;
        }
      
  
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = 2;
      
      <FONT COLOR="#228B22"><B>int</FONT></B> nr = 2;
      <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <FONT COLOR="#BC8F8F"><B>&quot;-r&quot;</FONT></B>) )
        nr = command_line.next(nr);
      
      std::string order = <FONT COLOR="#BC8F8F"><B>&quot;FIRST&quot;</FONT></B>; 
      <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(2, <FONT COLOR="#BC8F8F"><B>&quot;-Order&quot;</FONT></B>, <FONT COLOR="#BC8F8F"><B>&quot;-o&quot;</FONT></B>) )
        order = command_line.next(order);
  
      std::string family = <FONT COLOR="#BC8F8F"><B>&quot;LAGRANGE&quot;</FONT></B>; 
      <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(2, <FONT COLOR="#BC8F8F"><B>&quot;-FEFamily&quot;</FONT></B>, <FONT COLOR="#BC8F8F"><B>&quot;-f&quot;</FONT></B>) )
        family = command_line.next(family);
      
      <B><FONT COLOR="#A020F0">if</FONT></B> ((family == <FONT COLOR="#BC8F8F"><B>&quot;MONOMIAL&quot;</FONT></B>) || (family == <FONT COLOR="#BC8F8F"><B>&quot;XYZ&quot;</FONT></B>))
        {
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;ex19 currently requires a C^0 (or higher) FE basis.&quot;</FONT></B> &lt;&lt; std::endl;
  	error();
        }
  
      <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <FONT COLOR="#BC8F8F"><B>&quot;-pre&quot;</FONT></B>) )
        {
  #ifdef LIBMESH_HAVE_PETSC
          PetscOptionsSetValue(<FONT COLOR="#BC8F8F"><B>&quot;-snes_mf_operator&quot;</FONT></B>,PETSC_NULL);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
          std::cerr&lt;&lt;<FONT COLOR="#BC8F8F"><B>&quot;Must be using PetsC to use jacobian based preconditioning&quot;</FONT></B>&lt;&lt;std::endl;
  
          <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  #endif <I><FONT COLOR="#B22222">//LIBMESH_HAVE_PETSC
</FONT></I>        }  
        
      Mesh mesh (dim);    
      mesh.read (<FONT COLOR="#BC8F8F"><B>&quot;lshaped.xda&quot;</FONT></B>);
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (order != <FONT COLOR="#BC8F8F"><B>&quot;FIRST&quot;</FONT></B>)
        mesh.all_second_order();
  
      MeshRefinement(mesh).uniformly_refine(nr);
  
      mesh.print_info();    
      
      EquationSystems equation_systems (mesh);
      _equation_system = &amp;equation_systems;
      
      {
        NonlinearImplicitSystem&amp; system =
  	equation_systems.add_system&lt;NonlinearImplicitSystem&gt; (<FONT COLOR="#BC8F8F"><B>&quot;Laplace-Young&quot;</FONT></B>);
  
  
        equation_systems.parameters.set&lt;Real&gt;         (<FONT COLOR="#BC8F8F"><B>&quot;nonlinear solver tolerance&quot;</FONT></B>)          = 1.e-12;
        equation_systems.parameters.set&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt; (<FONT COLOR="#BC8F8F"><B>&quot;nonlinear solver maximum iterations&quot;</FONT></B>) = 50;
  
        
        system.add_variable(<FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>,
  			  Utility::string_to_enum&lt;Order&gt;   (order),
  			  Utility::string_to_enum&lt;FEFamily&gt;(family));
  
        
        system.nonlinear_solver-&gt;residual = compute_residual;
        system.nonlinear_solver-&gt;jacobian = compute_jacobian;
  
        equation_systems.init();
  
        equation_systems.print_info();
      }
      
      equation_systems.get_system(<FONT COLOR="#BC8F8F"><B>&quot;Laplace-Young&quot;</FONT></B>).solve();
  
      GMVIO (mesh).write_equation_systems (<FONT COLOR="#BC8F8F"><B>&quot;out.gmv&quot;</FONT></B>, 
  					 equation_systems);
    }
    
    <B><FONT COLOR="#A020F0">return</FONT></B> 0; 
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example  ./ex19-dbg
***************************************************************
 
Running ./ex19-dbg -r 3 -o FIRST

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=225
    n_local_nodes()=225
  n_elem()=255
    n_local_elem()=255
    n_active_elem()=192
  n_subdomains()=1
  n_processors()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Laplace-Young"
    Type "NonlinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=225
    n_local_dofs()=225
    n_constrained_dofs()=0
    n_vectors()=1

  NL step 0, |residual|_2 = 2.000000e-01
  NL step 1, |residual|_2 = 4.432961e-03
  NL step 2, |residual|_2 = 2.163781e-04
  NL step 3, |residual|_2 = 1.157690e-05
  NL step 4, |residual|_2 = 6.567452e-07
  NL step 5, |residual|_2 = 3.849499e-08
  NL step 6, |residual|_2 = 2.293601e-09

 ---------------------------------------------------------------------------- 
| Reference count information                                                |
 ---------------------------------------------------------------------------- 
| 12SparseMatrixIdE reference count information:
|  Creations:    13
|  Destructions: 13
| 13NumericVectorIdE reference count information:
|  Creations:    23
|  Destructions: 23
| 15NonlinearSolverIdE reference count information:
|  Creations:    1
|  Destructions: 1
| 4Elem reference count information:
|  Creations:    1615
|  Destructions: 1615
| 4Node reference count information:
|  Creations:    225
|  Destructions: 225
| 5QBase reference count information:
|  Creations:    33
|  Destructions: 33
| 6DofMap reference count information:
|  Creations:    1
|  Destructions: 1
| 6FEBase reference count information:
|  Creations:    21
|  Destructions: 21
| 6System reference count information:
|  Creations:    1
|  Destructions: 1
| 9DofObject reference count information:
|  Creations:    1840
|  Destructions: 1840
| N10Parameters5ValueE reference count information:
|  Creations:    10
|  Destructions: 10
 ---------------------------------------------------------------------------- 

----------------------------------------------------------------------------------------------------------------
| Time:           Fri Oct 17 18:21:04 2008                                                                      |
| OS:             Darwin                                                                                        |
| HostName:       benjamin-kirks-macbook.local                                                                  |
| OS Release:     9.5.0                                                                                         |
| OS Version:     Darwin Kernel Version 9.5.0: Wed Sep  3 11:29:43 PDT 2008; root:xnu-1228.7.58~1/RELEASE_I386  |
| Machine:        i386                                                                                          |
| Username:       benkirk                                                                                       |
| Configuration:  ./configure run on Fri Oct 17 17:29:52 CDT 2008                                               |
----------------------------------------------------------------------------------------------------------------
 -------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.633672, Active time=0.503336                |
 -------------------------------------------------------------------------------
| Event                          nCalls    Total       Avg         Percent of   |
|                                          Time        Time        Active Time  |
|-------------------------------------------------------------------------------|
|                                                                               |
|                                                                               |
| DofMap                                                                        |
|   add_neighbors_to_send_list() 1         0.0021      0.002114    0.42         |
|   compute_sparsity()           1         0.0643      0.064275    12.77        |
|   create_dof_constraints()     1         0.0018      0.001836    0.36         |
|   distribute_dofs()            1         0.0030      0.003020    0.60         |
|   dof_indices()                2880      0.0166      0.000006    3.29         |
|   reinit()                     1         0.0035      0.003548    0.70         |
|   sort_send_list()             1         0.0015      0.001458    0.29         |
|                                                                               |
| FE                                                                            |
|   compute_affine_map()         2944      0.0562      0.000019    11.17        |
|   compute_face_map()           448       0.0146      0.000033    2.90         |
|   compute_shape_functions()    2944      0.0653      0.000022    12.97        |
|   init_face_shape_functions()  161       0.0016      0.000010    0.31         |
|   init_shape_functions()       461       0.0159      0.000034    3.16         |
|   inverse_map()                2688      0.0399      0.000015    7.93         |
|                                                                               |
| GMVIO                                                                         |
|   write_nodal_data()           1         0.0032      0.003171    0.63         |
|                                                                               |
| LocationMap                                                                   |
|   find()                       756       0.0032      0.000004    0.64         |
|   init()                       3         0.0010      0.000332    0.20         |
|                                                                               |
| Mesh                                                                          |
|   find_neighbors()             2         0.0251      0.012527    4.98         |
|   renumber_nodes_and_elem()    2         0.0014      0.000716    0.28         |
|                                                                               |
| MeshRefinement                                                                |
|   _refine_elements()           3         0.0077      0.002561    1.53         |
|   add_point()                  756       0.0046      0.000006    0.91         |
|                                                                               |
| Parallel                                                                      |
|   allgather()                  1         0.0000      0.000008    0.00         |
|                                                                               |
| Partitioner                                                                   |
|   single_partition()           2         0.0016      0.000804    0.32         |
|                                                                               |
| System                                                                        |
|   assemble()                   1         0.0000      0.000003    0.00         |
|   solve()                      1         0.1693      0.169319    33.64        |
 -------------------------------------------------------------------------------
| Totals:                        14060     0.5033                  100.00       |
 -------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./ex19-dbg
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
