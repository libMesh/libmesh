<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("miscellaneous_ex3",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file miscellaneous_ex3.C with comments: </h1> 
<div class = "comment">
<h1>Miscellaneous Example 3 - 2D Laplace-Young Problem Using Nonlinear Solvers</h1>

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
        #include "libmesh/libmesh.h"
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_refinement.h"
        #include "libmesh/exodusII_io.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/fe.h"
        #include "libmesh/quadrature_gauss.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/dense_vector.h"
        #include "libmesh/elem.h"
        #include "libmesh/string_to_enum.h"
        #include "libmesh/getpot.h"
        
</pre>
</div>
<div class = "comment">
The nonlinear solver and system we will be using
</div>

<div class ="fragment">
<pre>
        #include "libmesh/nonlinear_solver.h"
        #include "libmesh/nonlinear_implicit_system.h"
        
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
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
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
                               SparseMatrix&lt;Number&gt;&  jacobian,
                               NonlinearImplicitSystem& sys)
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
          std::vector&lt;dof_id_type&gt; dof_indices;
        
</pre>
</div>
<div class = "comment">
Now we will loop over all the active elements in the mesh which
are local to this processor.
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
                  Gradient grad_u;
            
                  for (unsigned int i=0; i&lt;phi.size(); i++)
                    grad_u += dphi[i][qp]*soln(dof_indices[i]);
                  
                  const Number K = 1./std::sqrt(1. + grad_u*grad_u);
                  
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
                               NumericVector&lt;Number&gt;& residual,
                               NonlinearImplicitSystem& sys)
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
          libmesh_assert_equal_to (dim, 2);
        
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
          std::vector&lt;dof_id_type&gt; dof_indices;
        
</pre>
</div>
<div class = "comment">
Now we will loop over all the active elements in the mesh which
are local to this processor.
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
                  Number u = 0;
                  Gradient grad_u;
                  
                  for (unsigned int j=0; j&lt;phi.size(); j++)
                    {
                      u      += phi[j][qp]*soln(dof_indices[j]);
                      grad_u += dphi[j][qp]*soln(dof_indices[j]);
                    }
                  
                  const Number K = 1./std::sqrt(1. + grad_u*grad_u);
                  
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
        
        #if !defined(LIBMESH_HAVE_PETSC) && !defined(LIBMESH_HAVE_TRILINOS)
          if (libMesh::processor_id() == 0)
            std::cerr &lt;&lt; "ERROR: This example requires libMesh to be\n"
                      &lt;&lt; "compiled with nonlinear solver support from\n"
                      &lt;&lt; "PETSc or Trilinos!"
                      &lt;&lt; std::endl;
          return 0;
        #endif
        
        #ifndef LIBMESH_ENABLE_AMR
          if (libMesh::processor_id() == 0)
            std::cerr &lt;&lt; "ERROR: This example requires libMesh to be\n"
                      &lt;&lt; "compiled with AMR support!"
                      &lt;&lt; std::endl;
          return 0;
        #else
        
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
and then abort.
</div>

<div class ="fragment">
<pre>
              libmesh_error();
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
              libmesh_error();
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
Skip this 2D example if libMesh was compiled as 1D-only.
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(2 &lt;= LIBMESH_DIM, "2D support");
          
</pre>
</div>
<div class = "comment">
Create a mesh from file.
</div>

<div class ="fragment">
<pre>
          Mesh mesh;    
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
  

<br><br>Creates a system named "Laplace-Young"
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
</div>

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
          
</pre>
</div>
<div class = "comment">
Solve the system "Laplace-Young", print the number of iterations
and final residual
</div>

<div class ="fragment">
<pre>
          equation_systems.get_system("Laplace-Young").solve();
        
</pre>
</div>
<div class = "comment">
Print out final convergence information.  This duplicates some
output from during the solve itself, but demonstrates another way
to get this information after the solve is complete.
</div>

<div class ="fragment">
<pre>
          std::cout &lt;&lt; "Laplace-Young system solved at nonlinear iteration "
        	    &lt;&lt; system.n_nonlinear_iterations()
        	    &lt;&lt; " , final nonlinear residual norm: "
        	    &lt;&lt; system.final_nonlinear_residual()
        	    &lt;&lt; std::endl;
        
        #ifdef LIBMESH_HAVE_EXODUS_API
</pre>
</div>
<div class = "comment">
After solving the system write the solution
</div>

<div class ="fragment">
<pre>
          ExodusII_IO (mesh).write_equation_systems ("out.e", 
                                               equation_systems);
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
        #endif // #ifndef LIBMESH_ENABLE_AMR
        
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
<br><br><br> <h1> The source file miscellaneous_ex3.C without comments: </h1> 
<pre> 
   
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;cmath&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/string_to_enum.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/nonlinear_solver.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/nonlinear_implicit_system.h&quot;</FONT></B>
  
  #ifdef LIBMESH_HAVE_PETSC
  #include &lt;petsc.h&gt;
  #endif
  
  using namespace libMesh;
  
  EquationSystems *_equation_system = NULL;
  
  <B><FONT COLOR="#228B22">const</FONT></B> Real kappa = 1.;
  <B><FONT COLOR="#228B22">const</FONT></B> Real sigma = 0.2;
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> compute_jacobian (<B><FONT COLOR="#228B22">const</FONT></B> NumericVector&lt;Number&gt;&amp; soln,
                         SparseMatrix&lt;Number&gt;&amp;  jacobian,
                         NonlinearImplicitSystem&amp; sys)
  {
    EquationSystems &amp;es = *_equation_system;
  
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    NonlinearImplicitSystem&amp; system = 
      es.get_system&lt;NonlinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Laplace-Young&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
  
    FEType fe_type = dof_map.variable_type(0);
  
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    QGauss qrule (dim, FIFTH);
  
    fe-&gt;attach_quadrature_rule (&amp;qrule);
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
    
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    DenseMatrix&lt;Number&gt; Ke;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
  
        fe-&gt;reinit (elem);
  
        Ke.resize (dof_indices.size(),
                   dof_indices.size());
             
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
          {
            Gradient grad_u;
      
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
              grad_u += dphi[i][qp]*soln(dof_indices[i]);
            
            <B><FONT COLOR="#228B22">const</FONT></B> Number K = 1./std::sqrt(1. + grad_u*grad_u);
            
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
                Ke(i,j) += JxW[qp]*(
                                    K*(dphi[i][qp]*dphi[j][qp]) +
                                    kappa*phi[i][qp]*phi[j][qp]
                                    );
          }
        
        dof_map.constrain_element_matrix (Ke, dof_indices);
        
        jacobian.add_matrix (Ke, dof_indices);
      }
  
  }
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> compute_residual (<B><FONT COLOR="#228B22">const</FONT></B> NumericVector&lt;Number&gt;&amp; soln,
                         NumericVector&lt;Number&gt;&amp; residual,
                         NonlinearImplicitSystem&amp; sys)
  {
    EquationSystems &amp;es = *_equation_system;
  
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
    libmesh_assert_equal_to (dim, 2);
  
    NonlinearImplicitSystem&amp; system = 
      es.get_system&lt;NonlinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Laplace-Young&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
  
    FEType fe_type = dof_map.variable_type(0);
  
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    QGauss qrule (dim, FIFTH);
  
    fe-&gt;attach_quadrature_rule (&amp;qrule);
  
    AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
                
    QGauss qface(dim-1, FIFTH);
    
    fe_face-&gt;attach_quadrature_rule (&amp;qface);
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
    
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    DenseVector&lt;Number&gt; Re;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
  
    residual.zero();
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
  
        fe-&gt;reinit (elem);
  
        Re.resize (dof_indices.size());
        
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
          {
            Number u = 0;
            Gradient grad_u;
            
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
              {
                u      += phi[j][qp]*soln(dof_indices[j]);
                grad_u += dphi[j][qp]*soln(dof_indices[j]);
              }
            
            <B><FONT COLOR="#228B22">const</FONT></B> Number K = 1./std::sqrt(1. + grad_u*grad_u);
            
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
              Re(i) += JxW[qp]*(
                                K*(dphi[i][qp]*grad_u) +
                                kappa*phi[i][qp]*u
                                );
          }
  
        
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> side=0; side&lt;elem-&gt;n_sides(); side++)
          <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(side) == NULL)
            {
              <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp;  phi_face = fe_face-&gt;get_phi();
  
              <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
  
              fe_face-&gt;reinit(elem, side);
  
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qface.n_points(); qp++)
                {
                  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi_face.size(); i++)
                    Re(i) -= JxW_face[qp]*sigma*phi_face[i][qp];
                } 
            }
        
        dof_map.constrain_element_vector (Re, dof_indices);
        residual.add_vector (Re, dof_indices);
      }
  
  }
  
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_HAVE_PETSC) &amp;&amp; !defined(LIBMESH_HAVE_TRILINOS)
    <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR: This example requires libMesh to be\n&quot;</FONT></B>
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;compiled with nonlinear solver support from\n&quot;</FONT></B>
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;PETSc or Trilinos!&quot;</FONT></B>
                &lt;&lt; std::endl;
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  #endif
  
  #ifndef LIBMESH_ENABLE_AMR
    <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR: This example requires libMesh to be\n&quot;</FONT></B>
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;compiled with AMR support!&quot;</FONT></B>
                &lt;&lt; std::endl;
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    GetPot command_line (argc, argv);
    
    <B><FONT COLOR="#A020F0">if</FONT></B> (argc &lt; 3)
      {
        <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Usage:\n&quot;</FONT></B>
                    &lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;\t &quot;</FONT></B> &lt;&lt; argv[0] &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; -r 2&quot;</FONT></B>
                    &lt;&lt; std::endl;
  
        libmesh_error();
      }
    
    <B><FONT COLOR="#A020F0">else</FONT></B> 
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Running &quot;</FONT></B> &lt;&lt; argv[0];
        
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">int</FONT></B> i=1; i&lt;argc; i++)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B> &lt;&lt; argv[i];
        
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl &lt;&lt; std::endl;
      }
    
  
    <B><FONT COLOR="#228B22">int</FONT></B> nr = 2;
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <B><FONT COLOR="#BC8F8F">&quot;-r&quot;</FONT></B>) )
      nr = command_line.next(nr);
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string order = <B><FONT COLOR="#BC8F8F">&quot;FIRST&quot;</FONT></B>; 
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(2, <B><FONT COLOR="#BC8F8F">&quot;-Order&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;-o&quot;</FONT></B>) )
      order = command_line.next(order);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string family = <B><FONT COLOR="#BC8F8F">&quot;LAGRANGE&quot;</FONT></B>; 
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(2, <B><FONT COLOR="#BC8F8F">&quot;-FEFamily&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;-f&quot;</FONT></B>) )
      family = command_line.next(family);
    
    <B><FONT COLOR="#A020F0">if</FONT></B> ((family == <B><FONT COLOR="#BC8F8F">&quot;MONOMIAL&quot;</FONT></B>) || (family == <B><FONT COLOR="#BC8F8F">&quot;XYZ&quot;</FONT></B>))
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ex19 currently requires a C^0 (or higher) FE basis.&quot;</FONT></B> &lt;&lt; std::endl;
        libmesh_error();
      }
  
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <B><FONT COLOR="#BC8F8F">&quot;-pre&quot;</FONT></B>) )
      {
  #ifdef LIBMESH_HAVE_PETSC
        PetscOptionsSetValue(<B><FONT COLOR="#BC8F8F">&quot;-snes_mf_operator&quot;</FONT></B>,PETSC_NULL);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;Must be using PetsC to use jacobian based preconditioning&quot;</FONT></B>&lt;&lt;std::endl;
  
        <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  #endif <I><FONT COLOR="#B22222">//LIBMESH_HAVE_PETSC
</FONT></I>      }  
      
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
    
    Mesh mesh;    
    mesh.read (<B><FONT COLOR="#BC8F8F">&quot;lshaped.xda&quot;</FONT></B>);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (order != <B><FONT COLOR="#BC8F8F">&quot;FIRST&quot;</FONT></B>)
      mesh.all_second_order();
  
    MeshRefinement(mesh).uniformly_refine(nr);
  
    mesh.print_info();    
    
    EquationSystems equation_systems (mesh);
    _equation_system = &amp;equation_systems;
    
    
    NonlinearImplicitSystem&amp; system =
      equation_systems.add_system&lt;NonlinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Laplace-Young&quot;</FONT></B>);
  
  
    equation_systems.parameters.set&lt;Real&gt;         (<B><FONT COLOR="#BC8F8F">&quot;nonlinear solver tolerance&quot;</FONT></B>)          = 1.e-12;
    equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; (<B><FONT COLOR="#BC8F8F">&quot;nonlinear solver maximum iterations&quot;</FONT></B>) = 50;
  
      
    system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>,
  		      <B><FONT COLOR="#5F9EA0">Utility</FONT></B>::string_to_enum&lt;Order&gt;   (order),
  		      <B><FONT COLOR="#5F9EA0">Utility</FONT></B>::string_to_enum&lt;FEFamily&gt;(family));
  
    system.nonlinear_solver-&gt;residual = compute_residual;
    system.nonlinear_solver-&gt;jacobian = compute_jacobian;
  
    equation_systems.init();
  
    equation_systems.print_info();
    
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Laplace-Young&quot;</FONT></B>).solve();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Laplace-Young system solved at nonlinear iteration &quot;</FONT></B>
  	    &lt;&lt; system.n_nonlinear_iterations()
  	    &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; , final nonlinear residual norm: &quot;</FONT></B>
  	    &lt;&lt; system.final_nonlinear_residual()
  	    &lt;&lt; std::endl;
  
  #ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO (mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;out.e&quot;</FONT></B>, 
                                         equation_systems);
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>  #endif <I><FONT COLOR="#B22222">// #ifndef LIBMESH_ENABLE_AMR
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0; 
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example miscellaneous_ex3:
*  mpirun -np 12 example-devel -r 3 -o FIRST -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/miscellaneous/miscellaneous_ex3/.libs/lt-example-devel -r 3 -o FIRST -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=225
    n_local_nodes()=26
  n_elem()=255
    n_local_elem()=24
    n_active_elem()=192
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Laplace-Young"
    Type "NonlinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=225
    n_local_dofs()=26
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.12889
      Average Off-Processor Bandwidth <= 2.09778
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 9
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

  NL step  0, |residual|_2 = 2.000000e-01
  NL step  1, |residual|_2 = 4.432603e-03
  NL step  2, |residual|_2 = 2.163693e-04
  NL step  3, |residual|_2 = 1.157686e-05
  NL step  4, |residual|_2 = 6.567530e-07
  NL step  5, |residual|_2 = 3.849507e-08
  NL step  6, |residual|_2 = 2.293588e-09
Laplace-Young system solved at nonlinear iteration 6 , final nonlinear residual norm: 2.293588e-09
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/miscellaneous/miscellaneous_ex3/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:08:58 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           2.877e-01      1.00164   2.876e-01
Objects:              7.600e+01      1.00000   7.600e+01
Flops:                4.194e+05      3.34309   2.801e+05  3.361e+06
Flops/sec:            1.460e+06      3.34663   9.740e+05  1.169e+07
MPI Messages:         9.980e+02      2.52020   6.627e+02  7.952e+03
MPI Message Lengths:  4.084e+04      2.07414   4.612e+01  3.667e+05
MPI Reductions:       5.650e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.8755e-01 100.0%  3.3611e+06 100.0%  7.952e+03 100.0%  4.612e+01      100.0%  5.640e+02  99.8% 

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

VecDot                 6 1.0 6.1989e-05 1.4 3.18e+02 3.1 0.0e+00 0.0e+00 6.0e+00  0  0  0  0  1   0  0  0  0  1    42
VecMDot              158 1.0 1.5593e-03 1.4 1.15e+05 3.1 0.0e+00 0.0e+00 1.6e+02  0 28  0  0 28   0 28  0  0 28   608
VecNorm              176 1.0 1.9641e-03 2.0 9.50e+03 3.0 0.0e+00 0.0e+00 1.8e+02  1  2  0  0 31   1  2  0  0 31    40
VecScale             164 1.0 9.4414e-05 1.3 4.43e+03 3.0 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0   391
VecCopy               32 1.0 3.1233e-05 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               194 1.0 7.7963e-05 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                6 1.0 3.8862e-05 1.8 3.24e+02 3.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0    69
VecWAXPY               6 1.0 2.4796e-05 2.3 1.62e+02 3.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0    54
VecMAXPY             164 1.0 1.3041e-04 1.8 1.25e+05 3.0 0.0e+00 0.0e+00 0.0e+00  0 31  0  0  0   0 31  0  0  0  8009
VecAssemblyBegin      54 1.0 9.7177e-03 5.0 0.00e+00 0.0 2.8e+02 5.2e+01 1.6e+02  2  0  4  4 29   2  0  4  4 29     0
VecAssemblyEnd        54 1.0 7.6294e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      178 1.0 5.6839e-04 1.5 0.00e+00 0.0 7.1e+03 4.0e+01 0.0e+00  0  0 90 79  0   0  0 90 79  0     0
VecScatterEnd        178 1.0 3.1528e-0311.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
VecReduceArith        13 1.0 1.5273e-02 2.1 6.89e+02 3.1 0.0e+00 0.0e+00 0.0e+00  3  0  0  0  0   3  0  0  0  0     0
VecReduceComm          7 1.0 8.1971e-0381.7 0.00e+00 0.0 0.0e+00 0.0e+00 7.0e+00  2  0  0  0  1   2  0  0  0  1     0
VecNormalize         164 1.0 2.0597e-03 1.9 1.33e+04 3.0 0.0e+00 0.0e+00 1.6e+02  1  3  0  0 29   1  3  0  0 29    54
MatMult              164 1.0 1.0922e-03 1.5 6.48e+04 2.7 6.6e+03 3.9e+01 0.0e+00  0 17 82 69  0   0 17 82 69  0   514
MatSolve             164 1.0 2.4509e-04 3.7 8.97e+04 5.4 0.0e+00 0.0e+00 0.0e+00  0 18  0  0  0   0 18  0  0  0  2533
MatLUFactorNum         6 1.0 8.5831e-05 2.5 9.56e+03 8.6 0.0e+00 0.0e+00 0.0e+00  0  2  0  0  0   0  2  0  0  0   680
MatILUFactorSym        1 1.0 1.0085e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  0  0  0  0  1   0  0  0  0  1     0
MatAssemblyBegin      12 1.0 4.8513e-03 9.4 0.00e+00 0.0 3.6e+02 1.7e+02 2.4e+01  1  0  5 17  4   1  0  5 17  4     0
MatAssemblyEnd        12 1.0 1.1139e-03 1.8 0.00e+00 0.0 8.0e+01 1.2e+01 8.0e+00  0  0  1  0  1   0  0  1  0  1     0
MatGetRowIJ            1 1.0 2.1458e-06 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 7.2002e-05 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         8 1.0 1.7166e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
SNESSolve              1 1.0 7.5193e-02 1.0 4.19e+05 3.3 7.8e+03 4.6e+01 5.4e+02 26100 98 99 96  26100 98 99 97    45
SNESFunctionEval       7 1.0 2.8230e-02 1.0 0.00e+00 0.0 5.7e+02 5.5e+01 1.0e+02 10  0  7  9 19  10  0  7  9 19     0
SNESJacobianEval       6 1.0 2.3197e-02 1.0 0.00e+00 0.0 6.9e+02 1.1e+02 8.6e+01  8  0  9 21 15   8  0  9 21 15     0
SNESLineSearch         6 1.0 2.3658e-02 1.0 4.13e+03 2.8 7.3e+02 5.0e+01 1.1e+02  8  1  9 10 20   8  1  9 10 20     1
KSPGMRESOrthog       158 1.0 1.7982e-03 1.3 2.31e+05 3.1 0.0e+00 0.0e+00 1.6e+02  1 57  0  0 28   1 57  0  0 28  1068
KSPSetUp              12 1.0 1.1396e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               6 1.0 7.4975e-03 1.0 4.15e+05 3.3 6.3e+03 3.9e+01 3.3e+02  3 99 79 67 58   3 99 79 67 58   444
PCSetUp               12 1.0 8.1897e-04 1.3 9.56e+03 8.6 0.0e+00 0.0e+00 7.0e+00  0  2  0  0  1   0  2  0  0  1    71
PCSetUpOnBlocks        6 1.0 4.6396e-04 1.3 9.56e+03 8.6 0.0e+00 0.0e+00 5.0e+00  0  2  0  0  1   0  2  0  0  1   126
PCApply              164 1.0 1.7033e-03 1.2 8.97e+04 5.4 0.0e+00 0.0e+00 0.0e+00  1 18  0  0  0   1 18  0  0  0   364
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

           Container     2              2         1096     0
              Vector    47             47        79000     0
      Vector Scatter     2              2         2072     0
           Index Set     7              7         5400     0
   IS L to G Mapping     1              1          564     0
              Matrix     4              4        16992     0
                SNES     1              1         1268     0
      SNESLineSearch     1              1          840     0
    Distributed Mesh     2              2         8552     0
     Bipartite Graph     4              4         2736     0
       Krylov Solver     2              2        19360     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 5.00679e-06
Average time for zero size MPI_Send(): 1.32521e-05
#PETSc Option Table entries:
-ksp_right_pc
-log_summary
-o FIRST
-pc_type bjacobi
-r 3
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
| Time:           Thu Jan 31 22:08:58 2013                                                                             |
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
| libMesh Performance: Alive time=0.560824, Active time=0.271314                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0011      0.001141    0.0023      0.002310    0.42     0.85     |
|   build_sparsity()                 1         0.0014      0.001389    0.0037      0.003704    0.51     1.37     |
|   create_dof_constraints()         1         0.0009      0.000877    0.0009      0.000877    0.32     0.32     |
|   distribute_dofs()                1         0.0043      0.004348    0.0135      0.013462    1.60     4.96     |
|   dof_indices()                    261       0.0145      0.000056    0.0145      0.000056    5.36     5.36     |
|   prepare_send_list()              1         0.0000      0.000035    0.0000      0.000035    0.01     0.01     |
|   reinit()                         1         0.0078      0.007755    0.0078      0.007755    2.86     2.86     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0005      0.000492    0.0018      0.001810    0.18     0.67     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0030      0.002951    0.0030      0.002951    1.09     1.09     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        264       0.0055      0.000021    0.0055      0.000021    2.02     2.02     |
|   init_shape_functions()           69        0.0012      0.000017    0.0012      0.000017    0.44     0.44     |
|   inverse_map()                    168       0.0016      0.000010    0.0016      0.000010    0.60     0.60     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             264       0.0034      0.000013    0.0034      0.000013    1.25     1.25     |
|   compute_face_map()               56        0.0014      0.000025    0.0031      0.000055    0.52     1.14     |
|   init_face_shape_functions()      7         0.0001      0.000015    0.0001      0.000015    0.04     0.04     |
|   init_reference_to_physical_map() 69        0.0019      0.000028    0.0019      0.000028    0.71     0.71     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           756       0.0032      0.000004    0.0032      0.000004    1.17     1.17     |
|   init()                           3         0.0003      0.000096    0.0003      0.000096    0.11     0.11     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 2         0.0062      0.003095    0.0077      0.003872    2.28     2.85     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   broadcast()                      1         0.0009      0.000905    0.0032      0.003229    0.33     1.19     |
|   compute_hilbert_indices()        3         0.0017      0.000581    0.0017      0.000581    0.64     0.64     |
|   find_global_indices()            3         0.0013      0.000441    0.0081      0.002712    0.49     3.00     |
|   parallel_sort()                  3         0.0030      0.000987    0.0033      0.001098    1.09     1.21     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0001      0.000139    0.0050      0.005023    0.05     1.85     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _refine_elements()               3         0.0054      0.001783    0.0123      0.004086    1.97     4.52     |
|   add_point()                      756       0.0033      0.000004    0.0066      0.000009    1.20     2.42     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      2         0.0274      0.013690    0.0330      0.016514    10.09    12.17    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      11        0.0004      0.000038    0.0005      0.000046    0.15     0.19     |
|   broadcast()                      9         0.0006      0.000065    0.0005      0.000058    0.22     0.19     |
|   max(bool)                        4         0.0002      0.000047    0.0002      0.000047    0.07     0.07     |
|   max(scalar)                      261       0.0041      0.000016    0.0041      0.000016    1.50     1.50     |
|   max(vector)                      62        0.0013      0.000021    0.0040      0.000065    0.49     1.48     |
|   min(bool)                        312       0.0044      0.000014    0.0044      0.000014    1.63     1.63     |
|   min(scalar)                      253       0.0068      0.000027    0.0068      0.000027    2.49     2.49     |
|   min(vector)                      62        0.0014      0.000022    0.0041      0.000067    0.51     1.52     |
|   probe()                          176       0.0010      0.000006    0.0010      0.000006    0.38     0.38     |
|   receive()                        176       0.0010      0.000006    0.0021      0.000012    0.36     0.76     |
|   send()                           176       0.0005      0.000003    0.0005      0.000003    0.19     0.19     |
|   send_receive()                   182       0.0014      0.000008    0.0044      0.000024    0.51     1.62     |
|   sum()                            36        0.0005      0.000015    0.0007      0.000021    0.20     0.27     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           176       0.0003      0.000002    0.0003      0.000002    0.12     0.12     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         2         0.0011      0.000556    0.0024      0.001182    0.41     0.87     |
|   set_parent_processor_ids()       2         0.0006      0.000316    0.0006      0.000316    0.23     0.23     |
|                                                                                                                |
| PetscNonlinearSolver                                                                                           |
|   jacobian()                       6         0.0105      0.001746    0.0228      0.003793    3.86     8.39     |
|   residual()                       7         0.0092      0.001308    0.0281      0.004017    3.38     10.36    |
|   solve()                          1         0.1246      0.124573    0.1755      0.175452    45.91    64.67    |
|                                                                                                                |
| System                                                                                                         |
|   solve()                          1         0.0000      0.000050    0.1755      0.175503    0.02     64.69    |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            4614      0.2713                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example miscellaneous_ex3:
*  mpirun -np 12 example-devel -r 3 -o FIRST -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
