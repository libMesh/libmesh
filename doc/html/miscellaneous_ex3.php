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
                               NonlinearImplicitSystem& /*sys*/)
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
                               NonlinearImplicitSystem& /*sys*/)
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
Create a mesh, with dimension to be overridden by the file,
distributed across the default MPI communicator.
</div>

<div class ="fragment">
<pre>
          Mesh mesh(init.comm());
        
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
                         NonlinearImplicitSystem&amp; <I><FONT COLOR="#B22222">/*sys*/</FONT></I>)
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
                         NonlinearImplicitSystem&amp; <I><FONT COLOR="#B22222">/*sys*/</FONT></I>)
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
  
    Mesh mesh(init.comm());
  
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
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/miscellaneous/miscellaneous_ex3'
***************************************************************
* Running Example miscellaneous_ex3:
*  mpirun -np 4 example-devel -r 3 -o FIRST -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
Running /net/spark/workspace/roystgnr/libmesh/git/devel/examples/miscellaneous/miscellaneous_ex3/.libs/lt-example-devel -r 3 -o FIRST -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=225
    n_local_nodes()=66
  n_elem()=255
    n_local_elem()=70
    n_active_elem()=192
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
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
    n_local_dofs()=66
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.68444
      Average Off-Processor Bandwidth <= 0.853333
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 6
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

  NL step  0, |residual|_2 = 2.000000e-01
  NL step  1, |residual|_2 = 4.433526e-03
  NL step  2, |residual|_2 = 2.163914e-04
  NL step  3, |residual|_2 = 1.157687e-05
  NL step  4, |residual|_2 = 6.567375e-07
  NL step  5, |residual|_2 = 3.849456e-08
  NL step  6, |residual|_2 = 2.293554e-09
Laplace-Young system solved at nonlinear iteration 6 , final nonlinear residual norm: 2.293554e-09

 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:50:53 2013                                                                          |
| OS:             Linux                                                                                             |
| HostName:       spark.ices.utexas.edu                                                                             |
| OS Release:     2.6.32-279.22.1.el6.x86_64                                                                        |
| OS Version:     #1 SMP Tue Feb 5 14:33:39 CST 2013                                                                |
| Machine:        x86_64                                                                                            |
| Username:       roystgnr                                                                                          |
| Configuration:  ../configure  '--enable-everything'                                                               |
|  'METHODS=devel'                                                                                                  |
|  '--prefix=/h2/roystgnr/libmesh-test'                                                                             |
|  'CXX=distcc /usr/bin/g++'                                                                                        |
|  'CC=distcc /usr/bin/gcc'                                                                                         |
|  'FC=distcc /usr/bin/gfortran'                                                                                    |
|  'F77=distcc /usr/bin/gfortran'                                                                                   |
|  'PETSC_DIR=/opt/apps/ossw/libraries/petsc/petsc-3.3-p2'                                                          |
|  'PETSC_ARCH=gcc-system-mkl-gf-10.3.12.361-mpich2-1.4.1p1-cxx-opt'                                                |
|  'SLEPC_DIR=/opt/apps/ossw/libraries/slepc/slepc-3.3-p2-petsc-3.3-p2-cxx-opt'                                     |
|  'TRILINOS_DIR=/opt/apps/ossw/libraries/trilinos/trilinos-10.12.2/sl6/gcc-system/mpich2-1.4.1p1/mkl-gf-10.3.12.361'|
|  'VTK_DIR=/opt/apps/ossw/libraries/vtk/vtk-5.10.0/sl6/gcc-system'                                                 |
|  'HDF5_DIR=/opt/apps/ossw/libraries/hdf5/hdf5-1.8.9/sl6/gcc-system'                                               |
 -------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.434914, Active time=0.237467                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0003      0.000330    0.0006      0.000557    0.14     0.23     |
|   build_sparsity()                 1         0.0004      0.000359    0.0009      0.000946    0.15     0.40     |
|   create_dof_constraints()         1         0.0002      0.000229    0.0002      0.000229    0.10     0.10     |
|   distribute_dofs()                1         0.0007      0.000735    0.0020      0.002036    0.31     0.86     |
|   dof_indices()                    761       0.0044      0.000006    0.0044      0.000006    1.86     1.86     |
|   prepare_send_list()              1         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   reinit()                         1         0.0010      0.000952    0.0010      0.000952    0.40     0.40     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0002      0.000236    0.0008      0.000777    0.10     0.33     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.1547      0.154657    0.1547      0.154657    65.13    65.13    |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        722       0.0036      0.000005    0.0036      0.000005    1.51     1.51     |
|   init_shape_functions()           111       0.0003      0.000003    0.0003      0.000003    0.14     0.14     |
|   inverse_map()                    294       0.0011      0.000004    0.0011      0.000004    0.45     0.45     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             722       0.0025      0.000003    0.0025      0.000003    1.05     1.05     |
|   compute_face_map()               98        0.0008      0.000008    0.0019      0.000019    0.34     0.80     |
|   init_face_shape_functions()      7         0.0000      0.000004    0.0000      0.000004    0.01     0.01     |
|   init_reference_to_physical_map() 111       0.0008      0.000007    0.0008      0.000007    0.32     0.32     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           756       0.0010      0.000001    0.0010      0.000001    0.40     0.40     |
|   init()                           3         0.0001      0.000023    0.0001      0.000023    0.03     0.03     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 2         0.0011      0.000531    0.0012      0.000615    0.45     0.52     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   broadcast()                      1         0.0002      0.000241    0.0005      0.000548    0.10     0.23     |
|   compute_hilbert_indices()        3         0.0008      0.000277    0.0008      0.000277    0.35     0.35     |
|   find_global_indices()            3         0.0002      0.000083    0.0021      0.000702    0.10     0.89     |
|   parallel_sort()                  3         0.0004      0.000150    0.0006      0.000203    0.19     0.26     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0001      0.000068    0.1556      0.155581    0.03     65.52    |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _refine_elements()               3         0.0013      0.000445    0.0041      0.001358    0.56     1.72     |
|   add_point()                      756       0.0014      0.000002    0.0025      0.000003    0.61     1.04     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      2         0.0026      0.001279    0.0043      0.002138    1.08     1.80     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      11        0.0001      0.000013    0.0002      0.000018    0.06     0.08     |
|   broadcast()                      27        0.0001      0.000005    0.0001      0.000003    0.05     0.03     |
|   max(bool)                        4         0.0001      0.000035    0.0001      0.000035    0.06     0.06     |
|   max(scalar)                      261       0.0013      0.000005    0.0013      0.000005    0.56     0.56     |
|   max(vector)                      62        0.0005      0.000008    0.0015      0.000024    0.21     0.62     |
|   min(bool)                        312       0.0014      0.000004    0.0014      0.000004    0.57     0.57     |
|   min(scalar)                      253       0.0019      0.000007    0.0019      0.000007    0.79     0.79     |
|   min(vector)                      62        0.0021      0.000034    0.0030      0.000048    0.89     1.26     |
|   probe()                          48        0.0001      0.000003    0.0001      0.000003    0.06     0.06     |
|   receive()                        48        0.0002      0.000003    0.0003      0.000006    0.06     0.13     |
|   send()                           48        0.0001      0.000002    0.0001      0.000002    0.04     0.04     |
|   send_receive()                   54        0.0003      0.000005    0.0008      0.000014    0.11     0.32     |
|   sum()                            36        0.0003      0.000007    0.0004      0.000011    0.11     0.17     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           48        0.0001      0.000001    0.0001      0.000001    0.03     0.03     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         2         0.0003      0.000168    0.0007      0.000329    0.14     0.28     |
|   set_parent_processor_ids()       2         0.0001      0.000061    0.0001      0.000061    0.05     0.05     |
|                                                                                                                |
| PetscNonlinearSolver                                                                                           |
|   jacobian()                       6         0.0096      0.001603    0.0166      0.002770    4.05     7.00     |
|   residual()                       7         0.0096      0.001366    0.0191      0.002734    4.03     8.06     |
|   solve()                          1         0.0290      0.028998    0.0648      0.064763    12.21    27.27    |
|                                                                                                                |
| System                                                                                                         |
|   solve()                          1         0.0000      0.000024    0.0648      0.064787    0.01     27.28    |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            5660      0.2375                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example miscellaneous_ex3:
*  mpirun -np 4 example-devel -r 3 -o FIRST -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/miscellaneous/miscellaneous_ex3'
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
