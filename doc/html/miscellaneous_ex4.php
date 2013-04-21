<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("miscellaneous_ex4",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file miscellaneous_ex4.C with comments: </h1> 
<div class = "comment">
<h1>Miscellaneous Example 4 - Using a shell matrix</h1>

<br><br>This example solves the equation

<br><br>\f$-\Delta u+\int u = 1\f$

<br><br>with homogeneous Dirichlet boundary conditions.  This system has
a full system matrix which can be written as the sum of of sparse
matrix and a rank 1 matrix.  The shell matrix concept is used to
solve this problem.

<br><br>The problem is solved in parallel on a non-uniform grid in order
to demonstrate all the techniques that are required for this.
The grid is fixed, however, i.e. no adaptive mesh refinement is
used, so that the example remains simple.

<br><br>The example is 2d; extension to 3d is straight forward.


<br><br>C++ include files that we need
</div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
        #include &lt;algorithm&gt;
        #include &lt;cstdlib&gt; // *must* precede &lt;cmath&gt; for proper std:abs() on PGI, Sun Studio CC
        #include &lt;cmath&gt;
        
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
        #include "libmesh/vtk_io.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/fe.h"
        #include "libmesh/quadrature_gauss.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/dense_vector.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/sum_shell_matrix.h"
        #include "libmesh/tensor_shell_matrix.h"
        #include "libmesh/sparse_shell_matrix.h"
        #include "libmesh/mesh_refinement.h"
        
        #include "libmesh/getpot.h"
        
</pre>
</div>
<div class = "comment">
This example will solve a linear transient system,
so we need to include the \p TransientLinearImplicitSystem definition.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/transient_system.h"
        #include "libmesh/linear_implicit_system.h"
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
Function prototype.  This function will assemble the system matrix
and right-hand-side.
</div>

<div class ="fragment">
<pre>
        void assemble (EquationSystems& es,
        	       const std::string& system_name);
        
</pre>
</div>
<div class = "comment">
Begin the main program.  Note that the first
statement in the program throws an error if
you are in complex number mode, since this
example is only intended to work with real
numbers.
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
        
        #if !defined(LIBMESH_ENABLE_AMR)
          libmesh_example_assert(false, "--enable-amr");
        #else
          libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");
        
</pre>
</div>
<div class = "comment">
Brief message to the user regarding the program name
and command line arguments.


<br><br></div>

<div class ="fragment">
<pre>
          std::cout &lt;&lt; "Running: " &lt;&lt; argv[0];
        
          for (int i=1; i&lt;argc; i++)
            std::cout &lt;&lt; " " &lt;&lt; argv[i];
        
          std::cout &lt;&lt; std::endl &lt;&lt; std::endl;
        
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
Create a mesh, with dimension to be overridden later, distributed
across the default MPI communicator.
</div>

<div class ="fragment">
<pre>
          Mesh mesh(init.comm());
        
</pre>
</div>
<div class = "comment">
Create an equation systems object.
</div>

<div class ="fragment">
<pre>
          EquationSystems equation_systems (mesh);
        
          MeshTools::Generation::build_square (mesh,
        				       16,
        				       16,
        				       -1., 1.,
        				       -1., 1.,
        				       QUAD4);
        
          LinearImplicitSystem & system =
            equation_systems.add_system&lt;LinearImplicitSystem&gt;
            ("System");
        
</pre>
</div>
<div class = "comment">
Adds the variable "u" to "System".  "u"
will be approximated using first-order approximation.
</div>

<div class ="fragment">
<pre>
          system.add_variable ("u", FIRST);
        
</pre>
</div>
<div class = "comment">
Also, we need to add two vectors.  The tensor matrix v*w^T of
these two vectors will be part of the system matrix.
</div>

<div class ="fragment">
<pre>
          system.add_vector("v",false);
          system.add_vector("w",false);
        
</pre>
</div>
<div class = "comment">
We need an additional matrix to be used for preconditioning since
a shell matrix is not suitable for that.
</div>

<div class ="fragment">
<pre>
          system.add_matrix("Preconditioner");
        
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly function.
</div>

<div class ="fragment">
<pre>
          system.attach_assemble_function (assemble);
        
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
        
          equation_systems.parameters.set&lt;unsigned int&gt;
            ("linear solver maximum iterations") = 250;
          equation_systems.parameters.set&lt;Real&gt;
            ("linear solver tolerance") = TOLERANCE;
        
</pre>
</div>
<div class = "comment">
Refine arbitrarily some elements.
</div>

<div class ="fragment">
<pre>
          for(unsigned int i=0; i&lt;2; i++)
            {
              MeshRefinement mesh_refinement(mesh);
              MeshBase::element_iterator       elem_it  = mesh.elements_begin();
              const MeshBase::element_iterator elem_end = mesh.elements_end();
              for (; elem_it != elem_end; ++elem_it)
        	{
        	  Elem* elem = *elem_it;
        	  if(elem-&gt;active())
        	    {
        	      if((elem-&gt;id()%20)&gt;8)
        		{
        		  elem-&gt;set_refinement_flag(Elem::REFINE);
        		}
        	      else
        		{
        		  elem-&gt;set_refinement_flag(Elem::DO_NOTHING);
        		}
        	    }
        	  else
        	    {
        	      elem-&gt;set_refinement_flag(Elem::INACTIVE);
        	    }
        	}
              mesh_refinement.refine_elements();
              equation_systems.reinit();
            }
        
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
Before the assemblation of the matrix, we have to clear the two
vectors that form the tensor matrix (since this is not performed
automatically).
</div>

<div class ="fragment">
<pre>
          system.get_vector("v").init(system.n_dofs(), system.n_local_dofs());
          system.get_vector("w").init(system.n_dofs(), system.n_local_dofs());
        
</pre>
</div>
<div class = "comment">
We need a shell matrix to solve.  There is currently no way to
store the shell matrix in the system.  We just create it locally
here (a shell matrix does not occupy much memory).
</div>

<div class ="fragment">
<pre>
          SumShellMatrix&lt;Number&gt; shellMatrix(system.comm());
          TensorShellMatrix&lt;Number&gt; shellMatrix0(system.get_vector("v"),system.get_vector("w"));
          shellMatrix.matrices.push_back(&shellMatrix0);
          SparseShellMatrix&lt;Number&gt; shellMatrix1(*system.matrix);
          shellMatrix.matrices.push_back(&shellMatrix1);
        
</pre>
</div>
<div class = "comment">
Attach that to the system.
</div>

<div class ="fragment">
<pre>
          system.attach_shell_matrix(&shellMatrix);
        
</pre>
</div>
<div class = "comment">
Reset the preconditioning matrix to zero (for the system matrix,
the same thing is done automatically).
</div>

<div class ="fragment">
<pre>
          system.get_matrix("Preconditioner").zero();
        
</pre>
</div>
<div class = "comment">
Assemble & solve the linear system
</div>

<div class ="fragment">
<pre>
          system.solve();
        
</pre>
</div>
<div class = "comment">
Detach the shell matrix from the system since it will go out of
scope.  Nobody should solve the system outside this function.
</div>

<div class ="fragment">
<pre>
          system.detach_shell_matrix();
        
</pre>
</div>
<div class = "comment">
Print a nice message.
</div>

<div class ="fragment">
<pre>
          std::cout &lt;&lt; "Solved linear system in " &lt;&lt; system.n_linear_iterations() &lt;&lt; " iterations, residual norm is " &lt;&lt; system.final_linear_residual() &lt;&lt; "." &lt;&lt; std::endl;
        
        #if defined(LIBMESH_HAVE_VTK) && !defined(LIBMESH_ENABLE_PARMESH)
</pre>
</div>
<div class = "comment">
Write result to file.
</div>

<div class ="fragment">
<pre>
          VTKIO(mesh).write_equation_systems ("out.pvtu", equation_systems);
        #endif // #ifdef LIBMESH_HAVE_VTK
        
        #endif // #ifndef LIBMESH_ENABLE_AMR
        
          return 0;
        }
        
        
        
</pre>
</div>
<div class = "comment">
This function defines the assembly routine.  It is responsible for
computing the proper matrix entries for the element stiffness
matrices and right-hand sides.
</div>

<div class ="fragment">
<pre>
        void assemble (EquationSystems& es,
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
          libmesh_assert_equal_to (system_name, "System");
        
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
          LinearImplicitSystem & system =
            es.get_system&lt;LinearImplicitSystem&gt; ("System");
        
</pre>
</div>
<div class = "comment">
Get the Finite Element type for the first (and only)
variable in the system.
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
const std::vector<Point>& qface_points = fe_face->get_xyz();


<br><br>A reference to the \p DofMap object for this system.  The \p DofMap
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
Analogous data structures for thw two vectors v and w that form
the tensor shell matrix.
</div>

<div class ="fragment">
<pre>
          DenseVector&lt;Number&gt; Ve;
          DenseVector&lt;Number&gt; We;
        
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
              Ve.resize (dof_indices.size());
              We.resize (dof_indices.size());
        
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
                                        phi[i][qp]
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
Stiffness matrix
</div>

<div class ="fragment">
<pre>
                                              (dphi[i][qp]*dphi[j][qp])
                                              );
                        }
        
</pre>
</div>
<div class = "comment">
V and W are the same for this example.
</div>

<div class ="fragment">
<pre>
                      Ve(i) += JxW[qp]*(
                                        phi[i][qp]
                                        );
                      We(i) += JxW[qp]*(
                                        phi[i][qp]
                                        );
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
We have now built the element matrix and RHS vector in terms
of the element degrees of freedom.  However, it is possible
that some of the element DOFs are constrained to enforce
solution continuity, i.e. they are not really "free".  We need
to constrain those DOFs in terms of non-constrained DOFs to
ensure a continuous solution.  The
\p DofMap::constrain_element_matrix_and_vector() method does
just that.


<br><br>However, constraining both the sparse matrix (and right hand
side) plus the rank 1 matrix is tricky.  The dof_indices
vector has to be backuped for that because the constraining
functions modify it.


<br><br></div>

<div class ="fragment">
<pre>
              std::vector&lt;dof_id_type&gt; dof_indices_backup(dof_indices);
              dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
              dof_indices = dof_indices_backup;
              dof_map.constrain_element_dyad_matrix(Ve,We,dof_indices);
        
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
              system.get_matrix("Preconditioner").add_matrix (Ke, dof_indices);
              system.rhs-&gt;add_vector    (Fe, dof_indices);
              system.get_vector("v").add_vector(Ve,dof_indices);
              system.get_vector("w").add_vector(We,dof_indices);
            }
</pre>
</div>
<div class = "comment">
Finished computing the sytem matrix and right-hand side.


<br><br>Matrices and vectors must be closed manually.  This is necessary
because the matrix is not directly used as the system matrix (in
which case the solver closes it) but as a part of a shell matrix.
</div>

<div class ="fragment">
<pre>
          system.matrix-&gt;close();
          system.get_matrix("Preconditioner").close();
          system.rhs-&gt;close();
          system.get_vector("v").close();
          system.get_vector("w").close();
        
        #endif // #ifdef LIBMESH_ENABLE_AMR
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file miscellaneous_ex4.C without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;cstdlib&gt; <I><FONT COLOR="#B22222">// *must* precede &lt;cmath&gt; for proper std:abs() on PGI, Sun Studio CC
</FONT></I>  #include &lt;cmath&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/vtk_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sum_shell_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/tensor_shell_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_shell_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_refinement.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/transient_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/vector_value.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble (EquationSystems&amp; es,
  	       <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> !defined(LIBMESH_ENABLE_AMR)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-amr&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
    libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Running: &quot;</FONT></B> &lt;&lt; argv[0];
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">int</FONT></B> i=1; i&lt;argc; i++)
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B> &lt;&lt; argv[i];
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl &lt;&lt; std::endl;
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
  
    Mesh mesh(init.comm());
  
    EquationSystems equation_systems (mesh);
  
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square (mesh,
  				       16,
  				       16,
  				       -1., 1.,
  				       -1., 1.,
  				       QUAD4);
  
    LinearImplicitSystem &amp; system =
      equation_systems.add_system&lt;LinearImplicitSystem&gt;
      (<B><FONT COLOR="#BC8F8F">&quot;System&quot;</FONT></B>);
  
    system.add_variable (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, FIRST);
  
    system.add_vector(<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>,false);
    system.add_vector(<B><FONT COLOR="#BC8F8F">&quot;w&quot;</FONT></B>,false);
  
    system.add_matrix(<B><FONT COLOR="#BC8F8F">&quot;Preconditioner&quot;</FONT></B>);
  
    system.attach_assemble_function (assemble);
  
    equation_systems.init ();
  
    equation_systems.print_info();
  
    equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;
      (<B><FONT COLOR="#BC8F8F">&quot;linear solver maximum iterations&quot;</FONT></B>) = 250;
    equation_systems.parameters.set&lt;Real&gt;
      (<B><FONT COLOR="#BC8F8F">&quot;linear solver tolerance&quot;</FONT></B>) = TOLERANCE;
  
    <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;2; i++)
      {
        MeshRefinement mesh_refinement(mesh);
        <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::element_iterator       elem_it  = mesh.elements_begin();
        <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::element_iterator elem_end = mesh.elements_end();
        <B><FONT COLOR="#A020F0">for</FONT></B> (; elem_it != elem_end; ++elem_it)
  	{
  	  Elem* elem = *elem_it;
  	  <B><FONT COLOR="#A020F0">if</FONT></B>(elem-&gt;active())
  	    {
  	      <B><FONT COLOR="#A020F0">if</FONT></B>((elem-&gt;id()%20)&gt;8)
  		{
  		  elem-&gt;set_refinement_flag(Elem::REFINE);
  		}
  	      <B><FONT COLOR="#A020F0">else</FONT></B>
  		{
  		  elem-&gt;set_refinement_flag(Elem::DO_NOTHING);
  		}
  	    }
  	  <B><FONT COLOR="#A020F0">else</FONT></B>
  	    {
  	      elem-&gt;set_refinement_flag(Elem::INACTIVE);
  	    }
  	}
        mesh_refinement.refine_elements();
        equation_systems.reinit();
      }
  
    equation_systems.print_info();
  
    system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>).init(system.n_dofs(), system.n_local_dofs());
    system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;w&quot;</FONT></B>).init(system.n_dofs(), system.n_local_dofs());
  
    SumShellMatrix&lt;Number&gt; shellMatrix(system.comm());
    TensorShellMatrix&lt;Number&gt; shellMatrix0(system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>),system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;w&quot;</FONT></B>));
    shellMatrix.matrices.push_back(&amp;shellMatrix0);
    SparseShellMatrix&lt;Number&gt; shellMatrix1(*system.matrix);
    shellMatrix.matrices.push_back(&amp;shellMatrix1);
  
    system.attach_shell_matrix(&amp;shellMatrix);
  
    system.get_matrix(<B><FONT COLOR="#BC8F8F">&quot;Preconditioner&quot;</FONT></B>).zero();
  
    system.solve();
  
    system.detach_shell_matrix();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Solved linear system in &quot;</FONT></B> &lt;&lt; system.n_linear_iterations() &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; iterations, residual norm is &quot;</FONT></B> &lt;&lt; system.final_linear_residual() &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.&quot;</FONT></B> &lt;&lt; std::endl;
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> defined(LIBMESH_HAVE_VTK) &amp;&amp; !defined(LIBMESH_ENABLE_PARMESH)
    VTKIO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;out.pvtu&quot;</FONT></B>, equation_systems);
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_VTK
</FONT></I>  
  #endif <I><FONT COLOR="#B22222">// #ifndef LIBMESH_ENABLE_AMR
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble (EquationSystems&amp; es,
  	       <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
  #ifdef LIBMESH_ENABLE_AMR
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;System&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    LinearImplicitSystem &amp; system =
      es.get_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;System&quot;</FONT></B>);
  
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
  
  
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
  
    DenseMatrix&lt;Number&gt; Ke;
    DenseVector&lt;Number&gt; Fe;
  
    DenseVector&lt;Number&gt; Ve;
    DenseVector&lt;Number&gt; We;
  
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
  
        Fe.resize (dof_indices.size());
        Ve.resize (dof_indices.size());
        We.resize (dof_indices.size());
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
          {
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
              {
                Fe(i) += JxW[qp]*(
                                  phi[i][qp]
                                  );
  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
                  {
                    Ke(i,j) += JxW[qp]*(
  				      (dphi[i][qp]*dphi[j][qp])
                                        );
                  }
  
                Ve(i) += JxW[qp]*(
                                  phi[i][qp]
                                  );
                We(i) += JxW[qp]*(
                                  phi[i][qp]
                                  );
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
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;psi.size(); i++)
                      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;psi.size(); j++)
                        Ke(i,j) += penalty*JxW_face[qp]*psi[i][qp]*psi[j][qp];
                  }
              }
        }
  
  
  
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices_backup(dof_indices);
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        dof_indices = dof_indices_backup;
        dof_map.constrain_element_dyad_matrix(Ve,We,dof_indices);
  
        system.matrix-&gt;add_matrix (Ke, dof_indices);
        system.get_matrix(<B><FONT COLOR="#BC8F8F">&quot;Preconditioner&quot;</FONT></B>).add_matrix (Ke, dof_indices);
        system.rhs-&gt;add_vector    (Fe, dof_indices);
        system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>).add_vector(Ve,dof_indices);
        system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;w&quot;</FONT></B>).add_vector(We,dof_indices);
      }
  
    system.matrix-&gt;close();
    system.get_matrix(<B><FONT COLOR="#BC8F8F">&quot;Preconditioner&quot;</FONT></B>).close();
    system.rhs-&gt;close();
    system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>).close();
    system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;w&quot;</FONT></B>).close();
  
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_ENABLE_AMR
</FONT></I>  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/miscellaneous/miscellaneous_ex4'
***************************************************************
* Running Example miscellaneous_ex4:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
Running: /net/spark/workspace/roystgnr/libmesh/git/devel/examples/miscellaneous/miscellaneous_ex4/.libs/lt-example-devel -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc

 EquationSystems
  n_systems()=1
   System #0, "System"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=289
    n_local_dofs()=82
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=3
    n_matrices()=2
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.99308
      Average Off-Processor Bandwidth <= 0.719723
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 7
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

 EquationSystems
  n_systems()=1
   System #0, "System"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=2093
    n_local_dofs()=546
    n_constrained_dofs()=458
    n_local_constrained_dofs()=109
    n_vectors()=3
    n_matrices()=2
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 9.41806
      Average Off-Processor Bandwidth <= 0.300048
      Maximum  On-Processor Bandwidth <= 21
      Maximum Off-Processor Bandwidth <= 8
    DofMap Constraints
      Number of DoF Constraints = 458
      Average DoF Constraint Length= 2
      Number of Node Constraints = 909
      Maximum Node Constraint Length= 5
      Average Node Constraint Length= 2.51155

Solved linear system in 26 iterations, residual norm is 1.11367e-06.
*** Warning, This code is untested, experimental, or likely to see future API changes: ../src/mesh/vtk_io.C, line 358, compiled Apr 19 2013 at 11:42:41 ***

 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:51:02 2013                                                                          |
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
| libMesh Performance: Alive time=0.342003, Active time=0.256288                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     3         0.0038      0.001262    0.0047      0.001565    1.48     1.83     |
|   build_constraint_matrix()        896       0.0023      0.000003    0.0023      0.000003    0.89     0.89     |
|   build_sparsity()                 3         0.0035      0.001182    0.0081      0.002686    1.38     3.14     |
|   cnstrn_elem_dyad_mat()           448       0.0008      0.000002    0.0008      0.000002    0.33     0.33     |
|   cnstrn_elem_mat_vec()            448       0.0081      0.000018    0.0081      0.000018    3.17     3.17     |
|   create_dof_constraints()         3         0.0129      0.004295    0.0253      0.008441    5.03     9.88     |
|   distribute_dofs()                3         0.0082      0.002721    0.0222      0.007392    3.19     8.65     |
|   dof_indices()                    4672      0.0200      0.000004    0.0200      0.000004    7.79     7.79     |
|   enforce_constraints_exactly()    2         0.0005      0.000237    0.0005      0.000237    0.18     0.18     |
|   old_dof_indices()                1234      0.0062      0.000005    0.0062      0.000005    2.42     2.42     |
|   prepare_send_list()              3         0.0000      0.000006    0.0000      0.000006    0.01     0.01     |
|   reinit()                         3         0.0125      0.004178    0.0125      0.004178    4.89     4.89     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0014      0.001442    0.0042      0.004216    0.56     1.65     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        488       0.0016      0.000003    0.0016      0.000003    0.64     0.64     |
|   init_shape_functions()           41        0.0001      0.000002    0.0001      0.000002    0.03     0.03     |
|   inverse_map()                    3532      0.0065      0.000002    0.0065      0.000002    2.54     2.54     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             488       0.0014      0.000003    0.0014      0.000003    0.55     0.55     |
|   compute_face_map()               40        0.0003      0.000007    0.0006      0.000014    0.10     0.23     |
|   init_face_shape_functions()      1         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   init_reference_to_physical_map() 41        0.0002      0.000005    0.0002      0.000005    0.08     0.08     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           6156      0.0078      0.000001    0.0078      0.000001    3.06     3.06     |
|   init()                           4         0.0018      0.000446    0.0018      0.000446    0.70     0.70     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       2         0.0005      0.000230    0.0008      0.000411    0.18     0.32     |
|   find_neighbors()                 3         0.0109      0.003629    0.0114      0.003801    4.25     4.45     |
|   renumber_nodes_and_elem()        8         0.0012      0.000144    0.0012      0.000144    0.45     0.45     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        4         0.0122      0.003057    0.0122      0.003057    4.77     4.77     |
|   find_global_indices()            4         0.0020      0.000498    0.0159      0.003983    0.78     6.22     |
|   parallel_sort()                  4         0.0006      0.000146    0.0011      0.000284    0.23     0.44     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0076      0.007561    0.0119      0.011867    2.95     4.63     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              2         0.0004      0.000212    0.0004      0.000219    0.17     0.17     |
|   _refine_elements()               4         0.0116      0.002905    0.0325      0.008130    4.53     12.69    |
|   add_point()                      6156      0.0116      0.000002    0.0201      0.000003    4.51     7.83     |
|   make_coarsening_compatible()     4         0.0044      0.001109    0.0044      0.001109    1.73     1.73     |
|   make_flags_parallel_consistent() 6         0.0050      0.000833    0.0069      0.001154    1.95     2.70     |
|   make_refinement_compatible()     9         0.0010      0.000107    0.0010      0.000113    0.38     0.40     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0003      0.000274    0.0003      0.000274    0.11     0.11     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      3         0.0287      0.009573    0.0432      0.014393    11.21    16.85    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      17        0.0005      0.000030    0.0006      0.000033    0.20     0.22     |
|   max(bool)                        19        0.0003      0.000016    0.0003      0.000016    0.12     0.12     |
|   max(scalar)                      465       0.0020      0.000004    0.0020      0.000004    0.78     0.78     |
|   max(vector)                      112       0.0008      0.000007    0.0023      0.000021    0.32     0.91     |
|   min(bool)                        578       0.0023      0.000004    0.0023      0.000004    0.91     0.91     |
|   min(scalar)                      451       0.0038      0.000009    0.0038      0.000009    1.50     1.50     |
|   min(vector)                      112       0.0009      0.000008    0.0024      0.000021    0.36     0.92     |
|   probe()                          168       0.0007      0.000004    0.0007      0.000004    0.29     0.29     |
|   receive()                        168       0.0006      0.000003    0.0013      0.000008    0.23     0.52     |
|   send()                           168       0.0003      0.000002    0.0003      0.000002    0.13     0.13     |
|   send_receive()                   176       0.0008      0.000005    0.0027      0.000015    0.31     1.06     |
|   sum()                            39        0.0008      0.000019    0.0010      0.000025    0.30     0.38     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           168       0.0002      0.000001    0.0002      0.000001    0.07     0.07     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         3         0.0030      0.001002    0.0036      0.001214    1.17     1.42     |
|   set_parent_processor_ids()       3         0.0013      0.000426    0.0013      0.000426    0.50     0.50     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.0185      0.018520    0.0185      0.018520    7.23     7.23     |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       2         0.0029      0.001444    0.0112      0.005609    1.13     4.38     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0162      0.016213    0.0344      0.034434    6.33     13.44    |
|   project_vector()                 2         0.0024      0.001191    0.0175      0.008746    0.93     6.83     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            27374     0.2563                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example miscellaneous_ex4:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/miscellaneous/miscellaneous_ex4'
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
