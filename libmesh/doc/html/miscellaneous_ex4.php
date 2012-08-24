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
        #include &lt;cmath&gt;
        
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
        #include "vtk_io.h"
        #include "equation_systems.h"
        #include "fe.h"
        #include "quadrature_gauss.h"
        #include "dof_map.h"
        #include "sparse_matrix.h"
        #include "numeric_vector.h"
        #include "dense_matrix.h"
        #include "dense_vector.h"
        #include "mesh_generation.h"
        #include "sum_shell_matrix.h"
        #include "tensor_shell_matrix.h"
        #include "sparse_shell_matrix.h"
        #include "mesh_refinement.h"
        
        #include "getpot.h"
        
</pre>
</div>
<div class = "comment">
Some (older) compilers do not offer full stream 
functionality, \p OStringStream works around this.
Check example 9 for details.
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
        #include "transient_system.h"
        #include "linear_implicit_system.h"
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
Create a mesh
</div>

<div class ="fragment">
<pre>
          Mesh mesh;
        
</pre>
</div>
<div class = "comment">
Create an equation systems object.
</div>

<div class ="fragment">
<pre>
          EquationSystems equation_systems (mesh);
          MeshRefinement mesh_refinement (mesh);
        
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
          SumShellMatrix&lt;Number&gt; shellMatrix;
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
          libmesh_assert (system_name == "System");
          
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
          std::vector&lt;unsigned int&gt; dof_indices;
        
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
              std::vector&lt;unsigned int&gt; dof_indices_backup(dof_indices);
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
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
   
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;cmath&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;vtk_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sum_shell_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;tensor_shell_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_shell_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_refinement.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;o_string_stream.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;transient_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;vector_value.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  
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
  
    Mesh mesh;
  
    EquationSystems equation_systems (mesh);
    MeshRefinement mesh_refinement (mesh);
  
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
  
    SumShellMatrix&lt;Number&gt; shellMatrix;
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
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;System&quot;</FONT></B>);
    
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
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
  
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
  
        
  
  
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices_backup(dof_indices);
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
Linking miscellaneous_ex4-opt...
***************************************************************
* Running Example  mpirun -np 6 ./miscellaneous_ex4-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running: ./miscellaneous_ex4-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 EquationSystems
  n_systems()=1
   System #0, "System"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=289
    n_local_dofs()=57
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=3
    n_matrices()=2
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.45614
      Average Off-Processor Bandwidth <= 0.824561
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 5
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
    n_local_dofs()=377
    n_constrained_dofs()=458
    n_local_constrained_dofs()=80
    n_vectors()=3
    n_matrices()=2
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 9.15385
      Average Off-Processor Bandwidth <= 0.437666
      Maximum  On-Processor Bandwidth <= 18
      Maximum Off-Processor Bandwidth <= 6
    DofMap Constraints
      Number of DoF Constraints = 458
      Average DoF Constraint Length= 2
      Number of Node Constraints = 909
      Maximum Node Constraint Length= 5
      Average Node Constraint Length= 2.51155

Solved linear system in 34 iterations, residual norm is 1.05719e-07.
*** Warning, This code is untested, experimental, or likely to see future API changes: /h2/roystgnr/libmesh/svn/include/mesh/vtk_io.h, line 154, compiled Aug 24 2012 at 15:12:41 ***
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./miscellaneous_ex4-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:21:51 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           2.189e-01      1.26748   1.815e-01
Objects:              1.230e+02      1.00000   1.230e+02
Flops:                4.029e+06      1.68868   3.197e+06  1.918e+07
Flops/sec:            2.331e+07      1.69798   1.772e+07  1.063e+08
MPI Messages:         2.860e+02      1.67251   2.252e+02  1.351e+03
MPI Message Lengths:  7.117e+04      1.25079   2.774e+02  3.747e+05
MPI Reductions:       6.180e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.8149e-01 100.0%  1.9181e+07 100.0%  1.351e+03 100.0%  2.774e+02      100.0%  5.440e+02  88.0% 

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

VecDot                36 1.0 2.8708e-03 3.9 2.71e+04 1.2 0.0e+00 0.0e+00 3.6e+01  1  1  0  0  6   1  1  0  0  7    52
VecMDot               34 1.0 2.0499e-03 1.2 3.58e+05 1.2 0.0e+00 0.0e+00 3.4e+01  1 10  0  0  6   1 10  0  0  6   969
VecNorm               37 1.0 2.0254e-03 1.3 2.79e+04 1.2 0.0e+00 0.0e+00 3.7e+01  1  1  0  0  6   1  1  0  0  7    76
VecScale              36 1.0 1.7834e-04 5.2 1.36e+04 1.2 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   423
VecCopy               13 1.0 1.4544e-05 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                96 1.0 3.2830e-0410.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               40 1.0 2.1315e-04 3.4 2.94e+04 1.2 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0   766
VecMAXPY              36 1.0 4.2510e-04 3.5 3.84e+05 1.2 0.0e+00 0.0e+00 0.0e+00  0 11  0  0  0   0 11  0  0  0  5012
VecAssemblyBegin      96 1.0 2.0976e-02 3.2 0.00e+00 0.0 7.2e+01 2.1e+02 2.7e+02  8  0  5  4 44   8  0  5  4 50     0
VecAssemblyEnd        96 1.0 7.9870e-05 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       45 1.0 3.1018e-04 2.3 0.00e+00 0.0 8.4e+02 2.9e+02 0.0e+00  0  0 62 65  0   0  0 62 65  0     0
VecScatterEnd         45 1.0 3.1514e-03 3.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
VecNormalize          36 1.0 2.0847e-03 1.1 4.07e+04 1.2 0.0e+00 0.0e+00 3.6e+01  1  1  0  0  6   1  1  0  0  7   108
MatMult               36 1.0 2.2015e-02 1.6 3.16e+05 1.2 6.5e+02 1.6e+02 3.6e+02 11  9 48 28 58  11  9 48 28 66    79
MatMultAdd            36 1.0 1.8768e-03 4.4 2.63e+05 1.2 6.5e+02 1.6e+02 0.0e+00  1  8 48 28  0   1  8 48 28  0   769
MatSolve              36 1.0 5.0514e-03 4.5 1.96e+06 1.8 0.0e+00 0.0e+00 0.0e+00  1 47  0  0  0   1 47  0  0  0  1785
MatLUFactorNum         1 1.0 9.7322e-04 2.1 9.78e+05 2.4 0.0e+00 0.0e+00 0.0e+00  0 21  0  0  0   0 21  0  0  0  4172
MatILUFactorSym        1 1.0 2.3172e-03 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+00  1  0  0  0  0   1  0  0  0  0     0
MatAssemblyBegin      38 1.0 6.7410e-03 2.1 0.00e+00 0.0 7.2e+01 6.4e+02 7.6e+01  3  0  5 12 12   3  0  5 12 14     0
MatAssemblyEnd        38 1.0 5.7118e-03 1.0 0.00e+00 0.0 7.2e+01 4.3e+01 5.0e+01  3  0  5  1  8   3  0  5  1  9     0
MatGetRowIJ            1 1.0 1.1921e-0512.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 5.7936e-05 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries        14 1.0 2.3127e-05 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        34 1.0 2.3520e-03 1.2 7.16e+05 1.2 0.0e+00 0.0e+00 3.4e+01  1 21  0  0  6   1 21  0  0  6  1690
KSPSetup               2 1.0 5.0068e-05 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 3.0088e-02 1.0 4.03e+06 1.7 6.5e+02 1.6e+02 4.3e+02 16100 48 28 70  17100 48 28 80   637
PCSetUp                2 1.0 4.7731e-03 2.4 9.78e+05 2.4 0.0e+00 0.0e+00 3.0e+00  2 21  0  0  0   2 21  0  0  1   851
PCSetUpOnBlocks        1 1.0 3.4149e-03 2.0 9.78e+05 2.4 0.0e+00 0.0e+00 3.0e+00  1 21  0  0  0   1 21  0  0  1  1189
PCApply               36 1.0 6.1200e-03 4.2 1.96e+06 1.8 0.0e+00 0.0e+00 0.0e+00  1 47  0  0  0   1 47  0  0  0  1473
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    70             70       287200     0
         Vec Scatter     9              9         7812     0
           Index Set    17             17        14320     0
   IS L to G Mapping     3              3         3952     0
              Matrix    20             20       478028     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.282e-05
Average time for zero size MPI_Send(): 5.45184e-05
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
Configure run at: Sat May 19 03:47:23 2012
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared=1 --with-shared-libraries=1 --with-mpi-dir=/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid --with-mumps=true --download-mumps=1 --with-parmetis=true --download-parmetis=1 --with-superlu=true --download-superlu=1 --with-superludir=true --download-superlu_dist=1 --with-blacs=true --download-blacs=1 --with-scalapack=true --download-scalapack=1 --with-hypre=true --download-hypre=1 --with-blas-lib="[/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_intel_lp64.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_sequential.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_core.so]" --with-lapack-lib=/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t/libmkl_solver_lp64_sequential.a
-----------------------------------------
Libraries compiled on Sat May 19 03:47:23 CDT 2012 on daedalus 
Machine characteristics: Linux daedalus 2.6.32-34-generic #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011 x86_64 GNU/Linux 
Using PETSc directory: /org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5
Using PETSc arch: intel-11.1-lucid-mpich2-1.4.1-cxx-opt
-----------------------------------------
Using C compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpicxx -O3   -fPIC   
Using Fortran compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpif90 -fPIC -O3    
-----------------------------------------
Using include paths: -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/include  
------------------------------------------
Using C linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpicxx -O3 
Using Fortran linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/bin/mpif90 -fPIC -O3  
Using libraries: -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -lpetsc       -lX11 -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/intel-11.1-lucid-mpich2-1.4.1-cxx-opt/lib -lHYPRE -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lblacs -lsuperlu_dist_2.4 -lparmetis -lmetis -lsuperlu_4.0 -Wl,-rpath,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t -L/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-intel-11.1-lucid/lib/em64t -lmkl_solver_lp64_sequential -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -ldl -Wl,-rpath,/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/lib -L/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.4.1-intel-11.1-lucid/lib -lmpich -lopa -lmpl -lrt -lpthread -Wl,-rpath,/opt/intel/Compiler/11.1/073/lib/intel64 -L/opt/intel/Compiler/11.1/073/lib/intel64 -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -L/usr/lib/gcc/x86_64-linux-gnu/4.4.3 -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -lmpichf90 -lifport -lifcore -lm -lm -lmpichcxx -lstdc++ -lmpichcxx -lstdc++ -ldl -lmpich -lopa -lmpl -lrt -lpthread -limf -lsvml -lipgo -ldecimal -lgcc_s -lirc -lirc_s -ldl  
------------------------------------------

-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 6                                                |
| Time:           Fri Aug 24 15:21:51 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.323955, Active time=0.15315                                                  |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     3         0.0004      0.000126    0.0004      0.000148    0.25     0.29     |
|   build_constraint_matrix()        598       0.0009      0.000001    0.0009      0.000001    0.57     0.57     |
|   build_sparsity()                 3         0.0039      0.001298    0.0048      0.001612    2.54     3.16     |
|   cnstrn_elem_dyad_mat()           299       0.0002      0.000001    0.0002      0.000001    0.13     0.13     |
|   cnstrn_elem_mat_vec()            299       0.0079      0.000027    0.0079      0.000027    5.18     5.18     |
|   create_dof_constraints()         3         0.0052      0.001743    0.0070      0.002317    3.41     4.54     |
|   distribute_dofs()                3         0.0010      0.000344    0.0151      0.005046    0.67     9.89     |
|   dof_indices()                    5918      0.0017      0.000000    0.0017      0.000000    1.14     1.14     |
|   enforce_constraints_exactly()    2         0.0010      0.000515    0.0010      0.000515    0.67     0.67     |
|   old_dof_indices()                822       0.0002      0.000000    0.0002      0.000000    0.16     0.16     |
|   prepare_send_list()              3         0.0000      0.000004    0.0000      0.000004    0.01     0.01     |
|   reinit()                         3         0.0018      0.000602    0.0018      0.000602    1.18     1.18     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        325       0.0005      0.000002    0.0005      0.000002    0.35     0.35     |
|   init_shape_functions()           27        0.0000      0.000001    0.0000      0.000001    0.02     0.02     |
|   inverse_map()                    3322      0.0012      0.000000    0.0012      0.000000    0.75     0.75     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             325       0.0003      0.000001    0.0003      0.000001    0.19     0.19     |
|   compute_face_map()               26        0.0001      0.000003    0.0002      0.000006    0.05     0.10     |
|   init_face_shape_functions()      1         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   init_reference_to_physical_map() 27        0.0001      0.000003    0.0001      0.000003    0.05     0.05     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           6156      0.0019      0.000000    0.0019      0.000000    1.26     1.26     |
|   init()                           4         0.0005      0.000120    0.0005      0.000120    0.31     0.31     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       2         0.0001      0.000034    0.0002      0.000100    0.05     0.13     |
|   find_neighbors()                 3         0.0030      0.000984    0.0073      0.002425    1.93     4.75     |
|   renumber_nodes_and_elem()        8         0.0005      0.000060    0.0005      0.000060    0.31     0.31     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        4         0.0051      0.001284    0.0051      0.001284    3.35     3.35     |
|   find_global_indices()            4         0.0006      0.000146    0.0172      0.004289    0.38     11.20    |
|   parallel_sort()                  4         0.0023      0.000571    0.0082      0.002041    1.49     5.33     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              2         0.0001      0.000027    0.0002      0.000087    0.04     0.11     |
|   _refine_elements()               4         0.0053      0.001313    0.0286      0.007156    3.43     18.69    |
|   add_point()                      6156      0.0043      0.000001    0.0067      0.000001    2.82     4.35     |
|   make_coarsening_compatible()     2         0.0005      0.000236    0.0005      0.000236    0.31     0.31     |
|   make_refinement_compatible()     5         0.0001      0.000027    0.0008      0.000169    0.09     0.55     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0001      0.000114    0.0001      0.000114    0.07     0.07     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      3         0.0052      0.001733    0.0209      0.006951    3.39     13.62    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      17        0.0035      0.000205    0.0035      0.000205    2.28     2.28     |
|   max(bool)                        13        0.0170      0.001311    0.0170      0.001311    11.13    11.13    |
|   max(scalar)                      6         0.0044      0.000736    0.0044      0.000736    2.88     2.88     |
|   max(vector)                      4         0.0003      0.000086    0.0003      0.000086    0.22     0.22     |
|   min(bool)                        7         0.0007      0.000105    0.0007      0.000105    0.48     0.48     |
|   min(vector)                      4         0.0012      0.000296    0.0012      0.000296    0.77     0.77     |
|   probe()                          130       0.0185      0.000142    0.0185      0.000142    12.05    12.05    |
|   receive()                        130       0.0002      0.000002    0.0187      0.000144    0.14     12.20    |
|   send()                           130       0.0001      0.000001    0.0001      0.000001    0.07     0.07     |
|   send_receive()                   138       0.0003      0.000002    0.0191      0.000139    0.17     12.49    |
|   sum()                            16        0.0062      0.000386    0.0062      0.000386    4.03     4.03     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           130       0.0001      0.000000    0.0001      0.000000    0.04     0.04     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         3         0.0005      0.000160    0.0075      0.002514    0.31     4.93     |
|   set_parent_processor_ids()       3         0.0002      0.000078    0.0002      0.000078    0.15     0.15     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.0312      0.031176    0.0312      0.031176    20.36    20.36    |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       2         0.0009      0.000464    0.0016      0.000781    0.61     1.02     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0074      0.007408    0.0180      0.018008    4.84     11.76    |
|   project_vector()                 2         0.0044      0.002216    0.0072      0.003614    2.89     4.72     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            25104     0.1532                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
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
