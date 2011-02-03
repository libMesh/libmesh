<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex20",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 20 - Using a shell matrix</h1>

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
        #include "exodusII_io.h"
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
        #elif !defined(LIBMESH_HAVE_PETSC)
          libmesh_example_assert(false, "--enable-petsc");
        #else
        
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
          
        #ifdef LIBMESH_HAVE_EXODUS_API
</pre>
</div>
<div class = "comment">
Write result to file.
</div>

<div class ="fragment">
<pre>
          {
            OStringStream file_name;
            
            file_name &lt;&lt; "out.exd";
            ExodusII_IO(mesh).write_equation_systems (file_name.str(),
        					equation_systems);
          }
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
        
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
right-hand-side vector.  The \p PetscMatrix::add_matrix()
and \p PetscVector::add_vector() members do this for us.
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
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
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
  #elif !defined(LIBMESH_HAVE_PETSC)
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
  
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
    
  #ifdef LIBMESH_HAVE_EXODUS_API
    {
      OStringStream file_name;
      
      file_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out.exd&quot;</FONT></B>;
      ExodusII_IO(mesh).write_equation_systems (file_name.str(),
  					equation_systems);
    }
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
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
***************************************************************
* Running Example  mpirun -np 2 ./ex20-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running: ./ex20-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 EquationSystems
  n_systems()=1
   System "System"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=289
    n_local_dofs()=153
    n_constrained_dofs()=0
    n_vectors()=3

 EquationSystems
  n_systems()=1
   System "System"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=2093
    n_local_dofs()=1061
    n_constrained_dofs()=458
    n_vectors()=3

Solved linear system in 23 iterations, residual norm is 1.13645e-07.
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex20-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Thu Feb  3 12:11:24 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           9.458e-02      1.02004   9.365e-02
Objects:              1.330e+02      1.00000   1.330e+02
Flops:                8.686e+06      1.04826   8.486e+06  1.697e+07
Flops/sec:            9.185e+07      1.02766   9.061e+07  1.812e+08
MPI Messages:         7.700e+01      1.00000   7.700e+01  1.540e+02
MPI Message Lengths:  7.142e+04      1.01674   9.199e+02  1.417e+05
MPI Reductions:       3.410e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 9.3621e-02 100.0%  1.6973e+07 100.0%  1.540e+02 100.0%  9.199e+02      100.0%  2.670e+02  78.3% 

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

VecDot                24 1.0 1.3447e-04 1.2 5.09e+04 1.0 0.0e+00 0.0e+00 2.4e+01  0  1  0  0  7   0  1  0  0  9   747
VecMDot               23 1.0 2.0862e-04 1.0 5.85e+05 1.0 0.0e+00 0.0e+00 2.3e+01  0  7  0  0  7   0  7  0  0  9  5535
VecNorm               25 1.0 2.3031e-04 1.0 5.30e+04 1.0 0.0e+00 0.0e+00 2.5e+01  0  1  0  0  7   0  1  0  0  9   454
VecScale              24 1.0 2.3603e-05 1.3 2.55e+04 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2128
VecCopy                5 1.0 4.7684e-06 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                71 1.0 3.0041e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               26 1.0 4.3392e-05 1.1 5.30e+04 1.0 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  2412
VecMAXPY              24 1.0 1.5926e-04 1.0 6.34e+05 1.0 0.0e+00 0.0e+00 0.0e+00  0  7  0  0  0   0  7  0  0  0  7859
VecAssemblyBegin      22 1.0 1.9670e-04 1.1 0.00e+00 0.0 6.0e+00 5.5e+02 4.8e+01  0  0  4  2 14   0  0  4  2 18     0
VecAssemblyEnd        22 1.0 2.1219e-05 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       38 1.0 1.3590e-04 1.0 0.00e+00 0.0 6.6e+01 8.3e+02 0.0e+00  0  0 43 39  0   0  0 43 39  0     0
VecScatterEnd         38 1.0 5.6982e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize          24 1.0 2.6345e-04 1.0 7.64e+04 1.0 0.0e+00 0.0e+00 2.4e+01  0  1  0  0  7   0  1  0  0  9   572
MatMult               24 1.0 2.1708e-03 1.4 5.92e+05 1.0 4.8e+01 4.1e+02 9.6e+01  2  7 31 14 28   2  7 31 14 36   534
MatMultAdd            24 1.0 5.4908e-04 1.0 4.92e+05 1.0 4.8e+01 4.1e+02 0.0e+00  1  6 31 14  0   1  6 31 14  0  1752
MatSolve              24 1.0 3.0849e-03 1.1 3.75e+06 1.1 0.0e+00 0.0e+00 0.0e+00  3 43  0  0  0   3 43  0  0  0  2361
MatLUFactorNum         1 1.0 2.8989e-03 1.0 3.05e+06 1.0 0.0e+00 0.0e+00 0.0e+00  3 35  0  0  0   3 35  0  0  0  2057
MatILUFactorSym        1 1.0 6.2120e-03 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+00  6  0  0  0  0   6  0  0  0  0     0
MatAssemblyBegin      26 1.0 9.7370e-04 1.2 0.00e+00 0.0 6.0e+00 2.2e+03 5.2e+01  1  0  4  9 15   1  0  4  9 19     0
MatAssemblyEnd        26 1.0 8.6880e-04 1.0 0.00e+00 0.0 8.0e+00 1.0e+02 3.8e+01  1  0  5  1 11   1  0  5  1 14     0
MatGetRowIJ            1 1.0 0.0000e+00 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 3.9101e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  1   0  0  0  0  1     0
MatZeroEntries        14 1.0 8.4162e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        23 1.0 3.8290e-04 1.0 1.17e+06 1.0 0.0e+00 0.0e+00 2.3e+01  0 14  0  0  7   0 14  0  0  9  6033
KSPSetup               2 1.0 3.5048e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 1.5227e-02 1.0 8.69e+06 1.0 4.8e+01 4.1e+02 1.5e+02 16100 31 14 43  16100 31 14 55  1115
PCSetUp                2 1.0 9.3770e-03 1.1 3.05e+06 1.0 0.0e+00 0.0e+00 3.0e+00 10 35  0  0  1  10 35  0  0  1   636
PCSetUpOnBlocks        1 1.0 9.2041e-03 1.1 3.05e+06 1.0 0.0e+00 0.0e+00 3.0e+00 10 35  0  0  1  10 35  0  0  1   648
PCApply               24 1.0 3.2909e-03 1.0 3.75e+06 1.1 0.0e+00 0.0e+00 0.0e+00  3 43  0  0  0   3 43  0  0  0  2213
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    70             70       536312     0
         Vec Scatter    14             14        12152     0
           Index Set    22             22        25076     0
   IS L to G Mapping     3              3         8132     0
              Matrix    20             20      1447452     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 1.38283e-06
Average time for zero size MPI_Send(): 5.48363e-06
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
| Time:           Thu Feb  3 12:11:24 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 -------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.102748, Active time=0.083697                                              |
 -------------------------------------------------------------------------------------------------------------
| Event                           nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                           w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-------------------------------------------------------------------------------------------------------------|
|                                                                                                             |
|                                                                                                             |
| DofMap                                                                                                      |
|   add_neighbors_to_send_list()  3         0.0007      0.000223    0.0007      0.000238    0.80     0.85     |
|   build_constraint_matrix()     1794      0.0015      0.000001    0.0015      0.000001    1.80     1.80     |
|   cnstrn_elem_dyad_mat()        897       0.0002      0.000000    0.0002      0.000000    0.29     0.29     |
|   cnstrn_elem_mat_vec()         897       0.0044      0.000005    0.0044      0.000005    5.30     5.30     |
|   compute_sparsity()            3         0.0043      0.001439    0.0053      0.001759    5.16     6.30     |
|   create_dof_constraints()      3         0.0036      0.001204    0.0052      0.001733    4.31     6.21     |
|   distribute_dofs()             3         0.0019      0.000632    0.0049      0.001618    2.26     5.80     |
|   dof_indices()                 8234      0.0025      0.000000    0.0025      0.000000    3.00     3.00     |
|   enforce_constraints_exactly() 2         0.0002      0.000109    0.0002      0.000109    0.26     0.26     |
|   old_dof_indices()             2466      0.0006      0.000000    0.0006      0.000000    0.75     0.75     |
|   prepare_send_list()           3         0.0000      0.000004    0.0000      0.000004    0.01     0.01     |
|   reinit()                      3         0.0028      0.000945    0.0028      0.000945    3.39     3.39     |
|                                                                                                             |
| FE                                                                                                          |
|   compute_affine_map()          949       0.0005      0.000001    0.0005      0.000001    0.63     0.63     |
|   compute_face_map()            52        0.0001      0.000002    0.0002      0.000004    0.13     0.27     |
|   compute_shape_functions()     949       0.0002      0.000000    0.0002      0.000000    0.29     0.29     |
|   init_face_shape_functions()   12        0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|   init_shape_functions()        53        0.0001      0.000002    0.0001      0.000002    0.14     0.14     |
|   inverse_map()                 3614      0.0017      0.000000    0.0017      0.000000    2.07     2.07     |
|                                                                                                             |
| LocationMap                                                                                                 |
|   find()                        6156      0.0021      0.000000    0.0021      0.000000    2.48     2.48     |
|   init()                        4         0.0007      0.000175    0.0007      0.000175    0.84     0.84     |
|                                                                                                             |
| Mesh                                                                                                        |
|   contract()                    2         0.0002      0.000084    0.0003      0.000138    0.20     0.33     |
|   find_neighbors()              3         0.0046      0.001550    0.0047      0.001559    5.56     5.59     |
|   renumber_nodes_and_elem()     8         0.0003      0.000043    0.0003      0.000043    0.41     0.41     |
|                                                                                                             |
| MeshCommunication                                                                                           |
|   compute_hilbert_indices()     4         0.0049      0.001213    0.0049      0.001213    5.80     5.80     |
|   find_global_indices()         4         0.0008      0.000210    0.0063      0.001582    1.00     7.56     |
|   parallel_sort()               4         0.0004      0.000107    0.0005      0.000126    0.51     0.60     |
|                                                                                                             |
| MeshRefinement                                                                                              |
|   _coarsen_elements()           2         0.0001      0.000075    0.0002      0.000079    0.18     0.19     |
|   _refine_elements()            4         0.0057      0.001433    0.0130      0.003249    6.85     15.53    |
|   add_point()                   6156      0.0041      0.000001    0.0066      0.000001    4.90     7.84     |
|   make_coarsening_compatible()  2         0.0008      0.000387    0.0008      0.000387    0.92     0.92     |
|   make_refinement_compatible()  5         0.0002      0.000040    0.0002      0.000044    0.24     0.26     |
|                                                                                                             |
| MeshTools::Generation                                                                                       |
|   build_cube()                  1         0.0001      0.000117    0.0001      0.000117    0.14     0.14     |
|                                                                                                             |
| MetisPartitioner                                                                                            |
|   partition()                   3         0.0035      0.001168    0.0090      0.002993    4.19     10.73    |
|                                                                                                             |
| Parallel                                                                                                    |
|   allgather()                   13        0.0001      0.000004    0.0001      0.000004    0.07     0.07     |
|   max(bool)                     13        0.0003      0.000024    0.0003      0.000024    0.37     0.37     |
|   max(scalar)                   8         0.0001      0.000009    0.0001      0.000009    0.08     0.08     |
|   max(vector)                   4         0.0000      0.000002    0.0000      0.000002    0.01     0.01     |
|   min(bool)                     7         0.0000      0.000005    0.0000      0.000005    0.04     0.04     |
|   min(vector)                   4         0.0000      0.000005    0.0000      0.000005    0.02     0.02     |
|   probe()                       26        0.0001      0.000002    0.0001      0.000002    0.08     0.08     |
|   receive()                     26        0.0001      0.000003    0.0001      0.000005    0.08     0.16     |
|   send()                        26        0.0000      0.000001    0.0000      0.000001    0.03     0.03     |
|   send_receive()                34        0.0001      0.000002    0.0002      0.000007    0.08     0.30     |
|   sum()                         15        0.0001      0.000010    0.0001      0.000010    0.18     0.18     |
|   wait()                        26        0.0000      0.000001    0.0000      0.000001    0.02     0.02     |
|                                                                                                             |
| Partitioner                                                                                                 |
|   set_node_processor_ids()      3         0.0010      0.000339    0.0011      0.000355    1.21     1.27     |
|   set_parent_processor_ids()    3         0.0002      0.000070    0.0002      0.000070    0.25     0.25     |
|                                                                                                             |
| PetscLinearSolver                                                                                           |
|   solve()                       1         0.0156      0.015564    0.0156      0.015564    18.60    18.60    |
|                                                                                                             |
| ProjectVector                                                                                               |
|   operator()                    2         0.0023      0.001154    0.0040      0.002020    2.76     4.83     |
|                                                                                                             |
| System                                                                                                      |
|   assemble()                    1         0.0081      0.008087    0.0164      0.016418    9.66     19.62    |
|   project_vector()              2         0.0014      0.000684    0.0060      0.003014    1.63     7.20     |
 -------------------------------------------------------------------------------------------------------------
| Totals:                         32509     0.0837                                          100.00            |
 -------------------------------------------------------------------------------------------------------------

 
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
