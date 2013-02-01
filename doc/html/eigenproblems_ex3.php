<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("eigenproblems_ex3",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file eigenproblems_ex3.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include &lt;fstream&gt;
        
</pre>
</div>
<div class = "comment">
libMesh include files.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh.h"
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/exodusII_io.h"
        #include "libmesh/condensed_eigen_system.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/fe.h"
        #include "libmesh/quadrature_gauss.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/fe_interface.h"
        #include "libmesh/getpot.h"
        
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
Function prototype.  This is the function that will assemble
the eigen system. Here, we will simply assemble a mass matrix.
</div>

<div class ="fragment">
<pre>
        void assemble_matrices(EquationSystems& es,
                               const std::string& system_name);
        
</pre>
</div>
<div class = "comment">
We store the Dirichlet dofs in a set in order to impose the boundary conditions
</div>

<div class ="fragment">
<pre>
        void get_dirichlet_dofs(EquationSystems& es,
                                const std::string& system_name,
                                std::set&lt;unsigned int&gt;& global_dirichlet_dofs_set);
        
        
        int main (int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Initialize libMesh and the dependent libraries.
</div>

<div class ="fragment">
<pre>
          LibMeshInit init (argc, argv);
        
</pre>
</div>
<div class = "comment">
This example uses an ExodusII input file
</div>

<div class ="fragment">
<pre>
        #ifndef LIBMESH_HAVE_EXODUS_API
          libmesh_example_assert(false, "--enable-exodus");
        #endif
        
</pre>
</div>
<div class = "comment">
This example is designed for the SLEPc eigen solver interface.
</div>

<div class ="fragment">
<pre>
        #ifndef LIBMESH_HAVE_SLEPC
          if (libMesh::processor_id() == 0)
            std::cerr &lt;&lt; "ERROR: This example requires libMesh to be\n"
                      &lt;&lt; "compiled with SLEPc eigen solvers support!"
                      &lt;&lt; std::endl;
        
          return 0;
        #else
        
        #ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
</pre>
</div>
<div class = "comment">
SLEPc currently gives us a nasty crash with Real==float
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(false, "--disable-singleprecision");
        #endif
        
        #ifdef LIBMESH_USE_COMPLEX_NUMBERS
</pre>
</div>
<div class = "comment">
SLEPc currently gives us an "inner product not well defined" with
Number==complex
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(false, "--disable-complex");
        #endif
        
</pre>
</div>
<div class = "comment">
Tell the user what we are doing.
</div>

<div class ="fragment">
<pre>
          {
            std::cout &lt;&lt; "Running " &lt;&lt; argv[0];
              
            for (int i=1; i&lt;argc; i++)
              std::cout &lt;&lt; " " &lt;&lt; argv[i];
              
            std::cout &lt;&lt; std::endl &lt;&lt; std::endl;
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
Use GetPot to parse the command line arguments
</div>

<div class ="fragment">
<pre>
          GetPot command_line (argc, argv);
        
</pre>
</div>
<div class = "comment">
Read the mesh name from the command line
</div>

<div class ="fragment">
<pre>
          std::string mesh_name = "";
          if ( command_line.search(1, "-mesh_name") )
            mesh_name = command_line.next(mesh_name);
            
</pre>
</div>
<div class = "comment">
Also, read in the index of the eigenvector that we should plot
(zero-based indexing, as usual!)
</div>

<div class ="fragment">
<pre>
          unsigned int plotting_index = 0;
          if ( command_line.search(1, "-plotting_index") )
            plotting_index = command_line.next(plotting_index);
        
</pre>
</div>
<div class = "comment">
Finally, read in the number of eigenpairs we want to compute!
</div>

<div class ="fragment">
<pre>
          unsigned int n_evals = 0;
          if ( command_line.search(1, "-n_evals") )
            n_evals = command_line.next(n_evals);
        
</pre>
</div>
<div class = "comment">
Append the .e to mesh_name
</div>

<div class ="fragment">
<pre>
          std::ostringstream mesh_name_exodus;
          mesh_name_exodus &lt;&lt; mesh_name &lt;&lt; "_mesh.e";
        
          Mesh mesh;
          mesh.read(mesh_name_exodus.str());
        
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
Create a CondensedEigenSystem named "Eigensystem" and (for convenience)
use a reference to the system we create.
</div>

<div class ="fragment">
<pre>
          CondensedEigenSystem & eigen_system =
            equation_systems.add_system&lt;CondensedEigenSystem&gt; ("Eigensystem");
        
</pre>
</div>
<div class = "comment">
Declare the system variables.
Adds the variable "p" to "Eigensystem".   "p"
will be approximated using second-order approximation.
</div>

<div class ="fragment">
<pre>
          eigen_system.add_variable("p", SECOND);
        
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly
function defined below.
</div>

<div class ="fragment">
<pre>
          eigen_system.attach_assemble_function (assemble_matrices);
        
</pre>
</div>
<div class = "comment">
Set the number of requested eigenpairs \p n_evals and the number
of basis vectors used in the solution algorithm.
</div>

<div class ="fragment">
<pre>
          equation_systems.parameters.set&lt;unsigned int&gt;("eigenpairs")    = n_evals;
          equation_systems.parameters.set&lt;unsigned int&gt;("basis vectors") = n_evals*3;
        
</pre>
</div>
<div class = "comment">
Set the solver tolerance and the maximum number of iterations. 
</div>

<div class ="fragment">
<pre>
          equation_systems.parameters.set&lt;Real&gt;("linear solver tolerance") = pow(TOLERANCE, 5./3.);
          equation_systems.parameters.set&lt;unsigned int&gt;
            ("linear solver maximum iterations") = 1000;
        
</pre>
</div>
<div class = "comment">
Set the type of the problem, here we deal with
a generalized Hermitian problem.
</div>

<div class ="fragment">
<pre>
          eigen_system.set_eigenproblem_type(GHEP);
        
</pre>
</div>
<div class = "comment">
Order the eigenvalues "smallest first"
</div>

<div class ="fragment">
<pre>
          eigen_system.eigen_solver-&gt;set_position_of_spectrum(SMALLEST_MAGNITUDE);
        
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
Pass the Dirichlet dof IDs to the CondensedEigenSystem
</div>

<div class ="fragment">
<pre>
          std::set&lt;unsigned int&gt; dirichlet_dof_ids;
          get_dirichlet_dofs(equation_systems, "Eigensystem", dirichlet_dof_ids);
          eigen_system.initialize_condensed_dofs(dirichlet_dof_ids);
        
</pre>
</div>
<div class = "comment">
Solve the system "Eigensystem".
</div>

<div class ="fragment">
<pre>
          eigen_system.solve();
        
</pre>
</div>
<div class = "comment">
Get the number of converged eigen pairs.
</div>

<div class ="fragment">
<pre>
          unsigned int nconv = eigen_system.get_n_converged();
        
          std::cout &lt;&lt; "Number of converged eigenpairs: " &lt;&lt; nconv
                    &lt;&lt; "\n" &lt;&lt; std::endl;
        
          if (plotting_index &gt; n_evals)
            {
              std::cout &lt;&lt; "WARNING: Solver did not converge for the requested eigenvector!" &lt;&lt; std::endl;
            }
        
</pre>
</div>
<div class = "comment">
write out all of the computed eigenvalues and plot the specified eigenvector
</div>

<div class ="fragment">
<pre>
          std::ostringstream eigenvalue_output_name;
          eigenvalue_output_name &lt;&lt; mesh_name &lt;&lt; "_evals.txt";
          std::ofstream evals_file(eigenvalue_output_name.str().c_str());
        
          for(unsigned int i=0; i&lt;nconv; i++)
          {
            std::pair&lt;Real,Real&gt; eval = eigen_system.get_eigenpair(i);
            
</pre>
</div>
<div class = "comment">
The eigenvalues should be real!
</div>

<div class ="fragment">
<pre>
            libmesh_assert_less (eval.second, TOLERANCE);
            evals_file &lt;&lt; eval.first &lt;&lt; std::endl;
            
</pre>
</div>
<div class = "comment">
plot the specified eigenvector
</div>

<div class ="fragment">
<pre>
            if(i == plotting_index)
            {
        #ifdef LIBMESH_HAVE_EXODUS_API
</pre>
</div>
<div class = "comment">
Write the eigen vector to file.
</div>

<div class ="fragment">
<pre>
              std::ostringstream eigenvector_output_name;
              eigenvector_output_name &lt;&lt; mesh_name &lt;&lt; "_evec.e";
              ExodusII_IO (mesh).write_equation_systems (eigenvector_output_name.str(), equation_systems);
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
            }
          }
          
          evals_file.close();
        
        #endif // LIBMESH_HAVE_SLEPC
        
</pre>
</div>
<div class = "comment">
All done.  
</div>

<div class ="fragment">
<pre>
          return 0;
        }
        
        
        
        void assemble_matrices(EquationSystems& es,
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
          libmesh_assert_equal_to (system_name, "Eigensystem");
        
        #ifdef LIBMESH_HAVE_SLEPC
        
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
The dimension that we are running.
</div>

<div class ="fragment">
<pre>
          const unsigned int dim = mesh.mesh_dimension();
        
</pre>
</div>
<div class = "comment">
Get a reference to our system.
</div>

<div class ="fragment">
<pre>
          EigenSystem & eigen_system = es.get_system&lt;EigenSystem&gt; (system_name);
        
</pre>
</div>
<div class = "comment">
Get a constant reference to the Finite Element type
for the first (and only) variable in the system.
</div>

<div class ="fragment">
<pre>
          FEType fe_type = eigen_system.get_dof_map().variable_type(0);
        
</pre>
</div>
<div class = "comment">
A reference to the two system matrices
</div>

<div class ="fragment">
<pre>
          SparseMatrix&lt;Number&gt;&  matrix_A = *eigen_system.matrix_A;
          SparseMatrix&lt;Number&gt;&  matrix_B = *eigen_system.matrix_B;
        
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
A  Gauss quadrature rule for numerical integration.
Use the default quadrature order.
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
The element Jacobian * quadrature weight at each integration point.   
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
to degree of freedom numbers.
</div>

<div class ="fragment">
<pre>
          const DofMap& dof_map = eigen_system.get_dof_map();
        
</pre>
</div>
<div class = "comment">
The element mass and stiffness matrices.
</div>

<div class ="fragment">
<pre>
          DenseMatrix&lt;Number&gt;   Me;
          DenseMatrix&lt;Number&gt;   Ke;
        
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
matrix and right-hand-side contribution.  In case users
later modify this program to include refinement, we will
be safe and will only consider the active elements;
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
Zero the element matrices before
summing them.  We use the resize member here because
the number of degrees of freedom might have changed from
the last element.  Note that this will be the case if the
element type is different (i.e. the last element was a
triangle, now we are on a quadrilateral).
</div>

<div class ="fragment">
<pre>
              Ke.resize (dof_indices.size(), dof_indices.size());
              Me.resize (dof_indices.size(), dof_indices.size());
        
</pre>
</div>
<div class = "comment">
Now loop over the quadrature points.  This handles
the numeric integration.

<br><br>We will build the element matrix.  This involves
a double loop to integrate the test funcions (i) against
the trial functions (j).
</div>

<div class ="fragment">
<pre>
              for (unsigned int qp=0; qp&lt;qrule.n_points(); qp++)
                for (unsigned int i=0; i&lt;phi.size(); i++)
                  for (unsigned int j=0; j&lt;phi.size(); j++)
                    {
                      Me(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
                      Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
                    }
        
</pre>
</div>
<div class = "comment">
On an unrefined mesh, constrain_element_matrix does
nothing.  If this assembly function is ever repurposed to
run on a refined mesh, getting the hanging node constraints
right will be important.  Note that, even with
asymmetric_constraint_rows = false, the constrained dof
diagonals still exist in the matrix, with diagonal entries
that are there to ensure non-singular matrices for linear
solves but which would generate positive non-physical
eigenvalues for eigensolves.
dof_map.constrain_element_matrix(Ke, dof_indices, false);
dof_map.constrain_element_matrix(Me, dof_indices, false);


<br><br>Finally, simply add the element contribution to the
overall matrices A and B.
</div>

<div class ="fragment">
<pre>
              matrix_A.add_matrix (Ke, dof_indices);
              matrix_B.add_matrix (Me, dof_indices);
        
            } // end of element loop
        
        
        #endif // LIBMESH_HAVE_SLEPC
        
          /**
           * All done!
           */
          return;
        
        }
        
        void get_dirichlet_dofs(EquationSystems& es,
                                const std::string& system_name,
                                std::set&lt;unsigned int&gt;& dirichlet_dof_ids)
        {
        #ifdef LIBMESH_HAVE_SLEPC
        
          dirichlet_dof_ids.clear();
        
</pre>
</div>
<div class = "comment">
It is a good idea to make sure we are assembling
the proper system.
</div>

<div class ="fragment">
<pre>
          libmesh_assert_equal_to (system_name, "Eigensystem");
        
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
The dimension that we are running.
</div>

<div class ="fragment">
<pre>
          const unsigned int dim = mesh.mesh_dimension();
        
</pre>
</div>
<div class = "comment">
Get a reference to our system.
</div>

<div class ="fragment">
<pre>
          EigenSystem & eigen_system = es.get_system&lt;EigenSystem&gt; (system_name);
        
</pre>
</div>
<div class = "comment">
Get a constant reference to the Finite Element type
for the first (and only) variable in the system.
</div>

<div class ="fragment">
<pre>
          FEType fe_type = eigen_system.get_dof_map().variable_type(0);
          
          const DofMap& dof_map = eigen_system.get_dof_map();
        
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
matrix and right-hand-side contribution.  In case users
later modify this program to include refinement, we will
be safe and will only consider the active elements;
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
        
             {
</pre>
</div>
<div class = "comment">
All boundary dofs are Dirichlet dofs in this case
</div>

<div class ="fragment">
<pre>
                for (unsigned int s=0; s&lt;elem-&gt;n_sides(); s++)
                  if (elem-&gt;neighbor(s) == NULL)
                    {
                      std::vector&lt;unsigned int&gt; side_dofs;
                      FEInterface::dofs_on_side(elem, dim, fe_type,
                                                s, side_dofs);
        
                      for(unsigned int ii=0; ii&lt;side_dofs.size(); ii++)
                        dirichlet_dof_ids.insert(dof_indices[side_dofs[ii]]);
                    }
              }
        
            } // end of element loop
        
        #endif // LIBMESH_HAVE_SLEPC
        
          /**
           * All done!
           */
          return;
        
        }
        
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file eigenproblems_ex3.C without comments: </h1> 
<pre> 
  #include &lt;fstream&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/condensed_eigen_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_interface.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_matrices(EquationSystems&amp; es,
                         <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  <B><FONT COLOR="#228B22">void</FONT></B> get_dirichlet_dofs(EquationSystems&amp; es,
                          <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name,
                          <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;&amp; global_dirichlet_dofs_set);
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #ifndef LIBMESH_HAVE_EXODUS_API
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-exodus&quot;</FONT></B>);
  #endif
  
  #ifndef LIBMESH_HAVE_SLEPC
    <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR: This example requires libMesh to be\n&quot;</FONT></B>
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;compiled with SLEPc eigen solvers support!&quot;</FONT></B>
                &lt;&lt; std::endl;
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
  #ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-singleprecision&quot;</FONT></B>);
  #endif
  
  #ifdef LIBMESH_USE_COMPLEX_NUMBERS
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-complex&quot;</FONT></B>);
  #endif
  
    {
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Running &quot;</FONT></B> &lt;&lt; argv[0];
        
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">int</FONT></B> i=1; i&lt;argc; i++)
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B> &lt;&lt; argv[i];
        
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl &lt;&lt; std::endl;
    }
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
  
    GetPot command_line (argc, argv);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string mesh_name = <B><FONT COLOR="#BC8F8F">&quot;&quot;</FONT></B>;
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <B><FONT COLOR="#BC8F8F">&quot;-mesh_name&quot;</FONT></B>) )
      mesh_name = command_line.next(mesh_name);
      
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> plotting_index = 0;
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <B><FONT COLOR="#BC8F8F">&quot;-plotting_index&quot;</FONT></B>) )
      plotting_index = command_line.next(plotting_index);
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_evals = 0;
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <B><FONT COLOR="#BC8F8F">&quot;-n_evals&quot;</FONT></B>) )
      n_evals = command_line.next(n_evals);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream mesh_name_exodus;
    mesh_name_exodus &lt;&lt; mesh_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;_mesh.e&quot;</FONT></B>;
  
    Mesh mesh;
    mesh.read(mesh_name_exodus.str());
  
    mesh.print_info();
    
    EquationSystems equation_systems (mesh);
  
    CondensedEigenSystem &amp; eigen_system =
      equation_systems.add_system&lt;CondensedEigenSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Eigensystem&quot;</FONT></B>);
  
    eigen_system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>, SECOND);
  
    eigen_system.attach_assemble_function (assemble_matrices);
  
    equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;eigenpairs&quot;</FONT></B>)    = n_evals;
    equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;basis vectors&quot;</FONT></B>) = n_evals*3;
  
    equation_systems.parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;linear solver tolerance&quot;</FONT></B>) = pow(TOLERANCE, 5./3.);
    equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;
      (<B><FONT COLOR="#BC8F8F">&quot;linear solver maximum iterations&quot;</FONT></B>) = 1000;
  
    eigen_system.set_eigenproblem_type(GHEP);
  
    eigen_system.eigen_solver-&gt;set_position_of_spectrum(SMALLEST_MAGNITUDE);
  
    equation_systems.init();
  
    equation_systems.print_info();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dirichlet_dof_ids;
    get_dirichlet_dofs(equation_systems, <B><FONT COLOR="#BC8F8F">&quot;Eigensystem&quot;</FONT></B>, dirichlet_dof_ids);
    eigen_system.initialize_condensed_dofs(dirichlet_dof_ids);
  
    eigen_system.solve();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> nconv = eigen_system.get_n_converged();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Number of converged eigenpairs: &quot;</FONT></B> &lt;&lt; nconv
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B> &lt;&lt; std::endl;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (plotting_index &gt; n_evals)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;WARNING: Solver did not converge for the requested eigenvector!&quot;</FONT></B> &lt;&lt; std::endl;
      }
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream eigenvalue_output_name;
    eigenvalue_output_name &lt;&lt; mesh_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;_evals.txt&quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream evals_file(eigenvalue_output_name.str().c_str());
  
    <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;nconv; i++)
    {
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::pair&lt;Real,Real&gt; eval = eigen_system.get_eigenpair(i);
      
      libmesh_assert_less (eval.second, TOLERANCE);
      evals_file &lt;&lt; eval.first &lt;&lt; std::endl;
      
      <B><FONT COLOR="#A020F0">if</FONT></B>(i == plotting_index)
      {
  #ifdef LIBMESH_HAVE_EXODUS_API
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream eigenvector_output_name;
        eigenvector_output_name &lt;&lt; mesh_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;_evec.e&quot;</FONT></B>;
        ExodusII_IO (mesh).write_equation_systems (eigenvector_output_name.str(), equation_systems);
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>      }
    }
    
    evals_file.close();
  
  #endif <I><FONT COLOR="#B22222">// LIBMESH_HAVE_SLEPC
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_matrices(EquationSystems&amp; es,
                         <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Eigensystem&quot;</FONT></B>);
  
  #ifdef LIBMESH_HAVE_SLEPC
  
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    EigenSystem &amp; eigen_system = es.get_system&lt;EigenSystem&gt; (system_name);
  
    FEType fe_type = eigen_system.get_dof_map().variable_type(0);
  
    SparseMatrix&lt;Number&gt;&amp;  matrix_A = *eigen_system.matrix_A;
    SparseMatrix&lt;Number&gt;&amp;  matrix_B = *eigen_system.matrix_B;
  
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    QGauss qrule (dim, fe_type.default_quadrature_order());
  
    fe-&gt;attach_quadrature_rule (&amp;qrule);
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = eigen_system.get_dof_map();
  
    DenseMatrix&lt;Number&gt;   Me;
    DenseMatrix&lt;Number&gt;   Ke;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
  
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
   
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
  
        fe-&gt;reinit (elem);
  
        Ke.resize (dof_indices.size(), dof_indices.size());
        Me.resize (dof_indices.size(), dof_indices.size());
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
              {
                Me(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
                Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
              }
  
  
        matrix_A.add_matrix (Ke, dof_indices);
        matrix_B.add_matrix (Me, dof_indices);
  
      } <I><FONT COLOR="#B22222">// end of element loop
</FONT></I>  
  
  #endif <I><FONT COLOR="#B22222">// LIBMESH_HAVE_SLEPC
</FONT></I>  
    <I><FONT COLOR="#B22222">/**
     * All done!
     */</FONT></I>
    <B><FONT COLOR="#A020F0">return</FONT></B>;
  
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> get_dirichlet_dofs(EquationSystems&amp; es,
                          <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name,
                          <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;&amp; dirichlet_dof_ids)
  {
  #ifdef LIBMESH_HAVE_SLEPC
  
    dirichlet_dof_ids.clear();
  
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Eigensystem&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    EigenSystem &amp; eigen_system = es.get_system&lt;EigenSystem&gt; (system_name);
  
    FEType fe_type = eigen_system.get_dof_map().variable_type(0);
    
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = eigen_system.get_dof_map();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
  
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
   
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
  
       {
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> s=0; s&lt;elem-&gt;n_sides(); s++)
            <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(s) == NULL)
              {
                <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; side_dofs;
                <B><FONT COLOR="#5F9EA0">FEInterface</FONT></B>::dofs_on_side(elem, dim, fe_type,
                                          s, side_dofs);
  
                <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> ii=0; ii&lt;side_dofs.size(); ii++)
                  dirichlet_dof_ids.insert(dof_indices[side_dofs[ii]]);
              }
        }
  
      } <I><FONT COLOR="#B22222">// end of element loop
</FONT></I>  
  #endif <I><FONT COLOR="#B22222">// LIBMESH_HAVE_SLEPC
</FONT></I>  
    <I><FONT COLOR="#B22222">/**
     * All done!
     */</FONT></I>
    <B><FONT COLOR="#A020F0">return</FONT></B>;
  
  }
  
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example eigenproblems_ex3:
*  mpirun -np 12 example-devel -n_evals 5 -mesh_name drum1 -plotting_index 2 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/eigenproblems/eigenproblems_ex3/.libs/lt-example-devel -n_evals 5 -mesh_name drum1 -plotting_index 2 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1157
    n_local_nodes()=105
  n_elem()=536
    n_local_elem()=44
    n_active_elem()=536
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Eigensystem"
    Type "Eigen"
    Variables="p" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" 
    n_dofs()=1157
    n_local_dofs()=105
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=0
    n_matrices()=2
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 10.1106
      Average Off-Processor Bandwidth <= 1.36387
      Maximum  On-Processor Bandwidth <= 22
      Maximum Off-Processor Bandwidth <= 18
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

Number of converged eigenpairs: 5

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/eigenproblems/eigenproblems_ex3/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:05:57 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           1.782e+00      1.00000   1.782e+00
Objects:              1.020e+02      1.00000   1.020e+02
Flops:                6.824e+08      1.06796   6.653e+08  7.983e+09
Flops/sec:            3.829e+08      1.06796   3.733e+08  4.480e+09
MPI Messages:         2.791e+05      2.33252   1.994e+05  2.393e+06
MPI Message Lengths:  7.590e+07      1.21965   3.407e+02  8.153e+08
MPI Reductions:       3.208e+04      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.7820e+00 100.0%  7.9831e+09 100.0%  2.393e+06 100.0%  3.407e+02      100.0%  3.207e+04 100.0% 

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

STSetUp                1 1.0 6.0861e-03 1.0 9.02e+05 1.0 7.1e+02 2.3e+03 2.9e+01  0  0  0  0  0   0  0  0  0  0  1779
STApply             7974 1.0 9.4450e-01 1.0 5.75e+08 1.0 1.4e+06 5.2e+02 0.0e+00 52 86 58 89  0  52 86 58 89  0  7275
EPSSetUp               1 1.0 6.9010e-03 1.0 9.02e+05 1.0 7.1e+02 2.3e+03 6.5e+01  0  0  0  0  0   0  0  0  0  0  1569
EPSSolve               1 1.0 1.5585e+00 1.0 6.81e+08 1.1 2.4e+06 3.4e+02 3.2e+04 87100100100 99  87100100100 99  5112
IPOrthogonalize     7969 1.0 5.7067e-01 1.1 1.04e+08 1.5 1.0e+06 8.9e+01 3.2e+04 31 13 42 11 99  31 13 42 11 99  1885
IPInnerProduct     79648 1.0 5.2473e-01 1.1 7.77e+07 1.6 1.0e+06 8.9e+01 3.2e+04 29 10 42 11 99  29 10 42 11 99  1518
IPApplyMatrix      23899 1.0 2.3976e-01 1.6 4.72e+07 1.7 1.0e+06 8.9e+01 0.0e+00 11  6 42 11  0  11  6 42 11  0  1969
DSSolve              607 1.0 3.6283e-02 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  2  0  0  0  0   2  0  0  0  0     0
DSVectors            612 1.0 3.5334e-03 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
DSOther             1214 1.0 1.2735e-02 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
UpdateVectors        607 1.0 7.5262e-03 1.6 1.71e+06 1.4 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2413
VecMAXPBY          15925 1.0 2.4782e-02 1.2 2.62e+07 1.4 0.0e+00 0.0e+00 0.0e+00  1  3  0  0  0   1  3  0  0  0 11264
VecScale            7367 1.0 5.2018e-03 1.4 6.85e+05 1.4 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1401
VecCopy               15 1.0 3.0041e-05 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet             15964 1.0 1.1081e-02 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
VecAssemblyBegin      15 1.0 4.1223e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 4.5e+01  0  0  0  0  0   0  0  0  0  0     0
VecAssemblyEnd        15 1.0 4.6015e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin    47826 1.0 1.7823e-01 1.6 0.00e+00 0.0 2.4e+06 3.4e+02 0.0e+00  8  0100100  0   8  0100100  0     0
VecScatterEnd      47826 1.0 1.7680e-01 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  8  0  0  0  0   8  0  0  0  0     0
VecReduceArith     39824 1.0 4.3311e-02 1.2 3.05e+07 1.4 0.0e+00 0.0e+00 0.0e+00  2  4  0  0  0   2  4  0  0  0  7491
VecReduceComm      31856 1.0 2.9741e-01 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 3.2e+04 14  0  0  0 99  14  0  0  0 99     0
MatMult            31873 1.0 3.1583e-01 1.7 6.29e+07 1.7 1.3e+06 8.9e+01 0.0e+00 14  8 56 15  0  14  8 56 15  0  1994
MatSolve            7974 1.0 5.6140e-01 1.0 5.59e+08 1.0 0.0e+00 0.0e+00 0.0e+00 31 84  0  0  0  31 84  0  0  0 11959
MatLUFactorSym         1 1.0 1.1790e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatLUFactorNum         1 1.0 1.6680e-03 1.1 9.02e+05 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  6492
MatAssemblyBegin      11 1.0 7.7701e-04 2.9 0.00e+00 0.0 1.2e+02 4.8e+02 1.6e+01  0  0  0  0  0   0  0  0  0  0     0
MatAssemblyEnd        11 1.0 1.1303e-03 1.0 0.00e+00 0.0 3.4e+02 2.5e+01 3.2e+01  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ            1 1.0 1.3304e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 1.0359e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         4 1.0 2.0981e-05 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRedundant        1 1.0 7.7486e-04 1.0 0.00e+00 0.0 4.0e+02 3.7e+03 4.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSetUp               2 1.0 2.1458e-06 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve            7974 1.0 8.0519e-01 1.0 5.59e+08 1.0 1.1e+06 6.6e+02 0.0e+00 45 84 44 85  0  45 84 44 85  0  8338
PCSetUp                1 1.0 5.9872e-03 1.0 9.02e+05 1.0 7.1e+02 2.3e+03 2.7e+01  0  0  0  0  0   0  0  0  0  0  1809
PCApply             7974 1.0 7.8652e-01 1.0 5.59e+08 1.0 1.1e+06 6.6e+02 0.0e+00 44 84 44 85  0  44 84 44 85  0  8536
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

           Container     3              3         1644     0
  Spectral Transform     1              1          736     0
 Eigenproblem Solver     1              1         1680     0
       Inner product     1              1          624     0
       Direct solver     1              1         9408     0
              Vector    42             42        89296     0
      Vector Scatter     7              7         7252     0
           Index Set    23             23        34184     0
   IS L to G Mapping     1              1          564     0
              Matrix    16             16       834072     0
         PetscRandom     1              1          608     0
       Krylov Solver     2              2         2144     0
      Preconditioner     2              2         1800     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 1.90735e-07
Average time for MPI_Barrier(): 5.62668e-06
Average time for zero size MPI_Send(): 1.31726e-05
#PETSc Option Table entries:
-ksp_right_pc
-log_summary
-mesh_name drum1
-n_evals 5
-pc_type bjacobi
-plotting_index 2
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
| Time:           Thu Jan 31 22:05:57 2013                                                                             |
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
| libMesh Performance: Alive time=1.86087, Active time=1.76021                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| CondensedEigenSystem                                                                                           |
|   get_eigenpair()                  5         0.0021      0.000411    0.0022      0.000447    0.12     0.13     |
|   solve()                          1         0.0039      0.003921    1.5790      1.579035    0.22     89.71    |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0038      0.003825    0.0059      0.005913    0.22     0.34     |
|   build_sparsity()                 1         0.0029      0.002861    0.0082      0.008204    0.16     0.47     |
|   create_dof_constraints()         1         0.0011      0.001105    0.0011      0.001105    0.06     0.06     |
|   distribute_dofs()                1         0.0149      0.014909    0.0466      0.046610    0.85     2.65     |
|   dof_indices()                    202       0.0161      0.000080    0.0161      0.000080    0.91     0.91     |
|   prepare_send_list()              1         0.0001      0.000056    0.0001      0.000056    0.00     0.00     |
|   reinit()                         1         0.0288      0.028767    0.0288      0.028767    1.63     1.63     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0009      0.000909    0.0051      0.005104    0.05     0.29     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0055      0.005498    0.0055      0.005498    0.31     0.31     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        44        0.0011      0.000024    0.0011      0.000024    0.06     0.06     |
|   init_shape_functions()           1         0.0001      0.000136    0.0001      0.000136    0.01     0.01     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             44        0.0006      0.000014    0.0006      0.000014    0.04     0.04     |
|   init_reference_to_physical_map() 1         0.0001      0.000089    0.0001      0.000089    0.01     0.01     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0087      0.008733    0.0092      0.009221    0.50     0.52     |
|   read()                           1         0.0143      0.014344    0.0143      0.014344    0.81     0.81     |
|   renumber_nodes_and_elem()        2         0.0006      0.000323    0.0006      0.000323    0.04     0.04     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   broadcast()                      1         0.0085      0.008545    0.0158      0.015760    0.49     0.90     |
|   compute_hilbert_indices()        2         0.0094      0.004681    0.0094      0.004681    0.53     0.53     |
|   find_global_indices()            2         0.0039      0.001969    0.0193      0.009651    0.22     1.10     |
|   parallel_sort()                  2         0.0029      0.001458    0.0046      0.002277    0.17     0.26     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0002      0.000152    0.0109      0.010936    0.01     0.62     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0333      0.033343    0.0416      0.041606    1.89     2.36     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      9         0.0014      0.000154    0.0015      0.000164    0.08     0.08     |
|   broadcast()                      6         0.0004      0.000067    0.0004      0.000067    0.02     0.02     |
|   max(bool)                        2         0.0000      0.000007    0.0000      0.000007    0.00     0.00     |
|   max(scalar)                      120       0.0009      0.000008    0.0009      0.000008    0.05     0.05     |
|   max(vector)                      27        0.0004      0.000013    0.0009      0.000032    0.02     0.05     |
|   min(bool)                        138       0.0008      0.000006    0.0008      0.000006    0.05     0.05     |
|   min(scalar)                      113       0.0137      0.000121    0.0137      0.000121    0.78     0.78     |
|   min(vector)                      27        0.0005      0.000018    0.0012      0.000045    0.03     0.07     |
|   probe()                          132       0.0016      0.000012    0.0016      0.000012    0.09     0.09     |
|   receive()                        132       0.0008      0.000006    0.0025      0.000019    0.05     0.14     |
|   send()                           132       0.0004      0.000003    0.0004      0.000003    0.02     0.02     |
|   send_receive()                   136       0.0011      0.000008    0.0044      0.000032    0.06     0.25     |
|   sum()                            25        0.0018      0.000072    0.0023      0.000093    0.10     0.13     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           132       0.0003      0.000002    0.0003      0.000002    0.02     0.02     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0015      0.001528    0.0026      0.002577    0.09     0.15     |
|   set_parent_processor_ids()       1         0.0010      0.000986    0.0010      0.000986    0.06     0.06     |
|                                                                                                                |
| SlepcEigenSolver                                                                                               |
|   solve_generalized()              1         1.5672      1.567212    1.5672      1.567212    89.04    89.04    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0024      0.002400    0.0079      0.007901    0.14     0.45     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1454      1.7602                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example eigenproblems_ex3:
*  mpirun -np 12 example-devel -n_evals 5 -mesh_name drum1 -plotting_index 2 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
***************************************************************
* Running Example eigenproblems_ex3:
*  mpirun -np 12 example-devel -n_evals 5 -mesh_name drum2 -plotting_index 2 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/eigenproblems/eigenproblems_ex3/.libs/lt-example-devel -n_evals 5 -mesh_name drum2 -plotting_index 2 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1189
    n_local_nodes()=113
  n_elem()=552
    n_local_elem()=46
    n_active_elem()=552
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Eigensystem"
    Type "Eigen"
    Variables="p" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" 
    n_dofs()=1189
    n_local_dofs()=113
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=0
    n_matrices()=2
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 10.1396
      Average Off-Processor Bandwidth <= 1.37595
      Maximum  On-Processor Bandwidth <= 23
      Maximum Off-Processor Bandwidth <= 15
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

Number of converged eigenpairs: 5

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/eigenproblems/eigenproblems_ex3/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:05:59 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           1.798e+00      1.00000   1.798e+00
Objects:              1.020e+02      1.00000   1.020e+02
Flops:                7.164e+08      1.06211   6.946e+08  8.336e+09
Flops/sec:            3.985e+08      1.06211   3.864e+08  4.637e+09
MPI Messages:         2.783e+05      1.84170   1.935e+05  2.322e+06
MPI Message Lengths:  7.957e+07      1.25776   3.614e+02  8.392e+08
MPI Reductions:       3.198e+04      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.7977e+00 100.0%  8.3355e+09 100.0%  2.322e+06 100.0%  3.614e+02      100.0%  3.198e+04 100.0% 

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

STSetUp                1 1.0 6.3040e-03 1.0 9.66e+05 1.0 7.1e+02 2.3e+03 2.9e+01  0  0  0  0  0   0  0  0  0  0  1839
STApply             7951 1.0 9.6566e-01 1.0 6.02e+08 1.0 1.4e+06 5.4e+02 0.0e+00 53 86 59 89  0  53 86 59 89  0  7443
EPSSetUp               1 1.0 7.0720e-03 1.0 9.66e+05 1.0 7.1e+02 2.3e+03 6.5e+01  0  0  0  0  0   0  0  0  0  0  1639
EPSSolve               1 1.0 1.5781e+00 1.0 7.15e+08 1.1 2.3e+06 3.6e+02 3.2e+04 88100100100 99  88100100100 99  5272
IPOrthogonalize     7946 1.0 5.5866e-01 1.0 1.11e+08 1.5 9.5e+05 9.7e+01 3.2e+04 31 13 41 11 99  31 13 41 11 99  1987
IPInnerProduct     79418 1.0 5.1794e-01 1.1 8.28e+07 1.5 9.5e+05 9.7e+01 3.2e+04 28 10 41 11 99  28 10 41 11 99  1587
IPApplyMatrix      23830 1.0 2.4141e-01 1.5 5.03e+07 1.6 9.5e+05 9.7e+01 0.0e+00 10  6 41 11  0  10  6 41 11  0  2017
DSSolve              608 1.0 3.6442e-02 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  2  0  0  0  0   2  0  0  0  0     0
DSVectors            613 1.0 3.7684e-03 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
DSOther             1216 1.0 1.1626e-02 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
UpdateVectors        608 1.0 6.5277e-03 1.9 1.82e+06 1.4 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2876
VecMAXPBY          15879 1.0 2.5368e-02 1.4 2.80e+07 1.4 0.0e+00 0.0e+00 0.0e+00  1  3  0  0  0   1  3  0  0  0 11363
VecScale            7343 1.0 5.4152e-03 1.6 7.27e+05 1.4 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1384
VecCopy               15 1.0 2.5988e-05 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet             15918 1.0 1.1646e-02 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
VecAssemblyBegin      15 1.0 3.8218e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 4.5e+01  0  0  0  0  0   0  0  0  0  0     0
VecAssemblyEnd        15 1.0 4.4107e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin    47688 1.0 1.7448e-01 1.3 0.00e+00 0.0 2.3e+06 3.6e+02 0.0e+00  8  0100100  0   8  0100100  0     0
VecScatterEnd      47688 1.0 1.7442e-01 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  8  0  0  0  0   8  0  0  0  0     0
VecReduceArith     39709 1.0 4.3162e-02 1.2 3.25e+07 1.4 0.0e+00 0.0e+00 0.0e+00  2  4  0  0  0   2  4  0  0  0  7760
VecReduceComm      31764 1.0 2.7426e-01 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 3.2e+04 13  0  0  0 99  13  0  0  0 99     0
MatMult            31781 1.0 3.1578e-01 1.5 6.70e+07 1.6 1.3e+06 9.7e+01 0.0e+00 14  8 55 15  0  14  8 55 15  0  2056
MatSolve            7951 1.0 5.8667e-01 1.0 5.85e+08 1.0 0.0e+00 0.0e+00 0.0e+00 32 84  0  0  0  32 84  0  0  0 11974
MatLUFactorSym         1 1.0 1.2400e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatLUFactorNum         1 1.0 1.7231e-03 1.0 9.66e+05 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  6728
MatAssemblyBegin      11 1.0 6.4802e-04 2.5 0.00e+00 0.0 1.2e+02 5.1e+02 1.6e+01  0  0  0  0  0   0  0  0  0  0     0
MatAssemblyEnd        11 1.0 1.0960e-03 1.0 0.00e+00 0.0 3.2e+02 2.7e+01 3.2e+01  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ            1 1.0 1.2517e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 1.1110e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         4 1.0 2.3127e-05 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRedundant        1 1.0 7.9393e-04 1.0 0.00e+00 0.0 4.0e+02 3.8e+03 4.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSetUp               2 1.0 9.5367e-07 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve            7951 1.0 8.3778e-01 1.0 5.85e+08 1.0 1.0e+06 6.8e+02 0.0e+00 46 84 45 85  0  46 84 45 85  0  8385
PCSetUp                1 1.0 6.1941e-03 1.0 9.66e+05 1.0 7.1e+02 2.3e+03 2.7e+01  0  0  0  0  0   0  0  0  0  0  1872
PCApply             7951 1.0 8.1767e-01 1.0 5.85e+08 1.0 1.0e+06 6.8e+02 0.0e+00 45 84 45 85  0  45 84 45 85  0  8592
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

           Container     3              3         1644     0
  Spectral Transform     1              1          736     0
 Eigenproblem Solver     1              1         1680     0
       Inner product     1              1          624     0
       Direct solver     1              1         9408     0
              Vector    42             42        89360     0
      Vector Scatter     7              7         7252     0
           Index Set    23             23        34048     0
   IS L to G Mapping     1              1          564     0
              Matrix    16             16       861480     0
         PetscRandom     1              1          608     0
       Krylov Solver     2              2         2144     0
      Preconditioner     2              2         1800     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 8.2016e-06
Average time for zero size MPI_Send(): 1.33316e-05
#PETSc Option Table entries:
-ksp_right_pc
-log_summary
-mesh_name drum2
-n_evals 5
-pc_type bjacobi
-plotting_index 2
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
| Time:           Thu Jan 31 22:06:00 2013                                                                             |
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
| libMesh Performance: Alive time=1.82159, Active time=1.77593                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| CondensedEigenSystem                                                                                           |
|   get_eigenpair()                  5         0.0020      0.000404    0.0022      0.000449    0.11     0.13     |
|   solve()                          1         0.0035      0.003471    1.5985      1.598503    0.20     90.01    |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0039      0.003921    0.0056      0.005604    0.22     0.32     |
|   build_sparsity()                 1         0.0029      0.002904    0.0082      0.008186    0.16     0.46     |
|   create_dof_constraints()         1         0.0010      0.000977    0.0010      0.000977    0.06     0.06     |
|   distribute_dofs()                1         0.0157      0.015708    0.0472      0.047166    0.88     2.66     |
|   dof_indices()                    205       0.0162      0.000079    0.0162      0.000079    0.91     0.91     |
|   prepare_send_list()              1         0.0000      0.000047    0.0000      0.000047    0.00     0.00     |
|   reinit()                         1         0.0297      0.029743    0.0297      0.029743    1.67     1.67     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0009      0.000941    0.0052      0.005189    0.05     0.29     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0056      0.005605    0.0056      0.005605    0.32     0.32     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        46        0.0011      0.000023    0.0011      0.000023    0.06     0.06     |
|   init_shape_functions()           1         0.0001      0.000141    0.0001      0.000141    0.01     0.01     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             46        0.0007      0.000015    0.0007      0.000015    0.04     0.04     |
|   init_reference_to_physical_map() 1         0.0001      0.000090    0.0001      0.000090    0.01     0.01     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0089      0.008911    0.0090      0.009044    0.50     0.51     |
|   read()                           1         0.0069      0.006939    0.0069      0.006939    0.39     0.39     |
|   renumber_nodes_and_elem()        2         0.0007      0.000329    0.0007      0.000329    0.04     0.04     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   broadcast()                      1         0.0088      0.008798    0.0159      0.015931    0.50     0.90     |
|   compute_hilbert_indices()        2         0.0095      0.004757    0.0095      0.004757    0.54     0.54     |
|   find_global_indices()            2         0.0040      0.002017    0.0182      0.009121    0.23     1.03     |
|   parallel_sort()                  2         0.0029      0.001425    0.0033      0.001656    0.16     0.19     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0001      0.000132    0.0111      0.011140    0.01     0.63     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0376      0.037642    0.0459      0.045935    2.12     2.59     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      9         0.0005      0.000056    0.0005      0.000061    0.03     0.03     |
|   broadcast()                      6         0.0004      0.000068    0.0004      0.000068    0.02     0.02     |
|   max(bool)                        2         0.0000      0.000007    0.0000      0.000007    0.00     0.00     |
|   max(scalar)                      120       0.0008      0.000007    0.0008      0.000007    0.05     0.05     |
|   max(vector)                      27        0.0004      0.000013    0.0009      0.000033    0.02     0.05     |
|   min(bool)                        138       0.0009      0.000006    0.0009      0.000006    0.05     0.05     |
|   min(scalar)                      113       0.0129      0.000114    0.0129      0.000114    0.73     0.73     |
|   min(vector)                      27        0.0005      0.000017    0.0010      0.000038    0.03     0.06     |
|   probe()                          132       0.0013      0.000010    0.0013      0.000010    0.07     0.07     |
|   receive()                        132       0.0008      0.000006    0.0021      0.000016    0.05     0.12     |
|   send()                           132       0.0004      0.000003    0.0004      0.000003    0.02     0.02     |
|   send_receive()                   136       0.0011      0.000008    0.0040      0.000029    0.06     0.23     |
|   sum()                            25        0.0007      0.000028    0.0012      0.000047    0.04     0.07     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           132       0.0003      0.000002    0.0003      0.000002    0.02     0.02     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0016      0.001555    0.0022      0.002193    0.09     0.12     |
|   set_parent_processor_ids()       1         0.0011      0.001074    0.0011      0.001074    0.06     0.06     |
|                                                                                                                |
| SlepcEigenSolver                                                                                               |
|   solve_generalized()              1         1.5870      1.586954    1.5870      1.586954    89.36    89.36    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0024      0.002430    0.0081      0.008078    0.14     0.45     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1461      1.7759                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example eigenproblems_ex3:
*  mpirun -np 12 example-devel -n_evals 5 -mesh_name drum2 -plotting_index 2 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
