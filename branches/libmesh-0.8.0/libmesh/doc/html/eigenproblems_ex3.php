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
        #include "libmesh.h"
        #include "mesh.h"
        #include "mesh_generation.h"
        #include "exodusII_io.h"
        #include "condensed_eigen_system.h"
        #include "equation_systems.h"
        #include "fe.h"
        #include "quadrature_gauss.h"
        #include "dense_matrix.h"
        #include "sparse_matrix.h"
        #include "numeric_vector.h"
        #include "dof_map.h"
        #include "fe_interface.h"
        #include "getpot.h"
        
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
            libmesh_assert(eval.second &lt; TOLERANCE);
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
          libmesh_assert (system_name == "Eigensystem");
        
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
          libmesh_assert (system_name == "Eigensystem");
        
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
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  #include &lt;fstream&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;condensed_eigen_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe_interface.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  
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
      
      libmesh_assert(eval.second &lt; TOLERANCE);
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
    
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Eigensystem&quot;</FONT></B>);
  
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
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
  
  
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
  
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Eigensystem&quot;</FONT></B>);
  
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
Linking eigenproblems_ex3-opt...
***************************************************************
* Running Example  mpirun -np 6 ./eigenproblems_ex3-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running ./eigenproblems_ex3-opt -n_evals 5 -mesh_name drum1 -plotting_index 2 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1157
    n_local_nodes()=213
  n_elem()=536
    n_local_elem()=90
    n_active_elem()=536
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
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
    n_local_dofs()=213
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=0
    n_matrices()=2
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 10.3615
      Average Off-Processor Bandwidth <= 0.394366
      Maximum  On-Processor Bandwidth <= 22
      Maximum Off-Processor Bandwidth <= 10
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

Number of converged eigenpairs: 5

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./eigenproblems_ex3-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:21:09 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           3.438e+00      1.00022   3.437e+00
Objects:              9.700e+01      1.00000   9.700e+01
Flops:                7.211e+08      1.03477   7.113e+08  4.268e+09
Flops/sec:            2.098e+08      1.03454   2.069e+08  1.242e+09
MPI Messages:         1.447e+05      1.33315   1.206e+05  7.237e+05
MPI Message Lengths:  6.544e+07      1.19131   4.987e+02  3.609e+08
MPI Reductions:       2.913e+04      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 3.4371e+00 100.0%  4.2675e+09 100.0%  7.237e+05 100.0%  4.987e+02      100.0%  2.906e+04  99.8% 

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

STSetUp                1 1.0 7.8549e-03 1.5 8.09e+05 1.0 1.7e+02 4.2e+03 1.7e+01  0  0  0  0  0   0  0  0  0  0   618
STApply             7239 1.0 9.0502e-01 1.1 5.21e+08 1.0 3.2e+05 9.5e+02 0.0e+00 25 73 44 83  0  25 73 44 83  0  3444
EPSSetUp               1 1.0 9.9211e-03 1.3 8.09e+05 1.0 1.7e+02 4.2e+03 3.2e+01  0  0  0  0  0   0  0  0  0  0   489
EPSSolve               1 1.0 3.3619e+00 1.0 7.20e+08 1.0 7.2e+05 5.0e+02 2.9e+04 98100100100 99  98100100100 99  1267
EPSDense            1108 1.0 3.5072e-02 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
IPOrthogonalize     7234 1.0 2.4776e+00 1.0 1.95e+08 1.1 4.0e+05 1.5e+02 2.9e+04 71 26 56 16 99  71 26 56 16 99   453
IPInnerProduct     57825 1.0 2.4399e+00 1.0 1.52e+08 1.1 4.0e+05 1.5e+02 2.9e+04 70 20 56 16 99  70 20 56 16 99   356
IPApplyMatrix      28908 1.0 1.0896e+00 2.8 1.01e+08 1.1 4.0e+05 1.5e+02 0.0e+00 19 13 56 16  0  19 13 56 16  0   524
UpdateVectors        554 1.0 2.7199e-03 1.2 2.85e+06 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  6098
VecMAXPBY          14451 1.0 2.3681e-02 1.3 4.38e+07 1.1 0.0e+00 0.0e+00 0.0e+00  1  6  0  0  0   1  6  0  0  0 10749
VecDot             14457 1.0 7.2897e-01 2.2 4.90e+06 1.1 0.0e+00 0.0e+00 1.4e+04 15  1  0  0 50  15  1  0  0 50    39
VecScale            6685 1.0 4.4646e-03 1.2 1.14e+06 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1481
VecCopy               15 1.0 9.2983e-06 2.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet              7251 1.0 4.8170e-03 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAssemblyBegin      15 1.0 1.5845e-03 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 4.5e+01  0  0  0  0  0   0  0  0  0  0     0
VecAssemblyEnd        15 1.0 3.4332e-05 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin    50630 1.0 1.0066e-01 1.4 0.00e+00 0.0 7.2e+05 5.0e+02 0.0e+00  2  0100100  0   2  0100100  0     0
VecScatterEnd      50630 1.0 1.1677e+00 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 24  0  0  0  0  24  0  0  0  0     0
VecReduceArith     21684 1.0 2.5478e-02 1.3 4.61e+07 1.1 0.0e+00 0.0e+00 0.0e+00  1  6  0  0  0   1  6  0  0  0 10520
VecReduceComm      14451 1.0 1.4937e+00 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 1.4e+04 35  0  0  0 50  35  0  0  0 50     0
MatMult            36147 1.0 1.2679e+00 3.0 1.26e+08 1.1 5.1e+05 1.5e+02 0.0e+00 22 17 70 20  0  22 17 70 20  0   563
MatSolve            7239 1.0 4.4975e-01 1.2 4.96e+08 1.0 0.0e+00 0.0e+00 0.0e+00 11 70  0  0  0  11 70  0  0  0  6613
MatLUFactorSym         1 1.0 3.3629e-03 4.3 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatLUFactorNum         1 1.0 1.4100e-03 1.4 8.09e+05 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  3441
MatAssemblyBegin      11 1.0 2.2397e-03 2.4 0.00e+00 0.0 4.2e+01 8.3e+02 1.6e+01  0  0  0  0  0   0  0  0  0  0     0
MatAssemblyEnd        11 1.0 1.4291e-03 1.2 0.00e+00 0.0 1.1e+02 4.1e+01 3.2e+01  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ            1 1.0 5.4693e-04 7.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetSubMatrice       2 1.0 4.5300e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+01  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 1.3349e-03 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         4 1.0 2.0981e-05 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRedundant        1 1.0 1.6651e-03 2.9 0.00e+00 0.0 9.0e+01 7.4e+03 2.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSetup               2 1.0 9.5367e-07 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve            7239 1.0 7.9428e-01 1.3 4.96e+08 1.0 2.2e+05 1.3e+03 0.0e+00 20 70 30 79  0  20 70 30 79  0  3745
PCSetUp                1 1.0 7.8001e-03 1.5 8.09e+05 1.0 1.7e+02 4.2e+03 1.7e+01  0  0  0  0  0   0  0  0  0  0   622
PCApply             7239 1.0 7.7690e-01 1.3 4.96e+08 1.0 2.2e+05 1.3e+03 0.0e+00 20 70 30 79  0  20 70 30 79  0  3828
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

           Container     1              1          388     0
  Spectral Transform     1              1          568     0
 Eigenproblem Solver     1              1          924     0
       Inner product     1              1          428     0
                 Vec    41             41        89944     0
         Vec Scatter     7              7         6076     0
           Index Set    23             23        40344     0
   IS L to G Mapping     1              1         1396     0
              Matrix    16             16       893544     0
       Krylov Solver     2              2         1664     0
      Preconditioner     2              2         1416     0
         PetscRandom     1              1          448     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.03881e-05
Average time for zero size MPI_Send(): 6.75122e-05
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
| Time:           Fri Aug 24 15:21:09 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=3.64366, Active time=3.3998                                                    |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| CondensedEigenSystem                                                                                           |
|   get_eigenpair()                  5         0.0038      0.000759    0.0044      0.000886    0.11     0.13     |
|   solve()                          1         0.0064      0.006383    3.3771      3.377132    0.19     99.33    |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0001      0.000071    0.0001      0.000081    0.00     0.00     |
|   build_sparsity()                 1         0.0007      0.000726    0.0009      0.000938    0.02     0.03     |
|   create_dof_constraints()         1         0.0001      0.000066    0.0001      0.000066    0.00     0.00     |
|   distribute_dofs()                1         0.0003      0.000258    0.0024      0.002365    0.01     0.07     |
|   dof_indices()                    832       0.0003      0.000000    0.0003      0.000000    0.01     0.01     |
|   prepare_send_list()              1         0.0000      0.000004    0.0000      0.000004    0.00     0.00     |
|   reinit()                         1         0.0005      0.000493    0.0005      0.000493    0.01     0.01     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0002      0.000214    0.0008      0.000750    0.01     0.02     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0010      0.000966    0.0010      0.000966    0.03     0.03     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        90        0.0001      0.000002    0.0001      0.000002    0.00     0.00     |
|   init_shape_functions()           1         0.0000      0.000014    0.0000      0.000014    0.00     0.00     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             90        0.0001      0.000001    0.0001      0.000001    0.00     0.00     |
|   init_reference_to_physical_map() 1         0.0000      0.000018    0.0000      0.000018    0.00     0.00     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0003      0.000327    0.0011      0.001137    0.01     0.03     |
|   read()                           1         0.0010      0.000954    0.0010      0.000954    0.03     0.03     |
|   renumber_nodes_and_elem()        2         0.0001      0.000035    0.0001      0.000035    0.00     0.00     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   broadcast()                      1         0.0003      0.000277    0.0008      0.000801    0.01     0.02     |
|   compute_hilbert_indices()        2         0.0020      0.000987    0.0020      0.000987    0.06     0.06     |
|   find_global_indices()            2         0.0002      0.000118    0.0066      0.003294    0.01     0.19     |
|   parallel_sort()                  2         0.0012      0.000599    0.0033      0.001661    0.04     0.10     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0000      0.000027    0.0017      0.001743    0.00     0.05     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0012      0.001170    0.0041      0.004115    0.03     0.12     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      8         0.0013      0.000163    0.0013      0.000163    0.04     0.04     |
|   broadcast()                      7         0.0001      0.000014    0.0001      0.000014    0.00     0.00     |
|   gather()                         1         0.0001      0.000052    0.0001      0.000052    0.00     0.00     |
|   max(scalar)                      3         0.0013      0.000425    0.0013      0.000425    0.04     0.04     |
|   max(vector)                      2         0.0001      0.000038    0.0001      0.000038    0.00     0.00     |
|   min(vector)                      2         0.0002      0.000118    0.0002      0.000118    0.01     0.01     |
|   probe()                          50        0.0028      0.000055    0.0028      0.000055    0.08     0.08     |
|   receive()                        50        0.0001      0.000002    0.0029      0.000057    0.00     0.08     |
|   send()                           50        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   send_receive()                   54        0.0001      0.000002    0.0030      0.000056    0.00     0.09     |
|   sum()                            15        0.0034      0.000227    0.0034      0.000227    0.10     0.10     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           50        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0001      0.000105    0.0017      0.001750    0.00     0.05     |
|   set_parent_processor_ids()       1         0.0000      0.000035    0.0000      0.000035    0.00     0.00     |
|                                                                                                                |
| SlepcEigenSolver                                                                                               |
|   solve_generalized()              1         3.3699      3.369850    3.3699      3.369850    99.12    99.12    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0006      0.000589    0.0009      0.000899    0.02     0.03     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1337      3.3998                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

Running ./eigenproblems_ex3-opt -n_evals 5 -mesh_name drum2 -plotting_index 2 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1189
    n_local_nodes()=211
  n_elem()=552
    n_local_elem()=92
    n_active_elem()=552
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
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
    n_local_dofs()=211
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=0
    n_matrices()=2
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 10.5261
      Average Off-Processor Bandwidth <= 0.966825
      Maximum  On-Processor Bandwidth <= 22
      Maximum Off-Processor Bandwidth <= 16
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

Number of converged eigenpairs: 5

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./eigenproblems_ex3-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:21:14 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           4.799e+00      1.00017   4.799e+00
Objects:              9.700e+01      1.00000   9.700e+01
Flops:                1.042e+09      1.08629   9.871e+08  5.923e+09
Flops/sec:            2.172e+08      1.08643   2.057e+08  1.234e+09
MPI Messages:         2.234e+05      1.66642   1.787e+05  1.072e+06
MPI Message Lengths:  8.738e+07      1.27453   4.271e+02  4.581e+08
MPI Reductions:       3.594e+04      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 4.7986e+00 100.0%  5.9228e+09 100.0%  1.072e+06 100.0%  4.271e+02      100.0%  3.587e+04  99.8% 

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

STSetUp                1 1.0 6.0461e-03 1.1 1.20e+06 1.0 1.7e+02 4.4e+03 1.7e+01  0  0  0  0  0   0  0  0  0  0  1186
STApply             8939 1.0 1.7078e+00 1.6 7.49e+08 1.0 4.3e+05 8.9e+02 0.0e+00 28 75 40 84  0  28 75 40 84  0  2605
EPSSetUp               1 1.0 7.8928e-03 1.1 1.20e+06 1.0 1.7e+02 4.4e+03 3.2e+01  0  0  0  0  0   0  0  0  0  0   909
EPSSolve               1 1.0 4.6835e+00 1.0 1.04e+09 1.1 1.1e+06 4.3e+02 3.6e+04 98100100100 99  98100100100100  1263
EPSDense            1382 1.0 3.9797e-02 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
IPOrthogonalize     8934 1.0 3.5967e+00 1.2 2.87e+08 1.3 6.4e+05 1.1e+02 3.6e+04 68 24 60 16 99  68 24 60 16100   400
IPInnerProduct     71440 1.0 3.5505e+00 1.2 2.24e+08 1.4 6.4e+05 1.1e+02 3.6e+04 67 19 60 16 99  67 19 60 16100   313
IPApplyMatrix      35718 1.0 1.6428e+00 3.3 1.51e+08 1.4 6.4e+05 1.1e+02 0.0e+00 24 12 60 16  0  24 12 60 16  0   444
UpdateVectors        691 1.0 3.6130e-03 1.3 4.09e+06 1.2 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  5901
VecMAXPBY          17856 1.0 2.9634e-02 1.2 6.27e+07 1.2 0.0e+00 0.0e+00 0.0e+00  1  6  0  0  0   1  6  0  0  0 11030
VecDot             17862 1.0 1.1174e+00 1.9 6.98e+06 1.2 0.0e+00 0.0e+00 1.8e+04 18  1  0  0 50  18  1  0  0 50    33
VecScale            8248 1.0 5.3411e-03 1.2 1.62e+06 1.2 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1577
VecCopy               15 1.0 8.8215e-06 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet              8951 1.0 5.4462e-03 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAssemblyBegin      15 1.0 1.8647e-03 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 4.5e+01  0  0  0  0  0   0  0  0  0  0     0
VecAssemblyEnd        15 1.0 3.6478e-05 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin    62540 1.0 1.3357e-01 1.4 0.00e+00 0.0 1.1e+06 4.3e+02 0.0e+00  2  0100100  0   2  0100100  0     0
VecScatterEnd      62540 1.0 2.0532e+00 2.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 32  0  0  0  0  32  0  0  0  0     0
VecReduceArith     26789 1.0 2.8459e-02 1.1 6.61e+07 1.2 0.0e+00 0.0e+00 0.0e+00  1  6  0  0  0   1  6  0  0  0 12091
VecReduceComm      17856 1.0 1.6820e+00 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 1.8e+04 24  0  0  0 50  24  0  0  0 50     0
MatMult            44657 1.0 2.0915e+00 3.7 1.89e+08 1.4 8.0e+05 1.1e+02 0.0e+00 28 15 75 20  0  28 15 75 20  0   436
MatSolve            8939 1.0 5.7215e-01 1.1 7.11e+08 1.0 0.0e+00 0.0e+00 0.0e+00 11 72  0  0  0  11 72  0  0  0  7456
MatLUFactorSym         1 1.0 1.0331e-03 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatLUFactorNum         1 1.0 1.5850e-03 1.2 1.20e+06 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  4525
MatAssemblyBegin      11 1.0 2.3642e-03 1.7 0.00e+00 0.0 5.4e+01 6.2e+02 1.6e+01  0  0  0  0  0   0  0  0  0  0     0
MatAssemblyEnd        11 1.0 3.7580e-03 1.1 0.00e+00 0.0 1.4e+02 3.2e+01 3.2e+01  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ            1 1.0 9.4891e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetSubMatrice       2 1.0 6.4421e-04 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+01  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 1.5609e-03 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatZeroEntries         4 1.0 1.8120e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRedundant        1 1.0 5.8007e-04 1.1 0.00e+00 0.0 9.0e+01 7.6e+03 2.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSetup               2 1.0 9.5367e-07 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve            8939 1.0 1.5797e+00 2.0 7.11e+08 1.0 2.7e+05 1.4e+03 0.0e+00 22 72 25 80  0  22 72 25 80  0  2701
PCSetUp                1 1.0 6.0072e-03 1.1 1.20e+06 1.0 1.7e+02 4.4e+03 1.7e+01  0  0  0  0  0   0  0  0  0  0  1194
PCApply             8939 1.0 1.5578e+00 2.0 7.11e+08 1.0 2.7e+05 1.4e+03 0.0e+00 22 72 25 80  0  22 72 25 80  0  2738
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

           Container     1              1          388     0
  Spectral Transform     1              1          568     0
 Eigenproblem Solver     1              1          924     0
       Inner product     1              1          428     0
                 Vec    41             41        95264     0
         Vec Scatter     7              7         6076     0
           Index Set    23             23        43524     0
   IS L to G Mapping     1              1         1608     0
              Matrix    16             16      1002520     0
       Krylov Solver     2              2         1664     0
      Preconditioner     2              2         1416     0
         PetscRandom     1              1          448     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 4.14371e-05
Average time for zero size MPI_Send(): 7.7486e-05
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
| Time:           Fri Aug 24 15:21:14 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=5.02253, Active time=4.75412                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| CondensedEigenSystem                                                                                           |
|   get_eigenpair()                  5         0.0064      0.001286    0.0345      0.006893    0.14     0.72     |
|   solve()                          1         0.0073      0.007292    4.7008      4.700764    0.15     98.88    |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0001      0.000086    0.0001      0.000119    0.00     0.00     |
|   build_sparsity()                 1         0.0007      0.000733    0.0010      0.000956    0.02     0.02     |
|   create_dof_constraints()         1         0.0001      0.000077    0.0001      0.000077    0.00     0.00     |
|   distribute_dofs()                1         0.0003      0.000267    0.0022      0.002215    0.01     0.05     |
|   dof_indices()                    889       0.0003      0.000000    0.0003      0.000000    0.01     0.01     |
|   prepare_send_list()              1         0.0000      0.000009    0.0000      0.000009    0.00     0.00     |
|   reinit()                         1         0.0005      0.000508    0.0005      0.000508    0.01     0.01     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0002      0.000209    0.0008      0.000841    0.00     0.02     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0010      0.000980    0.0010      0.000980    0.02     0.02     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        92        0.0002      0.000002    0.0002      0.000002    0.00     0.00     |
|   init_shape_functions()           1         0.0000      0.000012    0.0000      0.000012    0.00     0.00     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             92        0.0001      0.000001    0.0001      0.000001    0.00     0.00     |
|   init_reference_to_physical_map() 1         0.0000      0.000019    0.0000      0.000019    0.00     0.00     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0003      0.000337    0.0012      0.001179    0.01     0.02     |
|   read()                           1         0.0010      0.000957    0.0010      0.000957    0.02     0.02     |
|   renumber_nodes_and_elem()        2         0.0001      0.000036    0.0001      0.000036    0.00     0.00     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   broadcast()                      1         0.0003      0.000286    0.0009      0.000875    0.01     0.02     |
|   compute_hilbert_indices()        2         0.0020      0.001013    0.0020      0.001013    0.04     0.04     |
|   find_global_indices()            2         0.0002      0.000120    0.0067      0.003331    0.01     0.14     |
|   parallel_sort()                  2         0.0012      0.000613    0.0035      0.001773    0.03     0.07     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0000      0.000027    0.0018      0.001849    0.00     0.04     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0013      0.001313    0.0044      0.004400    0.03     0.09     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      8         0.0016      0.000199    0.0016      0.000199    0.03     0.03     |
|   broadcast()                      7         0.0001      0.000016    0.0001      0.000016    0.00     0.00     |
|   gather()                         1         0.0000      0.000032    0.0000      0.000032    0.00     0.00     |
|   max(scalar)                      3         0.0015      0.000513    0.0015      0.000513    0.03     0.03     |
|   max(vector)                      2         0.0001      0.000047    0.0001      0.000047    0.00     0.00     |
|   min(vector)                      2         0.0003      0.000174    0.0003      0.000174    0.01     0.01     |
|   probe()                          50        0.0022      0.000044    0.0022      0.000044    0.05     0.05     |
|   receive()                        50        0.0001      0.000002    0.0023      0.000046    0.00     0.05     |
|   send()                           50        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   send_receive()                   54        0.0001      0.000002    0.0025      0.000046    0.00     0.05     |
|   sum()                            15        0.0310      0.002064    0.0310      0.002064    0.65     0.65     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           50        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0001      0.000111    0.0016      0.001645    0.00     0.03     |
|   set_parent_processor_ids()       1         0.0000      0.000037    0.0000      0.000037    0.00     0.00     |
|                                                                                                                |
| SlepcEigenSolver                                                                                               |
|   solve_generalized()              1         4.6913      4.691315    4.6913      4.691315    98.68    98.68    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0018      0.001825    0.0022      0.002157    0.04     0.05     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1398      4.7541                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 6 ./eigenproblems_ex3-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
