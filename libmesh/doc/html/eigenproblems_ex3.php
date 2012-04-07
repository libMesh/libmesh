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
Compiling C++ (in optimized mode) eigenproblems_ex3.C...
Linking eigenproblems_ex3-opt...
***************************************************************
* Running Example  ./eigenproblems_ex3-opt
***************************************************************
 
Running ./eigenproblems_ex3-opt -n_evals 5 -mesh_name drum1 -plotting_index 2

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1157
    n_local_nodes()=1157
  n_elem()=536
    n_local_elem()=536
    n_active_elem()=536
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
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
    n_local_dofs()=1157
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=0
    n_matrices()=2
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 10.9464
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 22
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

Number of converged eigenpairs: 5


-------------------------------------------------------------------
| Time:           Sat Apr  7 15:58:15 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=3.01797, Active time=2.94393                                               |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| CondensedEigenSystem                                                                                       |
|   get_eigenpair()              5         0.0017      0.000334    0.0017      0.000334    0.06     0.06     |
|   solve()                      1         0.0014      0.001439    2.9287      2.928706    0.05     99.48    |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 1         0.0004      0.000357    0.0004      0.000357    0.01     0.01     |
|   compute_sparsity()           1         0.0019      0.001939    0.0024      0.002397    0.07     0.08     |
|   create_dof_constraints()     1         0.0004      0.000409    0.0004      0.000409    0.01     0.01     |
|   distribute_dofs()            1         0.0007      0.000744    0.0026      0.002554    0.03     0.09     |
|   dof_indices()                2144      0.0017      0.000001    0.0017      0.000001    0.06     0.06     |
|   prepare_send_list()          1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                     1         0.0018      0.001808    0.0018      0.001808    0.06     0.06     |
|                                                                                                            |
| EquationSystems                                                                                            |
|   build_solution_vector()      1         0.0008      0.000841    0.0014      0.001383    0.03     0.05     |
|                                                                                                            |
| ExodusII_IO                                                                                                |
|   write_nodal_data()           1         0.0011      0.001124    0.0011      0.001124    0.04     0.04     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         536       0.0006      0.000001    0.0006      0.000001    0.02     0.02     |
|   compute_shape_functions()    536       0.0004      0.000001    0.0004      0.000001    0.01     0.01     |
|   init_shape_functions()       1         0.0001      0.000063    0.0001      0.000063    0.00     0.00     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             1         0.0014      0.001449    0.0014      0.001449    0.05     0.05     |
|   read()                       1         0.0034      0.003401    0.0034      0.003401    0.12     0.12     |
|   renumber_nodes_and_elem()    2         0.0002      0.000097    0.0002      0.000097    0.01     0.01     |
|                                                                                                            |
| MeshOutput                                                                                                 |
|   write_equation_systems()     1         0.0000      0.000032    0.0025      0.002540    0.00     0.09     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  1         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|                                                                                                            |
| Partitioner                                                                                                |
|   single_partition()           1         0.0001      0.000135    0.0001      0.000135    0.00     0.00     |
|                                                                                                            |
| SlepcEigenSolver                                                                                           |
|   solve_generalized()          1         2.9213      2.921316    2.9213      2.921316    99.23    99.23    |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   1         0.0043      0.004326    0.0060      0.005951    0.15     0.20     |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        3240      2.9439                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

Running ./eigenproblems_ex3-opt -n_evals 5 -mesh_name drum2 -plotting_index 2

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1189
    n_local_nodes()=1189
  n_elem()=552
    n_local_elem()=552
    n_active_elem()=552
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
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
    n_local_dofs()=1189
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=0
    n_matrices()=2
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 10.9613
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 22
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

Number of converged eigenpairs: 5


-------------------------------------------------------------------
| Time:           Sat Apr  7 15:58:18 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=3.00282, Active time=2.92409                                               |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| CondensedEigenSystem                                                                                       |
|   get_eigenpair()              5         0.0017      0.000346    0.0017      0.000346    0.06     0.06     |
|   solve()                      1         0.0029      0.002938    2.9047      2.904662    0.10     99.34    |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 1         0.0005      0.000482    0.0005      0.000482    0.02     0.02     |
|   compute_sparsity()           1         0.0039      0.003877    0.0058      0.005791    0.13     0.20     |
|   create_dof_constraints()     1         0.0003      0.000278    0.0003      0.000278    0.01     0.01     |
|   distribute_dofs()            1         0.0008      0.000753    0.0026      0.002612    0.03     0.09     |
|   dof_indices()                2208      0.0044      0.000002    0.0044      0.000002    0.15     0.15     |
|   prepare_send_list()          1         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   reinit()                     1         0.0019      0.001858    0.0019      0.001858    0.06     0.06     |
|                                                                                                            |
| EquationSystems                                                                                            |
|   build_solution_vector()      1         0.0009      0.000866    0.0015      0.001466    0.03     0.05     |
|                                                                                                            |
| ExodusII_IO                                                                                                |
|   write_nodal_data()           1         0.0011      0.001134    0.0011      0.001134    0.04     0.04     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         552       0.0013      0.000002    0.0013      0.000002    0.04     0.04     |
|   compute_shape_functions()    552       0.0009      0.000002    0.0009      0.000002    0.03     0.03     |
|   init_shape_functions()       1         0.0001      0.000103    0.0001      0.000103    0.00     0.00     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             1         0.0015      0.001507    0.0015      0.001507    0.05     0.05     |
|   read()                       1         0.0035      0.003489    0.0035      0.003489    0.12     0.12     |
|   renumber_nodes_and_elem()    2         0.0002      0.000104    0.0002      0.000104    0.01     0.01     |
|                                                                                                            |
| MeshOutput                                                                                                 |
|   write_equation_systems()     1         0.0000      0.000035    0.0026      0.002635    0.00     0.09     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                            |
| Partitioner                                                                                                |
|   single_partition()           1         0.0001      0.000138    0.0001      0.000138    0.00     0.00     |
|                                                                                                            |
| SlepcEigenSolver                                                                                           |
|   solve_generalized()          1         2.8884      2.888437    2.8884      2.888437    98.78    98.78    |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   1         0.0096      0.009646    0.0133      0.013285    0.33     0.45     |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        3336      2.9241                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./eigenproblems_ex3-opt
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
