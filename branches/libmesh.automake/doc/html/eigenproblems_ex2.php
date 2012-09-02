<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("eigenproblems_ex2",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh.h"
        #include "mesh.h"
        #include "mesh_generation.h"
        #include "exodusII_io.h"
        #include "eigen_system.h"
        #include "equation_systems.h"
        #include "fe.h"
        #include "quadrature_gauss.h"
        #include "dense_matrix.h"
        #include "sparse_matrix.h"
        #include "numeric_vector.h"
        #include "dof_map.h"
        
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
        void assemble_mass(EquationSystems& es,
                           const std::string& system_name);
        
        
        
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
        
</pre>
</div>
<div class = "comment">
Check for proper usage.
</div>

<div class ="fragment">
<pre>
          if (argc &lt; 3)
            {
              if (libMesh::processor_id() == 0)
                std::cerr &lt;&lt; "\nUsage: " &lt;&lt; argv[0]
                          &lt;&lt; " -n &lt;number of eigen values&gt;"
                          &lt;&lt; std::endl;
              libmesh_error();
            }
          
</pre>
</div>
<div class = "comment">
Tell the user what we are doing.
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
Get the number of eigen values to be computed from argv[2]
</div>

<div class ="fragment">
<pre>
          const unsigned int nev = std::atoi(argv[2]);
        
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
Create a mesh.
</div>

<div class ="fragment">
<pre>
          Mesh mesh;
        
</pre>
</div>
<div class = "comment">
Use the internal mesh generator to create a uniform
2D grid on a square.
</div>

<div class ="fragment">
<pre>
          MeshTools::Generation::build_square (mesh, 
                                               20, 20,
                                               -1., 1.,
                                               -1., 1.,
                                               QUAD4);
        
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
Create a EigenSystem named "Eigensystem" and (for convenience)
use a reference to the system we create.
</div>

<div class ="fragment">
<pre>
          EigenSystem & eigen_system =
            equation_systems.add_system&lt;EigenSystem&gt; ("Eigensystem");
        
</pre>
</div>
<div class = "comment">
Declare the system variables.
Adds the variable "p" to "Eigensystem".   "p"
will be approximated using second-order approximation.
</div>

<div class ="fragment">
<pre>
          eigen_system.add_variable("p", FIRST);
        
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly
function defined below.
</div>

<div class ="fragment">
<pre>
          eigen_system.attach_assemble_function (assemble_mass);
        
</pre>
</div>
<div class = "comment">
Set necessary parametrs used in EigenSystem::solve(),
i.e. the number of requested eigenpairs \p nev and the number
of basis vectors \p ncv used in the solution algorithm. Note that
ncv >= nev must hold and ncv >= 2*nev is recommended.
</div>

<div class ="fragment">
<pre>
          equation_systems.parameters.set&lt;unsigned int&gt;("eigenpairs")    = nev;
          equation_systems.parameters.set&lt;unsigned int&gt;("basis vectors") = nev*3;
        
</pre>
</div>
<div class = "comment">
You may optionally change the default eigensolver used by SLEPc. 
The Krylov-Schur method is mathematically equivalent to implicitly
restarted Arnoldi, the method of Arpack, so there is currently no
point in using SLEPc with Arpack.
ARNOLDI     = default in SLEPc 2.3.1 and earlier
KRYLOVSCHUR default in SLEPc 2.3.2 and later
eigen_system.eigen_solver->set_eigensolver_type(KRYLOVSCHUR); 


<br><br>Set the solver tolerance and the maximum number of iterations. 
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
Set the eigenvalues to be computed. Note that not
all solvers support this.
eigen_system.eigen_solver->set_position_of_spectrum(SMALLEST_MAGNITUDE);


<br><br>Initialize the data structures for the equation system.
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
        
</pre>
</div>
<div class = "comment">
Get the last converged eigenpair
</div>

<div class ="fragment">
<pre>
          if (nconv != 0)
            {
              eigen_system.get_eigenpair(nconv-1);
              
        #ifdef LIBMESH_HAVE_EXODUS_API
</pre>
</div>
<div class = "comment">
Write the eigen vector to file.
</div>

<div class ="fragment">
<pre>
              ExodusII_IO (mesh).write_equation_systems ("out.e", equation_systems);
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
            }
          else
            {
              std::cout &lt;&lt; "WARNING: Solver did not converge!\n" &lt;&lt; nconv &lt;&lt; std::endl;
            }
        
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
        
        
        
        void assemble_mass(EquationSystems& es,
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
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;eigen_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_mass(EquationSystems&amp; es,
                     <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
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
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (argc &lt; 3)
      {
        <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\nUsage: &quot;</FONT></B> &lt;&lt; argv[0]
                    &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; -n &lt;number of eigen values&gt;&quot;</FONT></B>
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
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> nev = std::atoi(argv[2]);
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
    
    Mesh mesh;
  
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square (mesh, 
                                         20, 20,
                                         -1., 1.,
                                         -1., 1.,
                                         QUAD4);
  
    mesh.print_info();
    
    EquationSystems equation_systems (mesh);
  
    EigenSystem &amp; eigen_system =
      equation_systems.add_system&lt;EigenSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Eigensystem&quot;</FONT></B>);
  
    eigen_system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>, FIRST);
  
    eigen_system.attach_assemble_function (assemble_mass);
  
    equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;eigenpairs&quot;</FONT></B>)    = nev;
    equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;basis vectors&quot;</FONT></B>) = nev*3;
  
  
    equation_systems.parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;linear solver tolerance&quot;</FONT></B>) = pow(TOLERANCE, 5./3.);
    equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;
      (<B><FONT COLOR="#BC8F8F">&quot;linear solver maximum iterations&quot;</FONT></B>) = 1000;
  
    eigen_system.set_eigenproblem_type(GHEP);
  
  
    equation_systems.init();
  
    equation_systems.print_info();
       
    eigen_system.solve();
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> nconv = eigen_system.get_n_converged();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Number of converged eigenpairs: &quot;</FONT></B> &lt;&lt; nconv
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B> &lt;&lt; std::endl;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (nconv != 0)
      {
        eigen_system.get_eigenpair(nconv-1);
        
  #ifdef LIBMESH_HAVE_EXODUS_API
        ExodusII_IO (mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;out.e&quot;</FONT></B>, equation_systems);
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>      }
    <B><FONT COLOR="#A020F0">else</FONT></B>
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;WARNING: Solver did not converge!\n&quot;</FONT></B> &lt;&lt; nconv &lt;&lt; std::endl;
      }
  
  #endif <I><FONT COLOR="#B22222">// LIBMESH_HAVE_SLEPC
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_mass(EquationSystems&amp; es,
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
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Linking eigenproblems_ex2-opt...
***************************************************************
* Running Example  mpirun -np 6 ./eigenproblems_ex2-opt -n 5 -eps_type lapack -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running ./eigenproblems_ex2-opt -n 5 -eps_type lapack -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=441
    n_local_nodes()=86
  n_elem()=400
    n_local_elem()=67
    n_active_elem()=400
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
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=441
    n_local_dofs()=86
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=0
    n_matrices()=2
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.69767
      Average Off-Processor Bandwidth <= 0.686047
      Maximum  On-Processor Bandwidth <= 9
      Maximum Off-Processor Bandwidth <= 5
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

Number of converged eigenpairs: 441

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./eigenproblems_ex2-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:21:03 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           2.334e-01      1.13060   2.110e-01
Objects:              4.820e+02      1.00000   4.820e+02
Flops:                0.000e+00      0.00000   0.000e+00  0.000e+00
Flops/sec:            0.000e+00      0.00000   0.000e+00  0.000e+00
MPI Messages:         4.500e+01      2.50000   2.900e+01  1.740e+02
MPI Message Lengths:  4.932e+03      1.83482   1.189e+02  2.068e+04
MPI Reductions:       9.570e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.1094e-01 100.0%  0.0000e+00   0.0%  1.740e+02 100.0%  1.189e+02      100.0%  4.790e+02  50.1% 

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

STSetUp                1 1.0 9.0837e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
EPSSetUp               1 1.0 1.3653e-02 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 4.5e+02  6  0  0  0 47   6  0  0  0 94     0
EPSSolve               1 1.0 1.6237e-01 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 51  0  0  0  0  51  0  0  0  0     0
EPSDense               1 1.0 1.6211e-01 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 51  0  0  0  0  51  0  0  0  0     0
VecCopy                1 1.0 3.0041e-05 9.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                 2 1.0 8.1062e-06 3.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAssemblyBegin       1 1.0 7.3297e-02632.6 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00 26  0  0  0  0  26  0  0  0  1     0
VecAssemblyEnd         1 1.0 1.1802e-04 3.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatConvert             2 1.0 2.3010e-03 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
MatAssemblyBegin       6 1.0 1.1261e-03 3.6 0.00e+00 0.0 5.4e+01 3.4e+02 4.0e+00  0  0 31 88  0   0  0 31 88  1     0
MatAssemblyEnd         6 1.0 1.4260e-03 1.2 0.00e+00 0.0 8.0e+01 1.7e+01 1.4e+01  1  0 46  7  1   1  0 46  7  3     0
MatGetRow            882 1.0 2.0576e-04 4.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetSubMatrice       2 1.0 5.8062e-03 2.4 0.00e+00 0.0 0.0e+00 0.0e+00 6.0e+00  2  0  0  0  1   2  0  0  0  1     0
MatZeroEntries         4 1.0 1.1921e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSetup               1 1.0 0.0000e+00 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
PCSetUp                1 1.0 1.1921e-06 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

  Spectral Transform     1              1          568     0
 Eigenproblem Solver     1              1          948     0
       Inner product     1              1          428     0
                 Vec   450            450       607912     0
         Vec Scatter     3              3         2604     0
           Index Set    10             10         5436     0
   IS L to G Mapping     1              1          840     0
              Matrix    12             10      6267208     0
       Krylov Solver     1              1          832     0
      Preconditioner     1              1          664     0
         PetscRandom     1              1          448     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.32014e-05
Average time for zero size MPI_Send(): 6.09954e-05
#PETSc Option Table entries:
-eps_type lapack
-ksp_right_pc
-log_summary
-n 5
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
| Time:           Fri Aug 24 15:21:03 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.33037, Active time=0.196947                                                  |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0001      0.000057    0.0001      0.000077    0.03     0.04     |
|   build_sparsity()                 1         0.0004      0.000411    0.0005      0.000547    0.21     0.28     |
|   create_dof_constraints()         1         0.0001      0.000057    0.0001      0.000057    0.03     0.03     |
|   distribute_dofs()                1         0.0001      0.000148    0.0012      0.001160    0.08     0.59     |
|   dof_indices()                    575       0.0002      0.000000    0.0002      0.000000    0.10     0.10     |
|   prepare_send_list()              1         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   reinit()                         1         0.0003      0.000257    0.0003      0.000257    0.13     0.13     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0002      0.000184    0.0016      0.001637    0.09     0.83     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0015      0.001460    0.0015      0.001460    0.74     0.74     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        67        0.0001      0.000001    0.0001      0.000001    0.04     0.04     |
|   init_shape_functions()           1         0.0000      0.000012    0.0000      0.000012    0.01     0.01     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             67        0.0000      0.000001    0.0000      0.000001    0.02     0.02     |
|   init_reference_to_physical_map() 1         0.0000      0.000027    0.0000      0.000027    0.01     0.01     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0004      0.000361    0.0014      0.001380    0.18     0.70     |
|   renumber_nodes_and_elem()        2         0.0000      0.000019    0.0000      0.000019    0.02     0.02     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0015      0.000752    0.0015      0.000752    0.76     0.76     |
|   find_global_indices()            2         0.0002      0.000117    0.0058      0.002909    0.12     2.95     |
|   parallel_sort()                  2         0.0023      0.001154    0.0037      0.001861    1.17     1.89     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0000      0.000032    0.0031      0.003129    0.02     1.59     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0002      0.000199    0.0002      0.000199    0.10     0.10     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0021      0.002077    0.0050      0.004952    1.05     2.51     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      8         0.0006      0.000071    0.0006      0.000071    0.29     0.29     |
|   broadcast()                      1         0.0000      0.000011    0.0000      0.000011    0.01     0.01     |
|   gather()                         1         0.0000      0.000005    0.0000      0.000005    0.00     0.00     |
|   max(scalar)                      2         0.0010      0.000512    0.0010      0.000512    0.52     0.52     |
|   max(vector)                      2         0.0000      0.000018    0.0000      0.000018    0.02     0.02     |
|   min(vector)                      2         0.0018      0.000918    0.0018      0.000918    0.93     0.93     |
|   probe()                          50        0.0021      0.000042    0.0021      0.000042    1.07     1.07     |
|   receive()                        50        0.0001      0.000002    0.0022      0.000044    0.04     1.11     |
|   send()                           50        0.0000      0.000001    0.0000      0.000001    0.02     0.02     |
|   send_receive()                   54        0.0001      0.000002    0.0024      0.000044    0.05     1.21     |
|   sum()                            10        0.0030      0.000299    0.0030      0.000299    1.52     1.52     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           50        0.0000      0.000001    0.0000      0.000001    0.02     0.02     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0001      0.000063    0.0016      0.001615    0.03     0.82     |
|   set_parent_processor_ids()       1         0.0000      0.000027    0.0000      0.000027    0.01     0.01     |
|                                                                                                                |
| SlepcEigenSolver                                                                                               |
|   solve_generalized()              1         0.1780      0.178015    0.1780      0.178015    90.39    90.39    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0003      0.000315    0.0005      0.000522    0.16     0.27     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1015      0.1969                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 6 ./eigenproblems_ex2-opt -n 5 -eps_type lapack -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
