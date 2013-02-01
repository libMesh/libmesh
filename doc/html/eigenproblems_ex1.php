<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("eigenproblems_ex1",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file eigenproblems_ex1.C with comments: </h1> 
<div class = "comment">
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh.h"
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/exodusII_io.h"
        #include "libmesh/eigen_system.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/fe.h"
        #include "libmesh/quadrature_gauss.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/dof_map.h"
        
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
Skip SLEPc examples on a non-SLEPc libMesh build
</div>

<div class ="fragment">
<pre>
        #ifndef LIBMESH_HAVE_SLEPC
          libmesh_example_assert(false, "--enable-slepc");
        }
        
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
          equation_systems.parameters.set&lt;Real&gt;
            ("linear solver tolerance") = pow(TOLERANCE, 5./3.);
          equation_systems.parameters.set&lt;unsigned int&gt;
            ("linear solver maximum iterations") = 1000;
        
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
        
</pre>
</div>
<div class = "comment">
All done.  
</div>

<div class ="fragment">
<pre>
          return 0;
        }
        
        #endif // LIBMESH_HAVE_SLEPC
        
        
        
        
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
A reference to the system matrix
</div>

<div class ="fragment">
<pre>
          SparseMatrix&lt;Number&gt;&  matrix_A = *eigen_system.matrix_A;
        
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
const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();


<br><br>A reference to the \p DofMap object for this system.  The \p DofMap
object handles the index translation from node and element numbers
to degree of freedom numbers.
</div>

<div class ="fragment">
<pre>
          const DofMap& dof_map = eigen_system.get_dof_map();
        
</pre>
</div>
<div class = "comment">
The element mass matrix.
</div>

<div class ="fragment">
<pre>
          DenseMatrix&lt;Number&gt;   Me;
        
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
Zero the element matrices and rhs before
summing them.  We use the resize member here because
the number of degrees of freedom might have changed from
the last element.  Note that this will be the case if the
element type is different (i.e. the last element was a
triangle, now we are on a quadrilateral).
</div>

<div class ="fragment">
<pre>
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
                      Me(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
                  
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
</div>

<div class ="fragment">
<pre>
              dof_map.constrain_element_matrix(Me, dof_indices, false);
        
</pre>
</div>
<div class = "comment">
Finally, simply add the element contribution to the
overall matrix.
</div>

<div class ="fragment">
<pre>
              matrix_A.add_matrix (Me, dof_indices);
        
        
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
<br><br><br> <h1> The source file eigenproblems_ex1.C without comments: </h1> 
<pre> 
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/eigen_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_mass(EquationSystems&amp; es,
                     <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #ifndef LIBMESH_HAVE_SLEPC
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-slepc&quot;</FONT></B>);
  }
  
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
  
  
    equation_systems.parameters.set&lt;Real&gt;
      (<B><FONT COLOR="#BC8F8F">&quot;linear solver tolerance&quot;</FONT></B>) = pow(TOLERANCE, 5./3.);
    equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;
      (<B><FONT COLOR="#BC8F8F">&quot;linear solver maximum iterations&quot;</FONT></B>) = 1000;
  
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
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  #endif <I><FONT COLOR="#B22222">// LIBMESH_HAVE_SLEPC
</FONT></I>  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_mass(EquationSystems&amp; es,
                     <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Eigensystem&quot;</FONT></B>);
  
  #ifdef LIBMESH_HAVE_SLEPC
  
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    EigenSystem &amp; eigen_system = es.get_system&lt;EigenSystem&gt; (system_name);
  
    FEType fe_type = eigen_system.get_dof_map().variable_type(0);
  
    SparseMatrix&lt;Number&gt;&amp;  matrix_A = *eigen_system.matrix_A;
  
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    QGauss qrule (dim, fe_type.default_quadrature_order());
  
    fe-&gt;attach_quadrature_rule (&amp;qrule);
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
  
  
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = eigen_system.get_dof_map();
  
    DenseMatrix&lt;Number&gt;   Me;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
  
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
   
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
  
        fe-&gt;reinit (elem);
  
        Me.resize (dof_indices.size(), dof_indices.size());
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
                Me(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
            
        dof_map.constrain_element_matrix(Me, dof_indices, false);
  
        matrix_A.add_matrix (Me, dof_indices);
  
  
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
* Running Example eigenproblems_ex1:
*  mpirun -np 12 example-devel -n 5 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/eigenproblems/eigenproblems_ex1/.libs/lt-example-devel -n 5 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=441
    n_local_nodes()=49
  n_elem()=400
    n_local_elem()=34
    n_active_elem()=400
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
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=441
    n_local_dofs()=49
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=0
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 7.76417
      Average Off-Processor Bandwidth <= 1.51474
      Maximum  On-Processor Bandwidth <= 11
      Maximum Off-Processor Bandwidth <= 7
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

Number of converged eigenpairs: 5

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/eigenproblems/eigenproblems_ex1/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:05:02 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           2.029e-01      1.00221   2.028e-01
Objects:              4.400e+01      1.00000   4.400e+01
Flops:                2.565e+06      1.80568   1.937e+06  2.325e+07
Flops/sec:            1.265e+07      1.80565   9.554e+06  1.146e+08
MPI Messages:         3.936e+03      3.49867   2.343e+03  2.812e+04
MPI Message Lengths:  1.379e+05      2.24604   4.312e+01  1.212e+06
MPI Reductions:       1.192e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.0275e-01 100.0%  2.3249e+07 100.0%  2.812e+04 100.0%  4.312e+01      100.0%  1.191e+03  99.9% 

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

STSetUp                1 1.0 4.9114e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  0   0  0  0  0  0     0
STApply              557 1.0 1.9524e-02 6.1 4.12e+05 1.7 2.8e+04 4.3e+01 0.0e+00  5 17 99 98  0   5 17 99 98  0   200
EPSSetUp               1 1.0 8.1396e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 3.6e+01  0  0  0  0  3   0  0  0  0  3     0
EPSSolve               1 1.0 5.9568e-02 1.0 2.56e+06 1.8 2.8e+04 4.3e+01 1.1e+03 29100 99 98 95  29100 99 98 95   390
IPOrthogonalize      558 1.0 4.9822e-02 1.5 2.06e+06 1.8 0.0e+00 0.0e+00 1.1e+03 21 79  0  0 95  21 79  0  0 95   371
IPInnerProduct      4482 1.0 4.7614e-02 1.5 1.08e+06 1.8 0.0e+00 0.0e+00 1.1e+03 20 42  0  0 95  20 42  0  0 95   203
DSSolve               43 1.0 3.7587e-03 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  2  0  0  0  0   2  0  0  0  0     0
DSVectors             49 1.0 3.6597e-04 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
DSOther               86 1.0 9.8753e-04 2.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
UpdateVectors         44 1.0 2.2769e-04 1.7 7.30e+04 1.8 0.0e+00 0.0e+00 0.0e+00  0  3  0  0  0   0  3  0  0  0  2886
VecMAXPBY           1114 1.0 1.2579e-03 1.4 9.78e+05 1.8 0.0e+00 0.0e+00 0.0e+00  1 38  0  0  0   1 38  0  0  0  6998
VecScale             515 1.0 2.8205e-04 1.6 2.52e+04 1.8 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0   805
VecCopy                2 1.0 1.7166e-05 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                 3 1.0 1.1921e-05 2.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAssemblyBegin       1 1.0 2.8610e-04 2.6 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAssemblyEnd         1 1.0 2.7895e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      557 1.0 2.0373e-03 1.7 0.00e+00 0.0 2.8e+04 4.3e+01 0.0e+00  1  0 99 98  0   1  0 99 98  0     0
VecScatterEnd        557 1.0 1.6396e-0218.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  4  0  0  0  0   4  0  0  0  0     0
VecReduceArith      2241 1.0 1.1910e-02 1.0 1.08e+06 1.8 0.0e+00 0.0e+00 0.0e+00  6 42  0  0  0   6 42  0  0  0   811
VecReduceComm       1127 1.0 3.4664e-02 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.1e+03 14  0  0  0 95  14  0  0  0 95     0
MatMult              557 1.0 1.9291e-02 6.4 4.12e+05 1.7 2.8e+04 4.3e+01 0.0e+00  5 17 99 98  0   5 17 99 98  0   202
MatAssemblyBegin       1 1.0 1.4114e-04 1.1 0.00e+00 0.0 6.9e+01 2.2e+02 2.0e+00  0  0  0  1  0   0  0  0  1  0     0
MatAssemblyEnd         1 1.0 5.3811e-04 1.0 0.00e+00 0.0 1.0e+02 1.3e+01 8.0e+00  0  0  0  0  1   0  0  0  0  1     0
MatZeroEntries         2 1.0 1.4067e-05 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

           Container     2              2         1096     0
  Spectral Transform     1              1          736     0
 Eigenproblem Solver     1              1         1680     0
       Inner product     1              1          624     0
       Direct solver     1              1         9408     0
              Vector    24             24        38008     0
      Vector Scatter     2              2         2072     0
           Index Set     4              4         3096     0
   IS L to G Mapping     1              1          564     0
              Matrix     3              3        14500     0
         PetscRandom     1              1          608     0
       Krylov Solver     1              1         1072     0
      Preconditioner     1              1          856     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 1.19209e-07
Average time for MPI_Barrier(): 4.57764e-06
Average time for zero size MPI_Send(): 1.34309e-05
#PETSc Option Table entries:
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
| Time:           Thu Jan 31 22:05:02 2013                                                                             |
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
| libMesh Performance: Alive time=0.282838, Active time=0.186116                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0021      0.002126    0.0037      0.003688    1.14     1.98     |
|   build_sparsity()                 1         0.0019      0.001877    0.0053      0.005329    1.01     2.86     |
|   create_dof_constraints()         1         0.0008      0.000809    0.0008      0.000809    0.43     0.43     |
|   distribute_dofs()                1         0.0075      0.007491    0.0245      0.024454    4.02     13.14    |
|   dof_indices()                    131       0.0071      0.000054    0.0071      0.000054    3.80     3.80     |
|   prepare_send_list()              1         0.0000      0.000037    0.0000      0.000037    0.02     0.02     |
|   reinit()                         1         0.0153      0.015280    0.0153      0.015280    8.21     8.21     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0008      0.000772    0.0031      0.003062    0.41     1.65     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0054      0.005431    0.0054      0.005431    2.92     2.92     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        34        0.0003      0.000008    0.0003      0.000008    0.15     0.15     |
|   init_shape_functions()           1         0.0001      0.000077    0.0001      0.000077    0.04     0.04     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             34        0.0004      0.000010    0.0004      0.000010    0.19     0.19     |
|   init_reference_to_physical_map() 1         0.0001      0.000096    0.0001      0.000096    0.05     0.05     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0093      0.009268    0.0095      0.009453    4.98     5.08     |
|   renumber_nodes_and_elem()        2         0.0004      0.000192    0.0004      0.000192    0.21     0.21     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0068      0.003413    0.0068      0.003413    3.67     3.67     |
|   find_global_indices()            2         0.0030      0.001497    0.0145      0.007229    1.61     7.77     |
|   parallel_sort()                  2         0.0028      0.001418    0.0033      0.001668    1.52     1.79     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0002      0.000159    0.0088      0.008778    0.09     4.72     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0022      0.002206    0.0022      0.002206    1.19     1.19     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0405      0.040546    0.0469      0.046942    21.79    25.22    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      9         0.0006      0.000063    0.0006      0.000068    0.30     0.33     |
|   max(bool)                        1         0.0000      0.000008    0.0000      0.000008    0.00     0.00     |
|   max(scalar)                      105       0.0011      0.000010    0.0011      0.000010    0.59     0.59     |
|   max(vector)                      24        0.0004      0.000017    0.0012      0.000048    0.22     0.62     |
|   min(bool)                        121       0.0012      0.000010    0.0012      0.000010    0.65     0.65     |
|   min(scalar)                      99        0.0065      0.000066    0.0065      0.000066    3.49     3.49     |
|   min(vector)                      24        0.0005      0.000020    0.0013      0.000054    0.26     0.69     |
|   probe()                          132       0.0010      0.000008    0.0010      0.000008    0.54     0.54     |
|   receive()                        132       0.0008      0.000006    0.0019      0.000014    0.45     1.00     |
|   send()                           132       0.0004      0.000003    0.0004      0.000003    0.22     0.22     |
|   send_receive()                   136       0.0012      0.000009    0.0038      0.000028    0.62     2.04     |
|   sum()                            20        0.0004      0.000022    0.0008      0.000042    0.23     0.45     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           132       0.0003      0.000002    0.0003      0.000002    0.15     0.15     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0012      0.001192    0.0020      0.001981    0.64     1.06     |
|   set_parent_processor_ids()       1         0.0008      0.000788    0.0008      0.000788    0.42     0.42     |
|                                                                                                                |
| SlepcEigenSolver                                                                                               |
|   solve_standard()                 1         0.0620      0.062022    0.0620      0.062022    33.32    33.32    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0008      0.000795    0.0035      0.003507    0.43     1.88     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1292      0.1861                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example eigenproblems_ex1:
*  mpirun -np 12 example-devel -n 5 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
