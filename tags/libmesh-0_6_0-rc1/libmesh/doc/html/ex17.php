<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex17",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
  

<br><br>This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
  

<br><br>This library is distributd in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
  

<br><br>You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


<br><br>

<br><br>C++ include files that we need


<br><br>

<br><br>libMesh include files.
</div>

<div class ="fragment">
<pre>
        #include "libmesh.h"
        #include "mesh.h"
        #include "mesh_generation.h"
        #include "gmv_io.h"
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
<h1>Example 17 - Solving a generalized Eigen Problem</h1>

<br><br>This example shows how the previous EigenSolver example 
can be adapted to solve generailzed eigenvalue problems.

<br><br>For solving eigen problems, libMesh interfaces
SLEPc (www.grycap.upv.es/slepc/) which again is based on PETSc.
Hence, this example will only work if the library is compiled
with SLEPc support enabled.

<br><br>In this example some eigenvalues for a generalized symmetric
eigenvalue problem A*x=lambda*B*x are computed, where the
matrices A and B are assembled according to stiffness and
mass matrix, respectively.


<br><br>

<br><br>Function prototype.  This is the function that will assemble
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
          libMesh::init (argc, argv);
        
</pre>
</div>
<div class = "comment">
This example is designed for the SLEPc eigen solver interface.
</div>

<div class ="fragment">
<pre>
        #ifndef HAVE_SLEPC
        
          std::cerr &lt;&lt; "ERROR: This example requires libMesh to be\n"
        	    &lt;&lt; "compiled with SLEPc eigen solvers support!"
        	    &lt;&lt; std::endl;
        
          return 0;
        #else
        
        
</pre>
</div>
<div class = "comment">
Braces are used to force object scope.  
</div>

<div class ="fragment">
<pre>
          {
</pre>
</div>
<div class = "comment">
Check for proper usage.
</div>

<div class ="fragment">
<pre>
            if (argc &lt; 3)
              {
        	std::cerr &lt;&lt; "\nUsage: " &lt;&lt; argv[0]
        		  &lt;&lt; " -n &lt;number of eigen values&gt;"
        		  &lt;&lt; std::endl;
        	error();
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
Set the dimensionality.
</div>

<div class ="fragment">
<pre>
            const unsigned int dim = 2;
        
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
Create a dim-dimensional mesh.
</div>

<div class ="fragment">
<pre>
            Mesh mesh (dim);
        
</pre>
</div>
<div class = "comment">
Use the internal mesh generator to create a uniform
grid on a square.
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
</div>

<div class ="fragment">
<pre>
            {
</pre>
</div>
<div class = "comment">
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
Set the eigen solver type. SLEPc offers various solvers such as
the Arnoldi and subspace method. It
also offers interfaces to other solver packages (e.g. ARPACK).
</div>

<div class ="fragment">
<pre>
              eigen_system.eigen_solver-&gt;set_eigensolver_type(ARNOLDI);
        
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
        
            }
               
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
        	
</pre>
</div>
<div class = "comment">
Write the eigen vector to file.
</div>

<div class ="fragment">
<pre>
                char buf[14];
        	sprintf (buf, "out.gmv");
        	GMVIO (mesh).write_equation_systems (buf, equation_systems);
              }
            else
              {
        	std::cout &lt;&lt; "WARNING: Solver did not converge!\n" &lt;&lt; nconv &lt;&lt; std::endl;
              }
          }
        
        #endif // HAVE_SLEPC
        
</pre>
</div>
<div class = "comment">
All done.  
</div>

<div class ="fragment">
<pre>
          return libMesh::close ();
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
          assert (system_name == "Eigensystem");
        
        #ifdef HAVE_SLEPC
        
</pre>
</div>
<div class = "comment">
Get a constant reference to the mesh object.
</div>

<div class ="fragment">
<pre>
          const Mesh& mesh = es.get_mesh();
        
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
Now we will loop over all the elements in the mesh.
We will compute the element matrix contribution.


<br><br></div>

<div class ="fragment">
<pre>
          MeshBase::const_element_iterator       el     = mesh.elements_begin();
          const MeshBase::const_element_iterator end_el = mesh.elements_end();
         
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
Finally, simply add the element contribution to the
overall matrices A and B.
</div>

<div class ="fragment">
<pre>
              matrix_A.add_matrix (Ke, dof_indices);
              matrix_B.add_matrix (Me, dof_indices);
        
            } // end of element loop
        
        
        #endif // HAVE_SLEPC
        
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
    
    
    
  
  
  
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;libmesh.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;mesh.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;mesh_generation.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;gmv_io.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;eigen_system.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;equation_systems.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;fe.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_matrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;sparse_matrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;numeric_vector.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dof_map.h&quot;</FONT></B>
  
  
  
  
  
  <FONT COLOR="#228B22"><B>void</FONT></B> assemble_mass(EquationSystems&amp; es,
  		   <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name);
  
  
  
  <FONT COLOR="#228B22"><B>int</FONT></B> main (<FONT COLOR="#228B22"><B>int</FONT></B> argc, <FONT COLOR="#228B22"><B>char</FONT></B>** argv)
  {
    libMesh::init (argc, argv);
  
  #ifndef HAVE_SLEPC
  
    std::cerr &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;ERROR: This example requires libMesh to be\n&quot;</FONT></B>
  	    &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;compiled with SLEPc eigen solvers support!&quot;</FONT></B>
  	    &lt;&lt; std::endl;
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
  
    {
      <B><FONT COLOR="#A020F0">if</FONT></B> (argc &lt; 3)
        {
  	std::cerr &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;\nUsage: &quot;</FONT></B> &lt;&lt; argv[0]
  		  &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; -n &lt;number of eigen values&gt;&quot;</FONT></B>
  		  &lt;&lt; std::endl;
  	error();
        }
      
      <B><FONT COLOR="#A020F0">else</FONT></B> 
        {
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Running &quot;</FONT></B> &lt;&lt; argv[0];
  	
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>int</FONT></B> i=1; i&lt;argc; i++)
  	  std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; &quot;</FONT></B> &lt;&lt; argv[i];
  	
  	std::cout &lt;&lt; std::endl &lt;&lt; std::endl;
        }
  
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = 2;
  
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> nev = std::atoi(argv[2]);
  
      Mesh mesh (dim);
  
      MeshTools::Generation::build_square (mesh, 
  					 20, 20,
  					 -1., 1.,
  					 -1., 1.,
  					 QUAD4);
  
      mesh.print_info();
      
      EquationSystems equation_systems (mesh);
  
      EigenSystem &amp; eigen_system =
        equation_systems.add_system&lt;EigenSystem&gt; (<FONT COLOR="#BC8F8F"><B>&quot;Eigensystem&quot;</FONT></B>);
  
      {
        eigen_system.add_variable(<FONT COLOR="#BC8F8F"><B>&quot;p&quot;</FONT></B>, FIRST);
  
        eigen_system.attach_assemble_function (assemble_mass);
  
        equation_systems.parameters.set&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt;(<FONT COLOR="#BC8F8F"><B>&quot;eigenpairs&quot;</FONT></B>)    = nev;
        equation_systems.parameters.set&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt;(<FONT COLOR="#BC8F8F"><B>&quot;basis vectors&quot;</FONT></B>) = nev*3;
  
        eigen_system.eigen_solver-&gt;set_eigensolver_type(ARNOLDI);
  
        equation_systems.parameters.set&lt;Real&gt;(<FONT COLOR="#BC8F8F"><B>&quot;linear solver tolerance&quot;</FONT></B>) = pow(TOLERANCE, 5./3.);
        equation_systems.parameters.set&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt;
  	(<FONT COLOR="#BC8F8F"><B>&quot;linear solver maximum iterations&quot;</FONT></B>) = 1000;
  
        eigen_system.set_eigenproblem_type(GHEP);
  
  
        equation_systems.init();
  
        equation_systems.print_info();
  
      }
         
      eigen_system.solve();
  
      <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> nconv = eigen_system.get_n_converged();
  
      std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Number of converged eigenpairs: &quot;</FONT></B> &lt;&lt; nconv
  	      &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;\n&quot;</FONT></B> &lt;&lt; std::endl;
  
      <B><FONT COLOR="#A020F0">if</FONT></B> (nconv != 0)
        {
  	eigen_system.get_eigenpair(nconv-1);
  	
  	<FONT COLOR="#228B22"><B>char</FONT></B> buf[14];
  	sprintf (buf, <FONT COLOR="#BC8F8F"><B>&quot;out.gmv&quot;</FONT></B>);
  	GMVIO (mesh).write_equation_systems (buf, equation_systems);
        }
      <B><FONT COLOR="#A020F0">else</FONT></B>
        {
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;WARNING: Solver did not converge!\n&quot;</FONT></B> &lt;&lt; nconv &lt;&lt; std::endl;
        }
    }
  
  #endif <I><FONT COLOR="#B22222">// HAVE_SLEPC
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close ();
  }
  
  
  
  
  <FONT COLOR="#228B22"><B>void</FONT></B> assemble_mass(EquationSystems&amp; es,
  		   <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name)
  {
    
    assert (system_name == <FONT COLOR="#BC8F8F"><B>&quot;Eigensystem&quot;</FONT></B>);
  
  #ifdef HAVE_SLEPC
  
    <FONT COLOR="#228B22"><B>const</FONT></B> Mesh&amp; mesh = es.get_mesh();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = mesh.mesh_dimension();
  
    EigenSystem &amp; eigen_system = es.get_system&lt;EigenSystem&gt; (system_name);
  
    FEType fe_type = eigen_system.get_dof_map().variable_type(0);
  
    SparseMatrix&lt;Number&gt;&amp;  matrix_A = *eigen_system.matrix_A;
    SparseMatrix&lt;Number&gt;&amp;  matrix_B = *eigen_system.matrix_B;
  
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    QGauss qrule (dim, fe_type.default_quadrature_order());
  
    fe-&gt;attach_quadrature_rule (&amp;qrule);
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> DofMap&amp; dof_map = eigen_system.get_dof_map();
  
    DenseMatrix&lt;Number&gt;   Me;
    DenseMatrix&lt;Number&gt;   Ke;
  
    std::vector&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt; dof_indices;
  
  
  
    MeshBase::const_element_iterator       el     = mesh.elements_begin();
    <FONT COLOR="#228B22"><B>const</FONT></B> MeshBase::const_element_iterator end_el = mesh.elements_end();
   
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        <FONT COLOR="#228B22"><B>const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
  
        fe-&gt;reinit (elem);
  
        Ke.resize (dof_indices.size(), dof_indices.size());
        Me.resize (dof_indices.size(), dof_indices.size());
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi.size(); i++)
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j&lt;phi.size(); j++)
  	    {
  	      Me(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp];
  	      Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
  	    }
  
        matrix_A.add_matrix (Ke, dof_indices);
        matrix_B.add_matrix (Me, dof_indices);
  
      } <I><FONT COLOR="#B22222">// end of element loop
</FONT></I>  
  
  #endif <I><FONT COLOR="#B22222">// HAVE_SLEPC
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
make: Entering directory `/aerolab/benkirk/libmesh/examples/ex17'
***************************************************************
* Running Example  ./ex17
***************************************************************
 
ERROR: This example requires libMesh to be
compiled with SLEPc eigen solvers support!
 
***************************************************************
* Done Running Example  ./ex17
***************************************************************
make: Leaving directory `/aerolab/benkirk/libmesh/examples/ex17'
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
