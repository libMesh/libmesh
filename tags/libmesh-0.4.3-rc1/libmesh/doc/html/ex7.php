<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex7",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 7 - Introduction to Complex Numbers and the "FrequencySystem"</h1>

<br><br>This is the seventh example program.  It builds on
the previous example programs, introduces complex
numbers and the FrequencySystem class to solve a 
simple Helmholtz equation grad(p)*grad(p)+(omega/c)^2*p=0,
for multiple frequencies rather efficiently.

<br><br>The FrequencySystem class offers two solution styles,
namely to solve large systems, or to solve
moderately-sized systems fast, for multiple frequencies.
The latter approach is implemented here.

<br><br>For this example the library has to be compiled with
complex numbers enabled. 
 

<br><br>C++ include files that we need
</div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
        #include &lt;algorithm&gt;
        #include &lt;stdio.h&gt;
        
</pre>
</div>
<div class = "comment">
Basic include files needed for overall functionality.
</div>

<div class ="fragment">
<pre>
        #include "libmesh.h"
        #include "libmesh_logging.h"
        #include "mesh.h"
        #include "gmv_io.h"
        #include "equation_systems.h"
        
</pre>
</div>
<div class = "comment">
Include FrequencySystem.  Compared to GeneralSystem,
this class offers added functionality for the solution of 
frequency-dependent systems.
</div>

<div class ="fragment">
<pre>
        #include "frequency_system.h"
        
</pre>
</div>
<div class = "comment">
Define the Finite Element object.
</div>

<div class ="fragment">
<pre>
        #include "fe.h"
        
</pre>
</div>
<div class = "comment">
Define Gauss quadrature rules.
</div>

<div class ="fragment">
<pre>
        #include "quadrature_gauss.h"
        
</pre>
</div>
<div class = "comment">
Define useful datatypes for finite element
matrix and vector components.
</div>

<div class ="fragment">
<pre>
        #include "dense_matrix.h"
        #include "dense_vector.h"
        
</pre>
</div>
<div class = "comment">
Define matrix and vector data types for the global 
equation system.  These are base classes,
from which specific implementations, like
the PETSc or LASPACK implementations, are derived.
</div>

<div class ="fragment">
<pre>
        #include "sparse_matrix.h"
        #include "numeric_vector.h"
        
</pre>
</div>
<div class = "comment">
Define the DofMap, which handles degree of freedom
indexing.
</div>

<div class ="fragment">
<pre>
        #include "dof_map.h"
        
</pre>
</div>
<div class = "comment">
Function prototype.  This is the function that will assemble
the mass, damping and stiffness matrices.  It will <i>not</i>
form an overall system matrix ready for solution.
</div>

<div class ="fragment">
<pre>
        void assemble_helmholtz(EquationSystems& es,
        			const std::string& system_name);
        
</pre>
</div>
<div class = "comment">
Function prototype.  This is the function that will combine
the previously-assembled mass, damping and stiffness matrices
to the overall matrix, which then renders ready for solution.
</div>

<div class ="fragment">
<pre>
        void add_M_C_K_helmholtz(EquationSystems& es,
        			 const std::string& system_name);
        
</pre>
</div>
<div class = "comment">
Begin the main program.  Note that this example only
works correctly if complex numbers have been enabled
in the library.  In order to link against the complex
PETSc libraries, you must have built PETSc with the same
C++ compiler that you used to build libMesh.  This is
so that the name mangling will be the same for the
routines in both libraries.
</div>

<div class ="fragment">
<pre>
        int main (int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Initialize Petsc, like in example 2.
</div>

<div class ="fragment">
<pre>
          libMesh::init (argc, argv);
          
</pre>
</div>
<div class = "comment">
This example is designed for complex numbers.   
</div>

<div class ="fragment">
<pre>
        #ifndef USE_COMPLEX_NUMBERS
        
          std::cerr &lt;&lt; "ERROR: This example is intended for " &lt;&lt; std::endl
        	    &lt;&lt; " use with complex numbers." &lt;&lt; std::endl;
          here();
        
          return 0;
        
        #else
          
</pre>
</div>
<div class = "comment">
Braces are used to force object scope, like in example 2
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
        	std::cerr &lt;&lt; "Usage: " &lt;&lt; argv[0] &lt;&lt; " -f [frequency]"
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
For now, restrict to dim=2, though this
may easily be changed, see example 4
</div>

<div class ="fragment">
<pre>
            const unsigned int dim = 2;
            
</pre>
</div>
<div class = "comment">
Get the frequency from argv[2] as a <i>float</i>,
currently, solve for 1/3rd, 2/3rd and 1/1th of the given frequency
</div>

<div class ="fragment">
<pre>
            const Real frequency_in = atof(argv[2]);
            const unsigned int n_frequencies = 3;
            
</pre>
</div>
<div class = "comment">
mesh discretization depends on frequency (badly guessed estimate...?)
</div>

<div class ="fragment">
<pre>
            const unsigned int n_el_per_dim =
              static_cast&lt;unsigned int&gt;(frequency_in*40.);
            
</pre>
</div>
<div class = "comment">
Tell the user the number of elements
</div>

<div class ="fragment">
<pre>
            std::cout &lt;&lt; " Using " &lt;&lt; n_el_per_dim &lt;&lt; " x " 
        	      &lt;&lt; n_el_per_dim &lt;&lt; " = " 
        	      &lt;&lt; n_el_per_dim*n_el_per_dim
        	      &lt;&lt; " QUAD9 elements"
        	      &lt;&lt; std::endl &lt;&lt; std::endl;
            
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
grid on the square [-1,1]^2.  We instruct the mesh generator
to build a mesh of n x n Quad9 elements.
</div>

<div class ="fragment">
<pre>
            mesh.build_square (n_el_per_dim, n_el_per_dim,
        		       -1., 1.,
        		       -1., 1.,
        		       QUAD9);
            
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
Create an equation systems object, which now handles
a frequency system, as opposed to previous examples.
</div>

<div class ="fragment">
<pre>
            EquationSystems equation_systems (mesh);
            
</pre>
</div>
<div class = "comment">
Create a FrequencySystem named "Helmholtz" & store a
reference to it.
</div>

<div class ="fragment">
<pre>
            FrequencySystem & f_system =      
              equation_systems.add_system&lt;FrequencySystem&gt; ("Helmholtz");
            
</pre>
</div>
<div class = "comment">
Add the variable "p" to "Helmholtz".  "p"
will be approximated using second-order approximation.
</div>

<div class ="fragment">
<pre>
            f_system.add_variable("p", SECOND);
            
</pre>
</div>
<div class = "comment">
Tell the frequency system about the two user-provided
functions.  In other circumstances, at least the
solve function has to be attached.
</div>

<div class ="fragment">
<pre>
            f_system.attach_assemble_function (assemble_helmholtz);
            f_system.attach_solve_function    (add_M_C_K_helmholtz);
            
</pre>
</div>
<div class = "comment">
To enable the fast solution scheme, additional
<i>global</i> matrices and one global vector, all appropriately sized,
have to be added.  The system object takes care of the
appropriate size, but the user should better fill explicitly
the sparsity structure of the overall matrix, so that the
fast matrix addition method can be used, as will be shown later.
</div>

<div class ="fragment">
<pre>
            f_system.add_matrix ("stiffness");
            f_system.add_matrix ("damping");
            f_system.add_matrix ("mass");
            f_system.add_vector ("rhs");
            
</pre>
</div>
<div class = "comment">
Communicates the frequencies to the system.  Note that
the frequency system stores the frequencies as parameters
in the equation systems object, so that our assemble and solve
functions may directly access them.
Will solve for 1/3rd, 2/3rd and 1/1th of the given frequency
</div>

<div class ="fragment">
<pre>
            f_system.set_frequencies_by_steps (frequency_in/n_frequencies,
        				       frequency_in,
        				       n_frequencies);
            
</pre>
</div>
<div class = "comment">
Use the parameters of the equation systems object to
tell the frequency system about the wave velocity and fluid
density.  The frequency system provides default values, but
these may be overridden, as shown here.
</div>

<div class ="fragment">
<pre>
            equation_systems.set_parameter ("wave speed") = 1.;
            equation_systems.set_parameter ("rho")        = 1.;
            
</pre>
</div>
<div class = "comment">
Initialize the data structures for the equation system.  <i>Always</i>
prior to this, the frequencies have to be communicated to the system.
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
            equation_systems.print_info ();
        
            for (unsigned int n=0; n &lt; n_frequencies; n++)
              {
</pre>
</div>
<div class = "comment">
Solve the system "Helmholtz" for the n-th frequency.  
Since we attached an assemble() function to the system,
the mass, damping and stiffness contributions will only
be assembled once.  Then, the system is solved for the
given frequencies.  Note that solve() may also solve 
the system only for specific frequencies.
</div>

<div class ="fragment">
<pre>
                f_system.solve (n,n);
        	
</pre>
</div>
<div class = "comment">
After solving the system, write the solution
to a GMV-formatted plot file, for every frequency.  
Now this is nice ;-) : we have the <i>identical</i> 
interface to the mesh write method as in the real-only 
case, but we output the real and imaginary 
part, and the magnitude, where the variable 
"p" is prepended with "r_", "i_", and "a_", 
respectively.
</div>

<div class ="fragment">
<pre>
                char buf[14];
        	sprintf (buf, "out%04d.gmv", n);
        	GMVIO(mesh).write_equation_systems (buf,
        					    equation_systems);
              }
            
</pre>
</div>
<div class = "comment">
Alternatively, the whole EquationSystems object can be
written to disk.  By default, the additional vectors are also
saved.
</div>

<div class ="fragment">
<pre>
            equation_systems.write ("eqn_sys.dat", libMeshEnums::WRITE);
          }
          
</pre>
</div>
<div class = "comment">
All done.  
</div>

<div class ="fragment">
<pre>
          return libMesh::close ();
        
        #endif 
        }
        
        
</pre>
</div>
<div class = "comment">
Here we define the matrix assembly routine for
the Helmholtz system.  This function will be
called to form the stiffness matrix and right-hand side.
</div>

<div class ="fragment">
<pre>
        void assemble_helmholtz(EquationSystems& es,
        			const std::string& system_name)
        {
        #ifdef USE_COMPLEX_NUMBERS
            
</pre>
</div>
<div class = "comment">
It is a good idea to make sure we are assembling
the proper system.
</div>

<div class ="fragment">
<pre>
          assert (system_name == "Helmholtz");
          
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
The dimension that we are in
</div>

<div class ="fragment">
<pre>
          const unsigned int dim = mesh.mesh_dimension();
          
</pre>
</div>
<div class = "comment">
Get a reference to our system, as before
</div>

<div class ="fragment">
<pre>
          FrequencySystem & f_system =
            es.get_system&lt;FrequencySystem&gt; (system_name);
          
</pre>
</div>
<div class = "comment">
A const reference to the DofMap object for this system.  The DofMap
object handles the index translation from node and element numbers
to degree of freedom numbers.
</div>

<div class ="fragment">
<pre>
          const DofMap& dof_map = f_system.get_dof_map();
          
</pre>
</div>
<div class = "comment">
Get a constant reference to the Finite Element type
for the first (and only) variable in the system.
</div>

<div class ="fragment">
<pre>
          const FEType fe_type = dof_map.variable_type(0);
        
</pre>
</div>
<div class = "comment">
For the admittance boundary condition,
get the fluid density
</div>

<div class ="fragment">
<pre>
          const Real rho = es.parameter("rho");
          
</pre>
</div>
<div class = "comment">
In here, we will add the element matrices to the
<i>additional</i> matrices "stiffness_mass", "damping",
and the additional vector "rhs", not to the members 
"matrix" and "rhs".  Therefore, get writable
references to them
</div>

<div class ="fragment">
<pre>
          SparseMatrix&lt;Number&gt;&   stiffness      = f_system.get_matrix("stiffness");
          SparseMatrix&lt;Number&gt;&   damping        = f_system.get_matrix("damping");
          SparseMatrix&lt;Number&gt;&   mass           = f_system.get_matrix("mass");
          NumericVector&lt;Number&gt;&  freq_indep_rhs = f_system.get_vector("rhs");
          
</pre>
</div>
<div class = "comment">
Some solver packages (PETSc) are especially picky about
allocating sparsity structure and truly assigning values
to this structure.  Namely, matrix additions, as performed
later, exhibit acceptable performance only for identical
sparsity structures.  Therefore, explicitly zero the
values in the collective matrix, so that matrix additions
encounter identical sparsity structures.
</div>

<div class ="fragment">
<pre>
          SparseMatrix&lt;Number&gt;&  matrix           = *f_system.matrix;
          
</pre>
</div>
<div class = "comment">
------------------------------------------------------------------
Finite Element related stuff

<br><br>Build a Finite Element object of the specified type.  Since the
FEBase::build() member dynamically creates memory we will
store the object as an AutoPtr<FEBase>.  This can be thought
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
The element Jacobian// quadrature weight at each integration point.   
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
Now this is slightly different from example 4.
We will not add directly to the overall (PETSc/LASPACK) matrix,
but to the additional matrices "stiffness_mass" and "damping".
The same holds for the right-hand-side vector Fe, which we will
later on store in the additional vector "rhs". 
The zero_matrix is used to explicitly induce the same sparsity
structure in the overall matrix.
see later on. (At least) the mass, and stiffness matrices, however, 
are inherently real.  Therefore, store these as one complex
matrix.  This will definitely save memory.
</div>

<div class ="fragment">
<pre>
          DenseMatrix&lt;Number&gt; Ke, Ce, Me, zero_matrix;
          DenseVector&lt;Number&gt; Fe;
          
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
We will compute the element matrix and right-hand-side
contribution.
</div>

<div class ="fragment">
<pre>
          const_local_elem_iterator           el (mesh.elements_begin());
          const const_local_elem_iterator end_el (mesh.elements_end());
          
          for ( ; el != end_el; ++el)
            {
</pre>
</div>
<div class = "comment">
Start logging the element initialization.
</div>

<div class ="fragment">
<pre>
              START_LOG("elem init","assemble_helmholtz");
              
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
Zero & resize the element matrix and right-hand side before
summing them, with different element types in the mesh this
is quite necessary.
</div>

<div class ="fragment">
<pre>
              {
                const unsigned int n_dof_indices = dof_indices.size();
        
        	Ke.resize          (n_dof_indices, n_dof_indices);
        	Ce.resize          (n_dof_indices, n_dof_indices);
        	Me.resize          (n_dof_indices, n_dof_indices);
        	zero_matrix.resize (n_dof_indices, n_dof_indices);
        	Fe.resize          (n_dof_indices);
              }
              
</pre>
</div>
<div class = "comment">
Stop logging the element initialization.
</div>

<div class ="fragment">
<pre>
              STOP_LOG("elem init","assemble_helmholtz");
        
</pre>
</div>
<div class = "comment">
Now loop over the quadrature points.  This handles
the numeric integration.
</div>

<div class ="fragment">
<pre>
              START_LOG("stiffness & mass","assemble_helmholtz");
        
              for (unsigned int qp=0; qp&lt;qrule.n_points(); qp++)
        	{
</pre>
</div>
<div class = "comment">
Now we will build the element matrix.  This involves
a double loop to integrate the test funcions (i) against
the trial functions (j).  Note the braces on the rhs
of Ke(i,j): these are quite necessary to finally compute
Real*(Point*Point) = Real, and not something else...
</div>

<div class ="fragment">
<pre>
                  for (unsigned int i=0; i&lt;phi.size(); i++)
        	    for (unsigned int j=0; j&lt;phi.size(); j++)
        	      {
        		Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
        		Me(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
        	      }	  
        	}
        
              STOP_LOG("stiffness & mass","assemble_helmholtz");
        
</pre>
</div>
<div class = "comment">
Now compute the contribution to the element matrix and the
right-hand-side vector if the current element lies on the
boundary. 

<br><br>The following loops over the sides of the element.
If the element has no neighbor on a side then that
side MUST live on a boundary of the domain.
	 

<br><br></div>

<div class ="fragment">
<pre>
              for (unsigned int side=0; side&lt;elem-&gt;n_sides(); side++)
        	if (elem-&gt;neighbor(side) == NULL)
        	  {
        	    START_LOG("damping & rhs","assemble_helmholtz");
        	      
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
                    QGauss qface(dim-1, SECOND);
        	      
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
The value of the shape functions at the quadrature
points.
</div>

<div class ="fragment">
<pre>
                    const std::vector&lt;std::vector&lt;Real&gt; &gt;&  phi_face =
        	      fe_face-&gt;get_phi();
        	      
</pre>
</div>
<div class = "comment">
The Jacobian// Quadrature Weight at the quadrature
points on the face.
</div>

<div class ="fragment">
<pre>
                    const std::vector&lt;Real&gt;& JxW_face = fe_face-&gt;get_JxW();
        	      
</pre>
</div>
<div class = "comment">
Compute the shape function values on the element
face.
</div>

<div class ="fragment">
<pre>
                    fe_face-&gt;reinit(elem, side);
        
</pre>
</div>
<div class = "comment">
Here we consider a normal velocity vn=1 applied to
the whole boundary of our mesh.
</div>

<div class ="fragment">
<pre>
                    const Real vn_value = 1.;
        	      
</pre>
</div>
<div class = "comment">
Consider a normal admittance an=1
at some parts of the bounfdary
</div>

<div class ="fragment">
<pre>
                    const Real an_value = 1.;
        	      
</pre>
</div>
<div class = "comment">
Loop over the face quagrature points for integration.
</div>

<div class ="fragment">
<pre>
                    for (unsigned int qp=0; qp&lt;qface.n_points(); qp++)
        	      {
</pre>
</div>
<div class = "comment">
Right-hand-side contribution due to prescribed
normal velocity.
</div>

<div class ="fragment">
<pre>
                        for (unsigned int i=0; i&lt;phi_face.size(); i++)
        		  Fe(i) += vn_value*phi_face[i][qp]*JxW_face[qp];
        		
</pre>
</div>
<div class = "comment">
Element matrix contributrion due to precribed
admittance boundary conditions.
</div>

<div class ="fragment">
<pre>
                        for (unsigned int i=0; i&lt;phi_face.size(); i++)
        		  for (unsigned int j=0; j&lt;phi_face.size(); j++)
        		    Ce(i,j) += rho*an_value*JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
        	      }
        
        	    STOP_LOG("damping & rhs","assemble_helmholtz");
        	  }
              
</pre>
</div>
<div class = "comment">
Finally, simply add the contributions to the additional
matrices and vector.
</div>

<div class ="fragment">
<pre>
              stiffness.add_matrix      (Ke, dof_indices);
              damping.add_matrix        (Ce, dof_indices);
              mass.add_matrix           (Me, dof_indices);
              freq_indep_rhs.add_vector (Fe, dof_indices);
              
</pre>
</div>
<div class = "comment">
For the overall matrix, explicitly zero the entries where
we added values in the other ones, so that we have 
identical sparsity footprints.
</div>

<div class ="fragment">
<pre>
              matrix.add_matrix(zero_matrix, dof_indices);
            } 
            
</pre>
</div>
<div class = "comment">
All done!
</div>

<div class ="fragment">
<pre>
        #endif
        }
        
        
</pre>
</div>
<div class = "comment">
We now define the function which will combine
the previously-assembled mass, stiffness, and
damping matrices into a single system matrix.
</div>

<div class ="fragment">
<pre>
        void add_M_C_K_helmholtz(EquationSystems& es,
        			 const std::string& system_name)
        {
        #ifdef USE_COMPLEX_NUMBERS
        
          START_LOG("init phase","add_M_C_K_helmholtz");
          
</pre>
</div>
<div class = "comment">
It is a good idea to make sure we are assembling
the proper system.
</div>

<div class ="fragment">
<pre>
          assert (system_name == "Helmholtz");
          
</pre>
</div>
<div class = "comment">
Get a reference to our system, as before
</div>

<div class ="fragment">
<pre>
          FrequencySystem & f_system =
            es.get_system&lt;FrequencySystem&gt; (system_name);
          
</pre>
</div>
<div class = "comment">
Get the frequency, fluid density, and speed of sound
for which we should currently solve
</div>

<div class ="fragment">
<pre>
          const Real frequency = es.parameter ("current frequency");
          const Real rho       = es.parameter ("rho");
          const Real speed     = es.parameter ("wave speed");
          
</pre>
</div>
<div class = "comment">
Compute angular frequency omega and wave number k
</div>

<div class ="fragment">
<pre>
          const Real omega = 2.0*libMesh::pi*frequency;
          const Real k     = omega / speed;
          
</pre>
</div>
<div class = "comment">
Get writable references to the overall matrix and vector, where the 
frequency-dependent system is to be collected
</div>

<div class ="fragment">
<pre>
          SparseMatrix&lt;Number&gt;&  matrix          = *f_system.matrix;
          NumericVector&lt;Number&gt;& rhs             = *f_system.rhs;
          
</pre>
</div>
<div class = "comment">
Get writable references to the frequency-independent matrices
and rhs, though we only need to extract values.  This write access
is necessary, since solver packages have to close the data structure 
before they can extract values for computation.
</div>

<div class ="fragment">
<pre>
          SparseMatrix&lt;Number&gt;&   stiffness      = f_system.get_matrix("stiffness");
          SparseMatrix&lt;Number&gt;&   damping        = f_system.get_matrix("damping");
          SparseMatrix&lt;Number&gt;&   mass           = f_system.get_matrix("mass");
          NumericVector&lt;Number&gt;&  freq_indep_rhs = f_system.get_vector("rhs");
          
</pre>
</div>
<div class = "comment">
form the scaling values for the coming matrix and vector axpy's
</div>

<div class ="fragment">
<pre>
          const Number scale_stiffness (  1., 0.   );
          const Number scale_damping   (  0., omega);
          const Number scale_mass      (-k*k, 0.   );
          const Number scale_rhs       (  0., -(rho*omega));
          
</pre>
</div>
<div class = "comment">
Now simply add the matrices together, store the result
in matrix and rhs.  Clear them first.
</div>

<div class ="fragment">
<pre>
          matrix.zero ();
          rhs.zero    ();
          
</pre>
</div>
<div class = "comment">
The matrices from which values are added to another matrix
have to be closed.  The add() method does take care of 
that, but let us do it explicitly.
</div>

<div class ="fragment">
<pre>
          stiffness.close ();
          damping.close   ();
          mass.close      ();
        
          STOP_LOG("init phase","add_M_C_K_helmholtz");
        
          START_LOG("global matrix & vector additions","add_M_C_K_helmholtz");
          
</pre>
</div>
<div class = "comment">
add the stiffness and mass with the proper frequency to the
overall system.  For this to work properly, matrix has
to be not only initialized, but filled with the identical
sparsity structure as the matrix added to it, otherwise
solver packages like PETSc crash.

<br><br>Note that we have to add the mass and stiffness contributions
one at a time; otherwise, the real part of matrix would
be fine, but the imaginary part cluttered with unwanted products.
</div>

<div class ="fragment">
<pre>
          matrix.add (scale_stiffness, stiffness);
          matrix.add (scale_mass,      mass);
          matrix.add (scale_damping,   damping);
          rhs.add    (scale_rhs,       freq_indep_rhs);
        
          STOP_LOG("global matrix & vector additions","add_M_C_K_helmholtz");
          
</pre>
</div>
<div class = "comment">
The "matrix" and "rhs" are now ready for solution   
</div>

<div class ="fragment">
<pre>
        #endif
        }
        
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
   
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;stdio.h&gt;
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;libmesh.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;libmesh_logging.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;mesh.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;gmv_io.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;equation_systems.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;frequency_system.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;fe.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;quadrature_gauss.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_matrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_vector.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;sparse_matrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;numeric_vector.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;dof_map.h&quot;</FONT></B>
  
  <FONT COLOR="#228B22"><B>void</FONT></B> assemble_helmholtz(EquationSystems&amp; es,
  			<FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name);
  
  <FONT COLOR="#228B22"><B>void</FONT></B> add_M_C_K_helmholtz(EquationSystems&amp; es,
  			 <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name);
  
  <FONT COLOR="#228B22"><B>int</FONT></B> main (<FONT COLOR="#228B22"><B>int</FONT></B> argc, <FONT COLOR="#228B22"><B>char</FONT></B>** argv)
  {
    libMesh::init (argc, argv);
    
  #ifndef USE_COMPLEX_NUMBERS
  
    std::cerr &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;ERROR: This example is intended for &quot;</FONT></B> &lt;&lt; std::endl
  	    &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; use with complex numbers.&quot;</FONT></B> &lt;&lt; std::endl;
    here();
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  
  #<B><FONT COLOR="#A020F0">else</FONT></B>
    
    {
      <B><FONT COLOR="#A020F0">if</FONT></B> (argc &lt; 3)
        {
  	std::cerr &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Usage: &quot;</FONT></B> &lt;&lt; argv[0] &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; -f [frequency]&quot;</FONT></B>
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
      
      <FONT COLOR="#228B22"><B>const</FONT></B> Real frequency_in = atof(argv[2]);
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_frequencies = 3;
      
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_el_per_dim =
        static_cast&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt;(frequency_in*40.);
      
      std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; Using &quot;</FONT></B> &lt;&lt; n_el_per_dim &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; x &quot;</FONT></B> 
  	      &lt;&lt; n_el_per_dim &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; = &quot;</FONT></B> 
  	      &lt;&lt; n_el_per_dim*n_el_per_dim
  	      &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; QUAD9 elements&quot;</FONT></B>
  	      &lt;&lt; std::endl &lt;&lt; std::endl;
      
      Mesh mesh (dim);
      
      mesh.build_square (n_el_per_dim, n_el_per_dim,
  		       -1., 1.,
  		       -1., 1.,
  		       QUAD9);
      
      mesh.print_info();
      
      EquationSystems equation_systems (mesh);
      
      FrequencySystem &amp; f_system =      
        equation_systems.add_system&lt;FrequencySystem&gt; (<FONT COLOR="#BC8F8F"><B>&quot;Helmholtz&quot;</FONT></B>);
      
      f_system.add_variable(<FONT COLOR="#BC8F8F"><B>&quot;p&quot;</FONT></B>, SECOND);
      
      f_system.attach_assemble_function (assemble_helmholtz);
      f_system.attach_solve_function    (add_M_C_K_helmholtz);
      
      f_system.add_matrix (<FONT COLOR="#BC8F8F"><B>&quot;stiffness&quot;</FONT></B>);
      f_system.add_matrix (<FONT COLOR="#BC8F8F"><B>&quot;damping&quot;</FONT></B>);
      f_system.add_matrix (<FONT COLOR="#BC8F8F"><B>&quot;mass&quot;</FONT></B>);
      f_system.add_vector (<FONT COLOR="#BC8F8F"><B>&quot;rhs&quot;</FONT></B>);
      
      f_system.set_frequencies_by_steps (frequency_in/n_frequencies,
  				       frequency_in,
  				       n_frequencies);
      
      equation_systems.set_parameter (<FONT COLOR="#BC8F8F"><B>&quot;wave speed&quot;</FONT></B>) = 1.;
      equation_systems.set_parameter (<FONT COLOR="#BC8F8F"><B>&quot;rho&quot;</FONT></B>)        = 1.;
      
      equation_systems.init ();
      
      equation_systems.print_info ();
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n=0; n &lt; n_frequencies; n++)
        {
  	f_system.solve (n,n);
  	
  	<FONT COLOR="#228B22"><B>char</FONT></B> buf[14];
  	sprintf (buf, <FONT COLOR="#BC8F8F"><B>&quot;out%04d.gmv&quot;</FONT></B>, n);
  	GMVIO(mesh).write_equation_systems (buf,
  					    equation_systems);
        }
      
      equation_systems.write (<FONT COLOR="#BC8F8F"><B>&quot;eqn_sys.dat&quot;</FONT></B>, libMeshEnums::WRITE);
    }
    
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close ();
  
  #endif 
  }
  
  
  <FONT COLOR="#228B22"><B>void</FONT></B> assemble_helmholtz(EquationSystems&amp; es,
  			<FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name)
  {
  #ifdef USE_COMPLEX_NUMBERS
      
    assert (system_name == <FONT COLOR="#BC8F8F"><B>&quot;Helmholtz&quot;</FONT></B>);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> Mesh&amp; mesh = es.get_mesh();
    
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = mesh.mesh_dimension();
    
    FrequencySystem &amp; f_system =
      es.get_system&lt;FrequencySystem&gt; (system_name);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> DofMap&amp; dof_map = f_system.get_dof_map();
    
    <FONT COLOR="#228B22"><B>const</FONT></B> FEType fe_type = dof_map.variable_type(0);
  
    <FONT COLOR="#228B22"><B>const</FONT></B> Real rho = es.parameter(<FONT COLOR="#BC8F8F"><B>&quot;rho&quot;</FONT></B>);
    
    SparseMatrix&lt;Number&gt;&amp;   stiffness      = f_system.get_matrix(<FONT COLOR="#BC8F8F"><B>&quot;stiffness&quot;</FONT></B>);
    SparseMatrix&lt;Number&gt;&amp;   damping        = f_system.get_matrix(<FONT COLOR="#BC8F8F"><B>&quot;damping&quot;</FONT></B>);
    SparseMatrix&lt;Number&gt;&amp;   mass           = f_system.get_matrix(<FONT COLOR="#BC8F8F"><B>&quot;mass&quot;</FONT></B>);
    NumericVector&lt;Number&gt;&amp;  freq_indep_rhs = f_system.get_vector(<FONT COLOR="#BC8F8F"><B>&quot;rhs&quot;</FONT></B>);
    
    SparseMatrix&lt;Number&gt;&amp;  matrix           = *f_system.matrix;
    
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    QGauss qrule (dim, FIFTH);
  
    fe-&gt;attach_quadrature_rule (&amp;qrule);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
    
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
    
    DenseMatrix&lt;Number&gt; Ke, Ce, Me, zero_matrix;
    DenseVector&lt;Number&gt; Fe;
    
    std::vector&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt; dof_indices;
  
    const_local_elem_iterator           el (mesh.elements_begin());
    <FONT COLOR="#228B22"><B>const</FONT></B> const_local_elem_iterator end_el (mesh.elements_end());
    
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        START_LOG(<FONT COLOR="#BC8F8F"><B>&quot;elem init&quot;</FONT></B>,<FONT COLOR="#BC8F8F"><B>&quot;assemble_helmholtz&quot;</FONT></B>);
        
        <FONT COLOR="#228B22"><B>const</FONT></B> Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
  
        fe-&gt;reinit (elem);
        
        {
          <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_dof_indices = dof_indices.size();
  
  	Ke.resize          (n_dof_indices, n_dof_indices);
  	Ce.resize          (n_dof_indices, n_dof_indices);
  	Me.resize          (n_dof_indices, n_dof_indices);
  	zero_matrix.resize (n_dof_indices, n_dof_indices);
  	Fe.resize          (n_dof_indices);
        }
        
        STOP_LOG(<FONT COLOR="#BC8F8F"><B>&quot;elem init&quot;</FONT></B>,<FONT COLOR="#BC8F8F"><B>&quot;assemble_helmholtz&quot;</FONT></B>);
  
        START_LOG(<FONT COLOR="#BC8F8F"><B>&quot;stiffness &amp; mass&quot;</FONT></B>,<FONT COLOR="#BC8F8F"><B>&quot;assemble_helmholtz&quot;</FONT></B>);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
  	{
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi.size(); i++)
  	    <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j&lt;phi.size(); j++)
  	      {
  		Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
  		Me(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
  	      }	  
  	}
  
        STOP_LOG(<FONT COLOR="#BC8F8F"><B>&quot;stiffness &amp; mass&quot;</FONT></B>,<FONT COLOR="#BC8F8F"><B>&quot;assemble_helmholtz&quot;</FONT></B>);
  
  	 
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> side=0; side&lt;elem-&gt;n_sides(); side++)
  	<B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(side) == NULL)
  	  {
  	    START_LOG(<FONT COLOR="#BC8F8F"><B>&quot;damping &amp; rhs&quot;</FONT></B>,<FONT COLOR="#BC8F8F"><B>&quot;assemble_helmholtz&quot;</FONT></B>);
  	      
  	    AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
  	      
  	    QGauss qface(dim-1, SECOND);
  	      
  	    fe_face-&gt;attach_quadrature_rule (&amp;qface);
  	      
  	    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp;  phi_face =
  	      fe_face-&gt;get_phi();
  	      
  	    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
  	      
  	    fe_face-&gt;reinit(elem, side);
  
  	    <FONT COLOR="#228B22"><B>const</FONT></B> Real vn_value = 1.;
  	      
  	    <FONT COLOR="#228B22"><B>const</FONT></B> Real an_value = 1.;
  	      
  	    <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp&lt;qface.n_points(); qp++)
  	      {
  		<B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi_face.size(); i++)
  		  Fe(i) += vn_value*phi_face[i][qp]*JxW_face[qp];
  		
  		<B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi_face.size(); i++)
  		  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j&lt;phi_face.size(); j++)
  		    Ce(i,j) += rho*an_value*JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
  	      }
  
  	    STOP_LOG(<FONT COLOR="#BC8F8F"><B>&quot;damping &amp; rhs&quot;</FONT></B>,<FONT COLOR="#BC8F8F"><B>&quot;assemble_helmholtz&quot;</FONT></B>);
  	  }
        
        stiffness.add_matrix      (Ke, dof_indices);
        damping.add_matrix        (Ce, dof_indices);
        mass.add_matrix           (Me, dof_indices);
        freq_indep_rhs.add_vector (Fe, dof_indices);
        
        matrix.add_matrix(zero_matrix, dof_indices);
      } 
      
  #endif
  }
  
  
  <FONT COLOR="#228B22"><B>void</FONT></B> add_M_C_K_helmholtz(EquationSystems&amp; es,
  			 <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name)
  {
  #ifdef USE_COMPLEX_NUMBERS
  
    START_LOG(<FONT COLOR="#BC8F8F"><B>&quot;init phase&quot;</FONT></B>,<FONT COLOR="#BC8F8F"><B>&quot;add_M_C_K_helmholtz&quot;</FONT></B>);
    
    assert (system_name == <FONT COLOR="#BC8F8F"><B>&quot;Helmholtz&quot;</FONT></B>);
    
    FrequencySystem &amp; f_system =
      es.get_system&lt;FrequencySystem&gt; (system_name);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> Real frequency = es.parameter (<FONT COLOR="#BC8F8F"><B>&quot;current frequency&quot;</FONT></B>);
    <FONT COLOR="#228B22"><B>const</FONT></B> Real rho       = es.parameter (<FONT COLOR="#BC8F8F"><B>&quot;rho&quot;</FONT></B>);
    <FONT COLOR="#228B22"><B>const</FONT></B> Real speed     = es.parameter (<FONT COLOR="#BC8F8F"><B>&quot;wave speed&quot;</FONT></B>);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> Real omega = 2.0*libMesh::pi*frequency;
    <FONT COLOR="#228B22"><B>const</FONT></B> Real k     = omega / speed;
    
    SparseMatrix&lt;Number&gt;&amp;  matrix          = *f_system.matrix;
    NumericVector&lt;Number&gt;&amp; rhs             = *f_system.rhs;
    
    SparseMatrix&lt;Number&gt;&amp;   stiffness      = f_system.get_matrix(<FONT COLOR="#BC8F8F"><B>&quot;stiffness&quot;</FONT></B>);
    SparseMatrix&lt;Number&gt;&amp;   damping        = f_system.get_matrix(<FONT COLOR="#BC8F8F"><B>&quot;damping&quot;</FONT></B>);
    SparseMatrix&lt;Number&gt;&amp;   mass           = f_system.get_matrix(<FONT COLOR="#BC8F8F"><B>&quot;mass&quot;</FONT></B>);
    NumericVector&lt;Number&gt;&amp;  freq_indep_rhs = f_system.get_vector(<FONT COLOR="#BC8F8F"><B>&quot;rhs&quot;</FONT></B>);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> Number scale_stiffness (  1., 0.   );
    <FONT COLOR="#228B22"><B>const</FONT></B> Number scale_damping   (  0., omega);
    <FONT COLOR="#228B22"><B>const</FONT></B> Number scale_mass      (-k*k, 0.   );
    <FONT COLOR="#228B22"><B>const</FONT></B> Number scale_rhs       (  0., -(rho*omega));
    
    matrix.zero ();
    rhs.zero    ();
    
    stiffness.close ();
    damping.close   ();
    mass.close      ();
  
    STOP_LOG(<FONT COLOR="#BC8F8F"><B>&quot;init phase&quot;</FONT></B>,<FONT COLOR="#BC8F8F"><B>&quot;add_M_C_K_helmholtz&quot;</FONT></B>);
  
    START_LOG(<FONT COLOR="#BC8F8F"><B>&quot;global matrix &amp; vector additions&quot;</FONT></B>,<FONT COLOR="#BC8F8F"><B>&quot;add_M_C_K_helmholtz&quot;</FONT></B>);
    
    matrix.add (scale_stiffness, stiffness);
    matrix.add (scale_mass,      mass);
    matrix.add (scale_damping,   damping);
    rhs.add    (scale_rhs,       freq_indep_rhs);
  
    STOP_LOG(<FONT COLOR="#BC8F8F"><B>&quot;global matrix &amp; vector additions&quot;</FONT></B>,<FONT COLOR="#BC8F8F"><B>&quot;add_M_C_K_helmholtz&quot;</FONT></B>);
    
  #endif
  }
  
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Linking ex7...
/home/peterson/code/libmesh/contrib/tecplot/lib/i686-pc-linux-gnu/tecio.a(tecxxx.o)(.text+0x1a7): In function `tecini':
: the use of `mktemp' is dangerous, better use `mkstemp'
***************************************************************
*** Skipping Example  ./ex7 , only good with --enable-complex
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
