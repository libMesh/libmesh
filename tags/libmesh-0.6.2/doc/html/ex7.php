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

<br><br>This example uses an L--shaped mesh and nodal boundary data
given in the files lshape.un and lshape_data.unv

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
        #include "mesh_generation.h"
        #include "gmv_io.h"
        #include "equation_systems.h"
        #include "elem.h"
        
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
Defines the MeshData class, which allows you to store
data about the mesh when reading in files, etc.
</div>

<div class ="fragment">
<pre>
        #include "mesh_data.h"
        
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
Create a dim-dimensional mesh.
</div>

<div class ="fragment">
<pre>
            Mesh mesh (dim);
        
</pre>
</div>
<div class = "comment">
Create a corresponding MeshData
and activate it. For more information on this object
cf. example 12.
</div>

<div class ="fragment">
<pre>
            MeshData mesh_data(mesh);
            mesh_data.activate();
        
</pre>
</div>
<div class = "comment">
Read the mesh file. Here the file lshape.unv contains
an L--shaped domain in .unv format.
</div>

<div class ="fragment">
<pre>
            mesh.read("lshape.unv", &mesh_data);
        
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
The load on the boundary of the domain is stored in
the .unv formated mesh data file lshape_data.unv.
At this, the data is given as complex valued normal
velocities.
</div>

<div class ="fragment">
<pre>
            mesh_data.read("lshape_data.unv");
        
</pre>
</div>
<div class = "comment">
Print information about the mesh to the screen.
</div>

<div class ="fragment">
<pre>
            mesh_data.print_info();
        
</pre>
</div>
<div class = "comment">
Create an equation systems object, which now handles
a frequency system, as opposed to previous examples.
Also pass a MeshData pointer so the data can be 
accessed in the matrix and rhs assembly.
</div>

<div class ="fragment">
<pre>
            EquationSystems equation_systems (mesh, &mesh_data);
            
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
            equation_systems.parameters.set&lt;Real&gt; ("wave speed") = 1.;
            equation_systems.parameters.set&lt;Real&gt; ("rho")        = 1.;
            
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
          const FEType& fe_type = dof_map.variable_type(0);
        
</pre>
</div>
<div class = "comment">
For the admittance boundary condition,
get the fluid density
</div>

<div class ="fragment">
<pre>
          const Real rho = es.parameters.get&lt;Real&gt;("rho");
          
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


<br><br></div>

<div class ="fragment">
<pre>
          MeshBase::const_element_iterator           el = mesh.local_elements_begin();
          const MeshBase::const_element_iterator end_el = mesh.local_elements_end();
          
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
Now compute the contribution to the element matrix
(due to mixed boundary conditions) if the current
element lies on the boundary. 

<br><br>The following loops over the sides of the element.
If the element has no neighbor on a side then that
side MUST live on a boundary of the domain.
	 

<br><br></div>

<div class ="fragment">
<pre>
              for (unsigned int side=0; side&lt;elem-&gt;n_sides(); side++)
        	if (elem-&gt;neighbor(side) == NULL)
        	  {
        	    START_LOG("damping","assemble_helmholtz");
        	      
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
For the Robin BCs consider a normal admittance an=1
at some parts of the bounfdary
</div>

<div class ="fragment">
<pre>
                    const Real an_value = 1.;
        	      
</pre>
</div>
<div class = "comment">
Loop over the face quadrature points for integration.
</div>

<div class ="fragment">
<pre>
                    for (unsigned int qp=0; qp&lt;qface.n_points(); qp++)
        	      {
        		
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
        
        	    STOP_LOG("damping","assemble_helmholtz");
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
        
            } // loop el
        
        
</pre>
</div>
<div class = "comment">
It now remains to compute the rhs. Here, we simply
get the normal velocities values on the boundary from the
mesh data.
</div>

<div class ="fragment">
<pre>
          {
            START_LOG("rhs","assemble_helmholtz");
            
</pre>
</div>
<div class = "comment">
get a reference to the mesh data.
</div>

<div class ="fragment">
<pre>
            const MeshData& mesh_data = es.get_mesh_data();
        
</pre>
</div>
<div class = "comment">
We will now loop over all nodes. In case nodal data 
for a certain node is available in the MeshData, we simply
adopt the corresponding value for the rhs vector.
Note that normal data was given in the mesh data file,
i.e. one value per node
</div>

<div class ="fragment">
<pre>
            assert(mesh_data.n_val_per_node() == 1);
        
            MeshBase::const_node_iterator       node_it  = mesh.nodes_begin();
            const MeshBase::const_node_iterator node_end = mesh.nodes_end();
            
            for ( ; node_it != node_end; ++node_it)
              {
</pre>
</div>
<div class = "comment">
the current node pointer
</div>

<div class ="fragment">
<pre>
                const Node* node = *node_it;
        
</pre>
</div>
<div class = "comment">
check if the mesh data has data for the current node
and do for all components
</div>

<div class ="fragment">
<pre>
                if (mesh_data.has_data(node))
        	  for (unsigned int comp=0; comp&lt;node-&gt;n_comp(0,0); comp++)
        	    {
</pre>
</div>
<div class = "comment">
the dof number
</div>

<div class ="fragment">
<pre>
                      unsigned int dn = node-&gt;dof_number(0,0,comp);
        
</pre>
</div>
<div class = "comment">
set the nodal value
</div>

<div class ="fragment">
<pre>
                      freq_indep_rhs.set(dn, mesh_data.get_data(node)[0]);
        	    }
              }
         
        
            STOP_LOG("rhs","assemble_helmholtz");
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
          const Real frequency = es.parameters.get&lt;Real&gt; ("current frequency");
          const Real rho       = es.parameters.get&lt;Real&gt; ("rho");
          const Real speed     = es.parameters.get&lt;Real&gt; ("wave speed");
          
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
          matrix.close(); matrix.zero ();
          rhs.close();    rhs.zero    ();
          
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
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh_logging.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;gmv_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;frequency_system.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_data.h&quot;</FONT></B>
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_helmholtz(EquationSystems&amp; es,
  			<B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  <B><FONT COLOR="#228B22">void</FONT></B> add_M_C_K_helmholtz(EquationSystems&amp; es,
  			 <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::init (argc, argv);
    
  #ifndef USE_COMPLEX_NUMBERS
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR: This example is intended for &quot;</FONT></B> &lt;&lt; std::endl
  	    &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; use with complex numbers.&quot;</FONT></B> &lt;&lt; std::endl;
    here();
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  
  #<B><FONT COLOR="#A020F0">else</FONT></B>
    
    {
      <B><FONT COLOR="#A020F0">if</FONT></B> (argc &lt; 3)
        {
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Usage: &quot;</FONT></B> &lt;&lt; argv[0] &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; -f [frequency]&quot;</FONT></B>
  		  &lt;&lt; std::endl;
  	
  	error();
        }
      
      <B><FONT COLOR="#A020F0">else</FONT></B> 
        {
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Running &quot;</FONT></B> &lt;&lt; argv[0];
  	
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">int</FONT></B> i=1; i&lt;argc; i++)
  	  <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B> &lt;&lt; argv[i];
  	
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl &lt;&lt; std::endl;
        }
      
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = 2;
      
      <B><FONT COLOR="#228B22">const</FONT></B> Real frequency_in = atof(argv[2]);
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_frequencies = 3;        
      
      Mesh mesh (dim);
  
      MeshData mesh_data(mesh);
      mesh_data.activate();
  
      mesh.read(<B><FONT COLOR="#BC8F8F">&quot;lshape.unv&quot;</FONT></B>, &amp;mesh_data);
  
      mesh.print_info();    
  
      mesh_data.read(<B><FONT COLOR="#BC8F8F">&quot;lshape_data.unv&quot;</FONT></B>);
  
      mesh_data.print_info();
  
      EquationSystems equation_systems (mesh, &amp;mesh_data);
      
      FrequencySystem &amp; f_system =      
        equation_systems.add_system&lt;FrequencySystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Helmholtz&quot;</FONT></B>);
      
      f_system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>, SECOND);
      
      f_system.attach_assemble_function (assemble_helmholtz);
      f_system.attach_solve_function    (add_M_C_K_helmholtz);
      
      f_system.add_matrix (<B><FONT COLOR="#BC8F8F">&quot;stiffness&quot;</FONT></B>);
      f_system.add_matrix (<B><FONT COLOR="#BC8F8F">&quot;damping&quot;</FONT></B>);
      f_system.add_matrix (<B><FONT COLOR="#BC8F8F">&quot;mass&quot;</FONT></B>);
      f_system.add_vector (<B><FONT COLOR="#BC8F8F">&quot;rhs&quot;</FONT></B>);
      
      f_system.set_frequencies_by_steps (frequency_in/n_frequencies,
  				       frequency_in,
  				       n_frequencies);
      
      equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;wave speed&quot;</FONT></B>) = 1.;
      equation_systems.parameters.set&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;rho&quot;</FONT></B>)        = 1.;
      
      equation_systems.init ();
      
      equation_systems.print_info ();
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n=0; n &lt; n_frequencies; n++)
        {
  	f_system.solve (n,n);
  	
  	<B><FONT COLOR="#228B22">char</FONT></B> buf[14];
  	sprintf (buf, <B><FONT COLOR="#BC8F8F">&quot;out%04d.gmv&quot;</FONT></B>, n);
  	GMVIO(mesh).write_equation_systems (buf,
  					    equation_systems);
        }
      
      equation_systems.write (<B><FONT COLOR="#BC8F8F">&quot;eqn_sys.dat&quot;</FONT></B>, libMeshEnums::WRITE);
    }
    
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close ();
  
  #endif 
  }
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_helmholtz(EquationSystems&amp; es,
  			<B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
  #ifdef USE_COMPLEX_NUMBERS
      
    assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Helmholtz&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> Mesh&amp; mesh = es.get_mesh();
    
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
    
    FrequencySystem &amp; f_system =
      es.get_system&lt;FrequencySystem&gt; (system_name);
    
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = f_system.get_dof_map();
    
    <B><FONT COLOR="#228B22">const</FONT></B> FEType&amp; fe_type = dof_map.variable_type(0);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real rho = es.parameters.get&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;rho&quot;</FONT></B>);
    
    SparseMatrix&lt;Number&gt;&amp;   stiffness      = f_system.get_matrix(<B><FONT COLOR="#BC8F8F">&quot;stiffness&quot;</FONT></B>);
    SparseMatrix&lt;Number&gt;&amp;   damping        = f_system.get_matrix(<B><FONT COLOR="#BC8F8F">&quot;damping&quot;</FONT></B>);
    SparseMatrix&lt;Number&gt;&amp;   mass           = f_system.get_matrix(<B><FONT COLOR="#BC8F8F">&quot;mass&quot;</FONT></B>);
    NumericVector&lt;Number&gt;&amp;  freq_indep_rhs = f_system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;rhs&quot;</FONT></B>);
    
    SparseMatrix&lt;Number&gt;&amp;  matrix           = *f_system.matrix;
    
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    QGauss qrule (dim, FIFTH);
  
    fe-&gt;attach_quadrature_rule (&amp;qrule);
    
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
    
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    DenseMatrix&lt;Number&gt; Ke, Ce, Me, zero_matrix;
    DenseVector&lt;Number&gt; Fe;
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
  
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator           el = mesh.local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.local_elements_end();
    
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        START_LOG(<B><FONT COLOR="#BC8F8F">&quot;elem init&quot;</FONT></B>,<B><FONT COLOR="#BC8F8F">&quot;assemble_helmholtz&quot;</FONT></B>);
        
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
  
        fe-&gt;reinit (elem);
        
        {
          <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_dof_indices = dof_indices.size();
  
  	Ke.resize          (n_dof_indices, n_dof_indices);
  	Ce.resize          (n_dof_indices, n_dof_indices);
  	Me.resize          (n_dof_indices, n_dof_indices);
  	zero_matrix.resize (n_dof_indices, n_dof_indices);
  	Fe.resize          (n_dof_indices);
        }
        
        STOP_LOG(<B><FONT COLOR="#BC8F8F">&quot;elem init&quot;</FONT></B>,<B><FONT COLOR="#BC8F8F">&quot;assemble_helmholtz&quot;</FONT></B>);
  
        START_LOG(<B><FONT COLOR="#BC8F8F">&quot;stiffness &amp; mass&quot;</FONT></B>,<B><FONT COLOR="#BC8F8F">&quot;assemble_helmholtz&quot;</FONT></B>);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
  	{
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
  	    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
  	      {
  		Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
  		Me(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp]);
  	      }	  
  	}
  
        STOP_LOG(<B><FONT COLOR="#BC8F8F">&quot;stiffness &amp; mass&quot;</FONT></B>,<B><FONT COLOR="#BC8F8F">&quot;assemble_helmholtz&quot;</FONT></B>);
  
  	 
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> side=0; side&lt;elem-&gt;n_sides(); side++)
  	<B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(side) == NULL)
  	  {
  	    START_LOG(<B><FONT COLOR="#BC8F8F">&quot;damping&quot;</FONT></B>,<B><FONT COLOR="#BC8F8F">&quot;assemble_helmholtz&quot;</FONT></B>);
  	      
  	    AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
  	      
  	    QGauss qface(dim-1, SECOND);
  	      
  	    fe_face-&gt;attach_quadrature_rule (&amp;qface);
  	      
  	    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp;  phi_face =
  	      fe_face-&gt;get_phi();
  	      
  	    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
  	      
  	    fe_face-&gt;reinit(elem, side);
  	      
  	    <B><FONT COLOR="#228B22">const</FONT></B> Real an_value = 1.;
  	      
  	    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qface.n_points(); qp++)
  	      {
  		
  		<B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi_face.size(); i++)
  		  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi_face.size(); j++)
  		    Ce(i,j) += rho*an_value*JxW_face[qp]*phi_face[i][qp]*phi_face[j][qp];
  	      }
  
  	    STOP_LOG(<B><FONT COLOR="#BC8F8F">&quot;damping&quot;</FONT></B>,<B><FONT COLOR="#BC8F8F">&quot;assemble_helmholtz&quot;</FONT></B>);
  	  }
  
       
        stiffness.add_matrix      (Ke, dof_indices);
        damping.add_matrix        (Ce, dof_indices);
        mass.add_matrix           (Me, dof_indices);
        
        matrix.add_matrix(zero_matrix, dof_indices);
  
      } <I><FONT COLOR="#B22222">// loop el
</FONT></I>  
  
    {
      START_LOG(<B><FONT COLOR="#BC8F8F">&quot;rhs&quot;</FONT></B>,<B><FONT COLOR="#BC8F8F">&quot;assemble_helmholtz&quot;</FONT></B>);
      
      <B><FONT COLOR="#228B22">const</FONT></B> MeshData&amp; mesh_data = es.get_mesh_data();
  
      assert(mesh_data.n_val_per_node() == 1);
  
      <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_node_iterator       node_it  = mesh.nodes_begin();
      <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_node_iterator node_end = mesh.nodes_end();
      
      <B><FONT COLOR="#A020F0">for</FONT></B> ( ; node_it != node_end; ++node_it)
        {
  	<B><FONT COLOR="#228B22">const</FONT></B> Node* node = *node_it;
  
  	<B><FONT COLOR="#A020F0">if</FONT></B> (mesh_data.has_data(node))
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> comp=0; comp&lt;node-&gt;n_comp(0,0); comp++)
  	    {
  	      <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dn = node-&gt;dof_number(0,0,comp);
  
  	      freq_indep_rhs.set(dn, mesh_data.get_data(node)[0]);
  	    }
        }
   
  
      STOP_LOG(<B><FONT COLOR="#BC8F8F">&quot;rhs&quot;</FONT></B>,<B><FONT COLOR="#BC8F8F">&quot;assemble_helmholtz&quot;</FONT></B>);
    }
  
  #endif
  }
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> add_M_C_K_helmholtz(EquationSystems&amp; es,
  			 <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
  #ifdef USE_COMPLEX_NUMBERS
  
    START_LOG(<B><FONT COLOR="#BC8F8F">&quot;init phase&quot;</FONT></B>,<B><FONT COLOR="#BC8F8F">&quot;add_M_C_K_helmholtz&quot;</FONT></B>);
    
    assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Helmholtz&quot;</FONT></B>);
    
    FrequencySystem &amp; f_system =
      es.get_system&lt;FrequencySystem&gt; (system_name);
    
    <B><FONT COLOR="#228B22">const</FONT></B> Real frequency = es.parameters.get&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;current frequency&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> Real rho       = es.parameters.get&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;rho&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> Real speed     = es.parameters.get&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;wave speed&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> Real omega = 2.0*libMesh::pi*frequency;
    <B><FONT COLOR="#228B22">const</FONT></B> Real k     = omega / speed;
    
    SparseMatrix&lt;Number&gt;&amp;  matrix          = *f_system.matrix;
    NumericVector&lt;Number&gt;&amp; rhs             = *f_system.rhs;
    
    SparseMatrix&lt;Number&gt;&amp;   stiffness      = f_system.get_matrix(<B><FONT COLOR="#BC8F8F">&quot;stiffness&quot;</FONT></B>);
    SparseMatrix&lt;Number&gt;&amp;   damping        = f_system.get_matrix(<B><FONT COLOR="#BC8F8F">&quot;damping&quot;</FONT></B>);
    SparseMatrix&lt;Number&gt;&amp;   mass           = f_system.get_matrix(<B><FONT COLOR="#BC8F8F">&quot;mass&quot;</FONT></B>);
    NumericVector&lt;Number&gt;&amp;  freq_indep_rhs = f_system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;rhs&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> Number scale_stiffness (  1., 0.   );
    <B><FONT COLOR="#228B22">const</FONT></B> Number scale_damping   (  0., omega);
    <B><FONT COLOR="#228B22">const</FONT></B> Number scale_mass      (-k*k, 0.   );
    <B><FONT COLOR="#228B22">const</FONT></B> Number scale_rhs       (  0., -(rho*omega));
    
    matrix.close(); matrix.zero ();
    rhs.close();    rhs.zero    ();
    
    stiffness.close ();
    damping.close   ();
    mass.close      ();
  
    STOP_LOG(<B><FONT COLOR="#BC8F8F">&quot;init phase&quot;</FONT></B>,<B><FONT COLOR="#BC8F8F">&quot;add_M_C_K_helmholtz&quot;</FONT></B>);
  
    START_LOG(<B><FONT COLOR="#BC8F8F">&quot;global matrix &amp; vector additions&quot;</FONT></B>,<B><FONT COLOR="#BC8F8F">&quot;add_M_C_K_helmholtz&quot;</FONT></B>);
    
    matrix.add (scale_stiffness, stiffness);
    matrix.add (scale_mass,      mass);
    matrix.add (scale_damping,   damping);
    rhs.add    (scale_rhs,       freq_indep_rhs);
  
    STOP_LOG(<B><FONT COLOR="#BC8F8F">&quot;global matrix &amp; vector additions&quot;</FONT></B>,<B><FONT COLOR="#BC8F8F">&quot;add_M_C_K_helmholtz&quot;</FONT></B>);
    
  #endif
  }
  
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
*** Skipping Example  ./ex7-devel , only good with --enable-complex
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
