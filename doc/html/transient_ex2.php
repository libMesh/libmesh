<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("transient_ex2",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file transient_ex2.C with comments: </h1> 
<div class = "comment">
<h1>Transient Example 2 - The Newmark System and the Wave Equation</h1>

<br><br>This is the eighth example program. It builds on
the previous example programs.  It introduces the
NewmarkSystem class.  In this example the wave equation
is solved using the time integration scheme provided
by the NewmarkSystem class.

<br><br>This example comes with a cylindrical mesh given in the
universal file pipe-mesh.unv.
The mesh contains HEX8 and PRISM6 elements.


<br><br>C++ include files that we need
</div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
        #include &lt;fstream&gt;
        #include &lt;algorithm&gt;
        #include &lt;stdio.h&gt;
        #include &lt;math.h&gt;
        
</pre>
</div>
<div class = "comment">
Basic include file needed for the mesh functionality.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh.h"
        #include "libmesh/serial_mesh.h"
        #include "libmesh/gmv_io.h"
        #include "libmesh/vtk_io.h"
        #include "libmesh/newmark_system.h"
        #include "libmesh/equation_systems.h"
        
</pre>
</div>
<div class = "comment">
Define the Finite Element object.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/fe.h"
        
</pre>
</div>
<div class = "comment">
Define Gauss quadrature rules.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/quadrature_gauss.h"
        
</pre>
</div>
<div class = "comment">
Define useful datatypes for finite element
</div>

<div class ="fragment">
<pre>
        #include "libmesh/dense_matrix.h"
        #include "libmesh/dense_vector.h"
        
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
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/numeric_vector.h"
        
</pre>
</div>
<div class = "comment">
Define the DofMap, which handles degree of freedom
indexing.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/dof_map.h"
        
</pre>
</div>
<div class = "comment">
The definition of a vertex associated with a Mesh.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/node.h"
        
</pre>
</div>
<div class = "comment">
The definition of a geometric element
</div>

<div class ="fragment">
<pre>
        #include "libmesh/elem.h"
        
</pre>
</div>
<div class = "comment">
Defines the MeshData class, which allows you to store
data about the mesh when reading in files, etc.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/mesh_data.h"
        
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
the linear system for our problem, governed by the linear
wave equation.
</div>

<div class ="fragment">
<pre>
        void assemble_wave(EquationSystems& es,
                           const std::string& system_name);
        
        
</pre>
</div>
<div class = "comment">
Function Prototype.  This function will be used to apply the
initial conditions.
</div>

<div class ="fragment">
<pre>
        void apply_initial(EquationSystems& es,
                           const std::string& system_name);
        
</pre>
</div>
<div class = "comment">
Function Prototype.  This function imposes
Dirichlet Boundary conditions via the penalty
method after the system is assembled.
</div>

<div class ="fragment">
<pre>
        void fill_dirichlet_bc(EquationSystems& es,
                               const std::string& system_name);
        
</pre>
</div>
<div class = "comment">
The main program
</div>

<div class ="fragment">
<pre>
        int main (int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Initialize libraries, like in example 2.
</div>

<div class ="fragment">
<pre>
          LibMeshInit init (argc, argv);
        
</pre>
</div>
<div class = "comment">
Check for proper usage.
</div>

<div class ="fragment">
<pre>
          if (argc &lt; 2)
            {
              if (libMesh::processor_id() == 0)
                std::cerr &lt;&lt; "Usage: " &lt;&lt; argv[0] &lt;&lt; " [meshfile]"
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
LasPack solvers don't work so well for this example
(not sure why), and Trilinos matrices don't work at all.
Print a warning to the user if PETSc is not in use.
</div>

<div class ="fragment">
<pre>
          if (libMesh::default_solver_package() == LASPACK_SOLVERS)
            {
              std::cout &lt;&lt; "WARNING! It appears you are using the\n"
                        &lt;&lt; "LasPack solvers.  This example may not converge\n"
                        &lt;&lt; "using LasPack, but should work OK with PETSc.\n"
                        &lt;&lt; "http://www.mcs.anl.gov/petsc/\n"
                        &lt;&lt; std::endl;
            }
          else if (libMesh::default_solver_package() == TRILINOS_SOLVERS)
            {
              std::cout &lt;&lt; "WARNING! It appears you are using the\n"
                        &lt;&lt; "Trilinos solvers.  The current libMesh Epetra\n"
                        &lt;&lt; "interface does not allow sparse matrix addition,\n"
                        &lt;&lt; "as is needed in this problem.  We recommend\n"
                        &lt;&lt; "using PETSc: http://www.mcs.anl.gov/petsc/\n"
                        &lt;&lt; std::endl;
              return 0;
            }
          
</pre>
</div>
<div class = "comment">
Get the name of the mesh file
from the command line.
</div>

<div class ="fragment">
<pre>
          std::string mesh_file = argv[1];
          std::cout &lt;&lt; "Mesh file is: " &lt;&lt; mesh_file &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Skip this 3D example if libMesh was compiled as 1D or 2D-only.
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(3 &lt;= LIBMESH_DIM, "3D support");
          
</pre>
</div>
<div class = "comment">
Create a mesh.
This example directly references all mesh nodes and is
incompatible with ParallelMesh use.


<br><br></div>

<div class ="fragment">
<pre>
          SerialMesh mesh;
          MeshData mesh_data(mesh);
          
</pre>
</div>
<div class = "comment">
Read the meshfile specified in the command line or
use the internal mesh generator to create a uniform
grid on an elongated cube.
</div>

<div class ="fragment">
<pre>
          mesh.read(mesh_file, &mesh_data);
           
</pre>
</div>
<div class = "comment">
mesh.build_cube (10, 10, 40,
-1., 1.,
-1., 1.,
0., 4.,
HEX8);


<br><br>Print information about the mesh to the screen.
</div>

<div class ="fragment">
<pre>
          mesh.print_info();
        
</pre>
</div>
<div class = "comment">
The node that should be monitored.
</div>

<div class ="fragment">
<pre>
          const unsigned int result_node = 274;
        
          
</pre>
</div>
<div class = "comment">
Time stepping issues

<br><br>Note that the total current time is stored as a parameter
in the \pEquationSystems object.

<br><br>the time step size
</div>

<div class ="fragment">
<pre>
          const Real delta_t = .0000625;
        
</pre>
</div>
<div class = "comment">
The number of time steps.
</div>

<div class ="fragment">
<pre>
          unsigned int n_time_steps = 300;
          
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
Declare the system and its variables.
Create a NewmarkSystem named "Wave"
</div>

<div class ="fragment">
<pre>
          equation_systems.add_system&lt;NewmarkSystem&gt; ("Wave");
        
</pre>
</div>
<div class = "comment">
Use a handy reference to this system
</div>

<div class ="fragment">
<pre>
          NewmarkSystem & t_system = equation_systems.get_system&lt;NewmarkSystem&gt; ("Wave");
          
</pre>
</div>
<div class = "comment">
Add the variable "p" to "Wave".   "p"
will be approximated using first-order approximation.
</div>

<div class ="fragment">
<pre>
          t_system.add_variable("p", FIRST);
        
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly
function and the initial condition function defined
below.
</div>

<div class ="fragment">
<pre>
          t_system.attach_assemble_function  (assemble_wave);
          t_system.attach_init_function      (apply_initial);
        
</pre>
</div>
<div class = "comment">
Set the time step size, and optionally the
Newmark parameters, so that \p NewmarkSystem can 
compute integration constants.  Here we simply use 
pass only the time step and use default values 
for \p alpha=.25  and \p delta=.5.
</div>

<div class ="fragment">
<pre>
          t_system.set_newmark_parameters(delta_t);
        
</pre>
</div>
<div class = "comment">
Set the speed of sound and fluid density
as \p EquationSystems parameter,
so that \p assemble_wave() can access it.
</div>

<div class ="fragment">
<pre>
          equation_systems.parameters.set&lt;Real&gt;("speed")          = 1000.;
          equation_systems.parameters.set&lt;Real&gt;("fluid density")  = 1000.;
        
</pre>
</div>
<div class = "comment">
Start time integration from t=0
</div>

<div class ="fragment">
<pre>
          t_system.time = 0.;
        
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
A file to store the results at certain nodes.
</div>

<div class ="fragment">
<pre>
          std::ofstream res_out("pressure_node.res");
        
</pre>
</div>
<div class = "comment">
get the dof_numbers for the nodes that
should be monitored.
</div>

<div class ="fragment">
<pre>
          const unsigned int res_node_no = result_node;
          const Node& res_node = mesh.node(res_node_no-1);
          unsigned int dof_no = res_node.dof_number(0,0,0);
        
</pre>
</div>
<div class = "comment">
Assemble the time independent system matrices and rhs.
This function will also compute the effective system matrix
K~=K+a_0*M+a_1*C and apply user specified initial
conditions. 
</div>

<div class ="fragment">
<pre>
          t_system.assemble();
        
</pre>
</div>
<div class = "comment">
Now solve for each time step.
For convenience, use a local buffer of the 
current time.  But once this time is updated,
also update the \p EquationSystems parameter
Start with t_time = 0 and write a short header
to the nodal result file
</div>

<div class ="fragment">
<pre>
          res_out &lt;&lt; "# pressure at node " &lt;&lt; res_node_no &lt;&lt; "\n"
                  &lt;&lt; "# time\tpressure\n"
                  &lt;&lt; t_system.time &lt;&lt; "\t" &lt;&lt; 0 &lt;&lt; std::endl;
        
        
          for (unsigned int time_step=0; time_step&lt;n_time_steps; time_step++)
            {
</pre>
</div>
<div class = "comment">
Update the time.  Both here and in the
\p EquationSystems object
</div>

<div class ="fragment">
<pre>
              t_system.time += delta_t;
        
</pre>
</div>
<div class = "comment">
Update the rhs.
</div>

<div class ="fragment">
<pre>
              t_system.update_rhs();
        
</pre>
</div>
<div class = "comment">
Impose essential boundary conditions.
Not that since the matrix is only assembled once,
the penalty parameter should be added to the matrix
only in the first time step.  The applied
boundary conditions may be time-dependent and hence
the rhs vector is considered in each time step. 
</div>

<div class ="fragment">
<pre>
              if (time_step == 0)
                {
</pre>
</div>
<div class = "comment">
The local function \p fill_dirichlet_bc()
may also set Dirichlet boundary conditions for the
matrix.  When you set the flag as shown below,
the flag will return true.  If you want it to return
false, simply do not set it.
</div>

<div class ="fragment">
<pre>
                  equation_systems.parameters.set&lt;bool&gt;("Newmark set BC for Matrix") = true;
        
                  fill_dirichlet_bc(equation_systems, "Wave");
        
</pre>
</div>
<div class = "comment">
unset the flag, so that it returns false
</div>

<div class ="fragment">
<pre>
                  equation_systems.parameters.set&lt;bool&gt;("Newmark set BC for Matrix") = false;
                }
              else
                fill_dirichlet_bc(equation_systems, "Wave");
        
</pre>
</div>
<div class = "comment">
Solve the system "Wave".
</div>

<div class ="fragment">
<pre>
              t_system.solve();
        
</pre>
</div>
<div class = "comment">
After solving the system, write the solution
to a GMV-formatted plot file.
Do only for a few time steps.
</div>

<div class ="fragment">
<pre>
              if (time_step == 30 || time_step == 60 ||
                  time_step == 90 || time_step == 120 )
                {
                  char buf[14];
        
        		  if (!libMesh::on_command_line("--vtk")){
        			  sprintf (buf, "out.%03d.gmv", time_step);
        
        	          GMVIO(mesh).write_equation_systems (buf,equation_systems);
        
        		  }else{
        #ifdef LIBMESH_HAVE_VTK
</pre>
</div>
<div class = "comment">
VTK viewers are generally not happy with two dots in a filename
</div>

<div class ="fragment">
<pre>
                                  sprintf (buf, "out_%03d.exd", time_step);
        
        	          VTKIO(mesh).write_equation_systems (buf,equation_systems);
        #endif // #ifdef LIBMESH_HAVE_VTK
        		  }
                }
        
</pre>
</div>
<div class = "comment">
Update the p, v and a.
</div>

<div class ="fragment">
<pre>
              t_system.update_u_v_a();
        
</pre>
</div>
<div class = "comment">
dof_no may not be local in parallel runs, so we may need a
global displacement vector
</div>

<div class ="fragment">
<pre>
              NumericVector&lt;Number&gt; &displacement
                = t_system.get_vector("displacement");
              std::vector&lt;Number&gt; global_displacement(displacement.size());
              displacement.localize(global_displacement);
        
</pre>
</div>
<div class = "comment">
Write nodal results to file.  The results can then
be viewed with e.g. gnuplot (run gnuplot and type
'plot "pressure_node.res" with lines' in the command line)
</div>

<div class ="fragment">
<pre>
              res_out &lt;&lt; t_system.time &lt;&lt; "\t"
                      &lt;&lt; global_displacement[dof_no]
                      &lt;&lt; std::endl;
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
        
</pre>
</div>
<div class = "comment">
This function assembles the system matrix and right-hand-side
for our wave equation.
</div>

<div class ="fragment">
<pre>
        void assemble_wave(EquationSystems& es,
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
          libmesh_assert_equal_to (system_name, "Wave");
        
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
Copy the speed of sound and fluid density
to a local variable.
</div>

<div class ="fragment">
<pre>
          const Real speed = es.parameters.get&lt;Real&gt;("speed");
          const Real rho   = es.parameters.get&lt;Real&gt;("fluid density");
        
</pre>
</div>
<div class = "comment">
Get a reference to our system, as before.
</div>

<div class ="fragment">
<pre>
          NewmarkSystem & t_system = es.get_system&lt;NewmarkSystem&gt; (system_name);
        
</pre>
</div>
<div class = "comment">
Get a constant reference to the Finite Element type
for the first (and only) variable in the system.
</div>

<div class ="fragment">
<pre>
          FEType fe_type = t_system.get_dof_map().variable_type(0);
        
</pre>
</div>
<div class = "comment">
In here, we will add the element matrices to the
@e additional matrices "stiffness_mass" and "damping"
and the additional vector "force", not to the members 
"matrix" and "rhs".  Therefore, get writable
references to them.
</div>

<div class ="fragment">
<pre>
          SparseMatrix&lt;Number&gt;&   stiffness = t_system.get_matrix("stiffness");
          SparseMatrix&lt;Number&gt;&   damping   = t_system.get_matrix("damping");
          SparseMatrix&lt;Number&gt;&   mass      = t_system.get_matrix("mass");
        
          NumericVector&lt;Number&gt;&  force     = t_system.get_vector("force");
        
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
          SparseMatrix&lt;Number&gt;&  matrix     = *t_system.matrix;
          DenseMatrix&lt;Number&gt;    zero_matrix;
        
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
A 2nd order Gauss quadrature rule for numerical integration.
</div>

<div class ="fragment">
<pre>
          QGauss qrule (dim, SECOND);
        
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
          const DofMap& dof_map = t_system.get_dof_map();
        
</pre>
</div>
<div class = "comment">
The element mass, damping and stiffness matrices
and the element contribution to the rhs.
</div>

<div class ="fragment">
<pre>
          DenseMatrix&lt;Number&gt;   Ke, Ce, Me;
          DenseVector&lt;Number&gt;   Fe;
        
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
Now we will loop over all the elements in the mesh.
We will compute the element matrix and right-hand-side
contribution.
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
element type is different (i.e. the last element was HEX8
and now have a PRISM6).
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
Now loop over the quadrature points.  This handles
the numeric integration.
</div>

<div class ="fragment">
<pre>
              for (unsigned int qp=0; qp&lt;qrule.n_points(); qp++)
                {
</pre>
</div>
<div class = "comment">
Now we will build the element matrix.  This involves
a double loop to integrate the test funcions (i) against
the trial functions (j).
</div>

<div class ="fragment">
<pre>
                  for (unsigned int i=0; i&lt;phi.size(); i++)
                    for (unsigned int j=0; j&lt;phi.size(); j++)
                      {
                        Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
                        Me(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp]
                                   *1./(speed*speed);
                      } // end of the matrix summation loop          
                } // end of quadrature point loop
        
</pre>
</div>
<div class = "comment">
Now compute the contribution to the element matrix and the
right-hand-side vector if the current element lies on the
boundary. 
</div>

<div class ="fragment">
<pre>
              {
</pre>
</div>
<div class = "comment">
In this example no natural boundary conditions will
be considered.  The code is left here so it can easily
be extended.

<br><br>don't do this for any side
</div>

<div class ="fragment">
<pre>
                for (unsigned int side=0; side&lt;elem-&gt;n_sides(); side++)
                  if (!true)
</pre>
</div>
<div class = "comment">
if (elem->neighbor(side) == NULL)
</div>

<div class ="fragment">
<pre>
                    {
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
                      const std::vector&lt;std::vector&lt;Real&gt; &gt;&  phi_face = fe_face-&gt;get_phi();
                      
</pre>
</div>
<div class = "comment">
The Jacobian * Quadrature Weight at the quadrature
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
Here we consider a normal acceleration acc_n=1 applied to
the whole boundary of our mesh.
</div>

<div class ="fragment">
<pre>
                      const Real acc_n_value = 1.0;
                      
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
Right-hand-side contribution due to prescribed
normal acceleration.
</div>

<div class ="fragment">
<pre>
                          for (unsigned int i=0; i&lt;phi_face.size(); i++)
                            {
                              Fe(i) += acc_n_value*rho
                                *phi_face[i][qp]*JxW_face[qp];
                            }
                        } // end face quadrature point loop          
                    } // end if (elem-&gt;neighbor(side) == NULL)
                
</pre>
</div>
<div class = "comment">
In this example the Dirichlet boundary conditions will be 
imposed via panalty method after the
system is assembled.
        

<br><br></div>

<div class ="fragment">
<pre>
              } // end boundary condition section          
        
</pre>
</div>
<div class = "comment">
If this assembly program were to be used on an adaptive mesh,
we would have to apply any hanging node constraint equations
by uncommenting the following lines:
std::vector<unsigned int> dof_indicesC = dof_indices;
std::vector<unsigned int> dof_indicesM = dof_indices;
dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
dof_map.constrain_element_matrix (Ce, dof_indicesC);
dof_map.constrain_element_matrix (Me, dof_indicesM);


<br><br>Finally, simply add the contributions to the additional
matrices and vector.
</div>

<div class ="fragment">
<pre>
              stiffness.add_matrix (Ke, dof_indices);
              damping.add_matrix   (Ce, dof_indices);
              mass.add_matrix      (Me, dof_indices);
              
              force.add_vector     (Fe, dof_indices);
            
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
            
            } // end of element loop
          
</pre>
</div>
<div class = "comment">
All done!
</div>

<div class ="fragment">
<pre>
          return;
        }
        
</pre>
</div>
<div class = "comment">
This function applies the initial conditions
</div>

<div class ="fragment">
<pre>
        void apply_initial(EquationSystems& es,
                           const std::string& system_name)
        {
</pre>
</div>
<div class = "comment">
Get a reference to our system, as before
</div>

<div class ="fragment">
<pre>
          NewmarkSystem & t_system = es.get_system&lt;NewmarkSystem&gt; (system_name);
          
</pre>
</div>
<div class = "comment">
Numeric vectors for the pressure, velocity and acceleration
values.
</div>

<div class ="fragment">
<pre>
          NumericVector&lt;Number&gt;&  pres_vec       = t_system.get_vector("displacement");
          NumericVector&lt;Number&gt;&  vel_vec        = t_system.get_vector("velocity");
          NumericVector&lt;Number&gt;&  acc_vec        = t_system.get_vector("acceleration");
        
</pre>
</div>
<div class = "comment">
Assume our fluid to be at rest, which would
also be the default conditions in class NewmarkSystem,
but let us do it explicetly here.
</div>

<div class ="fragment">
<pre>
          pres_vec.zero();
          vel_vec.zero();
          acc_vec.zero();
        }
        
</pre>
</div>
<div class = "comment">
This function applies the Dirichlet boundary conditions
</div>

<div class ="fragment">
<pre>
        void fill_dirichlet_bc(EquationSystems& es,
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
          libmesh_assert_equal_to (system_name, "Wave");
        
</pre>
</div>
<div class = "comment">
Get a reference to our system, as before.
</div>

<div class ="fragment">
<pre>
          NewmarkSystem & t_system = es.get_system&lt;NewmarkSystem&gt; (system_name);
        
</pre>
</div>
<div class = "comment">
Get writable references to the overall matrix and vector.
</div>

<div class ="fragment">
<pre>
          SparseMatrix&lt;Number&gt;&  matrix = *t_system.matrix;
          NumericVector&lt;Number&gt;& rhs    = *t_system.rhs;
        
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
Get \p libMesh's  \f$ \pi \f$
</div>

<div class ="fragment">
<pre>
          const Real pi = libMesh::pi;
        
</pre>
</div>
<div class = "comment">
Ask the \p EquationSystems flag whether
we should do this also for the matrix
</div>

<div class ="fragment">
<pre>
          const bool do_for_matrix =
            es.parameters.get&lt;bool&gt;("Newmark set BC for Matrix");
        
</pre>
</div>
<div class = "comment">
Number of nodes in the mesh.
</div>

<div class ="fragment">
<pre>
          unsigned int n_nodes = mesh.n_nodes();
        
          for (unsigned int n_cnt=0; n_cnt&lt;n_nodes; n_cnt++)
            {
</pre>
</div>
<div class = "comment">
Get a reference to the current node.
</div>

<div class ="fragment">
<pre>
              const Node& curr_node = mesh.node(n_cnt);
        
</pre>
</div>
<div class = "comment">
Check if Dirichlet BCs should be applied to this node.
Use the \p TOLERANCE from \p mesh_common.h as tolerance.
Here a pressure value is applied if the z-coord.
is equal to 4, which corresponds to one end of the
pipe-mesh in this directory.
</div>

<div class ="fragment">
<pre>
              const Real z_coo = 4.;
        
              if (fabs(curr_node(2)-z_coo) &lt; TOLERANCE)
                {
</pre>
</div>
<div class = "comment">
The global number of the respective degree of freedom.
</div>

<div class ="fragment">
<pre>
                  unsigned int dn = curr_node.dof_number(0,0,0);
        
</pre>
</div>
<div class = "comment">
The penalty parameter.
</div>

<div class ="fragment">
<pre>
                  const Real penalty = 1.e10;
        
</pre>
</div>
<div class = "comment">
Here we apply sinusoidal pressure values for 0<t<0.002
at one end of the pipe-mesh.
</div>

<div class ="fragment">
<pre>
                  Real p_value;
                  if (t_system.time &lt; .002 )
                    p_value = sin(2*pi*t_system.time/.002);
                  else
                    p_value = .0;
        
</pre>
</div>
<div class = "comment">
Now add the contributions to the matrix and the rhs.
</div>

<div class ="fragment">
<pre>
                  rhs.add(dn, p_value*penalty);
        
</pre>
</div>
<div class = "comment">
Add the panalty parameter to the global matrix
if desired.
</div>

<div class ="fragment">
<pre>
                  if (do_for_matrix)
                    matrix.add(dn, dn, penalty);
                }
            } // loop n_cnt
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file transient_ex2.C without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include &lt;fstream&gt;
  #include &lt;algorithm&gt;
  #include &lt;stdio.h&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/serial_mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/gmv_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/vtk_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/newmark_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/node.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_data.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_wave(EquationSystems&amp; es,
                     <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> apply_initial(EquationSystems&amp; es,
                     <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  <B><FONT COLOR="#228B22">void</FONT></B> fill_dirichlet_bc(EquationSystems&amp; es,
                         <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (argc &lt; 2)
      {
        <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Usage: &quot;</FONT></B> &lt;&lt; argv[0] &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; [meshfile]&quot;</FONT></B>
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
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::default_solver_package() == LASPACK_SOLVERS)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;WARNING! It appears you are using the\n&quot;</FONT></B>
                  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;LasPack solvers.  This example may not converge\n&quot;</FONT></B>
                  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;using LasPack, but should work OK with PETSc.\n&quot;</FONT></B>
                  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;http://www.mcs.anl.gov/petsc/\n&quot;</FONT></B>
                  &lt;&lt; std::endl;
      }
    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::default_solver_package() == TRILINOS_SOLVERS)
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;WARNING! It appears you are using the\n&quot;</FONT></B>
                  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Trilinos solvers.  The current libMesh Epetra\n&quot;</FONT></B>
                  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;interface does not allow sparse matrix addition,\n&quot;</FONT></B>
                  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;as is needed in this problem.  We recommend\n&quot;</FONT></B>
                  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;using PETSc: http://www.mcs.anl.gov/petsc/\n&quot;</FONT></B>
                  &lt;&lt; std::endl;
        <B><FONT COLOR="#A020F0">return</FONT></B> 0;
      }
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string mesh_file = argv[1];
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Mesh file is: &quot;</FONT></B> &lt;&lt; mesh_file &lt;&lt; std::endl;
  
    libmesh_example_assert(3 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;3D support&quot;</FONT></B>);
    
  
    SerialMesh mesh;
    MeshData mesh_data(mesh);
    
    mesh.read(mesh_file, &amp;mesh_data);
     
  
    mesh.print_info();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> result_node = 274;
  
    
    <B><FONT COLOR="#228B22">const</FONT></B> Real delta_t = .0000625;
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_time_steps = 300;
    
    EquationSystems equation_systems (mesh);
  
    equation_systems.add_system&lt;NewmarkSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Wave&quot;</FONT></B>);
  
    NewmarkSystem &amp; t_system = equation_systems.get_system&lt;NewmarkSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Wave&quot;</FONT></B>);
    
    t_system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>, FIRST);
  
    t_system.attach_assemble_function  (assemble_wave);
    t_system.attach_init_function      (apply_initial);
  
    t_system.set_newmark_parameters(delta_t);
  
    equation_systems.parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;speed&quot;</FONT></B>)          = 1000.;
    equation_systems.parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;fluid density&quot;</FONT></B>)  = 1000.;
  
    t_system.time = 0.;
  
    equation_systems.init();
  
    equation_systems.print_info();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream res_out(<B><FONT COLOR="#BC8F8F">&quot;pressure_node.res&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> res_node_no = result_node;
    <B><FONT COLOR="#228B22">const</FONT></B> Node&amp; res_node = mesh.node(res_node_no-1);
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dof_no = res_node.dof_number(0,0,0);
  
    t_system.assemble();
  
    res_out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;# pressure at node &quot;</FONT></B> &lt;&lt; res_node_no &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>
            &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;# time\tpressure\n&quot;</FONT></B>
            &lt;&lt; t_system.time &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\t&quot;</FONT></B> &lt;&lt; 0 &lt;&lt; std::endl;
  
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> time_step=0; time_step&lt;n_time_steps; time_step++)
      {
        t_system.time += delta_t;
  
        t_system.update_rhs();
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (time_step == 0)
          {
            equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">bool</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;Newmark set BC for Matrix&quot;</FONT></B>) = true;
  
            fill_dirichlet_bc(equation_systems, <B><FONT COLOR="#BC8F8F">&quot;Wave&quot;</FONT></B>);
  
            equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">bool</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;Newmark set BC for Matrix&quot;</FONT></B>) = false;
          }
        <B><FONT COLOR="#A020F0">else</FONT></B>
          fill_dirichlet_bc(equation_systems, <B><FONT COLOR="#BC8F8F">&quot;Wave&quot;</FONT></B>);
  
        t_system.solve();
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (time_step == 30 || time_step == 60 ||
            time_step == 90 || time_step == 120 )
          {
            <B><FONT COLOR="#228B22">char</FONT></B> buf[14];
  
  		  <B><FONT COLOR="#A020F0">if</FONT></B> (!libMesh::on_command_line(<B><FONT COLOR="#BC8F8F">&quot;--vtk&quot;</FONT></B>)){
  			  sprintf (buf, <B><FONT COLOR="#BC8F8F">&quot;out.%03d.gmv&quot;</FONT></B>, time_step);
  
  	          GMVIO(mesh).write_equation_systems (buf,equation_systems);
  
  		  }<B><FONT COLOR="#A020F0">else</FONT></B>{
  #ifdef LIBMESH_HAVE_VTK
  			  sprintf (buf, <B><FONT COLOR="#BC8F8F">&quot;out_%03d.exd&quot;</FONT></B>, time_step);
  
  	          VTKIO(mesh).write_equation_systems (buf,equation_systems);
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_VTK
</FONT></I>  		  }
          }
  
        t_system.update_u_v_a();
  
        NumericVector&lt;Number&gt; &amp;displacement
          = t_system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;displacement&quot;</FONT></B>);
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Number&gt; global_displacement(displacement.size());
        displacement.localize(global_displacement);
  
        res_out &lt;&lt; t_system.time &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\t&quot;</FONT></B>
                &lt;&lt; global_displacement[dof_no]
                &lt;&lt; std::endl;
      }
    
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_wave(EquationSystems&amp; es,
                     <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {  
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Wave&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real speed = es.parameters.get&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;speed&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> Real rho   = es.parameters.get&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;fluid density&quot;</FONT></B>);
  
    NewmarkSystem &amp; t_system = es.get_system&lt;NewmarkSystem&gt; (system_name);
  
    FEType fe_type = t_system.get_dof_map().variable_type(0);
  
    SparseMatrix&lt;Number&gt;&amp;   stiffness = t_system.get_matrix(<B><FONT COLOR="#BC8F8F">&quot;stiffness&quot;</FONT></B>);
    SparseMatrix&lt;Number&gt;&amp;   damping   = t_system.get_matrix(<B><FONT COLOR="#BC8F8F">&quot;damping&quot;</FONT></B>);
    SparseMatrix&lt;Number&gt;&amp;   mass      = t_system.get_matrix(<B><FONT COLOR="#BC8F8F">&quot;mass&quot;</FONT></B>);
  
    NumericVector&lt;Number&gt;&amp;  force     = t_system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;force&quot;</FONT></B>);
  
    SparseMatrix&lt;Number&gt;&amp;  matrix     = *t_system.matrix;
    DenseMatrix&lt;Number&gt;    zero_matrix;
  
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    QGauss qrule (dim, SECOND);
  
    fe-&gt;attach_quadrature_rule (&amp;qrule);
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = t_system.get_dof_map();
  
    DenseMatrix&lt;Number&gt;   Ke, Ce, Me;
    DenseVector&lt;Number&gt;   Fe;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
    
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
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
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
          {
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
                {
                  Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
                  Me(i,j) += JxW[qp]*phi[i][qp]*phi[j][qp]
                             *1./(speed*speed);
                } <I><FONT COLOR="#B22222">// end of the matrix summation loop          
</FONT></I>          } <I><FONT COLOR="#B22222">// end of quadrature point loop
</FONT></I>  
        {
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> side=0; side&lt;elem-&gt;n_sides(); side++)
            <B><FONT COLOR="#A020F0">if</FONT></B> (!true)
              {
                AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
                
                QGauss qface(dim-1, SECOND);
                
                fe_face-&gt;attach_quadrature_rule (&amp;qface);
                
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp;  phi_face = fe_face-&gt;get_phi();
                
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
                
                fe_face-&gt;reinit(elem, side);
  
                <B><FONT COLOR="#228B22">const</FONT></B> Real acc_n_value = 1.0;
                
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qface.n_points(); qp++)
                  {
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi_face.size(); i++)
                      {
                        Fe(i) += acc_n_value*rho
                          *phi_face[i][qp]*JxW_face[qp];
                      }
                  } <I><FONT COLOR="#B22222">// end face quadrature point loop          
</FONT></I>              } <I><FONT COLOR="#B22222">// end if (elem-&gt;neighbor(side) == NULL)
</FONT></I>          
          
        } <I><FONT COLOR="#B22222">// end boundary condition section          
</FONT></I>  
  
        stiffness.add_matrix (Ke, dof_indices);
        damping.add_matrix   (Ce, dof_indices);
        mass.add_matrix      (Me, dof_indices);
        
        force.add_vector     (Fe, dof_indices);
      
        matrix.add_matrix(zero_matrix, dof_indices);
      
      } <I><FONT COLOR="#B22222">// end of element loop
</FONT></I>    
    <B><FONT COLOR="#A020F0">return</FONT></B>;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> apply_initial(EquationSystems&amp; es,
                     <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    NewmarkSystem &amp; t_system = es.get_system&lt;NewmarkSystem&gt; (system_name);
    
    NumericVector&lt;Number&gt;&amp;  pres_vec       = t_system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;displacement&quot;</FONT></B>);
    NumericVector&lt;Number&gt;&amp;  vel_vec        = t_system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;velocity&quot;</FONT></B>);
    NumericVector&lt;Number&gt;&amp;  acc_vec        = t_system.get_vector(<B><FONT COLOR="#BC8F8F">&quot;acceleration&quot;</FONT></B>);
  
    pres_vec.zero();
    vel_vec.zero();
    acc_vec.zero();
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> fill_dirichlet_bc(EquationSystems&amp; es,
                         <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Wave&quot;</FONT></B>);
  
    NewmarkSystem &amp; t_system = es.get_system&lt;NewmarkSystem&gt; (system_name);
  
    SparseMatrix&lt;Number&gt;&amp;  matrix = *t_system.matrix;
    NumericVector&lt;Number&gt;&amp; rhs    = *t_system.rhs;
  
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real pi = libMesh::pi;
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> do_for_matrix =
      es.parameters.get&lt;<B><FONT COLOR="#228B22">bool</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;Newmark set BC for Matrix&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_nodes = mesh.n_nodes();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_cnt=0; n_cnt&lt;n_nodes; n_cnt++)
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Node&amp; curr_node = mesh.node(n_cnt);
  
        <B><FONT COLOR="#228B22">const</FONT></B> Real z_coo = 4.;
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (fabs(curr_node(2)-z_coo) &lt; TOLERANCE)
          {
            <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dn = curr_node.dof_number(0,0,0);
  
            <B><FONT COLOR="#228B22">const</FONT></B> Real penalty = 1.e10;
  
            Real p_value;
            <B><FONT COLOR="#A020F0">if</FONT></B> (t_system.time &lt; .002 )
              p_value = sin(2*pi*t_system.time/.002);
            <B><FONT COLOR="#A020F0">else</FONT></B>
              p_value = .0;
  
            rhs.add(dn, p_value*penalty);
  
            <B><FONT COLOR="#A020F0">if</FONT></B> (do_for_matrix)
              matrix.add(dn, dn, penalty);
          }
      } <I><FONT COLOR="#B22222">// loop n_cnt
</FONT></I>  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example transient_ex2:
*  mpirun -np 12 example-devel pipe-mesh.unv -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/transient/transient_ex2/.libs/lt-example-devel pipe-mesh.unv -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

Mesh file is: pipe-mesh.unv
*** Warning, This code is deprecated, and likely to be removed in future library versions! src/mesh/mesh_data.C, line 49, compiled Jan 31 2013 at 21:51:14 ***
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=3977
    n_local_nodes()=423
  n_elem()=3520
    n_local_elem()=286
    n_active_elem()=3520
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Wave"
    Type "Newmark"
    Variables="p" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=3977
    n_local_dofs()=423
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=9
    n_matrices()=4
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 21.5016
      Average Off-Processor Bandwidth <= 4.16344
      Maximum  On-Processor Bandwidth <= 36
      Maximum Off-Processor Bandwidth <= 27
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/transient/transient_ex2/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:21:46 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           9.369e+00      1.00063   9.367e+00
Objects:              1.604e+03      1.00000   1.604e+03
Flops:                2.470e+09      5.47389   1.221e+09  1.465e+10
Flops/sec:            2.637e+08      5.47164   1.303e+08  1.564e+09
MPI Messages:         3.472e+04      2.86319   2.200e+04  2.641e+05
MPI Message Lengths:  1.259e+07      1.99821   4.536e+02  1.198e+08
MPI Reductions:       1.896e+04      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 9.3674e+00 100.0%  1.4647e+10 100.0%  2.641e+05 100.0%  4.536e+02      100.0%  1.896e+04 100.0% 

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

VecMDot             4464 1.0 4.6387e-01 4.5 3.03e+07 2.0 0.0e+00 0.0e+00 4.5e+03  3  2  0  0 24   3  2  0  0 24   607
VecNorm             5064 1.0 3.4151e+0022.9 4.33e+06 2.0 0.0e+00 0.0e+00 5.1e+03 13  0  0  0 27  13  0  0  0 27    12
VecScale            5064 1.0 4.2763e-03 1.6 2.17e+06 2.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  4710
VecCopy             1500 1.0 2.5160e-03 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet              6292 1.0 5.3642e-03 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY             4200 1.0 1.7440e-02 1.1 3.34e+06 2.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1779
VecMAXPY            4764 1.0 2.0446e-02 1.8 3.42e+07 2.0 0.0e+00 0.0e+00 0.0e+00  0  2  0  0  0   0  2  0  0  0 15537
VecAssemblyBegin    1800 1.0 2.9381e-01 4.1 0.00e+00 0.0 1.3e+04 2.9e+02 5.4e+03  2  0  5  3 28   2  0  5  3 28     0
VecAssemblyEnd      1800 1.0 4.0882e-03 3.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     5664 1.0 3.0969e-02 1.7 0.00e+00 0.0 2.5e+05 4.6e+02 0.0e+00  0  0 95 95  0   0  0 95 95  0     0
VecScatterEnd       5664 1.0 4.4699e+00377.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 23  0  0  0  0  23  0  0  0  0     0
VecNormalize        4764 1.0 3.3697e+0026.7 6.12e+06 2.0 0.0e+00 0.0e+00 4.8e+03 12  0  0  0 25  12  0  0  0 25    17
MatMult             4764 1.0 4.5349e+0036.1 9.52e+07 2.0 2.1e+05 4.4e+02 0.0e+00 24  6 79 78  0  24  6 79 78  0   191
MatMultAdd           600 1.0 4.5132e-02 2.2 1.22e+07 2.0 2.6e+04 4.4e+02 0.0e+00  0  1 10 10  0   0  1 10 10  0  2474
MatSolve            5064 1.0 5.0835e-01 3.3 6.36e+08 4.0 0.0e+00 0.0e+00 0.0e+00  3 30  0  0  0   3 30  0  0  0  8584
MatLUFactorNum       300 1.0 1.8987e+00 7.5 1.66e+09 8.2 0.0e+00 0.0e+00 0.0e+00  9 59  0  0  0   9 59  0  0  0  4536
MatILUFactorSym      300 1.0 4.0512e+00 6.4 0.00e+00 0.0 0.0e+00 0.0e+00 9.0e+02 21  0  0  0  5  21  0  0  0  5     0
MatAssemblyBegin    1207 1.0 1.1069e-01 1.4 0.00e+00 0.0 3.3e+02 5.8e+03 2.4e+03  1  0  0  2 13   1  0  0  2 13     0
MatAssemblyEnd      1207 1.0 3.4465e-02 1.5 0.00e+00 0.0 6.2e+02 1.1e+02 5.6e+01  0  0  0  0  0   0  0  0  0  0     0
MatGetRow           2568 2.0 5.5647e-04 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ          300 1.0 9.4652e-05 2.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering       300 1.0 1.1343e-02 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 6.0e+02  0  0  0  0  3   0  0  0  0  3     0
MatZeroEntries        10 1.0 2.3627e-04 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatAXPY                3 1.0 4.9860e-03 1.0 0.00e+00 0.0 2.6e+02 1.1e+02 4.8e+01  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog      4464 1.0 4.7877e-01 3.8 6.07e+07 2.0 0.0e+00 0.0e+00 4.5e+03  3  4  0  0 24   3  4  0  0 24  1178
KSPSetUp             600 1.0 2.8026e-03 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve             300 1.0 7.0396e+00 1.0 2.45e+09 5.5 2.1e+05 4.4e+02 1.1e+04 75 99 79 78 58  75 99 79 78 58  2061
PCSetUp              600 1.0 6.0165e+00 6.4 1.66e+09 8.2 0.0e+00 0.0e+00 1.5e+03 30 59  0  0  8  30 59  0  0  8  1432
PCSetUpOnBlocks      300 1.0 6.0146e+00 6.4 1.66e+09 8.2 0.0e+00 0.0e+00 1.5e+03 30 59  0  0  8  30 59  0  0  8  1432
PCApply             5064 1.0 5.5878e-01 2.6 6.36e+08 4.0 0.0e+00 0.0e+00 0.0e+00  4 30  0  0  0   4 30  0  0  0  7809
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector   353            353      1680536     0
      Vector Scatter     8              8         8288     0
           Index Set   916            916      1197072     0
   IS L to G Mapping     1              1          564     0
              Matrix   321            321    226857488     0
       Krylov Solver     2              2        19360     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 1.19209e-07
Average time for MPI_Barrier(): 3.8147e-06
Average time for zero size MPI_Send(): 1.44839e-05
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
| Time:           Thu Jan 31 22:21:46 2013                                                                             |
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
| libMesh Performance: Alive time=9.50419, Active time=8.96498                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0284      0.028438    0.0640      0.064048    0.32     0.71     |
|   build_sparsity()                 1         0.0225      0.022546    0.0543      0.054350    0.25     0.61     |
|   create_dof_constraints()         1         0.0066      0.006610    0.0066      0.006610    0.07     0.07     |
|   distribute_dofs()                1         0.0593      0.059350    0.2180      0.218040    0.66     2.43     |
|   dof_indices()                    2105      0.1906      0.000091    0.1906      0.000091    2.13     2.13     |
|   prepare_send_list()              1         0.0006      0.000564    0.0006      0.000564    0.01     0.01     |
|   reinit()                         1         0.1542      0.154204    0.1542      0.154204    1.72     1.72     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          4         0.0177      0.004423    0.1292      0.032306    0.20     1.44     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        286       0.0158      0.000055    0.0158      0.000055    0.18     0.18     |
|   init_shape_functions()           5         0.0005      0.000108    0.0005      0.000108    0.01     0.01     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             232       0.0045      0.000019    0.0045      0.000019    0.05     0.05     |
|   compute_map()                    54        0.0025      0.000046    0.0025      0.000046    0.03     0.03     |
|   init_reference_to_physical_map() 5         0.0006      0.000111    0.0006      0.000111    0.01     0.01     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               4         0.1411      0.035276    0.1416      0.035388    1.57     1.58     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.1304      0.130368    0.2536      0.253587    1.45     2.83     |
|   read()                           1         0.0085      0.008512    0.0646      0.064638    0.09     0.72     |
|   renumber_nodes_and_elem()        2         0.0052      0.002583    0.0052      0.002583    0.06     0.06     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   broadcast()                      1         0.0202      0.020239    0.1373      0.137301    0.23     1.53     |
|   compute_hilbert_indices()        2         0.0709      0.035444    0.0709      0.035444    0.79     0.79     |
|   find_global_indices()            2         0.0235      0.011745    0.1147      0.057372    0.26     1.28     |
|   parallel_sort()                  2         0.0053      0.002655    0.0188      0.009399    0.06     0.21     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         4         0.0004      0.000100    0.2717      0.067936    0.00     3.03     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.2236      0.223584    0.2687      0.268712    2.49     3.00     |
|                                                                                                                |
| NewmarkSystem                                                                                                  |
|   initial_conditions ()            1         0.0000      0.000024    0.0000      0.000024    0.00     0.00     |
|   update_rhs ()                    300       0.1154      0.000385    0.1154      0.000385    1.29     1.29     |
|   update_u_v_a ()                  300       0.0490      0.000163    0.0490      0.000163    0.55     0.55     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      9         0.0024      0.000264    0.0024      0.000271    0.03     0.03     |
|   broadcast()                      6         0.0089      0.001480    0.0089      0.001480    0.10     0.10     |
|   max(bool)                        4         0.0000      0.000007    0.0000      0.000007    0.00     0.00     |
|   max(scalar)                      1689      0.0259      0.000015    0.0259      0.000015    0.29     0.29     |
|   max(vector)                      342       0.0071      0.000021    0.0228      0.000067    0.08     0.25     |
|   min(bool)                        2022      0.0260      0.000013    0.0260      0.000013    0.29     0.29     |
|   min(scalar)                      1682      0.2784      0.000165    0.2784      0.000165    3.11     3.11     |
|   min(vector)                      342       0.0078      0.000023    0.0259      0.000076    0.09     0.29     |
|   probe()                          132       0.0027      0.000021    0.0027      0.000021    0.03     0.03     |
|   receive()                        132       0.0010      0.000007    0.0037      0.000028    0.01     0.04     |
|   send()                           132       0.0004      0.000003    0.0004      0.000003    0.00     0.00     |
|   send_receive()                   136       0.0016      0.000011    0.0061      0.000045    0.02     0.07     |
|   sum()                            329       0.0570      0.000173    0.0756      0.000230    0.64     0.84     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           132       0.0003      0.000002    0.0003      0.000002    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0073      0.007343    0.0083      0.008266    0.08     0.09     |
|   set_parent_processor_ids()       1         0.0065      0.006466    0.0065      0.006466    0.07     0.07     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          300       7.1515      0.023838    7.1515      0.023838    79.77    79.77    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0268      0.026765    0.0770      0.077032    0.30     0.86     |
|                                                                                                                |
| UNVIO                                                                                                          |
|   count_elements()                 1         0.0040      0.003999    0.0040      0.003999    0.04     0.04     |
|   count_nodes()                    1         0.0057      0.005738    0.0057      0.005738    0.06     0.06     |
|   element_in()                     1         0.0356      0.035623    0.0356      0.035623    0.40     0.40     |
|   node_in()                        1         0.0108      0.010766    0.0108      0.010766    0.12     0.12     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            10714     8.9650                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example transient_ex2:
*  mpirun -np 12 example-devel pipe-mesh.unv -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
