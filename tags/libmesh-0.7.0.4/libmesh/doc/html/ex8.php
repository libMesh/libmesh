<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex8",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 8 - The Wave Equation</h1>

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
        #include "libmesh.h"
        #include "mesh.h"
        #include "gmv_io.h"
        #include "vtk_io.h"
        #include "newmark_system.h"
        #include "equation_systems.h"
        
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
The definition of a vertex associated with a Mesh.
</div>

<div class ="fragment">
<pre>
        #include "node.h"
        
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
Defines the MeshData class, which allows you to store
data about the mesh when reading in files, etc.
</div>

<div class ="fragment">
<pre>
        #include "mesh_data.h"
        
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
Initialize Petsc, like in example 2.
</div>

<div class ="fragment">
<pre>
          LibMeshInit init (argc, argv);
        
        #ifdef LIBMESH_ENABLE_PARMESH
          std::cout &lt;&lt; "ERROR: This example directly references\n"
                    &lt;&lt; "all mesh nodes and is incompatible with"
                    &lt;&lt; "ParallelMesh use."
                    &lt;&lt; std::endl;
          libmesh_example_assert(false, "--disable-parmesh");
        #else
        
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
                        &lt;&lt; "LasPack solvers.  ex8 may not converge\n"
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
</div>

<div class ="fragment">
<pre>
          Mesh mesh;
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
Store the current time as an
\p EquationSystems parameter, so that
\p fill_dirichlet_bc() can access it.
</div>

<div class ="fragment">
<pre>
          equation_systems.parameters.set&lt;Real&gt;("time")           = 0.;
        
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
          Real t_time = 0.;
          res_out &lt;&lt; "# pressure at node " &lt;&lt; res_node_no &lt;&lt; "\n"
                  &lt;&lt; "# time\tpressure\n"
                  &lt;&lt; t_time &lt;&lt; "\t" &lt;&lt; 0 &lt;&lt; std::endl;
        
        
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
              t_time += delta_t;
              equation_systems.parameters.set&lt;Real&gt;("time")  = t_time;
        
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
              res_out &lt;&lt; t_time &lt;&lt; "\t"
                      &lt;&lt; global_displacement[dof_no]
                      &lt;&lt; std::endl;
            }
        #endif
          
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
          libmesh_assert (system_name == "Wave");
        
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
          libmesh_assert (system_name == "Wave");
        
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
Ge the current time from \p EquationSystems
</div>

<div class ="fragment">
<pre>
          const Real current_time = es.parameters.get&lt;Real&gt;("time");
        
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
                  if (current_time &lt; .002 )
                    p_value = sin(2*pi*current_time/.002);
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
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include &lt;fstream&gt;
  #include &lt;algorithm&gt;
  #include &lt;stdio.h&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;gmv_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;vtk_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;newmark_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;node.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_data.h&quot;</FONT></B>
  
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
  
  #ifdef LIBMESH_ENABLE_PARMESH
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ERROR: This example directly references\n&quot;</FONT></B>
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;all mesh nodes and is incompatible with&quot;</FONT></B>
              &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ParallelMesh use.&quot;</FONT></B>
              &lt;&lt; std::endl;
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--disable-parmesh&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
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
                  &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;LasPack solvers.  ex8 may not converge\n&quot;</FONT></B>
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
    
    Mesh mesh;
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
  
    equation_systems.parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>)           = 0.;
  
    equation_systems.init();
  
    equation_systems.print_info();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::ofstream res_out(<B><FONT COLOR="#BC8F8F">&quot;pressure_node.res&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> res_node_no = result_node;
    <B><FONT COLOR="#228B22">const</FONT></B> Node&amp; res_node = mesh.node(res_node_no-1);
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dof_no = res_node.dof_number(0,0,0);
  
    t_system.assemble();
  
    Real t_time = 0.;
    res_out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;# pressure at node &quot;</FONT></B> &lt;&lt; res_node_no &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\n&quot;</FONT></B>
            &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;# time\tpressure\n&quot;</FONT></B>
            &lt;&lt; t_time &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\t&quot;</FONT></B> &lt;&lt; 0 &lt;&lt; std::endl;
  
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> time_step=0; time_step&lt;n_time_steps; time_step++)
      {
        t_time += delta_t;
        equation_systems.parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>)  = t_time;
  
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
  
        res_out &lt;&lt; t_time &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;\t&quot;</FONT></B>
                &lt;&lt; global_displacement[dof_no]
                &lt;&lt; std::endl;
      }
  #endif
    
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_wave(EquationSystems&amp; es,
                     <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {  
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Wave&quot;</FONT></B>);
  
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
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
  
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
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Wave&quot;</FONT></B>);
  
    NewmarkSystem &amp; t_system = es.get_system&lt;NewmarkSystem&gt; (system_name);
  
    SparseMatrix&lt;Number&gt;&amp;  matrix = *t_system.matrix;
    NumericVector&lt;Number&gt;&amp; rhs    = *t_system.rhs;
  
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real pi = libMesh::pi;
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> do_for_matrix =
      es.parameters.get&lt;<B><FONT COLOR="#228B22">bool</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;Newmark set BC for Matrix&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real current_time = es.parameters.get&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;time&quot;</FONT></B>);
  
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
            <B><FONT COLOR="#A020F0">if</FONT></B> (current_time &lt; .002 )
              p_value = sin(2*pi*current_time/.002);
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
* Running Example  mpirun -np 2 ./ex8-opt pipe-mesh.unv -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running ./ex8-opt pipe-mesh.unv -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

Mesh file is: pipe-mesh.unv
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=3977
    n_local_nodes()=2037
  n_elem()=3520
    n_local_elem()=1760
    n_active_elem()=3520
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Wave"
    Type "Newmark"
    Variables="p" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=3977
    n_local_dofs()=2037
    n_constrained_dofs()=0
    n_vectors()=9

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex8-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Thu Feb  3 12:14:11 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           3.713e+01      1.00001   3.712e+01
Objects:              1.859e+03      1.00000   1.859e+03
Flops:                1.652e+10      1.06440   1.602e+10  3.204e+10
Flops/sec:            4.450e+08      1.06441   4.315e+08  8.630e+08
MPI Messages:         2.748e+03      1.00000   2.748e+03  5.495e+03
MPI Message Lengths:  6.574e+06      1.03671   2.350e+03  1.291e+07
MPI Reductions:       9.761e+03      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 3.7125e+01 100.0%  3.2040e+10 100.0%  5.495e+03 100.0%  2.350e+03      100.0%  9.115e+03  93.4% 

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

VecMDot               30 1.0 4.5967e-04 1.9 1.22e+05 1.1 0.0e+00 0.0e+00 3.0e+01  0  0  0  0  0   0  0  0  0  0   519
VecNorm              630 1.0 2.8287e-02 2.6 2.57e+06 1.1 0.0e+00 0.0e+00 6.3e+02  0  0  0  0  6   0  0  0  0  7   177
VecScale             630 1.0 8.8358e-04 1.0 1.28e+06 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2836
VecCopy             1530 1.0 4.4994e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet              1011 1.0 1.5659e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY             3930 1.0 9.5201e-03 1.0 1.48e+07 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  3033
VecMAXPY              60 1.0 1.0920e-04 1.1 2.44e+05 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  4371
VecAssemblyBegin     900 1.0 6.6538e-03 1.0 0.00e+00 0.0 6.0e+02 5.8e+02 2.7e+03  0  0 11  3 28   0  0 11  3 30     0
VecAssemblyEnd       900 1.0 8.6761e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin     1530 1.0 6.2165e-03 1.3 0.00e+00 0.0 2.5e+03 8.7e+02 0.0e+00  0  0 45 17  0   0  0 45 17  0     0
VecScatterEnd       1530 1.0 2.1528e+00835.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  3  0  0  0  0   3  0  0  0  0     0
VecNormalize         330 1.0 2.4886e-02 3.2 2.02e+06 1.1 0.0e+00 0.0e+00 3.3e+02  0  0  0  0  3   0  0  0  0  4   158
MatMult              330 1.0 2.1744e+0087.1 3.08e+07 1.1 6.6e+02 7.8e+02 0.0e+00  3  0 12  4  0   3  0 12  4  0    28
MatMultAdd           600 1.0 5.0059e-02 1.0 5.72e+07 1.1 1.2e+03 7.8e+02 0.0e+00  0  0 22  7  0   0  0 22  7  0  2231
MatSolve              60 1.0 2.9930e-02 1.0 5.22e+07 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  3394
MatLUFactorNum       300 1.0 1.3236e+01 1.1 1.64e+10 1.1 0.0e+00 0.0e+00 0.0e+00 34 99  0  0  0  34 99  0  0  0  2397
MatILUFactorSym      300 1.0 2.3399e+01 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+02 61  0  0  0  3  61  0  0  0  3     0
MatAssemblyBegin    1207 1.0 1.0627e-02 1.1 0.00e+00 0.0 1.5e+01 1.1e+04 2.4e+03  0  0  0  1 25   0  0  0  1 26     0
MatAssemblyEnd      1207 1.0 7.0759e-02 1.0 0.00e+00 0.0 1.6e+01 2.0e+02 1.2e+03  0  0  0  0 13   0  0  0  0 14     0
MatGetRow           6111 1.1 7.6199e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetRowIJ          300 1.0 7.6056e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering       300 1.0 4.1714e-03 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 6.0e+02  0  0  0  0  6   0  0  0  0  7     0
MatZeroEntries        10 1.0 5.6219e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        30 1.0 5.7340e-04 1.6 2.44e+05 1.1 0.0e+00 0.0e+00 3.0e+01  0  0  0  0  0   0  0  0  0  0   832
KSPSetup             600 1.0 1.5330e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve             300 1.0 3.6726e+01 1.0 1.64e+10 1.1 6.6e+02 7.8e+02 1.6e+03 99100 12  4 16  99100 12  4 17   869
PCSetUp              600 1.0 3.6653e+01 1.1 1.64e+10 1.1 0.0e+00 0.0e+00 9.0e+02 96 99  0  0  9  96 99  0  0 10   866
PCSetUpOnBlocks      300 1.0 3.6651e+01 1.1 1.64e+10 1.1 0.0e+00 0.0e+00 9.0e+02 96 99  0  0  9  96 99  0  0 10   866
PCApply               60 1.0 3.0559e-02 1.0 5.22e+07 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  3324
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    27             27       300456     0
         Vec Scatter   305            305       264740     0
           Index Set  1210           1210      7964332     0
   IS L to G Mapping     1              1         8940     0
              Matrix   312            312   1573471792     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 1.38283e-06
Average time for zero size MPI_Send(): 6.55651e-06
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
| Time:           Thu Feb  3 12:14:11 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=37.1335, Active time=37.0319                                               |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 1         0.0008      0.000772    0.0008      0.000841    0.00     0.00     |
|   compute_sparsity()           1         0.0068      0.006839    0.0083      0.008270    0.02     0.02     |
|   create_dof_constraints()     1         0.0004      0.000388    0.0004      0.000388    0.00     0.00     |
|   distribute_dofs()            1         0.0021      0.002059    0.0057      0.005696    0.01     0.02     |
|   dof_indices()                12496     0.0045      0.000000    0.0045      0.000000    0.01     0.01     |
|   prepare_send_list()          1         0.0000      0.000021    0.0000      0.000021    0.00     0.00     |
|   reinit()                     1         0.0030      0.003030    0.0030      0.003030    0.01     0.01     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         1440      0.0016      0.000001    0.0016      0.000001    0.00     0.00     |
|   compute_map()                320       0.0007      0.000002    0.0007      0.000002    0.00     0.00     |
|   compute_shape_functions()    1760      0.0013      0.000001    0.0013      0.000001    0.00     0.00     |
|   init_shape_functions()       8         0.0002      0.000021    0.0002      0.000021    0.00     0.00     |
|                                                                                                            |
| GMVIO                                                                                                      |
|   write_nodal_data()           4         0.0482      0.012050    0.0482      0.012050    0.13     0.13     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             1         0.0069      0.006881    0.0092      0.009156    0.02     0.02     |
|   read()                       1         0.0001      0.000121    0.0134      0.013368    0.00     0.04     |
|   renumber_nodes_and_elem()    2         0.0003      0.000170    0.0003      0.000170    0.00     0.00     |
|                                                                                                            |
| MeshCommunication                                                                                          |
|   broadcast_bcs()              1         0.0000      0.000006    0.0000      0.000009    0.00     0.00     |
|   broadcast_mesh()             1         0.0013      0.001305    0.0015      0.001541    0.00     0.00     |
|   compute_hilbert_indices()    2         0.0134      0.006679    0.0134      0.006679    0.04     0.04     |
|   find_global_indices()        2         0.0016      0.000783    0.0157      0.007825    0.00     0.04     |
|   parallel_sort()              2         0.0005      0.000269    0.0006      0.000287    0.00     0.00     |
|                                                                                                            |
| MetisPartitioner                                                                                           |
|   partition()                  1         0.0047      0.004659    0.0123      0.012262    0.01     0.03     |
|                                                                                                            |
| NewmarkSystem                                                                                              |
|   initial_conditions ()        1         0.0000      0.000012    0.0000      0.000012    0.00     0.00     |
|   update_rhs ()                300       0.1060      0.000353    0.1060      0.000353    0.29     0.29     |
|   update_u_v_a ()              300       0.0075      0.000025    0.0075      0.000025    0.02     0.02     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  6         0.0006      0.000098    0.0006      0.000098    0.00     0.00     |
|   broadcast()                  6         0.0002      0.000037    0.0002      0.000037    0.00     0.00     |
|   gather()                     1         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   max(scalar)                  3         0.0024      0.000785    0.0024      0.000785    0.01     0.01     |
|   max(vector)                  2         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   min(vector)                  2         0.0000      0.000008    0.0000      0.000008    0.00     0.00     |
|   probe()                      10        0.0002      0.000022    0.0002      0.000022    0.00     0.00     |
|   receive()                    10        0.0000      0.000005    0.0003      0.000027    0.00     0.00     |
|   send()                       10        0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   send_receive()               14        0.0000      0.000003    0.0004      0.000025    0.00     0.00     |
|   sum()                        318       0.0078      0.000024    0.0078      0.000024    0.02     0.02     |
|   wait()                       10        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                            |
| Partitioner                                                                                                |
|   set_node_processor_ids()     1         0.0011      0.001088    0.0013      0.001268    0.00     0.00     |
|   set_parent_processor_ids()   1         0.0002      0.000155    0.0002      0.000155    0.00     0.00     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      300       36.7717     0.122572    36.7717     0.122572    99.30    99.30    |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   1         0.0225      0.022471    0.0276      0.027566    0.06     0.07     |
|                                                                                                            |
| UNVIO                                                                                                      |
|   count_elements()             1         0.0010      0.000961    0.0010      0.000961    0.00     0.00     |
|   count_nodes()                1         0.0007      0.000661    0.0007      0.000661    0.00     0.00     |
|   element_in()                 1         0.0063      0.006334    0.0063      0.006334    0.02     0.02     |
|   node_in()                    1         0.0053      0.005291    0.0053      0.005291    0.01     0.01     |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        17347     37.0319                                         100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 2 ./ex8-opt pipe-mesh.unv -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
