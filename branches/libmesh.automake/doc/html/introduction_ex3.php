<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("introduction_ex3",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Introduction Example 3 - Solving a Poisson Problem</h1>

<br><br>This is the third example program.  It builds on
the second example program by showing how to solve a simple
Poisson system.  This example also introduces the notion
of customized matrix assembly functions, working with an
exact solution, and using element iterators.
We will not comment on things that
were already explained in the second example.


<br><br>C++ include files that we need
</div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
        #include &lt;algorithm&gt;
        #include &lt;math.h&gt;
        
</pre>
</div>
<div class = "comment">
Basic include files needed for the mesh functionality.
</div>

<div class ="fragment">
<pre>
        #include "libmesh.h"
        #include "mesh.h"
        #include "mesh_generation.h"
        #include "vtk_io.h"
        #include "linear_implicit_system.h"
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
matrix and vector components.
</div>

<div class ="fragment">
<pre>
        #include "sparse_matrix.h"
        #include "numeric_vector.h"
        #include "dense_matrix.h"
        #include "dense_vector.h"
        #include "elem.h"
        
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
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
</pre>
</div>
<div class = "comment">
Function prototype.  This is the function that will assemble
the linear system for our Poisson problem.  Note that the
function will take the  EquationSystems object and the
name of the system we are assembling as input.  From the
EquationSystems object we have access to the  Mesh and
other objects we might need.
</div>

<div class ="fragment">
<pre>
        void assemble_poisson(EquationSystems& es,
                              const std::string& system_name);
        
</pre>
</div>
<div class = "comment">
Function prototype for the exact solution.
</div>

<div class ="fragment">
<pre>
        Real exact_solution (const Real x,
                             const Real y,
                             const Real z = 0.);
        
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
Brief message to the user regarding the program name
and command line arguments.
</div>

<div class ="fragment">
<pre>
          std::cout &lt;&lt; "Running " &lt;&lt; argv[0];
          
          for (int i=1; i&lt;argc; i++)
            std::cout &lt;&lt; " " &lt;&lt; argv[i];
          
          std::cout &lt;&lt; std::endl &lt;&lt; std::endl;
          
</pre>
</div>
<div class = "comment">
Skip this 2D example if libMesh was compiled as 1D-only.
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(2 &lt;= LIBMESH_DIM, "2D support");
          Mesh mesh;
          
          
</pre>
</div>
<div class = "comment">
Use the MeshTools::Generation mesh generator to create a uniform
2D grid on the square [-1,1]^2.  We instruct the mesh generator
to build a mesh of 15x15 QUAD9 elements.  Building QUAD9
elements instead of the default QUAD4's we used in example 2
allow us to use higher-order approximation.
</div>

<div class ="fragment">
<pre>
          MeshTools::Generation::build_square (mesh, 
                                               15, 15,
                                               -1., 1.,
                                               -1., 1.,
                                               QUAD9);
        
</pre>
</div>
<div class = "comment">
Print information about the mesh to the screen.
Note that 5x5 QUAD9 elements actually has 11x11 nodes,
so this mesh is significantly larger than the one in example 2.
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
Declare the Poisson system and its variables.
The Poisson system is another example of a steady system.
</div>

<div class ="fragment">
<pre>
          equation_systems.add_system&lt;LinearImplicitSystem&gt; ("Poisson");
        
</pre>
</div>
<div class = "comment">
Adds the variable "u" to "Poisson".  "u"
will be approximated using second-order approximation.
</div>

<div class ="fragment">
<pre>
          equation_systems.get_system("Poisson").add_variable("u", SECOND);
        
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly
function.  This will be called when needed by the
library.
</div>

<div class ="fragment">
<pre>
          equation_systems.get_system("Poisson").attach_assemble_function (assemble_poisson);
          
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
Solve the system "Poisson".  Note that calling this
member will assemble the linear system and invoke
the default numerical solver.  With PETSc the solver can be
controlled from the command line.  For example,
you can invoke conjugate gradient with:

<br><br>./introduction_ex3 -ksp_type cg

<br><br>You can also get a nice X-window that monitors the solver
convergence with:

<br><br>./introduction-ex3 -ksp_xmonitor

<br><br>if you linked against the appropriate X libraries when you
built PETSc.
</div>

<div class ="fragment">
<pre>
          equation_systems.get_system("Poisson").solve();
        
        #if defined(LIBMESH_HAVE_VTK) && !defined(LIBMESH_ENABLE_PARMESH)
        
</pre>
</div>
<div class = "comment">
After solving the system write the solution
to a VTK-formatted plot file.
</div>

<div class ="fragment">
<pre>
          VTKIO (mesh).write_equation_systems ("out.pvtu", equation_systems);
        
        #endif // #ifdef LIBMESH_HAVE_VTK
        
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
We now define the matrix assembly function for the
Poisson system.  We need to first compute element
matrices and right-hand sides, and then take into
account the boundary conditions, which will be handled
via a penalty method.
</div>

<div class ="fragment">
<pre>
        void assemble_poisson(EquationSystems& es,
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
          libmesh_assert (system_name == "Poisson");
        
          
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
The dimension that we are running
</div>

<div class ="fragment">
<pre>
          const unsigned int dim = mesh.mesh_dimension();
        
</pre>
</div>
<div class = "comment">
Get a reference to the LinearImplicitSystem we are solving
</div>

<div class ="fragment">
<pre>
          LinearImplicitSystem& system = es.get_system&lt;LinearImplicitSystem&gt; ("Poisson");
        
</pre>
</div>
<div class = "comment">
A reference to the  DofMap object for this system.  The  DofMap
object handles the index translation from node and element numbers
to degree of freedom numbers.  We will talk more about the  DofMap
in future examples.
</div>

<div class ="fragment">
<pre>
          const DofMap& dof_map = system.get_dof_map();
          
</pre>
</div>
<div class = "comment">
Get a constant reference to the Finite Element type
for the first (and only) variable in the system.
</div>

<div class ="fragment">
<pre>
          FEType fe_type = dof_map.variable_type(0);
          
</pre>
</div>
<div class = "comment">
Build a Finite Element object of the specified type.  Since the
FEBase::build() member dynamically creates memory we will
store the object as an AutoPtr<FEBase>.  This can be thought
of as a pointer that will clean up after itself.  Introduction Example 4
describes some advantages of  AutoPtr's in the context of
quadrature rules.
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
          QGauss qface(dim-1, FIFTH);
          
</pre>
</div>
<div class = "comment">
Tell the finite element object to use our
quadrature rule.
</div>

<div class ="fragment">
<pre>
          fe_face-&gt;attach_quadrature_rule (&qface);
        
</pre>
</div>
<div class = "comment">
Here we define some references to cell-specific data that
will be used to assemble the linear system.

<br><br>The element Jacobian * quadrature weight at each integration point.   
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt;& JxW = fe-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The physical XY locations of the quadrature points on the element.
These might be useful for evaluating spatially varying material
properties at the quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Point&gt;& q_point = fe-&gt;get_xyz();
        
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
Define data structures to contain the element matrix
and right-hand-side vector contribution.  Following
basic finite element terminology we will denote these
"Ke" and "Fe".  These datatypes are templated on
Number, which allows the same code to work for real
or complex numbers.
</div>

<div class ="fragment">
<pre>
          DenseMatrix&lt;Number&gt; Ke;
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

<br><br>Element iterators are a nice way to iterate through all the
elements, or all the elements that have some property.  The
iterator el will iterate from the first to the last element on
the local processor.  The iterator end_el tells us when to stop.
It is smart to make this one const so that we don't accidentally
mess it up!  In case users later modify this program to include
refinement, we will be safe and will only consider the active
elements; hence we use a variant of the \p active_elem_iterator.
</div>

<div class ="fragment">
<pre>
          MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
          const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
         
</pre>
</div>
<div class = "comment">
Loop over the elements.  Note that  ++el is preferred to
el++ since the latter requires an unnecessary temporary
object.
</div>

<div class ="fragment">
<pre>
          for ( ; el != end_el ; ++el)
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
Zero the element matrix and right-hand side before
summing them.  We use the resize member here because
the number of degrees of freedom might have changed from
the last element.  Note that this will be the case if the
element type is different (i.e. the last element was a
triangle, now we are on a quadrilateral).


<br><br>The  DenseMatrix::resize() and the  DenseVector::resize()
members will automatically zero out the matrix  and vector.
</div>

<div class ="fragment">
<pre>
              Ke.resize (dof_indices.size(),
                         dof_indices.size());
        
              Fe.resize (dof_indices.size());
        
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
                      }
                  
</pre>
</div>
<div class = "comment">
This is the end of the matrix summation loop
Now we build the element right-hand-side contribution.
This involves a single loop in which we integrate the
"forcing function" in the PDE against the test functions.
</div>

<div class ="fragment">
<pre>
                  {
                    const Real x = q_point[qp](0);
                    const Real y = q_point[qp](1);
                    const Real eps = 1.e-3;
                    
        
</pre>
</div>
<div class = "comment">
"fxy" is the forcing function for the Poisson equation.
In this case we set fxy to be a finite difference
Laplacian approximation to the (known) exact solution.

<br><br>We will use the second-order accurate FD Laplacian
approximation, which in 2D is

<br><br>u_xx + u_yy = (u(i,j-1) + u(i,j+1) +
u(i-1,j) + u(i+1,j) +
-4*u(i,j))/h^2

<br><br>Since the value of the forcing function depends only
on the location of the quadrature point (q_point[qp])
we will compute it here, outside of the i-loop
</div>

<div class ="fragment">
<pre>
                    const Real fxy = -(exact_solution(x,y-eps) +
                                       exact_solution(x,y+eps) +
                                       exact_solution(x-eps,y) +
                                       exact_solution(x+eps,y) -
                                       4.*exact_solution(x,y))/eps/eps;
                    
                    for (unsigned int i=0; i&lt;phi.size(); i++)
                      Fe(i) += JxW[qp]*fxy*phi[i][qp];
                  } 
                } 
              
</pre>
</div>
<div class = "comment">
We have now reached the end of the RHS summation,
and the end of quadrature point loop, so
the interior element integration has
been completed.  However, we have not yet addressed
boundary conditions.  For this example we will only
consider simple Dirichlet boundary conditions.

<br><br>There are several ways Dirichlet boundary conditions
can be imposed.  A simple approach, which works for
interpolary bases like the standard Lagrange polynomials,
is to assign function values to the
degrees of freedom living on the domain boundary. This
works well for interpolary bases, but is more difficult
when non-interpolary (e.g Legendre or Hierarchic) bases
are used.

<br><br>Dirichlet boundary conditions can also be imposed with a
"penalty" method.  In this case essentially the L2 projection
of the boundary values are added to the matrix. The
projection is multiplied by some large factor so that, in
floating point arithmetic, the existing (smaller) entries
in the matrix and right-hand-side are effectively ignored.

<br><br>This amounts to adding a term of the form (in latex notation)

<br><br>\frac{1}{\epsilon} \int_{\delta \Omega} \phi_i \phi_j = \frac{1}{\epsilon} \int_{\delta \Omega} u \phi_i

<br><br>where

<br><br>\frac{1}{\epsilon} is the penalty parameter, defined such that \epsilon << 1
</div>

<div class ="fragment">
<pre>
              {
        
</pre>
</div>
<div class = "comment">
The following loop is over the sides of the element.
If the element has no neighbor on a side then that
side MUST live on a boundary of the domain.
</div>

<div class ="fragment">
<pre>
                for (unsigned int side=0; side&lt;elem-&gt;n_sides(); side++)
                  if (elem-&gt;neighbor(side) == NULL)
                    {
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
The XYZ locations (in physical space) of the
quadrature points on the face.  This is where
we will interpolate the boundary value function.
</div>

<div class ="fragment">
<pre>
                      const std::vector&lt;Point &gt;& qface_point = fe_face-&gt;get_xyz();
                      
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
Loop over the face quadrature points for integration.
</div>

<div class ="fragment">
<pre>
                      for (unsigned int qp=0; qp&lt;qface.n_points(); qp++)
                        {
        
</pre>
</div>
<div class = "comment">
The location on the boundary of the current
face quadrature point.
</div>

<div class ="fragment">
<pre>
                          const Real xf = qface_point[qp](0);
                          const Real yf = qface_point[qp](1);
        
</pre>
</div>
<div class = "comment">
The penalty value.  \frac{1}{\epsilon}
in the discussion above.
</div>

<div class ="fragment">
<pre>
                          const Real penalty = 1.e10;
        
</pre>
</div>
<div class = "comment">
The boundary value.
</div>

<div class ="fragment">
<pre>
                          const Real value = exact_solution(xf, yf);
                          
</pre>
</div>
<div class = "comment">
Matrix contribution of the L2 projection. 
</div>

<div class ="fragment">
<pre>
                          for (unsigned int i=0; i&lt;phi_face.size(); i++)
                            for (unsigned int j=0; j&lt;phi_face.size(); j++)
                              Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];
        
</pre>
</div>
<div class = "comment">
Right-hand-side contribution of the L2
projection.
</div>

<div class ="fragment">
<pre>
                          for (unsigned int i=0; i&lt;phi_face.size(); i++)
                            Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
                        } 
                    }
              }
              
</pre>
</div>
<div class = "comment">
We have now finished the quadrature point loop,
and have therefore applied all the boundary conditions.


<br><br>If this assembly program were to be used on an adaptive mesh,
we would have to apply any hanging node constraint equations
</div>

<div class ="fragment">
<pre>
              dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        
</pre>
</div>
<div class = "comment">
The element matrix and right-hand-side are now built
for this element.  Add them to the global matrix and
right-hand-side vector.  The  SparseMatrix::add_matrix()
and  NumericVector::add_vector() members do this for us.
</div>

<div class ="fragment">
<pre>
              system.matrix-&gt;add_matrix (Ke, dof_indices);
              system.rhs-&gt;add_vector    (Fe, dof_indices);
            }
          
</pre>
</div>
<div class = "comment">
All done!
</div>

<div class ="fragment">
<pre>
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;vtk_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_poisson(EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  Real exact_solution (<B><FONT COLOR="#228B22">const</FONT></B> Real x,
                       <B><FONT COLOR="#228B22">const</FONT></B> Real y,
                       <B><FONT COLOR="#228B22">const</FONT></B> Real z = 0.);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Running &quot;</FONT></B> &lt;&lt; argv[0];
    
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">int</FONT></B> i=1; i&lt;argc; i++)
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B> &lt;&lt; argv[i];
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl &lt;&lt; std::endl;
    
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
    Mesh mesh;
    
    
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square (mesh, 
                                         15, 15,
                                         -1., 1.,
                                         -1., 1.,
                                         QUAD9);
  
    mesh.print_info();
    
    EquationSystems equation_systems (mesh);
    
    equation_systems.add_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>);
  
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>).add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, SECOND);
  
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>).attach_assemble_function (assemble_poisson);
    
    equation_systems.init();
    
    equation_systems.print_info();
  
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>).solve();
  
  #<B><FONT COLOR="#A020F0">if</FONT></B> defined(LIBMESH_HAVE_VTK) &amp;&amp; !defined(LIBMESH_ENABLE_PARMESH)
  
    VTKIO (mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;out.pvtu&quot;</FONT></B>, equation_systems);
  
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_VTK
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_poisson(EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>);
  
    
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    LinearImplicitSystem&amp; system = es.get_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
    
    FEType fe_type = dof_map.variable_type(0);
    
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    QGauss qrule (dim, FIFTH);
    
    fe-&gt;attach_quadrature_rule (&amp;qrule);
    
    AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
    
    QGauss qface(dim-1, FIFTH);
    
    fe_face-&gt;attach_quadrature_rule (&amp;qface);
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; q_point = fe-&gt;get_xyz();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    DenseMatrix&lt;Number&gt; Ke;
    DenseVector&lt;Number&gt; Fe;
  
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
   
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el ; ++el)
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
  
        fe-&gt;reinit (elem);
  
  
  
        Ke.resize (dof_indices.size(),
                   dof_indices.size());
  
        Fe.resize (dof_indices.size());
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
          {
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
                {
                  Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
                }
            
            {
              <B><FONT COLOR="#228B22">const</FONT></B> Real x = q_point[qp](0);
              <B><FONT COLOR="#228B22">const</FONT></B> Real y = q_point[qp](1);
              <B><FONT COLOR="#228B22">const</FONT></B> Real eps = 1.e-3;
              
  
              <B><FONT COLOR="#228B22">const</FONT></B> Real fxy = -(exact_solution(x,y-eps) +
                                 exact_solution(x,y+eps) +
                                 exact_solution(x-eps,y) +
                                 exact_solution(x+eps,y) -
                                 4.*exact_solution(x,y))/eps/eps;
              
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
                Fe(i) += JxW[qp]*fxy*phi[i][qp];
            } 
          } 
        
        {
  
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> side=0; side&lt;elem-&gt;n_sides(); side++)
            <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(side) == NULL)
              {
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp;  phi_face = fe_face-&gt;get_phi();
                
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
                
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point &gt;&amp; qface_point = fe_face-&gt;get_xyz();
                
                fe_face-&gt;reinit(elem, side);
                
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qface.n_points(); qp++)
                  {
  
                    <B><FONT COLOR="#228B22">const</FONT></B> Real xf = qface_point[qp](0);
                    <B><FONT COLOR="#228B22">const</FONT></B> Real yf = qface_point[qp](1);
  
                    <B><FONT COLOR="#228B22">const</FONT></B> Real penalty = 1.e10;
  
                    <B><FONT COLOR="#228B22">const</FONT></B> Real value = exact_solution(xf, yf);
                    
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi_face.size(); i++)
                      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi_face.size(); j++)
                        Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];
  
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi_face.size(); i++)
                      Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
                  } 
              }
        }
        
  
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
  
        system.matrix-&gt;add_matrix (Ke, dof_indices);
        system.rhs-&gt;add_vector    (Fe, dof_indices);
      }
    
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Linking introduction_ex3-opt...
***************************************************************
* Running Example  mpirun -np 6 ./introduction_ex3-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running ./introduction_ex3-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=961
    n_local_nodes()=175
  n_elem()=225
    n_local_elem()=37
    n_active_elem()=225
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Poisson"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" 
    n_dofs()=961
    n_local_dofs()=175
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 14.1771
      Average Off-Processor Bandwidth <= 1.28
      Maximum  On-Processor Bandwidth <= 25
      Maximum Off-Processor Bandwidth <= 16
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

*** Warning, This code is untested, experimental, or likely to see future API changes: /h2/roystgnr/libmesh/svn/include/mesh/vtk_io.h, line 154, compiled Aug 24 2012 at 15:12:28 ***
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./introduction_ex3-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:21:35 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           8.568e-02      1.51926   6.279e-02
Objects:              5.300e+01      1.03922   5.200e+01
Flops:                7.449e+05      1.56738   6.355e+05  3.813e+06
Flops/sec:            1.321e+07      1.58747   1.025e+07  6.147e+07
MPI Messages:         1.125e+02      2.50000   6.800e+01  4.080e+02
MPI Message Lengths:  1.592e+04      1.64872   1.805e+02  7.365e+04
MPI Reductions:       7.800e+01      1.02632

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 6.2743e-02  99.9%  3.8128e+06 100.0%  4.080e+02 100.0%  1.805e+02      100.0%  6.100e+01  78.2% 

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

VecMDot               14 1.0 6.3467e-04 5.5 3.75e+04 1.3 0.0e+00 0.0e+00 1.4e+01  1  5  0  0 18   1  5  0  0 23   317
VecNorm               16 1.0 9.1290e-0316.5 5.73e+03 1.3 0.0e+00 0.0e+00 1.6e+01  8  1  0  0 21   8  1  0  0 26     3
VecScale              15 1.0 8.7738e-05 5.4 2.68e+03 1.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   164
VecCopy                4 1.0 4.2915e-06 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                20 1.0 1.2398e-05 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                2 1.0 1.3605e-02 2.4 7.16e+02 1.3 0.0e+00 0.0e+00 0.0e+00 14  0  0  0  0  14  0  0  0  0     0
VecMAXPY              15 1.0 2.0027e-05 1.5 4.26e+04 1.3 0.0e+00 0.0e+00 0.0e+00  0  6  0  0  0   0  6  0  0  0 11420
VecAssemblyBegin       3 1.0 9.3291e-03 1.2 0.00e+00 0.0 1.8e+01 1.2e+02 9.0e+00 14  0  4  3 12  14  0  4  3 15     0
VecAssemblyEnd         3 1.0 3.5286e-05 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       16 1.0 7.5340e-05 1.4 0.00e+00 0.0 2.9e+02 1.6e+02 0.0e+00  0  0 71 62  0   0  0 71 62  0     0
VecScatterEnd         16 1.0 6.5768e-03103.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  3  0  0  0  0   3  0  0  0  0     0
VecNormalize          15 1.0 8.6913e-0314.4 8.06e+03 1.3 0.0e+00 0.0e+00 1.5e+01  8  1  0  0 19   8  1  0  0 25     5
MatMult               15 1.0 6.7108e-0329.0 8.09e+04 1.4 2.7e+02 1.5e+02 0.0e+00  3 11 66 55  0   3 11 66 55  0    63
MatSolve              15 1.0 2.3365e-04 1.9 2.86e+05 1.5 0.0e+00 0.0e+00 0.0e+00  0 39  0  0  0   0 39  0  0  0  6360
MatLUFactorNum         1 1.0 2.9302e-04 1.5 2.90e+05 1.8 0.0e+00 0.0e+00 0.0e+00  0 37  0  0  0   0 37  0  0  0  4856
MatILUFactorSym        1 1.0 8.2421e-04 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+00  1  0  0  0  1   1  0  0  0  2     0
MatAssemblyBegin       2 1.0 2.0161e-03 3.0 0.00e+00 0.0 2.7e+01 8.2e+02 4.0e+00  3  0  7 30  5   3  0  7 30  7     0
MatAssemblyEnd         2 1.0 1.5681e-03 1.4 0.00e+00 0.0 3.6e+01 4.0e+01 8.0e+00  2  0  9  2 10   2  0  9  2 13     0
MatGetRowIJ            1 1.0 1.1590e-031215.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 1.2300e-0330.9 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  0  0  0  0  4   0  0  0  0  5     0
MatZeroEntries         3 1.0 3.4332e-05 2.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        14 1.0 6.7949e-04 3.9 7.51e+04 1.3 0.0e+00 0.0e+00 1.4e+01  1 11  0  0 18   1 11  0  0 23   593
KSPSetup               2 1.0 7.9870e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 1.8681e-02 1.0 7.45e+05 1.6 2.7e+02 1.5e+02 3.4e+01 30100 66 55 44  30100 66 55 56   204
PCSetUp                2 1.0 2.7039e-03 2.3 2.90e+05 1.8 0.0e+00 0.0e+00 4.0e+00  3 37  0  0  5   3 37  0  0  7   526
PCSetUpOnBlocks        1 1.0 2.4109e-03 2.9 2.90e+05 1.8 0.0e+00 0.0e+00 4.0e+00  2 37  0  0  5   2 37  0  0  7   590
PCApply               15 1.0 4.6277e-04 1.7 2.86e+05 1.5 0.0e+00 0.0e+00 0.0e+00  1 39  0  0  0   1 39  0  0  0  3211
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    33             33        84320     0
         Vec Scatter     2              2         1736     0
           Index Set     9              9         8516     0
   IS L to G Mapping     1              1         1360     0
              Matrix     4              4       160388     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 4.25816e-05
Average time for zero size MPI_Send(): 4.33524e-05
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
| Time:           Fri Aug 24 15:21:35 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.166934, Active time=0.047131                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0000      0.000049    0.0001      0.000066    0.10     0.14     |
|   build_sparsity()                 1         0.0004      0.000375    0.0005      0.000474    0.80     1.01     |
|   create_dof_constraints()         1         0.0000      0.000045    0.0000      0.000045    0.10     0.10     |
|   distribute_dofs()                1         0.0002      0.000198    0.0017      0.001721    0.42     3.65     |
|   dof_indices()                    291       0.0001      0.000000    0.0001      0.000000    0.27     0.27     |
|   prepare_send_list()              1         0.0000      0.000006    0.0000      0.000006    0.01     0.01     |
|   reinit()                         1         0.0003      0.000329    0.0003      0.000329    0.70     0.70     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        49        0.0002      0.000004    0.0002      0.000004    0.46     0.46     |
|   init_shape_functions()           13        0.0000      0.000002    0.0000      0.000002    0.05     0.05     |
|   inverse_map()                    36        0.0001      0.000002    0.0001      0.000002    0.16     0.16     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             49        0.0001      0.000001    0.0001      0.000001    0.14     0.14     |
|   compute_face_map()               12        0.0001      0.000006    0.0002      0.000013    0.16     0.33     |
|   init_face_shape_functions()      1         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   init_reference_to_physical_map() 13        0.0001      0.000008    0.0001      0.000008    0.22     0.22     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0002      0.000188    0.0002      0.000222    0.40     0.47     |
|   renumber_nodes_and_elem()        2         0.0000      0.000022    0.0000      0.000022    0.10     0.10     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0010      0.000476    0.0010      0.000476    2.02     2.02     |
|   find_global_indices()            2         0.0002      0.000078    0.0098      0.004876    0.33     20.69    |
|   parallel_sort()                  2         0.0017      0.000828    0.0032      0.001592    3.51     6.76     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0002      0.000203    0.0002      0.000203    0.43     0.43     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0008      0.000804    0.0030      0.003040    1.71     6.45     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      8         0.0058      0.000729    0.0058      0.000729    12.38    12.38    |
|   broadcast()                      1         0.0000      0.000012    0.0000      0.000012    0.03     0.03     |
|   gather()                         1         0.0000      0.000044    0.0000      0.000044    0.09     0.09     |
|   max(scalar)                      2         0.0002      0.000109    0.0002      0.000109    0.46     0.46     |
|   max(vector)                      2         0.0006      0.000289    0.0006      0.000289    1.23     1.23     |
|   min(vector)                      2         0.0004      0.000183    0.0004      0.000183    0.78     0.78     |
|   probe()                          50        0.0017      0.000034    0.0017      0.000034    3.56     3.56     |
|   receive()                        50        0.0001      0.000002    0.0018      0.000035    0.19     3.77     |
|   send()                           50        0.0000      0.000001    0.0000      0.000001    0.08     0.08     |
|   send_receive()                   54        0.0001      0.000002    0.0020      0.000036    0.19     4.16     |
|   sum()                            7         0.0016      0.000236    0.0016      0.000236    3.50     3.50     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           50        0.0000      0.000001    0.0000      0.000001    0.10     0.10     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0001      0.000087    0.0010      0.001048    0.18     2.22     |
|   set_parent_processor_ids()       1         0.0000      0.000017    0.0000      0.000017    0.04     0.04     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.0299      0.029888    0.0299      0.029888    63.41    63.41    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0008      0.000787    0.0014      0.001409    1.67     2.99     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            762       0.0471                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 6 ./introduction_ex3-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
