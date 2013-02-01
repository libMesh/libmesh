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
<br><br><br> <h1> The source file exact_solution.C with comments: </h1> 
<div class = "comment">
  

<br><br>This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
  

<br><br>This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
  

<br><br>You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


<br><br>

<br><br>

<br><br>C++ Includes
</div>

<div class ="fragment">
<pre>
        #include &lt;math.h&gt;
        
</pre>
</div>
<div class = "comment">
Mesh library includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh_common.h"
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        
        
        
        
        /**
         * This is the exact solution that
         * we are trying to obtain.  We will solve
         *
         * - (u_xx + u_yy) = f
         *
         * and take a finite difference approximation using this
         * function to get f.  This is the well-known "method of
         * manufactured solutions".
         */
        Real exact_solution (const Real x,
        		     const Real y,
        		     const Real z = 0.)
        {
          static const Real pi = acos(-1.);
        
          return cos(.5*pi*x)*sin(.5*pi*y)*cos(.5*pi*z);
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file introduction_ex3.C with comments: </h1> 
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
        #include "libmesh/libmesh.h"
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/vtk_io.h"
        #include "libmesh/linear_implicit_system.h"
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
matrix and vector components.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/dense_vector.h"
        #include "libmesh/elem.h"
        
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
          libmesh_assert_equal_to (system_name, "Poisson");
        
          
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
          std::vector&lt;dof_id_type&gt; dof_indices;
        
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
<br><br><br> <h1> The source file exact_solution.C without comments: </h1> 
<pre> 
    
    
    
  
  
  
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_common.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  
  
  
  <I><FONT COLOR="#B22222">/**
   * This is the exact solution that
   * we are trying to obtain.  We will solve
   *
   * - (u_xx + u_yy) = f
   *
   * and take a finite difference approximation using this
   * function to get f.  This is the well-known &quot;method of
   * manufactured solutions&quot;.
   */</FONT></I>
  Real exact_solution (<B><FONT COLOR="#228B22">const</FONT></B> Real x,
  		     <B><FONT COLOR="#228B22">const</FONT></B> Real y,
  		     <B><FONT COLOR="#228B22">const</FONT></B> Real z = 0.)
  {
    <B><FONT COLOR="#228B22">static</FONT></B> <B><FONT COLOR="#228B22">const</FONT></B> Real pi = acos(-1.);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> cos(.5*pi*x)*sin(.5*pi*y)*cos(.5*pi*z);
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file introduction_ex3.C without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/vtk_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  
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
    
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>);
  
    
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
  
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
  
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
***************************************************************
* Running Example introduction_ex3:
*  mpirun -np 12 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/introduction/introduction_ex3/.libs/lt-example-devel -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=961
    n_local_nodes()=97
  n_elem()=225
    n_local_elem()=19
    n_active_elem()=225
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
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
    n_local_dofs()=97
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 13.2112
      Average Off-Processor Bandwidth <= 2.7846
      Maximum  On-Processor Bandwidth <= 26
      Maximum Off-Processor Bandwidth <= 18
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

*** Warning, This code is untested, experimental, or likely to see future API changes: ../../../include/libmesh/vtk_io.h, line 178, compiled Jan 31 2013 at 21:56:01 ***
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/introduction/introduction_ex3/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 21:56:16 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           4.037e-01      1.00392   4.034e-01
Objects:              6.200e+01      1.00000   6.200e+01
Flops:                1.505e+06      1.72367   1.033e+06  1.240e+07
Flops/sec:            3.727e+06      1.72344   2.561e+06  3.074e+07
MPI Messages:         5.890e+02      3.48521   3.656e+02  4.387e+03
MPI Message Lengths:  5.099e+04      1.96516   9.952e+01  4.366e+05
MPI Reductions:       1.970e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 4.0334e-01 100.0%  1.2399e+07 100.0%  4.387e+03 100.0%  9.952e+01      100.0%  1.960e+02  99.5% 

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

VecMDot               74 1.0 1.0543e-03 1.7 2.00e+05 1.4 0.0e+00 0.0e+00 7.4e+01  0 16  0  0 38   0 16  0  0 38  1875
VecNorm               78 1.0 7.6747e-04 1.5 1.51e+04 1.4 0.0e+00 0.0e+00 7.8e+01  0  1  0  0 40   0  1  0  0 40   195
VecScale              77 1.0 7.6363e-03 1.0 7.47e+03 1.4 0.0e+00 0.0e+00 0.0e+00  2  1  0  0  0   2  1  0  0  0    10
VecCopy                4 1.0 8.8215e-06 4.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                86 1.0 4.7684e-05 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                6 1.0 5.5015e-02 1.0 1.16e+03 1.4 0.0e+00 0.0e+00 0.0e+00 14  0  0  0  0  14  0  0  0  0     0
VecMAXPY              77 1.0 1.4448e-04 1.3 2.15e+05 1.4 0.0e+00 0.0e+00 0.0e+00  0 17  0  0  0   0 17  0  0  0 14753
VecAssemblyBegin       3 1.0 7.9703e-04 2.1 0.00e+00 0.0 4.8e+01 7.8e+01 9.0e+00  0  0  1  1  5   0  0  1  1  5     0
VecAssemblyEnd         3 1.0 2.5988e-05 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       78 1.0 3.5858e-04 1.7 0.00e+00 0.0 4.1e+03 9.6e+01 0.0e+00  0  0 92 89  0   0  0 92 89  0     0
VecScatterEnd         78 1.0 5.1165e-04 3.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize          77 1.0 8.4391e-03 1.0 2.24e+04 1.4 0.0e+00 0.0e+00 7.7e+01  2  2  0  0 39   2  2  0  0 39    26
MatMult               77 1.0 1.1570e-03 1.7 2.18e+05 1.5 4.0e+03 9.5e+01 0.0e+00  0 18 91 87  0   0 18 91 87  0  1885
MatSolve              78 1.0 5.5480e-04 1.7 7.20e+05 2.0 0.0e+00 0.0e+00 0.0e+00  0 42  0  0  0   0 42  0  0  0  9369
MatLUFactorNum         1 1.0 8.0452e-02 1.0 1.28e+05 3.5 0.0e+00 0.0e+00 0.0e+00 20  5  0  0  0  20  5  0  0  0     8
MatILUFactorSym        1 1.0 5.2309e-04 3.0 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  0  0  0  0  2   0  0  0  0  2     0
MatAssemblyBegin       2 1.0 3.2520e-04 1.6 0.00e+00 0.0 7.2e+01 5.3e+02 4.0e+00  0  0  2  9  2   0  0  2  9  2     0
MatAssemblyEnd         2 1.0 2.5442e-03 1.0 0.00e+00 0.0 1.0e+02 2.6e+01 8.0e+00  1  0  2  1  4   1  0  2  1  4     0
MatGetRowIJ            1 1.0 2.1458e-06 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 7.2956e-05 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  1   0  0  0  0  1     0
MatZeroEntries         3 1.0 2.0981e-05 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        74 1.0 1.2789e-03 1.5 4.01e+05 1.4 0.0e+00 0.0e+00 7.4e+01  0 32  0  0 38   0 32  0  0 38  3101
KSPSetUp               2 1.0 9.7990e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 1.4957e-01 1.0 1.50e+06 1.7 4.0e+03 9.5e+01 1.6e+02 37100 91 87 81  37100 91 87 81    83
PCSetUp                2 1.0 8.2732e-02 1.0 1.28e+05 3.5 0.0e+00 0.0e+00 7.0e+00 20  5  0  0  4  20  5  0  0  4     8
PCSetUpOnBlocks        1 1.0 8.0977e-02 1.0 1.28e+05 3.5 0.0e+00 0.0e+00 5.0e+00 20  5  0  0  3  20  5  0  0  3     8
PCApply               78 1.0 1.3814e-03 1.3 7.20e+05 2.0 0.0e+00 0.0e+00 0.0e+00  0 42  0  0  0   0 42  0  0  0  3763
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    43             43        94328     0
      Vector Scatter     2              2         2072     0
           Index Set     7              7         5964     0
   IS L to G Mapping     1              1          564     0
              Matrix     4              4        87984     0
       Krylov Solver     2              2        19360     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 5.09739e-05
Average time for zero size MPI_Send(): 5.17567e-05
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
| Time:           Thu Jan 31 21:56:16 2013                                                                             |
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
| libMesh Performance: Alive time=0.534381, Active time=0.384137                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0025      0.002486    0.0050      0.004987    0.65     1.30     |
|   build_sparsity()                 1         0.0021      0.002071    0.0059      0.005945    0.54     1.55     |
|   create_dof_constraints()         1         0.0005      0.000531    0.0005      0.000531    0.14     0.14     |
|   distribute_dofs()                1         0.0117      0.011678    0.0652      0.065180    3.04     16.97    |
|   dof_indices()                    78        0.0092      0.000119    0.0092      0.000119    2.41     2.41     |
|   prepare_send_list()              1         0.0001      0.000062    0.0001      0.000062    0.02     0.02     |
|   reinit()                         1         0.0201      0.020128    0.0201      0.020128    5.24     5.24     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0006      0.000607    0.0034      0.003394    0.16     0.88     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        29        0.0011      0.000037    0.0011      0.000037    0.28     0.28     |
|   init_shape_functions()           11        0.0003      0.000026    0.0003      0.000026    0.08     0.08     |
|   inverse_map()                    30        0.0005      0.000018    0.0005      0.000018    0.14     0.14     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             29        0.0005      0.000018    0.0005      0.000018    0.13     0.13     |
|   compute_face_map()               10        0.0005      0.000045    0.0010      0.000101    0.12     0.26     |
|   init_face_shape_functions()      1         0.0000      0.000019    0.0000      0.000019    0.00     0.00     |
|   init_reference_to_physical_map() 11        0.0007      0.000061    0.0007      0.000061    0.18     0.18     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0055      0.005508    0.0058      0.005793    1.43     1.51     |
|   renumber_nodes_and_elem()        2         0.0006      0.000280    0.0006      0.000280    0.15     0.15     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0040      0.001976    0.0040      0.001976    1.03     1.03     |
|   find_global_indices()            2         0.0020      0.000977    0.0104      0.005191    0.51     2.70     |
|   parallel_sort()                  2         0.0026      0.001321    0.0032      0.001593    0.69     0.83     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0065      0.006522    0.0101      0.010078    1.70     2.62     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0028      0.002768    0.0028      0.002768    0.72     0.72     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0247      0.024688    0.0291      0.029103    6.43     7.58     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      9         0.0221      0.002454    0.0221      0.002460    5.75     5.76     |
|   max(bool)                        1         0.0000      0.000007    0.0000      0.000007    0.00     0.00     |
|   max(scalar)                      105       0.0011      0.000010    0.0011      0.000010    0.29     0.29     |
|   max(vector)                      24        0.0004      0.000018    0.0012      0.000049    0.11     0.31     |
|   min(bool)                        121       0.0012      0.000010    0.0012      0.000010    0.31     0.31     |
|   min(scalar)                      99        0.0154      0.000156    0.0154      0.000156    4.01     4.01     |
|   min(vector)                      24        0.0005      0.000020    0.0014      0.000058    0.13     0.36     |
|   probe()                          132       0.0048      0.000037    0.0048      0.000037    1.26     1.26     |
|   receive()                        132       0.0008      0.000006    0.0057      0.000043    0.22     1.49     |
|   send()                           132       0.0004      0.000003    0.0004      0.000003    0.11     0.11     |
|   send_receive()                   136       0.0012      0.000009    0.0076      0.000056    0.30     1.99     |
|   sum()                            20        0.0006      0.000029    0.0009      0.000045    0.15     0.23     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           132       0.0003      0.000002    0.0003      0.000002    0.07     0.07     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0012      0.001217    0.0019      0.001873    0.32     0.49     |
|   set_parent_processor_ids()       1         0.0005      0.000453    0.0005      0.000453    0.12     0.12     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.2320      0.232008    0.2320      0.232008    60.40    60.40    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0027      0.002659    0.0085      0.008503    0.69     2.21     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1289      0.3841                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example introduction_ex3:
*  mpirun -np 12 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
