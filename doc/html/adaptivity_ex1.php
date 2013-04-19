<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("adaptivity_ex1",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file adaptivity_ex1.C with comments: </h1> 
<div class = "comment">
<h1>Adaptivity Example 1 - Solving 1D PDE Using Adaptive Mesh Refinement</h1>

<br><br>This example demonstrates how to solve a simple 1D problem
using adaptive mesh refinement. The PDE that is solved is:
-epsilon*u''(x) + u(x) = 1, on the domain [0,1] with boundary conditions
u(0) = u(1) = 0 and where epsilon << 1.

<br><br>The approach used to solve 1D problems in libMesh is virtually identical to
solving 2D or 3D problems, so in this sense this example represents a good
starting point for new users. Note that many concepts are used in this
example which are explained more fully in subsequent examples.


<br><br>Libmesh includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/edge_edge3.h"
        #include "libmesh/gnuplot_io.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/linear_implicit_system.h"
        #include "libmesh/fe.h"
        #include "libmesh/getpot.h"
        #include "libmesh/quadrature_gauss.h"
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/dense_vector.h"
        #include "libmesh/error_vector.h"
        #include "libmesh/kelly_error_estimator.h"
        #include "libmesh/mesh_refinement.h"
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        void assemble_1D(EquationSystems& es, const std::string& system_name);
        
        int main(int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Initialize the library.  This is necessary because the library
may depend on a number of other libraries (i.e. MPI and PETSc)
that require initialization before use.  When the LibMeshInit
object goes out of scope, other libraries and resources are
finalized.
</div>

<div class ="fragment">
<pre>
          LibMeshInit init (argc, argv);
        
</pre>
</div>
<div class = "comment">
Skip adaptive examples on a non-adaptive libMesh build
</div>

<div class ="fragment">
<pre>
        #ifndef LIBMESH_ENABLE_AMR
          libmesh_example_assert(false, "--enable-amr");
        #else
        
</pre>
</div>
<div class = "comment">
Create a mesh, with dimension to be overridden later, on the
default MPI communicator.
</div>

<div class ="fragment">
<pre>
          Mesh mesh(init.comm());
        
          GetPot command_line (argc, argv);
        
          int n = 4;
          if ( command_line.search(1, "-n") )
            n = command_line.next(n);
        
</pre>
</div>
<div class = "comment">
Build a 1D mesh with 4 elements from x=0 to x=1, using
EDGE3 (i.e. quadratic) 1D elements. They are called EDGE3 elements
because a quadratic element contains 3 nodes.
</div>

<div class ="fragment">
<pre>
          MeshTools::Generation::build_line(mesh,n,0.,1.,EDGE3);
        
</pre>
</div>
<div class = "comment">
Define the equation systems object and the system we are going
to solve. See Introduction Example 2 for more details.
</div>

<div class ="fragment">
<pre>
          EquationSystems equation_systems(mesh);
          LinearImplicitSystem& system = equation_systems.add_system
            &lt;LinearImplicitSystem&gt;("1D");
        
</pre>
</div>
<div class = "comment">
Add a variable "u" to the system, using second-order approximation
</div>

<div class ="fragment">
<pre>
          system.add_variable("u",SECOND);
        
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly function. This
will be called when needed by the library.
</div>

<div class ="fragment">
<pre>
          system.attach_assemble_function(assemble_1D);
        
</pre>
</div>
<div class = "comment">
Define the mesh refinement object that takes care of adaptively
refining the mesh.
</div>

<div class ="fragment">
<pre>
          MeshRefinement mesh_refinement(mesh);
        
</pre>
</div>
<div class = "comment">
These parameters determine the proportion of elements that will
be refined and coarsened. Any element within 30% of the maximum
error on any element will be refined, and any element within 30%
of the minimum error on any element might be coarsened
</div>

<div class ="fragment">
<pre>
          mesh_refinement.refine_fraction()  = 0.7;
          mesh_refinement.coarsen_fraction() = 0.3;
</pre>
</div>
<div class = "comment">
We won't refine any element more than 5 times in total
</div>

<div class ="fragment">
<pre>
          mesh_refinement.max_h_level()      = 5;
        
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
Refinement parameters
</div>

<div class ="fragment">
<pre>
          const unsigned int max_r_steps = 5; // Refine the mesh 5 times
        
</pre>
</div>
<div class = "comment">
Define the refinement loop
</div>

<div class ="fragment">
<pre>
          for(unsigned int r_step=0; r_step&lt;=max_r_steps; r_step++)
            {
</pre>
</div>
<div class = "comment">
Solve the equation system
</div>

<div class ="fragment">
<pre>
              equation_systems.get_system("1D").solve();
        
</pre>
</div>
<div class = "comment">
We need to ensure that the mesh is not refined on the last iteration
of this loop, since we do not want to refine the mesh unless we are
going to solve the equation system for that refined mesh.
</div>

<div class ="fragment">
<pre>
              if(r_step != max_r_steps)
                {
</pre>
</div>
<div class = "comment">
Error estimation objects, see Adaptivity Example 2 for details
</div>

<div class ="fragment">
<pre>
                  ErrorVector error;
                  KellyErrorEstimator error_estimator;
        
</pre>
</div>
<div class = "comment">
Compute the error for each active element
</div>

<div class ="fragment">
<pre>
                  error_estimator.estimate_error(system, error);
        
</pre>
</div>
<div class = "comment">
Output error estimate magnitude
</div>

<div class ="fragment">
<pre>
                  libMesh::out &lt;&lt; "Error estimate\nl2 norm = " &lt;&lt; error.l2_norm() &lt;&lt;
                    "\nmaximum = " &lt;&lt; error.maximum() &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Flag elements to be refined and coarsened
</div>

<div class ="fragment">
<pre>
                  mesh_refinement.flag_elements_by_error_fraction (error);
        
</pre>
</div>
<div class = "comment">
Perform refinement and coarsening
</div>

<div class ="fragment">
<pre>
                  mesh_refinement.refine_and_coarsen_elements();
        
</pre>
</div>
<div class = "comment">
Reinitialize the equation_systems object for the newly refined
mesh. One of the steps in this is project the solution onto the
new mesh
</div>

<div class ="fragment">
<pre>
                  equation_systems.reinit();
                }
            }
        
</pre>
</div>
<div class = "comment">
Construct gnuplot plotting object, pass in mesh, title of plot
and boolean to indicate use of grid in plot. The grid is used to
show the edges of each element in the mesh.
</div>

<div class ="fragment">
<pre>
          GnuPlotIO plot(mesh,"Adaptivity Example 1", GnuPlotIO::GRID_ON);
        
</pre>
</div>
<div class = "comment">
Write out script to be called from within gnuplot:
Load gnuplot, then type "call 'gnuplot_script'" from gnuplot prompt
</div>

<div class ="fragment">
<pre>
          plot.write_equation_systems("gnuplot_script",equation_systems);
        #endif // #ifndef LIBMESH_ENABLE_AMR
        
</pre>
</div>
<div class = "comment">
All done.  libMesh objects are destroyed here.  Because the
LibMeshInit object was created first, its destruction occurs
last, and it's destructor finalizes any external libraries and
checks for leaked memory.
</div>

<div class ="fragment">
<pre>
          return 0;
        }
        
        
        
        
</pre>
</div>
<div class = "comment">
Define the matrix assembly function for the 1D PDE we are solving
</div>

<div class ="fragment">
<pre>
        void assemble_1D(EquationSystems& es, const std::string& system_name)
        {
        
        #ifdef LIBMESH_ENABLE_AMR
        
</pre>
</div>
<div class = "comment">
It is a good idea to check we are solving the correct system
</div>

<div class ="fragment">
<pre>
          libmesh_assert_equal_to (system_name, "1D");
        
</pre>
</div>
<div class = "comment">
Get a reference to the mesh object
</div>

<div class ="fragment">
<pre>
          const MeshBase& mesh = es.get_mesh();
        
</pre>
</div>
<div class = "comment">
The dimension we are using, i.e. dim==1
</div>

<div class ="fragment">
<pre>
          const unsigned int dim = mesh.mesh_dimension();
        
</pre>
</div>
<div class = "comment">
Get a reference to the system we are solving
</div>

<div class ="fragment">
<pre>
          LinearImplicitSystem& system = es.get_system&lt;LinearImplicitSystem&gt;("1D");
        
</pre>
</div>
<div class = "comment">
Get a reference to the DofMap object for this system. The DofMap object
handles the index translation from node and element numbers to degree of
freedom numbers. DofMap's are discussed in more detail in future examples.
</div>

<div class ="fragment">
<pre>
          const DofMap& dof_map = system.get_dof_map();
        
</pre>
</div>
<div class = "comment">
Get a constant reference to the Finite Element type for the first
(and only) variable in the system.
</div>

<div class ="fragment">
<pre>
          FEType fe_type = dof_map.variable_type(0);
        
</pre>
</div>
<div class = "comment">
Build a finite element object of the specified type. The build
function dynamically allocates memory so we use an AutoPtr in this case.
An AutoPtr is a pointer that cleans up after itself. See examples 3 and 4
for more details on AutoPtr.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;FEBase&gt; fe(FEBase::build(dim, fe_type));
        
</pre>
</div>
<div class = "comment">
Tell the finite element object to use fifth order Gaussian quadrature
</div>

<div class ="fragment">
<pre>
          QGauss qrule(dim,FIFTH);
          fe-&gt;attach_quadrature_rule(&qrule);
        
</pre>
</div>
<div class = "comment">
Here we define some references to cell-specific data that will be used to
assemble the linear system.


<br><br>The element Jacobian * quadrature weight at each integration point.
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
The element shape function gradients evaluated at the quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi = fe-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
Declare a dense matrix and dense vector to hold the element matrix
and right-hand-side contribution
</div>

<div class ="fragment">
<pre>
          DenseMatrix&lt;Number&gt; Ke;
          DenseVector&lt;Number&gt; Fe;
        
</pre>
</div>
<div class = "comment">
This vector will hold the degree of freedom indices for the element.
These define where in the global system the element degrees of freedom
get mapped.
</div>

<div class ="fragment">
<pre>
          std::vector&lt;dof_id_type&gt; dof_indices;
        
</pre>
</div>
<div class = "comment">
We now loop over all the active elements in the mesh in order to calculate
the matrix and right-hand-side contribution from each element. Use a
const_element_iterator to loop over the elements. We make
el_end const as it is used only for the stopping condition of the loop.
</div>

<div class ="fragment">
<pre>
          MeshBase::const_element_iterator el     = mesh.active_local_elements_begin();
          const MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
        
</pre>
</div>
<div class = "comment">
Note that ++el is preferred to el++ when using loops with iterators
</div>

<div class ="fragment">
<pre>
          for( ; el != el_end; ++el)
          {
</pre>
</div>
<div class = "comment">
It is convenient to store a pointer to the current element
</div>

<div class ="fragment">
<pre>
            const Elem* elem = *el;
        
</pre>
</div>
<div class = "comment">
Get the degree of freedom indices for the current element.
These define where in the global matrix and right-hand-side this
element will contribute to.
</div>

<div class ="fragment">
<pre>
            dof_map.dof_indices(elem, dof_indices);
        
</pre>
</div>
<div class = "comment">
Compute the element-specific data for the current element. This
involves computing the location of the quadrature points (q_point)
and the shape functions (phi, dphi) for the current element.
</div>

<div class ="fragment">
<pre>
            fe-&gt;reinit(elem);
        
</pre>
</div>
<div class = "comment">
Store the number of local degrees of freedom contained in this element
</div>

<div class ="fragment">
<pre>
            const int n_dofs = dof_indices.size();
        
</pre>
</div>
<div class = "comment">
We resize and zero out Ke and Fe (resize() also clears the matrix and
vector). In this example, all elements in the mesh are EDGE3's, so
Ke will always be 3x3, and Fe will always be 3x1. If the mesh contained
different element types, then the size of Ke and Fe would change.
</div>

<div class ="fragment">
<pre>
            Ke.resize(n_dofs, n_dofs);
            Fe.resize(n_dofs);
        
        
</pre>
</div>
<div class = "comment">
Now loop over quadrature points to handle numerical integration
</div>

<div class ="fragment">
<pre>
            for(unsigned int qp=0; qp&lt;qrule.n_points(); qp++)
            {
</pre>
</div>
<div class = "comment">
Now build the element matrix and right-hand-side using loops to
integrate the test functions (i) against the trial functions (j).
</div>

<div class ="fragment">
<pre>
              for(unsigned int i=0; i&lt;phi.size(); i++)
              {
                Fe(i) += JxW[qp]*phi[i][qp];
        
                for(unsigned int j=0; j&lt;phi.size(); j++)
                {
                  Ke(i,j) += JxW[qp]*(1.e-3*dphi[i][qp]*dphi[j][qp] +
                                             phi[i][qp]*phi[j][qp]);
                }
              }
            }
        
        
</pre>
</div>
<div class = "comment">
At this point we have completed the matrix and RHS summation. The
final step is to apply boundary conditions, which in this case are
simple Dirichlet conditions with u(0) = u(1) = 0.


<br><br>Define the penalty parameter used to enforce the BC's
</div>

<div class ="fragment">
<pre>
            double penalty = 1.e10;
        
</pre>
</div>
<div class = "comment">
Loop over the sides of this element. For a 1D element, the "sides"
are defined as the nodes on each edge of the element, i.e. 1D elements
have 2 sides.
</div>

<div class ="fragment">
<pre>
            for(unsigned int s=0; s&lt;elem-&gt;n_sides(); s++)
            {
</pre>
</div>
<div class = "comment">
If this element has a NULL neighbor, then it is on the edge of the
mesh and we need to enforce a boundary condition using the penalty
method.
</div>

<div class ="fragment">
<pre>
              if(elem-&gt;neighbor(s) == NULL)
              {
                Ke(s,s) += penalty;
                Fe(s)   += 0*penalty;
              }
            }
        
</pre>
</div>
<div class = "comment">
This is a function call that is necessary when using adaptive
mesh refinement. See Adaptivity Example 2 for more details.
</div>

<div class ="fragment">
<pre>
            dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        
</pre>
</div>
<div class = "comment">
Add Ke and Fe to the global matrix and right-hand-side.
</div>

<div class ="fragment">
<pre>
            system.matrix-&gt;add_matrix(Ke, dof_indices);
            system.rhs-&gt;add_vector(Fe, dof_indices);
          }
        #endif // #ifdef LIBMESH_ENABLE_AMR
        }
</pre>
</div>

<a name="nocomments"></a> 
<br><br><br> <h1> The source file adaptivity_ex1.C without comments: </h1> 
<pre> 
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/edge_edge3.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/gnuplot_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/kelly_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_refinement.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_1D(EquationSystems&amp; es, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main(<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #ifndef LIBMESH_ENABLE_AMR
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-amr&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    Mesh mesh(init.comm());
  
    GetPot command_line (argc, argv);
  
    <B><FONT COLOR="#228B22">int</FONT></B> n = 4;
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <B><FONT COLOR="#BC8F8F">&quot;-n&quot;</FONT></B>) )
      n = command_line.next(n);
  
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_line(mesh,n,0.,1.,EDGE3);
  
    EquationSystems equation_systems(mesh);
    LinearImplicitSystem&amp; system = equation_systems.add_system
      &lt;LinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;1D&quot;</FONT></B>);
  
    system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>,SECOND);
  
    system.attach_assemble_function(assemble_1D);
  
    MeshRefinement mesh_refinement(mesh);
  
    mesh_refinement.refine_fraction()  = 0.7;
    mesh_refinement.coarsen_fraction() = 0.3;
    mesh_refinement.max_h_level()      = 5;
  
    equation_systems.init();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_r_steps = 5; <I><FONT COLOR="#B22222">// Refine the mesh 5 times
</FONT></I>  
    <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> r_step=0; r_step&lt;=max_r_steps; r_step++)
      {
        equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;1D&quot;</FONT></B>).solve();
  
        <B><FONT COLOR="#A020F0">if</FONT></B>(r_step != max_r_steps)
          {
            ErrorVector error;
            KellyErrorEstimator error_estimator;
  
            error_estimator.estimate_error(system, error);
  
  	  <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::out &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error estimate\nl2 norm = &quot;</FONT></B> &lt;&lt; error.l2_norm() &lt;&lt;
              <B><FONT COLOR="#BC8F8F">&quot;\nmaximum = &quot;</FONT></B> &lt;&lt; error.maximum() &lt;&lt; std::endl;
  
            mesh_refinement.flag_elements_by_error_fraction (error);
  
            mesh_refinement.refine_and_coarsen_elements();
  
            equation_systems.reinit();
          }
      }
  
    GnuPlotIO plot(mesh,<B><FONT COLOR="#BC8F8F">&quot;Adaptivity Example 1&quot;</FONT></B>, GnuPlotIO::GRID_ON);
  
    plot.write_equation_systems(<B><FONT COLOR="#BC8F8F">&quot;gnuplot_script&quot;</FONT></B>,equation_systems);
  #endif <I><FONT COLOR="#B22222">// #ifndef LIBMESH_ENABLE_AMR
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_1D(EquationSystems&amp; es, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
  
  #ifdef LIBMESH_ENABLE_AMR
  
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;1D&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    LinearImplicitSystem&amp; system = es.get_system&lt;LinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;1D&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
  
    FEType fe_type = dof_map.variable_type(0);
  
    AutoPtr&lt;FEBase&gt; fe(FEBase::build(dim, fe_type));
  
    QGauss qrule(dim,FIFTH);
    fe-&gt;attach_quadrature_rule(&amp;qrule);
  
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    DenseMatrix&lt;Number&gt; Ke;
    DenseVector&lt;Number&gt; Fe;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator el_end = mesh.active_local_elements_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B>( ; el != el_end; ++el)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
      dof_map.dof_indices(elem, dof_indices);
  
      fe-&gt;reinit(elem);
  
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_dofs = dof_indices.size();
  
      Ke.resize(n_dofs, n_dofs);
      Fe.resize(n_dofs);
  
  
      <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
      {
        <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
        {
          Fe(i) += JxW[qp]*phi[i][qp];
  
          <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
          {
            Ke(i,j) += JxW[qp]*(1.e-3*dphi[i][qp]*dphi[j][qp] +
                                       phi[i][qp]*phi[j][qp]);
          }
        }
      }
  
  
  
      <B><FONT COLOR="#228B22">double</FONT></B> penalty = 1.e10;
  
      <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> s=0; s&lt;elem-&gt;n_sides(); s++)
      {
        <B><FONT COLOR="#A020F0">if</FONT></B>(elem-&gt;neighbor(s) == NULL)
        {
          Ke(s,s) += penalty;
          Fe(s)   += 0*penalty;
        }
      }
  
      dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
  
      system.matrix-&gt;add_matrix(Ke, dof_indices);
      system.rhs-&gt;add_vector(Fe, dof_indices);
    }
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_ENABLE_AMR
</FONT></I>  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex1'
***************************************************************
* Running Example adaptivity_ex1:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
Error estimate
l2 norm = 6.92141
maximum = 3.48539
Error estimate
l2 norm = 4.26669
maximum = 2.13334
Error estimate
l2 norm = 1.4544
maximum = 0.727055
Error estimate
l2 norm = 0.308552
maximum = 0.152269
Error estimate
l2 norm = 0.0599993
maximum = 0.0230905

 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:45:06 2013                                                                          |
| OS:             Linux                                                                                             |
| HostName:       spark.ices.utexas.edu                                                                             |
| OS Release:     2.6.32-279.22.1.el6.x86_64                                                                        |
| OS Version:     #1 SMP Tue Feb 5 14:33:39 CST 2013                                                                |
| Machine:        x86_64                                                                                            |
| Username:       roystgnr                                                                                          |
| Configuration:  ../configure  '--enable-everything'                                                               |
|  'METHODS=devel'                                                                                                  |
|  '--prefix=/h2/roystgnr/libmesh-test'                                                                             |
|  'CXX=distcc /usr/bin/g++'                                                                                        |
|  'CC=distcc /usr/bin/gcc'                                                                                         |
|  'FC=distcc /usr/bin/gfortran'                                                                                    |
|  'F77=distcc /usr/bin/gfortran'                                                                                   |
|  'PETSC_DIR=/opt/apps/ossw/libraries/petsc/petsc-3.3-p2'                                                          |
|  'PETSC_ARCH=gcc-system-mkl-gf-10.3.12.361-mpich2-1.4.1p1-cxx-opt'                                                |
|  'SLEPC_DIR=/opt/apps/ossw/libraries/slepc/slepc-3.3-p2-petsc-3.3-p2-cxx-opt'                                     |
|  'TRILINOS_DIR=/opt/apps/ossw/libraries/trilinos/trilinos-10.12.2/sl6/gcc-system/mpich2-1.4.1p1/mkl-gf-10.3.12.361'|
|  'VTK_DIR=/opt/apps/ossw/libraries/vtk/vtk-5.10.0/sl6/gcc-system'                                                 |
|  'HDF5_DIR=/opt/apps/ossw/libraries/hdf5/hdf5-1.8.9/sl6/gcc-system'                                               |
 -------------------------------------------------------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.172276, Active time=0.050797                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     6         0.0001      0.000025    0.0002      0.000033    0.29     0.39     |
|   build_sparsity()                 6         0.0004      0.000061    0.0015      0.000251    0.72     2.96     |
|   create_dof_constraints()         6         0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|   distribute_dofs()                6         0.0007      0.000122    0.0036      0.000598    1.44     7.06     |
|   dof_indices()                    128       0.0004      0.000003    0.0004      0.000003    0.83     0.83     |
|   old_dof_indices()                48        0.0001      0.000003    0.0001      0.000003    0.26     0.26     |
|   prepare_send_list()              6         0.0000      0.000002    0.0000      0.000002    0.02     0.02     |
|   reinit()                         6         0.0005      0.000088    0.0005      0.000088    1.04     1.04     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0001      0.000079    0.0004      0.000353    0.16     0.69     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        55        0.0001      0.000001    0.0001      0.000001    0.13     0.13     |
|   init_shape_functions()           36        0.0001      0.000003    0.0001      0.000003    0.20     0.20     |
|   inverse_map()                    68        0.0001      0.000002    0.0001      0.000002    0.21     0.21     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             55        0.0001      0.000002    0.0001      0.000002    0.19     0.19     |
|   compute_face_map()               15        0.0000      0.000002    0.0000      0.000002    0.06     0.06     |
|   init_face_shape_functions()      5         0.0000      0.000003    0.0000      0.000003    0.03     0.03     |
|   init_reference_to_physical_map() 36        0.0001      0.000003    0.0001      0.000003    0.18     0.18     |
|                                                                                                                |
| GnuPlotIO                                                                                                      |
|   write_nodal_data()               1         0.0006      0.000556    0.0006      0.000556    1.09     1.09     |
|                                                                                                                |
| JumpErrorEstimator                                                                                             |
|   estimate_error()                 5         0.0004      0.000087    0.0013      0.000250    0.86     2.46     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           76        0.0001      0.000001    0.0001      0.000001    0.22     0.22     |
|   init()                           10        0.0002      0.000016    0.0002      0.000016    0.31     0.31     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       5         0.0000      0.000007    0.0001      0.000011    0.07     0.11     |
|   find_neighbors()                 6         0.0005      0.000084    0.0011      0.000177    0.99     2.09     |
|   renumber_nodes_and_elem()        17        0.0001      0.000005    0.0001      0.000005    0.15     0.15     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        7         0.0004      0.000051    0.0004      0.000051    0.70     0.70     |
|   find_global_indices()            7         0.0002      0.000033    0.0024      0.000345    0.46     4.76     |
|   parallel_sort()                  7         0.0006      0.000081    0.0010      0.000150    1.11     2.06     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0001      0.000076    0.0011      0.001068    0.15     2.10     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              10        0.0000      0.000004    0.0001      0.000011    0.09     0.21     |
|   _refine_elements()               10        0.0003      0.000029    0.0009      0.000087    0.57     1.72     |
|   add_point()                      76        0.0001      0.000002    0.0003      0.000003    0.26     0.50     |
|   make_coarsening_compatible()     21        0.0007      0.000031    0.0007      0.000031    1.30     1.30     |
|   make_flags_parallel_consistent() 15        0.0007      0.000048    0.0051      0.000342    1.42     10.10    |
|   make_refinement_compatible()     21        0.0001      0.000004    0.0002      0.000011    0.18     0.44     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0001      0.000082    0.0001      0.000082    0.16     0.16     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      6         0.0013      0.000224    0.0040      0.000669    2.65     7.90     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      32        0.0006      0.000018    0.0007      0.000021    1.16     1.34     |
|   max(bool)                        52        0.0007      0.000013    0.0007      0.000013    1.32     1.32     |
|   max(scalar)                      1081      0.0044      0.000004    0.0044      0.000004    8.69     8.69     |
|   max(vector)                      264       0.0016      0.000006    0.0049      0.000019    3.18     9.62     |
|   min(bool)                        1378      0.0060      0.000004    0.0060      0.000004    11.82    11.82    |
|   min(scalar)                      1069      0.0069      0.000006    0.0069      0.000006    13.55    13.55    |
|   min(vector)                      264       0.0018      0.000007    0.0050      0.000019    3.54     9.91     |
|   probe()                          366       0.0019      0.000005    0.0019      0.000005    3.75     3.75     |
|   receive()                        366       0.0008      0.000002    0.0027      0.000007    1.49     5.31     |
|   send()                           366       0.0004      0.000001    0.0004      0.000001    0.88     0.88     |
|   send_receive()                   380       0.0012      0.000003    0.0047      0.000012    2.27     9.17     |
|   sum()                            34        0.0005      0.000013    0.0010      0.000030    0.90     2.04     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           366       0.0003      0.000001    0.0003      0.000001    0.51     0.51     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         6         0.0003      0.000052    0.0016      0.000260    0.61     3.07     |
|   set_parent_processor_ids()       6         0.0001      0.000009    0.0001      0.000009    0.11     0.11     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          6         0.0122      0.002030    0.0122      0.002030    23.98    23.98    |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       5         0.0002      0.000042    0.0004      0.000086    0.41     0.85     |
|                                                                                                                |
| StatisticsVector                                                                                               |
|   maximum()                        5         0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       6         0.0003      0.000052    0.0006      0.000097    0.61     1.15     |
|   project_vector()                 5         0.0014      0.000274    0.0026      0.000517    2.70     5.09     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            6842      0.0508                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example adaptivity_ex1:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/adaptivity/adaptivity_ex1'
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
