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
Create a new mesh
</div>

<div class ="fragment">
<pre>
          Mesh mesh;
        
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
  
    Mesh mesh;
  
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
***************************************************************
* Running Example adaptivity_ex1:
*  mpirun -np 12 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/adaptivity/adaptivity_ex1/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 21:57:46 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           1.687e-01      1.00000   1.687e-01
Objects:              3.200e+02      1.02564   3.158e+02
Flops:                2.584e+04      2.23070   1.874e+04  2.248e+05
Flops/sec:            1.531e+05      2.23070   1.111e+05  1.333e+06
MPI Messages:         3.225e+02      3.02817   2.087e+02  2.504e+03
MPI Message Lengths:  3.616e+03      2.86529   1.127e+01  2.822e+04
MPI Reductions:       5.620e+02      1.01444

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.6863e-01 100.0%  2.2484e+05 100.0%  2.504e+03 100.0%  1.127e+01      100.0%  5.568e+02  99.1% 

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

VecMDot               81 1.0 1.2305e-03 1.2 7.53e+03 2.1 0.0e+00 0.0e+00 8.1e+01  1 30  0  0 14   1 30  0  0 15    55
VecNorm               93 1.0 9.9349e-03 8.5 1.01e+03 2.3 0.0e+00 0.0e+00 9.3e+01  4  4  0  0 17   4  4  0  0 17     1
VecScale              87 1.0 5.2691e-05 1.5 4.75e+02 2.3 0.0e+00 0.0e+00 0.0e+00  0  2  0  0  0   0  2  0  0  0    75
VecCopy               27 1.0 2.6703e-05 2.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               149 1.0 6.7234e-05 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               12 1.0 1.6253e-02 2.2 1.16e+02 2.6 0.0e+00 0.0e+00 0.0e+00  6  0  0  0  0   6  0  0  0  0     0
VecMAXPY              87 1.0 3.1710e-05 1.8 9.11e+03 2.1 0.0e+00 0.0e+00 0.0e+00  0 37  0  0  0   0 37  0  0  0  2599
VecAssemblyBegin      53 1.0 3.1378e-03 1.4 0.00e+00 0.0 1.6e+02 6.0e+00 1.3e+02  2  0  6  3 23   2  0  6  3 23     0
VecAssemblyEnd        53 1.0 8.7261e-05 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      108 1.0 2.7657e-04 1.8 0.00e+00 0.0 1.7e+03 1.3e+01 0.0e+00  0  0 69 78  0   0  0 69 78  0     0
VecScatterEnd        108 1.0 4.7207e-04 5.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize          87 1.0 9.9187e-03 8.5 1.42e+03 2.3 0.0e+00 0.0e+00 8.7e+01  4  5  0  0 15   4  5  0  0 16     1
MatMult               87 1.0 6.1297e-04 2.4 3.28e+03 2.4 1.4e+03 1.2e+01 0.0e+00  0 12 57 58  0   0 12 57 58  0    44
MatSolve              93 2.0 5.5552e-05 1.7 4.01e+03 2.6 0.0e+00 0.0e+00 0.0e+00  0 14  0  0  0   0 14  0  0  0   576
MatLUFactorNum         6 1.0 6.9141e-05 2.0 3.31e+02 3.3 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0    31
MatILUFactorSym        6 1.0 1.9312e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 1.8e+01  0  0  0  0  3   0  0  0  0  3     0
MatAssemblyBegin      12 1.0 2.3067e-03 1.5 0.00e+00 0.0 1.2e+02 1.7e+01 2.4e+01  1  0  5  8  4   1  0  5  8  4     0
MatAssemblyEnd        12 1.0 1.2372e-03 1.0 0.00e+00 0.0 1.7e+02 4.9e+00 4.8e+01  1  0  7  3  9   1  0  7  3  9     0
MatGetRowIJ            6 3.0 1.8835e-0519.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         6 3.0 3.4404e-04 3.1 0.00e+00 0.0 0.0e+00 0.0e+00 1.2e+01  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries        18 3.0 2.0504e-05 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        81 1.0 1.3509e-03 1.2 1.57e+04 2.1 0.0e+00 0.0e+00 8.1e+01  1 64  0  0 14   1 64  0  0 15   106
KSPSetUp              12 1.0 2.2554e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               6 1.0 2.3280e-02 1.0 2.58e+04 2.2 1.4e+03 1.2e+01 2.2e+02 14100 57 58 39  14100 57 58 39    10
PCSetUp               12 1.0 2.2600e-03 1.1 3.31e+02 3.3 0.0e+00 0.0e+00 4.6e+01  1  1  0  0  8   1  1  0  0  8     1
PCSetUpOnBlocks        6 1.0 1.2162e-03 1.1 3.31e+02 3.3 0.0e+00 0.0e+00 3.4e+01  1  1  0  0  6   1  1  0  0  6     2
PCApply               93 1.0 8.8787e-04 1.1 4.01e+03 2.6 0.0e+00 0.0e+00 0.0e+00  1 14  0  0  0   1 14  0  0  0    36
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector   193            193       295224     0
      Vector Scatter    17             17        17612     0
           Index Set    55             55        41144     0
   IS L to G Mapping     6              6         3384     0
              Matrix    24             24        64344     0
       Krylov Solver    12             12       116160     0
      Preconditioner    12             12        10704     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 4.62532e-06
Average time for zero size MPI_Send(): 1.31528e-05
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
| Time:           Thu Jan 31 21:57:46 2013                                                                             |
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
| libMesh Performance: Alive time=0.314226, Active time=0.145099                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     6         0.0010      0.000161    0.0017      0.000277    0.67     1.14     |
|   build_sparsity()                 6         0.0020      0.000336    0.0081      0.001357    1.39     5.61     |
|   create_dof_constraints()         6         0.0000      0.000007    0.0000      0.000007    0.03     0.03     |
|   distribute_dofs()                6         0.0058      0.000973    0.0195      0.003255    4.02     13.46    |
|   dof_indices()                    62        0.0028      0.000044    0.0028      0.000044    1.90     1.90     |
|   old_dof_indices()                18        0.0008      0.000044    0.0008      0.000044    0.54     0.54     |
|   prepare_send_list()              6         0.0001      0.000009    0.0001      0.000009    0.04     0.04     |
|   reinit()                         6         0.0062      0.001040    0.0062      0.001040    4.30     4.30     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0004      0.000404    0.0010      0.000995    0.28     0.69     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        24        0.0002      0.000006    0.0002      0.000006    0.11     0.11     |
|   init_shape_functions()           19        0.0005      0.000029    0.0005      0.000029    0.38     0.38     |
|   inverse_map()                    22        0.0001      0.000006    0.0001      0.000006    0.09     0.09     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             24        0.0002      0.000010    0.0002      0.000010    0.17     0.17     |
|   compute_face_map()               7         0.0000      0.000006    0.0000      0.000006    0.03     0.03     |
|   init_face_shape_functions()      4         0.0000      0.000008    0.0000      0.000008    0.02     0.02     |
|   init_reference_to_physical_map() 19        0.0003      0.000015    0.0003      0.000015    0.20     0.20     |
|                                                                                                                |
| GnuPlotIO                                                                                                      |
|   write_nodal_data()               1         0.0019      0.001909    0.0019      0.001909    1.32     1.32     |
|                                                                                                                |
| JumpErrorEstimator                                                                                             |
|   estimate_error()                 5         0.0012      0.000246    0.0059      0.001182    0.85     4.07     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           76        0.0007      0.000009    0.0007      0.000009    0.46     0.46     |
|   init()                           10        0.0009      0.000092    0.0009      0.000092    0.63     0.63     |
|                                                                                                                |
| Mesh                                                                                                           |
|   contract()                       5         0.0002      0.000031    0.0003      0.000050    0.11     0.17     |
|   find_neighbors()                 6         0.0041      0.000678    0.0048      0.000806    2.80     3.33     |
|   renumber_nodes_and_elem()        17        0.0004      0.000024    0.0004      0.000024    0.28     0.28     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        7         0.0011      0.000161    0.0011      0.000161    0.78     0.78     |
|   find_global_indices()            7         0.0018      0.000263    0.0127      0.001811    1.27     8.73     |
|   parallel_sort()                  7         0.0047      0.000666    0.0059      0.000848    3.21     4.09     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0002      0.000214    0.0032      0.003242    0.15     2.23     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _coarsen_elements()              10        0.0002      0.000022    0.0003      0.000030    0.15     0.21     |
|   _refine_elements()               10        0.0023      0.000229    0.0056      0.000558    1.58     3.85     |
|   add_point()                      76        0.0006      0.000009    0.0013      0.000018    0.45     0.92     |
|   make_coarsening_compatible()     21        0.0039      0.000184    0.0039      0.000184    2.67     2.67     |
|   make_refinement_compatible()     21        0.0004      0.000019    0.0006      0.000027    0.28     0.39     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0005      0.000495    0.0005      0.000495    0.34     0.34     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      6         0.0102      0.001698    0.0209      0.003489    7.02     14.43    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      32        0.0005      0.000015    0.0007      0.000021    0.34     0.46     |
|   max(bool)                        52        0.0023      0.000044    0.0023      0.000044    1.57     1.57     |
|   max(scalar)                      921       0.0071      0.000008    0.0071      0.000008    4.91     4.91     |
|   max(vector)                      224       0.0032      0.000014    0.0084      0.000037    2.22     5.78     |
|   min(bool)                        1163      0.0085      0.000007    0.0085      0.000007    5.87     5.87     |
|   min(scalar)                      909       0.0133      0.000015    0.0133      0.000015    9.14     9.14     |
|   min(vector)                      224       0.0034      0.000015    0.0088      0.000039    2.32     6.04     |
|   probe()                          682       0.0028      0.000004    0.0028      0.000004    1.91     1.91     |
|   receive()                        682       0.0038      0.000006    0.0067      0.000010    2.61     4.61     |
|   send()                           682       0.0020      0.000003    0.0020      0.000003    1.37     1.37     |
|   send_receive()                   696       0.0050      0.000007    0.0153      0.000022    3.42     10.55    |
|   sum()                            34        0.0008      0.000023    0.0020      0.000057    0.54     1.35     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           682       0.0012      0.000002    0.0012      0.000002    0.85     0.85     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         6         0.0019      0.000314    0.0053      0.000884    1.30     3.66     |
|   set_parent_processor_ids()       6         0.0004      0.000073    0.0004      0.000073    0.30     0.30     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          6         0.0284      0.004732    0.0284      0.004732    19.57    19.57    |
|                                                                                                                |
| ProjectVector                                                                                                  |
|   operator()                       5         0.0005      0.000104    0.0014      0.000275    0.36     0.95     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       6         0.0008      0.000137    0.0020      0.000328    0.57     1.35     |
|   project_vector()                 5         0.0034      0.000684    0.0065      0.001300    2.36     4.48     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            7540      0.1451                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example adaptivity_ex1:
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
