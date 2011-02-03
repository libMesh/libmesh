<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex0",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 0 - Solving 1D PDE Using Adaptive Mesh Refinement</h1>

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
        #include "mesh.h"
        #include "mesh_generation.h"
        #include "edge_edge3.h"
        #include "gnuplot_io.h"
        #include "equation_systems.h"
        #include "linear_implicit_system.h"
        #include "fe.h"
        #include "quadrature_gauss.h"
        #include "sparse_matrix.h"
        #include "dof_map.h"
        #include "numeric_vector.h"
        #include "dense_matrix.h"
        #include "dense_vector.h"
        #include "error_vector.h"
        #include "kelly_error_estimator.h"
        #include "mesh_refinement.h"
        
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
may depend on a number of other libraries (i.e. MPI  and Petsc)
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
        
</pre>
</div>
<div class = "comment">
Build a 1D mesh with 4 elements from x=0 to x=1, using 
EDGE3 (i.e. quadratic) 1D elements. They are called EDGE3 elements
because a quadratic element contains 3 nodes.
</div>

<div class ="fragment">
<pre>
          MeshTools::Generation::build_line(mesh,4,0.,1.,EDGE3);
        
</pre>
</div>
<div class = "comment">
Define the equation systems object and the system we are going
to solve. See Example 2 for more details.
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
Objects for error estimation, see Example 10 for more details.
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
          GnuPlotIO plot(mesh,"Example 0", GnuPlotIO::GRID_ON);
        
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
          libmesh_assert(system_name == "1D");
        
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
          std::vector&lt;unsigned int&gt; dof_indices;
        
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
mesh refinement. See Example 10 for more details.
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
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;edge_edge3.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;gnuplot_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;kelly_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_refinement.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_1D(EquationSystems&amp; es, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main(<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {   
    LibMeshInit init (argc, argv);
  
  #ifndef LIBMESH_ENABLE_AMR
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-amr&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    Mesh mesh;
  
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_line(mesh,4,0.,1.,EDGE3);
  
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
  
    GnuPlotIO plot(mesh,<B><FONT COLOR="#BC8F8F">&quot;Example 0&quot;</FONT></B>, GnuPlotIO::GRID_ON);
  
    plot.write_equation_systems(<B><FONT COLOR="#BC8F8F">&quot;gnuplot_script&quot;</FONT></B>,equation_systems);
  #endif <I><FONT COLOR="#B22222">// #ifndef LIBMESH_ENABLE_AMR
</FONT></I>    
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_1D(EquationSystems&amp; es, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
  
  #ifdef LIBMESH_ENABLE_AMR
  
    libmesh_assert(system_name == <B><FONT COLOR="#BC8F8F">&quot;1D&quot;</FONT></B>);
  
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
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
  
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
* Running Example  mpirun -np 2 ./ex0-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex0-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Thu Feb  3 12:09:28 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.357e-01      1.00438   1.355e-01
Objects:              3.120e+02      1.00000   3.120e+02
Flops:                1.496e+04      1.11489   1.419e+04  2.839e+04
Flops/sec:            1.102e+05      1.11002   1.048e+05  2.095e+05
MPI Messages:         1.585e+02      1.00000   1.585e+02  3.170e+02
MPI Message Lengths:  4.178e+03      1.03160   2.596e+01  8.228e+03
MPI Reductions:       4.770e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.3543e-01 100.0%  2.8386e+04 100.0%  3.170e+02 100.0%  2.596e+01      100.0%  3.660e+02  76.7% 

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

VecMDot               18 1.0 5.5313e-05 1.2 1.28e+03 1.1 0.0e+00 0.0e+00 1.8e+01  0  9  0  0  4   0  9  0  0  5    45
VecNorm               30 1.0 8.5115e-05 1.1 1.10e+03 1.1 0.0e+00 0.0e+00 3.0e+01  0  8  0  0  6   0  8  0  0  8    25
VecScale              24 1.0 1.9312e-05 1.1 4.40e+02 1.1 0.0e+00 0.0e+00 0.0e+00  0  3  0  0  0   0  3  0  0  0    44
VecCopy               23 1.0 9.0599e-06 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                64 1.0 1.5736e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               12 1.0 4.1027e-03 1.0 4.40e+02 1.1 0.0e+00 0.0e+00 0.0e+00  3  3  0  0  0   3  3  0  0  0     0
VecMAXPY              24 1.0 7.1526e-06 1.4 1.98e+03 1.1 0.0e+00 0.0e+00 0.0e+00  0 14  0  0  0   0 14  0  0  0   539
VecAssemblyBegin      48 1.0 2.6488e-04 1.0 0.00e+00 0.0 2.2e+01 6.0e+00 1.1e+02  0  0  7  2 24   0  0  7  2 31     0
VecAssemblyEnd        48 1.0 3.1948e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       61 1.0 7.4387e-05 1.1 0.00e+00 0.0 8.7e+01 2.1e+01 0.0e+00  0  0 27 23  0   0  0 27 23  0     0
VecScatterEnd         61 1.0 4.6015e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize          24 1.0 1.1277e-04 1.1 1.32e+03 1.1 0.0e+00 0.0e+00 2.4e+01  0  9  0  0  5   0  9  0  0  7    23
MatMult               24 1.0 1.0777e-04 1.0 3.03e+03 1.1 4.8e+01 1.2e+01 0.0e+00  0 21 15  7  0   0 21 15  7  0    54
MatSolve              24 1.0 2.2411e-05 1.4 4.90e+03 1.1 0.0e+00 0.0e+00 0.0e+00  0 32  0  0  0   0 32  0  0  0   409
MatLUFactorNum         6 1.0 3.6001e-05 1.0 1.78e+03 1.3 0.0e+00 0.0e+00 0.0e+00  0 11  0  0  0   0 11  0  0  0    88
MatILUFactorSym        6 1.0 9.1076e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 6.0e+00  0  0  0  0  1   0  0  0  0  2     0
MatAssemblyBegin      12 1.0 2.0242e-04 1.0 0.00e+00 0.0 1.8e+01 1.7e+01 2.4e+01  0  0  6  4  5   0  0  6  4  7     0
MatAssemblyEnd        12 1.0 4.1795e-04 1.0 0.00e+00 0.0 2.4e+01 5.0e+00 4.8e+01  0  0  8  1 10   0  0  8  1 13     0
MatGetRowIJ            6 1.0 2.1458e-06 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         6 1.0 7.8440e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.6e+01  0  0  0  0  3   0  0  0  0  4     0
MatZeroEntries        18 1.0 1.1444e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        18 1.0 8.6069e-05 1.2 2.60e+03 1.1 0.0e+00 0.0e+00 1.8e+01  0 18  0  0  4   0 18  0  0  5    59
KSPSetup              12 1.0 9.2030e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               6 1.0 5.5852e-03 1.0 1.50e+04 1.1 4.8e+01 1.2e+01 7.0e+01  4100 15  7 15   4100 15  7 19     5
PCSetUp               12 1.0 6.6161e-04 1.0 1.78e+03 1.3 0.0e+00 0.0e+00 2.2e+01  0 11  0  0  5   0 11  0  0  6     5
PCSetUpOnBlocks        6 1.0 3.3689e-04 1.0 1.78e+03 1.3 0.0e+00 0.0e+00 2.2e+01  0 11  0  0  5   0 11  0  0  6     9
PCApply               24 1.0 1.5950e-04 1.0 4.90e+03 1.1 0.0e+00 0.0e+00 0.0e+00  0 32  0  0  0   0 32  0  0  0    58
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec   158            158       230256     0
         Vec Scatter    33             33        28644     0
           Index Set    67             67        36304     0
   IS L to G Mapping     6              6         2912     0
              Matrix    24             24        69472     0
       Krylov Solver    12             12       113280     0
      Preconditioner    12             12         8448     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 1.38283e-06
Average time for zero size MPI_Send(): 5.00679e-06
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
| Time:           Thu Feb  3 12:09:28 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.288514, Active time=0.012526                                             |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 6         0.0001      0.000011    0.0001      0.000012    0.53     0.59     |
|   compute_sparsity()           6         0.0002      0.000039    0.0003      0.000045    1.88     2.16     |
|   create_dof_constraints()     6         0.0000      0.000000    0.0000      0.000000    0.01     0.01     |
|   distribute_dofs()            6         0.0002      0.000039    0.0006      0.000106    1.86     5.06     |
|   dof_indices()                299       0.0001      0.000000    0.0001      0.000000    0.89     0.89     |
|   old_dof_indices()            100       0.0000      0.000000    0.0000      0.000000    0.27     0.27     |
|   prepare_send_list()          6         0.0000      0.000001    0.0000      0.000001    0.03     0.03     |
|   reinit()                     6         0.0003      0.000043    0.0003      0.000043    2.05     2.05     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         112       0.0000      0.000000    0.0000      0.000000    0.30     0.30     |
|   compute_face_map()           30        0.0000      0.000000    0.0000      0.000000    0.10     0.10     |
|   compute_shape_functions()    112       0.0000      0.000000    0.0000      0.000000    0.18     0.18     |
|   init_face_shape_functions()  9         0.0000      0.000001    0.0000      0.000001    0.07     0.07     |
|   init_shape_functions()       66        0.0002      0.000003    0.0002      0.000003    1.48     1.48     |
|   inverse_map()                137       0.0001      0.000000    0.0001      0.000000    0.48     0.48     |
|                                                                                                            |
| JumpErrorEstimator                                                                                         |
|   estimate_error()             5         0.0004      0.000071    0.0006      0.000113    2.83     4.52     |
|                                                                                                            |
| LocationMap                                                                                                |
|   find()                       76        0.0001      0.000001    0.0001      0.000001    0.45     0.45     |
|   init()                       10        0.0001      0.000009    0.0001      0.000009    0.74     0.74     |
|                                                                                                            |
| Mesh                                                                                                       |
|   contract()                   5         0.0000      0.000004    0.0000      0.000005    0.15     0.20     |
|   find_neighbors()             6         0.0003      0.000050    0.0003      0.000054    2.39     2.59     |
|   renumber_nodes_and_elem()    17        0.0000      0.000001    0.0000      0.000001    0.12     0.12     |
|                                                                                                            |
| MeshCommunication                                                                                          |
|   compute_hilbert_indices()    7         0.0003      0.000045    0.0003      0.000045    2.49     2.49     |
|   find_global_indices()        7         0.0002      0.000023    0.0011      0.000163    1.30     9.09     |
|   parallel_sort()              7         0.0004      0.000054    0.0005      0.000067    2.99     3.73     |
|                                                                                                            |
| MeshRefinement                                                                                             |
|   _coarsen_elements()          10        0.0000      0.000003    0.0000      0.000004    0.22     0.32     |
|   _refine_elements()           10        0.0002      0.000019    0.0003      0.000034    1.52     2.75     |
|   add_point()                  76        0.0001      0.000001    0.0001      0.000002    0.57     1.05     |
|   make_coarsening_compatible() 11        0.0002      0.000022    0.0002      0.000022    1.89     1.89     |
|   make_refinement_compatible() 11        0.0000      0.000003    0.0000      0.000004    0.26     0.38     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0000      0.000037    0.0000      0.000037    0.30     0.30     |
|                                                                                                            |
| MetisPartitioner                                                                                           |
|   partition()                  6         0.0004      0.000066    0.0011      0.000175    3.17     8.40     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  25        0.0001      0.000003    0.0001      0.000003    0.55     0.55     |
|   max(bool)                    36        0.0001      0.000001    0.0001      0.000001    0.43     0.43     |
|   max(scalar)                  22        0.0001      0.000003    0.0001      0.000003    0.53     0.53     |
|   max(vector)                  7         0.0000      0.000002    0.0000      0.000002    0.10     0.10     |
|   min(bool)                    22        0.0000      0.000002    0.0000      0.000002    0.31     0.31     |
|   min(scalar)                  5         0.0000      0.000001    0.0000      0.000001    0.06     0.06     |
|   min(vector)                  7         0.0000      0.000003    0.0000      0.000003    0.18     0.18     |
|   probe()                      50        0.0001      0.000002    0.0001      0.000002    0.61     0.61     |
|   receive()                    50        0.0001      0.000002    0.0002      0.000004    0.81     1.44     |
|   send()                       50        0.0000      0.000001    0.0000      0.000001    0.21     0.21     |
|   send_receive()               64        0.0001      0.000002    0.0003      0.000005    0.85     2.79     |
|   sum()                        29        0.0001      0.000005    0.0001      0.000005    1.06     1.06     |
|   wait()                       50        0.0000      0.000001    0.0000      0.000001    0.20     0.20     |
|                                                                                                            |
| Partitioner                                                                                                |
|   set_node_processor_ids()     6         0.0001      0.000023    0.0002      0.000038    1.13     1.83     |
|   set_parent_processor_ids()   6         0.0000      0.000003    0.0000      0.000003    0.15     0.15     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      6         0.0068      0.001133    0.0068      0.001133    54.27    54.27    |
|                                                                                                            |
| ProjectVector                                                                                              |
|   operator()                   5         0.0002      0.000037    0.0003      0.000051    1.48     2.04     |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   6         0.0003      0.000046    0.0005      0.000078    2.20     3.75     |
|   project_vector()             5         0.0004      0.000084    0.0007      0.000142    3.35     5.65     |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        1617      0.0125                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 2 ./ex0-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
