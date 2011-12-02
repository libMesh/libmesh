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
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example  ./ex0-dbg
***************************************************************
 

 ---------------------------------------------------------------------------- 
| Reference count information                                                |
 ---------------------------------------------------------------------------- 
| N7libMesh10Parameters5ValueE reference count information:
|  Creations:    2
|  Destructions: 2
| N7libMesh12LinearSolverIdEE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh12SparseMatrixIdEE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh13NumericVectorIdEE reference count information:
|  Creations:    18
|  Destructions: 18
| N7libMesh15EquationSystemsE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh4ElemE reference count information:
|  Creations:    544
|  Destructions: 544
| N7libMesh4NodeE reference count information:
|  Creations:    85
|  Destructions: 85
| N7libMesh5QBaseE reference count information:
|  Creations:    26
|  Destructions: 26
| N7libMesh6DofMapE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh6FEBaseE reference count information:
|  Creations:    32
|  Destructions: 32
| N7libMesh6SystemE reference count information:
|  Creations:    1
|  Destructions: 1
| N7libMesh9DofObjectE reference count information:
|  Creations:    862
|  Destructions: 862
 ---------------------------------------------------------------------------- 

-------------------------------------------------------------------
| Time:           Thu Dec  1 21:20:45 2011                         |
| OS:             Linux                                            |
| HostName:       dknez-laptop                                     |
| OS Release:     3.0.0-13-generic                                 |
| OS Version:     #22-Ubuntu SMP Wed Nov 2 13:27:26 UTC 2011       |
| Machine:        x86_64                                           |
| Username:       dknez                                            |
| Configuration:  ./configure run on Tue Nov 22 16:02:41 EST 2011  |
-------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.073559, Active time=0.059951                                             |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 6         0.0009      0.000148    0.0009      0.000148    1.48     1.48     |
|   compute_sparsity()           6         0.0040      0.000665    0.0052      0.000862    6.66     8.63     |
|   create_dof_constraints()     6         0.0000      0.000003    0.0000      0.000003    0.03     0.03     |
|   distribute_dofs()            6         0.0051      0.000846    0.0110      0.001829    8.47     18.30    |
|   dof_indices()                464       0.0045      0.000010    0.0045      0.000010    7.58     7.58     |
|   old_dof_indices()            200       0.0018      0.000009    0.0018      0.000009    3.01     3.01     |
|   prepare_send_list()          6         0.0000      0.000006    0.0000      0.000006    0.06     0.06     |
|   reinit()                     6         0.0059      0.000980    0.0059      0.000980    9.80     9.80     |
|                                                                                                            |
| EquationSystems                                                                                            |
|   build_solution_vector()      1         0.0004      0.000350    0.0007      0.000727    0.58     1.21     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         218       0.0011      0.000005    0.0011      0.000005    1.85     1.85     |
|   compute_face_map()           57        0.0002      0.000004    0.0002      0.000004    0.37     0.37     |
|   compute_shape_functions()    218       0.0009      0.000004    0.0009      0.000004    1.56     1.56     |
|   init_face_shape_functions()  5         0.0000      0.000006    0.0000      0.000006    0.05     0.05     |
|   init_shape_functions()       120       0.0014      0.000012    0.0014      0.000012    2.36     2.36     |
|   inverse_map()                210       0.0009      0.000004    0.0009      0.000004    1.51     1.51     |
|                                                                                                            |
| GnuPlotIO                                                                                                  |
|   write_nodal_data()           1         0.0009      0.000943    0.0009      0.000943    1.57     1.57     |
|                                                                                                            |
| JumpErrorEstimator                                                                                         |
|   estimate_error()             5         0.0034      0.000689    0.0074      0.001473    5.75     12.29    |
|                                                                                                            |
| LocationMap                                                                                                |
|   find()                       76        0.0004      0.000005    0.0004      0.000005    0.68     0.68     |
|   init()                       10        0.0006      0.000057    0.0006      0.000057    0.95     0.95     |
|                                                                                                            |
| Mesh                                                                                                       |
|   contract()                   5         0.0002      0.000045    0.0003      0.000067    0.37     0.56     |
|   find_neighbors()             6         0.0037      0.000622    0.0037      0.000622    6.23     6.23     |
|   renumber_nodes_and_elem()    17        0.0003      0.000020    0.0003      0.000020    0.58     0.58     |
|                                                                                                            |
| MeshOutput                                                                                                 |
|   write_equation_systems()     1         0.0000      0.000042    0.0017      0.001713    0.07     2.86     |
|                                                                                                            |
| MeshRefinement                                                                                             |
|   _coarsen_elements()          10        0.0002      0.000015    0.0002      0.000015    0.25     0.25     |
|   _refine_elements()           10        0.0015      0.000149    0.0026      0.000256    2.49     4.27     |
|   add_point()                  76        0.0006      0.000008    0.0010      0.000014    1.01     1.74     |
|   make_coarsening_compatible() 21        0.0033      0.000158    0.0033      0.000158    5.54     5.54     |
|   make_refinement_compatible() 21        0.0003      0.000014    0.0003      0.000014    0.48     0.48     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0001      0.000095    0.0001      0.000095    0.16     0.16     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  6         0.0000      0.000002    0.0000      0.000002    0.03     0.03     |
|                                                                                                            |
| Partitioner                                                                                                |
|   single_partition()           6         0.0002      0.000034    0.0002      0.000034    0.35     0.35     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      6         0.0085      0.001419    0.0085      0.001419    14.20    14.20    |
|                                                                                                            |
| ProjectVector                                                                                              |
|   operator()                   5         0.0021      0.000414    0.0047      0.000935    3.45     7.80     |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   6         0.0026      0.000426    0.0050      0.000841    4.26     8.42     |
|   project_vector()             5         0.0037      0.000746    0.0093      0.001866    6.23     15.56    |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        1823      0.0600                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./ex0-dbg
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
