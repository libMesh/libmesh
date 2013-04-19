<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("miscellaneous_ex1",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file miscellaneous_ex1.C with comments: </h1> 
<div class = "comment">
<h1>Miscellaneous Example 1 - Infinite Elements for the Wave Equation</h1>

<br><br>This is the sixth example program.  It builds on
the previous examples, and introduces the Infinite
Element class.  Note that the library must be compiled
with Infinite Elements enabled.  Otherwise, this
example will abort.
This example intends to demonstrate the similarities
between the \p FE and the \p InfFE classes in libMesh.
The matrices are assembled according to the wave equation.
However, for practical applications a time integration
scheme (as introduced in subsequent examples) should be
used.


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
Basic include file needed for the mesh functionality.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/exodusII_io.h"
        #include "libmesh/libmesh.h"
        #include "libmesh/serial_mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/linear_implicit_system.h"
        #include "libmesh/equation_systems.h"
        
</pre>
</div>
<div class = "comment">
Define the Finite and Infinite Element object.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/fe.h"
        #include "libmesh/inf_fe.h"
        #include "libmesh/inf_elem_builder.h"
        
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
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
</pre>
</div>
<div class = "comment">
Function prototype.  This is similar to the Poisson
assemble function of example 4.
</div>

<div class ="fragment">
<pre>
        void assemble_wave (EquationSystems& es,
                            const std::string& system_name);
        
</pre>
</div>
<div class = "comment">
Begin the main program.
</div>

<div class ="fragment">
<pre>
        int main (int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Initialize libMesh, like in example 2.
</div>

<div class ="fragment">
<pre>
          LibMeshInit init (argc, argv);
        
</pre>
</div>
<div class = "comment">
This example requires Infinite Elements
</div>

<div class ="fragment">
<pre>
        #ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS
          libmesh_example_assert(false, "--enable-ifem");
        #else
        
</pre>
</div>
<div class = "comment">
Skip this 3D example if libMesh was compiled as 1D/2D-only.
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(3 &lt;= LIBMESH_DIM, "3D support");
        
</pre>
</div>
<div class = "comment">
Tell the user what we are doing.
</div>

<div class ="fragment">
<pre>
          std::cout &lt;&lt; "Running ex6 with dim = 3" &lt;&lt; std::endl &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
Create a serialized mesh, distributed across the default MPI
communicator.
InfElemBuilder still requires some updates to be ParallelMesh
compatible


<br><br></div>

<div class ="fragment">
<pre>
          SerialMesh mesh(init.comm());
        
</pre>
</div>
<div class = "comment">
Use the internal mesh generator to create elements
on the square [-1,1]^3, of type Hex8.
</div>

<div class ="fragment">
<pre>
          MeshTools::Generation::build_cube (mesh,
                                             4, 4, 4,
                                             -1., 1.,
                                             -1., 1.,
                                             -1., 1.,
                                             HEX8);
        
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
Write the mesh before the infinite elements are added
</div>

<div class ="fragment">
<pre>
        #ifdef LIBMESH_HAVE_EXODUS_API
          ExodusII_IO(mesh).write ("orig_mesh.e");
        #endif
        
</pre>
</div>
<div class = "comment">
Normally, when a mesh is imported or created in
libMesh, only conventional elements exist.  The infinite
elements used here, however, require prescribed
nodal locations (with specified distances from an imaginary
origin) and configurations that a conventional mesh creator
in general does not offer.  Therefore, an efficient method
for building infinite elements is offered.  It can account
for symmetry planes and creates infinite elements in a fully
automatic way.

<br><br>Right now, the simplified interface is used, automatically
determining the origin.  Check \p MeshBase for a generalized
method that can even return the element faces of interior
vibrating surfaces.  The \p bool determines whether to be
verbose.
</div>

<div class ="fragment">
<pre>
          InfElemBuilder builder(mesh);
          builder.build_inf_elem(true);
        
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
Write the mesh with the infinite elements added.
Compare this to the original mesh.
</div>

<div class ="fragment">
<pre>
        #ifdef LIBMESH_HAVE_EXODUS_API
          ExodusII_IO(mesh).write ("ifems_added.e");
        #endif
        
</pre>
</div>
<div class = "comment">
After building infinite elements, we have to let
the elements find their neighbors again.
</div>

<div class ="fragment">
<pre>
          mesh.find_neighbors();
        
</pre>
</div>
<div class = "comment">
Create an equation systems object, where \p ThinSystem
offers only the crucial functionality for solving a
system.  Use \p ThinSystem when you want the sleekest
system possible.
</div>

<div class ="fragment">
<pre>
          EquationSystems equation_systems (mesh);
        
</pre>
</div>
<div class = "comment">
Declare the system and its variables.
Create a system named "Wave".  This can
be a simple, steady system
</div>

<div class ="fragment">
<pre>
          equation_systems.add_system&lt;LinearImplicitSystem&gt; ("Wave");
        
</pre>
</div>
<div class = "comment">
Create an FEType describing the approximation
characteristics of the InfFE object.  Note that
the constructor automatically defaults to some
sensible values.  But use \p FIRST order
approximation.
</div>

<div class ="fragment">
<pre>
          FEType fe_type(FIRST);
        
</pre>
</div>
<div class = "comment">
Add the variable "p" to "Wave".  Note that there exist
various approaches in adding variables.  In example 3,
\p add_variable took the order of approximation and used
default values for the \p FEFamily, while here the \p FEType
is used.
</div>

<div class ="fragment">
<pre>
          equation_systems.get_system("Wave").add_variable("p", fe_type);
        
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly
function.
</div>

<div class ="fragment">
<pre>
          equation_systems.get_system("Wave").attach_assemble_function (assemble_wave);
        
</pre>
</div>
<div class = "comment">
Set the speed of sound and fluid density
as \p EquationSystems parameter,
so that \p assemble_wave() can access it.
</div>

<div class ="fragment">
<pre>
          equation_systems.parameters.set&lt;Real&gt;("speed")          = 1.;
          equation_systems.parameters.set&lt;Real&gt;("fluid density")  = 1.;
        
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
Solve the system "Wave".
</div>

<div class ="fragment">
<pre>
          equation_systems.get_system("Wave").solve();
        
</pre>
</div>
<div class = "comment">
Write the whole EquationSystems object to file.
For infinite elements, the concept of nodal_soln()
is not applicable. Therefore, writing the mesh in
some format @e always gives all-zero results at
the nodes of the infinite elements.  Instead,
use the FEInterface::compute_data() methods to
determine physically correct results within an
infinite element.
</div>

<div class ="fragment">
<pre>
          equation_systems.write ("eqn_sys.dat", libMeshEnums::WRITE);
        
</pre>
</div>
<div class = "comment">
All done.
</div>

<div class ="fragment">
<pre>
          return 0;
        
        #endif // else part of ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS
        }
        
</pre>
</div>
<div class = "comment">
This function assembles the system matrix and right-hand-side
for the discrete form of our wave equation.
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
        
        
        #ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
        
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
Get a reference to the system we are solving.
</div>

<div class ="fragment">
<pre>
          LinearImplicitSystem & system = es.get_system&lt;LinearImplicitSystem&gt;("Wave");
        
</pre>
</div>
<div class = "comment">
A reference to the \p DofMap object for this system.  The \p DofMap
object handles the index translation from node and element numbers
to degree of freedom numbers.
</div>

<div class ="fragment">
<pre>
          const DofMap& dof_map = system.get_dof_map();
        
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
Copy the speed of sound to a local variable.
</div>

<div class ="fragment">
<pre>
          const Real speed = es.parameters.get&lt;Real&gt;("speed");
        
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
Build a Finite Element object of the specified type.  Since the
\p FEBase::build() member dynamically creates memory we will
store the object as an \p AutoPtr<FEBase>.  Check ex5 for details.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
        
</pre>
</div>
<div class = "comment">
Do the same for an infinite element.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;FEBase&gt; inf_fe (FEBase::build_InfFE(dim, fe_type));
        
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
Due to its internal structure, the infinite element handles
quadrature rules differently.  It takes the quadrature
rule which has been initialized for the FE object, but
creates suitable quadrature rules by @e itself.  The user
need not worry about this.
</div>

<div class ="fragment">
<pre>
          inf_fe-&gt;attach_quadrature_rule (&qrule);
        
</pre>
</div>
<div class = "comment">
Define data structures to contain the element matrix
and right-hand-side vector contribution.  Following
basic finite element terminology we will denote these
"Ke",  "Ce", "Me", and "Fe" for the stiffness, damping
and mass matrices, and the load vector.  Note that in
Acoustics, these descriptors though do @e not match the
true physical meaning of the projectors.  The final
overall system, however, resembles the conventional
notation again.
</div>

<div class ="fragment">
<pre>
          DenseMatrix&lt;Number&gt; Ke;
          DenseMatrix&lt;Number&gt; Ce;
          DenseMatrix&lt;Number&gt; Me;
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
</div>

<div class ="fragment">
<pre>
          MeshBase::const_element_iterator           el = mesh.active_local_elements_begin();
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
The mesh contains both finite and infinite elements.  These
elements are handled through different classes, namely
\p FE and \p InfFE, respectively.  However, since both
are derived from \p FEBase, they share the same interface,
and overall burden of coding is @e greatly reduced through
using a pointer, which is adjusted appropriately to the
current element type.
</div>

<div class ="fragment">
<pre>
              FEBase* cfe=NULL;
        
</pre>
</div>
<div class = "comment">
This here is almost the only place where we need to
distinguish between finite and infinite elements.
For faster computation, however, different approaches
may be feasible.

<br><br>Up to now, we do not know what kind of element we
have.  Aske the element of what type it is:
</div>

<div class ="fragment">
<pre>
              if (elem-&gt;infinite())
                {
</pre>
</div>
<div class = "comment">
We have an infinite element.  Let \p cfe point
to our \p InfFE object.  This is handled through
an AutoPtr.  Through the \p AutoPtr::get() we "borrow"
the pointer, while the \p  AutoPtr \p inf_fe is
still in charge of memory management.
</div>

<div class ="fragment">
<pre>
                  cfe = inf_fe.get();
                }
              else
                {
</pre>
</div>
<div class = "comment">
This is a conventional finite element.  Let \p fe handle it.
</div>

<div class ="fragment">
<pre>
                    cfe = fe.get();
        
</pre>
</div>
<div class = "comment">
Boundary conditions.
Here we just zero the rhs-vector. For natural boundary
conditions check e.g. previous examples.
</div>

<div class ="fragment">
<pre>
                  {
</pre>
</div>
<div class = "comment">
Zero the RHS for this element.
</div>

<div class ="fragment">
<pre>
                    Fe.resize (dof_indices.size());
        
                    system.rhs-&gt;add_vector (Fe, dof_indices);
                  } // end boundary condition section
                } // else ( if (elem-&gt;infinite())) )
        
</pre>
</div>
<div class = "comment">
This is slightly different from the Poisson solver:
Since the finite element object may change, we have to
initialize the constant references to the data fields
each time again, when a new element is processed.

<br><br>The element Jacobian * quadrature weight at each integration point.
</div>

<div class ="fragment">
<pre>
              const std::vector&lt;Real&gt;& JxW = cfe-&gt;get_JxW();
        
</pre>
</div>
<div class = "comment">
The element shape functions evaluated at the quadrature points.
</div>

<div class ="fragment">
<pre>
              const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi = cfe-&gt;get_phi();
        
</pre>
</div>
<div class = "comment">
The element shape function gradients evaluated at the quadrature
points.
</div>

<div class ="fragment">
<pre>
              const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi = cfe-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
The infinite elements need more data fields than conventional FE.
These are the gradients of the phase term \p dphase, an additional
radial weight for the test functions \p Sobolev_weight, and its
gradient.

<br><br>Note that these data fields are also initialized appropriately by
the \p FE method, so that the weak form (below) is valid for @e both
finite and infinite elements.
</div>

<div class ="fragment">
<pre>
              const std::vector&lt;RealGradient&gt;& dphase  = cfe-&gt;get_dphase();
              const std::vector&lt;Real&gt;&         weight  = cfe-&gt;get_Sobolev_weight();
              const std::vector&lt;RealGradient&gt;& dweight = cfe-&gt;get_Sobolev_dweight();
        
</pre>
</div>
<div class = "comment">
Now this is all independent of whether we use an \p FE
or an \p InfFE.  Nice, hm? ;-)

<br><br>Compute the element-specific data, as described
in previous examples.
</div>

<div class ="fragment">
<pre>
              cfe-&gt;reinit (elem);
        
</pre>
</div>
<div class = "comment">
Zero the element matrices.  Boundary conditions were already
processed in the \p FE-only section, see above.
</div>

<div class ="fragment">
<pre>
              Ke.resize (dof_indices.size(), dof_indices.size());
              Ce.resize (dof_indices.size(), dof_indices.size());
              Me.resize (dof_indices.size(), dof_indices.size());
        
</pre>
</div>
<div class = "comment">
The total number of quadrature points for infinite elements
@e has to be determined in a different way, compared to
conventional finite elements.  This type of access is also
valid for finite elements, so this can safely be used
anytime, instead of asking the quadrature rule, as
seen in previous examples.
</div>

<div class ="fragment">
<pre>
              unsigned int max_qp = cfe-&gt;n_quadrature_points();
        
</pre>
</div>
<div class = "comment">
Loop over the quadrature points.
</div>

<div class ="fragment">
<pre>
              for (unsigned int qp=0; qp&lt;max_qp; qp++)
                {
</pre>
</div>
<div class = "comment">
Similar to the modified access to the number of quadrature
points, the number of shape functions may also be obtained
in a different manner.  This offers the great advantage
of being valid for both finite and infinite elements.
</div>

<div class ="fragment">
<pre>
                  const unsigned int n_sf = cfe-&gt;n_shape_functions();
        
</pre>
</div>
<div class = "comment">
Now we will build the element matrices.  Since the infinite
elements are based on a Petrov-Galerkin scheme, the
resulting system matrices are non-symmetric. The additional
weight, described before, is part of the trial space.

<br><br>For the finite elements, though, these matrices are symmetric
just as we know them, since the additional fields \p dphase,
\p weight, and \p dweight are initialized appropriately.

<br><br>test functions:    weight[qp]*phi[i][qp]
trial functions:   phi[j][qp]
phase term:        phase[qp]

<br><br>derivatives are similar, but note that these are of type
Point, not of type Real.
</div>

<div class ="fragment">
<pre>
                  for (unsigned int i=0; i&lt;n_sf; i++)
                    for (unsigned int j=0; j&lt;n_sf; j++)
                      {
</pre>
</div>
<div class = "comment">
(ndt*Ht + nHt*d) * nH
</div>

<div class ="fragment">
<pre>
                        Ke(i,j) +=
                          (                            //    (
                           (                           //      (
                            dweight[qp] * phi[i][qp]   //        Point * Real  = Point
                            +                          //        +
                            dphi[i][qp] * weight[qp]   //        Point * Real  = Point
                            ) * dphi[j][qp]            //      )       * Point = Real
                           ) * JxW[qp];                //    )         * Real  = Real
        
</pre>
</div>
<div class = "comment">
(d*Ht*nmut*nH - ndt*nmu*Ht*H - d*nHt*nmu*H)
</div>

<div class ="fragment">
<pre>
                        Ce(i,j) +=
                          (                                //    (
                           (dphase[qp] * dphi[j][qp])      //      (Point * Point) = Real
                           * weight[qp] * phi[i][qp]       //      * Real * Real   = Real
                           -                               //      -
                           (dweight[qp] * dphase[qp])      //      (Point * Point) = Real
                           * phi[i][qp] * phi[j][qp]       //      * Real * Real   = Real
                           -                               //      -
                           (dphi[i][qp] * dphase[qp])      //      (Point * Point) = Real
                           * weight[qp] * phi[j][qp]       //      * Real * Real   = Real
                           ) * JxW[qp];                    //    )         * Real  = Real
        
</pre>
</div>
<div class = "comment">
(d*Ht*H * (1 - nmut*nmu))
</div>

<div class ="fragment">
<pre>
                        Me(i,j) +=
                          (                                       //    (
                           (1. - (dphase[qp] * dphase[qp]))       //      (Real  - (Point * Point)) = Real
                           * phi[i][qp] * phi[j][qp] * weight[qp] //      * Real *  Real  * Real    = Real
                           ) * JxW[qp];                           //    ) * Real                    = Real
        
                      } // end of the matrix summation loop
                } // end of quadrature point loop
        
</pre>
</div>
<div class = "comment">
The element matrices are now built for this element.
Collect them in Ke, and then add them to the global matrix.
The \p SparseMatrix::add_matrix() member does this for us.
</div>

<div class ="fragment">
<pre>
              Ke.add(1./speed        , Ce);
              Ke.add(1./(speed*speed), Me);
        
</pre>
</div>
<div class = "comment">
If this assembly program were to be used on an adaptive mesh,
we would have to apply any hanging node constraint equations
</div>

<div class ="fragment">
<pre>
              dof_map.constrain_element_matrix(Ke, dof_indices);
        
              system.matrix-&gt;add_matrix (Ke, dof_indices);
            } // end of element loop
        
</pre>
</div>
<div class = "comment">
Note that we have not applied any boundary conditions so far.
Here we apply a unit load at the node located at (0,0,0).
</div>

<div class ="fragment">
<pre>
          {
</pre>
</div>
<div class = "comment">
Iterate over local nodes
</div>

<div class ="fragment">
<pre>
            MeshBase::const_node_iterator           nd = mesh.local_nodes_begin();
            const MeshBase::const_node_iterator nd_end = mesh.local_nodes_end();
        
            for (; nd != nd_end; ++nd)
              {
</pre>
</div>
<div class = "comment">
Get a reference to the current node.
</div>

<div class ="fragment">
<pre>
                const Node& node = **nd;
        
</pre>
</div>
<div class = "comment">
Check the location of the current node.
</div>

<div class ="fragment">
<pre>
                if (fabs(node(0)) &lt; TOLERANCE &&
                    fabs(node(1)) &lt; TOLERANCE &&
                    fabs(node(2)) &lt; TOLERANCE)
                  {
</pre>
</div>
<div class = "comment">
The global number of the respective degree of freedom.
</div>

<div class ="fragment">
<pre>
                    unsigned int dn = node.dof_number(0,0,0);
        
                    system.rhs-&gt;add (dn, 1.);
                  }
              }
          }
        
        #else
        
</pre>
</div>
<div class = "comment">
dummy assert
</div>

<div class ="fragment">
<pre>
          libmesh_assert_not_equal_to (es.get_mesh().mesh_dimension(), 1);
        
        #endif //ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
        
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

<a name="nocomments"></a> 
<br><br><br> <h1> The source file miscellaneous_ex1.C without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/serial_mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/inf_fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/inf_elem_builder.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/node.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_wave (EquationSystems&amp; es,
                      <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
  #ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-ifem&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    libmesh_example_assert(3 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;3D support&quot;</FONT></B>);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Running ex6 with dim = 3&quot;</FONT></B> &lt;&lt; std::endl &lt;&lt; std::endl;
  
  
    SerialMesh mesh(init.comm());
  
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube (mesh,
                                       4, 4, 4,
                                       -1., 1.,
                                       -1., 1.,
                                       -1., 1.,
                                       HEX8);
  
    mesh.print_info();
  
  #ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO(mesh).write (<B><FONT COLOR="#BC8F8F">&quot;orig_mesh.e&quot;</FONT></B>);
  #endif
  
    InfElemBuilder builder(mesh);
    builder.build_inf_elem(true);
  
    mesh.print_info();
  
  #ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO(mesh).write (<B><FONT COLOR="#BC8F8F">&quot;ifems_added.e&quot;</FONT></B>);
  #endif
  
    mesh.find_neighbors();
  
    EquationSystems equation_systems (mesh);
  
    equation_systems.add_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Wave&quot;</FONT></B>);
  
    FEType fe_type(FIRST);
  
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Wave&quot;</FONT></B>).add_variable(<B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>, fe_type);
  
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Wave&quot;</FONT></B>).attach_assemble_function (assemble_wave);
  
    equation_systems.parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;speed&quot;</FONT></B>)          = 1.;
    equation_systems.parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;fluid density&quot;</FONT></B>)  = 1.;
  
    equation_systems.init();
  
    equation_systems.print_info();
  
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Wave&quot;</FONT></B>).solve();
  
    equation_systems.write (<B><FONT COLOR="#BC8F8F">&quot;eqn_sys.dat&quot;</FONT></B>, libMeshEnums::WRITE);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  
  #endif <I><FONT COLOR="#B22222">// else part of ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS
</FONT></I>  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_wave(EquationSystems&amp; es,
                     <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Wave&quot;</FONT></B>);
  
  
  #ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
  
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    LinearImplicitSystem &amp; system = es.get_system&lt;LinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Wave&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    <B><FONT COLOR="#228B22">const</FONT></B> Real speed = es.parameters.get&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;speed&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> FEType&amp; fe_type = dof_map.variable_type(0);
  
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
  
    AutoPtr&lt;FEBase&gt; inf_fe (FEBase::build_InfFE(dim, fe_type));
  
    QGauss qrule (dim, SECOND);
  
    fe-&gt;attach_quadrature_rule (&amp;qrule);
  
    inf_fe-&gt;attach_quadrature_rule (&amp;qrule);
  
    DenseMatrix&lt;Number&gt; Ke;
    DenseMatrix&lt;Number&gt; Ce;
    DenseMatrix&lt;Number&gt; Me;
    DenseVector&lt;Number&gt; Fe;
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator           el = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
  
  
        FEBase* cfe=NULL;
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;infinite())
          {
            cfe = inf_fe.get();
          }
        <B><FONT COLOR="#A020F0">else</FONT></B>
          {
              cfe = fe.get();
  
            {
              Fe.resize (dof_indices.size());
  
              system.rhs-&gt;add_vector (Fe, dof_indices);
            } <I><FONT COLOR="#B22222">// end boundary condition section
</FONT></I>          } <I><FONT COLOR="#B22222">// else ( if (elem-&gt;infinite())) )
</FONT></I>  
        <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = cfe-&gt;get_JxW();
  
        <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = cfe-&gt;get_phi();
  
        <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = cfe-&gt;get_dphi();
  
        <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;RealGradient&gt;&amp; dphase  = cfe-&gt;get_dphase();
        <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp;         weight  = cfe-&gt;get_Sobolev_weight();
        <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;RealGradient&gt;&amp; dweight = cfe-&gt;get_Sobolev_dweight();
  
        cfe-&gt;reinit (elem);
  
        Ke.resize (dof_indices.size(), dof_indices.size());
        Ce.resize (dof_indices.size(), dof_indices.size());
        Me.resize (dof_indices.size(), dof_indices.size());
  
        <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_qp = cfe-&gt;n_quadrature_points();
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;max_qp; qp++)
          {
            <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_sf = cfe-&gt;n_shape_functions();
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_sf; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_sf; j++)
                {
                  Ke(i,j) +=
                    (                            <I><FONT COLOR="#B22222">//    (
</FONT></I>                     (                           <I><FONT COLOR="#B22222">//      (
</FONT></I>                      dweight[qp] * phi[i][qp]   <I><FONT COLOR="#B22222">//        Point * Real  = Point
</FONT></I>                      +                          <I><FONT COLOR="#B22222">//        +
</FONT></I>                      dphi[i][qp] * weight[qp]   <I><FONT COLOR="#B22222">//        Point * Real  = Point
</FONT></I>                      ) * dphi[j][qp]            <I><FONT COLOR="#B22222">//      )       * Point = Real
</FONT></I>                     ) * JxW[qp];                <I><FONT COLOR="#B22222">//    )         * Real  = Real
</FONT></I>  
                  Ce(i,j) +=
                    (                                <I><FONT COLOR="#B22222">//    (
</FONT></I>                     (dphase[qp] * dphi[j][qp])      <I><FONT COLOR="#B22222">//      (Point * Point) = Real
</FONT></I>                     * weight[qp] * phi[i][qp]       <I><FONT COLOR="#B22222">//      * Real * Real   = Real
</FONT></I>                     -                               <I><FONT COLOR="#B22222">//      -
</FONT></I>                     (dweight[qp] * dphase[qp])      <I><FONT COLOR="#B22222">//      (Point * Point) = Real
</FONT></I>                     * phi[i][qp] * phi[j][qp]       <I><FONT COLOR="#B22222">//      * Real * Real   = Real
</FONT></I>                     -                               <I><FONT COLOR="#B22222">//      -
</FONT></I>                     (dphi[i][qp] * dphase[qp])      <I><FONT COLOR="#B22222">//      (Point * Point) = Real
</FONT></I>                     * weight[qp] * phi[j][qp]       <I><FONT COLOR="#B22222">//      * Real * Real   = Real
</FONT></I>                     ) * JxW[qp];                    <I><FONT COLOR="#B22222">//    )         * Real  = Real
</FONT></I>  
                  Me(i,j) +=
                    (                                       <I><FONT COLOR="#B22222">//    (
</FONT></I>                     (1. - (dphase[qp] * dphase[qp]))       <I><FONT COLOR="#B22222">//      (Real  - (Point * Point)) = Real
</FONT></I>                     * phi[i][qp] * phi[j][qp] * weight[qp] <I><FONT COLOR="#B22222">//      * Real *  Real  * Real    = Real
</FONT></I>                     ) * JxW[qp];                           <I><FONT COLOR="#B22222">//    ) * Real                    = Real
</FONT></I>  
                } <I><FONT COLOR="#B22222">// end of the matrix summation loop
</FONT></I>          } <I><FONT COLOR="#B22222">// end of quadrature point loop
</FONT></I>  
        Ke.add(1./speed        , Ce);
        Ke.add(1./(speed*speed), Me);
  
        dof_map.constrain_element_matrix(Ke, dof_indices);
  
        system.matrix-&gt;add_matrix (Ke, dof_indices);
      } <I><FONT COLOR="#B22222">// end of element loop
</FONT></I>  
    {
      <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_node_iterator           nd = mesh.local_nodes_begin();
      <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_node_iterator nd_end = mesh.local_nodes_end();
  
      <B><FONT COLOR="#A020F0">for</FONT></B> (; nd != nd_end; ++nd)
        {
          <B><FONT COLOR="#228B22">const</FONT></B> Node&amp; node = **nd;
  
          <B><FONT COLOR="#A020F0">if</FONT></B> (fabs(node(0)) &lt; TOLERANCE &amp;&amp;
              fabs(node(1)) &lt; TOLERANCE &amp;&amp;
              fabs(node(2)) &lt; TOLERANCE)
            {
              <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dn = node.dof_number(0,0,0);
  
              system.rhs-&gt;add (dn, 1.);
            }
        }
    }
  
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    libmesh_assert_not_equal_to (es.get_mesh().mesh_dimension(), 1);
  
  #endif <I><FONT COLOR="#B22222">//ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B>;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
make[4]: Entering directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/miscellaneous/miscellaneous_ex1'
***************************************************************
* Running Example miscellaneous_ex1:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
 
Running ex6 with dim = 3

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=125
    n_local_nodes()=45
  n_elem()=64
    n_local_elem()=16
    n_active_elem()=64
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
  n_threads()=1
  processor_id()=0

 Verbose mode disabled in non-debug mode.
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=125
    n_local_nodes()=45
  n_elem()=64
    n_local_elem()=16
    n_active_elem()=64
  n_subdomains()=1
  n_partitions()=4
  n_processors()=4
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Wave"
    Type "LinearImplicit"
    Variables="p" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=125
    n_local_dofs()=45
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 14.512
      Average Off-Processor Bandwidth <= 6.544
      Maximum  On-Processor Bandwidth <= 33
      Maximum Off-Processor Bandwidth <= 24
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0


 -------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                 |
| Num Processors: 4                                                                                                 |
| Time:           Fri Apr 19 11:50:37 2013                                                                          |
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
| libMesh Performance: Alive time=0.176819, Active time=0.026243                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0002      0.000214    0.0006      0.000597    0.82     2.27     |
|   build_sparsity()                 1         0.0002      0.000222    0.0006      0.000589    0.85     2.24     |
|   create_dof_constraints()         1         0.0001      0.000058    0.0001      0.000058    0.22     0.22     |
|   distribute_dofs()                1         0.0004      0.000364    0.0011      0.001130    1.39     4.31     |
|   dof_indices()                    81        0.0006      0.000008    0.0006      0.000008    2.43     2.43     |
|   prepare_send_list()              1         0.0000      0.000011    0.0000      0.000011    0.04     0.04     |
|   reinit()                         1         0.0005      0.000478    0.0005      0.000478    1.82     1.82     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   write()                          1         0.0009      0.000901    0.0012      0.001185    3.43     4.52     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        16        0.0002      0.000013    0.0002      0.000013    0.77     0.77     |
|   init_shape_functions()           1         0.0000      0.000041    0.0000      0.000041    0.16     0.16     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             16        0.0001      0.000007    0.0001      0.000007    0.44     0.44     |
|   init_reference_to_physical_map() 1         0.0001      0.000054    0.0001      0.000054    0.21     0.21     |
|                                                                                                                |
| InfElemBuilder                                                                                                 |
|   build_inf_elem()                 1         0.0001      0.000078    0.0001      0.000078    0.30     0.30     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 4         0.0009      0.000229    0.0012      0.000291    3.48     4.44     |
|   renumber_nodes_and_elem()        4         0.0001      0.000016    0.0001      0.000016    0.25     0.25     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0021      0.002119    0.0025      0.002489    8.07     9.48     |
|   compute_hilbert_indices()        3         0.0008      0.000251    0.0008      0.000251    2.87     2.87     |
|   find_global_indices()            3         0.0002      0.000070    0.0018      0.000616    0.80     7.04     |
|   parallel_sort()                  3         0.0004      0.000137    0.0005      0.000172    1.56     1.97     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0002      0.000185    0.0002      0.000185    0.70     0.70     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      2         0.0023      0.001148    0.0034      0.001696    8.75     12.93    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      18        0.0001      0.000008    0.0002      0.000011    0.54     0.77     |
|   broadcast()                      2         0.0000      0.000006    0.0000      0.000006    0.05     0.05     |
|   max(bool)                        1         0.0000      0.000008    0.0000      0.000008    0.03     0.03     |
|   max(scalar)                      176       0.0007      0.000004    0.0007      0.000004    2.59     2.59     |
|   max(vector)                      42        0.0003      0.000006    0.0007      0.000016    1.00     2.61     |
|   min(bool)                        204       0.0007      0.000003    0.0007      0.000003    2.51     2.51     |
|   min(scalar)                      167       0.0009      0.000005    0.0009      0.000005    3.46     3.46     |
|   min(vector)                      42        0.0003      0.000008    0.0008      0.000019    1.27     2.97     |
|   probe()                          92        0.0002      0.000003    0.0002      0.000003    0.88     0.88     |
|   receive()                        76        0.0002      0.000003    0.0004      0.000006    0.87     1.64     |
|   send()                           64        0.0001      0.000002    0.0001      0.000002    0.50     0.50     |
|   send_receive()                   70        0.0003      0.000004    0.0009      0.000013    1.19     3.33     |
|   sum()                            26        0.0002      0.000006    0.0002      0.000008    0.60     0.83     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           64        0.0001      0.000001    0.0001      0.000001    0.32     0.32     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         2         0.0002      0.000121    0.0005      0.000262    0.92     2.00     |
|   set_parent_processor_ids()       2         0.0000      0.000025    0.0000      0.000025    0.19     0.19     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.0100      0.010037    0.0100      0.010037    38.25    38.25    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0014      0.001433    0.0020      0.001983    5.46     7.56     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1194      0.0262                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example miscellaneous_ex1:
*  mpirun -np 4 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc
***************************************************************
make[4]: Leaving directory `/net/spark/workspace/roystgnr/libmesh/git/devel/examples/miscellaneous/miscellaneous_ex1'
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
