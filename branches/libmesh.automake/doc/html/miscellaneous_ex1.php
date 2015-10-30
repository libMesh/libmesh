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
        #include "exodusII_io.h"
        #include "libmesh.h"
        #include "serial_mesh.h"
        #include "mesh_generation.h"
        #include "linear_implicit_system.h"
        #include "equation_systems.h"
        
</pre>
</div>
<div class = "comment">
Define the Finite and Infinite Element object.
</div>

<div class ="fragment">
<pre>
        #include "fe.h"
        #include "inf_fe.h"
        #include "inf_elem_builder.h"
        
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
Create a serialized mesh.
InfElemBuilder still requires some updates to be ParallelMesh
compatible
  

<br><br></div>

<div class ="fragment">
<pre>
          SerialMesh mesh;
        
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
          libmesh_assert (system_name == "Wave");
        
        
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
          libmesh_assert(es.get_mesh().mesh_dimension() != 1);
        
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
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;serial_mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;inf_fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;inf_elem_builder.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;node.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  
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
    
    
    SerialMesh mesh;
  
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
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Wave&quot;</FONT></B>);
  
  
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
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
    
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
  
    libmesh_assert(es.get_mesh().mesh_dimension() != 1);
  
  #endif <I><FONT COLOR="#B22222">//ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
</FONT></I>    
    <B><FONT COLOR="#A020F0">return</FONT></B>;
  }
  
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Linking miscellaneous_ex1-opt...
***************************************************************
* Running Example  mpirun -np 6 ./miscellaneous_ex1-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running ex6 with dim = 3

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=125
    n_local_nodes()=33
  n_elem()=64
    n_local_elem()=10
    n_active_elem()=64
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0

 Verbose mode disabled in non-debug mode.
 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=125
    n_local_nodes()=33
  n_elem()=64
    n_local_elem()=10
    n_active_elem()=64
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
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
    n_local_dofs()=33
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 13.3636
      Average Off-Processor Bandwidth <= 5.78788
      Maximum  On-Processor Bandwidth <= 27
      Maximum Off-Processor Bandwidth <= 17
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./miscellaneous_ex1-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:21:44 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           2.141e-01      1.01365   2.129e-01
Objects:              5.300e+01      1.03922   5.233e+01
Flops:                1.082e+05      6.21516   5.704e+04  3.423e+05
Flops/sec:            5.051e+05      6.19254   2.676e+05  1.606e+06
MPI Messages:         1.525e+02      1.54822   1.345e+02  8.070e+02
MPI Message Lengths:  1.458e+04      2.02726   8.158e+01  6.583e+04
MPI Reductions:       9.400e+01      1.02174

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.1283e-01 100.0%  3.4226e+05 100.0%  8.070e+02 100.0%  8.158e+01      100.0%  7.733e+01  84.1% 

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

VecMDot               22 1.0 1.4546e-03 7.0 1.64e+04 4.3 0.0e+00 0.0e+00 2.2e+01  0 18  0  0 24   0 18  0  0 28    42
VecNorm               24 1.0 4.1485e-03 5.7 1.58e+03 4.1 0.0e+00 0.0e+00 2.4e+01  1  2  0  0 26   1  2  0  0 31     1
VecScale              22 1.0 9.7990e-05 3.3 7.26e+02 4.1 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0    28
VecCopy                4 1.0 1.1206e-05 5.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                28 1.0 1.6952e-0428.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                2 1.0 9.3660e-03 1.4 1.32e+02 4.1 0.0e+00 0.0e+00 0.0e+00  4  0  0  0  0   4  0  0  0  0     0
VecMAXPY              23 1.0 1.3351e-05 2.2 1.82e+04 4.1 0.0e+00 0.0e+00 0.0e+00  0 20  0  0  0   0 20  0  0  0  5149
VecAssemblyBegin       3 1.0 6.4206e-04 1.2 0.00e+00 0.0 2.8e+01 9.8e+01 9.0e+00  0  0  3  4 10   0  0  3  4 12     0
VecAssemblyEnd         3 1.0 3.2902e-05 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       24 1.0 8.8453e-05 1.2 0.00e+00 0.0 6.3e+02 6.1e+01 0.0e+00  0  0 78 58  0   0  0 78 58  0     0
VecScatterEnd         24 1.0 1.7660e-03 2.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
VecNormalize          23 1.0 4.1869e-03 5.9 2.24e+03 4.1 0.0e+00 0.0e+00 2.3e+01  1  2  0  0 25   1  2  0  0 30     2
MatMult               23 1.0 1.8537e-03 2.1 2.83e+04 5.1 6.0e+02 5.9e+01 0.0e+00  1 29 74 54  0   1 29 74 54  0    53
MatSolve              23 1.0 4.1962e-05 2.3 3.33e+0412.1 0.0e+00 0.0e+00 0.0e+00  0 25  0  0  0   0 25  0  0  0  2057
MatLUFactorNum         1 1.0 2.9087e-05 1.7 9.47e+0332.9 0.0e+00 0.0e+00 0.0e+00  0  5  0  0  0   0  5  0  0  0   620
MatILUFactorSym        1 1.0 9.2983e-05 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+00  0  0  0  0  1   0  0  0  0  1     0
MatAssemblyBegin       2 1.0 1.0548e-03 2.6 0.00e+00 0.0 4.2e+01 5.4e+02 4.0e+00  0  0  5 34  4   0  0  5 34  5     0
MatAssemblyEnd         2 1.0 1.0951e-03 1.1 0.00e+00 0.0 5.2e+01 1.7e+01 8.0e+00  1  0  6  1  9   1  0  6  1 10     0
MatGetRowIJ            1 1.0 2.7180e-0514.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 8.1062e-05 2.3 0.00e+00 0.0 0.0e+00 0.0e+00 3.3e+00  0  0  0  0  4   0  0  0  0  4     0
MatZeroEntries         3 1.0 1.2875e-05 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        22 1.0 1.4918e-03 6.1 3.31e+04 4.2 0.0e+00 0.0e+00 2.2e+01  0 37  0  0 24   0 37  0  0 28    84
KSPSetup               2 1.0 7.4148e-05 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 1.3409e-02 1.0 1.08e+05 6.2 6.0e+02 5.9e+01 5.0e+01  6100 74 54 55   6100 74 54 65    26
PCSetUp                2 1.0 1.0710e-03 2.5 9.47e+0332.9 0.0e+00 0.0e+00 4.3e+00  0  5  0  0  5   0  5  0  0  6    17
PCSetUpOnBlocks        1 1.0 2.2697e-04 1.2 9.47e+0332.9 0.0e+00 0.0e+00 4.3e+00  0  5  0  0  5   0  5  0  0  6    79
PCApply               23 1.0 4.0483e-04 2.1 3.33e+0412.1 0.0e+00 0.0e+00 0.0e+00  0 25  0  0  0   0 25  0  0  0   213
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    33             33        52176     0
         Vec Scatter     2              2         1736     0
           Index Set     7              7         4372     0
   IS L to G Mapping     1              1          708     0
              Matrix     4              4        26484     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.67777e-05
Average time for zero size MPI_Send(): 5.58297e-05
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
| Time:           Fri Aug 24 15:21:44 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.357349, Active time=0.039832                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0000      0.000046    0.0001      0.000060    0.12     0.15     |
|   build_sparsity()                 1         0.0002      0.000164    0.0002      0.000200    0.41     0.50     |
|   create_dof_constraints()         1         0.0002      0.000164    0.0002      0.000164    0.41     0.41     |
|   distribute_dofs()                1         0.0001      0.000070    0.0006      0.000555    0.18     1.39     |
|   dof_indices()                    112       0.0000      0.000000    0.0000      0.000000    0.12     0.12     |
|   prepare_send_list()              1         0.0000      0.000005    0.0000      0.000005    0.01     0.01     |
|   reinit()                         1         0.0001      0.000081    0.0001      0.000081    0.20     0.20     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   write()                          1         0.0004      0.000372    0.0027      0.002692    0.93     6.76     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        10        0.0000      0.000004    0.0000      0.000004    0.11     0.11     |
|   init_shape_functions()           1         0.0000      0.000024    0.0000      0.000024    0.06     0.06     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             10        0.0000      0.000002    0.0000      0.000002    0.05     0.05     |
|   init_reference_to_physical_map() 1         0.0000      0.000027    0.0000      0.000027    0.07     0.07     |
|                                                                                                                |
| InfElemBuilder                                                                                                 |
|   build_inf_elem()                 1         0.0001      0.000052    0.0001      0.000052    0.13     0.13     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 4         0.0005      0.000122    0.0014      0.000353    1.23     3.55     |
|   renumber_nodes_and_elem()        4         0.0000      0.000007    0.0000      0.000007    0.07     0.07     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   assign_global_indices()          1         0.0014      0.001440    0.0030      0.002972    3.62     7.46     |
|   compute_hilbert_indices()        3         0.0006      0.000192    0.0006      0.000192    1.45     1.45     |
|   find_global_indices()            3         0.0001      0.000046    0.0041      0.001375    0.35     10.35    |
|   parallel_sort()                  3         0.0018      0.000591    0.0024      0.000814    4.45     6.13     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0001      0.000115    0.0001      0.000115    0.29     0.29     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      2         0.0015      0.000754    0.0039      0.001951    3.79     9.80     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      16        0.0012      0.000076    0.0012      0.000076    3.07     3.07     |
|   broadcast()                      2         0.0000      0.000006    0.0000      0.000006    0.03     0.03     |
|   gather()                         6         0.0024      0.000392    0.0024      0.000392    5.90     5.90     |
|   max(scalar)                      5         0.0009      0.000185    0.0009      0.000185    2.33     2.33     |
|   max(vector)                      5         0.0076      0.001527    0.0076      0.001527    19.17    19.17    |
|   min(vector)                      5         0.0002      0.000045    0.0002      0.000045    0.56     0.56     |
|   probe()                          90        0.0018      0.000020    0.0018      0.000020    4.43     4.43     |
|   receive()                        114       0.0001      0.000001    0.0019      0.000017    0.35     4.81     |
|   send()                           94        0.0001      0.000002    0.0001      0.000002    0.38     0.38     |
|   send_receive()                   100       0.0002      0.000002    0.0023      0.000023    0.43     5.66     |
|   sum()                            16        0.0011      0.000069    0.0011      0.000069    2.78     2.78     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           114       0.0001      0.000001    0.0001      0.000001    0.16     0.16     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         2         0.0001      0.000028    0.0007      0.000350    0.14     1.76     |
|   set_parent_processor_ids()       2         0.0000      0.000006    0.0000      0.000006    0.03     0.03     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.0165      0.016514    0.0165      0.016514    41.46    41.46    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0003      0.000293    0.0004      0.000415    0.74     1.04     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            736       0.0398                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 6 ./miscellaneous_ex1-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
