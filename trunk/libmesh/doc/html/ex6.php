<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("examples",$root)?>
 
<div class="content">
<div class = "comment">
Example 6 -- Infinite Elements

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
        #include "libmesh.h"
        #include "mesh.h"
        #include "steady_system.h"
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
Function prototype.  This is similar to the Poisson
assemble function of example 4.  
</div>

<div class ="fragment">
<pre>
        void assemble_wave(EquationSystems&amp; es,
                           const std::string&amp; system_name);
        
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
          libMesh::init (argc, argv);
          
</pre>
</div>
<div class = "comment">
This example requires Infinite Elements   
</div>

<div class ="fragment">
<pre>
        #ifndef ENABLE_INFINITE_ELEMENTS
        
          std::cerr &lt;&lt; "ERROR: This example requires the library to be compiled with Infinite Element support!"
                    &lt;&lt; std::endl;
          here();
        
          return 0;
        
        #else
          
</pre>
</div>
<div class = "comment">
This example is designed for real numbers only.   
</div>

<div class ="fragment">
<pre>
        #ifdef USE_COMPLEX_NUMBERS
          
          std::cerr &lt;&lt; "ERROR: Not intended for use with complex numbers."
                    &lt;&lt; std::endl;
          here();
          
          return 0;
          
        #endif
          
</pre>
</div>
<div class = "comment">
Braces are used to force object scope, like in example 2      
</div>

<div class ="fragment">
<pre>
          {        
</pre>
</div>
<div class = "comment">
For the moment, only allow 3D     
</div>

<div class ="fragment">
<pre>
            const unsigned int dim = 3; 
            
</pre>
</div>
<div class = "comment">
Tell the user what we are doing.
</div>

<div class ="fragment">
<pre>
            std::cout &lt;&lt; "Running ex6 with dim = " &lt;&lt; dim &lt;&lt; std::endl &lt;&lt; std::endl;        
            
</pre>
</div>
<div class = "comment">
Create a mesh with user-defined dimension 
</div>

<div class ="fragment">
<pre>
            Mesh mesh (dim);
        
</pre>
</div>
<div class = "comment">
Use the internal mesh generator to create elements
on the square [-1,1]^3, of type Hex8.
</div>

<div class ="fragment">
<pre>
            mesh.build_cube (4, 4, 4,
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
            mesh.build_inf_elem(true);
        
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
Save only the mesh, with infinite elements added.
</div>

<div class ="fragment">
<pre>
            mesh.write_gmv ("ifems_added.gmv");
            
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
</div>

<div class ="fragment">
<pre>
            {
</pre>
</div>
<div class = "comment">
Create a system named "Wave".  This can
be a simple, steady system
</div>

<div class ="fragment">
<pre>
              equation_systems.add_system&lt;SteadySystem&gt; ("Wave");
                    
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
              equation_systems("Wave").add_variable("p", fe_type);
              
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly
function.
</div>

<div class ="fragment">
<pre>
              equation_systems("Wave").attach_assemble_function (assemble_wave);
              
</pre>
</div>
<div class = "comment">
Set the speed of sound and fluid density
as \p EquationSystems parameter,
so that \p assemble_wave() can access it.
</div>

<div class ="fragment">
<pre>
              equation_systems.set_parameter("speed")          = 1.;
              equation_systems.set_parameter("fluid density")  = 1.;
              
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
            }
            
</pre>
</div>
<div class = "comment">
Solve the system "Wave".
</div>

<div class ="fragment">
<pre>
            equation_systems("Wave").solve();
            
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
          }
          
</pre>
</div>
<div class = "comment">
All done.  
</div>

<div class ="fragment">
<pre>
          return libMesh::close ();
        
        #endif // else part of ifndef ENABLE_INFINITE_ELEMENTS
        }
        
</pre>
</div>
<div class = "comment">
This function assembles the system matrix and right-hand-side
for the discrete form of our wave equation.
</div>

<div class ="fragment">
<pre>
        void assemble_wave(EquationSystems&amp; es,
                           const std::string&amp; system_name)
        {
</pre>
</div>
<div class = "comment">
It is a good idea to make sure we are assembling
the proper system.
</div>

<div class ="fragment">
<pre>
          assert (system_name == "Wave");
        
        
        #ifdef ENABLE_INFINITE_ELEMENTS
          
</pre>
</div>
<div class = "comment">
Get a constant reference to the mesh object.
</div>

<div class ="fragment">
<pre>
          const Mesh&amp; mesh = es.get_mesh();
        
</pre>
</div>
<div class = "comment">
A reference to the \p DofMap object for this system.  The \p DofMap
object handles the index translation from node and element numbers
to degree of freedom numbers.
</div>

<div class ="fragment">
<pre>
          const DofMap&amp; dof_map = es("Wave").get_dof_map();
          
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
          const Real speed = es.parameter("speed");
          
</pre>
</div>
<div class = "comment">
Get a constant reference to the Finite Element type
for the first (and only) variable in the system.
</div>

<div class ="fragment">
<pre>
          const FEType&amp; fe_type = dof_map.variable_type(0);
          
</pre>
</div>
<div class = "comment">
Build a Finite Element object of the specified type.  Since the
\p FEBase::build() member dynamically creates memory we will
store the object as an \p AutoPtr&lt;FEBase&gt;.  Check ex5 for details.
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
          fe-&gt;attach_quadrature_rule (&amp;qrule);
          
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
          inf_fe-&gt;attach_quadrature_rule (&amp;qrule);
          
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
          const_local_elem_iterator           el (mesh.elements_begin());
          const const_local_elem_iterator end_el (mesh.elements_end());
          
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
                    
                    es("Wave").rhs-&gt;add_vector (Fe, dof_indices);
                  } // end boundary condition section             
                } // else ( if (elem-&gt;infinite())) )
        
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
This is slightly different from the Poisson solver:
Since the finite element object may change, we have to
initialize the constant references to the data fields
each time again, when a new element is processed.

<br><br>The element Jacobian * quadrature weight at each integration point.          
</div>

<div class ="fragment">
<pre>
              const std::vector&lt;Real&gt;&amp; JxW = cfe-&gt;get_JxW();
              
</pre>
</div>
<div class = "comment">
The element shape functions evaluated at the quadrature points.       
</div>

<div class ="fragment">
<pre>
              const std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = cfe-&gt;get_phi();
              
</pre>
</div>
<div class = "comment">
The element shape function gradients evaluated at the quadrature
points.       
</div>

<div class ="fragment">
<pre>
              const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = cfe-&gt;get_dphi();
        
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
              const std::vector&lt;RealGradient&gt;&amp; dphase  = cfe-&gt;get_dphase();
              const std::vector&lt;Real&gt;&amp;         weight  = cfe-&gt;get_Sobolev_weight();
              const std::vector&lt;RealGradient&gt;&amp; dweight = cfe-&gt;get_Sobolev_dweight();
        
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
        
              es("Wave").matrix-&gt;add_matrix (Ke, dof_indices);
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
Number of nodes in the mesh.     
</div>

<div class ="fragment">
<pre>
            const unsigned int n_nodes = mesh.n_nodes();
            
            for (unsigned int n_cnt=0; n_cnt&lt;n_nodes; n_cnt++)
              {        
</pre>
</div>
<div class = "comment">
Get a reference to the current node.
</div>

<div class ="fragment">
<pre>
                const Node&amp; curr_node = mesh.node(n_cnt);
                
</pre>
</div>
<div class = "comment">
Check the location of the current node.
</div>

<div class ="fragment">
<pre>
                if (fabs(curr_node(0)) &lt; TOLERANCE &amp;&amp;
                    fabs(curr_node(1)) &lt; TOLERANCE &amp;&amp;
                    fabs(curr_node(2)) &lt; TOLERANCE)
                  {
</pre>
</div>
<div class = "comment">
The global number of the respective degree of freedom.
</div>

<div class ="fragment">
<pre>
                    unsigned int dn = curr_node.dof_number(0,0,0);
        
                    es("Wave").rhs-&gt;add (dn, 1.);
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
          assert(es.get_mesh().mesh_dimension() != 1);
        
        #endif //ifdef ENABLE_INFINITE_ELEMENTS
          
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

<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;libmesh.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;mesh.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;steady_system.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;equation_systems.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;fe.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;inf_fe.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;quadrature_gauss.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;sparse_matrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;numeric_vector.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_matrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_vector.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;dof_map.h&quot;</FONT></B>
  
  <FONT COLOR="#228B22"><B>void</FONT></B> assemble_wave(EquationSystems&amp; es,
  		   <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name);
  
  <FONT COLOR="#228B22"><B>int</FONT></B> main (<FONT COLOR="#228B22"><B>int</FONT></B> argc, <FONT COLOR="#228B22"><B>char</FONT></B>** argv)
  {
    libMesh::init (argc, argv);
    
  #ifndef ENABLE_INFINITE_ELEMENTS
  
    std::cerr &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;ERROR: This example requires the library to be compiled with Infinite Element support!&quot;</FONT></B>
  	    &lt;&lt; std::endl;
    here();
  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  
  #<B><FONT COLOR="#A020F0">else</FONT></B>
    
  #ifdef USE_COMPLEX_NUMBERS
    
    std::cerr &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;ERROR: Not intended for use with complex numbers.&quot;</FONT></B>
  	    &lt;&lt; std::endl;
    here();
    
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
    
  #endif
    
    {        
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = 3; 
      
      std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Running ex6 with dim = &quot;</FONT></B> &lt;&lt; dim &lt;&lt; std::endl &lt;&lt; std::endl;        
      
      Mesh mesh (dim);
  
      mesh.build_cube (4, 4, 4,
  		     -1., 1.,
  		     -1., 1.,
  		     -1., 1.,
  		     HEX8);
      
      mesh.print_info();
      
      mesh.build_inf_elem(true);
  
      mesh.print_info();
  
      mesh.write_gmv (<FONT COLOR="#BC8F8F"><B>&quot;ifems_added.gmv&quot;</FONT></B>);
      
      mesh.find_neighbors();
      
      EquationSystems equation_systems (mesh);
      
      {
        equation_systems.add_system&lt;SteadySystem&gt; (<FONT COLOR="#BC8F8F"><B>&quot;Wave&quot;</FONT></B>);
              
        FEType fe_type(FIRST);
        
        equation_systems(<FONT COLOR="#BC8F8F"><B>&quot;Wave&quot;</FONT></B>).add_variable(<FONT COLOR="#BC8F8F"><B>&quot;p&quot;</FONT></B>, fe_type);
        
        equation_systems(<FONT COLOR="#BC8F8F"><B>&quot;Wave&quot;</FONT></B>).attach_assemble_function (assemble_wave);
        
        equation_systems.set_parameter(<FONT COLOR="#BC8F8F"><B>&quot;speed&quot;</FONT></B>)          = 1.;
        equation_systems.set_parameter(<FONT COLOR="#BC8F8F"><B>&quot;fluid density&quot;</FONT></B>)  = 1.;
        
        equation_systems.init();
        
        equation_systems.print_info();
      }
      
      equation_systems(<FONT COLOR="#BC8F8F"><B>&quot;Wave&quot;</FONT></B>).solve();
      
      equation_systems.write (<FONT COLOR="#BC8F8F"><B>&quot;eqn_sys.dat&quot;</FONT></B>, libMeshEnums::WRITE);
    }
    
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close ();
  
  #endif <I><FONT COLOR="#B22222">// else part of ifndef ENABLE_INFINITE_ELEMENTS
</FONT></I>  }
  
  <FONT COLOR="#228B22"><B>void</FONT></B> assemble_wave(EquationSystems&amp; es,
  		   <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name)
  {
    assert (system_name == <FONT COLOR="#BC8F8F"><B>&quot;Wave&quot;</FONT></B>);
  
  
  #ifdef ENABLE_INFINITE_ELEMENTS
    
    <FONT COLOR="#228B22"><B>const</FONT></B> Mesh&amp; mesh = es.get_mesh();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> DofMap&amp; dof_map = es(<FONT COLOR="#BC8F8F"><B>&quot;Wave&quot;</FONT></B>).get_dof_map();
    
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = mesh.mesh_dimension();
    
    <FONT COLOR="#228B22"><B>const</FONT></B> Real speed = es.parameter(<FONT COLOR="#BC8F8F"><B>&quot;speed&quot;</FONT></B>);
    
    <FONT COLOR="#228B22"><B>const</FONT></B> FEType&amp; fe_type = dof_map.variable_type(0);
    
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    AutoPtr&lt;FEBase&gt; inf_fe (FEBase::build_InfFE(dim, fe_type));
    
    QGauss qrule (dim, SECOND);
    
    fe-&gt;attach_quadrature_rule (&amp;qrule);
    
    inf_fe-&gt;attach_quadrature_rule (&amp;qrule);
    
    DenseMatrix&lt;Number&gt; Ke;
    DenseMatrix&lt;Number&gt; Ce;
    DenseMatrix&lt;Number&gt; Me;
    DenseVector&lt;Number&gt; Fe;
    
    std::vector&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt; dof_indices;
    
    const_local_elem_iterator           el (mesh.elements_begin());
    <FONT COLOR="#228B22"><B>const</FONT></B> const_local_elem_iterator end_el (mesh.elements_end());
    
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {      
        <FONT COLOR="#228B22"><B>const</FONT></B> Elem* elem = *el;
        
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
  	    
  	    es(<FONT COLOR="#BC8F8F"><B>&quot;Wave&quot;</FONT></B>).rhs-&gt;add_vector (Fe, dof_indices);
  	  } <I><FONT COLOR="#B22222">// end boundary condition section	     
</FONT></I>  	} <I><FONT COLOR="#B22222">// else ( if (elem-&gt;infinite())) )
</FONT></I>  
        cfe-&gt;reinit (elem);
        
        <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = cfe-&gt;get_JxW();
        
        <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = cfe-&gt;get_phi();
        
        <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = cfe-&gt;get_dphi();
  
        <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;RealGradient&gt;&amp; dphase  = cfe-&gt;get_dphase();
        <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt;&amp;         weight  = cfe-&gt;get_Sobolev_weight();
        <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;RealGradient&gt;&amp; dweight = cfe-&gt;get_Sobolev_dweight();
  
        Ke.resize (dof_indices.size(), dof_indices.size());
        Ce.resize (dof_indices.size(), dof_indices.size());
        Me.resize (dof_indices.size(), dof_indices.size());
        
        <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> max_qp = cfe-&gt;n_quadrature_points();
        
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp&lt;max_qp; qp++)
          {	  
  	  <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_sf = cfe-&gt;n_shape_functions();
  
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;n_sf; i++)
  	    <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j&lt;n_sf; j++)
  	      {
  		Ke(i,j) +=
  		  (                            <I><FONT COLOR="#B22222">//    (                         
</FONT></I>  		   (                           <I><FONT COLOR="#B22222">//      (                       
</FONT></I>  		    dweight[qp] * phi[i][qp]   <I><FONT COLOR="#B22222">//        Point * Real  = Point 
</FONT></I>  		    +                          <I><FONT COLOR="#B22222">//        +                     
</FONT></I>  		    dphi[i][qp] * weight[qp]   <I><FONT COLOR="#B22222">//        Point * Real  = Point 
</FONT></I>  		    ) * dphi[j][qp]            <I><FONT COLOR="#B22222">//      )       * Point = Real  
</FONT></I>  		   ) * JxW[qp];                <I><FONT COLOR="#B22222">//    )         * Real  = Real  
</FONT></I>  
  		Ce(i,j) +=
  		  (                                <I><FONT COLOR="#B22222">//    (                         
</FONT></I>  		   (dphase[qp] * dphi[j][qp])      <I><FONT COLOR="#B22222">//      (Point * Point) = Real  
</FONT></I>  		   * weight[qp] * phi[i][qp]       <I><FONT COLOR="#B22222">//      * Real * Real   = Real  
</FONT></I>  		   -                               <I><FONT COLOR="#B22222">//      -                       
</FONT></I>  		   (dweight[qp] * dphase[qp])      <I><FONT COLOR="#B22222">//      (Point * Point) = Real  
</FONT></I>  		   * phi[i][qp] * phi[j][qp]       <I><FONT COLOR="#B22222">//      * Real * Real   = Real  
</FONT></I>  		   -                               <I><FONT COLOR="#B22222">//      -                       
</FONT></I>  		   (dphi[i][qp] * dphase[qp])      <I><FONT COLOR="#B22222">//      (Point * Point) = Real  
</FONT></I>  		   * weight[qp] * phi[j][qp]       <I><FONT COLOR="#B22222">//      * Real * Real   = Real  
</FONT></I>  		   ) * JxW[qp];                    <I><FONT COLOR="#B22222">//    )         * Real  = Real  
</FONT></I>  		
  		Me(i,j) +=
  		  (                                       <I><FONT COLOR="#B22222">//    (                                  
</FONT></I>  		   (1. - (dphase[qp] * dphase[qp]))       <I><FONT COLOR="#B22222">//      (Real  - (Point * Point)) = Real 
</FONT></I>  		   * phi[i][qp] * phi[j][qp] * weight[qp] <I><FONT COLOR="#B22222">//      * Real *  Real  * Real    = Real 
</FONT></I>  		   ) * JxW[qp];                           <I><FONT COLOR="#B22222">//    ) * Real                    = Real 
</FONT></I>  
  	      } <I><FONT COLOR="#B22222">// end of the matrix summation loop
</FONT></I>  	} <I><FONT COLOR="#B22222">// end of quadrature point loop
</FONT></I>  
        Ke.add(1./speed        , Ce);
        Ke.add(1./(speed*speed), Me);
  
        es(<FONT COLOR="#BC8F8F"><B>&quot;Wave&quot;</FONT></B>).matrix-&gt;add_matrix (Ke, dof_indices);
      } <I><FONT COLOR="#B22222">// end of element loop
</FONT></I>  
    {
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_nodes = mesh.n_nodes();
      
      <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n_cnt=0; n_cnt&lt;n_nodes; n_cnt++)
        {	
  	<FONT COLOR="#228B22"><B>const</FONT></B> Node&amp; curr_node = mesh.node(n_cnt);
  	
  	<B><FONT COLOR="#A020F0">if</FONT></B> (fabs(curr_node(0)) &lt; TOLERANCE &amp;&amp;
  	    fabs(curr_node(1)) &lt; TOLERANCE &amp;&amp;
  	    fabs(curr_node(2)) &lt; TOLERANCE)
  	  {
  	    <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dn = curr_node.dof_number(0,0,0);
  
  	    es(<FONT COLOR="#BC8F8F"><B>&quot;Wave&quot;</FONT></B>).rhs-&gt;add (dn, 1.);
  	  }
        }
    }
  
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    assert(es.get_mesh().mesh_dimension() != 1);
  
  #endif <I><FONT COLOR="#B22222">//ifdef ENABLE_INFINITE_ELEMENTS
</FONT></I>    
    <B><FONT COLOR="#A020F0">return</FONT></B>;
  }
  
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
