<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex6",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 6 - Infinite Elements for the Wave Equation</h1>

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
        #include <algorithm>
        #include <math.h>
        
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
          libMesh::init (argc, argv);
          
</pre>
</div>
<div class = "comment">
This example requires Infinite Elements   
</div>

<div class ="fragment">
<pre>
        #ifndef ENABLE_INFINITE_ELEMENTS
        
          std::cerr << "ERROR: This example requires the library to be compiled with Infinite Element support!"
        	    << std::endl;
          here();
        
          return 0;
        
        #else
          
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
Write the mesh before the infinite elements are added
</div>

<div class ="fragment">
<pre>
            mesh.write_gmv ("orig_mesh.gmv");
            
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
Write the mesh with the infinite elements added.
Compare this to the original mesh.
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
          DenseMatrix<Number> Ce;
          DenseMatrix<Number> Me;
          DenseVector<Number> Fe;
          
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
        	    
        	    es("Wave").rhs->add_vector (Fe, dof_indices);
        	  } // end boundary condition section	     
        	} // else ( if (elem->infinite())) )
        
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
              const std::vector<Real>&         weight  = cfe->get_Sobolev_weight();
              const std::vector<RealGradient>& dweight = cfe->get_Sobolev_dweight();
        
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
        	    for (unsigned int j=0; j<n_sf; j++)
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
        
              es("Wave").matrix->add_matrix (Ke, dof_indices);
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
            
            for (unsigned int n_cnt=0; n_cnt<n_nodes; n_cnt++)
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
        	    fabs(curr_node(1)) < TOLERANCE &&
        	    fabs(curr_node(2)) < TOLERANCE)
        	  {
</pre>
</div>
<div class = "comment">
The global number of the respective degree of freedom.
</div>

<div class ="fragment">
<pre>
                    unsigned int dn = curr_node.dof_number(0,0,0);
        
        	    es("Wave").rhs->add (dn, 1.);
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

<a name="nocomments"></a> 
<br><br><br> <h1> The program without comments: </h1> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Compiling C++ (in debug mode) ex6.C...
Linking ex6...
/home/benkirk/phd/code/libmesh/contrib/tecplot/lib/i686-pc-linux-gnu/tecio.a(tecxxx.o)(.text+0x1a7): In function `tecini':
: the use of `mktemp' is dangerous, better use `mkstemp'
***************************************************************
* Running Example  ./ex6
***************************************************************
 
Running ex6 with dim = 3

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=125
  n_elem()=64
   n_local_elem()=64
   n_active_elem()=64
  n_subdomains()=1
  n_processors()=1
  processor_id()=0

 Determined origin for Infinite Elements:
  0.00000 0.00000 0.00000 

 Building Infinite Elements:
  updating element neighbor tables...
  collecting boundary sides...
  found 0 inner and 96 outer boundary faces
  added 96 infinite elements and 98 nodes to the mesh

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=223
  n_elem()=160
   n_local_elem()=160
   n_active_elem()=160
  n_subdomains()=1
  n_processors()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Wave"
    Type "Steady"
    Variables="p" 
    Finite Element Types="0", "12" 
    Infinite Element Mapping="0" 
    Approximation Orders="1", "3" 
    n_dofs()=419
    n_local_dofs()=419
    n_constrained_dofs()=0
    n_additional_vectors()=0
    n_additional_matrices()=0
  n_parameters()=4
   Parameters:
    "fluid density"=1
    "linear solver maximum iterations"=5000
    "linear solver tolerance"=1e-12
    "speed"=1


 ---------------------------------------------------------------------------- 
| Reference count information                                                |
 ---------------------------------------------------------------------------- 
| 10SystemBase reference count information:
| Creations:    1
| Destructions: 1
| 12SparseMatrixISt7complexIdEE reference count information:
| Creations:    1
| Destructions: 1
| 13NumericVectorISt7complexIdEE reference count information:
| Creations:    3
| Destructions: 3
| 21LinearSolverInterfaceISt7complexIdEE reference count information:
| Creations:    1
| Destructions: 1
| 4Elem reference count information:
| Creations:    3207
| Destructions: 3207
| 4Node reference count information:
| Creations:    223
| Destructions: 223
| 5QBase reference count information:
| Creations:    5
| Destructions: 5
| 6DofMap reference count information:
| Creations:    1
| Destructions: 1
| 6FEBase reference count information:
| Creations:    3
| Destructions: 3
 ---------------------------------------------------------------------------- 

 ----------------------------------------------------------------------------
| Time:           Tue Nov 18 07:41:39 2003
| OS:             Linux
| HostName:       voyager.home.net
| OS Release      2.4.20-20.9
| OS Version:     #1 Mon Aug 18 11:45:58 EDT 2003
| Machine:        i686
| Username:       benkirk
 ----------------------------------------------------------------------------
 ----------------------------------------------------------------------------
| libMesh Performance: Alive time=0.436066, Active time=0.373535
 ----------------------------------------------------------------------------
| Event                         nCalls  Total       Avg         Percent of   |
|                                       Time        Time        Active Time  |
|----------------------------------------------------------------------------|
|                                                                            |
|                                                                            |
| DofMap                                                                     |
|   compute_sparsity()          1       0.0120      0.012025    3.22         |
|   create_dof_constraints()    1       0.0001      0.000088    0.02         |
|   distribute_dofs()           1       0.0005      0.000507    0.14         |
|   dof_indices()               320     0.0037      0.000012    0.99         |
|   reinit()                    1       0.0013      0.001250    0.33         |
|                                                                            |
| FE                                                                         |
|   compute_map()               160     0.0055      0.000034    1.47         |
|   compute_shape_functions()   64      0.0006      0.000010    0.17         |
|   init_shape_functions()      2       0.0004      0.000181    0.10         |
|                                                                            |
| InfFE                                                                      |
|   combine_base_radial()       96      0.0057      0.000060    1.53         |
|   compute_shape_functions()   96      0.0036      0.000038    0.97         |
|   init_radial_shape_functions()1       0.0001      0.000071    0.02         |
|   init_shape_functions()      1       0.0002      0.000248    0.07         |
|                                                                            |
| Mesh                                                                       |
|   build_cube()                1       0.0011      0.001093    0.29         |
|                                                                            |
| MeshBase                                                                   |
|   build_inf_elem()            1       0.0035      0.003458    0.93         |
|   find_neighbors()            4       0.0112      0.002809    3.01         |
|   renumber_nodes_and_elem()   2       0.0002      0.000085    0.05         |
|                                                                            |
| SystemBase                                                                 |
|   assemble()                  1       0.1991      0.199150    53.31        |
|   solve()                     1       0.1247      0.124697    33.38        |
 ----------------------------------------------------------------------------
| Totals:                       754     0.3735                  100.00       |
 ----------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./ex6
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
