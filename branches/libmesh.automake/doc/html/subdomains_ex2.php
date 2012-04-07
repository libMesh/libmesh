<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("subdomains_ex2",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Subdomains Example 2 - Subdomain-Restricted Variables</h1>

<br><br>This example builds on the fourth example program by showing how
to restrict solution fields to a subdomain (or union of
subdomains).


<br><br>

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
        #include "mesh_generation.h"
        #include "exodusII_io.h"
        #include "gnuplot_io.h"
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
Define the DofMap, which handles degree of freedom
indexing.
</div>

<div class ="fragment">
<pre>
        #include "dof_map.h"
        
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
Define the PerfLog, a performance logging utility.
It is useful for timing events in a code and giving
you an idea where bottlenecks lie.
</div>

<div class ="fragment">
<pre>
        #include "perf_log.h"
        
</pre>
</div>
<div class = "comment">
The definition of a geometric element
</div>

<div class ="fragment">
<pre>
        #include "elem.h"
        
        #include "string_to_enum.h"
        #include "getpot.h"
        
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
function will take the \p EquationSystems object and the
name of the system we are assembling as input.  From the
\p EquationSystems object we have acess to the \p Mesh and
other objects we might need.
</div>

<div class ="fragment">
<pre>
        void assemble_poisson(EquationSystems& es,
                              const std::string& system_name);
        
</pre>
</div>
<div class = "comment">
Exact solution function prototype.
</div>

<div class ="fragment">
<pre>
        Real exact_solution (const Real x,
                             const Real y = 0.,
                             const Real z = 0.);
        
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
Initialize libMesh and any dependent libaries, like in example 2.
</div>

<div class ="fragment">
<pre>
          LibMeshInit init (argc, argv);
        
</pre>
</div>
<div class = "comment">
Declare a performance log for the main program
PerfLog perf_main("Main Program");
  

<br><br>Create a GetPot object to parse the command line
</div>

<div class ="fragment">
<pre>
          GetPot command_line (argc, argv);
          
</pre>
</div>
<div class = "comment">
Check for proper calling arguments.
</div>

<div class ="fragment">
<pre>
          if (argc &lt; 3)
            {
              if (libMesh::processor_id() == 0)
                std::cerr &lt;&lt; "Usage:\n"
                          &lt;&lt;"\t " &lt;&lt; argv[0] &lt;&lt; " -d 2(3)" &lt;&lt; " -n 15"
                          &lt;&lt; std::endl;
        
</pre>
</div>
<div class = "comment">
This handy function will print the file name, line number,
and then abort.  Currrently the library does not use C++
exception handling.
</div>

<div class ="fragment">
<pre>
              libmesh_error();
            }
          
</pre>
</div>
<div class = "comment">
Brief message to the user regarding the program name
and command line arguments.
</div>

<div class ="fragment">
<pre>
          else 
            {
              std::cout &lt;&lt; "Running " &lt;&lt; argv[0];
              
              for (int i=1; i&lt;argc; i++)
                std::cout &lt;&lt; " " &lt;&lt; argv[i];
              
              std::cout &lt;&lt; std::endl &lt;&lt; std::endl;
            }
          
        
</pre>
</div>
<div class = "comment">
Read problem dimension from command line.  Use int
instead of unsigned since the GetPot overload is ambiguous
otherwise.
</div>

<div class ="fragment">
<pre>
          int dim = 2;
          if ( command_line.search(1, "-d") )
            dim = command_line.next(dim);
          
</pre>
</div>
<div class = "comment">
Skip higher-dimensional examples on a lower-dimensional libMesh build
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(dim &lt;= LIBMESH_DIM, "2D/3D support");
            
</pre>
</div>
<div class = "comment">
Create a mesh with user-defined dimension.
</div>

<div class ="fragment">
<pre>
          Mesh mesh (dim);
          
</pre>
</div>
<div class = "comment">
Read number of elements from command line
</div>

<div class ="fragment">
<pre>
          int ps = 15;
          if ( command_line.search(1, "-n") )
            ps = command_line.next(ps);
          
</pre>
</div>
<div class = "comment">
Read FE order from command line
</div>

<div class ="fragment">
<pre>
          std::string order = "SECOND"; 
          if ( command_line.search(2, "-Order", "-o") )
            order = command_line.next(order);
        
</pre>
</div>
<div class = "comment">
Read FE Family from command line
</div>

<div class ="fragment">
<pre>
          std::string family = "LAGRANGE"; 
          if ( command_line.search(2, "-FEFamily", "-f") )
            family = command_line.next(family);
          
</pre>
</div>
<div class = "comment">
Cannot use discontinuous basis.
</div>

<div class ="fragment">
<pre>
          if ((family == "MONOMIAL") || (family == "XYZ"))
            {
              if (libMesh::processor_id() == 0)
                std::cerr &lt;&lt; "ex28 currently requires a C^0 (or higher) FE basis." &lt;&lt; std::endl;
              libmesh_error();
            }
        
</pre>
</div>
<div class = "comment">
Use the MeshTools::Generation mesh generator to create a uniform
grid on the square [-1,1]^D.  We instruct the mesh generator
to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
elements in 3D.  Building these higher-order elements allows
us to use higher-order approximation, as in example 3.


<br><br></div>

<div class ="fragment">
<pre>
          Real halfwidth = dim &gt; 1 ? 1. : 0.;
          Real halfheight = dim &gt; 2 ? 1. : 0.;
        
          if ((family == "LAGRANGE") && (order == "FIRST"))
            {
</pre>
</div>
<div class = "comment">
No reason to use high-order geometric elements if we are
solving with low-order finite elements.
</div>

<div class ="fragment">
<pre>
              MeshTools::Generation::build_cube (mesh,
                                                 ps,
        					 (dim&gt;1) ? ps : 0,
        					 (dim&gt;2) ? ps : 0,
                                                 -1., 1.,
                                                 -halfwidth, halfwidth,
                                                 -halfheight, halfheight,
                                                 (dim==1)    ? EDGE2 : 
                                                 ((dim == 2) ? QUAD4 : HEX8));
            }
          
          else
            {
              MeshTools::Generation::build_cube (mesh,
        					 ps,
        					 (dim&gt;1) ? ps : 0,
        					 (dim&gt;2) ? ps : 0,
                                                 -1., 1.,
                                                 -halfwidth, halfwidth,
                                                 -halfheight, halfheight,
                                                 (dim==1)    ? EDGE3 : 
                                                 ((dim == 2) ? QUAD9 : HEX27));
            }
        
          {
            MeshBase::element_iterator       el     = mesh.elements_begin();
            const MeshBase::element_iterator end_el = mesh.elements_end();
            
            for ( ; el != end_el; ++el)
              {
        	Elem* elem = *el;
        	const Point cent = elem-&gt;centroid();
                if (dim &gt; 1)
                  {
        	    if ((cent(0) &gt; 0) == (cent(1) &gt; 0))
        	      elem-&gt;subdomain_id() = 1;	
                  }
                else
                  {
        	    if (cent(0) &gt; 0)
        	      elem-&gt;subdomain_id() = 1;	
                  }
              }
          }
        
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
Create an equation systems object.
</div>

<div class ="fragment">
<pre>
          EquationSystems equation_systems (mesh);
          
</pre>
</div>
<div class = "comment">
Declare the system and its variables.
Create a system named "Poisson"
</div>

<div class ="fragment">
<pre>
          LinearImplicitSystem& system =
            equation_systems.add_system&lt;LinearImplicitSystem&gt; ("Poisson");
        
          
          std::set&lt;subdomain_id_type&gt; active_subdomains;
        
          
</pre>
</div>
<div class = "comment">
Add the variable "u" to "Poisson".  "u"
will be approximated using second-order approximation.
</div>

<div class ="fragment">
<pre>
          active_subdomains.clear(); active_subdomains.insert(0);
          system.add_variable("u",
                              Utility::string_to_enum&lt;Order&gt;   (order),
                              Utility::string_to_enum&lt;FEFamily&gt;(family),
        		      &active_subdomains);
        
</pre>
</div>
<div class = "comment">
Add the variable "v" to "Poisson".  "v"
will be approximated using second-order approximation.
</div>

<div class ="fragment">
<pre>
          active_subdomains.clear(); active_subdomains.insert(1);
          system.add_variable("v",
                              Utility::string_to_enum&lt;Order&gt;   (order),
                              Utility::string_to_enum&lt;FEFamily&gt;(family),
        		      &active_subdomains);
        
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly
function.
</div>

<div class ="fragment">
<pre>
          system.attach_assemble_function (assemble_poisson);
          
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
Print information about the system to the screen.
</div>

<div class ="fragment">
<pre>
          equation_systems.print_info();
          mesh.print_info();
        
</pre>
</div>
<div class = "comment">
Solve the system "Poisson", just like example 2.
</div>

<div class ="fragment">
<pre>
          equation_systems.get_system("Poisson").solve();
        
</pre>
</div>
<div class = "comment">
After solving the system write the solution
to a GMV-formatted plot file.
</div>

<div class ="fragment">
<pre>
          if(dim == 1)
          {        
            GnuPlotIO plot(mesh,"Example 4, 1D",GnuPlotIO::GRID_ON);
            plot.write_equation_systems("out_1",equation_systems);
          }
          else
          {
        #ifdef LIBMESH_HAVE_EXODUS_API
            ExodusII_IO (mesh).write_equation_systems ((dim == 3) ? 
              "out_3.e" : "out_2.e",equation_systems);
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
          }
          
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

<br><br>
<br><br>
<br><br>We now define the matrix assembly function for the
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
Declare a performance log.  Give it a descriptive
string to identify what part of the code we are
logging, since there may be many PerfLogs in an
application.
</div>

<div class ="fragment">
<pre>
          PerfLog perf_log ("Matrix Assembly");
          
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
          LinearImplicitSystem& system = es.get_system&lt;LinearImplicitSystem&gt;("Poisson");
          
</pre>
</div>
<div class = "comment">
A reference to the \p DofMap object for this system.  The \p DofMap
object handles the index translation from node and element numbers
to degree of freedom numbers.  We will talk more about the \p DofMap
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
\p FEBase::build() member dynamically creates memory we will
store the object as an \p AutoPtr<FEBase>.  This can be thought
of as a pointer that will clean up after itself.
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
Tell the finte element object to use our
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
We begin with the element Jacobian * quadrature weight at each
integration point.   
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
"Ke" and "Fe". More detail is in example 3.
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
          std::vector&lt;unsigned int&gt; dof_indices, dof_indices2;
        
</pre>
</div>
<div class = "comment">
Now we will loop over all the elements in the mesh.
We will compute the element matrix and right-hand-side
contribution.  See example 3 for a discussion of the
element iterators.  Here we use the \p const_local_elem_iterator
to indicate we only want to loop over elements that are assigned
to the local processor.  This allows each processor to compute
its components of the global matrix.

<br><br>"PARALLEL CHANGE"
</div>

<div class ="fragment">
<pre>
          MeshBase::const_element_iterator       el     = mesh.local_elements_begin();
          const MeshBase::const_element_iterator end_el = mesh.local_elements_end();
        
          for ( ; el != end_el; ++el)
            {
</pre>
</div>
<div class = "comment">
Start logging the shape function initialization.
This is done through a simple function call with
the name of the event to log.
</div>

<div class ="fragment">
<pre>
              perf_log.push("elem init");      
        
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
              dof_map.dof_indices (elem, dof_indices,0);
              dof_map.dof_indices (elem, dof_indices2,1);
        
</pre>
</div>
<div class = "comment">
std::cout << "dof_indices.size()="
<< dof_indices.size() 
<< ", dof_indices2.size()="
<< dof_indices2.size()
<< std::endl;


<br><br>Compute the element-specific data for the current
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
</div>

<div class ="fragment">
<pre>
              Ke.resize (std::max(dof_indices.size(), dof_indices2.size()),
        		 std::max(dof_indices.size(), dof_indices2.size()));
        
              Fe.resize (std::max(dof_indices.size(), dof_indices2.size()));
        
</pre>
</div>
<div class = "comment">
Stop logging the shape function initialization.
If you forget to stop logging an event the PerfLog
object will probably catch the error and abort.
</div>

<div class ="fragment">
<pre>
              perf_log.pop("elem init");      
        
</pre>
</div>
<div class = "comment">
Now we will build the element matrix.  This involves
a double loop to integrate the test funcions (i) against
the trial functions (j).

<br><br>We have split the numeric integration into two loops
so that we can log the matrix and right-hand-side
computation seperately.

<br><br>Now start logging the element matrix computation
</div>

<div class ="fragment">
<pre>
              perf_log.push ("Ke");
        
              for (unsigned int qp=0; qp&lt;qrule.n_points(); qp++)
                for (unsigned int i=0; i&lt;phi.size(); i++)
                  for (unsigned int j=0; j&lt;phi.size(); j++)
                    Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
                    
        
</pre>
</div>
<div class = "comment">
Stop logging the matrix computation
</div>

<div class ="fragment">
<pre>
              perf_log.pop ("Ke");
        
</pre>
</div>
<div class = "comment">
Now we build the element right-hand-side contribution.
This involves a single loop in which we integrate the
"forcing function" in the PDE against the test functions.

<br><br>Start logging the right-hand-side computation
</div>

<div class ="fragment">
<pre>
              perf_log.push ("Fe");
              
              for (unsigned int qp=0; qp&lt;qrule.n_points(); qp++)
                {
</pre>
</div>
<div class = "comment">
fxy is the forcing function for the Poisson equation.
In this case we set fxy to be a finite difference
Laplacian approximation to the (known) exact solution.

<br><br>We will use the second-order accurate FD Laplacian
approximation, which in 2D on a structured grid is

<br><br>u_xx + u_yy = (u(i-1,j) + u(i+1,j) +
u(i,j-1) + u(i,j+1) +
-4*u(i,j))/h^2

<br><br>Since the value of the forcing function depends only
on the location of the quadrature point (q_point[qp])
we will compute it here, outside of the i-loop          
</div>

<div class ="fragment">
<pre>
                  const Real x = q_point[qp](0);
        #if LIBMESH_DIM &gt; 1
                  const Real y = q_point[qp](1);
        #else
                  const Real y = 0;
        #endif
        #if LIBMESH_DIM &gt; 2
                  const Real z = q_point[qp](2);
        #else
                  const Real z = 0;
        #endif
                  const Real eps = 1.e-3;
        
                  const Real uxx = (exact_solution(x-eps,y,z) +
                                    exact_solution(x+eps,y,z) +
                                    -2.*exact_solution(x,y,z))/eps/eps;
                      
                  const Real uyy = (exact_solution(x,y-eps,z) +
                                    exact_solution(x,y+eps,z) +
                                    -2.*exact_solution(x,y,z))/eps/eps;
                  
                  const Real uzz = (exact_solution(x,y,z-eps) +
                                    exact_solution(x,y,z+eps) +
                                    -2.*exact_solution(x,y,z))/eps/eps;
        
                  Real fxy;
                  if(dim==1)
                  {
</pre>
</div>
<div class = "comment">
In 1D, compute the rhs by differentiating the
exact solution twice.
</div>

<div class ="fragment">
<pre>
                    const Real pi = libMesh::pi;
                    fxy = (0.25*pi*pi)*sin(.5*pi*x);
                  }
                  else
                  {
                    fxy = - (uxx + uyy + ((dim==2) ? 0. : uzz));
                  } 
        
</pre>
</div>
<div class = "comment">
Add the RHS contribution
</div>

<div class ="fragment">
<pre>
                  for (unsigned int i=0; i&lt;phi.size(); i++)
                    Fe(i) += JxW[qp]*fxy*phi[i][qp];          
                }
              
</pre>
</div>
<div class = "comment">
Stop logging the right-hand-side computation
</div>

<div class ="fragment">
<pre>
              perf_log.pop ("Fe");
        
</pre>
</div>
<div class = "comment">
At this point the interior element integration has
been completed.  However, we have not yet addressed
boundary conditions.  For this example we will only
consider simple Dirichlet boundary conditions imposed
via the penalty method. This is discussed at length in
example 3.
</div>

<div class ="fragment">
<pre>
              {
                
</pre>
</div>
<div class = "comment">
Start logging the boundary condition computation
</div>

<div class ="fragment">
<pre>
                perf_log.push ("BCs");
        
</pre>
</div>
<div class = "comment">
The following loops over the sides of the element.
If the element has no neighbor on a side then that
side MUST live on a boundary of the domain.
</div>

<div class ="fragment">
<pre>
                for (unsigned int side=0; side&lt;elem-&gt;n_sides(); side++)
                  if ((elem-&gt;neighbor(side) == NULL) ||
        	      (elem-&gt;neighbor(side)-&gt;subdomain_id() != elem-&gt;subdomain_id()))
                    {
                    
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
        #if LIBMESH_DIM &gt; 1
                        const Real yf = qface_point[qp](1);
        #else
                        const Real yf = 0.;
        #endif
        #if LIBMESH_DIM &gt; 2
                        const Real zf = qface_point[qp](2);
        #else
                        const Real zf = 0.;
        #endif
        
        
</pre>
</div>
<div class = "comment">
The boundary value.
</div>

<div class ="fragment">
<pre>
                        const Real value = exact_solution(xf, yf, zf);
        
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
                    
                
</pre>
</div>
<div class = "comment">
Stop logging the boundary condition computation
</div>

<div class ="fragment">
<pre>
                perf_log.pop ("BCs");
              } 
              
        
</pre>
</div>
<div class = "comment">
The element matrix and right-hand-side are now built
for this element.  Add them to the global matrix and
right-hand-side vector.  The \p PetscMatrix::add_matrix()
and \p PetscVector::add_vector() members do this for us.
Start logging the insertion of the local (element)
matrix and vector into the global matrix and vector
</div>

<div class ="fragment">
<pre>
              perf_log.push ("matrix insertion");
              
              if (dof_indices.size())
        	{
        	  system.matrix-&gt;add_matrix (Ke, dof_indices);
        	  system.rhs-&gt;add_vector    (Fe, dof_indices);
        	}
              
              if (dof_indices2.size())
        	{
        	  system.matrix-&gt;add_matrix (Ke, dof_indices2);
        	  system.rhs-&gt;add_vector    (Fe, dof_indices2);
        	}
        
</pre>
</div>
<div class = "comment">
Start logging the insertion of the local (element)
matrix and vector into the global matrix and vector
</div>

<div class ="fragment">
<pre>
              perf_log.pop ("matrix insertion");
            }
        
</pre>
</div>
<div class = "comment">
That's it.  We don't need to do anything else to the
PerfLog.  When it goes out of scope (at this function return)
it will print its log to the screen. Pretty easy, huh?
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
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;gnuplot_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;perf_log.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;string_to_enum.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  
  using namespace libMesh;
   
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_poisson(EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  Real exact_solution (<B><FONT COLOR="#228B22">const</FONT></B> Real x,
                       <B><FONT COLOR="#228B22">const</FONT></B> Real y = 0.,
                       <B><FONT COLOR="#228B22">const</FONT></B> Real z = 0.);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
    
    GetPot command_line (argc, argv);
    
    <B><FONT COLOR="#A020F0">if</FONT></B> (argc &lt; 3)
      {
        <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Usage:\n&quot;</FONT></B>
                    &lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;\t &quot;</FONT></B> &lt;&lt; argv[0] &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; -d 2(3)&quot;</FONT></B> &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; -n 15&quot;</FONT></B>
                    &lt;&lt; std::endl;
  
        libmesh_error();
      }
    
    <B><FONT COLOR="#A020F0">else</FONT></B> 
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Running &quot;</FONT></B> &lt;&lt; argv[0];
        
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">int</FONT></B> i=1; i&lt;argc; i++)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B> &lt;&lt; argv[i];
        
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl &lt;&lt; std::endl;
      }
    
  
    <B><FONT COLOR="#228B22">int</FONT></B> dim = 2;
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <B><FONT COLOR="#BC8F8F">&quot;-d&quot;</FONT></B>) )
      dim = command_line.next(dim);
    
    libmesh_example_assert(dim &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D/3D support&quot;</FONT></B>);
      
    Mesh mesh (dim);
    
    <B><FONT COLOR="#228B22">int</FONT></B> ps = 15;
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(1, <B><FONT COLOR="#BC8F8F">&quot;-n&quot;</FONT></B>) )
      ps = command_line.next(ps);
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string order = <B><FONT COLOR="#BC8F8F">&quot;SECOND&quot;</FONT></B>; 
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(2, <B><FONT COLOR="#BC8F8F">&quot;-Order&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;-o&quot;</FONT></B>) )
      order = command_line.next(order);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string family = <B><FONT COLOR="#BC8F8F">&quot;LAGRANGE&quot;</FONT></B>; 
    <B><FONT COLOR="#A020F0">if</FONT></B> ( command_line.search(2, <B><FONT COLOR="#BC8F8F">&quot;-FEFamily&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;-f&quot;</FONT></B>) )
      family = command_line.next(family);
    
    <B><FONT COLOR="#A020F0">if</FONT></B> ((family == <B><FONT COLOR="#BC8F8F">&quot;MONOMIAL&quot;</FONT></B>) || (family == <B><FONT COLOR="#BC8F8F">&quot;XYZ&quot;</FONT></B>))
      {
        <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ex28 currently requires a C^0 (or higher) FE basis.&quot;</FONT></B> &lt;&lt; std::endl;
        libmesh_error();
      }
  
  
    Real halfwidth = dim &gt; 1 ? 1. : 0.;
    Real halfheight = dim &gt; 2 ? 1. : 0.;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> ((family == <B><FONT COLOR="#BC8F8F">&quot;LAGRANGE&quot;</FONT></B>) &amp;&amp; (order == <B><FONT COLOR="#BC8F8F">&quot;FIRST&quot;</FONT></B>))
      {
        <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube (mesh,
                                           ps,
  					 (dim&gt;1) ? ps : 0,
  					 (dim&gt;2) ? ps : 0,
                                           -1., 1.,
                                           -halfwidth, halfwidth,
                                           -halfheight, halfheight,
                                           (dim==1)    ? EDGE2 : 
                                           ((dim == 2) ? QUAD4 : HEX8));
      }
    
    <B><FONT COLOR="#A020F0">else</FONT></B>
      {
        <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube (mesh,
  					 ps,
  					 (dim&gt;1) ? ps : 0,
  					 (dim&gt;2) ? ps : 0,
                                           -1., 1.,
                                           -halfwidth, halfwidth,
                                           -halfheight, halfheight,
                                           (dim==1)    ? EDGE3 : 
                                           ((dim == 2) ? QUAD9 : HEX27));
      }
  
    {
      <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::element_iterator       el     = mesh.elements_begin();
      <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::element_iterator end_el = mesh.elements_end();
      
      <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
        {
  	Elem* elem = *el;
  	<B><FONT COLOR="#228B22">const</FONT></B> Point cent = elem-&gt;centroid();
          <B><FONT COLOR="#A020F0">if</FONT></B> (dim &gt; 1)
            {
  	    <B><FONT COLOR="#A020F0">if</FONT></B> ((cent(0) &gt; 0) == (cent(1) &gt; 0))
  	      elem-&gt;subdomain_id() = 1;	
            }
          <B><FONT COLOR="#A020F0">else</FONT></B>
            {
  	    <B><FONT COLOR="#A020F0">if</FONT></B> (cent(0) &gt; 0)
  	      elem-&gt;subdomain_id() = 1;	
            }
        }
    }
  
    mesh.print_info();
      
    EquationSystems equation_systems (mesh);
    
    LinearImplicitSystem&amp; system =
      equation_systems.add_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>);
  
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;subdomain_id_type&gt; active_subdomains;
  
    
    active_subdomains.clear(); active_subdomains.insert(0);
    system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>,
                        <B><FONT COLOR="#5F9EA0">Utility</FONT></B>::string_to_enum&lt;Order&gt;   (order),
                        <B><FONT COLOR="#5F9EA0">Utility</FONT></B>::string_to_enum&lt;FEFamily&gt;(family),
  		      &amp;active_subdomains);
  
    active_subdomains.clear(); active_subdomains.insert(1);
    system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>,
                        <B><FONT COLOR="#5F9EA0">Utility</FONT></B>::string_to_enum&lt;Order&gt;   (order),
                        <B><FONT COLOR="#5F9EA0">Utility</FONT></B>::string_to_enum&lt;FEFamily&gt;(family),
  		      &amp;active_subdomains);
  
    system.attach_assemble_function (assemble_poisson);
    
    equation_systems.init();
  
    equation_systems.print_info();
    mesh.print_info();
  
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>).solve();
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(dim == 1)
    {        
      GnuPlotIO plot(mesh,<B><FONT COLOR="#BC8F8F">&quot;Example 4, 1D&quot;</FONT></B>,GnuPlotIO::GRID_ON);
      plot.write_equation_systems(<B><FONT COLOR="#BC8F8F">&quot;out_1&quot;</FONT></B>,equation_systems);
    }
    <B><FONT COLOR="#A020F0">else</FONT></B>
    {
  #ifdef LIBMESH_HAVE_EXODUS_API
      ExodusII_IO (mesh).write_equation_systems ((dim == 3) ? 
        <B><FONT COLOR="#BC8F8F">&quot;out_3.e&quot;</FONT></B> : <B><FONT COLOR="#BC8F8F">&quot;out_2.e&quot;</FONT></B>,equation_systems);
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>    }
    
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_poisson(EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>);
  
    PerfLog perf_log (<B><FONT COLOR="#BC8F8F">&quot;Matrix Assembly&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    LinearImplicitSystem&amp; system = es.get_system&lt;LinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>);
    
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
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices, dof_indices2;
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.local_elements_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        perf_log.push(<B><FONT COLOR="#BC8F8F">&quot;elem init&quot;</FONT></B>);      
  
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices,0);
        dof_map.dof_indices (elem, dof_indices2,1);
  
  
        fe-&gt;reinit (elem);
  
        Ke.resize (std::max(dof_indices.size(), dof_indices2.size()),
  		 <B><FONT COLOR="#5F9EA0">std</FONT></B>::max(dof_indices.size(), dof_indices2.size()));
  
        Fe.resize (std::max(dof_indices.size(), dof_indices2.size()));
  
        perf_log.pop(<B><FONT COLOR="#BC8F8F">&quot;elem init&quot;</FONT></B>);      
  
        perf_log.push (<B><FONT COLOR="#BC8F8F">&quot;Ke&quot;</FONT></B>);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
              Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
              
  
        perf_log.pop (<B><FONT COLOR="#BC8F8F">&quot;Ke&quot;</FONT></B>);
  
        perf_log.push (<B><FONT COLOR="#BC8F8F">&quot;Fe&quot;</FONT></B>);
        
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
          {
            <B><FONT COLOR="#228B22">const</FONT></B> Real x = q_point[qp](0);
  #<B><FONT COLOR="#A020F0">if</FONT></B> LIBMESH_DIM &gt; 1
            <B><FONT COLOR="#228B22">const</FONT></B> Real y = q_point[qp](1);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
            <B><FONT COLOR="#228B22">const</FONT></B> Real y = 0;
  #endif
  #<B><FONT COLOR="#A020F0">if</FONT></B> LIBMESH_DIM &gt; 2
            <B><FONT COLOR="#228B22">const</FONT></B> Real z = q_point[qp](2);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
            <B><FONT COLOR="#228B22">const</FONT></B> Real z = 0;
  #endif
            <B><FONT COLOR="#228B22">const</FONT></B> Real eps = 1.e-3;
  
            <B><FONT COLOR="#228B22">const</FONT></B> Real uxx = (exact_solution(x-eps,y,z) +
                              exact_solution(x+eps,y,z) +
                              -2.*exact_solution(x,y,z))/eps/eps;
                
            <B><FONT COLOR="#228B22">const</FONT></B> Real uyy = (exact_solution(x,y-eps,z) +
                              exact_solution(x,y+eps,z) +
                              -2.*exact_solution(x,y,z))/eps/eps;
            
            <B><FONT COLOR="#228B22">const</FONT></B> Real uzz = (exact_solution(x,y,z-eps) +
                              exact_solution(x,y,z+eps) +
                              -2.*exact_solution(x,y,z))/eps/eps;
  
            Real fxy;
            <B><FONT COLOR="#A020F0">if</FONT></B>(dim==1)
            {
              <B><FONT COLOR="#228B22">const</FONT></B> Real pi = libMesh::pi;
              fxy = (0.25*pi*pi)*sin(.5*pi*x);
            }
            <B><FONT COLOR="#A020F0">else</FONT></B>
            {
              fxy = - (uxx + uyy + ((dim==2) ? 0. : uzz));
            } 
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
              Fe(i) += JxW[qp]*fxy*phi[i][qp];          
          }
        
        perf_log.pop (<B><FONT COLOR="#BC8F8F">&quot;Fe&quot;</FONT></B>);
  
        {
          
          perf_log.push (<B><FONT COLOR="#BC8F8F">&quot;BCs&quot;</FONT></B>);
  
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> side=0; side&lt;elem-&gt;n_sides(); side++)
            <B><FONT COLOR="#A020F0">if</FONT></B> ((elem-&gt;neighbor(side) == NULL) ||
  	      (elem-&gt;neighbor(side)-&gt;subdomain_id() != elem-&gt;subdomain_id()))
              {
              
                <B><FONT COLOR="#228B22">const</FONT></B> Real penalty = 1.e10;
  
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp;  phi_face = fe_face-&gt;get_phi();
  
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
  
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point &gt;&amp; qface_point = fe_face-&gt;get_xyz();
  
                fe_face-&gt;reinit(elem, side);
  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qface.n_points(); qp++)
                {
                  <B><FONT COLOR="#228B22">const</FONT></B> Real xf = qface_point[qp](0);
  #<B><FONT COLOR="#A020F0">if</FONT></B> LIBMESH_DIM &gt; 1
                  <B><FONT COLOR="#228B22">const</FONT></B> Real yf = qface_point[qp](1);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
                  <B><FONT COLOR="#228B22">const</FONT></B> Real yf = 0.;
  #endif
  #<B><FONT COLOR="#A020F0">if</FONT></B> LIBMESH_DIM &gt; 2
                  <B><FONT COLOR="#228B22">const</FONT></B> Real zf = qface_point[qp](2);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
                  <B><FONT COLOR="#228B22">const</FONT></B> Real zf = 0.;
  #endif
  
  
                  <B><FONT COLOR="#228B22">const</FONT></B> Real value = exact_solution(xf, yf, zf);
  
                  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi_face.size(); i++)
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi_face.size(); j++)
                      Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];
  
                  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi_face.size(); i++)
                    Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
                } 
              }
              
          
          perf_log.pop (<B><FONT COLOR="#BC8F8F">&quot;BCs&quot;</FONT></B>);
        } 
        
  
        perf_log.push (<B><FONT COLOR="#BC8F8F">&quot;matrix insertion&quot;</FONT></B>);
        
        <B><FONT COLOR="#A020F0">if</FONT></B> (dof_indices.size())
  	{
  	  system.matrix-&gt;add_matrix (Ke, dof_indices);
  	  system.rhs-&gt;add_vector    (Fe, dof_indices);
  	}
        
        <B><FONT COLOR="#A020F0">if</FONT></B> (dof_indices2.size())
  	{
  	  system.matrix-&gt;add_matrix (Ke, dof_indices2);
  	  system.rhs-&gt;add_vector    (Fe, dof_indices2);
  	}
  
        perf_log.pop (<B><FONT COLOR="#BC8F8F">&quot;matrix insertion&quot;</FONT></B>);
      }
  
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Compiling C++ (in optimized mode) subdomains_ex2.C...
Linking subdomains_ex2-opt...
***************************************************************
* Running Example  ./subdomains_ex2-opt
***************************************************************
 
Running ./subdomains_ex2-opt -d 1 -n 20

 Mesh Information:
  mesh_dimension()=1
  spatial_dimension()=3
  n_nodes()=41
    n_local_nodes()=41
  n_elem()=20
    n_local_elem()=20
    n_active_elem()=20
  n_subdomains()=2
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Poisson"
    Type "LinearImplicit"
    Variables="u" "v" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" "SECOND", "THIRD" 
    n_dofs()=42
    n_local_dofs()=42
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 3.85714
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 5
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=1
  spatial_dimension()=3
  n_nodes()=41
    n_local_nodes()=41
  n_elem()=20
    n_local_elem()=20
    n_active_elem()=20
  n_subdomains()=2
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0


-------------------------------------------------------------------
| Time:           Sat Apr  7 16:03:19 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.001311, Active time=0.000867                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           20        0.0002      0.000009    0.0002      0.000009    20.30    20.30    |
| Fe                            20        0.0002      0.000010    0.0002      0.000010    24.22    24.22    |
| Ke                            20        0.0000      0.000001    0.0000      0.000001    1.85     1.85     |
| elem init                     20        0.0004      0.000020    0.0004      0.000020    46.48    46.48    |
| matrix insertion              20        0.0001      0.000003    0.0001      0.000003    7.15     7.15     |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       100       0.0009                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.078748, Active time=0.005497                                             |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 1         0.0000      0.000034    0.0000      0.000034    0.62     0.62     |
|   compute_sparsity()           1         0.0002      0.000249    0.0003      0.000282    4.53     5.13     |
|   create_dof_constraints()     1         0.0000      0.000001    0.0000      0.000001    0.02     0.02     |
|   distribute_dofs()            1         0.0002      0.000161    0.0003      0.000347    2.93     6.31     |
|   dof_indices()                80        0.0001      0.000001    0.0001      0.000001    1.69     1.69     |
|   prepare_send_list()          1         0.0000      0.000002    0.0000      0.000002    0.04     0.04     |
|   reinit()                     1         0.0002      0.000184    0.0002      0.000184    3.35     3.35     |
|                                                                                                            |
| EquationSystems                                                                                            |
|   build_solution_vector()      1         0.0002      0.000152    0.0002      0.000185    2.77     3.37     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         24        0.0000      0.000002    0.0000      0.000002    0.71     0.71     |
|   compute_face_map()           4         0.0000      0.000002    0.0000      0.000002    0.15     0.15     |
|   compute_shape_functions()    24        0.0000      0.000001    0.0000      0.000001    0.29     0.29     |
|   init_face_shape_functions()  1         0.0000      0.000008    0.0000      0.000008    0.15     0.15     |
|   init_shape_functions()       5         0.0001      0.000020    0.0001      0.000020    1.82     1.82     |
|                                                                                                            |
| GnuPlotIO                                                                                                  |
|   write_nodal_data()           1         0.0007      0.000696    0.0007      0.000696    12.66    12.66    |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             1         0.0002      0.000202    0.0002      0.000202    3.67     3.67     |
|   renumber_nodes_and_elem()    2         0.0000      0.000005    0.0000      0.000005    0.16     0.16     |
|                                                                                                            |
| MeshOutput                                                                                                 |
|   write_equation_systems()     1         0.0000      0.000040    0.0009      0.000922    0.73     16.77    |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0001      0.000141    0.0001      0.000141    2.57     2.57     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  1         0.0000      0.000001    0.0000      0.000001    0.02     0.02     |
|                                                                                                            |
| Partitioner                                                                                                |
|   single_partition()           1         0.0000      0.000024    0.0000      0.000024    0.44     0.44     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      1         0.0020      0.002015    0.0020      0.002015    36.66    36.66    |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   1         0.0013      0.001322    0.0016      0.001563    24.05    28.43    |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        155       0.0055                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

Running ./subdomains_ex2-opt -d 2 -n 15

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=961
    n_local_nodes()=961
  n_elem()=225
    n_local_elem()=225
    n_active_elem()=225
  n_subdomains()=2
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Poisson"
    Type "LinearImplicit"
    Variables="u" "v" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" "SECOND", "THIRD" 
    n_dofs()=1022
    n_local_dofs()=1022
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 14.5616
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 25
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=961
    n_local_nodes()=961
  n_elem()=225
    n_local_elem()=225
    n_active_elem()=225
  n_subdomains()=2
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0


-------------------------------------------------------------------
| Time:           Sat Apr  7 16:03:19 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.021012, Active time=0.019377                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           225       0.0073      0.000032    0.0073      0.000032    37.68    37.68    |
| Fe                            225       0.0040      0.000018    0.0040      0.000018    20.57    20.57    |
| Ke                            225       0.0016      0.000007    0.0016      0.000007    8.24     8.24     |
| elem init                     225       0.0046      0.000021    0.0046      0.000021    23.97    23.97    |
| matrix insertion              225       0.0018      0.000008    0.0018      0.000008    9.54     9.54     |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       1125      0.0194                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.107728, Active time=0.040074                                             |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 1         0.0002      0.000236    0.0002      0.000236    0.59     0.59     |
|   compute_sparsity()           1         0.0026      0.002559    0.0030      0.002991    6.39     7.46     |
|   create_dof_constraints()     1         0.0002      0.000184    0.0002      0.000184    0.46     0.46     |
|   distribute_dofs()            1         0.0008      0.000772    0.0021      0.002117    1.93     5.28     |
|   dof_indices()                900       0.0014      0.000002    0.0014      0.000002    3.41     3.41     |
|   prepare_send_list()          1         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   reinit()                     1         0.0013      0.001343    0.0013      0.001343    3.35     3.35     |
|                                                                                                            |
| EquationSystems                                                                                            |
|   build_solution_vector()      1         0.0005      0.000459    0.0008      0.000779    1.15     1.94     |
|                                                                                                            |
| ExodusII_IO                                                                                                |
|   write_nodal_data()           1         0.0016      0.001648    0.0016      0.001648    4.11     4.11     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         345       0.0010      0.000003    0.0010      0.000003    2.54     2.54     |
|   compute_face_map()           120       0.0013      0.000011    0.0032      0.000027    3.36     7.95     |
|   compute_shape_functions()    345       0.0005      0.000002    0.0005      0.000002    1.36     1.36     |
|   init_face_shape_functions()  1         0.0000      0.000018    0.0000      0.000018    0.04     0.04     |
|   init_shape_functions()       121       0.0018      0.000014    0.0018      0.000014    4.37     4.37     |
|   inverse_map()                360       0.0017      0.000005    0.0017      0.000005    4.30     4.30     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             1         0.0010      0.000958    0.0010      0.000958    2.39     2.39     |
|   renumber_nodes_and_elem()    2         0.0001      0.000066    0.0001      0.000066    0.33     0.33     |
|                                                                                                            |
| MeshOutput                                                                                                 |
|   write_equation_systems()     1         0.0000      0.000039    0.0025      0.002467    0.10     6.16     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0010      0.001015    0.0010      0.001015    2.53     2.53     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                            |
| Partitioner                                                                                                |
|   single_partition()           1         0.0001      0.000107    0.0001      0.000107    0.27     0.27     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      1         0.0092      0.009198    0.0092      0.009198    22.95    22.95    |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   1         0.0137      0.013651    0.0213      0.021300    34.06    53.15    |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        2209      0.0401                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

Running ./subdomains_ex2-opt -d 3 -n 6

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=2197
    n_local_nodes()=2197
  n_elem()=216
    n_local_elem()=216
    n_active_elem()=216
  n_subdomains()=2
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Poisson"
    Type "LinearImplicit"
    Variables="u" "v" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" "SECOND", "THIRD" 
    n_dofs()=2522
    n_local_dofs()=2522
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 48.5337
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 125
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=2197
    n_local_nodes()=2197
  n_elem()=216
    n_local_elem()=216
    n_active_elem()=216
  n_subdomains()=2
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0


-------------------------------------------------------------------
| Time:           Sat Apr  7 16:03:19 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.105987, Active time=0.104937                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           216       0.0599      0.000277    0.0599      0.000277    57.04    57.04    |
| Fe                            216       0.0073      0.000034    0.0073      0.000034    6.98     6.98     |
| Ke                            216       0.0202      0.000093    0.0202      0.000093    19.24    19.24    |
| elem init                     216       0.0082      0.000038    0.0082      0.000038    7.81     7.81     |
| matrix insertion              216       0.0094      0.000043    0.0094      0.000043    8.93     8.93     |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       1080      0.1049                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.242877, Active time=0.161441                                             |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 1         0.0005      0.000463    0.0005      0.000463    0.29     0.29     |
|   compute_sparsity()           1         0.0100      0.010022    0.0115      0.011502    6.21     7.12     |
|   create_dof_constraints()     1         0.0002      0.000177    0.0002      0.000177    0.11     0.11     |
|   distribute_dofs()            1         0.0012      0.001190    0.0034      0.003390    0.74     2.10     |
|   dof_indices()                864       0.0025      0.000003    0.0025      0.000003    1.56     1.56     |
|   prepare_send_list()          1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                     1         0.0022      0.002198    0.0022      0.002198    1.36     1.36     |
|                                                                                                            |
| EquationSystems                                                                                            |
|   build_solution_vector()      1         0.0008      0.000845    0.0012      0.001180    0.52     0.73     |
|                                                                                                            |
| ExodusII_IO                                                                                                |
|   write_nodal_data()           1         0.0030      0.003007    0.0030      0.003007    1.86     1.86     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         576       0.0041      0.000007    0.0041      0.000007    2.51     2.51     |
|   compute_face_map()           360       0.0015      0.000004    0.0015      0.000004    0.91     0.91     |
|   compute_shape_functions()    576       0.0024      0.000004    0.0024      0.000004    1.49     1.49     |
|   init_face_shape_functions()  1         0.0001      0.000070    0.0001      0.000070    0.04     0.04     |
|   init_shape_functions()       361       0.0401      0.000111    0.0401      0.000111    24.83    24.83    |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             1         0.0011      0.001122    0.0011      0.001122    0.69     0.69     |
|   renumber_nodes_and_elem()    2         0.0003      0.000150    0.0003      0.000150    0.19     0.19     |
|                                                                                                            |
| MeshOutput                                                                                                 |
|   write_equation_systems()     1         0.0000      0.000038    0.0042      0.004226    0.02     2.62     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0018      0.001791    0.0018      0.001791    1.11     1.11     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                            |
| Partitioner                                                                                                |
|   single_partition()           1         0.0002      0.000153    0.0002      0.000153    0.09     0.09     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      1         0.0327      0.032663    0.0327      0.032663    20.23    20.23    |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   1         0.0569      0.056867    0.1062      0.106167    35.22    65.76    |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        2755      0.1614                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./subdomains_ex2-opt
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
