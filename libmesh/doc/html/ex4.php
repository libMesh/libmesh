<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex4",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 4 - Solving a 1D, 2D or 3D Poisson Problem in Parallel</h1>

<br><br>This is the fourth example program.  It builds on
the third example program by showing how to formulate
the code in a dimension-independent way.  Very minor
changes to the example will allow the problem to be
solved in one, two or three dimensions.

<br><br>This example will also introduce the PerfLog class
as a way to monitor your code's performance.  We will
use it to instrument the matrix assembly code and look
for bottlenecks where we should focus optimization efforts.

<br><br>This example also shows how to extend example 3 to run in
parallel.  Notice how litte has changed!  The significant
differences are marked with "PARALLEL CHANGE".


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
        #include "gmv_io.h"
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
          libMesh::init (argc, argv);
        
</pre>
</div>
<div class = "comment">
Declare a performance log for the main program
PerfLog perf_main("Main Program");
  

<br><br>Braces are used to force object scope, like in example 2
</div>

<div class ="fragment">
<pre>
          {
</pre>
</div>
<div class = "comment">
Create a GetPot object to parse the command line
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
                error();
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
Cannot use dicontinuous basis.
</div>

<div class ="fragment">
<pre>
            if ((family == "MONOMIAL") || (family == "XYZ"))
              {
        	std::cerr &lt;&lt; "ex4 currently requires a C^0 (or higher) FE basis." &lt;&lt; std::endl;
        	error();
              }
              
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
Use the MeshTools::Generation mesh generator to create a uniform
grid on the square [-1,1]^D.  We instruct the mesh generator
to build a mesh of 8x8 \p Quad9 elements in 2D, or \p Hex27
elements in 3D.  Building these higher-order elements allows
us to use higher-order approximation, as in example 3.
</div>

<div class ="fragment">
<pre>
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
        					   ps, ps, ps,
        					   -1., 1.,
        					   -1., 1.,
        					   -1., 1.,
        					   (dim==1)    ? EDGE2 : 
        					   ((dim == 2) ? QUAD4 : HEX8));
              }
            
            else
              {
        	MeshTools::Generation::build_cube (mesh,
        					   ps, ps, ps,
        					   -1., 1.,
        					   -1., 1.,
        					   -1., 1.,
        					   (dim==1)    ? EDGE3 : 
        					   ((dim == 2) ? QUAD9 : HEX27));
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
</div>

<div class ="fragment">
<pre>
            {
</pre>
</div>
<div class = "comment">
Creates a system named "Poisson"
</div>

<div class ="fragment">
<pre>
              LinearImplicitSystem& system =
        	equation_systems.add_system&lt;LinearImplicitSystem&gt; ("Poisson");
        
              
</pre>
</div>
<div class = "comment">
Adds the variable "u" to "Poisson".  "u"
will be approximated using second-order approximation.
</div>

<div class ="fragment">
<pre>
              system.add_variable("u",
        			  Utility::string_to_enum&lt;Order&gt;   (order),
        			  Utility::string_to_enum&lt;FEFamily&gt;(family));
        
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
Prints information about the system to the screen.
</div>

<div class ="fragment">
<pre>
              equation_systems.print_info();
            }
        
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
              GMVIO (mesh).write_equation_systems ((dim == 3) ? 
                "out_3.gmv" : "out_2.gmv",equation_systems);
            }
          }
          
</pre>
</div>
<div class = "comment">
All done.  
</div>

<div class ="fragment">
<pre>
          return libMesh::close ();
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
          assert (system_name == "Poisson");
        
        
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
          const Mesh& mesh = es.get_mesh();
        
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
          std::vector&lt;unsigned int&gt; dof_indices;
        
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
              perf_log.start_event("elem init");      
        
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
</div>

<div class ="fragment">
<pre>
              Ke.resize (dof_indices.size(),
        		 dof_indices.size());
        
              Fe.resize (dof_indices.size());
        
</pre>
</div>
<div class = "comment">
Stop logging the shape function initialization.
If you forget to stop logging an event the PerfLog
object will probably catch the error and abort.
</div>

<div class ="fragment">
<pre>
              perf_log.stop_event("elem init");      
        
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
              perf_log.start_event ("Ke");
        
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
              perf_log.stop_event ("Ke");
        
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
              perf_log.start_event ("Fe");
              
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
        	  const Real y = q_point[qp](1);
        	  const Real z = q_point[qp](2);
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
              perf_log.stop_event ("Fe");
        
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
                perf_log.start_event ("BCs");
        
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
        	  if (elem-&gt;neighbor(side) == NULL)
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
                        const Real yf = qface_point[qp](1);
                        const Real zf = qface_point[qp](2);
        
        
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
                perf_log.stop_event ("BCs");
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
              perf_log.start_event ("matrix insertion");
              
              system.matrix-&gt;add_matrix (Ke, dof_indices);
              system.rhs-&gt;add_vector    (Fe, dof_indices);
        
</pre>
</div>
<div class = "comment">
Start logging the insertion of the local (element)
matrix and vector into the global matrix and vector
</div>

<div class ="fragment">
<pre>
              perf_log.stop_event ("matrix insertion");
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
  #include <B><FONT COLOR="#BC8F8F">&quot;gmv_io.h&quot;</FONT></B>
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
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_poisson(EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  Real exact_solution (<B><FONT COLOR="#228B22">const</FONT></B> Real x,
  		     <B><FONT COLOR="#228B22">const</FONT></B> Real y = 0.,
  		     <B><FONT COLOR="#228B22">const</FONT></B> Real z = 0.);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    <B><FONT COLOR="#5F9EA0">libMesh</FONT></B>::init (argc, argv);
  
    
    {
      GetPot command_line (argc, argv);
      
      <B><FONT COLOR="#A020F0">if</FONT></B> (argc &lt; 3)
        {
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Usage:\n&quot;</FONT></B>
  		  &lt;&lt;<B><FONT COLOR="#BC8F8F">&quot;\t &quot;</FONT></B> &lt;&lt; argv[0] &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; -d 2(3)&quot;</FONT></B> &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; -n 15&quot;</FONT></B>
  		  &lt;&lt; std::endl;
  
  	error();
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
  	<B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ex4 currently requires a C^0 (or higher) FE basis.&quot;</FONT></B> &lt;&lt; std::endl;
  	error();
        }
        
      Mesh mesh (dim);
      
  
      <B><FONT COLOR="#A020F0">if</FONT></B> ((family == <B><FONT COLOR="#BC8F8F">&quot;LAGRANGE&quot;</FONT></B>) &amp;&amp; (order == <B><FONT COLOR="#BC8F8F">&quot;FIRST&quot;</FONT></B>))
        {
  	<B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube (mesh,
  					   ps, ps, ps,
  					   -1., 1.,
  					   -1., 1.,
  					   -1., 1.,
  					   (dim==1)    ? EDGE2 : 
  					   ((dim == 2) ? QUAD4 : HEX8));
        }
      
      <B><FONT COLOR="#A020F0">else</FONT></B>
        {
  	<B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube (mesh,
  					   ps, ps, ps,
  					   -1., 1.,
  					   -1., 1.,
  					   -1., 1.,
  					   (dim==1)    ? EDGE3 : 
  					   ((dim == 2) ? QUAD9 : HEX27));
        }
  
      
      mesh.print_info();
      
      
      EquationSystems equation_systems (mesh);
      
      {
        LinearImplicitSystem&amp; system =
  	equation_systems.add_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>);
  
        
        system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>,
  			  <B><FONT COLOR="#5F9EA0">Utility</FONT></B>::string_to_enum&lt;Order&gt;   (order),
  			  <B><FONT COLOR="#5F9EA0">Utility</FONT></B>::string_to_enum&lt;FEFamily&gt;(family));
  
        system.attach_assemble_function (assemble_poisson);
        
        equation_systems.init();
  
        equation_systems.print_info();
      }
  
      equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>).solve();
  
      <B><FONT COLOR="#A020F0">if</FONT></B>(dim == 1)
      {        
        GnuPlotIO plot(mesh,<B><FONT COLOR="#BC8F8F">&quot;Example 4, 1D&quot;</FONT></B>,GnuPlotIO::GRID_ON);
        plot.write_equation_systems(<B><FONT COLOR="#BC8F8F">&quot;out_1&quot;</FONT></B>,equation_systems);
      }
      <B><FONT COLOR="#A020F0">else</FONT></B>
      {
        GMVIO (mesh).write_equation_systems ((dim == 3) ? 
          <B><FONT COLOR="#BC8F8F">&quot;out_3.gmv&quot;</FONT></B> : <B><FONT COLOR="#BC8F8F">&quot;out_2.gmv&quot;</FONT></B>,equation_systems);
      }
    }
    
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close ();
  }
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_poisson(EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>);
  
  
    PerfLog perf_log (<B><FONT COLOR="#BC8F8F">&quot;Matrix Assembly&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> Mesh&amp; mesh = es.get_mesh();
  
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
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.local_elements_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        perf_log.start_event(<B><FONT COLOR="#BC8F8F">&quot;elem init&quot;</FONT></B>);      
  
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
  
        fe-&gt;reinit (elem);
  
        Ke.resize (dof_indices.size(),
  		 dof_indices.size());
  
        Fe.resize (dof_indices.size());
  
        perf_log.stop_event(<B><FONT COLOR="#BC8F8F">&quot;elem init&quot;</FONT></B>);      
  
        perf_log.start_event (<B><FONT COLOR="#BC8F8F">&quot;Ke&quot;</FONT></B>);
  
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
  	    Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
  	    
  
        perf_log.stop_event (<B><FONT COLOR="#BC8F8F">&quot;Ke&quot;</FONT></B>);
  
        perf_log.start_event (<B><FONT COLOR="#BC8F8F">&quot;Fe&quot;</FONT></B>);
        
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
  	{
  	  <B><FONT COLOR="#228B22">const</FONT></B> Real x = q_point[qp](0);
  	  <B><FONT COLOR="#228B22">const</FONT></B> Real y = q_point[qp](1);
  	  <B><FONT COLOR="#228B22">const</FONT></B> Real z = q_point[qp](2);
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
        
        perf_log.stop_event (<B><FONT COLOR="#BC8F8F">&quot;Fe&quot;</FONT></B>);
  
        {
  	
  	perf_log.start_event (<B><FONT COLOR="#BC8F8F">&quot;BCs&quot;</FONT></B>);
  
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> side=0; side&lt;elem-&gt;n_sides(); side++)
  	  <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(side) == NULL)
  	    {
              
  	      <B><FONT COLOR="#228B22">const</FONT></B> Real penalty = 1.e10;
  
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp;  phi_face = fe_face-&gt;get_phi();
  
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_face-&gt;get_JxW();
  
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point &gt;&amp; qface_point = fe_face-&gt;get_xyz();
  
                fe_face-&gt;reinit(elem, side);
  
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qface.n_points(); qp++)
                {
                  <B><FONT COLOR="#228B22">const</FONT></B> Real xf = qface_point[qp](0);
                  <B><FONT COLOR="#228B22">const</FONT></B> Real yf = qface_point[qp](1);
                  <B><FONT COLOR="#228B22">const</FONT></B> Real zf = qface_point[qp](2);
  
  
                  <B><FONT COLOR="#228B22">const</FONT></B> Real value = exact_solution(xf, yf, zf);
  
                  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi_face.size(); i++)
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi_face.size(); j++)
                      Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];
  
                  <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi_face.size(); i++)
                    Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
                } 
              }
              
  	
  	perf_log.stop_event (<B><FONT COLOR="#BC8F8F">&quot;BCs&quot;</FONT></B>);
        } 
        
  
        perf_log.start_event (<B><FONT COLOR="#BC8F8F">&quot;matrix insertion&quot;</FONT></B>);
        
        system.matrix-&gt;add_matrix (Ke, dof_indices);
        system.rhs-&gt;add_vector    (Fe, dof_indices);
  
        perf_log.stop_event (<B><FONT COLOR="#BC8F8F">&quot;matrix insertion&quot;</FONT></B>);
      }
  
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example  ./ex4-devel
***************************************************************
 
Running ./ex4-devel -d 1 -n 20

 Mesh Information:
  mesh_dimension()=1
  spatial_dimension()=3
  n_nodes()=41
  n_elem()=20
   n_local_elem()=20
   n_active_elem()=20
  n_subdomains()=1
  n_processors()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Poisson"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE" 
    Approximation Orders="SECOND" 
    n_dofs()=41
    n_local_dofs()=41
    n_constrained_dofs()=0
    n_vectors()=1


-------------------------------------------------------
| Time:           Wed Jun  6 12:18:07 2007             |
| OS:             Linux                                |
| HostName:       orville                              |
| OS Release:     2.6.21-1.3194.fc7PAE                 |
| OS Version:     #1 SMP Wed May 23 22:27:31 EDT 2007  |
| Machine:        i686                                 |
| Username:       peterson                             |
-------------------------------------------------------
 ------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.001533, Active time=0.000777       |
 ------------------------------------------------------------------------------
| Event                         nCalls    Total       Avg         Percent of   |
|                                         Time        Time        Active Time  |
|------------------------------------------------------------------------------|
|                                                                              |
| BCs                           20        0.0001      0.000005    13.26        |
| Fe                            20        0.0002      0.000010    26.25        |
| Ke                            20        0.0000      0.000002    5.53         |
| elem init                     20        0.0002      0.000011    28.44        |
| matrix insertion              20        0.0002      0.000010    26.51        |
 ------------------------------------------------------------------------------
| Totals:                       100       0.0008                  100.00       |
 ------------------------------------------------------------------------------

Running ./ex4-devel -d 2 -n 15

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=961
  n_elem()=225
   n_local_elem()=225
   n_active_elem()=225
  n_subdomains()=1
  n_processors()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Poisson"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE" 
    Approximation Orders="SECOND" 
    n_dofs()=961
    n_local_dofs()=961
    n_constrained_dofs()=0
    n_vectors()=1


-------------------------------------------------------
| Time:           Wed Jun  6 12:18:07 2007             |
| OS:             Linux                                |
| HostName:       orville                              |
| OS Release:     2.6.21-1.3194.fc7PAE                 |
| OS Version:     #1 SMP Wed May 23 22:27:31 EDT 2007  |
| Machine:        i686                                 |
| Username:       peterson                             |
-------------------------------------------------------
 ------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.02087, Active time=0.018112        |
 ------------------------------------------------------------------------------
| Event                         nCalls    Total       Avg         Percent of   |
|                                         Time        Time        Active Time  |
|------------------------------------------------------------------------------|
|                                                                              |
| BCs                           225       0.0042      0.000019    23.18        |
| Fe                            225       0.0046      0.000021    25.52        |
| Ke                            225       0.0038      0.000017    21.18        |
| elem init                     225       0.0037      0.000016    20.30        |
| matrix insertion              225       0.0018      0.000008    9.82         |
 ------------------------------------------------------------------------------
| Totals:                       1125      0.0181                  100.00       |
 ------------------------------------------------------------------------------

Running ./ex4-devel -d 3 -n 6

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=2197
  n_elem()=216
   n_local_elem()=216
   n_active_elem()=216
  n_subdomains()=1
  n_processors()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Poisson"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE" 
    Approximation Orders="SECOND" 
    n_dofs()=2197
    n_local_dofs()=2197
    n_constrained_dofs()=0
    n_vectors()=1


-------------------------------------------------------
| Time:           Wed Jun  6 12:18:07 2007             |
| OS:             Linux                                |
| HostName:       orville                              |
| OS Release:     2.6.21-1.3194.fc7PAE                 |
| OS Version:     #1 SMP Wed May 23 22:27:31 EDT 2007  |
| Machine:        i686                                 |
| Username:       peterson                             |
-------------------------------------------------------
 ------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.240115, Active time=0.237185       |
 ------------------------------------------------------------------------------
| Event                         nCalls    Total       Avg         Percent of   |
|                                         Time        Time        Active Time  |
|------------------------------------------------------------------------------|
|                                                                              |
| BCs                           216       0.1129      0.000523    47.59        |
| Fe                            216       0.0128      0.000059    5.41         |
| Ke                            216       0.0812      0.000376    34.22        |
| elem init                     216       0.0169      0.000078    7.12         |
| matrix insertion              216       0.0134      0.000062    5.66         |
 ------------------------------------------------------------------------------
| Totals:                       1080      0.2372                  100.00       |
 ------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./ex4-devel
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
