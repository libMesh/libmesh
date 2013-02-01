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
<br><br><br> <h1> The source file exact_solution.C with comments: </h1> 
<div class = "comment">
  

<br><br>This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.
  

<br><br>This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
  

<br><br>You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


<br><br>

<br><br>

<br><br>C++ Includes
</div>

<div class ="fragment">
<pre>
        #include &lt;math.h&gt;
        
</pre>
</div>
<div class = "comment">
Mesh library includes
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh_common.h"
        
</pre>
</div>
<div class = "comment">
Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        
        
        
        
        /**
         * This is the exact solution that
         * we are trying to obtain.  We will solve
         *
         * - (u_xx + u_yy) = f
         *
         * and take a finite difference approximation using this
         * function to get f.  This is the well-known "method of
         * manufactured solutions".
         */
        Real exact_solution (const Real x,
        		     const Real y,
        		     const Real z = 0.)
        {
          static const Real pi = acos(-1.);
        
          return cos(.5*pi*x)*sin(.5*pi*y)*cos(.5*pi*z);
        }
</pre>
</div>

<a name="comments"></a> 
<br><br><br> <h1> The source file subdomains_ex2.C with comments: </h1> 
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
        #include "libmesh/libmesh.h"
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/exodusII_io.h"
        #include "libmesh/gnuplot_io.h"
        #include "libmesh/linear_implicit_system.h"
        #include "libmesh/equation_systems.h"
        
</pre>
</div>
<div class = "comment">
Define the Finite Element object.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/fe.h"
        
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
Define the DofMap, which handles degree of freedom
indexing.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/dof_map.h"
        
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
Define the PerfLog, a performance logging utility.
It is useful for timing events in a code and giving
you an idea where bottlenecks lie.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/perf_log.h"
        
</pre>
</div>
<div class = "comment">
The definition of a geometric element
</div>

<div class ="fragment">
<pre>
        #include "libmesh/elem.h"
        
        #include "libmesh/string_to_enum.h"
        #include "libmesh/getpot.h"
        
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
            GnuPlotIO plot(mesh,"Subdomains Example 2, 1D",GnuPlotIO::GRID_ON);
            plot.write_equation_systems("gnuplot_script",equation_systems);
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
          libmesh_assert_equal_to (system_name, "Poisson");
        
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
          std::vector&lt;dof_id_type&gt; dof_indices, dof_indices2;
        
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
<br><br><br> <h1> The source file exact_solution.C without comments: </h1> 
<pre> 
    
    
    
  
  
  
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh_common.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  
  
  
  <I><FONT COLOR="#B22222">/**
   * This is the exact solution that
   * we are trying to obtain.  We will solve
   *
   * - (u_xx + u_yy) = f
   *
   * and take a finite difference approximation using this
   * function to get f.  This is the well-known &quot;method of
   * manufactured solutions&quot;.
   */</FONT></I>
  Real exact_solution (<B><FONT COLOR="#228B22">const</FONT></B> Real x,
  		     <B><FONT COLOR="#228B22">const</FONT></B> Real y,
  		     <B><FONT COLOR="#228B22">const</FONT></B> Real z = 0.)
  {
    <B><FONT COLOR="#228B22">static</FONT></B> <B><FONT COLOR="#228B22">const</FONT></B> Real pi = acos(-1.);
  
    <B><FONT COLOR="#A020F0">return</FONT></B> cos(.5*pi*x)*sin(.5*pi*y)*cos(.5*pi*z);
  }
</pre> 
<a name="nocomments"></a> 
<br><br><br> <h1> The source file subdomains_ex2.C without comments: </h1> 
<pre> 
  
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/gnuplot_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/perf_log.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/string_to_enum.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  
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
      GnuPlotIO plot(mesh,<B><FONT COLOR="#BC8F8F">&quot;Subdomains Example 2, 1D&quot;</FONT></B>,GnuPlotIO::GRID_ON);
      plot.write_equation_systems(<B><FONT COLOR="#BC8F8F">&quot;gnuplot_script&quot;</FONT></B>,equation_systems);
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
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>);
  
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
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices, dof_indices2;
  
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
***************************************************************
* Running Example subdomains_ex2:
*  mpirun -np 12 example-devel -d 1 -n 20 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/subdomains/subdomains_ex2/.libs/lt-example-devel -d 1 -n 20 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=1
  spatial_dimension()=3
  n_nodes()=41
    n_local_nodes()=5
  n_elem()=20
    n_local_elem()=2
    n_active_elem()=20
  n_subdomains()=2
  n_partitions()=12
  n_processors()=12
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
    n_local_dofs()=5
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 3.28571
      Average Off-Processor Bandwidth <= 0.761905
      Maximum  On-Processor Bandwidth <= 5
      Maximum Off-Processor Bandwidth <= 2
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=1
  spatial_dimension()=3
  n_nodes()=41
    n_local_nodes()=5
  n_elem()=20
    n_local_elem()=2
    n_active_elem()=20
  n_subdomains()=2
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 22:12:35 2013                                                                             |
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
 -----------------------------------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.001514, Active time=0.000726                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           2         0.0000      0.000001    0.0000      0.000001    0.28     0.28     |
| Fe                            2         0.0001      0.000044    0.0001      0.000044    12.12    12.12    |
| Ke                            2         0.0000      0.000005    0.0000      0.000005    1.52     1.52     |
| elem init                     2         0.0005      0.000275    0.0005      0.000275    75.62    75.62    |
| matrix insertion              2         0.0001      0.000038    0.0001      0.000038    10.47    10.47    |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       10        0.0007                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/subdomains/subdomains_ex2/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:12:35 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           5.235e-02      1.00082   5.233e-02
Objects:              5.500e+01      1.03774   5.350e+01
Flops:                6.910e+03      0.00000   3.256e+03  3.907e+04
Flops/sec:            1.320e+05      0.00000   6.220e+04  7.465e+05
MPI Messages:         5.100e+01      0.00000   3.400e+01  4.080e+02
MPI Message Lengths:  5.840e+02      0.00000   1.145e+01  4.672e+03
MPI Reductions:       8.300e+01      1.02469

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 5.2291e-02  99.9%  3.9067e+04 100.0%  4.080e+02 100.0%  1.145e+01      100.0%  8.050e+01  97.0% 

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

VecMDot               17 1.0 2.3842e-04 1.8 1.99e+03 0.0 0.0e+00 0.0e+00 1.7e+01  0 29  0  0 20   0 29  0  0 21    48
VecNorm               19 1.0 1.3924e-04 1.1 2.66e+02 0.0 0.0e+00 0.0e+00 1.9e+01  0  4  0  0 23   0  4  0  0 24    11
VecScale              18 1.0 3.0041e-05 1.9 1.26e+02 0.0 0.0e+00 0.0e+00 0.0e+00  0  2  0  0  0   0  2  0  0  0    25
VecCopy                2 1.0 3.0994e-06 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                25 1.0 1.3113e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                2 1.0 1.0701e-02 1.0 2.80e+01 0.0 0.0e+00 0.0e+00 0.0e+00 20  0  0  0  0  20  0  0  0  0     0
VecMAXPY              18 1.0 8.8215e-06 3.1 2.38e+03 0.0 0.0e+00 0.0e+00 0.0e+00  0 37  0  0  0   0 37  0  0  0  1619
VecAssemblyBegin       3 1.0 1.6594e-04 1.1 0.00e+00 0.0 1.6e+01 6.0e+00 9.0e+00  0  0  4  2 11   0  0  4  2 11     0
VecAssemblyEnd         3 1.0 3.1233e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       19 1.0 9.9659e-0513.9 0.00e+00 0.0 3.0e+02 1.2e+01 0.0e+00  0  0 75 81  0   0  0 75 81  0     0
VecScatterEnd         19 1.0 4.5061e-0515.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize          18 1.0 2.2507e-04 1.0 3.78e+02 0.0 0.0e+00 0.0e+00 1.8e+01  0  6  0  0 22   0  6  0  0 22    10
MatMult               18 1.0 1.7381e-04 3.7 9.18e+02 0.0 2.9e+02 1.2e+01 0.0e+00  0 13 71 74  0   0 13 71 74  0    29
MatSolve              19 0.0 2.0742e-05 0.0 1.12e+03 0.0 0.0e+00 0.0e+00 0.0e+00  0 14  0  0  0   0 14  0  0  0   258
MatLUFactorNum         1 1.0 5.8174e-05 3.3 8.20e+01 0.0 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0     6
MatILUFactorSym        1 1.0 9.8944e-05 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  0  0  0  0  4   0  0  0  0  4     0
MatAssemblyBegin       2 1.0 1.5759e-03 8.9 0.00e+00 0.0 2.4e+01 1.7e+01 4.0e+00  1  0  6  9  5   1  0  6  9  5     0
MatAssemblyEnd         2 1.0 4.5300e-04 1.0 0.00e+00 0.0 3.2e+01 5.0e+00 8.0e+00  1  0  8  3 10   1  0  8  3 10     0
MatGetRowIJ            1 0.0 1.0014e-05 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 0.0 1.2398e-04 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries         3 0.0 1.8835e-05 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        17 1.0 2.8563e-04 1.6 4.13e+03 0.0 0.0e+00 0.0e+00 1.7e+01  0 62  0  0 20   0 62  0  0 21    85
KSPSetUp               2 1.0 1.2589e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 1.2938e-02 1.0 6.91e+03 0.0 2.9e+02 1.2e+01 4.4e+01 25100 71 74 52  25100 71 74 54     3
PCSetUp                2 1.0 8.6808e-04 1.1 8.20e+01 0.0 0.0e+00 0.0e+00 7.5e+00  2  1  0  0  9   2  1  0  0  9     0
PCSetUpOnBlocks        1 1.0 3.7599e-04 1.1 8.20e+01 0.0 0.0e+00 0.0e+00 5.5e+00  1  1  0  0  7   1  1  0  0  7     1
PCApply               19 1.0 2.9063e-04 1.2 1.12e+03 0.0 0.0e+00 0.0e+00 0.0e+00  1 14  0  0  0   1 14  0  0  0    18
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    34             34        51832     0
      Vector Scatter     2              2         2072     0
           Index Set     9              9         6768     0
   IS L to G Mapping     1              1          564     0
              Matrix     4              4        10836     0
       Krylov Solver     2              2        19360     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 1.19209e-07
Average time for MPI_Barrier(): 5.19753e-06
Average time for zero size MPI_Send(): 1.33316e-05
#PETSc Option Table entries:
-d 1
-ksp_right_pc
-log_summary
-n 20
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

 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.200062, Active time=0.040676                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0003      0.000309    0.0006      0.000571    0.76     1.40     |
|   build_sparsity()                 1         0.0009      0.000935    0.0022      0.002179    2.30     5.36     |
|   create_dof_constraints()         1         0.0000      0.000030    0.0000      0.000030    0.07     0.07     |
|   distribute_dofs()                1         0.0014      0.001385    0.0040      0.003960    3.40     9.74     |
|   dof_indices()                    12        0.0007      0.000057    0.0007      0.000057    1.68     1.68     |
|   prepare_send_list()              1         0.0000      0.000017    0.0000      0.000017    0.04     0.04     |
|   reinit()                         1         0.0013      0.001291    0.0013      0.001291    3.17     3.17     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0003      0.000344    0.0012      0.001247    0.85     3.07     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        2         0.0000      0.000025    0.0000      0.000025    0.12     0.12     |
|   init_shape_functions()           1         0.0002      0.000169    0.0002      0.000169    0.42     0.42     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             2         0.0000      0.000020    0.0000      0.000020    0.10     0.10     |
|   init_reference_to_physical_map() 1         0.0001      0.000078    0.0001      0.000078    0.19     0.19     |
|                                                                                                                |
| GnuPlotIO                                                                                                      |
|   write_nodal_data()               1         0.0012      0.001240    0.0012      0.001240    3.05     3.05     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0011      0.001135    0.0013      0.001276    2.79     3.14     |
|   renumber_nodes_and_elem()        2         0.0001      0.000062    0.0001      0.000062    0.30     0.30     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0005      0.000254    0.0005      0.000254    1.25     1.25     |
|   find_global_indices()            2         0.0008      0.000386    0.0052      0.002586    1.90     12.72    |
|   parallel_sort()                  2         0.0023      0.001146    0.0026      0.001323    5.63     6.51     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0002      0.000166    0.0028      0.002806    0.41     6.90     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0006      0.000588    0.0006      0.000588    1.45     1.45     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0029      0.002874    0.0047      0.004676    7.07     11.50    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      11        0.0004      0.000033    0.0004      0.000038    0.90     1.02     |
|   max(bool)                        1         0.0000      0.000007    0.0000      0.000007    0.02     0.02     |
|   max(scalar)                      115       0.0010      0.000008    0.0010      0.000008    2.35     2.35     |
|   max(vector)                      26        0.0004      0.000015    0.0010      0.000040    0.96     2.54     |
|   min(bool)                        133       0.0011      0.000008    0.0011      0.000008    2.66     2.66     |
|   min(scalar)                      109       0.0021      0.000019    0.0021      0.000019    5.05     5.05     |
|   min(vector)                      26        0.0004      0.000017    0.0012      0.000048    1.10     3.04     |
|   probe()                          132       0.0006      0.000005    0.0006      0.000005    1.57     1.57     |
|   receive()                        132       0.0007      0.000006    0.0014      0.000011    1.83     3.48     |
|   send()                           132       0.0004      0.000003    0.0004      0.000003    0.99     0.99     |
|   send_receive()                   136       0.0011      0.000008    0.0033      0.000025    2.79     8.22     |
|   sum()                            22        0.0004      0.000020    0.0010      0.000045    1.07     2.44     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           132       0.0003      0.000002    0.0003      0.000002    0.72     0.72     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0005      0.000533    0.0011      0.001056    1.31     2.60     |
|   set_parent_processor_ids()       1         0.0001      0.000137    0.0001      0.000137    0.34     0.34     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.0146      0.014607    0.0146      0.014607    35.91    35.91    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0014      0.001409    0.0019      0.001894    3.46     4.66     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1149      0.0407                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example subdomains_ex2:
*  mpirun -np 12 example-devel -d 1 -n 20 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
***************************************************************
* Running Example subdomains_ex2:
*  mpirun -np 12 example-devel -d 2 -n 15 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/subdomains/subdomains_ex2/.libs/lt-example-devel -d 2 -n 15 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=961
    n_local_nodes()=97
  n_elem()=225
    n_local_elem()=19
    n_active_elem()=225
  n_subdomains()=2
  n_partitions()=12
  n_processors()=12
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
    n_local_dofs()=97
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 12.5294
      Average Off-Processor Bandwidth <= 2.63209
      Maximum  On-Processor Bandwidth <= 26
      Maximum Off-Processor Bandwidth <= 18
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=961
    n_local_nodes()=97
  n_elem()=225
    n_local_elem()=19
    n_active_elem()=225
  n_subdomains()=2
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 22:12:36 2013                                                                             |
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
 -----------------------------------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.010369, Active time=0.009646                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           19        0.0023      0.000123    0.0023      0.000123    24.25    24.25    |
| Fe                            19        0.0003      0.000016    0.0003      0.000016    3.11     3.11     |
| Ke                            19        0.0012      0.000062    0.0012      0.000062    12.13    12.13    |
| elem init                     19        0.0056      0.000297    0.0056      0.000297    58.55    58.55    |
| matrix insertion              19        0.0002      0.000010    0.0002      0.000010    1.96     1.96     |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       95        0.0096                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/subdomains/subdomains_ex2/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:12:36 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           1.785e-01      1.00034   1.785e-01
Objects:              6.400e+01      1.03226   6.217e+01
Flops:                1.051e+06      1.75249   7.208e+05  8.649e+06
Flops/sec:            5.891e+06      1.75190   4.039e+06  4.847e+07
MPI Messages:         4.140e+02      3.47899   2.572e+02  3.087e+03
MPI Message Lengths:  3.796e+04      2.03561   1.046e+02  3.228e+05
MPI Reductions:       1.500e+02      1.01351

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.7842e-01 100.0%  8.6494e+06 100.0%  3.087e+03 100.0%  1.046e+02      100.0%  1.472e+02  99.4% 

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

VecMDot               50 1.0 5.4178e-03 2.2 1.30e+05 1.4 0.0e+00 0.0e+00 5.0e+01  2 16  0  0 34   2 16  0  0 34   253
VecNorm               53 1.0 1.0948e-02 4.3 1.03e+04 1.4 0.0e+00 0.0e+00 5.3e+01  4  1  0  0 36   4  1  0  0 36    10
VecScale              52 1.0 4.3869e-05 1.6 5.04e+03 1.4 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  1211
VecCopy                3 1.0 8.8215e-06 4.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                60 1.0 3.8862e-05 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                4 1.0 1.5324e-02 2.3 7.76e+02 1.4 0.0e+00 0.0e+00 0.0e+00  5  0  0  0  0   5  0  0  0  0     1
VecMAXPY              52 1.0 9.9421e-05 2.2 1.41e+05 1.4 0.0e+00 0.0e+00 0.0e+00  0 17  0  0  0   0 17  0  0  0 14905
VecAssemblyBegin       3 1.0 2.1482e-04 1.3 0.00e+00 0.0 4.8e+01 7.8e+01 9.0e+00  0  0  2  1  6   0  0  2  1  6     0
VecAssemblyEnd         3 1.0 3.3140e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       53 1.0 2.5249e-04 1.6 0.00e+00 0.0 2.8e+03 9.9e+01 0.0e+00  0  0 89 85  0   0  0 89 85  0     0
VecScatterEnd         53 1.0 3.0153e-0330.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
VecNormalize          52 1.0 1.0890e-02 4.4 1.51e+04 1.4 0.0e+00 0.0e+00 5.2e+01  4  2  0  0 35   4  2  0  0 35    15
MatMult               52 1.0 3.4423e-03 7.3 1.47e+05 1.5 2.7e+03 9.8e+01 0.0e+00  1 17 88 82  0   1 17 88 82  0   434
MatSolve              53 1.0 3.8171e-04 2.3 4.89e+05 2.3 0.0e+00 0.0e+00 0.0e+00  0 40  0  0  0   0 40  0  0  0  9110
MatLUFactorNum         1 1.0 2.0480e-04 2.3 1.28e+05 4.7 0.0e+00 0.0e+00 0.0e+00  0  8  0  0  0   0  8  0  0  0  3194
MatILUFactorSym        1 1.0 5.6505e-04 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  0  0  0  0  2   0  0  0  0  2     0
MatAssemblyBegin       2 1.0 2.2840e-0314.9 0.00e+00 0.0 7.2e+01 5.3e+02 4.0e+00  1  0  2 12  3   1  0  2 12  3     0
MatAssemblyEnd         2 1.0 5.8198e-04 1.0 0.00e+00 0.0 1.0e+02 2.6e+01 8.0e+00  0  0  3  1  5   0  0  3  1  5     0
MatGetRowIJ            1 1.0 1.0967e-05 5.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 1.3518e-04 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 2.2e+00  0  0  0  0  1   0  0  0  0  1     0
MatZeroEntries         3 1.0 2.0027e-05 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        50 1.0 5.5857e-03 2.1 2.61e+05 1.4 0.0e+00 0.0e+00 5.0e+01  2 32  0  0 34   2 32  0  0 34   493
KSPSetUp               2 1.0 1.3208e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 2.6659e-02 1.0 1.05e+06 1.8 2.7e+03 9.8e+01 1.1e+02 15100 88 82 74  15100 88 82 75   324
PCSetUp                2 1.0 1.4973e-03 1.4 1.28e+05 4.7 0.0e+00 0.0e+00 7.2e+00  1  8  0  0  5   1  8  0  0  5   437
PCSetUpOnBlocks        1 1.0 1.0190e-03 1.6 1.28e+05 4.7 0.0e+00 0.0e+00 5.2e+00  0  8  0  0  3   0  8  0  0  4   642
PCApply               53 1.0 9.5224e-04 1.8 4.89e+05 2.3 0.0e+00 0.0e+00 0.0e+00  0 40  0  0  0   0 40  0  0  0  3652
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    43             43        94328     0
      Vector Scatter     2              2         2072     0
           Index Set     7              7         5964     0
   IS L to G Mapping     1              1          564     0
              Matrix     4              4        87984     0
       Krylov Solver     2              2        19360     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 1.19209e-07
Average time for MPI_Barrier(): 4.81606e-06
Average time for zero size MPI_Send(): 1.36693e-05
#PETSc Option Table entries:
-d 2
-ksp_right_pc
-log_summary
-n 15
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

 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.20184, Active time=0.160219                                                  |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0040      0.004021    0.0078      0.007757    2.51     4.84     |
|   build_sparsity()                 1         0.0021      0.002083    0.0071      0.007051    1.30     4.40     |
|   create_dof_constraints()         1         0.0007      0.000672    0.0007      0.000672    0.42     0.42     |
|   distribute_dofs()                1         0.0156      0.015638    0.0392      0.039239    9.76     24.49    |
|   dof_indices()                    97        0.0139      0.000144    0.0139      0.000144    8.70     8.70     |
|   prepare_send_list()              1         0.0001      0.000060    0.0001      0.000060    0.04     0.04     |
|   reinit()                         1         0.0216      0.021557    0.0216      0.021557    13.45    13.45    |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0007      0.000715    0.0078      0.007811    0.45     4.88     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0049      0.004919    0.0049      0.004919    3.07     3.07     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        29        0.0012      0.000043    0.0012      0.000043    0.78     0.78     |
|   init_shape_functions()           11        0.0003      0.000027    0.0003      0.000027    0.18     0.18     |
|   inverse_map()                    30        0.0005      0.000017    0.0005      0.000017    0.31     0.31     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             29        0.0005      0.000017    0.0005      0.000017    0.31     0.31     |
|   compute_face_map()               10        0.0004      0.000040    0.0009      0.000091    0.25     0.57     |
|   init_face_shape_functions()      1         0.0000      0.000018    0.0000      0.000018    0.01     0.01     |
|   init_reference_to_physical_map() 11        0.0006      0.000057    0.0006      0.000057    0.39     0.39     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0052      0.005179    0.0056      0.005599    3.23     3.49     |
|   renumber_nodes_and_elem()        2         0.0006      0.000283    0.0006      0.000283    0.35     0.35     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0039      0.001944    0.0039      0.001944    2.43     2.43     |
|   find_global_indices()            2         0.0019      0.000953    0.0102      0.005096    1.19     6.36     |
|   parallel_sort()                  2         0.0026      0.001313    0.0031      0.001549    1.64     1.93     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0002      0.000183    0.0130      0.013035    0.11     8.14     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0032      0.003172    0.0032      0.003172    1.98     1.98     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0242      0.024159    0.0285      0.028515    15.08    17.80    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      11        0.0010      0.000089    0.0011      0.000096    0.61     0.66     |
|   max(bool)                        1         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   max(scalar)                      115       0.0008      0.000007    0.0008      0.000007    0.49     0.49     |
|   max(vector)                      26        0.0003      0.000013    0.0008      0.000032    0.21     0.52     |
|   min(bool)                        133       0.0008      0.000006    0.0008      0.000006    0.51     0.51     |
|   min(scalar)                      109       0.0090      0.000083    0.0090      0.000083    5.64     5.64     |
|   min(vector)                      26        0.0004      0.000016    0.0010      0.000040    0.26     0.65     |
|   probe()                          132       0.0010      0.000007    0.0010      0.000007    0.61     0.61     |
|   receive()                        132       0.0008      0.000006    0.0018      0.000014    0.51     1.13     |
|   send()                           132       0.0004      0.000003    0.0004      0.000003    0.26     0.26     |
|   send_receive()                   136       0.0011      0.000008    0.0037      0.000027    0.70     2.33     |
|   sum()                            22        0.0006      0.000027    0.0040      0.000180    0.37     2.47     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           132       0.0003      0.000002    0.0003      0.000002    0.18     0.18     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0013      0.001303    0.0020      0.001959    0.81     1.22     |
|   set_parent_processor_ids()       1         0.0005      0.000482    0.0005      0.000482    0.30     0.30     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.0293      0.029251    0.0293      0.029251    18.26    18.26    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0037      0.003707    0.0107      0.010716    2.31     6.69     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1349      0.1602                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example subdomains_ex2:
*  mpirun -np 12 example-devel -d 2 -n 15 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
***************************************************************
* Running Example subdomains_ex2:
*  mpirun -np 12 example-devel -d 3 -n 6 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/subdomains/subdomains_ex2/.libs/lt-example-devel -d 3 -n 6 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=2197
    n_local_nodes()=245
  n_elem()=216
    n_local_elem()=18
    n_active_elem()=216
  n_subdomains()=2
  n_partitions()=12
  n_processors()=12
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
    n_local_dofs()=280
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 36.5321
      Average Off-Processor Bandwidth <= 14.5186
      Maximum  On-Processor Bandwidth <= 125
      Maximum Off-Processor Bandwidth <= 100
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=2197
    n_local_nodes()=245
  n_elem()=216
    n_local_elem()=18
    n_active_elem()=216
  n_subdomains()=2
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 22:12:37 2013                                                                             |
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
 -----------------------------------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.095968, Active time=0.095164                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           18        0.0397      0.002203    0.0397      0.002203    41.67    41.67    |
| Fe                            18        0.0010      0.000055    0.0010      0.000055    1.04     1.04     |
| Ke                            18        0.0283      0.001575    0.0283      0.001575    29.79    29.79    |
| elem init                     18        0.0255      0.001415    0.0255      0.001415    26.77    26.77    |
| matrix insertion              18        0.0007      0.000039    0.0007      0.000039    0.73     0.73     |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       90        0.0952                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/subdomains/subdomains_ex2/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:12:37 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           4.503e-01      1.00021   4.503e-01
Objects:              6.400e+01      1.03226   6.233e+01
Flops:                7.610e+06      5.19014   4.500e+06  5.400e+07
Flops/sec:            1.690e+07      5.19010   9.993e+06  1.199e+08
MPI Messages:         3.720e+02      2.61053   2.363e+02  2.836e+03
MPI Message Lengths:  1.561e+05      1.89868   4.439e+02  1.259e+06
MPI Reductions:       1.010e+02      1.02020

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 4.5027e-01 100.0%  5.4001e+07 100.0%  2.836e+03 100.0%  4.439e+02      100.0%  9.833e+01  99.3% 

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

VecMDot               26 1.0 5.1141e-03 3.3 1.96e+05 2.6 0.0e+00 0.0e+00 2.6e+01  1  3  0  0 26   1  3  0  0 26   345
VecNorm               28 1.0 1.6898e-02 8.5 1.57e+04 2.6 0.0e+00 0.0e+00 2.8e+01  2  0  0  0 28   2  0  0  0 28     8
VecScale              27 1.0 3.6240e-05 1.4 7.56e+03 2.6 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1879
VecCopy                2 1.0 5.2452e-06 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                34 1.0 2.9564e-05 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                2 1.0 1.4235e-02 2.5 1.12e+03 2.6 0.0e+00 0.0e+00 0.0e+00  2  0  0  0  0   2  0  0  0  0     1
VecMAXPY              27 1.0 1.2922e-04 2.2 2.11e+05 2.6 0.0e+00 0.0e+00 0.0e+00  0  4  0  0  0   0  4  0  0  0 14716
VecAssemblyBegin       3 1.0 5.4502e-04 1.2 0.00e+00 0.0 6.8e+01 3.1e+02 9.0e+00  0  0  2  2  9   0  0  2  2  9     0
VecAssemblyEnd         3 1.0 3.6716e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       28 1.0 2.3484e-04 1.7 0.00e+00 0.0 2.3e+03 2.9e+02 0.0e+00  0  0 81 54  0   0  0 81 54  0     0
VecScatterEnd         28 1.0 1.7843e-02130.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  2  0  0  0  0   2  0  0  0  0     0
VecNormalize          27 1.0 1.6800e-02 8.9 2.27e+04 2.6 0.0e+00 0.0e+00 2.7e+01  2  0  0  0 27   2  0  0  0 27    12
MatMult               27 1.0 1.8172e-0223.3 7.68e+05 2.7 2.2e+03 2.8e+02 0.0e+00  2 12 78 50  0   2 12 78 50  0   360
MatSolve              28 1.0 1.5588e-03 4.4 2.33e+06 4.5 0.0e+00 0.0e+00 0.0e+00  0 32  0  0  0   0 32  0  0  0 11180
MatLUFactorNum         1 1.0 4.5841e-03 7.8 4.08e+06 8.2 0.0e+00 0.0e+00 0.0e+00  1 48  0  0  0   1 48  0  0  0  5704
MatILUFactorSym        1 1.0 1.4646e-02 7.5 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  2  0  0  0  3   2  0  0  0  3     0
MatAssemblyBegin       2 1.0 2.3389e-0285.7 0.00e+00 0.0 1.0e+02 5.1e+03 4.0e+00  3  0  4 41  4   3  0  4 41  4     0
MatAssemblyEnd         2 1.0 1.0920e-03 1.2 0.00e+00 0.0 1.6e+02 7.3e+01 8.0e+00  0  0  6  1  8   0  0  6  1  8     0
MatGetRowIJ            1 1.0 1.1921e-05 3.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 1.4782e-04 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 2.3e+00  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries         3 1.0 7.5817e-05 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        26 1.0 5.2195e-03 3.1 3.93e+05 2.6 0.0e+00 0.0e+00 2.6e+01  1  7  0  0 26   1  7  0  0 26   678
KSPSetUp               2 1.0 1.1802e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 4.4281e-02 1.0 7.61e+06 5.2 2.2e+03 2.8e+02 6.1e+01 10100 78 50 62  10100 78 50 62  1220
PCSetUp                2 1.0 2.0015e-02 6.0 4.08e+06 8.2 0.0e+00 0.0e+00 7.3e+00  2 48  0  0  7   2 48  0  0  7  1306
PCSetUpOnBlocks        1 1.0 1.9500e-02 6.9 4.08e+06 8.2 0.0e+00 0.0e+00 5.3e+00  2 48  0  0  5   2 48  0  0  5  1341
PCApply               28 1.0 1.9162e-03 2.7 2.33e+06 4.5 0.0e+00 0.0e+00 0.0e+00  0 32  0  0  0   0 32  0  0  0  9095
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    43             43       154792     0
      Vector Scatter     2              2         2072     0
           Index Set     7              7         9112     0
   IS L to G Mapping     1              1          564     0
              Matrix     4              4       709268     0
       Krylov Solver     2              2        19360     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 4.3869e-06
Average time for zero size MPI_Send(): 1.35104e-05
#PETSc Option Table entries:
-d 3
-ksp_right_pc
-log_summary
-n 6
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

 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.474286, Active time=0.425871                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0113      0.011339    0.0365      0.036489    2.66     8.57     |
|   build_sparsity()                 1         0.0082      0.008224    0.0231      0.023132    1.93     5.43     |
|   create_dof_constraints()         1         0.0007      0.000690    0.0007      0.000690    0.16     0.16     |
|   distribute_dofs()                1         0.0346      0.034595    0.0863      0.086314    8.12     20.27    |
|   dof_indices()                    123       0.0637      0.000518    0.0637      0.000518    14.97    14.97    |
|   prepare_send_list()              1         0.0004      0.000368    0.0004      0.000368    0.09     0.09     |
|   reinit()                         1         0.0485      0.048480    0.0485      0.048480    11.38    11.38    |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0009      0.000874    0.0333      0.033321    0.21     7.82     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0071      0.007057    0.0071      0.007057    1.66     1.66     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        45        0.0134      0.000297    0.0134      0.000297    3.14     3.14     |
|   init_shape_functions()           28        0.0012      0.000042    0.0012      0.000042    0.28     0.28     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             45        0.0029      0.000065    0.0029      0.000065    0.69     0.69     |
|   compute_face_map()               27        0.0012      0.000045    0.0012      0.000045    0.29     0.29     |
|   init_face_shape_functions()      1         0.0002      0.000183    0.0002      0.000183    0.04     0.04     |
|   init_reference_to_physical_map() 28        0.0242      0.000864    0.0242      0.000864    5.68     5.68     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 1         0.0076      0.007633    0.0078      0.007760    1.79     1.82     |
|   renumber_nodes_and_elem()        2         0.0013      0.000647    0.0013      0.000647    0.30     0.30     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        2         0.0039      0.001951    0.0039      0.001951    0.92     0.92     |
|   find_global_indices()            2         0.0019      0.000936    0.0102      0.005104    0.44     2.40     |
|   parallel_sort()                  2         0.0027      0.001327    0.0031      0.001554    0.62     0.73     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0001      0.000139    0.0416      0.041647    0.03     9.78     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0062      0.006248    0.0062      0.006248    1.47     1.47     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      1         0.0379      0.037924    0.0423      0.042305    8.91     9.93     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      11        0.0009      0.000080    0.0009      0.000085    0.21     0.22     |
|   max(bool)                        1         0.0000      0.000007    0.0000      0.000007    0.00     0.00     |
|   max(scalar)                      115       0.0017      0.000015    0.0017      0.000015    0.41     0.41     |
|   max(vector)                      26        0.0005      0.000020    0.0015      0.000059    0.12     0.36     |
|   min(bool)                        133       0.0019      0.000015    0.0019      0.000015    0.46     0.46     |
|   min(scalar)                      109       0.0282      0.000259    0.0282      0.000259    6.63     6.63     |
|   min(vector)                      26        0.0006      0.000025    0.0020      0.000075    0.15     0.46     |
|   probe()                          132       0.0017      0.000013    0.0017      0.000013    0.41     0.41     |
|   receive()                        132       0.0009      0.000007    0.0026      0.000020    0.20     0.62     |
|   send()                           132       0.0005      0.000004    0.0005      0.000004    0.12     0.12     |
|   send_receive()                   136       0.0014      0.000010    0.0049      0.000036    0.33     1.16     |
|   sum()                            22        0.0019      0.000087    0.0113      0.000513    0.45     2.65     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           132       0.0003      0.000002    0.0003      0.000002    0.07     0.07     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         1         0.0020      0.001953    0.0027      0.002676    0.46     0.63     |
|   set_parent_processor_ids()       1         0.0004      0.000448    0.0004      0.000448    0.11     0.11     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.0589      0.058925    0.0589      0.058925    13.84    13.84    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0438      0.043766    0.0963      0.096346    10.28    22.62    |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1428      0.4259                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example subdomains_ex2:
*  mpirun -np 12 example-devel -d 3 -n 6 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
