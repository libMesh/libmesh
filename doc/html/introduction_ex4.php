<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("introduction_ex4",$root)?>
 
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
<br><br><br> <h1> The source file introduction_ex4.C with comments: </h1> 
<div class = "comment">
<h1>Introduction Example 4 - Solving a 1D, 2D or 3D Poisson Problem in Parallel</h1>

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
parallel.  Notice how little has changed!  The significant
differences are marked with "PARALLEL CHANGE".


<br><br>

<br><br>C++ include files that we need
</div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
        #include &lt;algorithm&gt;
        #include &lt;math.h&gt;
        #include &lt;set&gt;
        
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
        
</pre>
</div>
<div class = "comment">
To impose Dirichlet boundary conditions
</div>

<div class ="fragment">
<pre>
        #include "libmesh/dirichlet_boundaries.h"
        #include "libmesh/analytic_function.h"
        
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
        		     const Real y,
        		     const Real z = 0.);
        
</pre>
</div>
<div class = "comment">
Define a wrapper for exact_solution that will be needed below
</div>

<div class ="fragment">
<pre>
        void exact_solution_wrapper (DenseVector&lt;Number&gt;& output,
                                     const Point& p,
                                     const Real)
        {
          output(0) = exact_solution(p(0),
        			     (LIBMESH_DIM&gt;1)?p(1):0,
        			     (LIBMESH_DIM&gt;2)?p(2):0);
        }
        
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
                std::cerr &lt;&lt; "ex4 currently requires a C^0 (or higher) FE basis." &lt;&lt; std::endl;
              libmesh_error();
            }
            
          Mesh mesh;
          
        
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
        
          
</pre>
</div>
<div class = "comment">
Add the variable "u" to "Poisson".  "u"
will be approximated using second-order approximation.
</div>

<div class ="fragment">
<pre>
          unsigned int u_var = system.add_variable("u",
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
Construct a Dirichlet boundary condition object
  

<br><br>Indicate which boundary IDs we impose the BC on
We either build a line, a square or a cube, and
here we indicate the boundaries IDs in each case
</div>

<div class ="fragment">
<pre>
          std::set&lt;boundary_id_type&gt; boundary_ids;
</pre>
</div>
<div class = "comment">
the dim==1 mesh has two boundaries with IDs 0 and 1
</div>

<div class ="fragment">
<pre>
          boundary_ids.insert(0);
          boundary_ids.insert(1);
</pre>
</div>
<div class = "comment">
the dim==2 mesh has four boundaries with IDs 0, 1, 2 and 3
</div>

<div class ="fragment">
<pre>
          if(dim&gt;=2)
          {
            boundary_ids.insert(2);
            boundary_ids.insert(3);
          }
</pre>
</div>
<div class = "comment">
the dim==3 mesh has four boundaries with IDs 0, 1, 2, 3, 4 and 5
</div>

<div class ="fragment">
<pre>
          if(dim==3)
          {
            boundary_ids.insert(4);
            boundary_ids.insert(5);
          }
        
</pre>
</div>
<div class = "comment">
Create a vector storing the variable numbers which the BC applies to
</div>

<div class ="fragment">
<pre>
          std::vector&lt;unsigned int&gt; variables(1);
          variables[0] = u_var;
          
</pre>
</div>
<div class = "comment">
Create an AnalyticFunction object that we use to project the BC
This function just calls the function exact_solution via exact_solution_wrapper
</div>

<div class ="fragment">
<pre>
          AnalyticFunction&lt;&gt; exact_solution_object(exact_solution_wrapper);
          
          DirichletBoundary dirichlet_bc(boundary_ids,
                                         variables,
                                         &exact_solution_object);
        
</pre>
</div>
<div class = "comment">
We must add the Dirichlet boundary condition _before_ 
we call equation_systems.init()
</div>

<div class ="fragment">
<pre>
          system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
        
        
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
          system.solve();
        
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
            GnuPlotIO plot(mesh,"Introduction Example 4, 1D",GnuPlotIO::GRID_ON);
            plot.write_equation_systems("gnuplot_script",equation_systems);
          }
        #ifdef LIBMESH_HAVE_EXODUS_API
          else
          {
            ExodusII_IO (mesh).write_equation_systems ((dim == 3) ? 
              "out_3.e" : "out_2.e",equation_systems);
          }
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
          
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
We now define the matrix assembly function for the
Poisson system.  We need to first compute element
matrices and right-hand sides, and then take into
account the boundary conditions.
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
          std::vector&lt;dof_id_type&gt; dof_indices;
        
</pre>
</div>
<div class = "comment">
Now we will loop over all the elements in the mesh.
We will compute the element matrix and right-hand-side
contribution.  See example 3 for a discussion of the
element iterators.
</div>

<div class ="fragment">
<pre>
          MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
          const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
        
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
                  const Real y = 0.;
        #endif
        #if LIBMESH_DIM &gt; 2
                  const Real z = q_point[qp](2);
        #else
                  const Real z = 0.;
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
If this assembly program were to be used on an adaptive mesh,
we would have to apply any hanging node constraint equations
Also, note that here we call heterogenously_constrain_element_matrix_and_vector
to impose a inhomogeneous Dirichlet boundary conditions.
</div>

<div class ="fragment">
<pre>
              dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        
</pre>
</div>
<div class = "comment">
The element matrix and right-hand-side are now built
for this element.  Add them to the global matrix and
right-hand-side vector.  The \p SparseMatrix::add_matrix()
and \p NumericVector::add_vector() members do this for us.
Start logging the insertion of the local (element)
matrix and vector into the global matrix and vector
</div>

<div class ="fragment">
<pre>
              perf_log.push ("matrix insertion");
              
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
<br><br><br> <h1> The source file introduction_ex4.C without comments: </h1> 
<pre> 
  
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  #include &lt;set&gt;
  
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
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dirichlet_boundaries.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/analytic_function.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/string_to_enum.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_poisson(EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  Real exact_solution (<B><FONT COLOR="#228B22">const</FONT></B> Real x,
  		     <B><FONT COLOR="#228B22">const</FONT></B> Real y,
  		     <B><FONT COLOR="#228B22">const</FONT></B> Real z = 0.);
  
  <B><FONT COLOR="#228B22">void</FONT></B> exact_solution_wrapper (DenseVector&lt;Number&gt;&amp; output,
                               <B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                               <B><FONT COLOR="#228B22">const</FONT></B> Real)
  {
    output(0) = exact_solution(p(0),
  			     (LIBMESH_DIM&gt;1)?p(1):0,
  			     (LIBMESH_DIM&gt;2)?p(2):0);
  }
  
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
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;ex4 currently requires a C^0 (or higher) FE basis.&quot;</FONT></B> &lt;&lt; std::endl;
        libmesh_error();
      }
      
    Mesh mesh;
    
  
  
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
  
    
    mesh.print_info();
    
    
    EquationSystems equation_systems (mesh);
    
    LinearImplicitSystem&amp; system =
      equation_systems.add_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>);
  
    
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>,
                                             <B><FONT COLOR="#5F9EA0">Utility</FONT></B>::string_to_enum&lt;Order&gt;   (order),
                                             <B><FONT COLOR="#5F9EA0">Utility</FONT></B>::string_to_enum&lt;FEFamily&gt;(family));
  
    system.attach_assemble_function (assemble_poisson);
  
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;boundary_id_type&gt; boundary_ids;
    boundary_ids.insert(0);
    boundary_ids.insert(1);
    <B><FONT COLOR="#A020F0">if</FONT></B>(dim&gt;=2)
    {
      boundary_ids.insert(2);
      boundary_ids.insert(3);
    }
    <B><FONT COLOR="#A020F0">if</FONT></B>(dim==3)
    {
      boundary_ids.insert(4);
      boundary_ids.insert(5);
    }
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; variables(1);
    variables[0] = u_var;
    
    AnalyticFunction&lt;&gt; exact_solution_object(exact_solution_wrapper);
    
    DirichletBoundary dirichlet_bc(boundary_ids,
                                   variables,
                                   &amp;exact_solution_object);
  
    system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);
  
  
    equation_systems.init();
  
    equation_systems.print_info();
    mesh.print_info();
  
    system.solve();
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(dim == 1)
    {        
      GnuPlotIO plot(mesh,<B><FONT COLOR="#BC8F8F">&quot;Introduction Example 4, 1D&quot;</FONT></B>,GnuPlotIO::GRID_ON);
      plot.write_equation_systems(<B><FONT COLOR="#BC8F8F">&quot;gnuplot_script&quot;</FONT></B>,equation_systems);
    }
  #ifdef LIBMESH_HAVE_EXODUS_API
    <B><FONT COLOR="#A020F0">else</FONT></B>
    {
      ExodusII_IO (mesh).write_equation_systems ((dim == 3) ? 
        <B><FONT COLOR="#BC8F8F">&quot;out_3.e&quot;</FONT></B> : <B><FONT COLOR="#BC8F8F">&quot;out_2.e&quot;</FONT></B>,equation_systems);
    }
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>    
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
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        perf_log.push(<B><FONT COLOR="#BC8F8F">&quot;elem init&quot;</FONT></B>);      
  
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
        dof_map.dof_indices (elem, dof_indices);
  
        fe-&gt;reinit (elem);
  
        Ke.resize (dof_indices.size(),
                   dof_indices.size());
  
        Fe.resize (dof_indices.size());
  
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
            <B><FONT COLOR="#228B22">const</FONT></B> Real y = 0.;
  #endif
  #<B><FONT COLOR="#A020F0">if</FONT></B> LIBMESH_DIM &gt; 2
            <B><FONT COLOR="#228B22">const</FONT></B> Real z = q_point[qp](2);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
            <B><FONT COLOR="#228B22">const</FONT></B> Real z = 0.;
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
        
        dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
  
        perf_log.push (<B><FONT COLOR="#BC8F8F">&quot;matrix insertion&quot;</FONT></B>);
        
        system.matrix-&gt;add_matrix (Ke, dof_indices);
        system.rhs-&gt;add_vector    (Fe, dof_indices);
  
        perf_log.pop (<B><FONT COLOR="#BC8F8F">&quot;matrix insertion&quot;</FONT></B>);
      }
  
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example introduction_ex4:
*  mpirun -np 12 example-devel -d 1 -n 20 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/introduction/introduction_ex4/.libs/lt-example-devel -d 1 -n 20 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=1
  spatial_dimension()=3
  n_nodes()=41
    n_local_nodes()=5
  n_elem()=20
    n_local_elem()=2
    n_active_elem()=20
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Poisson"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" 
    n_dofs()=41
    n_local_dofs()=5
    n_constrained_dofs()=2
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 3.34146
      Average Off-Processor Bandwidth <= 0.780488
      Maximum  On-Processor Bandwidth <= 5
      Maximum Off-Processor Bandwidth <= 2
    DofMap Constraints
      Number of DoF Constraints = 2
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=1
  spatial_dimension()=3
  n_nodes()=41
    n_local_nodes()=5
  n_elem()=20
    n_local_elem()=2
    n_active_elem()=20
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 21:56:46 2013                                                                             |
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
| Matrix Assembly Performance: Alive time=0.001329, Active time=0.000643                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| Fe                            2         0.0000      0.000005    0.0000      0.000005    1.56     1.56     |
| Ke                            2         0.0000      0.000004    0.0000      0.000004    1.24     1.24     |
| elem init                     2         0.0005      0.000273    0.0005      0.000273    84.76    84.76    |
| matrix insertion              2         0.0001      0.000040    0.0001      0.000040    12.44    12.44    |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       8         0.0006                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/introduction/introduction_ex4/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 21:56:46 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           7.269e-02      1.00011   7.268e-02
Objects:              5.500e+01      1.03774   5.333e+01
Flops:                6.910e+03      0.00000   3.199e+03  3.839e+04
Flops/sec:            9.507e+04      0.00000   4.402e+04  5.282e+05
MPI Messages:         5.100e+01      0.00000   3.400e+01  4.080e+02
MPI Message Lengths:  5.840e+02      0.00000   1.145e+01  4.672e+03
MPI Reductions:       8.300e+01      1.02469

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 7.2619e-02  99.9%  3.8393e+04 100.0%  4.080e+02 100.0%  1.145e+01      100.0%  8.033e+01  96.8% 

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

VecMDot               17 1.0 2.4676e-04 1.8 1.99e+03 0.0 0.0e+00 0.0e+00 1.7e+01  0 29  0  0 20   0 29  0  0 21    45
VecNorm               19 1.0 9.2232e-0372.6 2.66e+02 0.0 0.0e+00 0.0e+00 1.9e+01  2  4  0  0 23   2  4  0  0 24     0
VecScale              18 1.0 3.1233e-05 1.3 1.26e+02 0.0 0.0e+00 0.0e+00 0.0e+00  0  2  0  0  0   0  2  0  0  0    24
VecCopy                2 1.0 7.8678e-06 8.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                25 1.0 1.3351e-05 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                2 1.0 9.1820e-03664.0 2.80e+01 0.0 0.0e+00 0.0e+00 0.0e+00 10  0  0  0  0  10  0  0  0  0     0
VecMAXPY              18 1.0 6.9141e-06 3.6 2.38e+03 0.0 0.0e+00 0.0e+00 0.0e+00  0 36  0  0  0   0 36  0  0  0  2016
VecAssemblyBegin       3 1.0 1.6499e-04 1.1 0.00e+00 0.0 1.6e+01 6.0e+00 9.0e+00  0  0  4  2 11   0  0  4  2 11     0
VecAssemblyEnd         3 1.0 3.0279e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       19 1.0 9.8705e-0513.8 0.00e+00 0.0 3.0e+02 1.2e+01 0.0e+00  0  0 75 81  0   0  0 75 81  0     0
VecScatterEnd         19 1.0 4.4346e-0511.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize          18 1.0 9.3262e-0342.8 3.78e+02 0.0 0.0e+00 0.0e+00 1.8e+01  3  6  0  0 22   3  6  0  0 22     0
MatMult               18 1.0 1.7524e-04 3.5 9.18e+02 0.0 2.9e+02 1.2e+01 0.0e+00  0 13 71 74  0   0 13 71 74  0    29
MatSolve              19 0.0 1.6212e-05 0.0 1.12e+03 0.0 0.0e+00 0.0e+00 0.0e+00  0 14  0  0  0   0 14  0  0  0   334
MatLUFactorNum         1 1.0 5.1975e-05 2.6 8.20e+01 0.0 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0     7
MatILUFactorSym        1 1.0 1.0014e-04 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  0  0  0  0  4   0  0  0  0  4     0
MatAssemblyBegin       2 1.0 2.3253e-02106.1 0.00e+00 0.0 2.4e+01 1.7e+01 4.0e+00 25  0  6  9  5  25  0  6  9  5     0
MatAssemblyEnd         2 1.0 4.4608e-04 1.0 0.00e+00 0.0 3.2e+01 5.0e+00 8.0e+00  1  0  8  3 10   1  0  8  3 10     0
MatGetRowIJ            1 0.0 5.0068e-06 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 0.0 1.2589e-04 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.8e+00  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries         3 0.0 1.3113e-05 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        17 1.0 2.9087e-04 1.5 4.13e+03 0.0 0.0e+00 0.0e+00 1.7e+01  0 62  0  0 20   0 62  0  0 21    82
KSPSetUp               2 1.0 1.1516e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 1.1387e-02 1.0 6.91e+03 0.0 2.9e+02 1.2e+01 4.3e+01 16100 71 74 52  16100 71 74 54     3
PCSetUp                2 1.0 8.5092e-04 1.0 8.20e+01 0.0 0.0e+00 0.0e+00 7.3e+00  1  1  0  0  9   1  1  0  0  9     0
PCSetUpOnBlocks        1 1.0 3.7289e-04 1.1 8.20e+01 0.0 0.0e+00 0.0e+00 5.3e+00  0  1  0  0  6   0  1  0  0  7     1
PCApply               19 1.0 2.7800e-04 1.2 1.12e+03 0.0 0.0e+00 0.0e+00 0.0e+00  0 14  0  0  0   0 14  0  0  0    19
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
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 8.58307e-06
Average time for zero size MPI_Send(): 1.37488e-05
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

 --------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.118531, Active time=0.060117                                                     |
 --------------------------------------------------------------------------------------------------------------------
| Event                                  nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                  w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------------|
|                                                                                                                    |
|                                                                                                                    |
| DofMap                                                                                                             |
|   add_neighbors_to_send_list()         1         0.0003      0.000265    0.0004      0.000449    0.44     0.75     |
|   build_constraint_matrix_and_vector() 2         0.0000      0.000005    0.0000      0.000005    0.02     0.02     |
|   build_sparsity()                     1         0.0006      0.000581    0.0017      0.001692    0.97     2.81     |
|   create_dof_constraints()             1         0.0012      0.001197    0.0021      0.002063    1.99     3.43     |
|   distribute_dofs()                    1         0.0011      0.001068    0.0034      0.003423    1.78     5.69     |
|   dof_indices()                        30        0.0013      0.000045    0.0013      0.000045    2.22     2.22     |
|   hetero_cnstrn_elem_mat_vec()         2         0.0000      0.000002    0.0000      0.000002    0.01     0.01     |
|   prepare_send_list()                  1         0.0000      0.000017    0.0000      0.000017    0.03     0.03     |
|   reinit()                             1         0.0011      0.001109    0.0011      0.001109    1.84     1.84     |
|                                                                                                                    |
| EquationSystems                                                                                                    |
|   build_solution_vector()              1         0.0003      0.000321    0.0009      0.000864    0.53     1.44     |
|                                                                                                                    |
| FE                                                                                                                 |
|   compute_shape_functions()            2         0.0001      0.000027    0.0001      0.000027    0.09     0.09     |
|   init_shape_functions()               1         0.0002      0.000181    0.0002      0.000181    0.30     0.30     |
|                                                                                                                    |
| FEMap                                                                                                              |
|   compute_affine_map()                 2         0.0000      0.000022    0.0000      0.000022    0.07     0.07     |
|   init_reference_to_physical_map()     1         0.0001      0.000084    0.0001      0.000084    0.14     0.14     |
|                                                                                                                    |
| GnuPlotIO                                                                                                          |
|   write_nodal_data()                   1         0.0012      0.001183    0.0012      0.001183    1.97     1.97     |
|                                                                                                                    |
| Mesh                                                                                                               |
|   find_neighbors()                     1         0.0011      0.001137    0.0013      0.001271    1.89     2.11     |
|   renumber_nodes_and_elem()            2         0.0001      0.000065    0.0001      0.000065    0.22     0.22     |
|                                                                                                                    |
| MeshCommunication                                                                                                  |
|   compute_hilbert_indices()            2         0.0005      0.000253    0.0005      0.000253    0.84     0.84     |
|   find_global_indices()                2         0.0008      0.000384    0.0051      0.002526    1.28     8.41     |
|   parallel_sort()                      2         0.0022      0.001112    0.0026      0.001282    3.70     4.26     |
|                                                                                                                    |
| MeshOutput                                                                                                         |
|   write_equation_systems()             1         0.0002      0.000183    0.0024      0.002357    0.30     3.92     |
|                                                                                                                    |
| MeshTools::Generation                                                                                              |
|   build_cube()                         1         0.0005      0.000498    0.0005      0.000498    0.83     0.83     |
|                                                                                                                    |
| MetisPartitioner                                                                                                   |
|   partition()                          1         0.0029      0.002895    0.0050      0.005018    4.82     8.35     |
|                                                                                                                    |
| Parallel                                                                                                           |
|   allgather()                          11        0.0004      0.000035    0.0005      0.000043    0.63     0.78     |
|   max(bool)                            1         0.0000      0.000007    0.0000      0.000007    0.01     0.01     |
|   max(scalar)                          113       0.0008      0.000007    0.0008      0.000007    1.34     1.34     |
|   max(vector)                          26        0.0004      0.000014    0.0009      0.000035    0.60     1.53     |
|   min(bool)                            131       0.0009      0.000007    0.0009      0.000007    1.42     1.42     |
|   min(scalar)                          107       0.0017      0.000016    0.0017      0.000016    2.78     2.78     |
|   min(vector)                          26        0.0004      0.000016    0.0013      0.000050    0.71     2.16     |
|   probe()                              132       0.0006      0.000004    0.0006      0.000004    0.97     0.97     |
|   receive()                            132       0.0008      0.000006    0.0014      0.000010    1.27     2.27     |
|   send()                               132       0.0004      0.000003    0.0004      0.000003    0.64     0.64     |
|   send_receive()                       136       0.0011      0.000008    0.0032      0.000023    1.76     5.24     |
|   sum()                                20        0.0004      0.000019    0.0006      0.000032    0.65     1.07     |
|                                                                                                                    |
| Parallel::Request                                                                                                  |
|   wait()                               132       0.0003      0.000002    0.0003      0.000002    0.44     0.44     |
|                                                                                                                    |
| Partitioner                                                                                                        |
|   set_node_processor_ids()             1         0.0005      0.000536    0.0011      0.001053    0.89     1.75     |
|   set_parent_processor_ids()           1         0.0001      0.000122    0.0001      0.000122    0.20     0.20     |
|                                                                                                                    |
| PetscLinearSolver                                                                                                  |
|   solve()                              1         0.0346      0.034556    0.0346      0.034556    57.48    57.48    |
|                                                                                                                    |
| System                                                                                                             |
|   assemble()                           1         0.0012      0.001159    0.0016      0.001639    1.93     2.73     |
 --------------------------------------------------------------------------------------------------------------------
| Totals:                                1163      0.0601                                          100.00            |
 --------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example introduction_ex4:
*  mpirun -np 12 example-devel -d 1 -n 20 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
***************************************************************
* Running Example introduction_ex4:
*  mpirun -np 12 example-devel -d 2 -n 15 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/introduction/introduction_ex4/.libs/lt-example-devel -d 2 -n 15 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=961
    n_local_nodes()=97
  n_elem()=225
    n_local_elem()=19
    n_active_elem()=225
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Poisson"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" 
    n_dofs()=961
    n_local_dofs()=97
    n_constrained_dofs()=120
    n_local_constrained_dofs()=21
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 13.2112
      Average Off-Processor Bandwidth <= 2.7846
      Maximum  On-Processor Bandwidth <= 26
      Maximum Off-Processor Bandwidth <= 18
    DofMap Constraints
      Number of DoF Constraints = 120
      Number of Heterogenous Constraints= 118
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=961
    n_local_nodes()=97
  n_elem()=225
    n_local_elem()=19
    n_active_elem()=225
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 21:56:46 2013                                                                             |
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
| Matrix Assembly Performance: Alive time=0.015138, Active time=0.00632                                     |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| Fe                            19        0.0002      0.000012    0.0002      0.000012    3.69     3.69     |
| Ke                            19        0.0013      0.000068    0.0013      0.000068    20.36    20.36    |
| elem init                     19        0.0046      0.000243    0.0046      0.000243    73.02    73.02    |
| matrix insertion              19        0.0002      0.000010    0.0002      0.000010    2.93     2.93     |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       76        0.0063                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/introduction/introduction_ex4/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 21:56:46 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           2.300e-01      1.00001   2.300e-01
Objects:              6.200e+01      1.00000   6.200e+01
Flops:                1.505e+06      1.69736   1.049e+06  1.259e+07
Flops/sec:            6.542e+06      1.69735   4.560e+06  5.472e+07
MPI Messages:         6.030e+02      3.48555   3.742e+02  4.491e+03
MPI Message Lengths:  5.214e+04      1.96532   9.941e+01  4.464e+05
MPI Reductions:       2.010e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.2998e-01 100.0%  1.2587e+07 100.0%  4.491e+03 100.0%  9.941e+01      100.0%  2.000e+02  99.5% 

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

VecMDot               76 1.0 8.3017e-04 1.5 2.06e+05 1.4 0.0e+00 0.0e+00 7.6e+01  0 16  0  0 38   0 16  0  0 38  2453
VecNorm               80 1.0 8.9016e-0317.1 1.55e+04 1.4 0.0e+00 0.0e+00 8.0e+01  3  1  0  0 40   3  1  0  0 40    17
VecScale              79 1.0 7.2718e-05 1.5 7.66e+03 1.4 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  1044
VecCopy                4 1.0 5.2452e-06 2.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                88 1.0 5.6982e-05 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                6 1.0 8.4088e-03641.3 1.16e+03 1.4 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     1
VecMAXPY              79 1.0 1.5283e-04 1.3 2.22e+05 1.4 0.0e+00 0.0e+00 0.0e+00  0 17  0  0  0   0 17  0  0  0 14362
VecAssemblyBegin       3 1.0 6.9189e-04 1.8 0.00e+00 0.0 4.8e+01 7.8e+01 9.0e+00  0  0  1  1  4   0  0  1  1  4     0
VecAssemblyEnd         3 1.0 2.5034e-05 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       80 1.0 3.4809e-04 1.7 0.00e+00 0.0 4.2e+03 9.5e+01 0.0e+00  0  0 93 89  0   0  0 93 89  0     0
VecScatterEnd         80 1.0 6.1369e-04 4.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize          79 1.0 9.0511e-0313.8 2.30e+04 1.4 0.0e+00 0.0e+00 7.9e+01  3  2  0  0 39   3  2  0  0 40    25
MatMult               79 1.0 1.2503e-03 1.9 2.24e+05 1.5 4.1e+03 9.5e+01 0.0e+00  0 18 91 87  0   0 18 91 87  0  1790
MatSolve              80 1.0 5.6934e-04 1.7 7.38e+05 2.0 0.0e+00 0.0e+00 0.0e+00  0 42  0  0  0   0 42  0  0  0  9364
MatLUFactorNum         1 1.0 1.8597e-04 2.1 9.16e+04 3.8 0.0e+00 0.0e+00 0.0e+00  0  4  0  0  0   0  4  0  0  0  2937
MatILUFactorSym        1 1.0 5.6314e-04 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  0  0  0  0  1   0  0  0  0  2     0
MatAssemblyBegin       2 1.0 2.5777e-0237.8 0.00e+00 0.0 7.2e+01 5.3e+02 4.0e+00  7  0  2  9  2   7  0  2  9  2     0
MatAssemblyEnd         2 1.0 1.1120e-03 1.8 0.00e+00 0.0 1.0e+02 2.6e+01 8.0e+00  0  0  2  1  4   0  0  2  1  4     0
MatGetRowIJ            1 1.0 2.1458e-06 2.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 8.1062e-05 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  1   0  0  0  0  1     0
MatZeroEntries         3 1.0 2.3842e-05 2.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        76 1.0 1.0660e-03 1.3 4.13e+05 1.4 0.0e+00 0.0e+00 7.6e+01  0 32  0  0 38   0 32  0  0 38  3832
KSPSetUp               2 1.0 1.0586e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 1.3735e-02 1.0 1.51e+06 1.7 4.1e+03 9.5e+01 1.6e+02  6100 91 87 81   6100 91 87 82   916
PCSetUp                2 1.0 1.3542e-03 1.4 9.16e+04 3.8 0.0e+00 0.0e+00 7.0e+00  0  4  0  0  3   0  4  0  0  4   403
PCSetUpOnBlocks        1 1.0 9.6917e-04 1.6 9.16e+04 3.8 0.0e+00 0.0e+00 5.0e+00  0  4  0  0  2   0  4  0  0  2   564
PCApply               80 1.0 1.3795e-03 1.3 7.38e+05 2.0 0.0e+00 0.0e+00 0.0e+00  1 42  0  0  0   1 42  0  0  0  3865
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
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.3869e-06
Average time for zero size MPI_Send(): 1.32521e-05
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

 --------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.279241, Active time=0.21264                                                      |
 --------------------------------------------------------------------------------------------------------------------
| Event                                  nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                  w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------------|
|                                                                                                                    |
|                                                                                                                    |
| DofMap                                                                                                             |
|   add_neighbors_to_send_list()         1         0.0025      0.002469    0.0049      0.004930    1.16     2.32     |
|   build_constraint_matrix_and_vector() 19        0.0004      0.000019    0.0004      0.000019    0.17     0.17     |
|   build_sparsity()                     1         0.0017      0.001713    0.0082      0.008157    0.81     3.84     |
|   create_dof_constraints()             1         0.0083      0.008293    0.0446      0.044580    3.90     20.97    |
|   distribute_dofs()                    1         0.0114      0.011378    0.0324      0.032366    5.35     15.22    |
|   dof_indices()                        303       0.0359      0.000119    0.0359      0.000119    16.91    16.91    |
|   hetero_cnstrn_elem_mat_vec()         19        0.0076      0.000401    0.0076      0.000401    3.59     3.59     |
|   prepare_send_list()                  1         0.0001      0.000057    0.0001      0.000057    0.03     0.03     |
|   reinit()                             1         0.0189      0.018900    0.0189      0.018900    8.89     8.89     |
|                                                                                                                    |
| EquationSystems                                                                                                    |
|   build_solution_vector()              1         0.0006      0.000591    0.0042      0.004189    0.28     1.97     |
|                                                                                                                    |
| ExodusII_IO                                                                                                        |
|   write_nodal_data()                   1         0.0046      0.004596    0.0046      0.004596    2.16     2.16     |
|                                                                                                                    |
| FE                                                                                                                 |
|   compute_shape_functions()            79        0.0019      0.000024    0.0019      0.000024    0.89     0.89     |
|   init_shape_functions()               61        0.0004      0.000007    0.0004      0.000007    0.20     0.20     |
|   inverse_map()                        180       0.0029      0.000016    0.0029      0.000016    1.35     1.35     |
|                                                                                                                    |
| FEMap                                                                                                              |
|   compute_affine_map()                 79        0.0010      0.000013    0.0010      0.000013    0.49     0.49     |
|   compute_face_map()                   60        0.0019      0.000032    0.0048      0.000080    0.90     2.26     |
|   init_face_shape_functions()          60        0.0004      0.000007    0.0004      0.000007    0.19     0.19     |
|   init_reference_to_physical_map()     61        0.0029      0.000047    0.0029      0.000047    1.35     1.35     |
|                                                                                                                    |
| Mesh                                                                                                               |
|   find_neighbors()                     1         0.0051      0.005123    0.0056      0.005563    2.41     2.62     |
|   renumber_nodes_and_elem()            2         0.0006      0.000279    0.0006      0.000279    0.26     0.26     |
|                                                                                                                    |
| MeshCommunication                                                                                                  |
|   compute_hilbert_indices()            2         0.0039      0.001956    0.0039      0.001956    1.84     1.84     |
|   find_global_indices()                2         0.0020      0.000976    0.0102      0.005076    0.92     4.77     |
|   parallel_sort()                      2         0.0026      0.001300    0.0030      0.001520    1.22     1.43     |
|                                                                                                                    |
| MeshOutput                                                                                                         |
|   write_equation_systems()             1         0.0002      0.000177    0.0096      0.009634    0.08     4.53     |
|                                                                                                                    |
| MeshTools::Generation                                                                                              |
|   build_cube()                         1         0.0029      0.002875    0.0029      0.002875    1.35     1.35     |
|                                                                                                                    |
| MetisPartitioner                                                                                                   |
|   partition()                          1         0.0243      0.024285    0.0286      0.028590    11.42    13.45    |
|                                                                                                                    |
| Parallel                                                                                                           |
|   allgather()                          11        0.0023      0.000212    0.0026      0.000234    1.10     1.21     |
|   max(bool)                            1         0.0000      0.000045    0.0000      0.000045    0.02     0.02     |
|   max(scalar)                          113       0.0020      0.000018    0.0020      0.000018    0.94     0.94     |
|   max(vector)                          26        0.0006      0.000024    0.0019      0.000075    0.29     0.91     |
|   min(bool)                            131       0.0022      0.000017    0.0022      0.000017    1.03     1.03     |
|   min(scalar)                          107       0.0179      0.000167    0.0179      0.000167    8.40     8.40     |
|   min(vector)                          26        0.0007      0.000027    0.0021      0.000082    0.33     1.00     |
|   probe()                              132       0.0012      0.000009    0.0012      0.000009    0.58     0.58     |
|   receive()                            132       0.0008      0.000006    0.0020      0.000015    0.36     0.96     |
|   send()                               132       0.0004      0.000003    0.0004      0.000003    0.19     0.19     |
|   send_receive()                       136       0.0011      0.000008    0.0039      0.000029    0.52     1.83     |
|   sum()                                20        0.0009      0.000045    0.0016      0.000078    0.42     0.73     |
|                                                                                                                    |
| Parallel::Request                                                                                                  |
|   wait()                               132       0.0003      0.000002    0.0003      0.000002    0.13     0.13     |
|                                                                                                                    |
| Partitioner                                                                                                        |
|   set_node_processor_ids()             1         0.0013      0.001275    0.0020      0.001981    0.60     0.93     |
|   set_parent_processor_ids()           1         0.0004      0.000443    0.0004      0.000443    0.21     0.21     |
|                                                                                                                    |
| PetscLinearSolver                                                                                                  |
|   solve()                              1         0.0326      0.032632    0.0326      0.032632    15.35    15.35    |
|                                                                                                                    |
| System                                                                                                             |
|   assemble()                           1         0.0030      0.003034    0.0154      0.015380    1.43     7.23     |
 --------------------------------------------------------------------------------------------------------------------
| Totals:                                2044      0.2126                                          100.00            |
 --------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example introduction_ex4:
*  mpirun -np 12 example-devel -d 2 -n 15 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
***************************************************************
* Running Example introduction_ex4:
*  mpirun -np 12 example-devel -d 3 -n 6 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/introduction/introduction_ex4/.libs/lt-example-devel -d 3 -n 6 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=2197
    n_local_nodes()=245
  n_elem()=216
    n_local_elem()=18
    n_active_elem()=216
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0

 EquationSystems
  n_systems()=1
   System #0, "Poisson"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" 
    n_dofs()=2197
    n_local_dofs()=245
    n_constrained_dofs()=866
    n_local_constrained_dofs()=101
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 41.4128
      Average Off-Processor Bandwidth <= 16.6208
      Maximum  On-Processor Bandwidth <= 125
      Maximum Off-Processor Bandwidth <= 126
    DofMap Constraints
      Number of DoF Constraints = 866
      Number of Heterogenous Constraints= 818
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=2197
    n_local_nodes()=245
  n_elem()=216
    n_local_elem()=18
    n_active_elem()=216
  n_subdomains()=1
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 21:56:48 2013                                                                             |
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
| Matrix Assembly Performance: Alive time=0.044919, Active time=0.033564                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| Fe                            18        0.0007      0.000036    0.0007      0.000036    1.96     1.96     |
| Ke                            18        0.0188      0.001045    0.0188      0.001045    56.03    56.03    |
| elem init                     18        0.0136      0.000754    0.0136      0.000754    40.41    40.41    |
| matrix insertion              18        0.0005      0.000030    0.0005      0.000030    1.60     1.60     |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       72        0.0336                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/introduction/introduction_ex4/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 21:56:48 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           1.262e+00      1.00000   1.262e+00
Objects:              6.200e+01      1.00000   6.200e+01
Flops:                6.179e+06      3.98700   3.925e+06  4.710e+07
Flops/sec:            4.896e+06      3.98700   3.109e+06  3.731e+07
MPI Messages:         4.710e+02      2.63866   2.978e+02  3.574e+03
MPI Message Lengths:  1.653e+05      1.83718   3.878e+02  1.386e+06
MPI Reductions:       1.160e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.2621e+00 100.0%  4.7095e+07 100.0%  3.574e+03 100.0%  3.878e+02      100.0%  1.150e+02  99.1% 

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

VecMDot               34 1.0 2.0378e-03 6.3 2.32e+05 2.3 0.0e+00 0.0e+00 3.4e+01  0  4  0  0 29   0  4  0  0 30  1021
VecNorm               37 1.0 7.3378e-0330.7 1.81e+04 2.3 0.0e+00 0.0e+00 3.7e+01  0  0  0  0 32   0  0  0  0 32    22
VecScale              36 1.0 4.7684e-05 1.2 8.82e+03 2.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1659
VecCopy                3 1.0 7.8678e-06 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                44 1.0 3.0279e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                4 1.0 2.9802e-05 2.0 1.96e+03 2.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   590
VecMAXPY              36 1.0 1.4710e-04 2.0 2.49e+05 2.3 0.0e+00 0.0e+00 0.0e+00  0  5  0  0  0   0  5  0  0  0 15204
VecAssemblyBegin       3 1.0 2.1195e-04 1.1 0.00e+00 0.0 6.8e+01 3.1e+02 9.0e+00  0  0  2  2  8   0  0  2  2  8     0
VecAssemblyEnd         3 1.0 3.4809e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       37 1.0 2.7871e-04 1.6 0.00e+00 0.0 3.0e+03 2.7e+02 0.0e+00  0  0 85 58  0   0  0 85 58  0     0
VecScatterEnd         37 1.0 1.2301e-0283.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize          36 1.0 7.4790e-0319.5 2.65e+04 2.3 0.0e+00 0.0e+00 3.6e+01  0  1  0  0 31   0  1  0  0 31    32
MatMult               36 1.0 1.2938e-0213.5 9.88e+05 2.6 3.0e+03 2.6e+02 0.0e+00  0 18 83 55  0   0 18 83 55  0   649
MatSolve              37 1.0 1.9011e-03 4.1 3.03e+06 4.4 0.0e+00 0.0e+00 0.0e+00  0 49  0  0  0   0 49  0  0  0 12235
MatLUFactorNum         1 1.0 2.0862e-03 5.1 1.65e+06 6.4 0.0e+00 0.0e+00 0.0e+00  0 23  0  0  0   0 23  0  0  0  5209
MatILUFactorSym        1 1.0 1.4885e-02 7.4 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  1  0  0  0  3   1  0  0  0  3     0
MatAssemblyBegin       2 1.0 6.9617e-0254.4 0.00e+00 0.0 1.0e+02 5.1e+03 4.0e+00  3  0  3 38  3   3  0  3 38  3     0
MatAssemblyEnd         2 1.0 1.1280e-03 1.1 0.00e+00 0.0 1.6e+02 6.7e+01 8.0e+00  0  0  5  1  7   0  0  5  1  7     0
MatGetRowIJ            1 1.0 7.8678e-06 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 9.5844e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries         3 1.0 8.4162e-05 2.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        34 1.0 2.1822e-03 4.1 4.65e+05 2.3 0.0e+00 0.0e+00 3.4e+01  0  9  0  0 29   0  9  0  0 30  1910
KSPSetUp               2 1.0 1.2517e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 2.2921e-02 1.0 6.18e+06 4.0 3.0e+03 2.6e+02 7.8e+01  2100 83 55 67   2100 83 55 68  2055
PCSetUp                2 1.0 1.7738e-02 5.6 1.65e+06 6.4 0.0e+00 0.0e+00 7.0e+00  1 23  0  0  6   1 23  0  0  6   613
PCSetUpOnBlocks        1 1.0 1.7235e-02 6.4 1.65e+06 6.4 0.0e+00 0.0e+00 5.0e+00  1 23  0  0  4   1 23  0  0  4   630
PCApply               37 1.0 2.3384e-03 2.6 3.03e+06 4.4 0.0e+00 0.0e+00 0.0e+00  0 49  0  0  0   0 49  0  0  0  9947
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    43             43       143704     0
      Vector Scatter     2              2         2072     0
           Index Set     7              7         8748     0
   IS L to G Mapping     1              1          564     0
              Matrix     4              4       698472     0
       Krylov Solver     2              2        19360     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 3.19481e-06
Average time for zero size MPI_Send(): 1.36693e-05
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

 --------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=1.30335, Active time=1.23348                                                       |
 --------------------------------------------------------------------------------------------------------------------
| Event                                  nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                  w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------------|
|                                                                                                                    |
|                                                                                                                    |
| DofMap                                                                                                             |
|   add_neighbors_to_send_list()         1         0.0066      0.006555    0.0242      0.024224    0.53     1.96     |
|   build_constraint_matrix_and_vector() 18        0.0009      0.000049    0.0009      0.000049    0.07     0.07     |
|   build_sparsity()                     1         0.0082      0.008152    0.0197      0.019701    0.66     1.60     |
|   create_dof_constraints()             1         0.0746      0.074622    0.7905      0.790477    6.05     64.09    |
|   distribute_dofs()                    1         0.0254      0.025360    0.0703      0.070316    2.06     5.70     |
|   dof_indices()                        321       0.1051      0.000327    0.1051      0.000327    8.52     8.52     |
|   hetero_cnstrn_elem_mat_vec()         18        0.0096      0.000532    0.0096      0.000532    0.78     0.78     |
|   prepare_send_list()                  1         0.0004      0.000384    0.0004      0.000384    0.03     0.03     |
|   reinit()                             1         0.0433      0.043323    0.0433      0.043323    3.51     3.51     |
|                                                                                                                    |
| EquationSystems                                                                                                    |
|   build_solution_vector()              1         0.0007      0.000716    0.0075      0.007540    0.06     0.61     |
|                                                                                                                    |
| ExodusII_IO                                                                                                        |
|   write_nodal_data()                   1         0.0069      0.006946    0.0069      0.006946    0.56     0.56     |
|                                                                                                                    |
| FE                                                                                                                 |
|   compute_shape_functions()            1026      0.0597      0.000058    0.0597      0.000058    4.84     4.84     |
|   init_shape_functions()               1009      0.0081      0.000008    0.0081      0.000008    0.66     0.66     |
|   inverse_map()                        2376      0.1597      0.000067    0.1597      0.000067    12.95    12.95    |
|                                                                                                                    |
| FEMap                                                                                                              |
|   compute_affine_map()                 1026      0.0318      0.000031    0.0318      0.000031    2.58     2.58     |
|   compute_edge_map()                   792       0.0080      0.000010    0.0080      0.000010    0.65     0.65     |
|   compute_face_map()                   216       0.0112      0.000052    0.0112      0.000052    0.90     0.90     |
|   init_edge_shape_functions()          792       0.0048      0.000006    0.0048      0.000006    0.39     0.39     |
|   init_face_shape_functions()          216       0.0237      0.000110    0.0237      0.000110    1.92     1.92     |
|   init_reference_to_physical_map()     1009      0.3449      0.000342    0.3449      0.000342    27.96    27.96    |
|                                                                                                                    |
| Mesh                                                                                                               |
|   find_neighbors()                     1         0.0075      0.007549    0.0077      0.007677    0.61     0.62     |
|   renumber_nodes_and_elem()            2         0.0013      0.000651    0.0013      0.000651    0.11     0.11     |
|                                                                                                                    |
| MeshCommunication                                                                                                  |
|   compute_hilbert_indices()            2         0.0039      0.001940    0.0039      0.001940    0.31     0.31     |
|   find_global_indices()                2         0.0019      0.000949    0.0101      0.005049    0.15     0.82     |
|   parallel_sort()                      2         0.0026      0.001321    0.0031      0.001538    0.21     0.25     |
|                                                                                                                    |
| MeshOutput                                                                                                         |
|   write_equation_systems()             1         0.0002      0.000182    0.0148      0.014805    0.01     1.20     |
|                                                                                                                    |
| MeshTools::Generation                                                                                              |
|   build_cube()                         1         0.0065      0.006518    0.0065      0.006518    0.53     0.53     |
|                                                                                                                    |
| MetisPartitioner                                                                                                   |
|   partition()                          1         0.0379      0.037882    0.0422      0.042176    3.07     3.42     |
|                                                                                                                    |
| Parallel                                                                                                           |
|   allgather()                          11        0.0004      0.000033    0.0004      0.000039    0.03     0.03     |
|   max(bool)                            1         0.0000      0.000008    0.0000      0.000008    0.00     0.00     |
|   max(scalar)                          113       0.0008      0.000007    0.0008      0.000007    0.06     0.06     |
|   max(vector)                          26        0.0003      0.000013    0.0008      0.000032    0.03     0.07     |
|   min(bool)                            131       0.0008      0.000006    0.0008      0.000006    0.06     0.06     |
|   min(scalar)                          107       0.1119      0.001045    0.1119      0.001045    9.07     9.07     |
|   min(vector)                          26        0.0004      0.000017    0.0011      0.000043    0.04     0.09     |
|   probe()                              132       0.0016      0.000012    0.0016      0.000012    0.13     0.13     |
|   receive()                            132       0.0008      0.000006    0.0025      0.000019    0.07     0.20     |
|   send()                               132       0.0004      0.000003    0.0004      0.000003    0.03     0.03     |
|   send_receive()                       136       0.0013      0.000009    0.0045      0.000033    0.10     0.37     |
|   sum()                                20        0.0006      0.000029    0.0010      0.000049    0.05     0.08     |
|                                                                                                                    |
| Parallel::Request                                                                                                  |
|   wait()                               132       0.0003      0.000002    0.0003      0.000002    0.02     0.02     |
|                                                                                                                    |
| Partitioner                                                                                                        |
|   set_node_processor_ids()             1         0.0019      0.001934    0.0027      0.002709    0.16     0.22     |
|   set_parent_processor_ids()           1         0.0004      0.000437    0.0004      0.000437    0.04     0.04     |
|                                                                                                                    |
| PetscLinearSolver                                                                                                  |
|   solve()                              1         0.0947      0.094730    0.0947      0.094730    7.68     7.68     |
|                                                                                                                    |
| System                                                                                                             |
|   assemble()                           1         0.0214      0.021423    0.0452      0.045172    1.74     3.66     |
 --------------------------------------------------------------------------------------------------------------------
| Totals:                                9942      1.2335                                          100.00            |
 --------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example introduction_ex4:
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
