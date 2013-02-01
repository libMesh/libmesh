<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("subdomains_ex1",$root)?>
 
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
<br><br><br> <h1> The source file subdomains_ex1.C with comments: </h1> 
<div class = "comment">
<h1>Subdomains Example 1 - Solving on a Subdomain</h1>

<br><br>This example builds on the example 4 by showing what to do in
order to solve an equation only on a subdomain.


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
        #include "libmesh/gmv_io.h"
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
        
        #include "libmesh/mesh_refinement.h"
        
</pre>
</div>
<div class = "comment">
Classes needed for subdomain computation.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/system_subset_by_subdomain.h"
        
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
Only our PETSc interface currently supports solves restricted to
subdomains
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, "--enable-petsc");
        
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
          std::string order = "FIRST"; 
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
To demonstate solving on a subdomain, we will solve only on the
interior of a circle (ball in 3d) with radius 0.8.  So show that
this also works well on locally refined meshes, we refine once
all elements that are located on the boundary of this circle (or
ball).
</div>

<div class ="fragment">
<pre>
          {
</pre>
</div>
<div class = "comment">
A MeshRefinement object is needed to refine meshes.
</div>

<div class ="fragment">
<pre>
            MeshRefinement meshRefinement(mesh);
        
</pre>
</div>
<div class = "comment">
Loop over all elements.
</div>

<div class ="fragment">
<pre>
            MeshBase::element_iterator       elem_it  = mesh.elements_begin();
            const MeshBase::element_iterator elem_end = mesh.elements_end(); 
            for (; elem_it != elem_end; ++elem_it)
              {
        	Elem* elem = *elem_it;
        	if(elem-&gt;active())
        	  {
</pre>
</div>
<div class = "comment">
Just check whether the current element has at least one
node inside and one node outside the circle.
</div>

<div class ="fragment">
<pre>
                    bool node_in = false;
        	    bool node_out = false;
        	    for(unsigned int i=0; i&lt;elem-&gt;n_nodes(); i++)
        	      {
        		double d = elem-&gt;point(i).size();
        		if(d&lt;0.8)
        		  {
        		    node_in = true;
        		  }
        		else
        		  {
        		    node_out = true;
        		  }
        	      }
        	    if(node_in && node_out)
        	      {
        		elem-&gt;set_refinement_flag(Elem::REFINE);
        	      }
        	    else
        	      {
        		elem-&gt;set_refinement_flag(Elem::DO_NOTHING);
        	      }
        	  }
        	else
        	  {
        	    elem-&gt;set_refinement_flag(Elem::INACTIVE);
        	  }
              }
        
</pre>
</div>
<div class = "comment">
Now actually refine.
</div>

<div class ="fragment">
<pre>
            meshRefinement.refine_elements();
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
Now set the subdomain_id of all elements whose centroid is inside
the circle to 1.
</div>

<div class ="fragment">
<pre>
          {
</pre>
</div>
<div class = "comment">
Loop over all elements.
</div>

<div class ="fragment">
<pre>
            MeshBase::element_iterator       elem_it  = mesh.elements_begin();
            const MeshBase::element_iterator elem_end = mesh.elements_end(); 
            for (; elem_it != elem_end; ++elem_it)
              {
        	Elem* elem = *elem_it;
        	double d = elem-&gt;centroid().size();
        	if(d&lt;0.8)
        	  {
        	    elem-&gt;subdomain_id() = 1;
        	  }
              }
          }
        
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
Print information about the system to the screen.
</div>

<div class ="fragment">
<pre>
          equation_systems.print_info();
          mesh.print_info();
        
</pre>
</div>
<div class = "comment">
Restrict solves to those elements that have subdomain_id set to 1.
</div>

<div class ="fragment">
<pre>
          std::set&lt;subdomain_id_type&gt; id_list;
          id_list.insert(1);
          SystemSubsetBySubdomain::SubdomainSelectionByList selection(id_list);
          SystemSubsetBySubdomain subset(system,selection);
          system.restrict_solve_to(&subset,SUBSET_ZERO);
        
</pre>
</div>
<div class = "comment">
Note that using \p SUBSET_ZERO will cause all dofs outside the
subdomain to be cleared.  This will, however, cause some hanging
nodes outside the subdomain to have inconsistent values.
  

<br><br>Solve the system "Poisson", just like example 2.
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
            GnuPlotIO plot(mesh,"Subdomains Example 1, 1D",GnuPlotIO::GRID_ON);
            plot.write_equation_systems("gnuplot_script",equation_systems);
          }
          else
          {
            GMVIO (mesh).write_equation_systems ((dim == 3) ? 
              "out_3.gmv" : "out_2.gmv",equation_systems);
        #ifdef LIBMESH_HAVE_EXODUS_API
            ExodusII_IO (mesh).write_equation_systems ((dim == 3) ? 
              "out_3.e" : "out_2.e",equation_systems);
        #endif // #ifdef LIBMESH_HAVE_EXODUS_API
          }
        
        #endif // #ifndef LIBMESH_ENABLE_AMR
          
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
Store a pointer to the element we are currently
working on.  This allows for nicer syntax later.
</div>

<div class ="fragment">
<pre>
              const Elem* elem = *el;
        
</pre>
</div>
<div class = "comment">
Elements with subdomain_id other than 1 are not in the active
subdomain.  We don't assemble anything for them.
</div>

<div class ="fragment">
<pre>
              if(elem-&gt;subdomain_id()==1)
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
The following loops over the sides of the element.  If
the element has no neighbor on a side then that side
MUST live on a boundary of the domain.  If there is a
neighbor, check that neighbor's subdomain_id; if that
is different from 1, the side is also located on the
boundary.
</div>

<div class ="fragment">
<pre>
                    for (unsigned int side=0; side&lt;elem-&gt;n_sides(); side++)
        	      if ((elem-&gt;neighbor(side) == NULL) ||
        		  (elem-&gt;neighbor(side)-&gt;subdomain_id()!=1))
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
If this assembly program were to be used on an adaptive mesh,
we would have to apply any hanging node constraint equations
</div>

<div class ="fragment">
<pre>
                  dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        	  
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
<br><br><br> <h1> The source file subdomains_ex1.C without comments: </h1> 
<pre> 
  
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/gmv_io.h&quot;</FONT></B>
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
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_refinement.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/system_subset_by_subdomain.h&quot;</FONT></B>
  
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
  
    libmesh_example_assert(libMesh::default_solver_package() == PETSC_SOLVERS, <B><FONT COLOR="#BC8F8F">&quot;--enable-petsc&quot;</FONT></B>);
  
  #ifndef LIBMESH_ENABLE_AMR
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-amr&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    
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
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string order = <B><FONT COLOR="#BC8F8F">&quot;FIRST&quot;</FONT></B>; 
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
  
  
    {
      MeshRefinement meshRefinement(mesh);
  
      <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::element_iterator       elem_it  = mesh.elements_begin();
      <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::element_iterator elem_end = mesh.elements_end(); 
      <B><FONT COLOR="#A020F0">for</FONT></B> (; elem_it != elem_end; ++elem_it)
        {
  	Elem* elem = *elem_it;
  	<B><FONT COLOR="#A020F0">if</FONT></B>(elem-&gt;active())
  	  {
  	    <B><FONT COLOR="#228B22">bool</FONT></B> node_in = false;
  	    <B><FONT COLOR="#228B22">bool</FONT></B> node_out = false;
  	    <B><FONT COLOR="#A020F0">for</FONT></B>(<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;elem-&gt;n_nodes(); i++)
  	      {
  		<B><FONT COLOR="#228B22">double</FONT></B> d = elem-&gt;point(i).size();
  		<B><FONT COLOR="#A020F0">if</FONT></B>(d&lt;0.8)
  		  {
  		    node_in = true;
  		  }
  		<B><FONT COLOR="#A020F0">else</FONT></B>
  		  {
  		    node_out = true;
  		  }
  	      }
  	    <B><FONT COLOR="#A020F0">if</FONT></B>(node_in &amp;&amp; node_out)
  	      {
  		elem-&gt;set_refinement_flag(Elem::REFINE);
  	      }
  	    <B><FONT COLOR="#A020F0">else</FONT></B>
  	      {
  		elem-&gt;set_refinement_flag(Elem::DO_NOTHING);
  	      }
  	  }
  	<B><FONT COLOR="#A020F0">else</FONT></B>
  	  {
  	    elem-&gt;set_refinement_flag(Elem::INACTIVE);
  	  }
        }
  
      meshRefinement.refine_elements();
    }
  
    mesh.print_info();
    
    {
      <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::element_iterator       elem_it  = mesh.elements_begin();
      <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::element_iterator elem_end = mesh.elements_end(); 
      <B><FONT COLOR="#A020F0">for</FONT></B> (; elem_it != elem_end; ++elem_it)
        {
  	Elem* elem = *elem_it;
  	<B><FONT COLOR="#228B22">double</FONT></B> d = elem-&gt;centroid().size();
  	<B><FONT COLOR="#A020F0">if</FONT></B>(d&lt;0.8)
  	  {
  	    elem-&gt;subdomain_id() = 1;
  	  }
        }
    }
  
    EquationSystems equation_systems (mesh);
    
    LinearImplicitSystem&amp; system =
      equation_systems.add_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>);
  
    
    system.add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>,
                        <B><FONT COLOR="#5F9EA0">Utility</FONT></B>::string_to_enum&lt;Order&gt;   (order),
                        <B><FONT COLOR="#5F9EA0">Utility</FONT></B>::string_to_enum&lt;FEFamily&gt;(family));
  
    system.attach_assemble_function (assemble_poisson);
    
    equation_systems.init();
  
    equation_systems.print_info();
    mesh.print_info();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;subdomain_id_type&gt; id_list;
    id_list.insert(1);
    <B><FONT COLOR="#5F9EA0">SystemSubsetBySubdomain</FONT></B>::SubdomainSelectionByList selection(id_list);
    SystemSubsetBySubdomain subset(system,selection);
    system.restrict_solve_to(&amp;subset,SUBSET_ZERO);
  
    
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>).solve();
  
    <B><FONT COLOR="#A020F0">if</FONT></B>(dim == 1)
    {        
      GnuPlotIO plot(mesh,<B><FONT COLOR="#BC8F8F">&quot;Subdomains Example 1, 1D&quot;</FONT></B>,GnuPlotIO::GRID_ON);
      plot.write_equation_systems(<B><FONT COLOR="#BC8F8F">&quot;gnuplot_script&quot;</FONT></B>,equation_systems);
    }
    <B><FONT COLOR="#A020F0">else</FONT></B>
    {
      GMVIO (mesh).write_equation_systems ((dim == 3) ? 
        <B><FONT COLOR="#BC8F8F">&quot;out_3.gmv&quot;</FONT></B> : <B><FONT COLOR="#BC8F8F">&quot;out_2.gmv&quot;</FONT></B>,equation_systems);
  #ifdef LIBMESH_HAVE_EXODUS_API
      ExodusII_IO (mesh).write_equation_systems ((dim == 3) ? 
        <B><FONT COLOR="#BC8F8F">&quot;out_3.e&quot;</FONT></B> : <B><FONT COLOR="#BC8F8F">&quot;out_2.e&quot;</FONT></B>,equation_systems);
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>    }
  
  #endif <I><FONT COLOR="#B22222">// #ifndef LIBMESH_ENABLE_AMR
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
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
        <B><FONT COLOR="#A020F0">if</FONT></B>(elem-&gt;subdomain_id()==1)
  	{
  	  perf_log.push(<B><FONT COLOR="#BC8F8F">&quot;elem init&quot;</FONT></B>);      
  	  
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
  	  
  	  {
  	    
  	    perf_log.push (<B><FONT COLOR="#BC8F8F">&quot;BCs&quot;</FONT></B>);
  	    
  	    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> side=0; side&lt;elem-&gt;n_sides(); side++)
  	      <B><FONT COLOR="#A020F0">if</FONT></B> ((elem-&gt;neighbor(side) == NULL) ||
  		  (elem-&gt;neighbor(side)-&gt;subdomain_id()!=1))
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
  	  
  	  dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
  	  
  	  perf_log.push (<B><FONT COLOR="#BC8F8F">&quot;matrix insertion&quot;</FONT></B>);
  	  
  	  system.matrix-&gt;add_matrix (Ke, dof_indices);
  	  system.rhs-&gt;add_vector    (Fe, dof_indices);
  	  
  	  perf_log.pop (<B><FONT COLOR="#BC8F8F">&quot;matrix insertion&quot;</FONT></B>);
  	}
      }
  
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example subdomains_ex1:
*  mpirun -np 12 example-devel -d 1 -n 40 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/subdomains/subdomains_ex1/.libs/lt-example-devel -d 1 -n 40 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=1
  spatial_dimension()=3
  n_nodes()=43
    n_local_nodes()=5
  n_elem()=44
    n_local_elem()=4
    n_active_elem()=42
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
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=43
    n_local_dofs()=5
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 2.69767
      Average Off-Processor Bandwidth <= 0.511628
      Maximum  On-Processor Bandwidth <= 3
      Maximum Off-Processor Bandwidth <= 1
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=1
  spatial_dimension()=3
  n_nodes()=43
    n_local_nodes()=5
  n_elem()=44
    n_local_elem()=4
    n_active_elem()=42
  n_subdomains()=2
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 22:11:59 2013                                                                             |
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
| Matrix Assembly Performance: Alive time=0.001512, Active time=0.000774                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           4         0.0000      0.000001    0.0000      0.000001    0.78     0.78     |
| Fe                            4         0.0001      0.000025    0.0001      0.000025    13.05    13.05    |
| Ke                            4         0.0000      0.000005    0.0000      0.000005    2.33     2.33     |
| elem init                     4         0.0006      0.000140    0.0006      0.000140    72.61    72.61    |
| matrix insertion              4         0.0001      0.000022    0.0001      0.000022    11.24    11.24    |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       20        0.0008                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/subdomains/subdomains_ex1/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:11:59 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           7.887e-02      1.00355   7.864e-02
Objects:              8.000e+01      1.02564   7.817e+01
Flops:                5.088e+03      0.00000   2.892e+03  3.470e+04
Flops/sec:            6.472e+04      0.00000   3.677e+04  4.412e+05
MPI Messages:         6.500e+01     21.66667   4.992e+01  5.990e+02
MPI Message Lengths:  5.040e+02     22.90909   7.593e+00  4.548e+03
MPI Reductions:       1.450e+02      1.01399

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 7.8587e-02  99.9%  3.4705e+04 100.0%  5.990e+02 100.0%  7.593e+00      100.0%  1.422e+02  99.4% 

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

VecMDot               19 1.0 1.8537e-0313.7 1.71e+03 0.0 0.0e+00 0.0e+00 1.9e+01  1 33  0  0 13   1 33  0  0 13     6
VecNorm               21 1.0 9.6173e-0320.1 2.10e+02 0.0 0.0e+00 0.0e+00 2.1e+01  9  4  0  0 15   9  4  0  0 15     0
VecScale              20 1.0 2.4080e-05 1.3 1.00e+02 0.0 0.0e+00 0.0e+00 0.0e+00  0  2  0  0  0   0  2  0  0  0    29
VecCopy                2 1.0 1.1206e-05 5.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                30 1.0 1.6928e-05 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                2 1.0 1.5836e-02 2.2 2.00e+01 0.0 0.0e+00 0.0e+00 0.0e+00 13  0  0  0  0  13  0  0  0  0     0
VecMAXPY              20 1.0 1.2875e-05 4.2 2.09e+03 0.0 0.0e+00 0.0e+00 0.0e+00  0 42  0  0  0   0 42  0  0  0  1136
VecAssemblyBegin       3 1.0 1.1921e-03 5.6 0.00e+00 0.0 1.8e+01 6.0e+00 9.0e+00  1  0  3  2  6   1  0  3  2  6     0
VecAssemblyEnd         3 1.0 3.1233e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       24 1.0 1.1134e-04 2.9 0.00e+00 0.0 3.8e+02 8.2e+00 0.0e+00  0  0 64 69  0   0  0 64 69  0     0
VecScatterEnd         24 1.0 1.1418e-03399.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize          20 1.0 9.5563e-0317.4 3.00e+02 0.0 0.0e+00 0.0e+00 2.0e+01  8  6  0  0 14   8  6  0  0 14     0
MatMult               20 1.0 1.2386e-0325.6 5.00e+02 0.0 3.6e+02 8.0e+00 0.0e+00  1 10 60 63  0   1 10 60 63  0     3
MatSolve              21 0.0 1.7405e-05 0.0 4.41e+02 0.0 0.0e+00 0.0e+00 0.0e+00  0  8  0  0  0   0  8  0  0  0   163
MatLUFactorNum         1 1.0 3.6955e-05 1.8 1.70e+01 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     3
MatILUFactorSym        1 1.0 8.5831e-05 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  0  0  0  0  2   0  0  0  0  2     0
MatAssemblyBegin       6 1.0 1.7691e-04 1.2 0.00e+00 0.0 2.7e+01 1.2e+01 8.0e+00  0  0  5  7  6   0  0  5  7  6     0
MatAssemblyEnd         6 1.0 6.6662e-04 1.1 0.00e+00 0.0 1.1e+02 4.0e+00 2.4e+01  1  0 18  9 17   1  0 18  9 17     0
MatGetRowIJ            1 0.0 1.2159e-05 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 0.0 1.2493e-04 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.8e+00  0  0  0  0  1   0  0  0  0  1     0
MatZeroEntries         3 1.0 1.4782e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        19 1.0 1.8992e-03 9.9 3.61e+03 0.0 0.0e+00 0.0e+00 1.9e+01  1 71  0  0 13   1 71  0  0 13    13
KSPSetUp               2 1.0 1.1802e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 2.0584e-02 1.0 5.09e+03 0.0 3.6e+02 8.0e+00 4.7e+01 25100 60 63 33  25100 60 63 33     2
PCSetUp                2 1.0 8.1182e-04 1.0 1.70e+01 0.0 0.0e+00 0.0e+00 7.2e+00  1  0  0  0  5   1  0  0  0  5     0
PCSetUpOnBlocks        1 1.0 3.4595e-04 1.1 1.70e+01 0.0 0.0e+00 0.0e+00 5.2e+00  0  0  0  0  4   0  0  0  0  4     0
PCApply               21 1.0 2.9349e-04 1.6 4.41e+02 0.0 0.0e+00 0.0e+00 0.0e+00  0  8  0  0  0   0  8  0  0  0    10
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    40             40        60760     0
      Vector Scatter     5              5         5180     0
           Index Set    15             15        11172     0
   IS L to G Mapping     1              1          564     0
              Matrix    12             12        32148     0
       Krylov Solver     2              2        19360     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 8.39233e-06
Average time for zero size MPI_Send(): 1.3411e-05
#PETSc Option Table entries:
-d 1
-ksp_right_pc
-log_summary
-n 40
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
| libMesh Performance: Alive time=0.198672, Active time=0.063772                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0003      0.000292    0.0004      0.000426    0.46     0.67     |
|   build_sparsity()                 1         0.0010      0.000967    0.0024      0.002363    1.52     3.71     |
|   create_dof_constraints()         1         0.0000      0.000029    0.0000      0.000029    0.05     0.05     |
|   distribute_dofs()                1         0.0012      0.001204    0.0039      0.003886    1.89     6.09     |
|   dof_indices()                    20        0.0006      0.000032    0.0006      0.000032    0.99     0.99     |
|   prepare_send_list()              1         0.0000      0.000011    0.0000      0.000011    0.02     0.02     |
|   reinit()                         1         0.0015      0.001522    0.0015      0.001522    2.39     2.39     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0003      0.000299    0.0036      0.003640    0.47     5.71     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        4         0.0001      0.000013    0.0001      0.000013    0.08     0.08     |
|   init_shape_functions()           1         0.0002      0.000174    0.0002      0.000174    0.27     0.27     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             4         0.0001      0.000014    0.0001      0.000014    0.09     0.09     |
|   init_reference_to_physical_map() 1         0.0001      0.000079    0.0001      0.000079    0.12     0.12     |
|                                                                                                                |
| GnuPlotIO                                                                                                      |
|   write_nodal_data()               1         0.0013      0.001255    0.0013      0.001255    1.97     1.97     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           4         0.0000      0.000010    0.0000      0.000010    0.07     0.07     |
|   init()                           1         0.0002      0.000165    0.0002      0.000165    0.26     0.26     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 2         0.0019      0.000958    0.0022      0.001090    3.01     3.42     |
|   renumber_nodes_and_elem()        4         0.0002      0.000043    0.0002      0.000043    0.27     0.27     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        3         0.0012      0.000390    0.0012      0.000390    1.83     1.83     |
|   find_global_indices()            3         0.0012      0.000410    0.0075      0.002496    1.93     11.74    |
|   parallel_sort()                  3         0.0028      0.000945    0.0034      0.001123    4.45     5.28     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0001      0.000138    0.0057      0.005682    0.22     8.91     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _refine_elements()               1         0.0005      0.000456    0.0006      0.000555    0.72     0.87     |
|   add_point()                      4         0.0000      0.000010    0.0001      0.000021    0.06     0.13     |
|   make_refinement_compatible()     2         0.0001      0.000035    0.0001      0.000044    0.11     0.14     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0006      0.000622    0.0006      0.000622    0.98     0.98     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      2         0.0068      0.003416    0.0109      0.005451    10.71    17.09    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      13        0.0004      0.000030    0.0005      0.000036    0.60     0.73     |
|   max(bool)                        5         0.0001      0.000011    0.0001      0.000011    0.08     0.08     |
|   max(scalar)                      170       0.0022      0.000013    0.0022      0.000013    3.38     3.38     |
|   max(vector)                      40        0.0008      0.000019    0.0022      0.000055    1.21     3.47     |
|   min(bool)                        202       0.0025      0.000012    0.0025      0.000012    3.87     3.87     |
|   min(scalar)                      163       0.0041      0.000025    0.0041      0.000025    6.41     6.41     |
|   min(vector)                      40        0.0008      0.000021    0.0023      0.000058    1.32     3.66     |
|   probe()                          187       0.0009      0.000005    0.0009      0.000005    1.44     1.44     |
|   receive()                        187       0.0010      0.000006    0.0020      0.000011    1.62     3.12     |
|   send()                           187       0.0006      0.000003    0.0006      0.000003    0.87     0.87     |
|   send_receive()                   182       0.0014      0.000008    0.0043      0.000023    2.16     6.69     |
|   sum()                            23        0.0006      0.000027    0.0018      0.000079    0.97     2.85     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           187       0.0004      0.000002    0.0004      0.000002    0.57     0.57     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         2         0.0008      0.000418    0.0020      0.001001    1.31     3.14     |
|   set_parent_processor_ids()       2         0.0002      0.000119    0.0002      0.000119    0.37     0.37     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.0234      0.023428    0.0234      0.023428    36.74    36.74    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0014      0.001385    0.0019      0.001874    2.17     2.94     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            1661      0.0638                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example subdomains_ex1:
*  mpirun -np 12 example-devel -d 1 -n 40 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
***************************************************************
* Running Example subdomains_ex1:
*  mpirun -np 12 example-devel -d 2 -n 30 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/subdomains/subdomains_ex1/.libs/lt-example-devel -d 2 -n 30 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1329
    n_local_nodes()=127
  n_elem()=1268
    n_local_elem()=107
    n_active_elem()=1176
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
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=1329
    n_local_dofs()=127
    n_constrained_dofs()=184
    n_local_constrained_dofs()=18
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.66667
      Average Off-Processor Bandwidth <= 1.02784
      Maximum  On-Processor Bandwidth <= 17
      Maximum Off-Processor Bandwidth <= 10
    DofMap Constraints
      Number of DoF Constraints = 184
      Average DoF Constraint Length= 2
      Number of Node Constraints = 368
      Maximum Node Constraint Length= 3
      Average Node Constraint Length= 2.5

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1329
    n_local_nodes()=127
  n_elem()=1268
    n_local_elem()=107
    n_active_elem()=1176
  n_subdomains()=2
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 22:12:00 2013                                                                             |
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
| Matrix Assembly Performance: Alive time=0.024301, Active time=0.012752                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           66        0.0030      0.000046    0.0030      0.000046    23.66    23.66    |
| Fe                            66        0.0007      0.000011    0.0007      0.000011    5.76     5.76     |
| Ke                            66        0.0010      0.000015    0.0010      0.000015    7.79     7.79     |
| elem init                     66        0.0076      0.000115    0.0076      0.000115    59.49    59.49    |
| matrix insertion              66        0.0004      0.000006    0.0004      0.000006    3.30     3.30     |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       330       0.0128                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/subdomains/subdomains_ex1/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:12:00 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           4.946e-01      1.00000   4.946e-01
Objects:              8.700e+01      1.00000   8.700e+01
Flops:                6.098e+05      4.29779   3.768e+05  4.522e+06
Flops/sec:            1.233e+06      4.29780   7.619e+05  9.142e+06
MPI Messages:         3.590e+02      2.96694   2.112e+02  2.534e+03
MPI Message Lengths:  2.500e+04      3.72638   7.016e+01  1.778e+05
MPI Reductions:       1.980e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 4.9455e-01 100.0%  4.5219e+06 100.0%  2.534e+03 100.0%  7.016e+01      100.0%  1.970e+02  99.5% 

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

VecMDot               46 1.0 6.6280e-04 1.7 1.21e+05 3.7 0.0e+00 0.0e+00 4.6e+01  0 20  0  0 23   0 20  0  0 23  1384
VecNorm               49 1.0 3.8433e-04 1.4 9.90e+03 3.6 0.0e+00 0.0e+00 4.9e+01  0  2  0  0 25   0  2  0  0 25   196
VecScale              48 1.0 4.8161e-05 1.4 4.85e+03 3.6 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0   766
VecCopy                3 1.0 9.0599e-06 4.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                59 1.0 3.0756e-05 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                4 1.0 2.8133e-05 2.0 8.08e+02 3.6 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   219
VecMAXPY              48 1.0 8.6784e-05 2.3 1.31e+05 3.6 0.0e+00 0.0e+00 0.0e+00  0 22  0  0  0   0 22  0  0  0 11466
VecAssemblyBegin       3 1.0 1.7309e-04 1.1 0.00e+00 0.0 4.4e+01 9.0e+01 9.0e+00  0  0  2  2  5   0  0  2  2  5     0
VecAssemblyEnd         3 1.0 3.0041e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       52 1.0 2.3961e-04 1.7 0.00e+00 0.0 2.1e+03 6.8e+01 0.0e+00  0  0 81 78  0   0  0 81 78  0     0
VecScatterEnd         52 1.0 3.0661e-04 3.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize          48 1.0 5.0068e-04 1.2 1.45e+04 3.6 0.0e+00 0.0e+00 4.8e+01  0  2  0  0 24   0  2  0  0 24   221
MatMult               48 1.0 6.5494e-04 1.8 8.34e+04 4.0 2.0e+03 6.6e+01 0.0e+00  0 14 80 75  0   0 14 80 75  0   932
MatSolve              49 1.0 2.5988e-04 4.7 2.56e+05 5.8 0.0e+00 0.0e+00 0.0e+00  0 37  0  0  0   0 37  0  0  0  6400
MatLUFactorNum         1 1.0 1.0300e-04 2.7 4.29e+0412.8 0.0e+00 0.0e+00 0.0e+00  0  5  0  0  0   0  5  0  0  0  2111
MatILUFactorSym        1 1.0 2.8110e-04 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  0  0  0  0  2   0  0  0  0  2     0
MatAssemblyBegin       6 1.0 2.1815e-04 1.3 0.00e+00 0.0 6.6e+01 3.1e+02 8.0e+00  0  0  3 11  4   0  0  3 11  4     0
MatAssemblyEnd         6 1.0 7.9298e-04 1.0 0.00e+00 0.0 2.5e+02 1.9e+01 2.4e+01  0  0 10  3 12   0  0 10  3 12     0
MatGetRowIJ            1 1.0 1.5020e-05 7.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 9.3937e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  1   0  0  0  0  1     0
MatZeroEntries         3 1.0 1.9789e-05 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        46 1.0 7.7319e-04 1.4 2.42e+05 3.6 0.0e+00 0.0e+00 4.6e+01  0 41  0  0 23   0 41  0  0 23  2382
KSPSetUp               2 1.0 1.1897e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 3.9020e-03 1.0 6.10e+05 4.3 2.0e+03 6.6e+01 1.0e+02  1100 80 75 52   1100 80 75 52  1159
PCSetUp                2 1.0 1.0920e-03 1.2 4.29e+0412.8 0.0e+00 0.0e+00 7.0e+00  0  5  0  0  4   0  5  0  0  4   199
PCSetUpOnBlocks        1 1.0 6.2990e-04 1.5 4.29e+0412.8 0.0e+00 0.0e+00 5.0e+00  0  5  0  0  3   0  5  0  0  3   345
PCApply               49 1.0 7.7128e-04 1.3 2.56e+05 5.8 0.0e+00 0.0e+00 0.0e+00  0 37  0  0  0   0 37  0  0  0  2157
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    49             49       104376     0
      Vector Scatter     5              5         5180     0
           Index Set    15             15        12076     0
   IS L to G Mapping     1              1          564     0
              Matrix    12             12       124832     0
       Krylov Solver     2              2        19360     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 3.00407e-06
Average time for zero size MPI_Send(): 1.33514e-05
#PETSc Option Table entries:
-d 2
-ksp_right_pc
-log_summary
-n 30
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
| libMesh Performance: Alive time=0.586205, Active time=0.459663                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0059      0.005886    0.0093      0.009287    1.28     2.02     |
|   build_constraint_matrix()        66        0.0011      0.000016    0.0011      0.000016    0.23     0.23     |
|   build_sparsity()                 1         0.0047      0.004743    0.0119      0.011916    1.03     2.59     |
|   cnstrn_elem_mat_vec()            66        0.0086      0.000131    0.0086      0.000131    1.88     1.88     |
|   create_dof_constraints()         1         0.0188      0.018762    0.0448      0.044777    4.08     9.74     |
|   distribute_dofs()                1         0.0212      0.021206    0.0678      0.067756    4.61     14.74    |
|   dof_indices()                    1225      0.0481      0.000039    0.0481      0.000039    10.47    10.47    |
|   prepare_send_list()              1         0.0001      0.000055    0.0001      0.000055    0.01     0.01     |
|   reinit()                         1         0.0445      0.044514    0.0445      0.044514    9.68     9.68     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          2         0.0031      0.001534    0.0148      0.007413    0.67     3.23     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0111      0.011140    0.0111      0.011140    2.42     2.42     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        86        0.0020      0.000024    0.0020      0.000024    0.44     0.44     |
|   init_shape_functions()           21        0.0003      0.000012    0.0003      0.000012    0.06     0.06     |
|   inverse_map()                    980       0.0045      0.000005    0.0045      0.000005    0.98     0.98     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             86        0.0012      0.000014    0.0012      0.000014    0.26     0.26     |
|   compute_face_map()               20        0.0005      0.000027    0.0011      0.000057    0.12     0.25     |
|   init_face_shape_functions()      1         0.0000      0.000021    0.0000      0.000021    0.00     0.00     |
|   init_reference_to_physical_map() 21        0.0004      0.000021    0.0004      0.000021    0.10     0.10     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               1         0.0104      0.010373    0.0105      0.010505    2.26     2.29     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           1104      0.0048      0.000004    0.0048      0.000004    1.04     1.04     |
|   init()                           1         0.0017      0.001691    0.0017      0.001691    0.37     0.37     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 2         0.0459      0.022926    0.0465      0.023230    9.98     10.11    |
|   renumber_nodes_and_elem()        4         0.0017      0.000426    0.0017      0.000426    0.37     0.37     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        3         0.0249      0.008304    0.0249      0.008304    5.42     5.42     |
|   find_global_indices()            3         0.0097      0.003229    0.0433      0.014426    2.11     9.41     |
|   parallel_sort()                  3         0.0039      0.001299    0.0067      0.002221    0.85     1.45     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         2         0.0002      0.000119    0.0369      0.018467    0.05     8.04     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _refine_elements()               1         0.0078      0.007776    0.0180      0.017980    1.69     3.91     |
|   add_point()                      1104      0.0047      0.000004    0.0097      0.000009    1.02     2.10     |
|   make_refinement_compatible()     2         0.0008      0.000397    0.0009      0.000434    0.17     0.19     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0044      0.004360    0.0044      0.004360    0.95     0.95     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      2         0.1166      0.058306    0.1463      0.073161    25.37    31.83    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      13        0.0010      0.000078    0.0012      0.000090    0.22     0.26     |
|   max(bool)                        5         0.0005      0.000096    0.0005      0.000096    0.10     0.10     |
|   max(scalar)                      189       0.0012      0.000007    0.0012      0.000007    0.27     0.27     |
|   max(vector)                      44        0.0006      0.000013    0.0014      0.000031    0.12     0.30     |
|   min(bool)                        225       0.0014      0.000006    0.0014      0.000006    0.29     0.29     |
|   min(scalar)                      182       0.0083      0.000046    0.0083      0.000046    1.82     1.82     |
|   min(vector)                      44        0.0007      0.000015    0.0015      0.000035    0.14     0.34     |
|   probe()                          187       0.0036      0.000019    0.0036      0.000019    0.79     0.79     |
|   receive()                        187       0.0011      0.000006    0.0048      0.000026    0.25     1.04     |
|   send()                           187       0.0006      0.000003    0.0006      0.000003    0.12     0.12     |
|   send_receive()                   182       0.0015      0.000008    0.0055      0.000030    0.32     1.20     |
|   sum()                            26        0.0026      0.000099    0.0037      0.000141    0.56     0.80     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           187       0.0003      0.000002    0.0003      0.000002    0.07     0.07     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         2         0.0050      0.002493    0.0080      0.003995    1.08     1.74     |
|   set_parent_processor_ids()       2         0.0038      0.001882    0.0038      0.001882    0.82     0.82     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.0079      0.007929    0.0079      0.007929    1.72     1.72     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0061      0.006080    0.0247      0.024666    1.32     5.37     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            6478      0.4597                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example subdomains_ex1:
*  mpirun -np 12 example-devel -d 2 -n 30 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
***************************************************************
* Running Example subdomains_ex1:
*  mpirun -np 12 example-devel -d 3 -n 15 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/subdomains/subdomains_ex1/.libs/lt-example-devel -d 3 -n 15 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=10854
    n_local_nodes()=1083
  n_elem()=8767
    n_local_elem()=738
    n_active_elem()=8093
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
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=10854
    n_local_dofs()=1083
    n_constrained_dofs()=4068
    n_local_constrained_dofs()=366
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 27.0504
      Average Off-Processor Bandwidth <= 4.52184
      Maximum  On-Processor Bandwidth <= 88
      Maximum Off-Processor Bandwidth <= 71
    DofMap Constraints
      Number of DoF Constraints = 4068
      Average DoF Constraint Length= 2.66667
      Number of Node Constraints = 5428
      Maximum Node Constraint Length= 13
      Average Node Constraint Length= 6.24613

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=10854
    n_local_nodes()=1083
  n_elem()=8767
    n_local_elem()=738
    n_active_elem()=8093
  n_subdomains()=2
  n_partitions()=12
  n_processors()=12
  n_threads()=1
  processor_id()=0


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 22:12:05 2013                                                                             |
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
| Matrix Assembly Performance: Alive time=0.214361, Active time=0.181901                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           228       0.0606      0.000266    0.0606      0.000266    33.32    33.32    |
| Fe                            228       0.0077      0.000034    0.0077      0.000034    4.24     4.24     |
| Ke                            228       0.0340      0.000149    0.0340      0.000149    18.68    18.68    |
| elem init                     228       0.0768      0.000337    0.0768      0.000337    42.20    42.20    |
| matrix insertion              228       0.0028      0.000012    0.0028      0.000012    1.56     1.56     |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       1140      0.1819                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/subdomains/subdomains_ex1/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:12:06 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           5.495e+00      1.00030   5.494e+00
Objects:              8.700e+01      1.00000   8.700e+01
Flops:                5.713e+07      5.65600   2.833e+07  3.399e+08
Flops/sec:            1.040e+07      5.65601   5.155e+06  6.186e+07
MPI Messages:         4.610e+02      2.74405   2.672e+02  3.206e+03
MPI Message Lengths:  2.271e+05      1.96708   5.840e+02  1.872e+06
MPI Reductions:       1.630e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 5.4944e+00 100.0%  3.3991e+08 100.0%  3.206e+03 100.0%  5.840e+02      100.0%  1.620e+02  99.4% 

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

VecMDot               29 1.0 1.0871e-0215.3 5.02e+05 2.0 0.0e+00 0.0e+00 2.9e+01  0  1  0  0 18   0  1  0  0 18   409
VecNorm               31 1.0 8.6665e-04 2.6 3.58e+04 2.0 0.0e+00 0.0e+00 3.1e+01  0  0  0  0 19   0  0  0  0 19   366
VecScale              30 1.0 7.4387e-05 1.1 1.73e+04 2.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2063
VecCopy                2 1.0 1.0967e-05 3.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                40 1.0 7.1049e-05 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                2 1.0 3.6955e-05 1.8 2.31e+03 2.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   554
VecMAXPY              30 1.0 3.0088e-04 1.5 5.36e+05 2.0 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0 15779
VecAssemblyBegin       3 1.0 5.2500e-04 2.7 0.00e+00 0.0 8.6e+01 9.1e+02 9.0e+00  0  0  3  4  6   0  0  3  4  6     0
VecAssemblyEnd         3 1.0 5.6744e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       34 1.0 4.9067e-04 1.4 0.00e+00 0.0 2.4e+03 4.0e+02 0.0e+00  0  0 74 50  0   0  0 74 50  0     0
VecScatterEnd         34 1.0 1.0377e-01360.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
VecNormalize          30 1.0 7.8392e-04 1.5 5.20e+04 2.0 0.0e+00 0.0e+00 3.0e+01  0  0  0  0 18   0  0  0  0 19   587
MatMult               30 1.0 1.0497e-0144.7 9.10e+05 2.3 2.3e+03 3.8e+02 0.0e+00  1  2 71 47  0   1  2 71 47  0    74
MatSolve              31 1.0 1.6212e-02 2.7 1.33e+07 3.5 0.0e+00 0.0e+00 0.0e+00  0 29  0  0  0   0 29  0  0  0  6092
MatLUFactorNum         1 1.0 4.3848e-02 7.1 4.19e+07 7.9 0.0e+00 0.0e+00 0.0e+00  0 66  0  0  0   0 66  0  0  0  5102
MatILUFactorSym        1 1.0 8.4404e-02 4.6 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  1  0  0  0  2   1  0  0  0  2     0
MatAssemblyBegin       6 1.0 9.7179e-04 3.3 0.00e+00 0.0 1.3e+02 5.7e+03 8.0e+00  0  0  4 39  5   0  0  4 39  5     0
MatAssemblyEnd         6 1.0 2.0180e-03 1.1 0.00e+00 0.0 4.6e+02 9.8e+01 2.4e+01  0  0 14  2 15   0  0 14  2 15     0
MatGetRowIJ            1 1.0 2.1935e-05 3.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 1.1206e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  1   0  0  0  0  1     0
MatZeroEntries         3 1.0 1.8001e-04 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        29 1.0 1.1257e-0210.3 1.01e+06 2.0 0.0e+00 0.0e+00 2.9e+01  0  3  0  0 18   0  3  0  0 18   790
KSPSetUp               2 1.0 1.6713e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 1.5095e-01 1.0 5.71e+07 5.7 2.3e+03 3.8e+02 6.7e+01  3100 71 47 41   3100 71 47 41  2252
PCSetUp                2 1.0 1.2912e-01 5.1 4.19e+07 7.9 0.0e+00 0.0e+00 7.0e+00  1 66  0  0  4   1 66  0  0  4  1732
PCSetUpOnBlocks        1 1.0 1.2854e-01 5.2 4.19e+07 7.9 0.0e+00 0.0e+00 5.0e+00  1 66  0  0  3   1 66  0  0  3  1740
PCApply               31 1.0 1.6892e-02 2.4 1.33e+07 3.5 0.0e+00 0.0e+00 0.0e+00  0 29  0  0  0   0 29  0  0  0  5847
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    49             49       258256     0
      Vector Scatter     5              5         5180     0
           Index Set    15             15        19512     0
   IS L to G Mapping     1              1          564     0
              Matrix    12             12      3382548     0
       Krylov Solver     2              2        19360     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.43459e-06
Average time for zero size MPI_Send(): 1.37488e-05
#PETSc Option Table entries:
-d 3
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
| libMesh Performance: Alive time=5.51913, Active time=5.23447                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0743      0.074305    0.1335      0.133454    1.42     2.55     |
|   build_constraint_matrix()        228       0.0102      0.000045    0.0102      0.000045    0.19     0.19     |
|   build_sparsity()                 1         0.0947      0.094666    0.1830      0.182996    1.81     3.50     |
|   cnstrn_elem_mat_vec()            228       0.0157      0.000069    0.0157      0.000069    0.30     0.30     |
|   create_dof_constraints()         1         0.4912      0.491150    1.3015      1.301506    9.38     24.86    |
|   distribute_dofs()                1         0.1693      0.169316    0.5815      0.581518    3.23     11.11    |
|   dof_indices()                    13923     0.8947      0.000064    0.8947      0.000064    17.09    17.09    |
|   prepare_send_list()              1         0.0009      0.000929    0.0009      0.000929    0.02     0.02     |
|   reinit()                         1         0.3978      0.397767    0.3978      0.397767    7.60     7.60     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          2         0.0220      0.010996    0.1678      0.083880    0.42     3.20     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0744      0.074370    0.0744      0.074370    1.42     1.42     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        441       0.0466      0.000106    0.0466      0.000106    0.89     0.89     |
|   init_shape_functions()           214       0.0015      0.000007    0.0015      0.000007    0.03     0.03     |
|   inverse_map()                    25764     0.2253      0.000009    0.2253      0.000009    4.30     4.30     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             441       0.0142      0.000032    0.0142      0.000032    0.27     0.27     |
|   compute_face_map()               213       0.0069      0.000032    0.0069      0.000032    0.13     0.13     |
|   init_face_shape_functions()      1         0.0001      0.000067    0.0001      0.000067    0.00     0.00     |
|   init_reference_to_physical_map() 214       0.0219      0.000102    0.0219      0.000102    0.42     0.42     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               1         0.0867      0.086678    0.0868      0.086816    1.66     1.66     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           37744     0.1306      0.000003    0.1306      0.000003    2.49     2.49     |
|   init()                           1         0.0070      0.006964    0.0070      0.006964    0.13     0.13     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 2         0.4473      0.223657    0.4550      0.227483    8.55     8.69     |
|   renumber_nodes_and_elem()        4         0.0180      0.004491    0.0180      0.004491    0.34     0.34     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        3         0.1287      0.042903    0.1287      0.042903    2.46     2.46     |
|   find_global_indices()            3         0.0470      0.015674    0.1900      0.063317    0.90     3.63     |
|   parallel_sort()                  3         0.0068      0.002281    0.0092      0.003070    0.13     0.18     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         2         0.0003      0.000140    0.3295      0.164731    0.01     6.29     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _refine_elements()               1         0.1724      0.172366    0.4594      0.459424    3.29     8.78     |
|   add_point()                      37744     0.1335      0.000004    0.2712      0.000007    2.55     5.18     |
|   make_refinement_compatible()     2         0.0055      0.002743    0.0055      0.002755    0.10     0.11     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0208      0.020820    0.0208      0.020820    0.40     0.40     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      2         1.0082      0.504078    1.1543      0.577174    19.26    22.05    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      13        0.0071      0.000544    0.0082      0.000632    0.14     0.16     |
|   max(bool)                        5         0.0093      0.001856    0.0093      0.001856    0.18     0.18     |
|   max(scalar)                      189       0.0018      0.000010    0.0018      0.000010    0.03     0.03     |
|   max(vector)                      44        0.0006      0.000014    0.0017      0.000038    0.01     0.03     |
|   min(bool)                        225       0.0017      0.000008    0.0017      0.000008    0.03     0.03     |
|   min(scalar)                      182       0.1150      0.000632    0.1150      0.000632    2.20     2.20     |
|   min(vector)                      44        0.0008      0.000017    0.0020      0.000046    0.01     0.04     |
|   probe()                          187       0.0263      0.000141    0.0263      0.000141    0.50     0.50     |
|   receive()                        187       0.0016      0.000008    0.0279      0.000149    0.03     0.53     |
|   send()                           187       0.0008      0.000004    0.0008      0.000004    0.01     0.01     |
|   send_receive()                   182       0.0020      0.000011    0.0192      0.000105    0.04     0.37     |
|   sum()                            26        0.0036      0.000137    0.0075      0.000288    0.07     0.14     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           187       0.0004      0.000002    0.0004      0.000002    0.01     0.01     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         2         0.0338      0.016912    0.0546      0.027318    0.65     1.04     |
|   set_parent_processor_ids()       2         0.0237      0.011851    0.0237      0.011851    0.45     0.45     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.1595      0.159485    0.1595      0.159485    3.05     3.05     |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0723      0.072318    0.2147      0.214729    1.38     4.10     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            118853    5.2345                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example subdomains_ex1:
*  mpirun -np 12 example-devel -d 3 -n 15 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
