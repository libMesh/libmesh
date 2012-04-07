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
        #include "libmesh.h"
        #include "mesh.h"
        #include "mesh_generation.h"
        #include "exodusII_io.h"
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
        
        #include "mesh_refinement.h"
        
</pre>
</div>
<div class = "comment">
Classes needed for subdomain computation.
</div>

<div class ="fragment">
<pre>
        #include "system_subset_by_subdomain.h"
        
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
            GnuPlotIO plot(mesh,"Example 4, 1D",GnuPlotIO::GRID_ON);
            plot.write_equation_systems("out_1",equation_systems);
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
          std::vector&lt;unsigned int&gt; dof_indices;
        
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
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  
  #include &lt;iostream&gt;
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
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
  
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_refinement.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;system_subset_by_subdomain.h&quot;</FONT></B>
  
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
      GnuPlotIO plot(mesh,<B><FONT COLOR="#BC8F8F">&quot;Example 4, 1D&quot;</FONT></B>,GnuPlotIO::GRID_ON);
      plot.write_equation_systems(<B><FONT COLOR="#BC8F8F">&quot;out_1&quot;</FONT></B>,equation_systems);
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
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
  
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
Compiling C++ (in optimized mode) subdomains_ex1.C...
Linking subdomains_ex1-opt...
***************************************************************
* Running Example  ./subdomains_ex1-opt
***************************************************************
 
Running ./subdomains_ex1-opt -d 1 -n 40

 Mesh Information:
  mesh_dimension()=1
  spatial_dimension()=3
  n_nodes()=43
    n_local_nodes()=43
  n_elem()=44
    n_local_elem()=44
    n_active_elem()=42
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
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
    n_local_dofs()=43
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 2.95349
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 3
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=1
  spatial_dimension()=3
  n_nodes()=43
    n_local_nodes()=43
  n_elem()=44
    n_local_elem()=44
    n_active_elem()=42
  n_subdomains()=2
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0


-------------------------------------------------------------------
| Time:           Sat Apr  7 16:03:06 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.001519, Active time=0.001018                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           34        0.0001      0.000003    0.0001      0.000003    11.59    11.59    |
| Fe                            34        0.0003      0.000010    0.0003      0.000010    32.71    32.71    |
| Ke                            34        0.0000      0.000001    0.0000      0.000001    2.36     2.36     |
| elem init                     34        0.0005      0.000014    0.0005      0.000014    45.58    45.58    |
| matrix insertion              34        0.0001      0.000002    0.0001      0.000002    7.76     7.76     |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       170       0.0010                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.079995, Active time=0.006212                                             |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 1         0.0001      0.000053    0.0001      0.000053    0.85     0.85     |
|   compute_sparsity()           1         0.0003      0.000296    0.0004      0.000358    4.76     5.76     |
|   create_dof_constraints()     1         0.0000      0.000002    0.0000      0.000002    0.03     0.03     |
|   distribute_dofs()            1         0.0001      0.000108    0.0003      0.000299    1.74     4.81     |
|   dof_indices()                152       0.0002      0.000001    0.0002      0.000001    2.43     2.43     |
|   prepare_send_list()          1         0.0000      0.000002    0.0000      0.000002    0.03     0.03     |
|   reinit()                     1         0.0002      0.000190    0.0002      0.000190    3.06     3.06     |
|                                                                                                            |
| EquationSystems                                                                                            |
|   build_solution_vector()      1         0.0002      0.000174    0.0002      0.000226    2.80     3.64     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         36        0.0000      0.000001    0.0000      0.000001    0.80     0.80     |
|   compute_face_map()           2         0.0000      0.000003    0.0000      0.000003    0.10     0.10     |
|   compute_shape_functions()    36        0.0000      0.000001    0.0000      0.000001    0.45     0.45     |
|   init_face_shape_functions()  1         0.0000      0.000009    0.0000      0.000009    0.14     0.14     |
|   init_shape_functions()       3         0.0001      0.000027    0.0001      0.000027    1.30     1.30     |
|                                                                                                            |
| GnuPlotIO                                                                                                  |
|   write_nodal_data()           1         0.0007      0.000730    0.0007      0.000730    11.75    11.75    |
|                                                                                                            |
| LocationMap                                                                                                |
|   find()                       4         0.0000      0.000003    0.0000      0.000003    0.18     0.18     |
|   init()                       1         0.0001      0.000053    0.0001      0.000053    0.85     0.85     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             2         0.0003      0.000163    0.0003      0.000163    5.26     5.26     |
|   renumber_nodes_and_elem()    4         0.0000      0.000003    0.0000      0.000003    0.23     0.23     |
|                                                                                                            |
| MeshOutput                                                                                                 |
|   write_equation_systems()     1         0.0000      0.000042    0.0010      0.001000    0.68     16.10    |
|                                                                                                            |
| MeshRefinement                                                                                             |
|   _refine_elements()           1         0.0001      0.000099    0.0001      0.000136    1.59     2.19     |
|   add_point()                  4         0.0000      0.000006    0.0000      0.000009    0.37     0.60     |
|   make_refinement_compatible() 1         0.0000      0.000024    0.0000      0.000024    0.39     0.39     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0002      0.000155    0.0002      0.000155    2.50     2.50     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  1         0.0000      0.000001    0.0000      0.000001    0.02     0.02     |
|                                                                                                            |
| Partitioner                                                                                                |
|   single_partition()           2         0.0000      0.000018    0.0000      0.000018    0.58     0.58     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      1         0.0020      0.002034    0.0020      0.002034    32.74    32.74    |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   1         0.0015      0.001513    0.0018      0.001758    24.36    28.30    |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        262       0.0062                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

Running ./subdomains_ex1-opt -d 2 -n 30

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1329
    n_local_nodes()=1329
  n_elem()=1268
    n_local_elem()=1268
    n_active_elem()=1176
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
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
    n_local_dofs()=1329
    n_constrained_dofs()=184
    n_local_constrained_dofs()=184
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 9.2167
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 15
      Maximum Off-Processor Bandwidth <= 0
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
    n_local_nodes()=1329
  n_elem()=1268
    n_local_elem()=1268
    n_active_elem()=1176
  n_subdomains()=2
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0


-------------------------------------------------------------------
| Time:           Sat Apr  7 16:03:06 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.022208, Active time=0.014648                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           628       0.0041      0.000007    0.0041      0.000007    28.12    28.12    |
| Fe                            628       0.0051      0.000008    0.0051      0.000008    34.71    34.71    |
| Ke                            628       0.0007      0.000001    0.0007      0.000001    4.49     4.49     |
| elem init                     628       0.0039      0.000006    0.0039      0.000006    26.65    26.65    |
| matrix insertion              628       0.0009      0.000001    0.0009      0.000001    6.03     6.03     |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       3140      0.0146                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.149911, Active time=0.074249                                             |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 1         0.0008      0.000768    0.0008      0.000768    1.03     1.03     |
|   build_constraint_matrix()    628       0.0011      0.000002    0.0011      0.000002    1.54     1.54     |
|   cnstrn_elem_mat_vec()        628       0.0030      0.000005    0.0030      0.000005    4.09     4.09     |
|   compute_sparsity()           1         0.0035      0.003460    0.0048      0.004785    4.66     6.44     |
|   create_dof_constraints()     1         0.0049      0.004853    0.0070      0.007009    6.54     9.44     |
|   distribute_dofs()            1         0.0010      0.001035    0.0040      0.004003    1.39     5.39     |
|   dof_indices()                5520      0.0041      0.000001    0.0041      0.000001    5.49     5.49     |
|   prepare_send_list()          1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                     1         0.0030      0.002967    0.0030      0.002967    4.00     4.00     |
|                                                                                                            |
| EquationSystems                                                                                            |
|   build_solution_vector()      2         0.0028      0.001404    0.0050      0.002489    3.78     6.70     |
|                                                                                                            |
| ExodusII_IO                                                                                                |
|   write_nodal_data()           1         0.0026      0.002579    0.0026      0.002579    3.47     3.47     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         828       0.0009      0.000001    0.0009      0.000001    1.16     1.16     |
|   compute_face_map()           200       0.0009      0.000004    0.0018      0.000009    1.17     2.44     |
|   compute_shape_functions()    828       0.0005      0.000001    0.0005      0.000001    0.61     0.61     |
|   init_face_shape_functions()  1         0.0000      0.000010    0.0000      0.000010    0.01     0.01     |
|   init_shape_functions()       201       0.0006      0.000003    0.0006      0.000003    0.86     0.86     |
|   inverse_map()                2072      0.0022      0.000001    0.0022      0.000001    2.92     2.92     |
|                                                                                                            |
| GMVIO                                                                                                      |
|   write_nodal_data()           1         0.0078      0.007781    0.0078      0.007781    10.48    10.48    |
|                                                                                                            |
| LocationMap                                                                                                |
|   find()                       1104      0.0016      0.000001    0.0016      0.000001    2.10     2.10     |
|   init()                       1         0.0006      0.000618    0.0006      0.000618    0.83     0.83     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             2         0.0075      0.003734    0.0075      0.003734    10.06    10.06    |
|   renumber_nodes_and_elem()    4         0.0004      0.000112    0.0004      0.000112    0.60     0.60     |
|                                                                                                            |
| MeshOutput                                                                                                 |
|   write_equation_systems()     2         0.0001      0.000027    0.0154      0.007697    0.07     20.73    |
|                                                                                                            |
| MeshRefinement                                                                                             |
|   _refine_elements()           1         0.0028      0.002788    0.0072      0.007191    3.75     9.68     |
|   add_point()                  1104      0.0021      0.000002    0.0041      0.000004    2.89     5.46     |
|   make_refinement_compatible() 1         0.0001      0.000106    0.0001      0.000106    0.14     0.14     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0016      0.001598    0.0016      0.001598    2.15     2.15     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                            |
| Partitioner                                                                                                |
|   single_partition()           2         0.0003      0.000154    0.0003      0.000154    0.42     0.42     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      1         0.0043      0.004284    0.0043      0.004284    5.77     5.77     |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   1         0.0134      0.013371    0.0224      0.022367    18.01    30.12    |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        13141     0.0742                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

Running ./subdomains_ex1-opt -d 3 -n 15

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=10854
    n_local_nodes()=10854
  n_elem()=8767
    n_local_elem()=8767
    n_active_elem()=8093
  n_subdomains()=1
  n_partitions()=1
  n_processors()=1
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
    n_local_dofs()=10854
    n_constrained_dofs()=4068
    n_local_constrained_dofs()=4068
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 29.7157
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 78
      Maximum Off-Processor Bandwidth <= 0
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
    n_local_nodes()=10854
  n_elem()=8767
    n_local_elem()=8767
    n_active_elem()=8093
  n_subdomains()=2
  n_partitions()=1
  n_processors()=1
  n_threads()=1
  processor_id()=0


-------------------------------------------------------------------
| Time:           Sat Apr  7 16:03:07 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.378895, Active time=0.297449                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           2931      0.1352      0.000046    0.1352      0.000046    45.44    45.44    |
| Fe                            2931      0.0807      0.000028    0.0807      0.000028    27.13    27.13    |
| Ke                            2931      0.0236      0.000008    0.0236      0.000008    7.94     7.94     |
| elem init                     2931      0.0389      0.000013    0.0389      0.000013    13.09    13.09    |
| matrix insertion              2931      0.0190      0.000006    0.0190      0.000006    6.40     6.40     |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       14655     0.2974                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=1.22689, Active time=1.08058                                               |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 1         0.0075      0.007451    0.0075      0.007451    0.69     0.69     |
|   build_constraint_matrix()    2931      0.0259      0.000009    0.0259      0.000009    2.40     2.40     |
|   cnstrn_elem_mat_vec()        2931      0.0379      0.000013    0.0379      0.000013    3.51     3.51     |
|   compute_sparsity()           1         0.0822      0.082197    0.1079      0.107898    7.61     9.99     |
|   create_dof_constraints()     1         0.1197      0.119702    0.2219      0.221950    11.08    20.54    |
|   distribute_dofs()            1         0.0113      0.011324    0.0324      0.032366    1.05     3.00     |
|   dof_indices()                40989     0.0586      0.000001    0.0586      0.000001    5.42     5.42     |
|   prepare_send_list()          1         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   reinit()                     1         0.0210      0.021041    0.0210      0.021041    1.95     1.95     |
|                                                                                                            |
| EquationSystems                                                                                            |
|   build_solution_vector()      2         0.0265      0.013267    0.0431      0.021559    2.46     3.99     |
|                                                                                                            |
| ExodusII_IO                                                                                                |
|   write_nodal_data()           1         0.0179      0.017940    0.0179      0.017940    1.66     1.66     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         5649      0.0168      0.000003    0.0168      0.000003    1.56     1.56     |
|   compute_face_map()           2718      0.0079      0.000003    0.0079      0.000003    0.73     0.73     |
|   compute_shape_functions()    5649      0.0097      0.000002    0.0097      0.000002    0.89     0.89     |
|   init_face_shape_functions()  1         0.0000      0.000020    0.0000      0.000020    0.00     0.00     |
|   init_shape_functions()       2719      0.0638      0.000023    0.0638      0.000023    5.90     5.90     |
|   inverse_map()                43392     0.0817      0.000002    0.0817      0.000002    7.56     7.56     |
|                                                                                                            |
| GMVIO                                                                                                      |
|   write_nodal_data()           1         0.0764      0.076440    0.0764      0.076440    7.07     7.07     |
|                                                                                                            |
| LocationMap                                                                                                |
|   find()                       37744     0.0236      0.000001    0.0236      0.000001    2.18     2.18     |
|   init()                       1         0.0018      0.001792    0.0018      0.001792    0.17     0.17     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             2         0.0464      0.023207    0.0464      0.023207    4.30     4.30     |
|   renumber_nodes_and_elem()    4         0.0063      0.001574    0.0063      0.001574    0.58     0.58     |
|                                                                                                            |
| MeshOutput                                                                                                 |
|   write_equation_systems()     2         0.0001      0.000037    0.1376      0.068786    0.01     12.73    |
|                                                                                                            |
| MeshRefinement                                                                                             |
|   _refine_elements()           1         0.0455      0.045523    0.1160      0.115999    4.21     10.73    |
|   add_point()                  37744     0.0337      0.000001    0.0639      0.000002    3.11     5.91     |
|   make_refinement_compatible() 1         0.0003      0.000284    0.0003      0.000284    0.03     0.03     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0054      0.005395    0.0054      0.005395    0.50     0.50     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                            |
| Partitioner                                                                                                |
|   single_partition()           2         0.0013      0.000635    0.0013      0.000635    0.12     0.12     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      1         0.0448      0.044756    0.0448      0.044756    4.14     4.14     |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   1         0.2066      0.206639    0.3791      0.379101    19.12    35.08    |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        182494    1.0806                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./subdomains_ex1-opt
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
