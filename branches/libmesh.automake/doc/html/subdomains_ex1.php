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
Linking subdomains_ex1-opt...
***************************************************************
* Running Example  mpirun -np 6 ./subdomains_ex1-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running ./subdomains_ex1-opt -d 1 -n 40 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=1
  spatial_dimension()=3
  n_nodes()=43
    n_local_nodes()=8
  n_elem()=44
    n_local_elem()=8
    n_active_elem()=42
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
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
    n_local_dofs()=8
    n_constrained_dofs()=0
    n_local_constrained_dofs()=0
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 2.75
      Average Off-Processor Bandwidth <= 0.125
      Maximum  On-Processor Bandwidth <= 3
      Maximum Off-Processor Bandwidth <= 1
    DofMap Constraints
      Number of DoF Constraints = 0
      Number of Node Constraints = 0

 Mesh Information:
  mesh_dimension()=1
  spatial_dimension()=3
  n_nodes()=43
    n_local_nodes()=8
  n_elem()=44
    n_local_elem()=8
    n_active_elem()=42
  n_subdomains()=2
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0


-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 6                                                |
| Time:           Fri Aug 24 15:26:12 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.000636, Active time=0.000407                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           3         0.0002      0.000074    0.0002      0.000074    54.55    54.55    |
| Fe                            3         0.0001      0.000020    0.0001      0.000020    14.74    14.74    |
| Ke                            3         0.0000      0.000000    0.0000      0.000000    0.25     0.25     |
| elem init                     3         0.0001      0.000027    0.0001      0.000027    19.90    19.90    |
| matrix insertion              3         0.0000      0.000014    0.0000      0.000014    10.57    10.57    |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       15        0.0004                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./subdomains_ex1-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:26:12 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           8.094e-02      2.87249   3.710e-02
Objects:              6.600e+01      1.00000   6.600e+01
Flops:                3.492e+03      2.87881   2.507e+03  1.504e+04
Flops/sec:            1.233e+05      5.99790   8.217e+04  4.930e+05
MPI Messages:         4.900e+01      1.92157   4.117e+01  2.470e+02
MPI Message Lengths:  4.000e+02      2.06186   7.951e+00  1.964e+03
MPI Reductions:       1.350e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 3.7057e-02  99.9%  1.5041e+04 100.0%  2.470e+02 100.0%  7.951e+00      100.0%  9.300e+01  68.9% 

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

VecMDot               11 1.0 1.3804e-04 2.4 9.90e+02 3.0 0.0e+00 0.0e+00 1.1e+01  0 28  0  0  8   0 28  0  0 12    31
VecNorm               13 1.0 9.2602e-04 7.8 2.08e+02 2.7 0.0e+00 0.0e+00 1.3e+01  1  6  0  0 10   1  6  0  0 14     1
VecScale              12 1.0 8.0109e-05 4.3 9.60e+01 2.7 0.0e+00 0.0e+00 0.0e+00  0  3  0  0  0   0  3  0  0  0     5
VecCopy                4 1.0 4.2915e-06 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                18 1.0 1.0252e-05 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                2 1.0 7.0598e-03 1.2 3.20e+01 2.7 0.0e+00 0.0e+00 0.0e+00 18  1  0  0  0  18  1  0  0  0     0
VecMAXPY              12 1.0 6.4373e-06 2.2 1.23e+03 2.7 0.0e+00 0.0e+00 0.0e+00  0 36  0  0  0   0 36  0  0  0   837
VecAssemblyBegin       3 1.0 1.1179e-03 1.3 0.00e+00 0.0 1.0e+01 6.0e+00 9.0e+00  3  0  4  3  7   3  0  4  3 10     0
VecAssemblyEnd         3 1.0 2.2650e-05 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       16 1.0 5.5742e-0411.7 0.00e+00 0.0 1.3e+02 8.3e+00 0.0e+00  0  0 53 55  0   0  0 53 55  0     0
VecScatterEnd         16 1.0 5.0998e-0420.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
VecNormalize          12 1.0 9.7609e-04 8.1 2.88e+02 2.7 0.0e+00 0.0e+00 1.2e+01  1  8  0  0  9   1  8  0  0 13     1
MatMult               12 1.0 6.0463e-04 4.6 4.80e+02 3.1 1.2e+02 8.0e+00 0.0e+00  1 14 49 49  0   1 14 49 49  0     3
MatSolve              12 1.0 2.7180e-05 6.7 4.32e+02 3.3 0.0e+00 0.0e+00 0.0e+00  0 12  0  0  0   0 12  0  0  0    67
MatLUFactorNum         1 1.0 2.0027e-05 1.7 2.20e+01 3.1 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0     5
MatILUFactorSym        1 1.0 7.0095e-05 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+00  0  0  0  0  1   0  0  0  0  1     0
MatAssemblyBegin       6 1.0 8.2803e-04 1.4 0.00e+00 0.0 1.5e+01 1.2e+01 8.0e+00  2  0  6  9  6   2  0  6  9  9     0
MatAssemblyEnd         6 1.0 9.6273e-04 1.2 0.00e+00 0.0 6.0e+01 4.0e+00 2.2e+01  2  0 24 12 16   2  0 24 12 24     0
MatGetRowIJ            1 1.0 1.2159e-0512.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetSubMatrice       2 1.0 3.7169e-04 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+01  1  0  0  0  7   1  0  0  0 11     0
MatGetOrdering         1 1.0 5.5075e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  1   0  0  0  0  2     0
MatZeroEntries         3 1.0 8.1062e-06 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        11 1.0 1.6499e-04 1.9 2.05e+03 2.8 0.0e+00 0.0e+00 1.1e+01  0 59  0  0  8   0 59  0  0 12    54
KSPSetup               2 1.0 7.2002e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 8.7020e-03 1.0 3.49e+03 2.9 1.2e+02 8.0e+00 2.7e+01 23100 49 49 20  23100 49 49 29     2
PCSetUp                2 1.0 7.5006e-04 2.1 2.20e+01 3.1 0.0e+00 0.0e+00 3.0e+00  1  1  0  0  2   1  1  0  0  3     0
PCSetUpOnBlocks        1 1.0 2.1482e-04 1.4 2.20e+01 3.1 0.0e+00 0.0e+00 3.0e+00  0  1  0  0  2   0  1  0  0  3     0
PCApply               12 1.0 1.9622e-04 1.2 4.32e+02 3.3 0.0e+00 0.0e+00 0.0e+00  0 12  0  0  0   0 12  0  0  0     9
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    29             29        39360     0
         Vec Scatter     5              5         4340     0
           Index Set    15             15         8128     0
   IS L to G Mapping     1              1          440     0
              Matrix    12             12        27924     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 4.85897e-05
Average time for zero size MPI_Send(): 4.33524e-05
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
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.272823, Active time=0.074086                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0000      0.000027    0.0000      0.000029    0.04     0.04     |
|   build_sparsity()                 1         0.0001      0.000072    0.0001      0.000089    0.10     0.12     |
|   create_dof_constraints()         1         0.0000      0.000000    0.0000      0.000000    0.00     0.00     |
|   distribute_dofs()                1         0.0001      0.000063    0.0005      0.000451    0.09     0.61     |
|   dof_indices()                    57        0.0000      0.000000    0.0000      0.000000    0.04     0.04     |
|   prepare_send_list()              1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                         1         0.0001      0.000054    0.0001      0.000054    0.07     0.07     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          1         0.0001      0.000131    0.0002      0.000239    0.18     0.32     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        4         0.0000      0.000002    0.0000      0.000002    0.01     0.01     |
|   init_shape_functions()           2         0.0000      0.000008    0.0000      0.000008    0.02     0.02     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             4         0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|   compute_face_map()               1         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   init_face_shape_functions()      1         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   init_reference_to_physical_map() 2         0.0002      0.000094    0.0002      0.000094    0.26     0.26     |
|                                                                                                                |
| GnuPlotIO                                                                                                      |
|   write_nodal_data()               1         0.0527      0.052705    0.0527      0.052705    71.14    71.14    |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           4         0.0000      0.000000    0.0000      0.000000    0.00     0.00     |
|   init()                           1         0.0000      0.000008    0.0000      0.000008    0.01     0.01     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 2         0.0001      0.000068    0.0003      0.000135    0.18     0.36     |
|   renumber_nodes_and_elem()        4         0.0000      0.000002    0.0000      0.000002    0.01     0.01     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        3         0.0003      0.000103    0.0003      0.000103    0.42     0.42     |
|   find_global_indices()            3         0.0001      0.000038    0.0035      0.001163    0.16     4.71     |
|   parallel_sort()                  3         0.0013      0.000418    0.0019      0.000642    1.69     2.60     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         1         0.0001      0.000107    0.0531      0.053051    0.14     71.61    |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _refine_elements()               1         0.0000      0.000048    0.0001      0.000147    0.06     0.20     |
|   add_point()                      4         0.0000      0.000001    0.0000      0.000002    0.01     0.01     |
|   make_refinement_compatible()     1         0.0000      0.000007    0.0001      0.000072    0.01     0.10     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0000      0.000049    0.0000      0.000049    0.07     0.07     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      2         0.0004      0.000215    0.0023      0.001136    0.58     3.07     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      11        0.0005      0.000042    0.0005      0.000042    0.63     0.63     |
|   broadcast()                      2         0.0000      0.000006    0.0000      0.000006    0.02     0.02     |
|   gather()                         2         0.0001      0.000030    0.0001      0.000030    0.08     0.08     |
|   max(bool)                        2         0.0001      0.000063    0.0001      0.000063    0.17     0.17     |
|   max(scalar)                      3         0.0004      0.000125    0.0004      0.000125    0.50     0.50     |
|   max(vector)                      3         0.0001      0.000020    0.0001      0.000020    0.08     0.08     |
|   min(bool)                        1         0.0001      0.000065    0.0001      0.000065    0.09     0.09     |
|   min(vector)                      3         0.0004      0.000120    0.0004      0.000120    0.49     0.49     |
|   probe()                          75        0.0015      0.000020    0.0015      0.000020    2.06     2.06     |
|   receive()                        75        0.0001      0.000001    0.0016      0.000022    0.13     2.21     |
|   send()                           75        0.0000      0.000001    0.0000      0.000001    0.06     0.06     |
|   send_receive()                   76        0.0001      0.000002    0.0018      0.000024    0.20     2.46     |
|   sum()                            13        0.0010      0.000075    0.0010      0.000075    1.31     1.31     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           75        0.0000      0.000000    0.0000      0.000000    0.05     0.05     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         2         0.0000      0.000022    0.0007      0.000374    0.06     1.01     |
|   set_parent_processor_ids()       2         0.0000      0.000005    0.0000      0.000005    0.01     0.01     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.0133      0.013301    0.0133      0.013301    17.95    17.95    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0006      0.000600    0.0008      0.000829    0.81     1.12     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            531       0.0741                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

Running ./subdomains_ex1-opt -d 2 -n 30 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1329
    n_local_nodes()=242
  n_elem()=1268
    n_local_elem()=214
    n_active_elem()=1176
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
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
    n_local_dofs()=242
    n_constrained_dofs()=184
    n_local_constrained_dofs()=36
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 8.82645
      Average Off-Processor Bandwidth <= 0.400826
      Maximum  On-Processor Bandwidth <= 15
      Maximum Off-Processor Bandwidth <= 7
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
    n_local_nodes()=242
  n_elem()=1268
    n_local_elem()=214
    n_active_elem()=1176
  n_subdomains()=2
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0


-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 6                                                |
| Time:           Fri Aug 24 15:26:13 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.009649, Active time=0.001734                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           82        0.0006      0.000007    0.0006      0.000007    32.87    32.87    |
| Fe                            82        0.0004      0.000005    0.0004      0.000005    22.03    22.03    |
| Ke                            82        0.0001      0.000001    0.0001      0.000001    4.56     4.56     |
| elem init                     82        0.0006      0.000007    0.0006      0.000007    33.04    33.04    |
| matrix insertion              82        0.0001      0.000002    0.0001      0.000002    7.50     7.50     |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       410       0.0017                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./subdomains_ex1-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:26:13 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.053e-01      1.02299   1.036e-01
Objects:              7.600e+01      1.00000   7.600e+01
Flops:                4.072e+05      1.54215   3.323e+05  1.994e+06
Flops/sec:            3.955e+06      1.54215   3.208e+06  1.925e+07
MPI Messages:         1.445e+02      2.44915   9.700e+01  5.820e+02
MPI Message Lengths:  1.283e+04      1.88507   9.534e+01  5.549e+04
MPI Reductions:       1.460e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.0358e-01 100.0%  1.9937e+06 100.0%  5.820e+02 100.0%  9.534e+01      100.0%  1.040e+02  71.2% 

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

VecMDot               16 1.0 1.6222e-03 9.7 4.72e+04 1.6 0.0e+00 0.0e+00 1.6e+01  1 10  0  0 11   1 10  0  0 15   128
VecNorm               18 1.0 3.2706e-03 3.3 6.26e+03 1.6 0.0e+00 0.0e+00 1.8e+01  3  1  0  0 12   3  1  0  0 17     8
VecScale              17 1.0 2.6703e-05 1.5 2.96e+03 1.6 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0   490
VecCopy                4 1.0 1.8835e-05 8.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                23 1.0 2.5988e-05 3.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                2 1.0 3.8862e-05 1.5 6.96e+02 1.6 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0    79
VecMAXPY              17 1.0 2.2650e-05 1.8 5.29e+04 1.6 0.0e+00 0.0e+00 0.0e+00  0 12  0  0  0   0 12  0  0  0 10321
VecAssemblyBegin       3 1.0 7.1120e-04 1.7 0.00e+00 0.0 2.0e+01 1.4e+02 9.0e+00  1  0  3  5  6   1  0  3  5  9     0
VecAssemblyEnd         3 1.0 2.2411e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       21 1.0 9.0837e-05 1.4 0.00e+00 0.0 3.6e+02 8.4e+01 0.0e+00  0  0 62 54  0   0  0 62 54  0     0
VecScatterEnd         21 1.0 5.1453e-03 5.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  2  0  0  0  0   2  0  0  0  0     0
VecNormalize          17 1.0 3.3307e-03 3.3 8.87e+03 1.6 0.0e+00 0.0e+00 1.7e+01  3  2  0  0 12   3  2  0  0 16    12
MatMult               17 1.0 5.0623e-03 4.9 4.98e+04 1.7 3.4e+02 8.0e+01 0.0e+00  2 11 58 49  0   2 11 58 49  0    43
MatSolve              17 1.0 1.5426e-04 1.8 1.65e+05 1.5 0.0e+00 0.0e+00 0.0e+00  0 42  0  0  0   0 42  0  0  0  5412
MatLUFactorNum         1 1.0 1.5497e-04 1.7 1.07e+05 2.7 0.0e+00 0.0e+00 0.0e+00  0 23  0  0  0   0 23  0  0  0  2946
MatILUFactorSym        1 1.0 3.9291e-04 1.7 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+00  0  0  0  0  1   0  0  0  0  1     0
MatAssemblyBegin       6 1.0 2.3081e-03 1.3 0.00e+00 0.0 3.0e+01 4.1e+02 8.0e+00  2  0  5 22  5   2  0  5 22  8     0
MatAssemblyEnd         6 1.0 1.7190e-03 1.1 0.00e+00 0.0 1.2e+02 2.2e+01 2.2e+01  2  0 21  5 15   2  0 21  5 21     0
MatGetRowIJ            1 1.0 1.2875e-0513.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetSubMatrice       2 1.0 8.9097e-04 3.2 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+01  0  0  0  0  7   0  0  0  0 10     0
MatGetOrdering         1 1.0 1.6940e-0347.1 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  1   0  0  0  0  2     0
MatZeroEntries         3 1.0 1.4067e-05 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        16 1.0 1.6649e-03 8.0 9.45e+04 1.6 0.0e+00 0.0e+00 1.6e+01  1 21  0  0 11   1 21  0  0 15   251
KSPSetup               2 1.0 9.4891e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 8.1079e-03 1.0 4.07e+05 1.5 3.4e+02 8.0e+01 3.7e+01  8100 58 49 25   8100 58 49 36   246
PCSetUp                2 1.0 2.5170e-03 3.8 1.07e+05 2.7 0.0e+00 0.0e+00 3.0e+00  1 23  0  0  2   1 23  0  0  3   181
PCSetUpOnBlocks        1 1.0 2.2840e-03 5.2 1.07e+05 2.7 0.0e+00 0.0e+00 3.0e+00  1 23  0  0  2   1 23  0  0  3   200
PCApply               17 1.0 3.3832e-04 1.2 1.65e+05 1.5 0.0e+00 0.0e+00 0.0e+00  0 42  0  0  0   0 42  0  0  0  2468
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    39             39        85000     0
         Vec Scatter     5              5         4340     0
           Index Set    15             15        15900     0
   IS L to G Mapping     1              1         1512     0
              Matrix    12             12       173064     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 4.39644e-05
Average time for zero size MPI_Send(): 7.43469e-05
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
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.280124, Active time=0.086886                                                 |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0001      0.000145    0.0002      0.000165    0.17     0.19     |
|   build_constraint_matrix()        82        0.0001      0.000001    0.0001      0.000001    0.10     0.10     |
|   build_sparsity()                 1         0.0016      0.001639    0.0021      0.002051    1.89     2.36     |
|   cnstrn_elem_mat_vec()            82        0.0072      0.000088    0.0072      0.000088    8.33     8.33     |
|   create_dof_constraints()         1         0.0016      0.001592    0.0022      0.002166    1.83     2.49     |
|   distribute_dofs()                1         0.0004      0.000382    0.0043      0.004323    0.44     4.98     |
|   dof_indices()                    2528      0.0008      0.000000    0.0008      0.000000    0.87     0.87     |
|   prepare_send_list()              1         0.0000      0.000004    0.0000      0.000004    0.00     0.00     |
|   reinit()                         1         0.0008      0.000767    0.0008      0.000767    0.88     0.88     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          2         0.0005      0.000260    0.0058      0.002919    0.60     6.72     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0018      0.001810    0.0018      0.001810    2.08     2.08     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        119       0.0002      0.000002    0.0002      0.000002    0.22     0.22     |
|   init_shape_functions()           38        0.0000      0.000001    0.0000      0.000001    0.04     0.04     |
|   inverse_map()                    1031      0.0004      0.000000    0.0004      0.000000    0.41     0.41     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             119       0.0001      0.000001    0.0001      0.000001    0.12     0.12     |
|   compute_face_map()               37        0.0001      0.000003    0.0002      0.000006    0.14     0.26     |
|   init_face_shape_functions()      1         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   init_reference_to_physical_map() 38        0.0001      0.000002    0.0001      0.000002    0.09     0.09     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               1         0.0041      0.004105    0.0041      0.004105    4.72     4.72     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           1104      0.0004      0.000000    0.0004      0.000000    0.46     0.46     |
|   init()                           1         0.0001      0.000109    0.0001      0.000109    0.13     0.13     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 2         0.0017      0.000861    0.0071      0.003553    1.98     8.18     |
|   renumber_nodes_and_elem()        4         0.0002      0.000047    0.0002      0.000047    0.22     0.22     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        3         0.0049      0.001632    0.0049      0.001632    5.63     5.63     |
|   find_global_indices()            3         0.0005      0.000178    0.0161      0.005351    0.61     18.47    |
|   parallel_sort()                  3         0.0074      0.002471    0.0094      0.003143    8.53     10.85    |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         2         0.0003      0.000151    0.0121      0.006027    0.35     13.87    |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _refine_elements()               1         0.0009      0.000941    0.0026      0.002602    1.08     2.99     |
|   add_point()                      1104      0.0007      0.000001    0.0012      0.000001    0.82     1.38     |
|   make_refinement_compatible()     1         0.0000      0.000024    0.0001      0.000075    0.03     0.09     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0003      0.000299    0.0003      0.000299    0.34     0.34     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      2         0.0038      0.001900    0.0175      0.008735    4.37     20.11    |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      11        0.0007      0.000067    0.0007      0.000067    0.85     0.85     |
|   broadcast()                      2         0.0000      0.000007    0.0000      0.000007    0.02     0.02     |
|   gather()                         2         0.0003      0.000149    0.0003      0.000149    0.34     0.34     |
|   max(bool)                        2         0.0004      0.000203    0.0004      0.000203    0.47     0.47     |
|   max(scalar)                      3         0.0059      0.001955    0.0059      0.001955    6.75     6.75     |
|   max(vector)                      3         0.0006      0.000195    0.0006      0.000195    0.67     0.67     |
|   min(bool)                        1         0.0001      0.000051    0.0001      0.000051    0.06     0.06     |
|   min(vector)                      3         0.0028      0.000942    0.0028      0.000942    3.25     3.25     |
|   probe()                          75        0.0097      0.000129    0.0097      0.000129    11.17    11.17    |
|   receive()                        75        0.0001      0.000002    0.0098      0.000131    0.14     11.32    |
|   send()                           75        0.0001      0.000001    0.0001      0.000001    0.07     0.07     |
|   send_receive()                   76        0.0001      0.000002    0.0099      0.000130    0.16     11.39    |
|   sum()                            16        0.0072      0.000451    0.0072      0.000451    8.30     8.30     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           75        0.0000      0.000001    0.0000      0.000001    0.05     0.05     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         2         0.0003      0.000144    0.0064      0.003184    0.33     7.33     |
|   set_parent_processor_ids()       2         0.0001      0.000066    0.0001      0.000066    0.15     0.15     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.0153      0.015335    0.0153      0.015335    17.65    17.65    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0018      0.001790    0.0098      0.009841    2.06     11.33    |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            6741      0.0869                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

Running ./subdomains_ex1-opt -d 3 -n 15 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=10854
    n_local_nodes()=2023
  n_elem()=8767
    n_local_elem()=1466
    n_active_elem()=8093
  n_subdomains()=1
  n_partitions()=6
  n_processors()=6
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
    n_local_dofs()=2023
    n_constrained_dofs()=4068
    n_local_constrained_dofs()=738
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 27.8156
      Average Off-Processor Bandwidth <= 2.60356
      Maximum  On-Processor Bandwidth <= 78
      Maximum Off-Processor Bandwidth <= 45
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
    n_local_nodes()=2023
  n_elem()=8767
    n_local_elem()=1466
    n_active_elem()=8093
  n_subdomains()=2
  n_partitions()=6
  n_processors()=6
  n_threads()=1
  processor_id()=0


-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 6                                                |
| Time:           Fri Aug 24 15:26:14 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.050237, Active time=0.036437                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           497       0.0117      0.000024    0.0117      0.000024    32.22    32.22    |
| Fe                            497       0.0067      0.000014    0.0067      0.000014    18.46    18.46    |
| Ke                            497       0.0049      0.000010    0.0049      0.000010    13.46    13.46    |
| elem init                     497       0.0105      0.000021    0.0105      0.000021    28.71    28.71    |
| matrix insertion              497       0.0026      0.000005    0.0026      0.000005    7.15     7.15     |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       2485      0.0364                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./subdomains_ex1-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:26:15 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.429e+00      1.00566   1.423e+00
Objects:              8.500e+01      1.00000   8.500e+01
Flops:                1.926e+08      2.42527   1.184e+08  7.102e+08
Flops/sec:            1.348e+08      2.41195   8.312e+07  4.987e+08
MPI Messages:         1.895e+02      1.26756   1.753e+02  1.052e+03
MPI Message Lengths:  2.107e+05      1.23076   1.092e+03  1.148e+06
MPI Reductions:       1.620e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.4232e+00 100.0%  7.1016e+08 100.0%  1.052e+03 100.0%  1.092e+03      100.0%  1.200e+02  74.1% 

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

VecMDot               24 1.0 5.3639e-02 7.3 5.98e+05 1.3 0.0e+00 0.0e+00 2.4e+01  3  0  0  0 15   3  0  0  0 20    57
VecNorm               26 1.0 1.4311e-01240.3 5.19e+04 1.3 0.0e+00 0.0e+00 2.6e+01  2  0  0  0 16   2  0  0  0 22     2
VecScale              25 1.0 6.0081e-05 1.2 2.50e+04 1.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2129
VecCopy                4 1.0 1.8597e-05 1.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                31 1.0 5.6505e-05 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                2 1.0 1.4877e-04 4.4 3.99e+03 1.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   138
VecMAXPY              25 1.0 2.5272e-04 1.3 6.47e+05 1.3 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0 13118
VecAssemblyBegin       3 1.0 7.2987e-03 1.1 0.00e+00 0.0 3.4e+01 1.4e+03 9.0e+00  0  0  3  4  6   0  0  3  4  8     0
VecAssemblyEnd         3 1.0 4.1246e-05 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       29 1.0 5.0473e-04 2.0 0.00e+00 0.0 7.3e+02 6.9e+02 0.0e+00  0  0 69 44  0   0  0 69 44  0     0
VecScatterEnd         29 1.0 3.1743e-01 7.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 14  0  0  0  0  14  0  0  0  0     0
VecNormalize          25 1.0 1.4324e-01230.1 7.48e+04 1.3 0.0e+00 0.0e+00 2.5e+01  2  0  0  0 15   2  0  0  0 21     3
MatMult               25 1.0 3.1949e-01 7.4 1.30e+06 1.4 7.0e+02 6.6e+02 0.0e+00 15  1 67 40  0  15  1 67 40  0    20
MatSolve              25 1.0 6.5646e-02 2.4 2.91e+07 1.7 0.0e+00 0.0e+00 0.0e+00  3 18  0  0  0   3 18  0  0  0  1953
MatLUFactorNum         1 1.0 1.1578e-01 2.2 1.61e+08 2.7 0.0e+00 0.0e+00 0.0e+00  5 80  0  0  0   5 80  0  0  0  4912
MatILUFactorSym        1 1.0 3.1307e-01 2.8 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+00 12  0  0  0  1  12  0  0  0  1     0
MatAssemblyBegin       6 1.0 2.5902e-03 3.2 0.00e+00 0.0 5.1e+01 9.9e+03 8.0e+00  0  0  5 44  5   0  0  5 44  7     0
MatAssemblyEnd         6 1.0 3.6426e-03 1.1 0.00e+00 0.0 1.7e+02 1.7e+02 2.2e+01  0  0 16  2 14   0  0 16  2 18     0
MatGetRowIJ            1 1.0 1.0395e-04109.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetSubMatrice       2 1.0 1.6730e-03 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+01  0  0  0  0  6   0  0  0  0  8     0
MatGetOrdering         1 1.0 1.5402e-04 3.6 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  1   0  0  0  0  2     0
MatZeroEntries         3 1.0 1.2360e-03 4.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        24 1.0 5.3897e-02 7.1 1.20e+06 1.3 0.0e+00 0.0e+00 2.4e+01  3  1  0  0 15   3  1  0  0 20   114
KSPSetup               2 1.0 1.1086e-04 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 5.6816e-01 1.0 1.93e+08 2.4 7.0e+02 6.6e+02 5.3e+01 40100 67 40 33  40100 67 40 44  1250
PCSetUp                2 1.0 4.2938e-01 2.6 1.61e+08 2.7 0.0e+00 0.0e+00 3.0e+00 17 80  0  0  2  17 80  0  0  2  1324
PCSetUpOnBlocks        1 1.0 4.2901e-01 2.6 1.61e+08 2.7 0.0e+00 0.0e+00 3.0e+00 17 80  0  0  2  17 80  0  0  2  1326
PCApply               25 1.0 6.6800e-02 2.4 2.91e+07 1.7 0.0e+00 0.0e+00 0.0e+00  3 18  0  0  0   3 18  0  0  0  1919
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    48             48       416616     0
         Vec Scatter     5              5         4340     0
           Index Set    15             15        68944     0
   IS L to G Mapping     1              1        10892     0
              Matrix    12             12      9191792     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 1.00136e-05
Average time for zero size MPI_Send(): 3.24647e-05
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
 ----------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=1.53425, Active time=1.24492                                                   |
 ----------------------------------------------------------------------------------------------------------------
| Event                              nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                              w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|----------------------------------------------------------------------------------------------------------------|
|                                                                                                                |
|                                                                                                                |
| DofMap                                                                                                         |
|   add_neighbors_to_send_list()     1         0.0024      0.002382    0.0029      0.002852    0.19     0.23     |
|   build_constraint_matrix()        497       0.0022      0.000004    0.0022      0.000004    0.17     0.17     |
|   build_sparsity()                 1         0.0307      0.030712    0.0351      0.035096    2.47     2.82     |
|   cnstrn_elem_mat_vec()            497       0.0087      0.000018    0.0087      0.000018    0.70     0.70     |
|   create_dof_constraints()         1         0.0532      0.053218    0.0848      0.084807    4.27     6.81     |
|   distribute_dofs()                1         0.0056      0.005623    0.0449      0.044907    0.45     3.61     |
|   dof_indices()                    23398     0.0114      0.000000    0.0114      0.000000    0.92     0.92     |
|   prepare_send_list()              1         0.0002      0.000153    0.0002      0.000153    0.01     0.01     |
|   reinit()                         1         0.0117      0.011666    0.0117      0.011666    0.94     0.94     |
|                                                                                                                |
| EquationSystems                                                                                                |
|   build_solution_vector()          2         0.0043      0.002129    0.0199      0.009971    0.34     1.60     |
|                                                                                                                |
| ExodusII_IO                                                                                                    |
|   write_nodal_data()               1         0.0104      0.010429    0.0104      0.010429    0.84     0.84     |
|                                                                                                                |
| FE                                                                                                             |
|   compute_shape_functions()        962       0.0075      0.000008    0.0075      0.000008    0.61     0.61     |
|   init_shape_functions()           466       0.0003      0.000001    0.0003      0.000001    0.02     0.02     |
|   inverse_map()                    25764     0.0238      0.000001    0.0238      0.000001    1.91     1.91     |
|                                                                                                                |
| FEMap                                                                                                          |
|   compute_affine_map()             962       0.0024      0.000002    0.0024      0.000002    0.19     0.19     |
|   compute_face_map()               465       0.0012      0.000003    0.0012      0.000003    0.09     0.09     |
|   init_face_shape_functions()      1         0.0000      0.000008    0.0000      0.000008    0.00     0.00     |
|   init_reference_to_physical_map() 466       0.0043      0.000009    0.0043      0.000009    0.35     0.35     |
|                                                                                                                |
| GMVIO                                                                                                          |
|   write_nodal_data()               1         0.0345      0.034539    0.0345      0.034539    2.77     2.77     |
|                                                                                                                |
| LocationMap                                                                                                    |
|   find()                           37744     0.0109      0.000000    0.0109      0.000000    0.87     0.87     |
|   init()                           1         0.0005      0.000467    0.0005      0.000467    0.04     0.04     |
|                                                                                                                |
| Mesh                                                                                                           |
|   find_neighbors()                 2         0.0218      0.010877    0.0532      0.026603    1.75     4.27     |
|   renumber_nodes_and_elem()        4         0.0035      0.000866    0.0035      0.000866    0.28     0.28     |
|                                                                                                                |
| MeshCommunication                                                                                              |
|   compute_hilbert_indices()        3         0.0273      0.009094    0.0273      0.009094    2.19     2.19     |
|   find_global_indices()            3         0.0035      0.001154    0.0751      0.025038    0.28     6.03     |
|   parallel_sort()                  3         0.0162      0.005406    0.0340      0.011318    1.30     2.73     |
|                                                                                                                |
| MeshOutput                                                                                                     |
|   write_equation_systems()         2         0.0000      0.000025    0.0650      0.032481    0.00     5.22     |
|                                                                                                                |
| MeshRefinement                                                                                                 |
|   _refine_elements()               1         0.0302      0.030183    0.1425      0.142513    2.42     11.45    |
|   add_point()                      37744     0.0244      0.000001    0.0379      0.000001    1.96     3.04     |
|   make_refinement_compatible()     1         0.0001      0.000104    0.0009      0.000947    0.01     0.08     |
|                                                                                                                |
| MeshTools::Generation                                                                                          |
|   build_cube()                     1         0.0014      0.001407    0.0014      0.001407    0.11     0.11     |
|                                                                                                                |
| MetisPartitioner                                                                                               |
|   partition()                      2         0.0300      0.015010    0.0856      0.042785    2.41     6.87     |
|                                                                                                                |
| Parallel                                                                                                       |
|   allgather()                      11        0.0207      0.001879    0.0207      0.001879    1.66     1.66     |
|   broadcast()                      2         0.0000      0.000011    0.0000      0.000011    0.00     0.00     |
|   gather()                         2         0.0001      0.000042    0.0001      0.000042    0.01     0.01     |
|   max(bool)                        2         0.0721      0.036049    0.0721      0.036049    5.79     5.79     |
|   max(scalar)                      3         0.0752      0.025066    0.0752      0.025066    6.04     6.04     |
|   max(vector)                      3         0.0001      0.000039    0.0001      0.000039    0.01     0.01     |
|   min(bool)                        1         0.0008      0.000843    0.0008      0.000843    0.07     0.07     |
|   min(vector)                      3         0.0028      0.000949    0.0028      0.000949    0.23     0.23     |
|   probe()                          75        0.0459      0.000612    0.0459      0.000612    3.68     3.68     |
|   receive()                        75        0.0003      0.000004    0.0461      0.000615    0.02     3.71     |
|   send()                           75        0.0002      0.000003    0.0002      0.000003    0.02     0.02     |
|   send_receive()                   76        0.0002      0.000003    0.0426      0.000560    0.02     3.42     |
|   sum()                            16        0.0360      0.002252    0.0360      0.002252    2.89     2.89     |
|                                                                                                                |
| Parallel::Request                                                                                              |
|   wait()                           75        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                                |
| Partitioner                                                                                                    |
|   set_node_processor_ids()         2         0.0029      0.001430    0.0227      0.011359    0.23     1.82     |
|   set_parent_processor_ids()       2         0.0010      0.000499    0.0010      0.000499    0.08     0.08     |
|                                                                                                                |
| PetscLinearSolver                                                                                              |
|   solve()                          1         0.5793      0.579312    0.5793      0.579312    46.53    46.53    |
|                                                                                                                |
| System                                                                                                         |
|   assemble()                       1         0.0228      0.022832    0.0504      0.050439    1.83     4.05     |
 ----------------------------------------------------------------------------------------------------------------
| Totals:                            129424    1.2449                                          100.00            |
 ----------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 6 ./subdomains_ex1-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
