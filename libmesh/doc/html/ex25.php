<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex25",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 25 - Solving a 1D, 2D, or 3D Poisson on a subdomain</h1>

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
              "out_3.exd" : "out_2.exd",equation_systems);
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
        <B><FONT COLOR="#BC8F8F">&quot;out_3.exd&quot;</FONT></B> : <B><FONT COLOR="#BC8F8F">&quot;out_2.exd&quot;</FONT></B>,equation_systems);
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
* Running Example  mpirun -np 2 ./ex25-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running ./ex25-opt -d 1 -n 40 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=1
  spatial_dimension()=3
  n_nodes()=43
    n_local_nodes()=22
  n_elem()=44
    n_local_elem()=22
    n_active_elem()=42
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Poisson"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=43
    n_local_dofs()=22
    n_constrained_dofs()=0
    n_vectors()=1

 Mesh Information:
  mesh_dimension()=1
  spatial_dimension()=3
  n_nodes()=43
    n_local_nodes()=22
  n_elem()=44
    n_local_elem()=22
    n_active_elem()=42
  n_subdomains()=2
  n_processors()=2
  processor_id()=0


-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 2                                                |
| Time:           Thu Feb  3 12:13:05 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.000383, Active time=0.000218                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           17        0.0000      0.000002    0.0000      0.000002    15.14    15.14    |
| Fe                            17        0.0000      0.000002    0.0000      0.000002    17.89    17.89    |
| Ke                            17        0.0000      0.000000    0.0000      0.000000    1.38     1.38     |
| elem init                     17        0.0001      0.000007    0.0001      0.000007    54.59    54.59    |
| matrix insertion              17        0.0000      0.000001    0.0000      0.000001    11.01    11.01    |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       85        0.0002                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex25-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Thu Feb  3 12:13:05 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.040e-02      1.03249   1.023e-02
Objects:              6.800e+01      1.00000   6.800e+01
Flops:                2.897e+03      1.10742   2.756e+03  5.513e+03
Flops/sec:            2.787e+05      1.07257   2.693e+05  5.385e+05
MPI Messages:         2.350e+01      1.00000   2.350e+01  4.700e+01
MPI Message Lengths:  4.620e+02      1.03587   1.932e+01  9.080e+02
MPI Reductions:       1.260e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.0205e-02  99.7%  5.5130e+03 100.0%  4.700e+01 100.0%  1.932e+01      100.0%  8.400e+01  66.7% 

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

VecMDot                5 1.0 2.0266e-05 1.1 5.25e+02 1.1 0.0e+00 0.0e+00 5.0e+00  0 19  0  0  4   0 19  0  0  6    50
VecNorm                7 1.0 3.6716e-05 1.1 2.52e+02 1.1 0.0e+00 0.0e+00 7.0e+00  0  9  0  0  6   0  9  0  0  8    13
VecScale               6 1.0 1.3113e-05 1.2 1.08e+02 1.1 0.0e+00 0.0e+00 0.0e+00  0  4  0  0  0   0  4  0  0  0    16
VecCopy                3 1.0 3.0994e-06 3.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                12 1.0 5.9605e-06 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                2 1.0 3.3269e-03 1.0 7.20e+01 1.1 0.0e+00 0.0e+00 0.0e+00 33  3  0  0  0  33  3  0  0  0     0
VecMAXPY               6 1.0 6.9141e-06 1.2 7.20e+02 1.1 0.0e+00 0.0e+00 0.0e+00  0 25  0  0  0   0 25  0  0  0   202
VecAssemblyBegin       3 1.0 3.9816e-05 1.0 0.00e+00 0.0 2.0e+00 6.0e+00 9.0e+00  0  0  4  1  7   0  0  4  1 11     0
VecAssemblyEnd         3 1.0 1.0967e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       11 1.0 2.9802e-05 1.1 0.00e+00 0.0 1.4e+01 8.6e+00 0.0e+00  0  0 30 13  0   0  0 30 13  0     0
VecScatterEnd         11 1.0 1.7881e-05 2.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize           6 1.0 5.6982e-05 1.1 3.24e+02 1.1 0.0e+00 0.0e+00 6.0e+00  1 11  0  0  5   1 11  0  0  7    11
MatMult                6 1.0 4.5538e-05 1.1 5.28e+02 1.1 1.2e+01 8.0e+00 0.0e+00  0 19 26 11  0   0 19 26 11  0    23
MatSolve               6 1.0 3.0994e-06 3.2 6.12e+02 1.3 0.0e+00 0.0e+00 0.0e+00  0 20  0  0  0   0 20  0  0  0   354
MatLUFactorNum         1 1.0 1.7881e-05 1.1 8.00e+01 1.6 0.0e+00 0.0e+00 0.0e+00  0  2  0  0  0   0  2  0  0  0     7
MatILUFactorSym        1 1.0 4.0054e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+00  0  0  0  0  1   0  0  0  0  1     0
MatAssemblyBegin       6 1.0 1.2302e-04 1.0 0.00e+00 0.0 3.0e+00 1.2e+01 8.0e+00  1  0  6  4  6   1  0  6  4 10     0
MatAssemblyEnd         6 1.0 2.9302e-04 1.0 0.00e+00 0.0 1.2e+01 4.0e+00 2.2e+01  3  0 26  5 17   3  0 26  5 26     0
MatGetRowIJ            1 1.0 9.5367e-07 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetSubMatrice       2 1.0 1.1015e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+01  1  0  0  0  8   1  0  0  0 12     0
MatGetOrdering         1 1.0 2.3842e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries         3 1.0 7.1526e-06 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog         5 1.0 4.1246e-05 1.0 1.06e+03 1.1 0.0e+00 0.0e+00 5.0e+00  0 38  0  0  4   0 38  0  0  6    50
KSPSetup               2 1.0 4.3154e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 4.0510e-03 1.0 2.90e+03 1.1 1.2e+01 8.0e+00 1.5e+01 40100 26 11 12  40100 26 11 18     1
PCSetUp                2 1.0 2.9564e-04 1.0 8.00e+01 1.6 0.0e+00 0.0e+00 3.0e+00  3  2  0  0  2   3  2  0  0  4     0
PCSetUpOnBlocks        1 1.0 1.3304e-04 1.0 8.00e+01 1.6 0.0e+00 0.0e+00 3.0e+00  1  2  0  0  2   1  2  0  0  4     1
PCApply                6 1.0 7.9155e-05 1.0 6.12e+02 1.3 0.0e+00 0.0e+00 0.0e+00  1 20  0  0  0   1 20  0  0  0    14
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    29             29        41600     0
         Vec Scatter     6              6         5208     0
           Index Set    16             16         8816     0
   IS L to G Mapping     1              1          496     0
              Matrix    12             12        32836     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 1.19209e-06
Average time for zero size MPI_Send(): 4.52995e-06
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
Configure run at: Fri Oct 15 13:01:23 2010
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared=1 --with-mpi-dir=/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid --with-mumps=true --download-mumps=ifneeded --with-parmetis=true --download-parmetis=ifneeded --with-superlu=true --download-superlu=ifneeded --with-superludir=true --download-superlu_dist=ifneeded --with-blacs=true --download-blacs=ifneeded --with-scalapack=true --download-scalapack=ifneeded --with-hypre=true --download-hypre=ifneeded --with-blas-lib="[/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_intel_lp64.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_sequential.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_core.so]" --with-lapack-lib=/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_solver_lp64_sequential.a
-----------------------------------------
Libraries compiled on Fri Oct 15 13:01:23 CDT 2010 on atreides 
Machine characteristics: Linux atreides 2.6.32-25-generic #44-Ubuntu SMP Fri Sep 17 20:05:27 UTC 2010 x86_64 GNU/Linux 
Using PETSc directory: /org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5
Using PETSc arch: gcc-4.5-lucid-mpich2-1.2.1-cxx-opt
-----------------------------------------
Using C compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpicxx -Wall -Wwrite-strings -Wno-strict-aliasing -O3   -fPIC   
Using Fortran compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpif90 -fPIC -Wall -Wno-unused-variable -O3    
-----------------------------------------
Using include paths: -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/include  
------------------------------------------
Using C linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpicxx -Wall -Wwrite-strings -Wno-strict-aliasing -O3 
Using Fortran linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpif90 -fPIC -Wall -Wno-unused-variable -O3  
Using libraries: -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -lpetsc       -lX11 -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -lHYPRE -lsuperlu_dist_2.4 -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -lscalapack -lblacs -lsuperlu_4.0 -Wl,-rpath,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t -L/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t -lmkl_solver_lp64_sequential -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -Wl,-rpath,/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/lib -L/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/lib -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib/gcc/x86_64-unknown-linux-gnu/4.5.1 -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib/gcc/x86_64-unknown-linux-gnu/4.5.1 -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib64 -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib64 -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib -ldl -lmpich -lopa -lpthread -lrt -lgcc_s -lmpichf90 -lgfortran -lm -lm -lmpichcxx -lstdc++ -ldl -lmpich -lopa -lpthread -lrt -lgcc_s -ldl  
------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.017591, Active time=0.00747                                              |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 1         0.0000      0.000016    0.0000      0.000018    0.21     0.24     |
|   compute_sparsity()           1         0.0001      0.000075    0.0001      0.000090    1.00     1.20     |
|   create_dof_constraints()     1         0.0000      0.000003    0.0000      0.000003    0.04     0.04     |
|   distribute_dofs()            1         0.0001      0.000057    0.0001      0.000141    0.76     1.89     |
|   dof_indices()                99        0.0000      0.000000    0.0000      0.000000    0.47     0.47     |
|   prepare_send_list()          1         0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|   reinit()                     1         0.0001      0.000056    0.0001      0.000056    0.75     0.75     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         18        0.0000      0.000001    0.0000      0.000001    0.15     0.15     |
|   compute_face_map()           1         0.0000      0.000003    0.0000      0.000003    0.04     0.04     |
|   compute_shape_functions()    18        0.0000      0.000000    0.0000      0.000000    0.08     0.08     |
|   init_face_shape_functions()  1         0.0000      0.000004    0.0000      0.000004    0.05     0.05     |
|   init_shape_functions()       2         0.0000      0.000020    0.0000      0.000020    0.54     0.54     |
|   inverse_map()                1         0.0000      0.000004    0.0000      0.000004    0.05     0.05     |
|                                                                                                            |
| LocationMap                                                                                                |
|   find()                       4         0.0000      0.000001    0.0000      0.000001    0.05     0.05     |
|   init()                       1         0.0000      0.000018    0.0000      0.000018    0.24     0.24     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             2         0.0001      0.000071    0.0002      0.000080    1.90     2.16     |
|   renumber_nodes_and_elem()    4         0.0000      0.000001    0.0000      0.000001    0.04     0.04     |
|                                                                                                            |
| MeshCommunication                                                                                          |
|   compute_hilbert_indices()    3         0.0003      0.000098    0.0003      0.000098    3.95     3.95     |
|   find_global_indices()        3         0.0001      0.000034    0.0008      0.000254    1.35     10.19    |
|   parallel_sort()              3         0.0002      0.000076    0.0003      0.000088    3.07     3.52     |
|                                                                                                            |
| MeshRefinement                                                                                             |
|   _refine_elements()           1         0.0000      0.000031    0.0001      0.000053    0.41     0.71     |
|   add_point()                  4         0.0000      0.000003    0.0000      0.000004    0.15     0.20     |
|   make_refinement_compatible() 1         0.0000      0.000010    0.0000      0.000016    0.13     0.21     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0000      0.000037    0.0000      0.000037    0.50     0.50     |
|                                                                                                            |
| MetisPartitioner                                                                                           |
|   partition()                  2         0.0002      0.000120    0.0006      0.000276    3.23     7.39     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  8         0.0000      0.000005    0.0000      0.000005    0.58     0.58     |
|   broadcast()                  2         0.0000      0.000002    0.0000      0.000002    0.07     0.07     |
|   gather()                     2         0.0000      0.000002    0.0000      0.000002    0.07     0.07     |
|   max(bool)                    2         0.0000      0.000005    0.0000      0.000005    0.15     0.15     |
|   max(scalar)                  3         0.0000      0.000008    0.0000      0.000008    0.31     0.31     |
|   max(vector)                  3         0.0000      0.000002    0.0000      0.000002    0.07     0.07     |
|   min(bool)                    1         0.0000      0.000006    0.0000      0.000006    0.08     0.08     |
|   min(vector)                  3         0.0000      0.000004    0.0000      0.000004    0.16     0.16     |
|   probe()                      15        0.0000      0.000002    0.0000      0.000002    0.31     0.31     |
|   receive()                    15        0.0000      0.000003    0.0001      0.000004    0.51     0.83     |
|   send()                       15        0.0000      0.000001    0.0000      0.000001    0.23     0.23     |
|   send_receive()               20        0.0000      0.000002    0.0001      0.000006    0.43     1.61     |
|   sum()                        12        0.0000      0.000004    0.0000      0.000004    0.59     0.59     |
|   wait()                       14        0.0000      0.000001    0.0000      0.000001    0.16     0.16     |
|                                                                                                            |
| Partitioner                                                                                                |
|   set_node_processor_ids()     2         0.0000      0.000023    0.0001      0.000035    0.63     0.94     |
|   set_parent_processor_ids()   2         0.0000      0.000003    0.0000      0.000003    0.08     0.08     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      1         0.0053      0.005282    0.0053      0.005282    70.71    70.71    |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   1         0.0004      0.000426    0.0005      0.000505    5.70     6.76     |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        296       0.0075                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

Running ./ex25-opt -d 2 -n 30 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1329
    n_local_nodes()=683
  n_elem()=1268
    n_local_elem()=633
    n_active_elem()=1176
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Poisson"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=1329
    n_local_dofs()=683
    n_constrained_dofs()=184
    n_vectors()=1

 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=1329
    n_local_nodes()=683
  n_elem()=1268
    n_local_elem()=633
    n_active_elem()=1176
  n_subdomains()=2
  n_processors()=2
  processor_id()=0


-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 2                                                |
| Time:           Thu Feb  3 12:13:06 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.009794, Active time=0.004693                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           318       0.0018      0.000006    0.0018      0.000006    37.72    37.72    |
| Fe                            318       0.0012      0.000004    0.0012      0.000004    25.46    25.46    |
| Ke                            318       0.0002      0.000001    0.0002      0.000001    3.77     3.77     |
| elem init                     318       0.0013      0.000004    0.0013      0.000004    26.66    26.66    |
| matrix insertion              318       0.0003      0.000001    0.0003      0.000001    6.39     6.39     |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       1590      0.0047                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex25-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Thu Feb  3 12:13:06 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           4.748e-02      1.03204   4.675e-02
Objects:              6.800e+01      1.00000   6.800e+01
Flops:                9.821e+05      1.04903   9.592e+05  1.918e+06
Flops/sec:            2.068e+07      1.01646   2.052e+07  4.103e+07
MPI Messages:         3.050e+01      1.00000   3.050e+01  6.100e+01
MPI Message Lengths:  1.668e+04      1.03679   5.373e+02  3.278e+04
MPI Reductions:       1.410e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 4.6718e-02  99.9%  1.9183e+06 100.0%  6.100e+01 100.0%  5.373e+02      100.0%  9.900e+01  70.2% 

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

VecMDot               12 1.0 7.9870e-05 1.3 6.28e+04 1.1 0.0e+00 0.0e+00 1.2e+01  0  6  0  0  9   0  6  0  0 12  1500
VecNorm               14 1.0 7.8917e-05 1.1 1.13e+04 1.1 0.0e+00 0.0e+00 1.4e+01  0  1  0  0 10   0  1  0  0 14   273
VecScale              13 1.0 1.3351e-05 1.0 5.24e+03 1.1 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0   749
VecCopy                3 1.0 2.8610e-06 3.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                19 1.0 9.7752e-06 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                2 1.0 1.8120e-05 1.2 1.61e+03 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   170
VecMAXPY              13 1.0 2.6703e-05 1.3 7.25e+04 1.1 0.0e+00 0.0e+00 0.0e+00  0  7  0  0  0   0  7  0  0  0  5184
VecAssemblyBegin       3 1.0 4.5061e-05 1.1 0.00e+00 0.0 2.0e+00 4.0e+02 9.0e+00  0  0  3  2  6   0  0  3  2  9     0
VecAssemblyEnd         3 1.0 1.0967e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       18 1.0 4.3392e-05 1.0 0.00e+00 0.0 2.8e+01 3.0e+02 0.0e+00  0  0 46 25  0   0  0 46 25  0     0
VecScatterEnd         18 1.0 9.3460e-05 4.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize          13 1.0 1.0705e-04 1.0 1.57e+04 1.1 0.0e+00 0.0e+00 1.3e+01  0  2  0  0  9   0  2  0  0 13   280
MatMult               13 1.0 2.1338e-04 1.5 8.70e+04 1.1 2.6e+01 2.8e+02 0.0e+00  0  9 43 22  0   0  9 43 22  0   775
MatSolve              13 1.0 2.6226e-04 1.1 4.05e+05 1.0 0.0e+00 0.0e+00 0.0e+00  1 41  0  0  0   1 41  0  0  0  3028
MatLUFactorNum         1 1.0 4.0102e-04 1.1 3.37e+05 1.0 0.0e+00 0.0e+00 0.0e+00  1 35  0  0  0   1 35  0  0  0  1661
MatILUFactorSym        1 1.0 9.8801e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+00  2  0  0  0  1   2  0  0  0  1     0
MatAssemblyBegin       6 1.0 1.3018e-04 1.1 0.00e+00 0.0 3.0e+00 1.5e+03 8.0e+00  0  0  5 14  6   0  0  5 14  8     0
MatAssemblyEnd         6 1.0 4.1699e-04 1.0 0.00e+00 0.0 1.2e+01 7.3e+01 2.2e+01  1  0 20  3 16   1  0 20  3 22     0
MatGetRowIJ            1 1.0 9.5367e-07 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetSubMatrice       2 1.0 2.8205e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+01  1  0  0  0  7   1  0  0  0 10     0
MatGetOrdering         1 1.0 2.9087e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  1   0  0  0  0  2     0
MatZeroEntries         3 1.0 2.2888e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        12 1.0 1.1873e-04 1.1 1.26e+05 1.1 0.0e+00 0.0e+00 1.2e+01  0 12  0  0  9   0 12  0  0 12  2019
KSPSetup               2 1.0 5.1975e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 2.5868e-03 1.0 9.82e+05 1.0 2.6e+01 2.8e+02 2.9e+01  6100 43 22 21   6100 43 22 29   742
PCSetUp                2 1.0 1.6260e-03 1.1 3.37e+05 1.0 0.0e+00 0.0e+00 3.0e+00  3 35  0  0  2   3 35  0  0  3   410
PCSetUpOnBlocks        1 1.0 1.4720e-03 1.1 3.37e+05 1.0 0.0e+00 0.0e+00 3.0e+00  3 35  0  0  2   3 35  0  0  3   452
PCApply               13 1.0 3.8505e-04 1.1 4.05e+05 1.0 0.0e+00 0.0e+00 0.0e+00  1 41  0  0  0   1 41  0  0  0  2062
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    29             29       111056     0
         Vec Scatter     6              6         5208     0
           Index Set    16             16        19924     0
   IS L to G Mapping     1              1         3312     0
              Matrix    12             12       522444     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 1.38283e-06
Average time for zero size MPI_Send(): 4.52995e-06
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
Configure run at: Fri Oct 15 13:01:23 2010
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared=1 --with-mpi-dir=/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid --with-mumps=true --download-mumps=ifneeded --with-parmetis=true --download-parmetis=ifneeded --with-superlu=true --download-superlu=ifneeded --with-superludir=true --download-superlu_dist=ifneeded --with-blacs=true --download-blacs=ifneeded --with-scalapack=true --download-scalapack=ifneeded --with-hypre=true --download-hypre=ifneeded --with-blas-lib="[/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_intel_lp64.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_sequential.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_core.so]" --with-lapack-lib=/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_solver_lp64_sequential.a
-----------------------------------------
Libraries compiled on Fri Oct 15 13:01:23 CDT 2010 on atreides 
Machine characteristics: Linux atreides 2.6.32-25-generic #44-Ubuntu SMP Fri Sep 17 20:05:27 UTC 2010 x86_64 GNU/Linux 
Using PETSc directory: /org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5
Using PETSc arch: gcc-4.5-lucid-mpich2-1.2.1-cxx-opt
-----------------------------------------
Using C compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpicxx -Wall -Wwrite-strings -Wno-strict-aliasing -O3   -fPIC   
Using Fortran compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpif90 -fPIC -Wall -Wno-unused-variable -O3    
-----------------------------------------
Using include paths: -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/include  
------------------------------------------
Using C linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpicxx -Wall -Wwrite-strings -Wno-strict-aliasing -O3 
Using Fortran linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpif90 -fPIC -Wall -Wno-unused-variable -O3  
Using libraries: -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -lpetsc       -lX11 -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -lHYPRE -lsuperlu_dist_2.4 -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -lscalapack -lblacs -lsuperlu_4.0 -Wl,-rpath,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t -L/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t -lmkl_solver_lp64_sequential -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -Wl,-rpath,/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/lib -L/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/lib -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib/gcc/x86_64-unknown-linux-gnu/4.5.1 -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib/gcc/x86_64-unknown-linux-gnu/4.5.1 -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib64 -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib64 -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib -ldl -lmpich -lopa -lpthread -lrt -lgcc_s -lmpichf90 -lgfortran -lm -lm -lmpichcxx -lstdc++ -ldl -lmpich -lopa -lpthread -lrt -lgcc_s -ldl  
------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.054996, Active time=0.039504                                             |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 1         0.0003      0.000251    0.0003      0.000280    0.64     0.71     |
|   build_constraint_matrix()    318       0.0003      0.000001    0.0003      0.000001    0.67     0.67     |
|   cnstrn_elem_mat_vec()        318       0.0036      0.000011    0.0036      0.000011    9.13     9.13     |
|   compute_sparsity()           1         0.0017      0.001700    0.0021      0.002109    4.30     5.34     |
|   create_dof_constraints()     1         0.0012      0.001246    0.0017      0.001737    3.15     4.40     |
|   distribute_dofs()            1         0.0007      0.000727    0.0017      0.001738    1.84     4.40     |
|   dof_indices()                3802      0.0010      0.000000    0.0010      0.000000    2.61     2.61     |
|   prepare_send_list()          1         0.0000      0.000008    0.0000      0.000008    0.02     0.02     |
|   reinit()                     1         0.0010      0.000954    0.0010      0.000954    2.41     2.41     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         415       0.0003      0.000001    0.0003      0.000001    0.74     0.74     |
|   compute_face_map()           97        0.0003      0.000003    0.0006      0.000006    0.74     1.53     |
|   compute_shape_functions()    415       0.0001      0.000000    0.0001      0.000000    0.35     0.35     |
|   init_face_shape_functions()  62        0.0000      0.000001    0.0000      0.000001    0.10     0.10     |
|   init_shape_functions()       98        0.0002      0.000002    0.0002      0.000002    0.56     0.56     |
|   inverse_map()                1318      0.0008      0.000001    0.0008      0.000001    1.94     1.94     |
|                                                                                                            |
| GMVIO                                                                                                      |
|   write_nodal_data()           1         0.0035      0.003505    0.0035      0.003505    8.87     8.87     |
|                                                                                                            |
| LocationMap                                                                                                |
|   find()                       1104      0.0004      0.000000    0.0004      0.000000    1.13     1.13     |
|   init()                       1         0.0002      0.000185    0.0002      0.000185    0.47     0.47     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             2         0.0029      0.001430    0.0029      0.001440    7.24     7.29     |
|   renumber_nodes_and_elem()    4         0.0001      0.000028    0.0001      0.000028    0.29     0.29     |
|                                                                                                            |
| MeshCommunication                                                                                          |
|   compute_hilbert_indices()    3         0.0053      0.001772    0.0053      0.001772    13.46    13.46    |
|   find_global_indices()        3         0.0008      0.000258    0.0066      0.002204    1.96     16.74    |
|   parallel_sort()              3         0.0004      0.000120    0.0004      0.000132    0.91     1.00     |
|                                                                                                            |
| MeshRefinement                                                                                             |
|   _refine_elements()           1         0.0010      0.001019    0.0023      0.002334    2.58     5.91     |
|   add_point()                  1104      0.0007      0.000001    0.0012      0.000001    1.78     3.11     |
|   make_refinement_compatible() 1         0.0000      0.000050    0.0001      0.000059    0.13     0.15     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0003      0.000319    0.0003      0.000319    0.81     0.81     |
|                                                                                                            |
| MetisPartitioner                                                                                           |
|   partition()                  2         0.0025      0.001259    0.0069      0.003458    6.37     17.51    |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  8         0.0001      0.000008    0.0001      0.000008    0.16     0.16     |
|   broadcast()                  2         0.0000      0.000003    0.0000      0.000003    0.02     0.02     |
|   gather()                     2         0.0000      0.000002    0.0000      0.000002    0.01     0.01     |
|   max(bool)                    2         0.0000      0.000004    0.0000      0.000004    0.02     0.02     |
|   max(scalar)                  3         0.0001      0.000027    0.0001      0.000027    0.21     0.21     |
|   max(vector)                  3         0.0000      0.000002    0.0000      0.000002    0.02     0.02     |
|   min(bool)                    1         0.0000      0.000009    0.0000      0.000009    0.02     0.02     |
|   min(vector)                  3         0.0000      0.000006    0.0000      0.000006    0.05     0.05     |
|   probe()                      15        0.0001      0.000005    0.0001      0.000005    0.18     0.18     |
|   receive()                    15        0.0000      0.000003    0.0001      0.000008    0.11     0.29     |
|   send()                       15        0.0000      0.000001    0.0000      0.000001    0.06     0.06     |
|   send_receive()               20        0.0000      0.000002    0.0002      0.000009    0.09     0.47     |
|   sum()                        15        0.0001      0.000008    0.0001      0.000008    0.32     0.32     |
|   wait()                       14        0.0000      0.000001    0.0000      0.000001    0.03     0.03     |
|                                                                                                            |
| Partitioner                                                                                                |
|   set_node_processor_ids()     2         0.0007      0.000330    0.0007      0.000357    1.67     1.81     |
|   set_parent_processor_ids()   2         0.0001      0.000055    0.0001      0.000055    0.28     0.28     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      1         0.0043      0.004286    0.0043      0.004286    10.85    10.85    |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   1         0.0042      0.004229    0.0099      0.009911    10.71    25.09    |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        9203      0.0395                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

Running ./ex25-opt -d 3 -n 15 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=10854
    n_local_nodes()=5663
  n_elem()=8767
    n_local_elem()=4383
    n_active_elem()=8093
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Poisson"
    Type "LinearImplicit"
    Variables="u" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" 
    Approximation Orders="FIRST", "THIRD" 
    n_dofs()=10854
    n_local_dofs()=5663
    n_constrained_dofs()=4068
    n_vectors()=1

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=10854
    n_local_nodes()=5663
  n_elem()=8767
    n_local_elem()=4383
    n_active_elem()=8093
  n_subdomains()=2
  n_processors()=2
  processor_id()=0


-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 2                                                |
| Time:           Thu Feb  3 12:13:06 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 -----------------------------------------------------------------------------------------------------------
| Matrix Assembly Performance: Alive time=0.137466, Active time=0.114544                                    |
 -----------------------------------------------------------------------------------------------------------
| Event                         nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                         w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|-----------------------------------------------------------------------------------------------------------|
|                                                                                                           |
| BCs                           1470      0.0676      0.000046    0.0676      0.000046    59.00    59.00    |
| Fe                            1470      0.0198      0.000013    0.0198      0.000013    17.25    17.25    |
| Ke                            1470      0.0072      0.000005    0.0072      0.000005    6.29     6.29     |
| elem init                     1470      0.0133      0.000009    0.0133      0.000009    11.59    11.59    |
| matrix insertion              1470      0.0067      0.000005    0.0067      0.000005    5.87     5.87     |
 -----------------------------------------------------------------------------------------------------------
| Totals:                       7350      0.1145                                          100.00            |
 -----------------------------------------------------------------------------------------------------------

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex25-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Thu Feb  3 12:13:08 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           2.007e+00      1.00365   2.004e+00
Objects:              7.800e+01      1.00000   7.800e+01
Flops:                8.973e+08      1.25174   8.071e+08  1.614e+09
Flops/sec:            4.471e+08      1.24720   4.028e+08  8.055e+08
MPI Messages:         3.450e+01      1.00000   3.450e+01  6.900e+01
MPI Message Lengths:  2.286e+05      1.02703   6.538e+03  4.511e+05
MPI Reductions:       1.490e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.0035e+00 100.0%  1.6142e+09 100.0%  6.900e+01 100.0%  6.538e+03      100.0%  1.070e+02  71.8% 

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

VecMDot               16 1.0 2.1791e-03 4.6 7.34e+05 1.1 0.0e+00 0.0e+00 1.6e+01  0  0  0  0 11   0  0  0  0 15   638
VecNorm               18 1.0 4.6539e-04 1.3 9.71e+04 1.1 0.0e+00 0.0e+00 1.8e+01  0  0  0  0 12   0  0  0  0 17   396
VecScale              17 1.0 3.9577e-05 1.0 4.59e+04 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2198
VecCopy                3 1.0 2.3127e-05 2.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                23 1.0 4.1723e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                2 1.0 4.3869e-05 1.3 1.08e+04 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   466
VecMAXPY              17 1.0 2.1172e-04 1.1 8.20e+05 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  7346
VecAssemblyBegin       3 1.0 7.6771e-05 1.1 0.00e+00 0.0 2.0e+00 6.4e+03 9.0e+00  0  0  3  3  6   0  0  3  3  8     0
VecAssemblyEnd         3 1.0 1.8120e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       22 1.0 2.2626e-04 1.5 0.00e+00 0.0 3.6e+01 3.1e+03 0.0e+00  0  0 52 25  0   0  0 52 25  0     0
VecScatterEnd         22 1.0 2.5539e-013867.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  6  0  0  0  0   6  0  0  0  0     0
VecNormalize          17 1.0 5.4026e-04 1.2 1.38e+05 1.1 0.0e+00 0.0e+00 1.7e+01  0  0  0  0 11   0  0  0  0 16   483
MatMult               17 1.0 2.5742e-01114.9 2.34e+06 1.1 3.4e+01 2.9e+03 0.0e+00  6  0 49 22  0   6  0 49 22  0    17
MatSolve              17 1.0 6.6412e-02 1.1 7.64e+07 1.2 0.0e+00 0.0e+00 0.0e+00  3  9  0  0  0   3  9  0  0  0  2137
MatLUFactorNum         1 1.0 6.2586e-01 1.3 8.17e+08 1.3 0.0e+00 0.0e+00 0.0e+00 28 91  0  0  0  28 91  0  0  0  2340
MatILUFactorSym        1 1.0 8.2042e-01 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+00 38  0  0  0  1  38  0  0  0  1     0
MatAssemblyBegin       6 1.0 6.0463e-04 3.7 0.00e+00 0.0 3.0e+00 5.7e+04 8.0e+00  0  0  4 38  5   0  0  4 38  7     0
MatAssemblyEnd         6 1.0 2.3627e-03 1.0 0.00e+00 0.0 1.2e+01 7.3e+02 2.2e+01  0  0 17  2 15   0  0 17  2 21     0
MatGetRowIJ            1 1.0 0.0000e+00 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetSubMatrice       2 1.0 2.1639e-03 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+01  0  0  0  0  7   0  0  0  0  9     0
MatGetOrdering         1 1.0 5.0068e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  1   0  0  0  0  2     0
MatZeroEntries         3 1.0 7.1788e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        16 1.0 2.3963e-03 3.4 1.47e+06 1.1 0.0e+00 0.0e+00 1.6e+01  0  0  0  0 11   0  0  0  0 15  1161
KSPSetup               2 1.0 8.4877e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 1.5173e+00 1.0 8.97e+08 1.3 3.4e+01 2.9e+03 3.7e+01 76100 49 22 25  76100 49 22 35  1064
PCSetUp                2 1.0 1.4466e+00 1.2 8.17e+08 1.3 0.0e+00 0.0e+00 3.0e+00 66 91  0  0  2  66 91  0  0  3  1012
PCSetUpOnBlocks        1 1.0 1.4464e+00 1.2 8.17e+08 1.3 0.0e+00 0.0e+00 3.0e+00 66 91  0  0  2  66 91  0  0  3  1013
PCApply               17 1.0 6.6727e-02 1.1 7.64e+07 1.2 0.0e+00 0.0e+00 0.0e+00  3  9  0  0  0   3  9  0  0  0  2127
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    39             39       786280     0
         Vec Scatter     6              6         5208     0
           Index Set    16             16        89392     0
   IS L to G Mapping     1              1        25360     0
              Matrix    12             12     32861396     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 1.38283e-06
Average time for zero size MPI_Send(): 6.55651e-06
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
Configure run at: Fri Oct 15 13:01:23 2010
Configure options: --with-debugging=false --COPTFLAGS=-O3 --CXXOPTFLAGS=-O3 --FOPTFLAGS=-O3 --with-clanguage=C++ --with-shared=1 --with-mpi-dir=/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid --with-mumps=true --download-mumps=ifneeded --with-parmetis=true --download-parmetis=ifneeded --with-superlu=true --download-superlu=ifneeded --with-superludir=true --download-superlu_dist=ifneeded --with-blacs=true --download-blacs=ifneeded --with-scalapack=true --download-scalapack=ifneeded --with-hypre=true --download-hypre=ifneeded --with-blas-lib="[/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_intel_lp64.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_sequential.so,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_core.so]" --with-lapack-lib=/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t/libmkl_solver_lp64_sequential.a
-----------------------------------------
Libraries compiled on Fri Oct 15 13:01:23 CDT 2010 on atreides 
Machine characteristics: Linux atreides 2.6.32-25-generic #44-Ubuntu SMP Fri Sep 17 20:05:27 UTC 2010 x86_64 GNU/Linux 
Using PETSc directory: /org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5
Using PETSc arch: gcc-4.5-lucid-mpich2-1.2.1-cxx-opt
-----------------------------------------
Using C compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpicxx -Wall -Wwrite-strings -Wno-strict-aliasing -O3   -fPIC   
Using Fortran compiler: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpif90 -fPIC -Wall -Wno-unused-variable -O3    
-----------------------------------------
Using include paths: -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/include -I/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/include -I/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/include  
------------------------------------------
Using C linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpicxx -Wall -Wwrite-strings -Wno-strict-aliasing -O3 
Using Fortran linker: /org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/bin/mpif90 -fPIC -Wall -Wno-unused-variable -O3  
Using libraries: -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -lpetsc       -lX11 -Wl,-rpath,/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -L/org/centers/pecos/LIBRARIES/PETSC3/petsc-3.1-p5/gcc-4.5-lucid-mpich2-1.2.1-cxx-opt/lib -lHYPRE -lsuperlu_dist_2.4 -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lparmetis -lmetis -lscalapack -lblacs -lsuperlu_4.0 -Wl,-rpath,/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t -L/org/centers/pecos/LIBRARIES/MKL/mkl-10.0.3.020-gcc-4.5-lucid/lib/em64t -lmkl_solver_lp64_sequential -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm -Wl,-rpath,/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/lib -L/org/centers/pecos/LIBRARIES/MPICH2/mpich2-1.2.1-gcc-4.5-lucid/lib -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib/gcc/x86_64-unknown-linux-gnu/4.5.1 -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib/gcc/x86_64-unknown-linux-gnu/4.5.1 -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib64 -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib64 -Wl,-rpath,/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib -L/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/lib -ldl -lmpich -lopa -lpthread -lrt -lgcc_s -lmpichf90 -lgfortran -lm -lm -lmpichcxx -lstdc++ -ldl -lmpich -lopa -lpthread -lrt -lgcc_s -ldl  
------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=2.01469, Active time=1.95602                                               |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 1         0.0028      0.002772    0.0031      0.003125    0.14     0.16     |
|   build_constraint_matrix()    1470      0.0057      0.000004    0.0057      0.000004    0.29     0.29     |
|   cnstrn_elem_mat_vec()        1470      0.0108      0.000007    0.0108      0.000007    0.55     0.55     |
|   compute_sparsity()           1         0.0372      0.037153    0.0415      0.041482    1.90     2.12     |
|   create_dof_constraints()     1         0.0292      0.029200    0.0522      0.052223    1.49     2.67     |
|   distribute_dofs()            1         0.0067      0.006744    0.0163      0.016315    0.34     0.83     |
|   dof_indices()                30763     0.0126      0.000000    0.0126      0.000000    0.64     0.64     |
|   prepare_send_list()          1         0.0001      0.000137    0.0001      0.000137    0.01     0.01     |
|   reinit()                     1         0.0095      0.009459    0.0095      0.009459    0.48     0.48     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         2829      0.0063      0.000002    0.0063      0.000002    0.32     0.32     |
|   compute_face_map()           1359      0.0027      0.000002    0.0027      0.000002    0.14     0.14     |
|   compute_shape_functions()    2829      0.0034      0.000001    0.0034      0.000001    0.18     0.18     |
|   init_face_shape_functions()  1111      0.0048      0.000004    0.0048      0.000004    0.24     0.24     |
|   init_shape_functions()       1360      0.0101      0.000007    0.0101      0.000007    0.52     0.52     |
|   inverse_map()                33927     0.0498      0.000001    0.0498      0.000001    2.55     2.55     |
|                                                                                                            |
| GMVIO                                                                                                      |
|   write_nodal_data()           1         0.0304      0.030432    0.0304      0.030432    1.56     1.56     |
|                                                                                                            |
| LocationMap                                                                                                |
|   find()                       37744     0.0117      0.000000    0.0117      0.000000    0.60     0.60     |
|   init()                       1         0.0007      0.000734    0.0007      0.000734    0.04     0.04     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             2         0.0251      0.012527    0.0251      0.012555    1.28     1.28     |
|   renumber_nodes_and_elem()    4         0.0015      0.000365    0.0015      0.000365    0.07     0.07     |
|                                                                                                            |
| MeshCommunication                                                                                          |
|   compute_hilbert_indices()    3         0.0264      0.008791    0.0264      0.008791    1.35     1.35     |
|   find_global_indices()        3         0.0036      0.001194    0.0313      0.010436    0.18     1.60     |
|   parallel_sort()              3         0.0009      0.000295    0.0010      0.000327    0.05     0.05     |
|                                                                                                            |
| MeshRefinement                                                                                             |
|   _refine_elements()           1         0.0294      0.029374    0.0690      0.068964    1.50     3.53     |
|   add_point()                  37744     0.0227      0.000001    0.0371      0.000001    1.16     1.90     |
|   make_refinement_compatible() 1         0.0002      0.000225    0.0002      0.000231    0.01     0.01     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0015      0.001481    0.0015      0.001481    0.08     0.08     |
|                                                                                                            |
| MetisPartitioner                                                                                           |
|   partition()                  2         0.0173      0.008658    0.0412      0.020609    0.89     2.11     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  8         0.0002      0.000029    0.0002      0.000029    0.01     0.01     |
|   broadcast()                  2         0.0000      0.000005    0.0000      0.000005    0.00     0.00     |
|   gather()                     2         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   max(bool)                    2         0.0000      0.000005    0.0000      0.000005    0.00     0.00     |
|   max(scalar)                  3         0.0018      0.000604    0.0018      0.000604    0.09     0.09     |
|   max(vector)                  3         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   min(bool)                    1         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   min(vector)                  3         0.0000      0.000007    0.0000      0.000007    0.00     0.00     |
|   probe()                      15        0.0002      0.000016    0.0002      0.000016    0.01     0.01     |
|   receive()                    15        0.0001      0.000006    0.0003      0.000022    0.00     0.02     |
|   send()                       15        0.0001      0.000004    0.0001      0.000004    0.00     0.00     |
|   send_receive()               20        0.0000      0.000002    0.0005      0.000023    0.00     0.02     |
|   sum()                        15        0.0007      0.000047    0.0007      0.000047    0.04     0.04     |
|   wait()                       14        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                            |
| Partitioner                                                                                                |
|   set_node_processor_ids()     2         0.0043      0.002146    0.0045      0.002242    0.22     0.23     |
|   set_parent_processor_ids()   2         0.0007      0.000363    0.0007      0.000363    0.04     0.04     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      1         1.5264      1.526411    1.5264      1.526411    78.04    78.04    |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   1         0.0584      0.058394    0.1376      0.137609    2.99     7.04     |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        152758    1.9560                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 2 ./ex25-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
