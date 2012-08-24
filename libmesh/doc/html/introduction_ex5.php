<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("introduction_ex5",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Introduction Example 5 - Run-Time Quadrature Rule Selection</h1>

<br><br>This is the fifth example program.  It builds on
the previous two examples, and extends the use
of the \p AutoPtr as a convenient build method to
determine the quadrature rule at run time.


<br><br>

<br><br>C++ include files that we need
</div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
        #include &lt;sstream&gt; 
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
Define the base quadrature class, with which
specialized quadrature rules will be built.
</div>

<div class ="fragment">
<pre>
        #include "quadrature.h"
        
</pre>
</div>
<div class = "comment">
Include the namespace \p QuadratureRules for
some handy descriptions.
</div>

<div class ="fragment">
<pre>
        #include "quadrature_rules.h"
        
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
To impose Dirichlet boundary conditions
</div>

<div class ="fragment">
<pre>
        #include "dirichlet_boundaries.h"
        #include "analytic_function.h"
        
</pre>
</div>
<div class = "comment">
The definition of a geometric element
</div>

<div class ="fragment">
<pre>
        #include "elem.h"
        
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
Function prototype, as before.
</div>

<div class ="fragment">
<pre>
        void assemble_poisson(EquationSystems& es,
                              const std::string& system_name);
        
        
        
</pre>
</div>
<div class = "comment">
Exact solution function prototype, as before.
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
          output(0) = exact_solution(p(0),p(1),p(2));
        }
        
        
</pre>
</div>
<div class = "comment">
The quadrature type the user requests.
</div>

<div class ="fragment">
<pre>
        QuadratureType quad_type=INVALID_Q_RULE;
        
        
        
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
Check for proper usage.  The quadrature rule
must be given at run time.
</div>

<div class ="fragment">
<pre>
          if (argc &lt; 3)
            {
              if (libMesh::processor_id() == 0)
                {
                  std::cerr &lt;&lt; "Usage: " &lt;&lt; argv[0] &lt;&lt; " -q n"
                            &lt;&lt; std::endl;
                  std::cerr &lt;&lt; "  where n stands for:" &lt;&lt; std::endl;
        
              
</pre>
</div>
<div class = "comment">
Note that only some of all quadrature rules are
valid choices.  For example, the Jacobi quadrature
is actually a "helper" for higher-order rules,
included in QGauss.
</div>

<div class ="fragment">
<pre>
                  for (unsigned int n=0; n&lt;QuadratureRules::num_valid_elem_rules; n++)
                    std::cerr &lt;&lt; "  " &lt;&lt; QuadratureRules::valid_elem_rules[n] &lt;&lt; "    " 
                              &lt;&lt; QuadratureRules::name(QuadratureRules::valid_elem_rules[n])
                              &lt;&lt; std::endl;
              
                  std::cerr &lt;&lt; std::endl;
                }
              
              libmesh_error();
            }
          
          
</pre>
</div>
<div class = "comment">
Tell the user what we are doing.
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
Set the quadrature rule type that the user wants from argv[2]
</div>

<div class ="fragment">
<pre>
          quad_type = static_cast&lt;QuadratureType&gt;(std::atoi(argv[2]));
        
</pre>
</div>
<div class = "comment">
Skip this 3D example if libMesh was compiled as 1D-only.
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(3 &lt;= LIBMESH_DIM, "3D support");
          
</pre>
</div>
<div class = "comment">
The following is identical to example 4, and therefore
not commented.  Differences are mentioned when present.
</div>

<div class ="fragment">
<pre>
          Mesh mesh;
        
</pre>
</div>
<div class = "comment">
We will use a linear approximation space in this example,
hence 8-noded hexahedral elements are sufficient.  This
is different than example 4 where we used 27-noded
hexahedral elements to support a second-order approximation
space.
</div>

<div class ="fragment">
<pre>
          MeshTools::Generation::build_cube (mesh,
                                             16, 16, 16,
                                             -1., 1.,
                                             -1., 1.,
                                             -1., 1.,
                                             HEX8);
          
          mesh.print_info();
          
          EquationSystems equation_systems (mesh);
          
          equation_systems.add_system&lt;LinearImplicitSystem&gt; ("Poisson");
          
          unsigned int u_var = equation_systems.get_system("Poisson").add_variable("u", FIRST);
        
          equation_systems.get_system("Poisson").attach_assemble_function (assemble_poisson);
        
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
          boundary_ids.insert(2);
          boundary_ids.insert(3);
          boundary_ids.insert(4);
          boundary_ids.insert(5);
        
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
          equation_systems.get_system("Poisson").get_dof_map().add_dirichlet_boundary(dirichlet_bc);
        
          equation_systems.init();
          
          equation_systems.print_info();
        
          equation_systems.get_system("Poisson").solve();
        
</pre>
</div>
<div class = "comment">
"Personalize" the output, with the
number of the quadrature rule appended.
</div>

<div class ="fragment">
<pre>
          std::ostringstream f_name;
          f_name &lt;&lt; "out_" &lt;&lt; quad_type &lt;&lt; ".e";
        
        #ifdef LIBMESH_HAVE_EXODUS_API
          ExodusII_IO(mesh).write_equation_systems (f_name.str(),
                                              equation_systems);
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
        
        
        
        
        void assemble_poisson(EquationSystems& es,
                              const std::string& system_name)
        {
          libmesh_assert (system_name == "Poisson");
        
          const MeshBase& mesh = es.get_mesh();
        
          const unsigned int dim = mesh.mesh_dimension();
        
          LinearImplicitSystem& system = es.get_system&lt;LinearImplicitSystem&gt;("Poisson");
          
          const DofMap& dof_map = system.get_dof_map();
          
          FEType fe_type = dof_map.variable_type(0);
        
          
</pre>
</div>
<div class = "comment">
Build a Finite Element object of the specified type.  Since the
\p FEBase::build() member dynamically creates memory we will
store the object as an \p AutoPtr<FEBase>.  Below, the
functionality of \p AutoPtr's is described more detailed in 
the context of building quadrature rules.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
          
          
</pre>
</div>
<div class = "comment">
Now this deviates from example 4.  we create a 
5th order quadrature rule of user-specified type
for numerical integration.  Note that not all
quadrature rules support this order.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;QBase&gt; qrule(QBase::build(quad_type, dim, THIRD));
        
        
          
</pre>
</div>
<div class = "comment">
Tell the finte element object to use our
quadrature rule.  Note that a \p AutoPtr<QBase> returns
a QBase* pointer to the object it handles with \p get().  
However, using \p get(), the \p AutoPtr<QBase> \p qrule is 
still in charge of this pointer. I.e., when \p qrule goes 
out of scope, it will safely delete the \p QBase object it 
points to.  This behavior may be overridden using
\p AutoPtr<Xyz>::release(), but is currently not
recommended.
</div>

<div class ="fragment">
<pre>
          fe-&gt;attach_quadrature_rule (qrule.get());
        
          
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
As already seen in example 3, boundary integration 
requires a quadrature rule.  Here, however,
we use the more convenient way of building this
rule at run-time using \p quad_type.  Note that one 
could also have initialized the face quadrature rules 
with the type directly determined from \p qrule, namely 
through:
\verbatim
AutoPtr<QBase>  qface (QBase::build(qrule->type(),
dim-1, 
THIRD));
\endverbatim
And again: using the \p AutoPtr<QBase> relaxes
the need to delete the object afterwards,
they clean up themselves.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;QBase&gt;  qface (QBase::build(quad_type,
                                              dim-1, 
                                              THIRD));
                      
          
</pre>
</div>
<div class = "comment">
Tell the finte element object to use our
quadrature rule.  Note that a \p AutoPtr<QBase> returns
a \p QBase* pointer to the object it handles with \p get().  
However, using \p get(), the \p AutoPtr<QBase> \p qface is 
still in charge of this pointer. I.e., when \p qface goes 
out of scope, it will safely delete the \p QBase object it 
points to.  This behavior may be overridden using
\p AutoPtr<Xyz>::release(), but is not recommended.
</div>

<div class ="fragment">
<pre>
          fe_face-&gt;attach_quadrature_rule (qface.get());
                      
        
          
</pre>
</div>
<div class = "comment">
This is again identical to example 4, and not commented.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt;& JxW = fe-&gt;get_JxW();
          
          const std::vector&lt;Point&gt;& q_point = fe-&gt;get_xyz();
          
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi = fe-&gt;get_phi();
          
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi = fe-&gt;get_dphi();
            
          DenseMatrix&lt;Number&gt; Ke;
          DenseVector&lt;Number&gt; Fe;
          
          std::vector&lt;unsigned int&gt; dof_indices;
          
          
          
          
          
</pre>
</div>
<div class = "comment">
Now we will loop over all the elements in the mesh.
See example 3 for details.
</div>

<div class ="fragment">
<pre>
          MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
          const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
          
          for ( ; el != end_el; ++el)
            {
              const Elem* elem = *el;
              
              dof_map.dof_indices (elem, dof_indices);
              
              fe-&gt;reinit (elem);
              
              Ke.resize (dof_indices.size(),
                         dof_indices.size());
              
              Fe.resize (dof_indices.size());
              
        
        
              
</pre>
</div>
<div class = "comment">
Now loop over the quadrature points.  This handles
the numeric integration.  Note the slightly different
access to the QBase members!
</div>

<div class ="fragment">
<pre>
              for (unsigned int qp=0; qp&lt;qrule-&gt;n_points(); qp++)
                {
</pre>
</div>
<div class = "comment">
Add the matrix contribution
</div>

<div class ="fragment">
<pre>
                  for (unsigned int i=0; i&lt;phi.size(); i++)
                    for (unsigned int j=0; j&lt;phi.size(); j++)
                      Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
                  
                  
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
        
                  const Real fxy = - (uxx + uyy + ((dim==2) ? 0. : uzz));
                  
        
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
If this assembly program were to be used on an adaptive mesh,
we would have to apply any hanging node constraint equations
Call heterogenously_constrain_element_matrix_and_vector to impose
non-homogeneous Dirichlet BCs
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
</div>

<div class ="fragment">
<pre>
              system.matrix-&gt;add_matrix (Ke, dof_indices);
              system.rhs-&gt;add_vector    (Fe, dof_indices);
              
            } // end of element loop
          
          
          
          
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
<pre> 
  
  
  #include &lt;iostream&gt;
  #include &lt;sstream&gt; 
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_rules.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;dirichlet_boundaries.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;analytic_function.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  
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
    output(0) = exact_solution(p(0),p(1),p(2));
  }
  
  
  QuadratureType quad_type=INVALID_Q_RULE;
  
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
    
    <B><FONT COLOR="#A020F0">if</FONT></B> (argc &lt; 3)
      {
        <B><FONT COLOR="#A020F0">if</FONT></B> (libMesh::processor_id() == 0)
          {
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Usage: &quot;</FONT></B> &lt;&lt; argv[0] &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; -q n&quot;</FONT></B>
                      &lt;&lt; std::endl;
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;  where n stands for:&quot;</FONT></B> &lt;&lt; std::endl;
  
        
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n=0; n&lt;QuadratureRules::num_valid_elem_rules; n++)
              <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;  &quot;</FONT></B> &lt;&lt; QuadratureRules::valid_elem_rules[n] &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;    &quot;</FONT></B> 
                        &lt;&lt; QuadratureRules::name(QuadratureRules::valid_elem_rules[n])
                        &lt;&lt; std::endl;
        
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; std::endl;
          }
        
        libmesh_error();
      }
    
    
    <B><FONT COLOR="#A020F0">else</FONT></B> 
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Running &quot;</FONT></B> &lt;&lt; argv[0];
        
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">int</FONT></B> i=1; i&lt;argc; i++)
          <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; &quot;</FONT></B> &lt;&lt; argv[i];
        
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; std::endl &lt;&lt; std::endl;
      }
    
  
    quad_type = static_cast&lt;QuadratureType&gt;(std::atoi(argv[2]));
  
    libmesh_example_assert(3 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;3D support&quot;</FONT></B>);
    
    Mesh mesh;
  
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_cube (mesh,
                                       16, 16, 16,
                                       -1., 1.,
                                       -1., 1.,
                                       -1., 1.,
                                       HEX8);
    
    mesh.print_info();
    
    EquationSystems equation_systems (mesh);
    
    equation_systems.add_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>).add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, FIRST);
  
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>).attach_assemble_function (assemble_poisson);
  
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::set&lt;boundary_id_type&gt; boundary_ids;
    boundary_ids.insert(0);
    boundary_ids.insert(1);
    boundary_ids.insert(2);
    boundary_ids.insert(3);
    boundary_ids.insert(4);
    boundary_ids.insert(5);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; variables(1);
    variables[0] = u_var;
    
    AnalyticFunction&lt;&gt; exact_solution_object(exact_solution_wrapper);
    
    DirichletBoundary dirichlet_bc(boundary_ids,
                                   variables,
                                   &amp;exact_solution_object);
  
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>).get_dof_map().add_dirichlet_boundary(dirichlet_bc);
  
    equation_systems.init();
    
    equation_systems.print_info();
  
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>).solve();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream f_name;
    f_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out_&quot;</FONT></B> &lt;&lt; quad_type &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.e&quot;</FONT></B>;
  
  #ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO(mesh).write_equation_systems (f_name.str(),
                                        equation_systems);
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_poisson(EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
  
    LinearImplicitSystem&amp; system = es.get_system&lt;LinearImplicitSystem&gt;(<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap&amp; dof_map = system.get_dof_map();
    
    FEType fe_type = dof_map.variable_type(0);
  
    
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    
    AutoPtr&lt;QBase&gt; qrule(QBase::build(quad_type, dim, THIRD));
  
  
    
    fe-&gt;attach_quadrature_rule (qrule.get());
  
    
    AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
    
    
    AutoPtr&lt;QBase&gt;  qface (QBase::build(quad_type,
                                        dim-1, 
                                        THIRD));
                
    
    fe_face-&gt;attach_quadrature_rule (qface.get());
                
  
    
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
        
        dof_map.dof_indices (elem, dof_indices);
        
        fe-&gt;reinit (elem);
        
        Ke.resize (dof_indices.size(),
                   dof_indices.size());
        
        Fe.resize (dof_indices.size());
        
  
  
        
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule-&gt;n_points(); qp++)
          {
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi.size(); j++)
                Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
            
            
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
  
            <B><FONT COLOR="#228B22">const</FONT></B> Real fxy = - (uxx + uyy + ((dim==2) ? 0. : uzz));
            
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi.size(); i++)
              Fe(i) += JxW[qp]*fxy*phi[i][qp];          
          }     
        
        dof_map.heterogenously_constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        
        system.matrix-&gt;add_matrix (Ke, dof_indices);
        system.rhs-&gt;add_vector    (Fe, dof_indices);
        
      } <I><FONT COLOR="#B22222">// end of element loop
</FONT></I>    
    
    
    
    <B><FONT COLOR="#A020F0">return</FONT></B>;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
Linking introduction_ex5-opt...
***************************************************************
* Running Example  mpirun -np 6 ./introduction_ex5-opt -q 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running ./introduction_ex5-opt -q 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=4913
    n_local_nodes()=958
  n_elem()=4096
    n_local_elem()=683
    n_active_elem()=4096
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
    n_dofs()=4913
    n_local_dofs()=958
    n_constrained_dofs()=1542
    n_local_constrained_dofs()=266
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 21.8685
      Average Off-Processor Bandwidth <= 2.46242
      Maximum  On-Processor Bandwidth <= 27
      Maximum Off-Processor Bandwidth <= 17
    DofMap Constraints
      Number of DoF Constraints = 1538
      Number of Heterogenous Constraints= 1474
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./introduction_ex5-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:21:41 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           3.108e-01      1.13999   2.799e-01
Objects:              6.000e+01      1.00000   6.000e+01
Flops:                1.647e+07      1.76243   1.235e+07  7.410e+07
Flops/sec:            5.299e+07      1.56762   4.395e+07  2.637e+08
MPI Messages:         1.925e+02      1.23794   1.793e+02  1.076e+03
MPI Message Lengths:  1.838e+05      1.22884   9.227e+02  9.928e+05
MPI Reductions:       1.090e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.7978e-01 100.0%  7.4103e+07 100.0%  1.076e+03 100.0%  9.227e+02      100.0%  9.300e+01  85.3% 

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

VecMDot               30 1.0 1.9555e-02 1.5 8.90e+05 1.5 0.0e+00 0.0e+00 3.0e+01  6  6  0  0 28   6  6  0  0 32   234
VecNorm               32 1.0 7.3779e-03 2.7 6.13e+04 1.5 0.0e+00 0.0e+00 3.2e+01  2  0  0  0 29   2  0  0  0 34    43
VecScale              31 1.0 4.8399e-05 1.4 2.97e+04 1.5 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  3147
VecCopy                4 1.0 8.8215e-06 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                36 1.0 1.2660e-04 5.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                2 1.0 1.1325e-04 3.9 3.83e+03 1.5 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   174
VecMAXPY              31 1.0 3.6240e-04 1.4 9.48e+05 1.5 0.0e+00 0.0e+00 0.0e+00  0  7  0  0  0   0  7  0  0  0 13421
VecAssemblyBegin       3 1.0 5.1432e-0313.2 0.00e+00 0.0 2.6e+01 1.8e+03 9.0e+00  2  0  2  5  8   2  0  2  5 10     0
VecAssemblyEnd         3 1.0 3.6001e-05 1.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       32 1.0 2.4533e-04 1.5 0.00e+00 0.0 9.0e+02 5.8e+02 0.0e+00  0  0 83 52  0   0  0 83 52  0     0
VecScatterEnd         32 1.0 2.7973e-02 4.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  7  0  0  0  0   7  0  0  0  0     0
VecNormalize          31 1.0 7.1006e-03 2.8 8.91e+04 1.5 0.0e+00 0.0e+00 3.1e+01  1  1  0  0 28   1  1  0  0 33    64
MatMult               31 1.0 2.9124e-02 4.3 1.42e+06 1.5 8.7e+02 5.7e+02 0.0e+00  7 10 81 50  0   7 10 81 50  0   245
MatSolve              31 1.0 4.5779e-03 1.6 7.44e+06 1.6 0.0e+00 0.0e+00 0.0e+00  1 47  0  0  0   1 47  0  0  0  7533
MatLUFactorNum         1 1.0 1.3123e-02 4.4 5.68e+06 2.8 0.0e+00 0.0e+00 0.0e+00  2 30  0  0  0   2 30  0  0  0  1719
MatILUFactorSym        1 1.0 3.0813e-02 2.4 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+00  7  0  0  0  1   7  0  0  0  1     0
MatAssemblyBegin       2 1.0 2.1869e-0224.4 0.00e+00 0.0 3.9e+01 1.0e+04 4.0e+00  3  0  4 41  4   3  0  4 41  4     0
MatAssemblyEnd         2 1.0 1.9920e-03 1.3 0.00e+00 0.0 5.6e+01 1.4e+02 8.0e+00  1  0  5  1  7   1  0  5  1  9     0
MatGetRowIJ            1 1.0 3.2902e-0517.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 9.7990e-05 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries         3 1.0 3.0923e-04 5.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        30 1.0 1.9917e-02 1.5 1.78e+06 1.5 0.0e+00 0.0e+00 3.0e+01  6 12  0  0 28   6 12  0  0 32   459
KSPSetup               2 1.0 7.8917e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 7.2039e-02 1.0 1.65e+07 1.8 8.7e+02 5.7e+02 6.5e+01 26100 81 50 60  26100 81 50 70  1029
PCSetUp                2 1.0 3.8911e-02 2.4 5.68e+06 2.8 0.0e+00 0.0e+00 3.0e+00  9 30  0  0  3   9 30  0  0  3   580
PCSetUpOnBlocks        1 1.0 3.8607e-02 2.4 5.68e+06 2.8 0.0e+00 0.0e+00 3.0e+00  9 30  0  0  3   9 30  0  0  3   584
PCApply               31 1.0 4.8978e-03 1.5 7.44e+06 1.6 0.0e+00 0.0e+00 0.0e+00  1 47  0  0  0   1 47  0  0  0  7041
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    42             42       345344     0
         Vec Scatter     2              2         1736     0
           Index Set     7              7        17920     0
   IS L to G Mapping     1              1         5632     0
              Matrix     4              4      1765508     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.282e-05
Average time for zero size MPI_Send(): 8.10226e-05
#PETSc Option Table entries:
-ksp_right_pc
-log_summary
-pc_type bjacobi
-q 0
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

-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 6                                                |
| Time:           Fri Aug 24 15:21:41 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.545511, Active time=0.239951                                                     |
 --------------------------------------------------------------------------------------------------------------------
| Event                                  nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                  w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------------|
|                                                                                                                    |
|                                                                                                                    |
| DofMap                                                                                                             |
|   add_neighbors_to_send_list()         1         0.0007      0.000726    0.0010      0.000954    0.30     0.40     |
|   build_constraint_matrix_and_vector() 683       0.0010      0.000002    0.0010      0.000002    0.43     0.43     |
|   build_sparsity()                     1         0.0102      0.010202    0.0123      0.012291    4.25     5.12     |
|   create_dof_constraints()             1         0.0086      0.008632    0.0114      0.011400    3.60     4.75     |
|   distribute_dofs()                    1         0.0013      0.001254    0.0187      0.018731    0.52     7.81     |
|   dof_indices()                        10100     0.0052      0.000001    0.0052      0.000001    2.16     2.16     |
|   hetero_cnstrn_elem_mat_vec()         683       0.0221      0.000032    0.0221      0.000032    9.22     9.22     |
|   prepare_send_list()                  1         0.0001      0.000079    0.0001      0.000079    0.03     0.03     |
|   reinit()                             1         0.0033      0.003348    0.0033      0.003348    1.40     1.40     |
|                                                                                                                    |
| EquationSystems                                                                                                    |
|   build_solution_vector()              1         0.0011      0.001108    0.0031      0.003107    0.46     1.29     |
|                                                                                                                    |
| ExodusII_IO                                                                                                        |
|   write_nodal_data()                   1         0.0046      0.004565    0.0046      0.004565    1.90     1.90     |
|                                                                                                                    |
| FE                                                                                                                 |
|   compute_shape_functions()            683       0.0026      0.000004    0.0026      0.000004    1.10     1.10     |
|   init_shape_functions()               1         0.0000      0.000020    0.0000      0.000020    0.01     0.01     |
|                                                                                                                    |
| FEMap                                                                                                              |
|   compute_affine_map()                 683       0.0009      0.000001    0.0009      0.000001    0.38     0.38     |
|   init_reference_to_physical_map()     1         0.0000      0.000029    0.0000      0.000029    0.01     0.01     |
|                                                                                                                    |
| Mesh                                                                                                               |
|   find_neighbors()                     1         0.0063      0.006284    0.0068      0.006790    2.62     2.83     |
|   renumber_nodes_and_elem()            2         0.0005      0.000268    0.0005      0.000268    0.22     0.22     |
|                                                                                                                    |
| MeshCommunication                                                                                                  |
|   compute_hilbert_indices()            2         0.0282      0.014125    0.0282      0.014125    11.77    11.77    |
|   find_global_indices()                2         0.0017      0.000853    0.0394      0.019688    0.71     16.41    |
|   parallel_sort()                      2         0.0028      0.001388    0.0035      0.001771    1.16     1.48     |
|                                                                                                                    |
| MeshOutput                                                                                                         |
|   write_equation_systems()             1         0.0000      0.000026    0.0077      0.007698    0.01     3.21     |
|                                                                                                                    |
| MeshTools::Generation                                                                                              |
|   build_cube()                         1         0.0017      0.001653    0.0017      0.001653    0.69     0.69     |
|                                                                                                                    |
| MetisPartitioner                                                                                                   |
|   partition()                          1         0.0088      0.008784    0.0315      0.031476    3.66     13.12    |
|                                                                                                                    |
| Parallel                                                                                                           |
|   allgather()                          8         0.0049      0.000617    0.0049      0.000617    2.06     2.06     |
|   broadcast()                          1         0.0000      0.000013    0.0000      0.000013    0.01     0.01     |
|   gather()                             1         0.0000      0.000005    0.0000      0.000005    0.00     0.00     |
|   max(scalar)                          2         0.0030      0.001501    0.0030      0.001501    1.25     1.25     |
|   max(vector)                          2         0.0002      0.000122    0.0002      0.000122    0.10     0.10     |
|   min(vector)                          2         0.0013      0.000655    0.0013      0.000655    0.55     0.55     |
|   probe()                              50        0.0225      0.000451    0.0225      0.000451    9.40     9.40     |
|   receive()                            50        0.0001      0.000003    0.0227      0.000454    0.06     9.46     |
|   send()                               50        0.0001      0.000002    0.0001      0.000002    0.03     0.03     |
|   send_receive()                       54        0.0001      0.000002    0.0229      0.000425    0.05     9.56     |
|   sum()                                10        0.0034      0.000345    0.0034      0.000345    1.44     1.44     |
|                                                                                                                    |
| Parallel::Request                                                                                                  |
|   wait()                               50        0.0000      0.000001    0.0000      0.000001    0.02     0.02     |
|                                                                                                                    |
| Partitioner                                                                                                        |
|   set_node_processor_ids()             1         0.0006      0.000578    0.0079      0.007895    0.24     3.29     |
|   set_parent_processor_ids()           1         0.0002      0.000235    0.0002      0.000235    0.10     0.10     |
|                                                                                                                    |
| PetscLinearSolver                                                                                                  |
|   solve()                              1         0.0825      0.082511    0.0825      0.082511    34.39    34.39    |
|                                                                                                                    |
| System                                                                                                             |
|   assemble()                           1         0.0089      0.008855    0.0363      0.036317    3.69     15.14    |
 --------------------------------------------------------------------------------------------------------------------
| Totals:                                13138     0.2400                                          100.00            |
 --------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 6 ./introduction_ex5-opt -q 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
