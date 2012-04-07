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
Compiling C++ (in optimized mode) introduction_ex5.C...
Linking introduction_ex5-opt...
***************************************************************
* Running Example  ./introduction_ex5-opt -q 0
***************************************************************
 
Running ./introduction_ex5-opt -q 0

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=4913
    n_local_nodes()=4913
  n_elem()=4096
    n_local_elem()=4096
    n_active_elem()=4096
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
    n_dofs()=4913
    n_local_dofs()=4913
    n_constrained_dofs()=1538
    n_local_constrained_dofs()=1538
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 23.9465
      Average Off-Processor Bandwidth <= 0
      Maximum  On-Processor Bandwidth <= 27
      Maximum Off-Processor Bandwidth <= 0
    DofMap Constraints
      Number of DoF Constraints = 1538
      Number of Heterogenous Constraints= 1474
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0


-------------------------------------------------------------------
| Time:           Sat Apr  7 15:59:38 2012                         |
| OS:             Linux                                            |
| HostName:       lkirk-home                                       |
| OS Release:     3.0.0-17-generic                                 |
| OS Version:     #30-Ubuntu SMP Thu Mar 8 20:45:39 UTC 2012       |
| Machine:        x86_64                                           |
| Username:       benkirk                                          |
| Configuration:  ./configure run on Sat Apr  7 15:49:27 CDT 2012  |
-------------------------------------------------------------------
 --------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.385387, Active time=0.300489                                                     |
 --------------------------------------------------------------------------------------------------------------------
| Event                                  nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                  w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------------|
|                                                                                                                    |
|                                                                                                                    |
| DofMap                                                                                                             |
|   add_neighbors_to_send_list()         1         0.0024      0.002396    0.0024      0.002396    0.80     0.80     |
|   build_constraint_matrix_and_vector() 4096      0.0118      0.000003    0.0118      0.000003    3.92     3.92     |
|   compute_sparsity()                   1         0.0174      0.017435    0.0215      0.021549    5.80     7.17     |
|   create_dof_constraints()             1         0.0199      0.019863    0.0291      0.029109    6.61     9.69     |
|   distribute_dofs()                    1         0.0031      0.003063    0.0111      0.011051    1.02     3.68     |
|   dof_indices()                        16384     0.0248      0.000002    0.0248      0.000002    8.26     8.26     |
|   hetero_cnstrn_elem_mat_vec()         4096      0.0198      0.000005    0.0198      0.000005    6.58     6.58     |
|   prepare_send_list()                  1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   reinit()                             1         0.0080      0.007987    0.0080      0.007987    2.66     2.66     |
|                                                                                                                    |
| EquationSystems                                                                                                    |
|   build_solution_vector()              1         0.0058      0.005774    0.0096      0.009612    1.92     3.20     |
|                                                                                                                    |
| ExodusII_IO                                                                                                        |
|   write_nodal_data()                   1         0.0090      0.008966    0.0090      0.008966    2.98     2.98     |
|                                                                                                                    |
| FE                                                                                                                 |
|   compute_affine_map()                 4096      0.0071      0.000002    0.0071      0.000002    2.37     2.37     |
|   compute_shape_functions()            4096      0.0048      0.000001    0.0048      0.000001    1.59     1.59     |
|   init_shape_functions()               1         0.0001      0.000119    0.0001      0.000119    0.04     0.04     |
|                                                                                                                    |
| Mesh                                                                                                               |
|   find_neighbors()                     1         0.0206      0.020576    0.0206      0.020576    6.85     6.85     |
|   renumber_nodes_and_elem()            2         0.0020      0.001018    0.0020      0.001018    0.68     0.68     |
|                                                                                                                    |
| MeshOutput                                                                                                         |
|   write_equation_systems()             1         0.0000      0.000042    0.0186      0.018620    0.01     6.20     |
|                                                                                                                    |
| MeshTools::Generation                                                                                              |
|   build_cube()                         1         0.0097      0.009686    0.0097      0.009686    3.22     3.22     |
|                                                                                                                    |
| Parallel                                                                                                           |
|   allgather()                          1         0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                                    |
| Partitioner                                                                                                        |
|   single_partition()                   1         0.0003      0.000303    0.0003      0.000303    0.10     0.10     |
|                                                                                                                    |
| PetscLinearSolver                                                                                                  |
|   solve()                              1         0.0338      0.033791    0.0338      0.033791    11.25    11.25    |
|                                                                                                                    |
| System                                                                                                             |
|   assemble()                           1         0.1002      0.100212    0.1569      0.156934    33.35    52.23    |
 --------------------------------------------------------------------------------------------------------------------
| Totals:                                32786     0.3005                                          100.00            |
 --------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  ./introduction_ex5-opt -q 0
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
