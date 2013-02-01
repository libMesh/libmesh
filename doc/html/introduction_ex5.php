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
<br><br><br> <h1> The source file introduction_ex5.C with comments: </h1> 
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
        #include "libmesh/libmesh.h"
        #include "libmesh/mesh.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/exodusII_io.h"
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
Define the base quadrature class, with which
specialized quadrature rules will be built.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/quadrature.h"
        
</pre>
</div>
<div class = "comment">
Include the namespace \p QuadratureRules for
some handy descriptions.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/quadrature_rules.h"
        
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
Define the DofMap, which handles degree of freedom
indexing.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/dof_map.h"
        
</pre>
</div>
<div class = "comment">
To impose Dirichlet boundary conditions
</div>

<div class ="fragment">
<pre>
        #include "libmesh/dirichlet_boundaries.h"
        #include "libmesh/analytic_function.h"
        
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
          libmesh_assert_equal_to (system_name, "Poisson");
        
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
          
          std::vector&lt;dof_id_type&gt; dof_indices;
          
          
          
          
          
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
<br><br><br> <h1> The source file introduction_ex5.C without comments: </h1> 
<pre> 
  
  
  #include &lt;iostream&gt;
  #include &lt;sstream&gt; 
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_rules.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dirichlet_boundaries.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/analytic_function.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  
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
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>);
  
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
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
    
    
    
    
    
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
***************************************************************
* Running Example introduction_ex5:
*  mpirun -np 12 example-devel -q 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running /workspace/libmesh/examples/introduction/introduction_ex5/.libs/lt-example-devel -q 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=4913
    n_local_nodes()=540
  n_elem()=4096
    n_local_elem()=343
    n_active_elem()=4096
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
    n_dofs()=4913
    n_local_dofs()=540
    n_constrained_dofs()=1538
    n_local_constrained_dofs()=184
    n_vectors()=1
    n_matrices()=1
    DofMap Sparsity
      Average  On-Processor Bandwidth <= 21.9516
      Average Off-Processor Bandwidth <= 4.4311
      Maximum  On-Processor Bandwidth <= 36
      Maximum Off-Processor Bandwidth <= 30
    DofMap Constraints
      Number of DoF Constraints = 1538
      Number of Heterogenous Constraints= 1474
      Average DoF Constraint Length= 0
      Number of Node Constraints = 0

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/introduction/introduction_ex5/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 21:57:17 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           2.151e+00      1.00000   2.151e+00
Objects:              6.200e+01      1.00000   6.200e+01
Flops:                7.624e+06      3.48639   5.216e+06  6.259e+07
Flops/sec:            3.544e+06      3.48639   2.425e+06  2.909e+07
MPI Messages:         4.525e+02      2.48626   2.873e+02  3.448e+03
MPI Message Lengths:  2.258e+05      1.93855   5.116e+02  1.764e+06
MPI Reductions:       1.200e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.1512e+00 100.0%  6.2589e+07 100.0%  3.448e+03 100.0%  5.116e+02      100.0%  1.190e+02  99.2% 

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

VecMDot               36 1.0 3.1800e-03 6.9 5.24e+05 2.2 0.0e+00 0.0e+00 3.6e+01  0  8  0  0 30   0  8  0  0 30  1500
VecNorm               39 1.0 1.9608e-03 8.3 4.21e+04 2.2 0.0e+00 0.0e+00 3.9e+01  0  1  0  0 32   0  1  0  0 33   195
VecScale              38 1.0 7.4387e-05 1.2 2.05e+04 2.2 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2510
VecCopy                3 1.0 7.1526e-06 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                46 1.0 4.9353e-05 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                4 1.0 4.5061e-05 2.4 4.32e+03 2.2 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   872
VecMAXPY              38 1.0 3.3164e-04 2.0 5.64e+05 2.2 0.0e+00 0.0e+00 0.0e+00  0  8  0  0  0   0  8  0  0  0 15466
VecAssemblyBegin       3 1.0 2.1005e-04 1.0 0.00e+00 0.0 7.2e+01 1.1e+03 9.0e+00  0  0  2  4  8   0  0  2  4  8     0
VecAssemblyEnd         3 1.0 4.1962e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       39 1.0 3.0351e-04 1.4 0.00e+00 0.0 3.0e+03 3.5e+02 0.0e+00  0  0 86 59  0   0  0 86 59  0     0
VecScatterEnd         39 1.0 8.2798e-0350.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecNormalize          38 1.0 2.0566e-03 5.0 6.16e+04 2.2 0.0e+00 0.0e+00 3.8e+01  0  1  0  0 32   0  1  0  0 32   272
MatMult               38 1.0 9.0966e-03 7.2 9.51e+05 2.1 2.9e+03 3.4e+02 0.0e+00  0 14 84 56  0   0 14 84 56  0   962
MatSolve              39 1.0 3.3357e-03 4.1 4.00e+06 4.3 0.0e+00 0.0e+00 0.0e+00  0 51  0  0  0   0 51  0  0  0  9543
MatLUFactorNum         1 1.0 2.1801e-03 5.4 1.61e+06 6.4 0.0e+00 0.0e+00 0.0e+00  0 18  0  0  0   0 18  0  0  0  5272
MatILUFactorSym        1 1.0 9.0692e-03 6.7 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  0  0  0  0  2   0  0  0  0  3     0
MatAssemblyBegin       2 1.0 1.0839e-0238.4 0.00e+00 0.0 1.1e+02 5.7e+03 4.0e+00  0  0  3 35  3   0  0  3 35  3     0
MatAssemblyEnd         2 1.0 1.1630e-03 1.2 0.00e+00 0.0 1.5e+02 8.8e+01 8.0e+00  0  0  4  1  7   0  0  4  1  7     0
MatGetRowIJ            1 1.0 1.0967e-05 5.8 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 9.9182e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries         3 1.0 8.2970e-05 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog        36 1.0 3.4204e-03 4.0 1.05e+06 2.2 0.0e+00 0.0e+00 3.6e+01  0 15  0  0 30   0 15  0  0 30  2791
KSPSetUp               2 1.0 1.3113e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 1.9294e-02 1.0 7.62e+06 3.5 2.9e+03 3.4e+02 8.2e+01  1100 84 56 68   1100 84 56 69  3244
PCSetUp                2 1.0 1.2045e-02 4.7 1.61e+06 6.4 0.0e+00 0.0e+00 7.0e+00  0 18  0  0  6   0 18  0  0  6   954
PCSetUpOnBlocks        1 1.0 1.1517e-02 5.7 1.61e+06 6.4 0.0e+00 0.0e+00 5.0e+00  0 18  0  0  4   0 18  0  0  4   998
PCApply               39 1.0 3.7823e-03 2.9 4.00e+06 4.3 0.0e+00 0.0e+00 0.0e+00  0 51  0  0  0   0 51  0  0  0  8416
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector    43             43       232584     0
      Vector Scatter     2              2         2072     0
           Index Set     7              7         9528     0
   IS L to G Mapping     1              1          564     0
              Matrix     4              4       824052     0
       Krylov Solver     2              2        19360     0
      Preconditioner     2              2         1784     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 4.3869e-06
Average time for zero size MPI_Send(): 1.357e-05
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


 ----------------------------------------------------------------------------------------------------------------------
| Processor id:   0                                                                                                    |
| Num Processors: 12                                                                                                   |
| Time:           Thu Jan 31 21:57:17 2013                                                                             |
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
 --------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=2.3202, Active time=2.08344                                                        |
 --------------------------------------------------------------------------------------------------------------------
| Event                                  nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                  w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|--------------------------------------------------------------------------------------------------------------------|
|                                                                                                                    |
|                                                                                                                    |
| DofMap                                                                                                             |
|   add_neighbors_to_send_list()         1         0.0398      0.039776    0.0834      0.083395    1.91     4.00     |
|   build_constraint_matrix_and_vector() 343       0.0063      0.000018    0.0063      0.000018    0.30     0.30     |
|   build_sparsity()                     1         0.0279      0.027930    0.0697      0.069729    1.34     3.35     |
|   create_dof_constraints()             1         0.1188      0.118824    0.5631      0.563056    5.70     27.03    |
|   distribute_dofs()                    1         0.0830      0.082988    0.2946      0.294642    3.98     14.14    |
|   dof_indices()                        5522      0.5988      0.000108    0.5988      0.000108    28.74    28.74    |
|   hetero_cnstrn_elem_mat_vec()         343       0.0148      0.000043    0.0148      0.000043    0.71     0.71     |
|   prepare_send_list()                  1         0.0006      0.000585    0.0006      0.000585    0.03     0.03     |
|   reinit()                             1         0.2084      0.208383    0.2084      0.208383    10.00    10.00    |
|                                                                                                                    |
| EquationSystems                                                                                                    |
|   build_solution_vector()              1         0.0055      0.005470    0.0425      0.042456    0.26     2.04     |
|                                                                                                                    |
| ExodusII_IO                                                                                                        |
|   write_nodal_data()                   1         0.0435      0.043484    0.0435      0.043484    2.09     2.09     |
|                                                                                                                    |
| FE                                                                                                                 |
|   compute_shape_functions()            343       0.0186      0.000054    0.0186      0.000054    0.89     0.89     |
|   init_shape_functions()               1         0.0003      0.000267    0.0003      0.000267    0.01     0.01     |
|                                                                                                                    |
| FEMap                                                                                                              |
|   compute_affine_map()                 343       0.0073      0.000021    0.0073      0.000021    0.35     0.35     |
|   init_reference_to_physical_map()     1         0.0002      0.000221    0.0002      0.000221    0.01     0.01     |
|                                                                                                                    |
| Mesh                                                                                                               |
|   find_neighbors()                     1         0.1529      0.152910    0.1589      0.158853    7.34     7.62     |
|   renumber_nodes_and_elem()            2         0.0061      0.003037    0.0061      0.003037    0.29     0.29     |
|                                                                                                                    |
| MeshCommunication                                                                                                  |
|   compute_hilbert_indices()            2         0.0681      0.034073    0.0681      0.034073    3.27     3.27     |
|   find_global_indices()                2         0.0277      0.013835    0.1043      0.052151    1.33     5.01     |
|   parallel_sort()                      2         0.0057      0.002861    0.0070      0.003502    0.27     0.34     |
|                                                                                                                    |
| MeshOutput                                                                                                         |
|   write_equation_systems()             1         0.0001      0.000144    0.0862      0.086223    0.01     4.14     |
|                                                                                                                    |
| MeshTools::Generation                                                                                              |
|   build_cube()                         1         0.0258      0.025752    0.0258      0.025752    1.24     1.24     |
|                                                                                                                    |
| MetisPartitioner                                                                                                   |
|   partition()                          1         0.5064      0.506430    0.5575      0.557514    24.31    26.76    |
|                                                                                                                    |
| Parallel                                                                                                           |
|   allgather()                          9         0.0006      0.000063    0.0006      0.000071    0.03     0.03     |
|   max(bool)                            1         0.0000      0.000007    0.0000      0.000007    0.00     0.00     |
|   max(scalar)                          105       0.0008      0.000008    0.0008      0.000008    0.04     0.04     |
|   max(vector)                          24        0.0003      0.000014    0.0008      0.000034    0.02     0.04     |
|   min(bool)                            121       0.0008      0.000007    0.0008      0.000007    0.04     0.04     |
|   min(scalar)                          99        0.0366      0.000370    0.0366      0.000370    1.76     1.76     |
|   min(vector)                          24        0.0005      0.000019    0.0016      0.000065    0.02     0.07     |
|   probe()                              132       0.0028      0.000021    0.0028      0.000021    0.14     0.14     |
|   receive()                            132       0.0010      0.000008    0.0039      0.000029    0.05     0.19     |
|   send()                               132       0.0005      0.000004    0.0005      0.000004    0.02     0.02     |
|   send_receive()                       136       0.0016      0.000012    0.0063      0.000047    0.08     0.30     |
|   sum()                                20        0.0014      0.000071    0.0019      0.000096    0.07     0.09     |
|                                                                                                                    |
| Parallel::Request                                                                                                  |
|   wait()                               132       0.0003      0.000002    0.0003      0.000002    0.01     0.01     |
|                                                                                                                    |
| Partitioner                                                                                                        |
|   set_node_processor_ids()             1         0.0093      0.009311    0.0134      0.013446    0.45     0.65     |
|   set_parent_processor_ids()           1         0.0074      0.007397    0.0074      0.007397    0.36     0.36     |
|                                                                                                                    |
| PetscLinearSolver                                                                                                  |
|   solve()                              1         0.0239      0.023890    0.0239      0.023890    1.15     1.15     |
|                                                                                                                    |
| System                                                                                                             |
|   assemble()                           1         0.0291      0.029073    0.1160      0.115991    1.40     5.57     |
 --------------------------------------------------------------------------------------------------------------------
| Totals:                                7987      2.0834                                          100.00            |
 --------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example introduction_ex5:
*  mpirun -np 12 example-devel -q 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
