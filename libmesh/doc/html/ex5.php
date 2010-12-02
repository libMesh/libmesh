<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex5",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 5 - Run-Time Quadrature Rule Selection</h1>

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
          
          equation_systems.get_system("Poisson").add_variable("u", FIRST);
        
          equation_systems.get_system("Poisson").attach_assemble_function (assemble_poisson);
        
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
          f_name &lt;&lt; "out_" &lt;&lt; quad_type &lt;&lt; ".exd";
        
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
Most of this has already been seen before, except
for the build routines of QBase, described below
</div>

<div class ="fragment">
<pre>
              {
                for (unsigned int side=0; side&lt;elem-&gt;n_sides(); side++)
                  if (elem-&gt;neighbor(side) == NULL)
                    {              
                      const std::vector&lt;std::vector&lt;Real&gt; &gt;& phi_face    = fe_face-&gt;get_phi();
                      const std::vector&lt;Real&gt;&               JxW_face    = fe_face-&gt;get_JxW();              
                      const std::vector&lt;Point &gt;&             qface_point = fe_face-&gt;get_xyz();
                      
                      
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
Note that the \p AutoPtr<QBase> overloaded the operator->,
so that QBase methods may safely be accessed.  It may
be said: accessing an \p AutoPtr<Xyz> through the
"." operator returns \p AutoPtr methods, while access
through the "->" operator returns Xyz methods.
This allows almost no change in syntax when switching
to "safe pointers".
</div>

<div class ="fragment">
<pre>
                      for (unsigned int qp=0; qp&lt;qface-&gt;n_points(); qp++)
                        {
                          const Real xf = qface_point[qp](0);
                          const Real yf = qface_point[qp](1);
                          const Real zf = qface_point[qp](2);
                          
                          const Real penalty = 1.e10;
                          
                          const Real value = exact_solution(xf, yf, zf);
                          
                          for (unsigned int i=0; i&lt;phi_face.size(); i++)
                            for (unsigned int j=0; j&lt;phi_face.size(); j++)
                              Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];
                          
                          
                          for (unsigned int i=0; i&lt;phi_face.size(); i++)
                            Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
                          
                        } // end face quadrature point loop          
                    } // end if (elem-&gt;neighbor(side) == NULL)
              } // end boundary condition section          
              
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
  
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  
  
  
  
  
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_poisson(EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  
  
  Real exact_solution (<B><FONT COLOR="#228B22">const</FONT></B> Real x,
                       <B><FONT COLOR="#228B22">const</FONT></B> Real y,
                       <B><FONT COLOR="#228B22">const</FONT></B> Real z = 0.);
  
  
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
    
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>).add_variable(<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, FIRST);
  
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>).attach_assemble_function (assemble_poisson);
  
    equation_systems.init();
    
    equation_systems.print_info();
  
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Poisson&quot;</FONT></B>).solve();
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::ostringstream f_name;
    f_name &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;out_&quot;</FONT></B> &lt;&lt; quad_type &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;.exd&quot;</FONT></B>;
  
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
  
  
  
  
        
        
        {
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> side=0; side&lt;elem-&gt;n_sides(); side++)
            <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(side) == NULL)
              {              
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi_face    = fe_face-&gt;get_phi();
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp;               JxW_face    = fe_face-&gt;get_JxW();              
                <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point &gt;&amp;             qface_point = fe_face-&gt;get_xyz();
                
                
                fe_face-&gt;reinit(elem, side);
                
                
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qface-&gt;n_points(); qp++)
                  {
                    <B><FONT COLOR="#228B22">const</FONT></B> Real xf = qface_point[qp](0);
                    <B><FONT COLOR="#228B22">const</FONT></B> Real yf = qface_point[qp](1);
                    <B><FONT COLOR="#228B22">const</FONT></B> Real zf = qface_point[qp](2);
                    
                    <B><FONT COLOR="#228B22">const</FONT></B> Real penalty = 1.e10;
                    
                    <B><FONT COLOR="#228B22">const</FONT></B> Real value = exact_solution(xf, yf, zf);
                    
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi_face.size(); i++)
                      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;phi_face.size(); j++)
                        Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];
                    
                    
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;phi_face.size(); i++)
                      Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
                    
                  } <I><FONT COLOR="#B22222">// end face quadrature point loop          
</FONT></I>              } <I><FONT COLOR="#B22222">// end if (elem-&gt;neighbor(side) == NULL)
</FONT></I>        } <I><FONT COLOR="#B22222">// end boundary condition section          
</FONT></I>        
        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        
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
Compiling C++ (in optimized mode) ex5.C...
/org/centers/pecos/LIBRARIES/GCC/gcc-4.5.1-lucid/libexec/gcc/x86_64-unknown-linux-gnu/4.5.1/cc1plus: error while loading shared libraries: libmpc.so.2: cannot open shared object file: No such file or directory
make[1]: *** [ex5.x86_64-unknown-linux-gnu.opt.o] Error 1
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
