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
***************************************************************
* Running Example  mpirun -np 2 ./ex5-opt -q 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
Running ./ex5-opt -q 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary

 Mesh Information:
  mesh_dimension()=3
  spatial_dimension()=3
  n_nodes()=4913
    n_local_nodes()=2601
  n_elem()=4096
    n_local_elem()=2048
    n_active_elem()=4096
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
    n_dofs()=4913
    n_local_dofs()=2601
    n_constrained_dofs()=0
    n_vectors()=1

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex5-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Thu Feb  3 12:13:33 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           2.102e-01      1.02292   2.079e-01
Objects:              4.300e+01      1.00000   4.300e+01
Flops:                3.992e+07      1.22999   3.619e+07  7.238e+07
Flops/sec:            1.899e+08      1.20243   1.739e+08  3.478e+08
MPI Messages:         1.950e+01      1.00000   1.950e+01  3.900e+01
MPI Message Lengths:  1.208e+05      1.01951   6.137e+03  2.393e+05
MPI Reductions:       6.600e+01      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 2.0786e-01 100.0%  7.2376e+07 100.0%  3.900e+01 100.0%  6.137e+03      100.0%  5.000e+01  75.8% 

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

VecMDot                7 1.0 5.8246e-04 5.3 1.46e+05 1.1 0.0e+00 0.0e+00 7.0e+00  0  0  0  0 11   0  0  0  0 14   472
VecNorm                9 1.0 1.9789e-04 1.0 4.68e+04 1.1 0.0e+00 0.0e+00 9.0e+00  0  0  0  0 14   0  0  0  0 18   447
VecScale               8 1.0 2.6226e-05 1.1 2.08e+04 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  1499
VecCopy                3 1.0 9.0599e-06 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                13 1.0 1.8835e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                2 1.0 3.9132e-03 1.0 1.04e+04 1.1 0.0e+00 0.0e+00 0.0e+00  2  0  0  0  0   2  0  0  0  0     5
VecMAXPY               8 1.0 5.0068e-05 1.2 1.82e+05 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  6869
VecAssemblyBegin       3 1.0 6.3896e-05 1.0 0.00e+00 0.0 2.0e+00 6.1e+03 9.0e+00  0  0  5  5 14   0  0  5  5 18     0
VecAssemblyEnd         3 1.0 1.5974e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin       10 1.0 5.7936e-05 1.2 0.00e+00 0.0 1.8e+01 2.4e+03 0.0e+00  0  0 46 18  0   0  0 46 18  0     0
VecScatterEnd         10 1.0 1.6134e-02704.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  4  0  0  0  0   4  0  0  0  0     0
VecNormalize           8 1.0 2.1720e-04 1.0 6.24e+04 1.1 0.0e+00 0.0e+00 8.0e+00  0  0  0  0 12   0  0  0  0 16   543
MatMult                8 1.0 1.6797e-0222.9 9.78e+05 1.1 1.6e+01 2.3e+03 0.0e+00  4  3 41 15  0   4  3 41 15  0   110
MatSolve               8 1.0 3.6678e-03 1.2 6.28e+06 1.2 0.0e+00 0.0e+00 0.0e+00  2 16  0  0  0   2 16  0  0  0  3165
MatLUFactorNum         1 1.0 2.4829e-02 1.2 3.23e+07 1.2 0.0e+00 0.0e+00 0.0e+00 11 80  0  0  0  11 80  0  0  0  2342
MatILUFactorSym        1 1.0 6.2803e-02 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+00 27  0  0  0  2  27  0  0  0  2     0
MatAssemblyBegin       2 1.0 5.6100e-04 4.6 0.00e+00 0.0 3.0e+00 4.4e+04 4.0e+00  0  0  8 55  6   0  0  8 55  8     0
MatAssemblyEnd         2 1.0 9.5510e-04 1.1 0.00e+00 0.0 4.0e+00 5.8e+02 8.0e+00  0  0 10  1 12   0  0 10  1 16     0
MatGetRowIJ            1 1.0 9.5367e-07 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 4.9114e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  3   0  0  0  0  4     0
MatZeroEntries         3 1.0 1.8907e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog         7 1.0 6.4516e-04 3.8 2.91e+05 1.1 0.0e+00 0.0e+00 7.0e+00  0  1  0  0 11   0  1  0  0 14   853
KSPSetup               2 1.0 7.2002e-05 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 9.7005e-02 1.0 3.99e+07 1.2 1.6e+01 2.3e+03 1.9e+01 47100 41 15 29  47100 41 15 38   746
PCSetUp                2 1.0 8.7902e-02 1.2 3.23e+07 1.2 0.0e+00 0.0e+00 3.0e+00 38 80  0  0  5  38 80  0  0  6   662
PCSetUpOnBlocks        1 1.0 8.7737e-02 1.2 3.23e+07 1.2 0.0e+00 0.0e+00 3.0e+00 38 80  0  0  5  38 80  0  0  6   663
PCApply                8 1.0 3.7963e-03 1.2 6.28e+06 1.2 0.0e+00 0.0e+00 0.0e+00  2 16  0  0  0   2 16  0  0  0  3058
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    23             23       409824     0
         Vec Scatter     3              3         2604     0
           Index Set     8              8        37676     0
   IS L to G Mapping     1              1        11964     0
              Matrix     4              4      5569500     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 1.43051e-06
Average time for zero size MPI_Send(): 5.60284e-06
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

-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 2                                                |
| Time:           Thu Feb  3 12:13:33 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.222039, Active time=0.195229                                             |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 1         0.0010      0.001014    0.0012      0.001223    0.52     0.63     |
|   compute_sparsity()           1         0.0079      0.007901    0.0096      0.009555    4.05     4.89     |
|   create_dof_constraints()     1         0.0005      0.000461    0.0005      0.000461    0.24     0.24     |
|   distribute_dofs()            1         0.0024      0.002436    0.0061      0.006056    1.25     3.10     |
|   dof_indices()                8704      0.0034      0.000000    0.0034      0.000000    1.76     1.76     |
|   prepare_send_list()          1         0.0001      0.000065    0.0001      0.000065    0.03     0.03     |
|   reinit()                     1         0.0035      0.003545    0.0035      0.003545    1.82     1.82     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         2816      0.0030      0.000001    0.0030      0.000001    1.52     1.52     |
|   compute_face_map()           768       0.0008      0.000001    0.0008      0.000001    0.40     0.40     |
|   compute_shape_functions()    2816      0.0017      0.000001    0.0017      0.000001    0.89     0.89     |
|   init_face_shape_functions()  376       0.0009      0.000002    0.0009      0.000002    0.47     0.47     |
|   init_shape_functions()       769       0.0037      0.000005    0.0037      0.000005    1.88     1.88     |
|   inverse_map()                3072      0.0081      0.000003    0.0081      0.000003    4.14     4.14     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             1         0.0079      0.007932    0.0084      0.008407    4.06     4.31     |
|   renumber_nodes_and_elem()    2         0.0004      0.000188    0.0004      0.000188    0.19     0.19     |
|                                                                                                            |
| MeshCommunication                                                                                          |
|   compute_hilbert_indices()    2         0.0130      0.006517    0.0130      0.006517    6.68     6.68     |
|   find_global_indices()        2         0.0018      0.000920    0.0156      0.007821    0.94     8.01     |
|   parallel_sort()              2         0.0006      0.000295    0.0006      0.000305    0.30     0.31     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0018      0.001830    0.0018      0.001830    0.94     0.94     |
|                                                                                                            |
| MetisPartitioner                                                                                           |
|   partition()                  1         0.0056      0.005607    0.0131      0.013127    2.87     6.72     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  6         0.0001      0.000010    0.0001      0.000010    0.03     0.03     |
|   broadcast()                  1         0.0000      0.000005    0.0000      0.000005    0.00     0.00     |
|   gather()                     1         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   max(scalar)                  2         0.0005      0.000248    0.0005      0.000248    0.25     0.25     |
|   max(vector)                  2         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   min(vector)                  2         0.0000      0.000008    0.0000      0.000008    0.01     0.01     |
|   probe()                      10        0.0001      0.000007    0.0001      0.000007    0.04     0.04     |
|   receive()                    10        0.0001      0.000006    0.0001      0.000013    0.03     0.07     |
|   send()                       10        0.0000      0.000002    0.0000      0.000002    0.01     0.01     |
|   send_receive()               14        0.0000      0.000002    0.0002      0.000014    0.02     0.10     |
|   sum()                        9         0.0002      0.000019    0.0002      0.000019    0.09     0.09     |
|   wait()                       10        0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|                                                                                                            |
| Partitioner                                                                                                |
|   set_node_processor_ids()     1         0.0013      0.001305    0.0013      0.001329    0.67     0.68     |
|   set_parent_processor_ids()   1         0.0002      0.000184    0.0002      0.000184    0.09     0.09     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      1         0.0985      0.098453    0.0985      0.098453    50.43    50.43    |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   1         0.0261      0.026113    0.0463      0.046341    13.38    23.74    |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        19419     0.1952                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 2 ./ex5-opt -q 0 -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
