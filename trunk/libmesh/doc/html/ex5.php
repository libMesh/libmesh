<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("examples",$root)?>
 
<div class="content">
<div class = "comment">
Example 5 -- Run-Time Quadrature Rule Selection

<br><br>This is the fifth example program.  It builds on
the previous two examples, and extends the use
of the \p AutoPtr&lt;&gt; as a convenient build method to
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
        #include "steady_system.h"
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
Function prototype, as before.
</div>

<div class ="fragment">
<pre>
        void assemble_poisson(EquationSystems&amp; es,
                              const std::string&amp; system_name);
        
        
        
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
          libMesh::init (argc, argv);
          
          
</pre>
</div>
<div class = "comment">
Braces are used to force object scope, like in example 2   
</div>

<div class ="fragment">
<pre>
          {
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
                
                error();
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
            quad_type = static_cast&lt;QuadratureType&gt;(atoi(argv[2]));
        
        
</pre>
</div>
<div class = "comment">
Independence of dimension has already been shown in
example 4.  For the time being, restrict to 3 dimensions.
</div>

<div class ="fragment">
<pre>
            const unsigned int dim=3;
            
</pre>
</div>
<div class = "comment">
The following is identical to example 4, and therefore
not commented.  Differences are mentioned when present.
</div>

<div class ="fragment">
<pre>
            Mesh mesh (dim);
        
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
            mesh.build_cube (16, 16, 16,
                             -1., 1.,
                             -1., 1.,
                             -1., 1.,
                             HEX8);
            
            mesh.print_info();
            
            EquationSystems equation_systems (mesh);
            
            {
              equation_systems.add_system&lt;SteadySystem&gt; ("Poisson");
              
              equation_systems("Poisson").add_variable("u", FIRST);
        
              equation_systems("Poisson").attach_assemble_function (assemble_poisson);
        
              equation_systems.init();
              
              equation_systems.print_info();
            }
        
            equation_systems("Poisson").solve();
        
</pre>
</div>
<div class = "comment">
"Personalize" the output, with the
number of the quadrature rule appended.
</div>

<div class ="fragment">
<pre>
            std::ostringstream f_name;
            f_name &lt;&lt; "out_" &lt;&lt; quad_type &lt;&lt; ".gmv";
        
            mesh.write_gmv (f_name.str(), equation_systems);
          }
        
        
</pre>
</div>
<div class = "comment">
All done.
</div>

<div class ="fragment">
<pre>
          return libMesh::close ();
        }
        
        
        
        
        void assemble_poisson(EquationSystems&amp; es,
                              const std::string&amp; system_name)
        {
          assert (system_name == "Poisson");
        
          const Mesh&amp; mesh = es.get_mesh();
        
          const unsigned int dim = mesh.mesh_dimension();
        
          FEType fe_type = es("Poisson").get_dof_map().variable_type(0);
        
          
</pre>
</div>
<div class = "comment">
Build a Finite Element object of the specified type.  Since the
\p FEBase::build() member dynamically creates memory we will
store the object as an \p AutoPtr&lt;FEBase&gt;.  Below, the
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
quadrature rule.  Note that a \p AutoPtr&lt;QBase&gt; returns
a QBase* pointer to the object it handles with \p get().  
However, using \p get(), the \p AutoPtr&lt;QBase&gt; \p qrule is 
still in charge of this pointer. I.e., when \p qrule goes 
out of scope, it will safely delete the \p QBase object it 
points to.  This behavior may be overridden using
\p AutoPtr&lt;Xyz&gt;::release(), but is currently not
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
AutoPtr&lt;QBase&gt;  qface (QBase::build(qrule-&gt;type(),
dim-1, 
THIRD));
\endverbatim
And again: using the \p AutoPtr&lt;QBase&gt; relaxes
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
quadrature rule.  Note that a \p AutoPtr&lt;QBase&gt; returns
a \p QBase* pointer to the object it handles with \p get().  
However, using \p get(), the \p AutoPtr&lt;QBase&gt; \p qface is 
still in charge of this pointer. I.e., when \p qface goes 
out of scope, it will safely delete the \p QBase object it 
points to.  This behavior may be overridden using
\p AutoPtr&lt;Xyz&gt;::release(), but is not recommended.
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
          const std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
          
          const std::vector&lt;Point&gt;&amp; q_point = fe-&gt;get_xyz();
          
          const std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
          
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
          
          const DofMap&amp; dof_map = es("Poisson").get_dof_map();
          
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
          const_elem_iterator           el (mesh.elements_begin());
          const const_elem_iterator end_el (mesh.elements_end());
          
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
                      const std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi_face    = fe_face-&gt;get_phi();
                      const std::vector&lt;Real&gt;&amp;               JxW_face    = fe_face-&gt;get_JxW();              
                      const std::vector&lt;Point &gt;&amp;             qface_point = fe_face-&gt;get_xyz();
                      
                      
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
Loop over the face quagrature points for integration.
Note that the \p AutoPtr&lt;QBase&gt; overloaded the operator-&gt;,
so that QBase methods may safely be accessed.  It may
be said: accessing an \p AutoPtr&lt;Xyz&gt; through the
"." operator returns \p AutoPtr methods, while access
through the "-&gt;" operator returns Xyz methods.
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
The element matrix and right-hand-side are now built
for this element.  Add them to the global matrix and
right-hand-side vector.  The \p PetscMatrix::add_matrix()
and \p PetscVector::add_vector() members do this for us.
</div>

<div class ="fragment">
<pre>
              es("Poisson").matrix-&gt;add_matrix (Ke, dof_indices);
              es("Poisson").rhs-&gt;add_vector    (Fe, dof_indices);
              
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

<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  
  #include &lt;iostream&gt;
  #include &lt;sstream&gt; 
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;libmesh.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;mesh.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;steady_system.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;equation_systems.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;fe.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;quadrature.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;quadrature_rules.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;sparse_matrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;numeric_vector.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_matrix.h&quot;</FONT></B>
  #include <FONT COLOR="#BC8F8F"><B>&quot;dense_vector.h&quot;</FONT></B>
  
  #include <FONT COLOR="#BC8F8F"><B>&quot;dof_map.h&quot;</FONT></B>
  
  
  
  
  
  
  
  
  <FONT COLOR="#228B22"><B>void</FONT></B> assemble_poisson(EquationSystems&amp; es,
                        <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name);
  
  
  
  Real exact_solution (<FONT COLOR="#228B22"><B>const</FONT></B> Real x,
  		     <FONT COLOR="#228B22"><B>const</FONT></B> Real y,
  		     <FONT COLOR="#228B22"><B>const</FONT></B> Real z = 0.);
  
  
  QuadratureType quad_type=INVALID_Q_RULE;
  
  
  
  <FONT COLOR="#228B22"><B>int</FONT></B> main (<FONT COLOR="#228B22"><B>int</FONT></B> argc, <FONT COLOR="#228B22"><B>char</FONT></B>** argv)
  {
    
    libMesh::init (argc, argv);
    
    
    {
      <B><FONT COLOR="#A020F0">if</FONT></B> (argc &lt; 3)
        {
  	std::cerr &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Usage: &quot;</FONT></B> &lt;&lt; argv[0] &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; -q n&quot;</FONT></B>
  		  &lt;&lt; std::endl;
  	std::cerr &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;  where n stands for:&quot;</FONT></B> &lt;&lt; std::endl;
  
  	
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> n=0; n&lt;QuadratureRules::num_valid_elem_rules; n++)
  	  std::cerr &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;  &quot;</FONT></B> &lt;&lt; QuadratureRules::valid_elem_rules[n] &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;    &quot;</FONT></B> 
  		    &lt;&lt; QuadratureRules::name(QuadratureRules::valid_elem_rules[n])
  		    &lt;&lt; std::endl;
  	
  	std::cerr &lt;&lt; std::endl;
  	
  	error();
        }
      
      
      <B><FONT COLOR="#A020F0">else</FONT></B> 
        {
  	std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;Running &quot;</FONT></B> &lt;&lt; argv[0];
  	
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>int</FONT></B> i=1; i&lt;argc; i++)
  	  std::cout &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot; &quot;</FONT></B> &lt;&lt; argv[i];
  	
  	std::cout &lt;&lt; std::endl &lt;&lt; std::endl;
        }
      
  
      quad_type = static_cast&lt;QuadratureType&gt;(atoi(argv[2]));
  
  
      <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim=3;
      
      Mesh mesh (dim);
  
      mesh.build_cube (16, 16, 16,
  		     -1., 1.,
  		     -1., 1.,
  		     -1., 1.,
  		     HEX8);
      
      mesh.print_info();
      
      EquationSystems equation_systems (mesh);
      
      {
        equation_systems.add_system&lt;SteadySystem&gt; (<FONT COLOR="#BC8F8F"><B>&quot;Poisson&quot;</FONT></B>);
        
        equation_systems(<FONT COLOR="#BC8F8F"><B>&quot;Poisson&quot;</FONT></B>).add_variable(<FONT COLOR="#BC8F8F"><B>&quot;u&quot;</FONT></B>, FIRST);
  
        equation_systems(<FONT COLOR="#BC8F8F"><B>&quot;Poisson&quot;</FONT></B>).attach_assemble_function (assemble_poisson);
  
        equation_systems.init();
        
        equation_systems.print_info();
      }
  
      equation_systems(<FONT COLOR="#BC8F8F"><B>&quot;Poisson&quot;</FONT></B>).solve();
  
      std::ostringstream f_name;
      f_name &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;out_&quot;</FONT></B> &lt;&lt; quad_type &lt;&lt; <FONT COLOR="#BC8F8F"><B>&quot;.gmv&quot;</FONT></B>;
  
      mesh.write_gmv (f_name.str(), equation_systems);
    }
  
  
    <B><FONT COLOR="#A020F0">return</FONT></B> libMesh::close ();
  }
  
  
  
  
  <FONT COLOR="#228B22"><B>void</FONT></B> assemble_poisson(EquationSystems&amp; es,
                        <FONT COLOR="#228B22"><B>const</FONT></B> std::string&amp; system_name)
  {
    assert (system_name == <FONT COLOR="#BC8F8F"><B>&quot;Poisson&quot;</FONT></B>);
  
    <FONT COLOR="#228B22"><B>const</FONT></B> Mesh&amp; mesh = es.get_mesh();
  
    <FONT COLOR="#228B22"><B>const</FONT></B> <FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> dim = mesh.mesh_dimension();
  
    FEType fe_type = es(<FONT COLOR="#BC8F8F"><B>&quot;Poisson&quot;</FONT></B>).get_dof_map().variable_type(0);
  
    
    AutoPtr&lt;FEBase&gt; fe (FEBase::build(dim, fe_type));
    
    
    AutoPtr&lt;QBase&gt; qrule(QBase::build(quad_type, dim, THIRD));
  
  
    
    fe-&gt;attach_quadrature_rule (qrule.get());
  
    
    AutoPtr&lt;FEBase&gt; fe_face (FEBase::build(dim, fe_type));
    
    
    AutoPtr&lt;QBase&gt;  qface (QBase::build(quad_type,
  				      dim-1, 
  				      THIRD));
  	      
    
    fe_face-&gt;attach_quadrature_rule (qface.get());
  	      
  
    
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
    
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Point&gt;&amp; q_point = fe-&gt;get_xyz();
    
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi = fe-&gt;get_phi();
    
    <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
    
    <FONT COLOR="#228B22"><B>const</FONT></B> DofMap&amp; dof_map = es(<FONT COLOR="#BC8F8F"><B>&quot;Poisson&quot;</FONT></B>).get_dof_map();
    
    DenseMatrix&lt;Number&gt; Ke;
    DenseVector&lt;Number&gt; Fe;
    
    std::vector&lt;<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B>&gt; dof_indices;
    
    
    
    
    
    const_elem_iterator           el (mesh.elements_begin());
    <FONT COLOR="#228B22"><B>const</FONT></B> const_elem_iterator end_el (mesh.elements_end());
    
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {
        <FONT COLOR="#228B22"><B>const</FONT></B> Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        
        fe-&gt;reinit (elem);
        
        Ke.resize (dof_indices.size(),
  		 dof_indices.size());
        
        Fe.resize (dof_indices.size());
        
  
  
        
        <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp&lt;qrule-&gt;n_points(); qp++)
  	{
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi.size(); i++)
  	    <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j&lt;phi.size(); j++)
  	      Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
  	  
  	  
  	  <FONT COLOR="#228B22"><B>const</FONT></B> Real x = q_point[qp](0);
  	  <FONT COLOR="#228B22"><B>const</FONT></B> Real y = q_point[qp](1);
  	  <FONT COLOR="#228B22"><B>const</FONT></B> Real z = q_point[qp](2);
  	  <FONT COLOR="#228B22"><B>const</FONT></B> Real eps = 1.e-3;
  
  	  <FONT COLOR="#228B22"><B>const</FONT></B> Real uxx = (exact_solution(x-eps,y,z) +
  			    exact_solution(x+eps,y,z) +
  			    -2.*exact_solution(x,y,z))/eps/eps;
  	      
  	  <FONT COLOR="#228B22"><B>const</FONT></B> Real uyy = (exact_solution(x,y-eps,z) +
  			    exact_solution(x,y+eps,z) +
  			    -2.*exact_solution(x,y,z))/eps/eps;
  	  
  	  <FONT COLOR="#228B22"><B>const</FONT></B> Real uzz = (exact_solution(x,y,z-eps) +
  			    exact_solution(x,y,z+eps) +
  			    -2.*exact_solution(x,y,z))/eps/eps;
  
  	  <FONT COLOR="#228B22"><B>const</FONT></B> Real fxy = - (uxx + uyy + ((dim==2) ? 0. : uzz));
  	  
  
  	  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi.size(); i++)
  	    Fe(i) += JxW[qp]*fxy*phi[i][qp];	  
  	}
  
  
  
  
        
        
        {
  	<B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> side=0; side&lt;elem-&gt;n_sides(); side++)
  	  <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(side) == NULL)
  	    {	      
  	      <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; phi_face    = fe_face-&gt;get_phi();
  	      <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Real&gt;&amp;               JxW_face    = fe_face-&gt;get_JxW();	      
  	      <FONT COLOR="#228B22"><B>const</FONT></B> std::vector&lt;Point &gt;&amp;             qface_point = fe_face-&gt;get_xyz();
  	      
  	      
  	      fe_face-&gt;reinit(elem, side);
  	      
  	      
  	      <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> qp=0; qp&lt;qface-&gt;n_points(); qp++)
  		{
  		  <FONT COLOR="#228B22"><B>const</FONT></B> Real xf = qface_point[qp](0);
  		  <FONT COLOR="#228B22"><B>const</FONT></B> Real yf = qface_point[qp](1);
  		  <FONT COLOR="#228B22"><B>const</FONT></B> Real zf = qface_point[qp](2);
  		  
  		  <FONT COLOR="#228B22"><B>const</FONT></B> Real penalty = 1.e10;
  		  
  		  <FONT COLOR="#228B22"><B>const</FONT></B> Real value = exact_solution(xf, yf, zf);
  		  
  		  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi_face.size(); i++)
  		    <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> j=0; j&lt;phi_face.size(); j++)
  		      Ke(i,j) += JxW_face[qp]*penalty*phi_face[i][qp]*phi_face[j][qp];
  		  
  		  
  		  <B><FONT COLOR="#A020F0">for</FONT></B> (<FONT COLOR="#228B22"><B>unsigned</FONT></B> <FONT COLOR="#228B22"><B>int</FONT></B> i=0; i&lt;phi_face.size(); i++)
  		    Fe(i) += JxW_face[qp]*penalty*value*phi_face[i][qp];
  		  
  		} <I><FONT COLOR="#B22222">// end face quadrature point loop	  
</FONT></I>  	    } <I><FONT COLOR="#B22222">// end if (elem-&gt;neighbor(side) == NULL)
</FONT></I>        } <I><FONT COLOR="#B22222">// end boundary condition section	  
</FONT></I>        
        
        
        
        
        es(<FONT COLOR="#BC8F8F"><B>&quot;Poisson&quot;</FONT></B>).matrix-&gt;add_matrix (Ke, dof_indices);
        es(<FONT COLOR="#BC8F8F"><B>&quot;Poisson&quot;</FONT></B>).rhs-&gt;add_vector    (Fe, dof_indices);
        
      } <I><FONT COLOR="#B22222">// end of element loop
</FONT></I>    
    
    
    
    <B><FONT COLOR="#A020F0">return</FONT></B>;
  }
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
