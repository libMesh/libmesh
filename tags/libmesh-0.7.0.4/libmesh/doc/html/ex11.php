<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex11",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 11 - Stokes Equations - Systems of Equations</h1>

<br><br>This example shows how a simple, linear system of equations
can be solved in parallel.  The system of equations are the familiar
Stokes equations for low-speed incompressible fluid flow.


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
        #include "equation_systems.h"
        #include "fe.h"
        #include "quadrature_gauss.h"
        #include "dof_map.h"
        #include "sparse_matrix.h"
        #include "numeric_vector.h"
        #include "dense_matrix.h"
        #include "dense_vector.h"
        #include "linear_implicit_system.h"
        
</pre>
</div>
<div class = "comment">
For systems of equations the \p DenseSubMatrix
and \p DenseSubVector provide convenient ways for
assembling the element matrix and vector on a
component-by-component basis.
</div>

<div class ="fragment">
<pre>
        #include "dense_submatrix.h"
        #include "dense_subvector.h"
        
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
Function prototype.  This function will assemble the system
matrix and right-hand-side.
</div>

<div class ="fragment">
<pre>
        void assemble_stokes (EquationSystems& es,
                              const std::string& system_name);
        
</pre>
</div>
<div class = "comment">
The main program.
</div>

<div class ="fragment">
<pre>
        int main (int argc, char** argv)
        {
</pre>
</div>
<div class = "comment">
Initialize libMesh.
</div>

<div class ="fragment">
<pre>
          LibMeshInit init (argc, argv);
        
</pre>
</div>
<div class = "comment">
Skip this 2D example if libMesh was compiled as 1D-only.
</div>

<div class ="fragment">
<pre>
          libmesh_example_assert(2 &lt;= LIBMESH_DIM, "2D support");
            
</pre>
</div>
<div class = "comment">
Create a mesh.
</div>

<div class ="fragment">
<pre>
          Mesh mesh;
            
</pre>
</div>
<div class = "comment">
Use the MeshTools::Generation mesh generator to create a uniform
2D grid on the square [-1,1]^2.  We instruct the mesh generator
to build a mesh of 8x8 \p Quad9 elements.  Building these
higher-order elements allows us to use higher-order
approximation, as in example 3.
</div>

<div class ="fragment">
<pre>
          MeshTools::Generation::build_square (mesh,
                                               15, 15,
                                               0., 1.,
                                               0., 1.,
                                               QUAD9);
          
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
Create an equation systems object.
</div>

<div class ="fragment">
<pre>
          EquationSystems equation_systems (mesh);
          
</pre>
</div>
<div class = "comment">
Declare the system and its variables.
Create a transient system named "Convection-Diffusion"
</div>

<div class ="fragment">
<pre>
          LinearImplicitSystem & system = 
            equation_systems.add_system&lt;LinearImplicitSystem&gt; ("Stokes");
          
</pre>
</div>
<div class = "comment">
Add the variables "u" & "v" to "Stokes".  They
will be approximated using second-order approximation.
</div>

<div class ="fragment">
<pre>
          system.add_variable ("u", SECOND);
          system.add_variable ("v", SECOND);
        
</pre>
</div>
<div class = "comment">
Add the variable "p" to "Stokes". This will
be approximated with a first-order basis,
providing an LBB-stable pressure-velocity pair.
</div>

<div class ="fragment">
<pre>
          system.add_variable ("p", FIRST);
        
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly
function.
</div>

<div class ="fragment">
<pre>
          system.attach_assemble_function (assemble_stokes);
          
</pre>
</div>
<div class = "comment">
Initialize the data structures for the equation system.
</div>

<div class ="fragment">
<pre>
          equation_systems.init ();
        
          equation_systems.parameters.set&lt;unsigned int&gt;("linear solver maximum iterations") = 250;
          equation_systems.parameters.set&lt;Real&gt;        ("linear solver tolerance") = TOLERANCE;
              
</pre>
</div>
<div class = "comment">
Prints information about the system to the screen.
</div>

<div class ="fragment">
<pre>
          equation_systems.print_info();
            
</pre>
</div>
<div class = "comment">
Assemble & solve the linear system,
then write the solution.
</div>

<div class ="fragment">
<pre>
          equation_systems.get_system("Stokes").solve();
        
        #ifdef LIBMESH_HAVE_EXODUS_API
          ExodusII_IO(mesh).write_equation_systems ("out.exd",
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
        
        void assemble_stokes (EquationSystems& es,
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
          libmesh_assert (system_name == "Stokes");
          
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
Get a reference to the Convection-Diffusion system object.
</div>

<div class ="fragment">
<pre>
          LinearImplicitSystem & system =
            es.get_system&lt;LinearImplicitSystem&gt; ("Stokes");
        
</pre>
</div>
<div class = "comment">
Numeric ids corresponding to each variable in the system
</div>

<div class ="fragment">
<pre>
          const unsigned int u_var = system.variable_number ("u");
          const unsigned int v_var = system.variable_number ("v");
          const unsigned int p_var = system.variable_number ("p");
          
</pre>
</div>
<div class = "comment">
Get the Finite Element type for "u".  Note this will be
the same as the type for "v".
</div>

<div class ="fragment">
<pre>
          FEType fe_vel_type = system.variable_type(u_var);
          
</pre>
</div>
<div class = "comment">
Get the Finite Element type for "p".
</div>

<div class ="fragment">
<pre>
          FEType fe_pres_type = system.variable_type(p_var);
        
</pre>
</div>
<div class = "comment">
Build a Finite Element object of the specified type for
the velocity variables.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;FEBase&gt; fe_vel  (FEBase::build(dim, fe_vel_type));
            
</pre>
</div>
<div class = "comment">
Build a Finite Element object of the specified type for
the pressure variables.
</div>

<div class ="fragment">
<pre>
          AutoPtr&lt;FEBase&gt; fe_pres (FEBase::build(dim, fe_pres_type));
          
</pre>
</div>
<div class = "comment">
A Gauss quadrature rule for numerical integration.
Let the \p FEType object decide what order rule is appropriate.
</div>

<div class ="fragment">
<pre>
          QGauss qrule (dim, fe_vel_type.default_quadrature_order());
        
</pre>
</div>
<div class = "comment">
Tell the finite element objects to use our quadrature rule.
</div>

<div class ="fragment">
<pre>
          fe_vel-&gt;attach_quadrature_rule (&qrule);
          fe_pres-&gt;attach_quadrature_rule (&qrule);
          
</pre>
</div>
<div class = "comment">
Here we define some references to cell-specific data that
will be used to assemble the linear system.

<br><br>The element Jacobian * quadrature weight at each integration point.   
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt;& JxW = fe_vel-&gt;get_JxW();
          
</pre>
</div>
<div class = "comment">
The element shape function gradients for the velocity
variables evaluated at the quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi = fe_vel-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
The element shape functions for the pressure variable
evaluated at the quadrature points.
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt;& psi = fe_pres-&gt;get_phi();
          
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
          const DofMap & dof_map = system.get_dof_map();
        
</pre>
</div>
<div class = "comment">
Define data structures to contain the element matrix
and right-hand-side vector contribution.  Following
basic finite element terminology we will denote these
"Ke" and "Fe".
</div>

<div class ="fragment">
<pre>
          DenseMatrix&lt;Number&gt; Ke;
          DenseVector&lt;Number&gt; Fe;
        
          DenseSubMatrix&lt;Number&gt;
            Kuu(Ke), Kuv(Ke), Kup(Ke),
            Kvu(Ke), Kvv(Ke), Kvp(Ke),
            Kpu(Ke), Kpv(Ke), Kpp(Ke);
        
          DenseSubVector&lt;Number&gt;
            Fu(Fe),
            Fv(Fe),
            Fp(Fe);
        
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
          std::vector&lt;unsigned int&gt; dof_indices_u;
          std::vector&lt;unsigned int&gt; dof_indices_v;
          std::vector&lt;unsigned int&gt; dof_indices_p;
          
</pre>
</div>
<div class = "comment">
Now we will loop over all the elements in the mesh that
live on the local processor. We will compute the element
matrix and right-hand-side contribution.  In case users later
modify this program to include refinement, we will be safe and
will only consider the active elements; hence we use a variant of
the \p active_elem_iterator.


<br><br></div>

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
Get the degree of freedom indices for the
current element.  These define where in the global
matrix and right-hand-side this element will
contribute to.
</div>

<div class ="fragment">
<pre>
              dof_map.dof_indices (elem, dof_indices);
              dof_map.dof_indices (elem, dof_indices_u, u_var);
              dof_map.dof_indices (elem, dof_indices_v, v_var);
              dof_map.dof_indices (elem, dof_indices_p, p_var);
        
              const unsigned int n_dofs   = dof_indices.size();
              const unsigned int n_u_dofs = dof_indices_u.size(); 
              const unsigned int n_v_dofs = dof_indices_v.size(); 
              const unsigned int n_p_dofs = dof_indices_p.size();
              
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
              fe_vel-&gt;reinit  (elem);
              fe_pres-&gt;reinit (elem);
        
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
              Ke.resize (n_dofs, n_dofs);
              Fe.resize (n_dofs);
        
</pre>
</div>
<div class = "comment">
Reposition the submatrices...  The idea is this:

<br><br>-           -          -  -
| Kuu Kuv Kup |        | Fu |
Ke = | Kvu Kvv Kvp |;  Fe = | Fv |
| Kpu Kpv Kpp |        | Fp |
-           -          -  -

<br><br>The \p DenseSubMatrix.repostition () member takes the
(row_offset, column_offset, row_size, column_size).

<br><br>Similarly, the \p DenseSubVector.reposition () member
takes the (row_offset, row_size)
</div>

<div class ="fragment">
<pre>
              Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
              Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
              Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);
              
              Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
              Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
              Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);
              
              Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
              Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
              Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);
        
              Fu.reposition (u_var*n_u_dofs, n_u_dofs);
              Fv.reposition (v_var*n_u_dofs, n_v_dofs);
              Fp.reposition (p_var*n_u_dofs, n_p_dofs);
              
</pre>
</div>
<div class = "comment">
Now we will build the element matrix.
</div>

<div class ="fragment">
<pre>
              for (unsigned int qp=0; qp&lt;qrule.n_points(); qp++)
                {
</pre>
</div>
<div class = "comment">
Assemble the u-velocity row
uu coupling
</div>

<div class ="fragment">
<pre>
                  for (unsigned int i=0; i&lt;n_u_dofs; i++)
                    for (unsigned int j=0; j&lt;n_u_dofs; j++)
                      Kuu(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
        
</pre>
</div>
<div class = "comment">
up coupling
</div>

<div class ="fragment">
<pre>
                  for (unsigned int i=0; i&lt;n_u_dofs; i++)
                    for (unsigned int j=0; j&lt;n_p_dofs; j++)
                      Kup(i,j) += -JxW[qp]*psi[j][qp]*dphi[i][qp](0);
        
        
</pre>
</div>
<div class = "comment">
Assemble the v-velocity row
vv coupling
</div>

<div class ="fragment">
<pre>
                  for (unsigned int i=0; i&lt;n_v_dofs; i++)
                    for (unsigned int j=0; j&lt;n_v_dofs; j++)
                      Kvv(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
        
</pre>
</div>
<div class = "comment">
vp coupling
</div>

<div class ="fragment">
<pre>
                  for (unsigned int i=0; i&lt;n_v_dofs; i++)
                    for (unsigned int j=0; j&lt;n_p_dofs; j++)
                      Kvp(i,j) += -JxW[qp]*psi[j][qp]*dphi[i][qp](1);
        
                  
</pre>
</div>
<div class = "comment">
Assemble the pressure row
pu coupling
</div>

<div class ="fragment">
<pre>
                  for (unsigned int i=0; i&lt;n_p_dofs; i++)
                    for (unsigned int j=0; j&lt;n_u_dofs; j++)
                      Kpu(i,j) += -JxW[qp]*psi[i][qp]*dphi[j][qp](0);
        
</pre>
</div>
<div class = "comment">
pv coupling
</div>

<div class ="fragment">
<pre>
                  for (unsigned int i=0; i&lt;n_p_dofs; i++)
                    for (unsigned int j=0; j&lt;n_v_dofs; j++)
                      Kpv(i,j) += -JxW[qp]*psi[i][qp]*dphi[j][qp](1);
                  
                } // end of the quadrature point qp-loop
        
</pre>
</div>
<div class = "comment">
At this point the interior element integration has
been completed.  However, we have not yet addressed
boundary conditions.  For this example we will only
consider simple Dirichlet boundary conditions imposed
via the penalty method. The penalty method used here
is equivalent (for Lagrange basis functions) to lumping
the matrix resulting from the L2 projection penalty
approach introduced in example 3.
</div>

<div class ="fragment">
<pre>
              {
</pre>
</div>
<div class = "comment">
The following loops over the sides of the element.
If the element has no neighbor on a side then that
side MUST live on a boundary of the domain.
</div>

<div class ="fragment">
<pre>
                for (unsigned int s=0; s&lt;elem-&gt;n_sides(); s++)
                  if (elem-&gt;neighbor(s) == NULL)
                    {
                      AutoPtr&lt;Elem&gt; side (elem-&gt;build_side(s));
                                    
</pre>
</div>
<div class = "comment">
Loop over the nodes on the side.
</div>

<div class ="fragment">
<pre>
                      for (unsigned int ns=0; ns&lt;side-&gt;n_nodes(); ns++)
                        {
</pre>
</div>
<div class = "comment">
The location on the boundary of the current
node.
                   

<br><br>const Real xf = side->point(ns)(0);
</div>

<div class ="fragment">
<pre>
                          const Real yf = side-&gt;point(ns)(1);
                          
</pre>
</div>
<div class = "comment">
The penalty value.  \f$ \frac{1}{\epsilon \f$
</div>

<div class ="fragment">
<pre>
                          const Real penalty = 1.e10;
                          
</pre>
</div>
<div class = "comment">
The boundary values.
                   

<br><br>Set u = 1 on the top boundary, 0 everywhere else
</div>

<div class ="fragment">
<pre>
                          const Real u_value = (yf &gt; .99) ? 1. : 0.;
                          
</pre>
</div>
<div class = "comment">
Set v = 0 everywhere
</div>

<div class ="fragment">
<pre>
                          const Real v_value = 0.;
                          
</pre>
</div>
<div class = "comment">
Find the node on the element matching this node on
the side.  That defined where in the element matrix
the boundary condition will be applied.
</div>

<div class ="fragment">
<pre>
                          for (unsigned int n=0; n&lt;elem-&gt;n_nodes(); n++)
                            if (elem-&gt;node(n) == side-&gt;node(ns))
                              {
</pre>
</div>
<div class = "comment">
Matrix contribution.
</div>

<div class ="fragment">
<pre>
                                Kuu(n,n) += penalty;
                                Kvv(n,n) += penalty;
                                            
</pre>
</div>
<div class = "comment">
Right-hand-side contribution.
</div>

<div class ="fragment">
<pre>
                                Fu(n) += penalty*u_value;
                                Fv(n) += penalty*v_value;
                              }
                        } // end face node loop          
                    } // end if (elem-&gt;neighbor(side) == NULL)
              } // end boundary condition section          
              
</pre>
</div>
<div class = "comment">
If this assembly program were to be used on an adaptive mesh,
we would have to apply any hanging node constraint equations.
</div>

<div class ="fragment">
<pre>
              dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);
        
</pre>
</div>
<div class = "comment">
The element matrix and right-hand-side are now built
for this element.  Add them to the global matrix and
right-hand-side vector.  The \p NumericMatrix::add_matrix()
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
That's it.
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
  #include &lt;algorithm&gt;
  #include &lt;math.h&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;linear_implicit_system.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_submatrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_subvector.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_stokes (EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name);
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init (argc, argv);
  
    libmesh_example_assert(2 &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D support&quot;</FONT></B>);
      
    Mesh mesh;
      
    <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_square (mesh,
                                         15, 15,
                                         0., 1.,
                                         0., 1.,
                                         QUAD9);
    
    mesh.print_info();
    
    EquationSystems equation_systems (mesh);
    
    LinearImplicitSystem &amp; system = 
      equation_systems.add_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Stokes&quot;</FONT></B>);
    
    system.add_variable (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, SECOND);
    system.add_variable (<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>, SECOND);
  
    system.add_variable (<B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>, FIRST);
  
    system.attach_assemble_function (assemble_stokes);
    
    equation_systems.init ();
  
    equation_systems.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;linear solver maximum iterations&quot;</FONT></B>) = 250;
    equation_systems.parameters.set&lt;Real&gt;        (<B><FONT COLOR="#BC8F8F">&quot;linear solver tolerance&quot;</FONT></B>) = TOLERANCE;
        
    equation_systems.print_info();
      
    equation_systems.get_system(<B><FONT COLOR="#BC8F8F">&quot;Stokes&quot;</FONT></B>).solve();
  
  #ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO(mesh).write_equation_systems (<B><FONT COLOR="#BC8F8F">&quot;out.exd&quot;</FONT></B>,
                                        equation_systems);
  #endif <I><FONT COLOR="#B22222">// #ifdef LIBMESH_HAVE_EXODUS_API
</FONT></I>  
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_stokes (EquationSystems&amp; es,
                        <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;Stokes&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
    
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
    
    LinearImplicitSystem &amp; system =
      es.get_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;Stokes&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> u_var = system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> v_var = system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;v&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> p_var = system.variable_number (<B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>);
    
    FEType fe_vel_type = system.variable_type(u_var);
    
    FEType fe_pres_type = system.variable_type(p_var);
  
    AutoPtr&lt;FEBase&gt; fe_vel  (FEBase::build(dim, fe_vel_type));
      
    AutoPtr&lt;FEBase&gt; fe_pres (FEBase::build(dim, fe_pres_type));
    
    QGauss qrule (dim, fe_vel_type.default_quadrature_order());
  
    fe_vel-&gt;attach_quadrature_rule (&amp;qrule);
    fe_pres-&gt;attach_quadrature_rule (&amp;qrule);
    
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe_vel-&gt;get_JxW();
    
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe_vel-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp; psi = fe_pres-&gt;get_phi();
    
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap &amp; dof_map = system.get_dof_map();
  
    DenseMatrix&lt;Number&gt; Ke;
    DenseVector&lt;Number&gt; Fe;
  
    DenseSubMatrix&lt;Number&gt;
      Kuu(Ke), Kuv(Ke), Kup(Ke),
      Kvu(Ke), Kvv(Ke), Kvp(Ke),
      Kpu(Ke), Kpv(Ke), Kpp(Ke);
  
    DenseSubVector&lt;Number&gt;
      Fu(Fe),
      Fv(Fe),
      Fp(Fe);
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices_u;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices_v;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices_p;
    
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator       el     = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
    
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
      {    
        <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
        
        dof_map.dof_indices (elem, dof_indices);
        dof_map.dof_indices (elem, dof_indices_u, u_var);
        dof_map.dof_indices (elem, dof_indices_v, v_var);
        dof_map.dof_indices (elem, dof_indices_p, p_var);
  
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_dofs   = dof_indices.size();
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_u_dofs = dof_indices_u.size(); 
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_v_dofs = dof_indices_v.size(); 
        <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_p_dofs = dof_indices_p.size();
        
        fe_vel-&gt;reinit  (elem);
        fe_pres-&gt;reinit (elem);
  
        Ke.resize (n_dofs, n_dofs);
        Fe.resize (n_dofs);
  
        Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
        Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
        Kup.reposition (u_var*n_u_dofs, p_var*n_u_dofs, n_u_dofs, n_p_dofs);
        
        Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
        Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
        Kvp.reposition (v_var*n_v_dofs, p_var*n_v_dofs, n_v_dofs, n_p_dofs);
        
        Kpu.reposition (p_var*n_u_dofs, u_var*n_u_dofs, n_p_dofs, n_u_dofs);
        Kpv.reposition (p_var*n_u_dofs, v_var*n_u_dofs, n_p_dofs, n_v_dofs);
        Kpp.reposition (p_var*n_u_dofs, p_var*n_u_dofs, n_p_dofs, n_p_dofs);
  
        Fu.reposition (u_var*n_u_dofs, n_u_dofs);
        Fv.reposition (v_var*n_u_dofs, n_v_dofs);
        Fp.reposition (p_var*n_u_dofs, n_p_dofs);
        
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
          {
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_u_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_u_dofs; j++)
                Kuu(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_u_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_p_dofs; j++)
                Kup(i,j) += -JxW[qp]*psi[j][qp]*dphi[i][qp](0);
  
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_v_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_v_dofs; j++)
                Kvv(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_v_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_p_dofs; j++)
                Kvp(i,j) += -JxW[qp]*psi[j][qp]*dphi[i][qp](1);
  
            
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_p_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_u_dofs; j++)
                Kpu(i,j) += -JxW[qp]*psi[i][qp]*dphi[j][qp](0);
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_p_dofs; i++)
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_v_dofs; j++)
                Kpv(i,j) += -JxW[qp]*psi[i][qp]*dphi[j][qp](1);
            
          } <I><FONT COLOR="#B22222">// end of the quadrature point qp-loop
</FONT></I>  
        {
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> s=0; s&lt;elem-&gt;n_sides(); s++)
            <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(s) == NULL)
              {
                AutoPtr&lt;Elem&gt; side (elem-&gt;build_side(s));
                              
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> ns=0; ns&lt;side-&gt;n_nodes(); ns++)
                  {
                     
                    <B><FONT COLOR="#228B22">const</FONT></B> Real yf = side-&gt;point(ns)(1);
                    
                    <B><FONT COLOR="#228B22">const</FONT></B> Real penalty = 1.e10;
                    
                     
                    <B><FONT COLOR="#228B22">const</FONT></B> Real u_value = (yf &gt; .99) ? 1. : 0.;
                    
                    <B><FONT COLOR="#228B22">const</FONT></B> Real v_value = 0.;
                    
                    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n=0; n&lt;elem-&gt;n_nodes(); n++)
                      <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;node(n) == side-&gt;node(ns))
                        {
                          Kuu(n,n) += penalty;
                          Kvv(n,n) += penalty;
                                      
                          Fu(n) += penalty*u_value;
                          Fv(n) += penalty*v_value;
                        }
                  } <I><FONT COLOR="#B22222">// end face node loop          
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
* Running Example  mpirun -np 2 ./ex11-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
 Mesh Information:
  mesh_dimension()=2
  spatial_dimension()=3
  n_nodes()=961
    n_local_nodes()=495
  n_elem()=225
    n_local_elem()=112
    n_active_elem()=225
  n_subdomains()=1
  n_processors()=2
  processor_id()=0

 EquationSystems
  n_systems()=1
   System "Stokes"
    Type "LinearImplicit"
    Variables="u" "v" "p" 
    Finite Element Types="LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" "LAGRANGE", "JACOBI_20_00" 
    Infinite Element Mapping="CARTESIAN" "CARTESIAN" "CARTESIAN" 
    Approximation Orders="SECOND", "THIRD" "SECOND", "THIRD" "FIRST", "THIRD" 
    n_dofs()=2178
    n_local_dofs()=1126
    n_constrained_dofs()=0
    n_vectors()=1

************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex11-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Thu Feb  3 12:09:33 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.641e-01      1.01145   1.631e-01
Objects:              3.300e+01      1.00000   3.300e+01
Flops:                3.593e+07      1.15870   3.347e+07  6.695e+07
Flops/sec:            2.190e+08      1.14557   2.051e+08  4.102e+08
MPI Messages:         1.350e+01      1.00000   1.350e+01  2.700e+01
MPI Message Lengths:  3.845e+04      1.01564   2.826e+03  7.631e+04
MPI Reductions:       5.400e+01      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.6313e-01 100.0%  6.6945e+07 100.0%  2.700e+01 100.0%  2.826e+03      100.0%  3.800e+01  70.4% 

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

VecMDot                1 1.0 2.6941e-05 2.3 2.25e+03 1.1 0.0e+00 0.0e+00 1.0e+00  0  0  0  0  2   0  0  0  0  3   162
VecNorm                3 1.0 4.1723e-05 1.0 6.76e+03 1.1 0.0e+00 0.0e+00 3.0e+00  0  0  0  0  6   0  0  0  0  8   313
VecScale               2 1.0 1.4067e-05 1.2 2.25e+03 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0   310
VecCopy                3 1.0 5.9605e-06 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet                 7 1.0 9.2983e-06 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY                2 1.0 3.8803e-03 1.0 4.50e+03 1.1 0.0e+00 0.0e+00 0.0e+00  2  0  0  0  0   2  0  0  0  0     2
VecMAXPY               2 1.0 4.0531e-06 1.0 4.50e+03 1.1 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2149
VecAssemblyBegin       3 1.0 4.0054e-05 1.0 0.00e+00 0.0 2.0e+00 7.7e+02 9.0e+00  0  0  7  2 17   0  0  7  2 24     0
VecAssemblyEnd         3 1.0 1.1921e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin        4 1.0 3.6001e-05 1.4 0.00e+00 0.0 6.0e+00 1.1e+03 0.0e+00  0  0 22  9  0   0  0 22  9  0     0
VecScatterEnd          4 1.0 1.7377e-022208.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  5  0  0  0  0   5  0  0  0  0     0
VecNormalize           2 1.0 6.4135e-05 1.0 6.76e+03 1.1 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  4   0  0  0  0  5   204
MatMult                2 1.0 1.7536e-02109.5 1.73e+05 1.1 4.0e+00 9.3e+02 0.0e+00  5  0 15  5  0   5  0 15  5  0    19
MatSolve               2 1.0 1.0653e-03 1.1 1.66e+06 1.1 0.0e+00 0.0e+00 0.0e+00  1  5  0  0  0   1  5  0  0  0  2965
MatLUFactorNum         1 1.0 2.7160e-02 1.2 3.41e+07 1.2 0.0e+00 0.0e+00 0.0e+00 16 95  0  0  0  16 95  0  0  0  2335
MatILUFactorSym        1 1.0 1.1455e-01 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 1.0e+00 66  0  0  0  2  66  0  0  0  3     0
MatAssemblyBegin       2 1.0 1.7190e-04 1.5 0.00e+00 0.0 3.0e+00 1.5e+04 4.0e+00  0  0 11 59  7   0  0 11 59 11     0
MatAssemblyEnd         2 1.0 4.8995e-04 1.1 0.00e+00 0.0 4.0e+00 2.3e+02 8.0e+00  0  0 15  1 15   0  0 15  1 21     0
MatGetRowIJ            1 1.0 9.5367e-07 0.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         1 1.0 3.6955e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 2.0e+00  0  0  0  0  4   0  0  0  0  5     0
MatZeroEntries         3 1.0 1.3614e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog         1 1.0 4.6968e-05 1.6 4.50e+03 1.1 0.0e+00 0.0e+00 1.0e+00  0  0  0  0  2   0  0  0  0  3   185
KSPSetup               2 1.0 4.7922e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               1 1.0 1.4742e-01 1.0 3.59e+07 1.2 4.0e+00 9.3e+02 7.0e+00 90100 15  5 13  90100 15  5 18   454
PCSetUp                2 1.0 1.4198e-01 1.1 3.41e+07 1.2 0.0e+00 0.0e+00 3.0e+00 82 95  0  0  6  82 95  0  0  8   447
PCSetUpOnBlocks        1 1.0 1.4181e-01 1.1 3.41e+07 1.2 0.0e+00 0.0e+00 3.0e+00 82 95  0  0  6  82 95  0  0  8   447
PCApply                2 1.0 1.1349e-03 1.1 1.66e+06 1.1 0.0e+00 0.0e+00 0.0e+00  1  5  0  0  0   1  5  0  0  0  2783
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec    13             13        91664     0
         Vec Scatter     3              3         2604     0
           Index Set     8              8        18856     0
   IS L to G Mapping     1              1         5504     0
              Matrix     4              4      5558796     0
       Krylov Solver     2              2        18880     0
      Preconditioner     2              2         1408     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 1.57356e-06
Average time for zero size MPI_Send(): 5.48363e-06
#PETSc Option Table entries:
-ksp_right_pc
-log_summary
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

-------------------------------------------------------------------
| Processor id:   0                                                |
| Num Processors: 2                                                |
| Time:           Thu Feb  3 12:09:33 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=0.172747, Active time=0.159178                                             |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 1         0.0001      0.000087    0.0001      0.000114    0.05     0.07     |
|   compute_sparsity()           1         0.0017      0.001717    0.0019      0.001891    1.08     1.19     |
|   create_dof_constraints()     1         0.0001      0.000059    0.0001      0.000059    0.04     0.04     |
|   distribute_dofs()            1         0.0006      0.000605    0.0013      0.001262    0.38     0.79     |
|   dof_indices()                1041      0.0005      0.000000    0.0005      0.000000    0.32     0.32     |
|   prepare_send_list()          1         0.0000      0.000017    0.0000      0.000017    0.01     0.01     |
|   reinit()                     1         0.0006      0.000593    0.0006      0.000593    0.37     0.37     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         224       0.0002      0.000001    0.0002      0.000001    0.15     0.15     |
|   compute_shape_functions()    224       0.0001      0.000000    0.0001      0.000000    0.06     0.06     |
|   init_shape_functions()       2         0.0001      0.000036    0.0001      0.000036    0.05     0.05     |
|                                                                                                            |
| Mesh                                                                                                       |
|   find_neighbors()             1         0.0004      0.000352    0.0010      0.001028    0.22     0.65     |
|   renumber_nodes_and_elem()    2         0.0000      0.000022    0.0000      0.000022    0.03     0.03     |
|                                                                                                            |
| MeshCommunication                                                                                          |
|   compute_hilbert_indices()    2         0.0009      0.000440    0.0009      0.000440    0.55     0.55     |
|   find_global_indices()        2         0.0002      0.000081    0.0014      0.000687    0.10     0.86     |
|   parallel_sort()              2         0.0002      0.000108    0.0002      0.000120    0.14     0.15     |
|                                                                                                            |
| MeshTools::Generation                                                                                      |
|   build_cube()                 1         0.0002      0.000206    0.0002      0.000206    0.13     0.13     |
|                                                                                                            |
| MetisPartitioner                                                                                           |
|   partition()                  1         0.0005      0.000474    0.0010      0.001014    0.30     0.64     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  6         0.0000      0.000007    0.0000      0.000007    0.03     0.03     |
|   broadcast()                  1         0.0000      0.000005    0.0000      0.000005    0.00     0.00     |
|   gather()                     1         0.0000      0.000003    0.0000      0.000003    0.00     0.00     |
|   max(scalar)                  2         0.0007      0.000342    0.0007      0.000342    0.43     0.43     |
|   max(vector)                  2         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   min(vector)                  2         0.0000      0.000007    0.0000      0.000007    0.01     0.01     |
|   probe()                      10        0.0000      0.000003    0.0000      0.000003    0.02     0.02     |
|   receive()                    10        0.0000      0.000004    0.0001      0.000007    0.02     0.05     |
|   send()                       10        0.0000      0.000002    0.0000      0.000002    0.01     0.01     |
|   send_receive()               14        0.0000      0.000003    0.0001      0.000010    0.02     0.09     |
|   sum()                        9         0.0001      0.000010    0.0001      0.000010    0.06     0.06     |
|   wait()                       10        0.0000      0.000001    0.0000      0.000001    0.01     0.01     |
|                                                                                                            |
| Partitioner                                                                                                |
|   set_node_processor_ids()     1         0.0002      0.000174    0.0002      0.000193    0.11     0.12     |
|   set_parent_processor_ids()   1         0.0000      0.000012    0.0000      0.000012    0.01     0.01     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      1         0.1483      0.148347    0.1483      0.148347    93.20    93.20    |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   1         0.0034      0.003350    0.0040      0.004017    2.10     2.52     |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        1589      0.1592                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 2 ./ex11-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
