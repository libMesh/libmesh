<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("miscellaneous_ex5",$root)?>
 
<div class="content">
<a name="comments"></a> 
<br><br><br> <h1> The source file miscellaneous_ex5.C with comments: </h1> 
<div class = "comment">
<h1>Miscellaneous Example 5 - Interior Penalty Discontinuous Galerkin</h1>

<br><br>By Lorenzo Botti

<br><br>This example is based on Adaptivity Example 3, but uses an
Interior Penalty Discontinuous Galerkin formulation.


<br><br></div>

<div class ="fragment">
<pre>
        #include &lt;iostream&gt;
        
</pre>
</div>
<div class = "comment">
LibMesh include files.
</div>

<div class ="fragment">
<pre>
        #include "libmesh/libmesh.h"
        #include "libmesh/mesh.h"
        #include "libmesh/equation_systems.h"
        #include "libmesh/mesh_data.h"
        #include "libmesh/mesh_generation.h"
        #include "libmesh/mesh_modification.h"
        #include "libmesh/elem.h"
        #include "libmesh/transient_system.h"
        #include "libmesh/fe.h"
        #include "libmesh/quadrature_gauss.h"
        #include "libmesh/dof_map.h"
        #include "libmesh/sparse_matrix.h"
        #include "libmesh/dense_matrix.h"
        #include "libmesh/dense_vector.h"
        #include "libmesh/dense_submatrix.h"
        #include "libmesh/dense_subvector.h"
        #include "libmesh/numeric_vector.h"
        #include "libmesh/linear_implicit_system.h"
        #include "libmesh/exodusII_io.h"
        #include "libmesh/fe_interface.h"
        #include "libmesh/getpot.h"
        #include "libmesh/mesh_refinement.h"
        #include "libmesh/error_vector.h"
        #include "libmesh/kelly_error_estimator.h"
        #include "libmesh/discontinuity_measure.h"
        #include "libmesh/string_to_enum.h"
        
        #include "libmesh/exact_solution.h"
</pre>
</div>
<div class = "comment">
#define QORDER TWENTYSIXTH 


<br><br>Bring in everything from the libMesh namespace
</div>

<div class ="fragment">
<pre>
        using namespace libMesh;
        
        Number exact_solution (const Point& p, const Parameters& parameters, const std::string&, const std::string&) 
        {
          const Real x = p(0);
          const Real y = p(1);
          const Real z = p(2);
        
          if (parameters.get&lt;bool&gt;("singularity"))
          {
              Real theta = atan2(y,x);
        
              if (theta &lt; 0)
                 theta += 2. * libMesh::pi;
                          
              return pow(x*x + y*y, 1./3.)*sin(2./3.*theta) + z;
          }
          else
          {
              return cos(x) * exp(y) * (1. - z);
          }
        }
        
</pre>
</div>
<div class = "comment">
We now define the gradient of the exact solution, again being careful
to obtain an angle from atan2 in the correct
quadrant.


<br><br></div>

<div class ="fragment">
<pre>
        Gradient exact_derivative(const Point& p,
                                  const Parameters& parameters,  // es parameters
                                  const std::string&,            // sys_name, not needed
                                  const std::string&)            // unk_name, not needed
        {
</pre>
</div>
<div class = "comment">
Gradient value to be returned.
</div>

<div class ="fragment">
<pre>
          Gradient gradu;
          
</pre>
</div>
<div class = "comment">
x and y coordinates in space
</div>

<div class ="fragment">
<pre>
          const Real x = p(0);
          const Real y = p(1);
          const Real z = p(2);
        
          if (parameters.get&lt;bool&gt;("singularity"))
          {
              libmesh_assert_not_equal_to (x, 0.);
        
</pre>
</div>
<div class = "comment">
For convenience...
</div>

<div class ="fragment">
<pre>
              const Real tt = 2./3.;
              const Real ot = 1./3.;
              
</pre>
</div>
<div class = "comment">
The value of the radius, squared
</div>

<div class ="fragment">
<pre>
              const Real r2 = x*x + y*y;
        
</pre>
</div>
<div class = "comment">
The boundary value, given by the exact solution,
u_exact = r^(2/3)*sin(2*theta/3).
</div>

<div class ="fragment">
<pre>
              Real theta = atan2(y,x);
              
</pre>
</div>
<div class = "comment">
Make sure 0 <= theta <= 2*pi
</div>

<div class ="fragment">
<pre>
              if (theta &lt; 0)
                theta += 2. * libMesh::pi;
        
</pre>
</div>
<div class = "comment">
du/dx
</div>

<div class ="fragment">
<pre>
              gradu(0) = tt*x*pow(r2,-tt)*sin(tt*theta) - pow(r2,ot)*cos(tt*theta)*tt/(1.+y*y/x/x)*y/x/x;
              gradu(1) = tt*y*pow(r2,-tt)*sin(tt*theta) + pow(r2,ot)*cos(tt*theta)*tt/(1.+y*y/x/x)*1./x;
              gradu(2) = 1.;
          }
          else
          {
              gradu(0) = -sin(x) * exp(y) * (1. - z);
              gradu(1) = cos(x) * exp(y) * (1. - z);
              gradu(2) = -cos(x) * exp(y);
          }
          return gradu;
        }
        
</pre>
</div>
<div class = "comment">
We now define the matrix assembly function for the
Laplace system.  We need to first compute element volume
matrices, and then take into account the boundary 
conditions and the flux integrals, which will be handled
via an interior penalty method.
</div>

<div class ="fragment">
<pre>
        void assemble_ellipticdg(EquationSystems& es, const std::string& system_name)
        {
          std::cout&lt;&lt;" assembling elliptic dg system... ";
          std::cout.flush();
        
</pre>
</div>
<div class = "comment">
It is a good idea to make sure we are assembling
the proper system.
</div>

<div class ="fragment">
<pre>
          libmesh_assert_equal_to (system_name, "EllipticDG");
          
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
          LinearImplicitSystem & ellipticdg_system = es.get_system&lt;LinearImplicitSystem&gt; ("EllipticDG");
</pre>
</div>
<div class = "comment">
Get some parameters that we need during assembly
</div>

<div class ="fragment">
<pre>
          const Real penalty = es.parameters.get&lt;Real&gt; ("penalty");
          std::string refinement_type = es.parameters.get&lt;std::string&gt; ("refinement");
        
</pre>
</div>
<div class = "comment">
A reference to the \p DofMap object for this system.  The \p DofMap
object handles the index translation from node and element numbers
to degree of freedom numbers.  We will talk more about the \p DofMap
</div>

<div class ="fragment">
<pre>
          const DofMap & dof_map = ellipticdg_system.get_dof_map();
          
</pre>
</div>
<div class = "comment">
Get a constant reference to the Finite Element type
for the first (and only) variable in the system.
</div>

<div class ="fragment">
<pre>
          FEType fe_type = ellipticdg_system.variable_type(0);
        
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
          AutoPtr&lt;FEBase&gt; fe  (FEBase::build(dim, fe_type));
          AutoPtr&lt;FEBase&gt; fe_elem_face(FEBase::build(dim, fe_type));
          AutoPtr&lt;FEBase&gt; fe_neighbor_face(FEBase::build(dim, fe_type));
        
</pre>
</div>
<div class = "comment">
Quadrature rules for numerical integration.
</div>

<div class ="fragment">
<pre>
        #ifdef QORDER 
          QGauss qrule (dim, QORDER);
        #else
          QGauss qrule (dim, fe_type.default_quadrature_order());
        #endif
          fe-&gt;attach_quadrature_rule (&qrule);
          
        #ifdef QORDER 
          QGauss qface(dim-1, QORDER);
        #else
          QGauss qface(dim-1, fe_type.default_quadrature_order());
        #endif
          
</pre>
</div>
<div class = "comment">
Tell the finite element object to use our quadrature rule.
</div>

<div class ="fragment">
<pre>
          fe_elem_face-&gt;attach_quadrature_rule(&qface);
          fe_neighbor_face-&gt;attach_quadrature_rule(&qface);
        
</pre>
</div>
<div class = "comment">
Here we define some references to cell-specific data that
will be used to assemble the linear system.
Data for interior volume integrals
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;Real&gt;& JxW = fe-&gt;get_JxW();
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi = fe-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
Data for surface integrals on the element boundary
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt;&  phi_face = fe_elem_face-&gt;get_phi();
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi_face = fe_elem_face-&gt;get_dphi();
          const std::vector&lt;Real&gt;& JxW_face = fe_elem_face-&gt;get_JxW();
          const std::vector&lt;Point&gt;& qface_normals = fe_elem_face-&gt;get_normals();
          const std::vector&lt;Point&gt;& qface_points = fe_elem_face-&gt;get_xyz();
            
</pre>
</div>
<div class = "comment">
Data for surface integrals on the neighbor boundary
</div>

<div class ="fragment">
<pre>
          const std::vector&lt;std::vector&lt;Real&gt; &gt;&  phi_neighbor_face = fe_neighbor_face-&gt;get_phi();
          const std::vector&lt;std::vector&lt;RealGradient&gt; &gt;& dphi_neighbor_face = fe_neighbor_face-&gt;get_dphi();
        
</pre>
</div>
<div class = "comment">
Define data structures to contain the element interior matrix
and right-hand-side vector contribution.  Following
basic finite element terminology we will denote these
"Ke" and "Fe". 
</div>

<div class ="fragment">
<pre>
          DenseMatrix&lt;Number&gt; Ke;
          DenseVector&lt;Number&gt; Fe;
        
</pre>
</div>
<div class = "comment">
Data structures to contain the element and neighbor boundary matrix
contribution. This matrices will do the coupling beetwen the dofs of
the element and those of his neighbors.
Ken: matrix coupling elem and neighbor dofs
</div>

<div class ="fragment">
<pre>
          DenseMatrix&lt;Number&gt; Kne;
          DenseMatrix&lt;Number&gt; Ken;
          DenseMatrix&lt;Number&gt; Kee;
          DenseMatrix&lt;Number&gt; Knn;
          
</pre>
</div>
<div class = "comment">
This vector will hold the degree of freedom indices for
the element.  These define where in the global system
the element degrees of freedom get mapped.
</div>

<div class ="fragment">
<pre>
          std::vector&lt;dof_id_type&gt; dof_indices;
        
</pre>
</div>
<div class = "comment">
Now we will loop over all the elements in the mesh.  We will
compute first the element interior matrix and right-hand-side contribution
and then the element and neighbors boundary matrix contributions.  
</div>

<div class ="fragment">
<pre>
          MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
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
            const unsigned int n_dofs   = dof_indices.size();
              
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
            fe-&gt;reinit  (elem);
        
</pre>
</div>
<div class = "comment">
Zero the element matrix and right-hand side before
summing them.  We use the resize member here because
the number of degrees of freedom might have changed from
the last element. 
</div>

<div class ="fragment">
<pre>
            Ke.resize (n_dofs, n_dofs);
            Fe.resize (n_dofs);
         
</pre>
</div>
<div class = "comment">
Now we will build the element interior matrix.  This involves
a double loop to integrate the test funcions (i) against
the trial functions (j).
</div>

<div class ="fragment">
<pre>
            for (unsigned int qp=0; qp&lt;qrule.n_points(); qp++)
            {
              for (unsigned int i=0; i&lt;n_dofs; i++)
              {
                for (unsigned int j=0; j&lt;n_dofs; j++)
                {
                   Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
                }
              }
            }
           
</pre>
</div>
<div class = "comment">
Now we adress boundary conditions. 
We consider Dirichlet bc imposed via the interior penalty method
The following loops over the sides of the element.
If the element has no neighbor on a side then that
side MUST live on a boundary of the domain.
</div>

<div class ="fragment">
<pre>
            for (unsigned int side=0; side&lt;elem-&gt;n_sides(); side++)
            {
              if (elem-&gt;neighbor(side) == NULL)
              {
</pre>
</div>
<div class = "comment">
Pointer to the element face
</div>

<div class ="fragment">
<pre>
                fe_elem_face-&gt;reinit(elem, side);
                
                AutoPtr&lt;Elem&gt; elem_side (elem-&gt;build_side(side));
</pre>
</div>
<div class = "comment">
h elemet dimension to compute the interior penalty penalty parameter
</div>

<div class ="fragment">
<pre>
                const unsigned int elem_b_order = static_cast&lt;unsigned int&gt; (fe_elem_face-&gt;get_order());
                const double h_elem = elem-&gt;volume()/elem_side-&gt;volume() * 1./pow(elem_b_order, 2.);
                
                for (unsigned int qp=0; qp&lt;qface.n_points(); qp++)
                {
                  Number bc_value = exact_solution(qface_points[qp], es.parameters,"null","void");
                  for (unsigned int i=0; i&lt;n_dofs; i++)
                  {
</pre>
</div>
<div class = "comment">
Matrix contribution
</div>

<div class ="fragment">
<pre>
                     for (unsigned int j=0; j&lt;n_dofs; j++)
                     { 
                        Ke(i,j) += JxW_face[qp] * penalty/h_elem * phi_face[i][qp] * phi_face[j][qp]; // stability 
                        Ke(i,j) -= JxW_face[qp] * (phi_face[i][qp] * (dphi_face[j][qp]*qface_normals[qp]) + phi_face[j][qp] * (dphi_face[i][qp]*qface_normals[qp])); // consistency
                     }
</pre>
</div>
<div class = "comment">
RHS contribution
</div>

<div class ="fragment">
<pre>
                     Fe(i) += JxW_face[qp] * bc_value * penalty/h_elem * phi_face[i][qp];      // stability
                     Fe(i) -= JxW_face[qp] * dphi_face[i][qp] * (bc_value*qface_normals[qp]);     // consistency
                  }
                }
              }
</pre>
</div>
<div class = "comment">
If the element is not on a boundary of the domain 
we loop over his neighbors to compute the element 
and neighbor boundary matrix contributions 
</div>

<div class ="fragment">
<pre>
              else 
              {
</pre>
</div>
<div class = "comment">
Store a pointer to the neighbor we are currently
working on.  
</div>

<div class ="fragment">
<pre>
                const Elem* neighbor = elem-&gt;neighbor(side);
                 
</pre>
</div>
<div class = "comment">
Get the global id of the element and the neighbor
</div>

<div class ="fragment">
<pre>
                const unsigned int elem_id = elem-&gt;id();
                const unsigned int neighbor_id = neighbor-&gt;id();
                
</pre>
</div>
<div class = "comment">
If the neighbor has the same h level and is active
perform integration only if our global id is bigger than our neighbor id.
We don't want to compute twice the same contributions.
If the neighbor has a different h level perform integration
only if the neighbor is at a lower level. 
</div>

<div class ="fragment">
<pre>
                if ((neighbor-&gt;active() && (neighbor-&gt;level() == elem-&gt;level()) && (elem_id &lt; neighbor_id)) || (neighbor-&gt;level() &lt; elem-&gt;level()))
                {
</pre>
</div>
<div class = "comment">
Pointer to the element side
</div>

<div class ="fragment">
<pre>
                  AutoPtr&lt;Elem&gt; elem_side (elem-&gt;build_side(side));
                  
</pre>
</div>
<div class = "comment">
h dimension to compute the interior penalty penalty parameter
</div>

<div class ="fragment">
<pre>
                  const unsigned int elem_b_order = static_cast&lt;unsigned int&gt;(fe_elem_face-&gt;get_order());
                  const unsigned int neighbor_b_order = static_cast&lt;unsigned int&gt;(fe_neighbor_face-&gt;get_order());
                  const double side_order = (elem_b_order + neighbor_b_order)/2.;
                  const double h_elem = (elem-&gt;volume()/elem_side-&gt;volume()) * 1./pow(side_order,2.);
        
</pre>
</div>
<div class = "comment">
The quadrature point locations on the neighbor side 
</div>

<div class ="fragment">
<pre>
                  std::vector&lt;Point&gt; qface_neighbor_point;
        
</pre>
</div>
<div class = "comment">
The quadrature point locations on the element side
</div>

<div class ="fragment">
<pre>
                  std::vector&lt;Point &gt; qface_point;
                 
</pre>
</div>
<div class = "comment">
Reinitialize shape functions on the element side
</div>

<div class ="fragment">
<pre>
                  fe_elem_face-&gt;reinit(elem, side);
        
</pre>
</div>
<div class = "comment">
Get the physical locations of the element quadrature points
</div>

<div class ="fragment">
<pre>
                  qface_point = fe_elem_face-&gt;get_xyz();
        
</pre>
</div>
<div class = "comment">
Find their locations on the neighbor 
</div>

<div class ="fragment">
<pre>
                  unsigned int side_neighbor = neighbor-&gt;which_neighbor_am_i(elem);
                  if (refinement_type == "p")
                    fe_neighbor_face-&gt;side_map (neighbor, elem_side.get(), side_neighbor, qface.get_points(), qface_neighbor_point);
        	  else
                    FEInterface::inverse_map (elem-&gt;dim(), fe-&gt;get_fe_type(), neighbor, qface_point, qface_neighbor_point);
</pre>
</div>
<div class = "comment">
Calculate the neighbor element shape functions at those locations
</div>

<div class ="fragment">
<pre>
                  fe_neighbor_face-&gt;reinit(neighbor, &qface_neighbor_point);
                  
</pre>
</div>
<div class = "comment">
Get the degree of freedom indices for the
neighbor.  These define where in the global
matrix this neighbor will contribute to.
</div>

<div class ="fragment">
<pre>
                  std::vector&lt;dof_id_type&gt; neighbor_dof_indices;
                  dof_map.dof_indices (neighbor, neighbor_dof_indices);
                  const unsigned int n_neighbor_dofs = neighbor_dof_indices.size();
        
</pre>
</div>
<div class = "comment">
Zero the element and neighbor side matrix before
summing them.  We use the resize member here because
the number of degrees of freedom might have changed from
the last element or neighbor. 
Note that Kne and Ken are not square matrices if neighbor 
and element have a different p level
</div>

<div class ="fragment">
<pre>
                  Kne.resize (n_neighbor_dofs, n_dofs);
                  Ken.resize (n_dofs, n_neighbor_dofs);
                  Kee.resize (n_dofs, n_dofs);
                  Knn.resize (n_neighbor_dofs, n_neighbor_dofs);
        
</pre>
</div>
<div class = "comment">
Now we will build the element and neighbor
boundary matrices.  This involves
a double loop to integrate the test funcions
(i) against the trial functions (j).
</div>

<div class ="fragment">
<pre>
                  for (unsigned int qp=0; qp&lt;qface.n_points(); qp++)
                  {
</pre>
</div>
<div class = "comment">
Kee Matrix. Integrate the element test function i
against the element test function j
</div>

<div class ="fragment">
<pre>
                    for (unsigned int i=0; i&lt;n_dofs; i++)          
                    {
                      for (unsigned int j=0; j&lt;n_dofs; j++)
                      {
                         Kee(i,j) -= 0.5 * JxW_face[qp] * (phi_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) + phi_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp])); // consistency
                         Kee(i,j) += JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*phi_face[i][qp];  // stability
                      }
                    }
        
</pre>
</div>
<div class = "comment">
Knn Matrix. Integrate the neighbor test function i
against the neighbor test function j
</div>

<div class ="fragment">
<pre>
                    for (unsigned int i=0; i&lt;n_neighbor_dofs; i++) 
                    {
                      for (unsigned int j=0; j&lt;n_neighbor_dofs; j++)
                      {
                         Knn(i,j) += 0.5 * JxW_face[qp] * (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]) + phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp])); // consistency
                         Knn(i,j) += JxW_face[qp] * penalty/h_elem * phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp]; // stability
                      }
                    }
        
</pre>
</div>
<div class = "comment">
Kne Matrix. Integrate the neighbor test function i
against the element test function j
</div>

<div class ="fragment">
<pre>
                    for (unsigned int i=0; i&lt;n_neighbor_dofs; i++) 
                    {
                      for (unsigned int j=0; j&lt;n_dofs; j++)
                      {
                         Kne(i,j) += 0.5 * JxW_face[qp] * (phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]) - phi_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp])); // consistency
                         Kne(i,j) -= JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*phi_neighbor_face[i][qp];  // stability
                      }
                    }
                    
</pre>
</div>
<div class = "comment">
Ken Matrix. Integrate the element test function i
against the neighbor test function j
</div>

<div class ="fragment">
<pre>
                    for (unsigned int i=0; i&lt;n_dofs; i++)         
                    {
                      for (unsigned int j=0; j&lt;n_neighbor_dofs; j++)
                      {
                         Ken(i,j) += 0.5 * JxW_face[qp] * (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) - phi_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));  // consistency
                         Ken(i,j) -= JxW_face[qp] * penalty/h_elem * phi_face[i][qp]*phi_neighbor_face[j][qp];  // stability
                      }
                    }
                  }
              
</pre>
</div>
<div class = "comment">
The element and neighbor boundary matrix are now built
for this side.  Add them to the global matrix 
The \p SparseMatrix::add_matrix() members do this for us.
</div>

<div class ="fragment">
<pre>
                  ellipticdg_system.matrix-&gt;add_matrix(Kne,neighbor_dof_indices,dof_indices);
                  ellipticdg_system.matrix-&gt;add_matrix(Ken,dof_indices,neighbor_dof_indices);
                  ellipticdg_system.matrix-&gt;add_matrix(Kee,dof_indices);
                  ellipticdg_system.matrix-&gt;add_matrix(Knn,neighbor_dof_indices);
                }
              }
            }
</pre>
</div>
<div class = "comment">
The element interior matrix and right-hand-side are now built
for this element.  Add them to the global matrix and
right-hand-side vector.  The \p SparseMatrix::add_matrix()
and \p NumericVector::add_vector() members do this for us.
</div>

<div class ="fragment">
<pre>
            ellipticdg_system.matrix-&gt;add_matrix(Ke, dof_indices);
            ellipticdg_system.rhs-&gt;add_vector(Fe, dof_indices);
          }
        
          std::cout &lt;&lt; "done" &lt;&lt; std::endl;
        }
        
        
        
        int main (int argc, char** argv)
        {
          LibMeshInit init(argc, argv);
        
</pre>
</div>
<div class = "comment">
Parse the input file
</div>

<div class ="fragment">
<pre>
          GetPot input_file("miscellaneous_ex5.in");
        
</pre>
</div>
<div class = "comment">
Read in parameters from the input file
</div>

<div class ="fragment">
<pre>
          const unsigned int adaptive_refinement_steps = input_file("max_adaptive_r_steps",3);
          const unsigned int uniform_refinement_steps  = input_file("uniform_h_r_steps",3);
          const Real refine_fraction                   = input_file("refine_fraction",0.5);
          const Real coarsen_fraction                  = input_file("coarsen_fraction",0.);
          const unsigned int max_h_level               = input_file("max_h_level", 10);
          const std::string refinement_type            = input_file("refinement_type","p");
          Order p_order                                = static_cast&lt;Order&gt;(input_file("p_order",1));
          const std::string element_type               = input_file("element_type", "tensor");
          const Real penalty                           = input_file("ip_penalty", 10.);
          const bool singularity                       = input_file("singularity", true);
          const unsigned int dim                       = input_file("dimension", 3);
        
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
Create or read the mesh
</div>

<div class ="fragment">
<pre>
          Mesh mesh;
        
          if (dim == 1)
            MeshTools::Generation::build_line(mesh,1,-1.,0.);
          else if (dim == 2)
            mesh.read("lshaped.xda");
          else
            mesh.read ("lshaped3D.xda");
        
</pre>
</div>
<div class = "comment">
Use triangles if the config file says so
</div>

<div class ="fragment">
<pre>
          if (element_type == "simplex")
            MeshTools::Modification::all_tri(mesh);
        
</pre>
</div>
<div class = "comment">
Mesh Refinement object 
</div>

<div class ="fragment">
<pre>
          MeshRefinement mesh_refinement(mesh);
          mesh_refinement.refine_fraction() = refine_fraction;
          mesh_refinement.coarsen_fraction() = coarsen_fraction;
          mesh_refinement.max_h_level() = max_h_level;
</pre>
</div>
<div class = "comment">
Do uniform refinement  
</div>

<div class ="fragment">
<pre>
          for (unsigned int rstep=0; rstep&lt;uniform_refinement_steps; rstep++)
          {
             mesh_refinement.uniformly_refine(1);
          }
         
</pre>
</div>
<div class = "comment">
Crate an equation system object 
</div>

<div class ="fragment">
<pre>
          EquationSystems equation_system (mesh);
        
</pre>
</div>
<div class = "comment">
Set parameters for the equation system and the solver
</div>

<div class ="fragment">
<pre>
          equation_system.parameters.set&lt;Real&gt;("linear solver tolerance") = TOLERANCE * TOLERANCE;
          equation_system.parameters.set&lt;unsigned int&gt;("linear solver maximum iterations") = 1000;
          equation_system.parameters.set&lt;Real&gt;("penalty") = penalty;
          equation_system.parameters.set&lt;bool&gt;("singularity") = singularity;
          equation_system.parameters.set&lt;std::string&gt;("refinement") = refinement_type;
        
</pre>
</div>
<div class = "comment">
Create a system named ellipticdg
</div>

<div class ="fragment">
<pre>
          LinearImplicitSystem& ellipticdg_system = equation_system.add_system&lt;LinearImplicitSystem&gt; ("EllipticDG");
         
</pre>
</div>
<div class = "comment">
Add a variable "u" to "ellipticdg" using the p_order specified in the config file
</div>

<div class ="fragment">
<pre>
          if( on_command_line( "element_type" ) )
            {
              std::string fe_str = command_line_value( std::string("element_type"),
        					       std::string("MONOMIAL") );
              if( fe_str != "MONOMIAL" || fe_str != "XYZ" )
        	{
        	  std::cerr &lt;&lt; "Error: This example must be run with MONOMIAL or XYZ element types." &lt;&lt; std::endl;
        	  libmesh_error();
        	}
              ellipticdg_system.add_variable ("u", p_order, Utility::string_to_enum&lt;FEFamily&gt;(fe_str) );
             } 
          else
            ellipticdg_system.add_variable ("u", p_order, MONOMIAL);
          
</pre>
</div>
<div class = "comment">
Give the system a pointer to the matrix assembly function 
</div>

<div class ="fragment">
<pre>
          ellipticdg_system.attach_assemble_function (assemble_ellipticdg);
        
</pre>
</div>
<div class = "comment">
Initialize the data structures for the equation system
</div>

<div class ="fragment">
<pre>
          equation_system.init();
        
</pre>
</div>
<div class = "comment">
Construct ExactSolution object and attach solution functions
</div>

<div class ="fragment">
<pre>
          ExactSolution exact_sol(equation_system);
          exact_sol.attach_exact_value(exact_solution);
          exact_sol.attach_exact_deriv(exact_derivative);
        
</pre>
</div>
<div class = "comment">
A refinement loop.
</div>

<div class ="fragment">
<pre>
          for (unsigned int rstep=0; rstep&lt;adaptive_refinement_steps; ++rstep)
          {
            std::cout &lt;&lt; "  Beginning Solve " &lt;&lt; rstep &lt;&lt; std::endl;
            std::cout &lt;&lt; "Number of elements: " &lt;&lt; mesh.n_elem() &lt;&lt; std::endl;
               
</pre>
</div>
<div class = "comment">
Solve the system
</div>

<div class ="fragment">
<pre>
            ellipticdg_system.solve();
                    
            std::cout &lt;&lt; "System has: " &lt;&lt; equation_system.n_active_dofs()
                      &lt;&lt; " degrees of freedom."
                      &lt;&lt; std::endl;
        
            std::cout &lt;&lt; "Linear solver converged at step: "
                      &lt;&lt; ellipticdg_system.n_linear_iterations()
                      &lt;&lt; ", final residual: "
                      &lt;&lt; ellipticdg_system.final_linear_residual()
                      &lt;&lt; std::endl;
         
</pre>
</div>
<div class = "comment">
Compute the error
</div>

<div class ="fragment">
<pre>
            exact_sol.compute_error("EllipticDG", "u");
        
</pre>
</div>
<div class = "comment">
Print out the error values
</div>

<div class ="fragment">
<pre>
            std::cout &lt;&lt; "L2-Error is: "
                      &lt;&lt; exact_sol.l2_error("EllipticDG", "u")
                      &lt;&lt; std::endl;
           
</pre>
</div>
<div class = "comment">
Possibly refine the mesh 
</div>

<div class ="fragment">
<pre>
            if (rstep+1 &lt; adaptive_refinement_steps)
            {
</pre>
</div>
<div class = "comment">
The ErrorVector is a particular StatisticsVector
for computing error information on a finite element mesh.
</div>

<div class ="fragment">
<pre>
              ErrorVector error;
        
</pre>
</div>
<div class = "comment">
The discontinuity error estimator 
evaluate the jump of the solution 
on elements faces 
</div>

<div class ="fragment">
<pre>
              DiscontinuityMeasure error_estimator;
              error_estimator.estimate_error(ellipticdg_system,error);
              
</pre>
</div>
<div class = "comment">
Take the error in error and decide which elements will be coarsened or refined  
</div>

<div class ="fragment">
<pre>
              mesh_refinement.flag_elements_by_error_fraction(error);
              if (refinement_type == "p")
                mesh_refinement.switch_h_to_p_refinement();
              if (refinement_type == "hp")
                mesh_refinement.add_p_to_h_refinement();
</pre>
</div>
<div class = "comment">
Refine and coarsen the flagged elements
</div>

<div class ="fragment">
<pre>
              mesh_refinement.refine_and_coarsen_elements();
              equation_system.reinit();
            }
          }
        
</pre>
</div>
<div class = "comment">
Write out the solution
After solving the system write the solution
to a ExodusII-formatted plot file.
</div>

<div class ="fragment">
<pre>
        #ifdef LIBMESH_HAVE_EXODUS_API
          ExodusII_IO (mesh).write_discontinuous_exodusII("lshaped_dg.e",equation_system);
        #endif
        
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

<a name="nocomments"></a> 
<br><br><br> <h1> The source file miscellaneous_ex5.C without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_data.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_modification.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/elem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/transient_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_submatrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/dense_subvector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/fe_interface.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/kelly_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/discontinuity_measure.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/string_to_enum.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh/exact_solution.h&quot;</FONT></B>
  
  using namespace libMesh;
  
  Number exact_solution (<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p, <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp; parameters, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;) 
  {
    <B><FONT COLOR="#228B22">const</FONT></B> Real x = p(0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real y = p(1);
    <B><FONT COLOR="#228B22">const</FONT></B> Real z = p(2);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (parameters.get&lt;<B><FONT COLOR="#228B22">bool</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;singularity&quot;</FONT></B>))
    {
        Real theta = atan2(y,x);
  
        <B><FONT COLOR="#A020F0">if</FONT></B> (theta &lt; 0)
           theta += 2. * libMesh::pi;
                    
        <B><FONT COLOR="#A020F0">return</FONT></B> pow(x*x + y*y, 1./3.)*sin(2./3.*theta) + z;
    }
    <B><FONT COLOR="#A020F0">else</FONT></B>
    {
        <B><FONT COLOR="#A020F0">return</FONT></B> cos(x) * exp(y) * (1. - z);
    }
  }
  
  
  Gradient exact_derivative(<B><FONT COLOR="#228B22">const</FONT></B> Point&amp; p,
                            <B><FONT COLOR="#228B22">const</FONT></B> Parameters&amp; parameters,  <I><FONT COLOR="#B22222">// es parameters
</FONT></I>                            <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;,            <I><FONT COLOR="#B22222">// sys_name, not needed
</FONT></I>                            <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp;)            <I><FONT COLOR="#B22222">// unk_name, not needed
</FONT></I>  {
    Gradient gradu;
    
    <B><FONT COLOR="#228B22">const</FONT></B> Real x = p(0);
    <B><FONT COLOR="#228B22">const</FONT></B> Real y = p(1);
    <B><FONT COLOR="#228B22">const</FONT></B> Real z = p(2);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (parameters.get&lt;<B><FONT COLOR="#228B22">bool</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;singularity&quot;</FONT></B>))
    {
        libmesh_assert_not_equal_to (x, 0.);
  
        <B><FONT COLOR="#228B22">const</FONT></B> Real tt = 2./3.;
        <B><FONT COLOR="#228B22">const</FONT></B> Real ot = 1./3.;
        
        <B><FONT COLOR="#228B22">const</FONT></B> Real r2 = x*x + y*y;
  
        Real theta = atan2(y,x);
        
        <B><FONT COLOR="#A020F0">if</FONT></B> (theta &lt; 0)
          theta += 2. * libMesh::pi;
  
        gradu(0) = tt*x*pow(r2,-tt)*sin(tt*theta) - pow(r2,ot)*cos(tt*theta)*tt/(1.+y*y/x/x)*y/x/x;
        gradu(1) = tt*y*pow(r2,-tt)*sin(tt*theta) + pow(r2,ot)*cos(tt*theta)*tt/(1.+y*y/x/x)*1./x;
        gradu(2) = 1.;
    }
    <B><FONT COLOR="#A020F0">else</FONT></B>
    {
        gradu(0) = -sin(x) * exp(y) * (1. - z);
        gradu(1) = cos(x) * exp(y) * (1. - z);
        gradu(2) = -cos(x) * exp(y);
    }
    <B><FONT COLOR="#A020F0">return</FONT></B> gradu;
  }
  
  <B><FONT COLOR="#228B22">void</FONT></B> assemble_ellipticdg(EquationSystems&amp; es, <B><FONT COLOR="#228B22">const</FONT></B> std::string&amp; system_name)
  {
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout&lt;&lt;<B><FONT COLOR="#BC8F8F">&quot; assembling elliptic dg system... &quot;</FONT></B>;
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout.flush();
  
    libmesh_assert_equal_to (system_name, <B><FONT COLOR="#BC8F8F">&quot;EllipticDG&quot;</FONT></B>);
    
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase&amp; mesh = es.get_mesh();
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim = mesh.mesh_dimension();
    
    LinearImplicitSystem &amp; ellipticdg_system = es.get_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;EllipticDG&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> Real penalty = es.parameters.get&lt;Real&gt; (<B><FONT COLOR="#BC8F8F">&quot;penalty&quot;</FONT></B>);
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::string refinement_type = es.parameters.get&lt;std::string&gt; (<B><FONT COLOR="#BC8F8F">&quot;refinement&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> DofMap &amp; dof_map = ellipticdg_system.get_dof_map();
    
    FEType fe_type = ellipticdg_system.variable_type(0);
  
    AutoPtr&lt;FEBase&gt; fe  (FEBase::build(dim, fe_type));
    AutoPtr&lt;FEBase&gt; fe_elem_face(FEBase::build(dim, fe_type));
    AutoPtr&lt;FEBase&gt; fe_neighbor_face(FEBase::build(dim, fe_type));
  
  #ifdef QORDER 
    QGauss qrule (dim, QORDER);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
    QGauss qrule (dim, fe_type.default_quadrature_order());
  #endif
    fe-&gt;attach_quadrature_rule (&amp;qrule);
    
  #ifdef QORDER 
    QGauss qface(dim-1, QORDER);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
    QGauss qface(dim-1, fe_type.default_quadrature_order());
  #endif
    
    fe_elem_face-&gt;attach_quadrature_rule(&amp;qface);
    fe_neighbor_face-&gt;attach_quadrature_rule(&amp;qface);
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW = fe-&gt;get_JxW();
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi = fe-&gt;get_dphi();
  
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp;  phi_face = fe_elem_face-&gt;get_phi();
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi_face = fe_elem_face-&gt;get_dphi();
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Real&gt;&amp; JxW_face = fe_elem_face-&gt;get_JxW();
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; qface_normals = fe_elem_face-&gt;get_normals();
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;Point&gt;&amp; qface_points = fe_elem_face-&gt;get_xyz();
      
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;Real&gt; &gt;&amp;  phi_neighbor_face = fe_neighbor_face-&gt;get_phi();
    <B><FONT COLOR="#228B22">const</FONT></B> std::vector&lt;std::vector&lt;RealGradient&gt; &gt;&amp; dphi_neighbor_face = fe_neighbor_face-&gt;get_dphi();
  
    DenseMatrix&lt;Number&gt; Ke;
    DenseVector&lt;Number&gt; Fe;
  
    DenseMatrix&lt;Number&gt; Kne;
    DenseMatrix&lt;Number&gt; Ken;
    DenseMatrix&lt;Number&gt; Kee;
    DenseMatrix&lt;Number&gt; Knn;
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; dof_indices;
  
    <B><FONT COLOR="#5F9EA0">MeshBase</FONT></B>::const_element_iterator el = mesh.active_local_elements_begin();
    <B><FONT COLOR="#228B22">const</FONT></B> MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
   
    <B><FONT COLOR="#A020F0">for</FONT></B> ( ; el != end_el; ++el)
    {
      <B><FONT COLOR="#228B22">const</FONT></B> Elem* elem = *el;
  
      dof_map.dof_indices (elem, dof_indices);
      <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_dofs   = dof_indices.size();
        
      fe-&gt;reinit  (elem);
  
      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);
   
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qrule.n_points(); qp++)
      {
        <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_dofs; i++)
        {
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_dofs; j++)
          {
             Ke(i,j) += JxW[qp]*(dphi[i][qp]*dphi[j][qp]);
          }
        }
      }
     
      <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> side=0; side&lt;elem-&gt;n_sides(); side++)
      {
        <B><FONT COLOR="#A020F0">if</FONT></B> (elem-&gt;neighbor(side) == NULL)
        {
          fe_elem_face-&gt;reinit(elem, side);
          
          AutoPtr&lt;Elem&gt; elem_side (elem-&gt;build_side(side));
          <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> elem_b_order = static_cast&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; (fe_elem_face-&gt;get_order());
          <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">double</FONT></B> h_elem = elem-&gt;volume()/elem_side-&gt;volume() * 1./pow(elem_b_order, 2.);
          
          <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qface.n_points(); qp++)
          {
            Number bc_value = exact_solution(qface_points[qp], es.parameters,<B><FONT COLOR="#BC8F8F">&quot;null&quot;</FONT></B>,<B><FONT COLOR="#BC8F8F">&quot;void&quot;</FONT></B>);
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_dofs; i++)
            {
               <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_dofs; j++)
               { 
                  Ke(i,j) += JxW_face[qp] * penalty/h_elem * phi_face[i][qp] * phi_face[j][qp]; <I><FONT COLOR="#B22222">// stability 
</FONT></I>                  Ke(i,j) -= JxW_face[qp] * (phi_face[i][qp] * (dphi_face[j][qp]*qface_normals[qp]) + phi_face[j][qp] * (dphi_face[i][qp]*qface_normals[qp])); <I><FONT COLOR="#B22222">// consistency
</FONT></I>               }
               Fe(i) += JxW_face[qp] * bc_value * penalty/h_elem * phi_face[i][qp];      <I><FONT COLOR="#B22222">// stability
</FONT></I>               Fe(i) -= JxW_face[qp] * dphi_face[i][qp] * (bc_value*qface_normals[qp]);     <I><FONT COLOR="#B22222">// consistency
</FONT></I>            }
          }
        }
        <B><FONT COLOR="#A020F0">else</FONT></B> 
        {
          <B><FONT COLOR="#228B22">const</FONT></B> Elem* neighbor = elem-&gt;neighbor(side);
           
          <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> elem_id = elem-&gt;id();
          <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> neighbor_id = neighbor-&gt;id();
          
          <B><FONT COLOR="#A020F0">if</FONT></B> ((neighbor-&gt;active() &amp;&amp; (neighbor-&gt;level() == elem-&gt;level()) &amp;&amp; (elem_id &lt; neighbor_id)) || (neighbor-&gt;level() &lt; elem-&gt;level()))
          {
            AutoPtr&lt;Elem&gt; elem_side (elem-&gt;build_side(side));
            
            <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> elem_b_order = static_cast&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;(fe_elem_face-&gt;get_order());
            <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> neighbor_b_order = static_cast&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;(fe_neighbor_face-&gt;get_order());
            <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">double</FONT></B> side_order = (elem_b_order + neighbor_b_order)/2.;
            <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">double</FONT></B> h_elem = (elem-&gt;volume()/elem_side-&gt;volume()) * 1./pow(side_order,2.);
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Point&gt; qface_neighbor_point;
  
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;Point &gt; qface_point;
           
            fe_elem_face-&gt;reinit(elem, side);
  
            qface_point = fe_elem_face-&gt;get_xyz();
  
  	  <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> side_neighbor = neighbor-&gt;which_neighbor_am_i(elem);
            <B><FONT COLOR="#A020F0">if</FONT></B> (refinement_type == <B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>)
              fe_neighbor_face-&gt;side_map (neighbor, elem_side.get(), side_neighbor, qface.get_points(), qface_neighbor_point);
  	  <B><FONT COLOR="#A020F0">else</FONT></B>
              <B><FONT COLOR="#5F9EA0">FEInterface</FONT></B>::inverse_map (elem-&gt;dim(), fe-&gt;get_fe_type(), neighbor, qface_point, qface_neighbor_point);
            fe_neighbor_face-&gt;reinit(neighbor, &amp;qface_neighbor_point);
            
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;dof_id_type&gt; neighbor_dof_indices;
            dof_map.dof_indices (neighbor, neighbor_dof_indices);
            <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> n_neighbor_dofs = neighbor_dof_indices.size();
  
            Kne.resize (n_neighbor_dofs, n_dofs);
            Ken.resize (n_dofs, n_neighbor_dofs);
            Kee.resize (n_dofs, n_dofs);
            Knn.resize (n_neighbor_dofs, n_neighbor_dofs);
  
            <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> qp=0; qp&lt;qface.n_points(); qp++)
            {
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_dofs; i++)          
              {
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_dofs; j++)
                {
                   Kee(i,j) -= 0.5 * JxW_face[qp] * (phi_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) + phi_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp])); <I><FONT COLOR="#B22222">// consistency
</FONT></I>                   Kee(i,j) += JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*phi_face[i][qp];  <I><FONT COLOR="#B22222">// stability
</FONT></I>                }
              }
  
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_neighbor_dofs; i++) 
              {
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_neighbor_dofs; j++)
                {
                   Knn(i,j) += 0.5 * JxW_face[qp] * (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp]) + phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp])); <I><FONT COLOR="#B22222">// consistency
</FONT></I>                   Knn(i,j) += JxW_face[qp] * penalty/h_elem * phi_neighbor_face[j][qp]*phi_neighbor_face[i][qp]; <I><FONT COLOR="#B22222">// stability
</FONT></I>                }
              }
  
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_neighbor_dofs; i++) 
              {
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_dofs; j++)
                {
                   Kne(i,j) += 0.5 * JxW_face[qp] * (phi_neighbor_face[i][qp]*(qface_normals[qp]*dphi_face[j][qp]) - phi_face[j][qp]*(qface_normals[qp]*dphi_neighbor_face[i][qp])); <I><FONT COLOR="#B22222">// consistency
</FONT></I>                   Kne(i,j) -= JxW_face[qp] * penalty/h_elem * phi_face[j][qp]*phi_neighbor_face[i][qp];  <I><FONT COLOR="#B22222">// stability
</FONT></I>                }
              }
              
              <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> i=0; i&lt;n_dofs; i++)         
              {
                <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> j=0; j&lt;n_neighbor_dofs; j++)
                {
                   Ken(i,j) += 0.5 * JxW_face[qp] * (phi_neighbor_face[j][qp]*(qface_normals[qp]*dphi_face[i][qp]) - phi_face[i][qp]*(qface_normals[qp]*dphi_neighbor_face[j][qp]));  <I><FONT COLOR="#B22222">// consistency
</FONT></I>                   Ken(i,j) -= JxW_face[qp] * penalty/h_elem * phi_face[i][qp]*phi_neighbor_face[j][qp];  <I><FONT COLOR="#B22222">// stability
</FONT></I>                }
              }
            }
        
            ellipticdg_system.matrix-&gt;add_matrix(Kne,neighbor_dof_indices,dof_indices);
            ellipticdg_system.matrix-&gt;add_matrix(Ken,dof_indices,neighbor_dof_indices);
            ellipticdg_system.matrix-&gt;add_matrix(Kee,dof_indices);
            ellipticdg_system.matrix-&gt;add_matrix(Knn,neighbor_dof_indices);
          }
        }
      }
      ellipticdg_system.matrix-&gt;add_matrix(Ke, dof_indices);
      ellipticdg_system.rhs-&gt;add_vector(Fe, dof_indices);
    }
  
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;done&quot;</FONT></B> &lt;&lt; std::endl;
  }
  
  
  
  <B><FONT COLOR="#228B22">int</FONT></B> main (<B><FONT COLOR="#228B22">int</FONT></B> argc, <B><FONT COLOR="#228B22">char</FONT></B>** argv)
  {
    LibMeshInit init(argc, argv);
  
    GetPot input_file(<B><FONT COLOR="#BC8F8F">&quot;miscellaneous_ex5.in&quot;</FONT></B>);
  
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> adaptive_refinement_steps = input_file(<B><FONT COLOR="#BC8F8F">&quot;max_adaptive_r_steps&quot;</FONT></B>,3);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> uniform_refinement_steps  = input_file(<B><FONT COLOR="#BC8F8F">&quot;uniform_h_r_steps&quot;</FONT></B>,3);
    <B><FONT COLOR="#228B22">const</FONT></B> Real refine_fraction                   = input_file(<B><FONT COLOR="#BC8F8F">&quot;refine_fraction&quot;</FONT></B>,0.5);
    <B><FONT COLOR="#228B22">const</FONT></B> Real coarsen_fraction                  = input_file(<B><FONT COLOR="#BC8F8F">&quot;coarsen_fraction&quot;</FONT></B>,0.);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> max_h_level               = input_file(<B><FONT COLOR="#BC8F8F">&quot;max_h_level&quot;</FONT></B>, 10);
    <B><FONT COLOR="#228B22">const</FONT></B> std::string refinement_type            = input_file(<B><FONT COLOR="#BC8F8F">&quot;refinement_type&quot;</FONT></B>,<B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>);
    Order p_order                                = static_cast&lt;Order&gt;(input_file(<B><FONT COLOR="#BC8F8F">&quot;p_order&quot;</FONT></B>,1));
    <B><FONT COLOR="#228B22">const</FONT></B> std::string element_type               = input_file(<B><FONT COLOR="#BC8F8F">&quot;element_type&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;tensor&quot;</FONT></B>);
    <B><FONT COLOR="#228B22">const</FONT></B> Real penalty                           = input_file(<B><FONT COLOR="#BC8F8F">&quot;ip_penalty&quot;</FONT></B>, 10.);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">bool</FONT></B> singularity                       = input_file(<B><FONT COLOR="#BC8F8F">&quot;singularity&quot;</FONT></B>, true);
    <B><FONT COLOR="#228B22">const</FONT></B> <B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> dim                       = input_file(<B><FONT COLOR="#BC8F8F">&quot;dimension&quot;</FONT></B>, 3);
  
    libmesh_example_assert(dim &lt;= LIBMESH_DIM, <B><FONT COLOR="#BC8F8F">&quot;2D/3D support&quot;</FONT></B>);
      
  #ifndef LIBMESH_ENABLE_AMR
    libmesh_example_assert(false, <B><FONT COLOR="#BC8F8F">&quot;--enable-amr&quot;</FONT></B>);
  #<B><FONT COLOR="#A020F0">else</FONT></B>
  
    Mesh mesh;
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 1)
      <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Generation::build_line(mesh,1,-1.,0.);
    <B><FONT COLOR="#A020F0">else</FONT></B> <B><FONT COLOR="#A020F0">if</FONT></B> (dim == 2)
      mesh.read(<B><FONT COLOR="#BC8F8F">&quot;lshaped.xda&quot;</FONT></B>);
    <B><FONT COLOR="#A020F0">else</FONT></B>
      mesh.read (<B><FONT COLOR="#BC8F8F">&quot;lshaped3D.xda&quot;</FONT></B>);
  
    <B><FONT COLOR="#A020F0">if</FONT></B> (element_type == <B><FONT COLOR="#BC8F8F">&quot;simplex&quot;</FONT></B>)
      <B><FONT COLOR="#5F9EA0">MeshTools</FONT></B>::Modification::all_tri(mesh);
  
    MeshRefinement mesh_refinement(mesh);
    mesh_refinement.refine_fraction() = refine_fraction;
    mesh_refinement.coarsen_fraction() = coarsen_fraction;
    mesh_refinement.max_h_level() = max_h_level;
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> rstep=0; rstep&lt;uniform_refinement_steps; rstep++)
    {
       mesh_refinement.uniformly_refine(1);
    }
   
    EquationSystems equation_system (mesh);
  
    equation_system.parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;linear solver tolerance&quot;</FONT></B>) = TOLERANCE * TOLERANCE;
    equation_system.parameters.set&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;linear solver maximum iterations&quot;</FONT></B>) = 1000;
    equation_system.parameters.set&lt;Real&gt;(<B><FONT COLOR="#BC8F8F">&quot;penalty&quot;</FONT></B>) = penalty;
    equation_system.parameters.set&lt;<B><FONT COLOR="#228B22">bool</FONT></B>&gt;(<B><FONT COLOR="#BC8F8F">&quot;singularity&quot;</FONT></B>) = singularity;
    equation_system.parameters.set&lt;std::string&gt;(<B><FONT COLOR="#BC8F8F">&quot;refinement&quot;</FONT></B>) = refinement_type;
  
    LinearImplicitSystem&amp; ellipticdg_system = equation_system.add_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;EllipticDG&quot;</FONT></B>);
   
    <B><FONT COLOR="#A020F0">if</FONT></B>( on_command_line( <B><FONT COLOR="#BC8F8F">&quot;element_type&quot;</FONT></B> ) )
      {
        <B><FONT COLOR="#5F9EA0">std</FONT></B>::string fe_str = command_line_value( std::string(<B><FONT COLOR="#BC8F8F">&quot;element_type&quot;</FONT></B>),
  					       <B><FONT COLOR="#5F9EA0">std</FONT></B>::string(<B><FONT COLOR="#BC8F8F">&quot;MONOMIAL&quot;</FONT></B>) );
        <B><FONT COLOR="#A020F0">if</FONT></B>( fe_str != <B><FONT COLOR="#BC8F8F">&quot;MONOMIAL&quot;</FONT></B> || fe_str != <B><FONT COLOR="#BC8F8F">&quot;XYZ&quot;</FONT></B> )
  	{
  	  <B><FONT COLOR="#5F9EA0">std</FONT></B>::cerr &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Error: This example must be run with MONOMIAL or XYZ element types.&quot;</FONT></B> &lt;&lt; std::endl;
  	  libmesh_error();
  	}
        ellipticdg_system.add_variable (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, p_order, Utility::string_to_enum&lt;FEFamily&gt;(fe_str) );
       } 
    <B><FONT COLOR="#A020F0">else</FONT></B>
      ellipticdg_system.add_variable (<B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>, p_order, MONOMIAL);
    
    ellipticdg_system.attach_assemble_function (assemble_ellipticdg);
  
    equation_system.init();
  
    ExactSolution exact_sol(equation_system);
    exact_sol.attach_exact_value(exact_solution);
    exact_sol.attach_exact_deriv(exact_derivative);
  
    <B><FONT COLOR="#A020F0">for</FONT></B> (<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B> rstep=0; rstep&lt;adaptive_refinement_steps; ++rstep)
    {
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;  Beginning Solve &quot;</FONT></B> &lt;&lt; rstep &lt;&lt; std::endl;
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Number of elements: &quot;</FONT></B> &lt;&lt; mesh.n_elem() &lt;&lt; std::endl;
         
      ellipticdg_system.solve();
              
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;System has: &quot;</FONT></B> &lt;&lt; equation_system.n_active_dofs()
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot; degrees of freedom.&quot;</FONT></B>
                &lt;&lt; std::endl;
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;Linear solver converged at step: &quot;</FONT></B>
                &lt;&lt; ellipticdg_system.n_linear_iterations()
                &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;, final residual: &quot;</FONT></B>
                &lt;&lt; ellipticdg_system.final_linear_residual()
                &lt;&lt; std::endl;
   
      exact_sol.compute_error(<B><FONT COLOR="#BC8F8F">&quot;EllipticDG&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>);
  
      <B><FONT COLOR="#5F9EA0">std</FONT></B>::cout &lt;&lt; <B><FONT COLOR="#BC8F8F">&quot;L2-Error is: &quot;</FONT></B>
                &lt;&lt; exact_sol.l2_error(<B><FONT COLOR="#BC8F8F">&quot;EllipticDG&quot;</FONT></B>, <B><FONT COLOR="#BC8F8F">&quot;u&quot;</FONT></B>)
                &lt;&lt; std::endl;
     
      <B><FONT COLOR="#A020F0">if</FONT></B> (rstep+1 &lt; adaptive_refinement_steps)
      {
        ErrorVector error;
  
        DiscontinuityMeasure error_estimator;
        error_estimator.estimate_error(ellipticdg_system,error);
        
        mesh_refinement.flag_elements_by_error_fraction(error);
        <B><FONT COLOR="#A020F0">if</FONT></B> (refinement_type == <B><FONT COLOR="#BC8F8F">&quot;p&quot;</FONT></B>)
          mesh_refinement.switch_h_to_p_refinement();
        <B><FONT COLOR="#A020F0">if</FONT></B> (refinement_type == <B><FONT COLOR="#BC8F8F">&quot;hp&quot;</FONT></B>)
          mesh_refinement.add_p_to_h_refinement();
        mesh_refinement.refine_and_coarsen_elements();
        equation_system.reinit();
      }
    }
  
  #ifdef LIBMESH_HAVE_EXODUS_API
    ExodusII_IO (mesh).write_discontinuous_exodusII(<B><FONT COLOR="#BC8F8F">&quot;lshaped_dg.e&quot;</FONT></B>,equation_system);
  #endif
  
  #endif <I><FONT COLOR="#B22222">// #ifndef LIBMESH_ENABLE_AMR
</FONT></I>    
    <B><FONT COLOR="#A020F0">return</FONT></B> 0;
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example miscellaneous_ex5:
*  mpirun -np 12 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
  Beginning Solve 0
Number of elements: 219
 assembling elliptic dg system... done
System has: 768 degrees of freedom.
Linear solver converged at step: 42, final residual: 8.61401e-12
L2-Error is: 0.00666744
  Beginning Solve 1
Number of elements: 827
 assembling elliptic dg system... done
System has: 2896 degrees of freedom.
Linear solver converged at step: 49, final residual: 9.80671e-12
L2-Error is: 0.00264921
  Beginning Solve 2
Number of elements: 3003
 assembling elliptic dg system... done
System has: 10512 degrees of freedom.
Linear solver converged at step: 68, final residual: 2.30969e-11
L2-Error is: 0.0016323
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

/workspace/libmesh/examples/miscellaneous/miscellaneous_ex5/.libs/lt-example-devel on a intel-12. named hbar.ices.utexas.edu with 12 processors, by benkirk Thu Jan 31 22:09:59 2013
Using Petsc Release Version 3.3.0, Patch 2, Fri Jul 13 15:42:00 CDT 2012 

                         Max       Max/Min        Avg      Total 
Time (sec):           3.569e+00      1.00040   3.568e+00
Objects:              2.020e+02      1.00000   2.020e+02
Flops:                7.765e+07      1.86013   5.768e+07  6.922e+08
Flops/sec:            2.176e+07      1.85941   1.616e+07  1.940e+08
MPI Messages:         1.386e+03      1.82075   1.017e+03  1.220e+04
MPI Message Lengths:  6.411e+05      1.97396   4.708e+02  5.745e+06
MPI Reductions:       5.130e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 3.5684e+00 100.0%  6.9216e+08 100.0%  1.220e+04 100.0%  4.708e+02      100.0%  5.120e+02  99.8% 

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

VecMDot              159 1.0 1.8934e-02 4.6 2.12e+06 1.0 0.0e+00 0.0e+00 1.6e+02  0  4  0  0 31   0  4  0  0 31  1316
VecNorm              169 1.0 2.7646e-02 4.9 1.61e+05 1.0 0.0e+00 0.0e+00 1.7e+02  0  0  0  0 33   0  0  0  0 33    68
VecScale             166 1.0 2.0695e-04 1.3 7.91e+04 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  4483
VecCopy               16 1.0 3.6001e-05 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               199 1.0 2.7466e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               14 1.0 1.3553e-02 1.9 1.32e+04 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0    11
VecMAXPY             166 1.0 1.3723e-03 1.1 2.28e+06 1.0 0.0e+00 0.0e+00 0.0e+00  0  4  0  0  0   0  4  0  0  0 19467
VecAssemblyBegin      23 1.0 3.5905e-0247.0 0.00e+00 0.0 0.0e+00 0.0e+00 5.7e+01  1  0  0  0 11   1  0  0  0 11     0
VecAssemblyEnd        23 1.0 5.0545e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      175 1.0 1.6456e-03 1.4 0.00e+00 0.0 1.1e+04 4.5e+02 0.0e+00  0  0 89 85  0   0  0 89 85  0     0
VecScatterEnd        175 1.0 4.4419e-0215.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0     0
VecNormalize         166 1.0 2.7815e-02 4.8 2.37e+05 1.0 0.0e+00 0.0e+00 1.7e+02  0  0  0  0 32   0  0  0  0 32   100
MatMult              166 1.0 5.2725e-02 4.7 4.36e+06 1.1 1.0e+04 4.5e+02 0.0e+00  1  7 85 80  0   1  7 85 80  0   945
MatSolve             169 1.0 4.4376e-02 1.5 3.27e+07 1.6 0.0e+00 0.0e+00 0.0e+00  1 46  0  0  0   1 46  0  0  0  7158
MatLUFactorNum         3 1.0 2.5683e-02 2.8 3.62e+07 2.9 0.0e+00 0.0e+00 0.0e+00  0 39  0  0  0   0 39  0  0  0 10517
MatILUFactorSym        3 1.0 4.8430e-02 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 9.0e+00  1  0  0  0  2   1  0  0  0  2     0
MatAssemblyBegin       6 1.0 8.3387e-02 3.7 0.00e+00 0.0 3.9e+02 1.9e+03 1.2e+01  2  0  3 13  2   2  0  3 13  2     0
MatAssemblyEnd         6 1.0 1.7877e-03 1.2 0.00e+00 0.0 3.6e+02 1.1e+02 2.4e+01  0  0  3  1  5   0  0  3  1  5     0
MatGetRowIJ            3 1.0 1.5020e-05 3.7 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         3 1.0 2.9111e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 1.2e+01  0  0  0  0  2   0  0  0  0  2     0
MatZeroEntries         9 1.0 1.6832e-04 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       159 1.0 2.0551e-02 3.7 4.25e+06 1.0 0.0e+00 0.0e+00 1.6e+02  0  7  0  0 31   0  7  0  0 31  2426
KSPSetUp               6 1.0 1.7929e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               3 1.0 1.6063e-01 1.0 7.76e+07 1.9 1.0e+04 4.5e+02 3.6e+02  4100 85 80 69   4100 85 80 69  4309
PCSetUp                6 1.0 7.5485e-02 2.2 3.62e+07 2.9 0.0e+00 0.0e+00 2.7e+01  1 39  0  0  5   1 39  0  0  5  3578
PCSetUpOnBlocks        3 1.0 7.4771e-02 2.3 3.62e+07 2.9 0.0e+00 0.0e+00 2.1e+01  1 39  0  0  4   1 39  0  0  4  3612
PCApply              169 1.0 4.6537e-02 1.4 3.27e+07 1.6 0.0e+00 0.0e+00 0.0e+00  1 46  0  0  0   1 46  0  0  0  6825
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

              Vector   137            137       722872     0
      Vector Scatter     8              8         8288     0
           Index Set    29             29        38848     0
   IS L to G Mapping     3              3         1692     0
              Matrix    12             12      2599580     0
       Krylov Solver     6              6        58080     0
      Preconditioner     6              6         5352     0
              Viewer     1              0            0     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 4.57764e-06
Average time for zero size MPI_Send(): 1.508e-05
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
| Time:           Thu Jan 31 22:09:59 2013                                                                             |
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
 ---------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=3.70542, Active time=3.44529                                                        |
 ---------------------------------------------------------------------------------------------------------------------
| Event                                   nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                   w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|---------------------------------------------------------------------------------------------------------------------|
|                                                                                                                     |
|                                                                                                                     |
| DofMap                                                                                                              |
|   add_neighbors_to_send_list()          3         0.0189      0.006286    0.0453      0.015112    0.55     1.32     |
|   build_sparsity()                      3         0.0913      0.030421    0.5473      0.182443    2.65     15.89    |
|   create_dof_constraints()              3         0.1001      0.033362    0.1617      0.053915    2.90     4.69     |
|   distribute_dofs()                     3         0.0710      0.023663    0.2168      0.072268    2.06     6.29     |
|   dof_indices()                         12370     0.7245      0.000059    0.7245      0.000059    21.03    21.03    |
|   old_dof_indices()                     556       0.0314      0.000057    0.0314      0.000057    0.91     0.91     |
|   prepare_send_list()                   3         0.0007      0.000232    0.0007      0.000232    0.02     0.02     |
|   reinit()                              3         0.1358      0.045277    0.1358      0.045277    3.94     3.94     |
|                                                                                                                     |
| EquationSystems                                                                                                     |
|   build_discontinuous_solution_vector() 1         0.0390      0.039014    0.1966      0.196598    1.13     5.71     |
|                                                                                                                     |
| ExodusII_IO                                                                                                         |
|   write_nodal_data_discontinuous()      1         0.0183      0.018300    0.0183      0.018300    0.53     0.53     |
|                                                                                                                     |
| FE                                                                                                                  |
|   compute_shape_functions()             3425      0.0679      0.000020    0.0679      0.000020    1.97     1.97     |
|   init_shape_functions()                2843      0.0321      0.000011    0.0321      0.000011    0.93     0.93     |
|   inverse_map()                         12324     0.2067      0.000017    0.2067      0.000017    6.00     6.00     |
|                                                                                                                     |
| FEMap                                                                                                               |
|   compute_affine_map()                  3425      0.0582      0.000017    0.0582      0.000017    1.69     1.69     |
|   compute_face_map()                    1228      0.0232      0.000019    0.0232      0.000019    0.67     0.67     |
|   init_face_shape_functions()           5         0.0002      0.000042    0.0002      0.000042    0.01     0.01     |
|   init_reference_to_physical_map()      2843      0.1614      0.000057    0.1614      0.000057    4.68     4.68     |
|                                                                                                                     |
| JumpErrorEstimator                                                                                                  |
|   estimate_error()                      2         0.0246      0.012309    0.1358      0.067923    0.71     3.94     |
|                                                                                                                     |
| LocationMap                                                                                                         |
|   find()                                22344     0.0939      0.000004    0.0939      0.000004    2.73     2.73     |
|   init()                                6         0.0111      0.001850    0.0111      0.001850    0.32     0.32     |
|                                                                                                                     |
| Mesh                                                                                                                |
|   contract()                            2         0.0022      0.001076    0.0057      0.002825    0.06     0.16     |
|   find_neighbors()                      5         0.1651      0.033026    0.1694      0.033885    4.79     4.92     |
|   renumber_nodes_and_elem()             12        0.0098      0.000814    0.0098      0.000814    0.28     0.28     |
|                                                                                                                     |
| MeshCommunication                                                                                                   |
|   compute_hilbert_indices()             5         0.0303      0.006061    0.0303      0.006061    0.88     0.88     |
|   find_global_indices()                 5         0.0130      0.002600    0.0525      0.010498    0.38     1.52     |
|   parallel_sort()                       5         0.0049      0.000973    0.0058      0.001155    0.14     0.17     |
|                                                                                                                     |
| MeshRefinement                                                                                                      |
|   _coarsen_elements()                   4         0.0023      0.000573    0.0024      0.000590    0.07     0.07     |
|   _refine_elements()                    6         0.1390      0.023161    0.3691      0.061512    4.03     10.71    |
|   add_point()                           22344     0.1129      0.000005    0.2131      0.000010    3.28     6.19     |
|   make_coarsening_compatible()          8         0.0423      0.005288    0.0423      0.005288    1.23     1.23     |
|   make_refinement_compatible()          8         0.0023      0.000282    0.0024      0.000302    0.07     0.07     |
|                                                                                                                     |
| MetisPartitioner                                                                                                    |
|   partition()                           5         0.3959      0.079184    0.4499      0.089974    11.49    13.06    |
|                                                                                                                     |
| Parallel                                                                                                            |
|   allgather()                           19        0.0031      0.000163    0.0032      0.000170    0.09     0.09     |
|   broadcast()                           20        0.0005      0.000024    0.0004      0.000019    0.01     0.01     |
|   max(bool)                             23        0.0114      0.000496    0.0114      0.000496    0.33     0.33     |
|   max(scalar)                           494       0.0160      0.000032    0.0160      0.000032    0.47     0.47     |
|   max(vector)                           119       0.0022      0.000019    0.0061      0.000051    0.06     0.18     |
|   min(bool)                             611       0.0089      0.000015    0.0089      0.000015    0.26     0.26     |
|   min(scalar)                           483       0.0927      0.000192    0.0927      0.000192    2.69     2.69     |
|   min(vector)                           119       0.0024      0.000020    0.0065      0.000054    0.07     0.19     |
|   probe()                               440       0.0043      0.000010    0.0043      0.000010    0.12     0.12     |
|   receive()                             440       0.0034      0.000008    0.0078      0.000018    0.10     0.23     |
|   send()                                440       0.0016      0.000004    0.0016      0.000004    0.05     0.05     |
|   send_receive()                        450       0.0042      0.000009    0.0152      0.000034    0.12     0.44     |
|   sum()                                 25        0.0008      0.000033    0.0018      0.000072    0.02     0.05     |
|                                                                                                                     |
| Parallel::Request                                                                                                   |
|   wait()                                440       0.0012      0.000003    0.0012      0.000003    0.04     0.04     |
|                                                                                                                     |
| Partitioner                                                                                                         |
|   set_node_processor_ids()              6         0.0134      0.002233    0.0212      0.003532    0.39     0.62     |
|   set_parent_processor_ids()            5         0.0122      0.002444    0.0122      0.002444    0.35     0.35     |
|                                                                                                                     |
| PetscLinearSolver                                                                                                   |
|   solve()                               3         0.2407      0.080244    0.2407      0.080244    6.99     6.99     |
|                                                                                                                     |
| ProjectVector                                                                                                       |
|   operator()                            2         0.0203      0.010127    0.1629      0.081443    0.59     4.73     |
|                                                                                                                     |
| System                                                                                                              |
|   assemble()                            3         0.1363      0.045446    0.4865      0.162183    3.96     14.12    |
|   project_vector()                      2         0.0385      0.019232    0.2179      0.108952    1.12     6.32     |
|                                                                                                                     |
| XdrIO                                                                                                               |
|   read()                                1         0.0011      0.001083    0.0012      0.001183    0.03     0.03     |
 ---------------------------------------------------------------------------------------------------------------------
| Totals:                                 87945     3.4453                                          100.00            |
 ---------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example miscellaneous_ex5:
*  mpirun -np 12 example-devel  -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
