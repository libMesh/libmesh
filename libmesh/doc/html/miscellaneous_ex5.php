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
        #include "libmesh.h"
        #include "mesh.h"
        #include "equation_systems.h"
        #include "mesh_data.h"
        #include "mesh_generation.h"
        #include "mesh_modification.h"
        #include "elem.h"
        #include "transient_system.h"
        #include "fe.h"
        #include "quadrature_gauss.h"
        #include "dof_map.h"
        #include "sparse_matrix.h"
        #include "dense_matrix.h"
        #include "dense_vector.h"
        #include "dense_submatrix.h"
        #include "dense_subvector.h"
        #include "numeric_vector.h"
        #include "linear_implicit_system.h"
        #include "exodusII_io.h"
        #include "fe_interface.h"
        #include "getpot.h"
        #include "mesh_refinement.h"
        #include "error_vector.h"
        #include "kelly_error_estimator.h"
        #include "discontinuity_measure.h"
        #include "string_to_enum.h"
        
        #include "exact_solution.h"
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
              libmesh_assert (x != 0.);
        
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
          libmesh_assert (system_name == "EllipticDG");
          
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
          std::vector&lt;unsigned int&gt; dof_indices;
        
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
                  std::vector&lt;unsigned int&gt; neighbor_dof_indices;
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
<br><br><br> <h1> The program without comments: </h1> 
<pre> 
  
  #include &lt;iostream&gt;
  
  #include <B><FONT COLOR="#BC8F8F">&quot;libmesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;equation_systems.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_data.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_generation.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_modification.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;elem.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;transient_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;quadrature_gauss.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dof_map.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;sparse_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_matrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_submatrix.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;dense_subvector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;numeric_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;linear_implicit_system.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;exodusII_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe_interface.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;kelly_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;discontinuity_measure.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;string_to_enum.h&quot;</FONT></B>
  
  #include <B><FONT COLOR="#BC8F8F">&quot;exact_solution.h&quot;</FONT></B>
  
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
        libmesh_assert (x != 0.);
  
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
  
    libmesh_assert (system_name == <B><FONT COLOR="#BC8F8F">&quot;EllipticDG&quot;</FONT></B>);
    
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
    
    <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; dof_indices;
  
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
            
            <B><FONT COLOR="#5F9EA0">std</FONT></B>::vector&lt;<B><FONT COLOR="#228B22">unsigned</FONT></B> <B><FONT COLOR="#228B22">int</FONT></B>&gt; neighbor_dof_indices;
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
Linking ./miscellaneous_ex5-opt...
***************************************************************
* Running Example  mpirun -np 6 ./miscellaneous_ex5-opt --element_type MONOMIAL -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
  Beginning Solve 0
Number of elements: 219
 assembling elliptic dg system... done
System has: 768 degrees of freedom.
Linear solver converged at step: 35, final residual: 1.7652e-11
L2-Error is: 0.00666744
  Beginning Solve 1
Number of elements: 827
 assembling elliptic dg system... done
System has: 2896 degrees of freedom.
Linear solver converged at step: 47, final residual: 2.41141e-11
L2-Error is: 0.00264921
  Beginning Solve 2
Number of elements: 3003
 assembling elliptic dg system... done
System has: 10512 degrees of freedom.
Linear solver converged at step: 60, final residual: 1.65543e-11
L2-Error is: 0.0016323
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./miscellaneous_ex5-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:21:54 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.152e+00      1.02372   1.139e+00
Objects:              1.980e+02      1.00000   1.980e+02
Flops:                2.239e+08      1.80081   1.720e+08  1.032e+09
Flops/sec:            1.943e+08      1.78436   1.510e+08  9.057e+08
MPI Messages:         7.125e+02      1.51596   6.490e+02  3.894e+03
MPI Message Lengths:  7.516e+05      1.92707   9.139e+02  3.559e+06
MPI Reductions:       4.830e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.1390e+00 100.0%  1.0321e+09 100.0%  3.894e+03 100.0%  9.139e+02      100.0%  4.290e+02  88.8% 

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

VecMDot              142 1.0 7.0386e-0223.5 3.98e+06 1.0 0.0e+00 0.0e+00 1.4e+02  4  2  0  0 29   4  2  0  0 33   339
VecNorm              151 1.0 1.1680e-01 3.5 2.79e+05 1.0 0.0e+00 0.0e+00 1.5e+02  5  0  0  0 31   5  0  0  0 35    14
VecScale             148 1.0 3.0613e-04 2.4 1.37e+05 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2685
VecCopy               27 1.0 5.4836e-05 2.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               170 1.0 2.4104e-04 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               12 1.0 9.8445e-03 1.6 1.89e+04 1.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0    12
VecMAXPY             148 1.0 1.8954e-03 1.2 4.25e+06 1.0 0.0e+00 0.0e+00 0.0e+00  0  2  0  0  0   0  2  0  0  0 13430
VecAssemblyBegin      23 1.0 3.8818e-02 3.5 0.00e+00 0.0 0.0e+00 0.0e+00 5.7e+01  3  0  0  0 12   3  0  0  0 13     0
VecAssemblyEnd        23 1.0 3.8147e-05 1.2 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      157 1.0 1.4141e-03 1.7 0.00e+00 0.0 3.4e+03 8.7e+02 0.0e+00  0  0 88 84  0   0  0 88 84  0     0
VecScatterEnd        157 1.0 2.4418e-01314.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 13  0  0  0  0  13  0  0  0  0     0
VecNormalize         148 1.0 1.1566e-01 6.5 4.12e+05 1.0 0.0e+00 0.0e+00 1.5e+02  4  0  0  0 31   4  0  0  0 34    21
MatMult              148 1.0 2.5525e-0128.7 7.45e+06 1.0 3.2e+03 8.7e+02 0.0e+00 14  4 82 78  0  14  4 82 78  0   173
MatSolve             148 1.0 9.3645e-02 2.1 7.35e+07 1.4 0.0e+00 0.0e+00 0.0e+00  6 37  0  0  0   6 37  0  0  0  4034
MatLUFactorNum         3 1.0 1.1166e-01 3.1 1.34e+08 2.5 0.0e+00 0.0e+00 0.0e+00  6 54  0  0  0   6 54  0  0  0  4999
MatILUFactorSym        3 1.0 1.6950e-01 2.4 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00 10  0  0  0  1  10  0  0  0  1     0
MatAssemblyBegin       6 1.0 1.0412e-01 9.1 0.00e+00 0.0 1.4e+02 3.6e+03 1.2e+01  5  0  4 14  2   5  0  4 14  3     0
MatAssemblyEnd         6 1.0 3.8126e-03 1.2 0.00e+00 0.0 1.3e+02 2.0e+02 2.4e+01  0  0  3  1  5   0  0  3  1  6     0
MatGetRowIJ            3 1.0 1.6928e-05 8.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         3 1.0 1.2064e-04 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 1.2e+01  0  0  0  0  2   0  0  0  0  3     0
MatZeroEntries         9 1.0 2.5988e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       142 1.0 7.2324e-0215.5 7.97e+06 1.0 0.0e+00 0.0e+00 1.4e+02  4  5  0  0 29   4  5  0  0 33   660
KSPSetup               6 1.0 1.2994e-04 1.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               3 1.0 5.2132e-01 1.0 2.24e+08 1.8 3.2e+03 8.7e+02 3.1e+02 46100 82 78 64  46100 82 78 72  1980
PCSetUp                6 1.0 2.8175e-01 2.6 1.34e+08 2.5 0.0e+00 0.0e+00 1.5e+01 16 54  0  0  3  16 54  0  0  3  1981
PCSetUpOnBlocks        3 1.0 2.8139e-01 2.6 1.34e+08 2.5 0.0e+00 0.0e+00 1.5e+01 16 54  0  0  3  16 54  0  0  3  1984
PCApply              148 1.0 9.5374e-02 2.1 7.35e+07 1.4 0.0e+00 0.0e+00 0.0e+00  6 37  0  0  0   6 37  0  0  0  3961
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec   134            134      1050304     0
         Vec Scatter     8              8         6944     0
           Index Set    29             29        54216     0
   IS L to G Mapping     3              3        14044     0
              Matrix    12             12      8216292     0
       Krylov Solver     6              6        56640     0
      Preconditioner     6              6         4224     0
========================================================================================================================
Average time to get PetscTime(): 0
Average time for MPI_Barrier(): 4.32014e-05
Average time for zero size MPI_Send(): 4.48227e-05
#PETSc Option Table entries:
--element_type MONOMIAL
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
| Time:           Fri Aug 24 15:21:54 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ---------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=1.34221, Active time=1.0546                                                         |
 ---------------------------------------------------------------------------------------------------------------------
| Event                                   nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                   w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|---------------------------------------------------------------------------------------------------------------------|
|                                                                                                                     |
|                                                                                                                     |
| DofMap                                                                                                              |
|   add_neighbors_to_send_list()          3         0.0007      0.000227    0.0009      0.000286    0.06     0.08     |
|   build_sparsity()                      3         0.0177      0.005891    0.0239      0.007981    1.68     2.27     |
|   create_dof_constraints()              3         0.0082      0.002719    0.0130      0.004342    0.77     1.24     |
|   distribute_dofs()                     3         0.0014      0.000466    0.0107      0.003555    0.13     1.01     |
|   dof_indices()                         24234     0.0087      0.000000    0.0087      0.000000    0.83     0.83     |
|   old_dof_indices()                     1214      0.0004      0.000000    0.0004      0.000000    0.04     0.04     |
|   prepare_send_list()                   3         0.0000      0.000016    0.0000      0.000016    0.00     0.00     |
|   reinit()                              3         0.0029      0.000971    0.0029      0.000971    0.28     0.28     |
|                                                                                                                     |
| EquationSystems                                                                                                     |
|   build_discontinuous_solution_vector() 1         0.0211      0.021083    0.0224      0.022445    2.00     2.13     |
|                                                                                                                     |
| ExodusII_IO                                                                                                         |
|   write_nodal_data_discontinuous()      1         0.0047      0.004715    0.0047      0.004715    0.45     0.45     |
|                                                                                                                     |
| FE                                                                                                                  |
|   compute_shape_functions()             6758      0.0076      0.000001    0.0076      0.000001    0.72     0.72     |
|   init_shape_functions()                5582      0.0056      0.000001    0.0056      0.000001    0.53     0.53     |
|   inverse_map()                         18340     0.0330      0.000002    0.0330      0.000002    3.13     3.13     |
|                                                                                                                     |
| FEMap                                                                                                               |
|   compute_affine_map()                  6758      0.0076      0.000001    0.0076      0.000001    0.72     0.72     |
|   compute_face_map()                    2463      0.0032      0.000001    0.0032      0.000001    0.31     0.31     |
|   init_face_shape_functions()           5         0.0000      0.000006    0.0000      0.000006    0.00     0.00     |
|   init_reference_to_physical_map()      5582      0.0307      0.000005    0.0307      0.000005    2.91     2.91     |
|                                                                                                                     |
| JumpErrorEstimator                                                                                                  |
|   estimate_error()                      2         0.0054      0.002691    0.0336      0.016801    0.51     3.19     |
|                                                                                                                     |
| LocationMap                                                                                                         |
|   find()                                22344     0.0057      0.000000    0.0057      0.000000    0.54     0.54     |
|   init()                                6         0.0008      0.000138    0.0008      0.000138    0.08     0.08     |
|                                                                                                                     |
| Mesh                                                                                                                |
|   contract()                            2         0.0002      0.000088    0.0005      0.000266    0.02     0.05     |
|   find_neighbors()                      5         0.0085      0.001703    0.0158      0.003161    0.81     1.50     |
|   renumber_nodes_and_elem()             12        0.0013      0.000105    0.0013      0.000105    0.12     0.12     |
|                                                                                                                     |
| MeshCommunication                                                                                                   |
|   compute_hilbert_indices()             5         0.0055      0.001094    0.0055      0.001094    0.52     0.52     |
|   find_global_indices()                 5         0.0007      0.000133    0.0172      0.003438    0.06     1.63     |
|   parallel_sort()                       5         0.0044      0.000881    0.0096      0.001921    0.42     0.91     |
|                                                                                                                     |
| MeshRefinement                                                                                                      |
|   _coarsen_elements()                   4         0.0002      0.000058    0.0007      0.000169    0.02     0.06     |
|   _refine_elements()                    6         0.0184      0.003073    0.0900      0.015002    1.75     8.54     |
|   add_point()                           22344     0.0155      0.000001    0.0228      0.000001    1.47     2.16     |
|   make_coarsening_compatible()          4         0.0018      0.000442    0.0018      0.000442    0.17     0.17     |
|   make_refinement_compatible()          4         0.0001      0.000020    0.0077      0.001935    0.01     0.73     |
|                                                                                                                     |
| MetisPartitioner                                                                                                    |
|   partition()                           5         0.0088      0.001759    0.0341      0.006826    0.83     3.24     |
|                                                                                                                     |
| Parallel                                                                                                            |
|   allgather()                           19        0.0054      0.000285    0.0054      0.000285    0.51     0.51     |
|   broadcast()                           20        0.0001      0.000004    0.0001      0.000003    0.01     0.01     |
|   max(bool)                             16        0.0520      0.003250    0.0520      0.003250    4.93     4.93     |
|   max(scalar)                           13        0.0340      0.002614    0.0340      0.002614    3.22     3.22     |
|   max(vector)                           5         0.0002      0.000038    0.0002      0.000038    0.02     0.02     |
|   min(bool)                             8         0.0097      0.001210    0.0097      0.001210    0.92     0.92     |
|   min(scalar)                           2         0.0001      0.000051    0.0001      0.000051    0.01     0.01     |
|   min(vector)                           5         0.0080      0.001590    0.0080      0.001590    0.75     0.75     |
|   probe()                               170       0.0171      0.000101    0.0171      0.000101    1.62     1.62     |
|   receive()                             170       0.0003      0.000002    0.0174      0.000102    0.02     1.65     |
|   send()                                170       0.0001      0.000001    0.0001      0.000001    0.01     0.01     |
|   send_receive()                        180       0.0003      0.000002    0.0180      0.000100    0.03     1.70     |
|   sum()                                 25        0.0228      0.000913    0.0228      0.000913    2.16     2.16     |
|                                                                                                                     |
| Parallel::Request                                                                                                   |
|   wait()                                170       0.0001      0.000000    0.0001      0.000000    0.01     0.01     |
|                                                                                                                     |
| Partitioner                                                                                                         |
|   set_node_processor_ids()              6         0.0010      0.000159    0.0166      0.002759    0.09     1.57     |
|   set_parent_processor_ids()            5         0.0006      0.000123    0.0006      0.000123    0.06     0.06     |
|                                                                                                                     |
| PetscLinearSolver                                                                                                   |
|   solve()                               3         0.5971      0.199046    0.5971      0.199046    56.62    56.62    |
|                                                                                                                     |
| ProjectVector                                                                                                       |
|   operator()                            2         0.0066      0.003298    0.0277      0.013872    0.63     2.63     |
|                                                                                                                     |
| System                                                                                                              |
|   assemble()                            3         0.0334      0.011133    0.0864      0.028804    3.17     8.19     |
|   project_vector()                      2         0.0349      0.017434    0.0629      0.031474    3.31     5.97     |
|                                                                                                                     |
| XdrIO                                                                                                               |
|   read()                                1         0.0001      0.000138    0.0002      0.000157    0.01     0.01     |
 ---------------------------------------------------------------------------------------------------------------------
| Totals:                                 116704    1.0546                                          100.00            |
 ---------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 6 ./miscellaneous_ex5-opt --element_type MONOMIAL -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
***************************************************************
* Running Example  mpirun -np 6 ./miscellaneous_ex5-opt --element_type XYZ -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
  Beginning Solve 0
Number of elements: 219
 assembling elliptic dg system... done
System has: 768 degrees of freedom.
Linear solver converged at step: 35, final residual: 1.76521e-11
L2-Error is: 0.00666744
  Beginning Solve 1
Number of elements: 827
 assembling elliptic dg system... done
System has: 2896 degrees of freedom.
Linear solver converged at step: 47, final residual: 2.41141e-11
L2-Error is: 0.00264921
  Beginning Solve 2
Number of elements: 3003
 assembling elliptic dg system... done
System has: 10512 degrees of freedom.
Linear solver converged at step: 60, final residual: 1.65543e-11
L2-Error is: 0.0016323
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./miscellaneous_ex5-opt on a intel-11. named daedalus with 6 processors, by roystgnr Fri Aug 24 15:21:56 2012
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.308e+00      1.07569   1.232e+00
Objects:              1.980e+02      1.00000   1.980e+02
Flops:                2.239e+08      1.80081   1.720e+08  1.032e+09
Flops/sec:            1.712e+08      1.67703   1.393e+08  8.356e+08
MPI Messages:         7.125e+02      1.51596   6.490e+02  3.894e+03
MPI Message Lengths:  7.516e+05      1.92707   9.139e+02  3.559e+06
MPI Reductions:       4.830e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.2317e+00 100.0%  1.0321e+09 100.0%  3.894e+03 100.0%  9.139e+02      100.0%  4.290e+02  88.8% 

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

VecMDot              142 1.0 4.4302e-02 3.1 3.98e+06 1.0 0.0e+00 0.0e+00 1.4e+02  2  2  0  0 29   2  2  0  0 33   539
VecNorm              151 1.0 1.9510e-0118.7 2.79e+05 1.0 0.0e+00 0.0e+00 1.5e+02  9  0  0  0 31   9  0  0  0 35     9
VecScale             148 1.0 4.0388e-04 2.9 1.37e+05 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  2035
VecCopy               27 1.0 4.7684e-05 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               170 1.0 2.9302e-04 1.5 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               12 1.0 6.7282e-03 1.1 1.89e+04 1.0 0.0e+00 0.0e+00 0.0e+00  1  0  0  0  0   1  0  0  0  0    17
VecMAXPY             148 1.0 2.3487e-03 1.4 4.25e+06 1.0 0.0e+00 0.0e+00 0.0e+00  0  2  0  0  0   0  2  0  0  0 10839
VecAssemblyBegin      23 1.0 5.4118e-02 4.6 0.00e+00 0.0 0.0e+00 0.0e+00 5.7e+01  3  0  0  0 12   3  0  0  0 13     0
VecAssemblyEnd        23 1.0 4.5538e-05 1.6 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      157 1.0 1.5037e-03 2.0 0.00e+00 0.0 3.4e+03 8.7e+02 0.0e+00  0  0 88 84  0   0  0 88 84  0     0
VecScatterEnd        157 1.0 2.8792e-01 7.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00 11  0  0  0  0  11  0  0  0  0     0
VecNormalize         148 1.0 1.9541e-0118.5 4.12e+05 1.0 0.0e+00 0.0e+00 1.5e+02  9  0  0  0 31   9  0  0  0 34    13
MatMult              148 1.0 2.9654e-01 8.0 7.45e+06 1.0 3.2e+03 8.7e+02 0.0e+00 12  4 82 78  0  12  4 82 78  0   149
MatSolve             148 1.0 8.1570e-02 1.6 7.35e+07 1.4 0.0e+00 0.0e+00 0.0e+00  5 37  0  0  0   5 37  0  0  0  4632
MatLUFactorNum         3 1.0 1.0142e-01 2.1 1.34e+08 2.5 0.0e+00 0.0e+00 0.0e+00  5 54  0  0  0   5 54  0  0  0  5504
MatILUFactorSym        3 1.0 2.0841e-01 2.6 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00  9  0  0  0  1   9  0  0  0  1     0
MatAssemblyBegin       6 1.0 1.5809e-0113.4 0.00e+00 0.0 1.4e+02 3.6e+03 1.2e+01  9  0  4 14  2   9  0  4 14  3     0
MatAssemblyEnd         6 1.0 4.6549e-03 1.2 0.00e+00 0.0 1.3e+02 2.0e+02 2.4e+01  0  0  3  1  5   0  0  3  1  6     0
MatGetRowIJ            3 1.0 8.6188e-04451.9 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         3 1.0 9.6798e-0410.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.2e+01  0  0  0  0  2   0  0  0  0  3     0
MatZeroEntries         9 1.0 3.6287e-04 2.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       142 1.0 4.6014e-02 2.8 7.97e+06 1.0 0.0e+00 0.0e+00 1.4e+02  2  5  0  0 29   2  5  0  0 33  1037
KSPSetup               6 1.0 1.1277e-04 1.3 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               3 1.0 5.4144e-01 1.0 2.24e+08 1.8 3.2e+03 8.7e+02 3.1e+02 44100 82 78 64  44100 82 78 72  1906
PCSetUp                6 1.0 3.1047e-01 2.4 1.34e+08 2.5 0.0e+00 0.0e+00 1.5e+01 15 54  0  0  3  15 54  0  0  3  1798
PCSetUpOnBlocks        3 1.0 3.1007e-01 2.4 1.34e+08 2.5 0.0e+00 0.0e+00 1.5e+01 15 54  0  0  3  15 54  0  0  3  1800
PCApply              148 1.0 8.3165e-02 1.5 7.35e+07 1.4 0.0e+00 0.0e+00 0.0e+00  6 37  0  0  0   6 37  0  0  0  4543
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec   134            134      1050304     0
         Vec Scatter     8              8         6944     0
           Index Set    29             29        54216     0
   IS L to G Mapping     3              3        14044     0
              Matrix    12             12      8216292     0
       Krylov Solver     6              6        56640     0
      Preconditioner     6              6         4224     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 9.01222e-06
Average time for zero size MPI_Send(): 1.13249e-05
#PETSc Option Table entries:
--element_type XYZ
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
| Time:           Fri Aug 24 15:21:56 2012                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-34-generic                                |
| OS Version:     #76-Ubuntu SMP Tue Aug 30 17:05:01 UTC 2011      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Wed Aug 22 12:44:06 CDT 2012  |
-------------------------------------------------------------------
 ---------------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=1.4253, Active time=1.14817                                                         |
 ---------------------------------------------------------------------------------------------------------------------
| Event                                   nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                                   w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|---------------------------------------------------------------------------------------------------------------------|
|                                                                                                                     |
|                                                                                                                     |
| DofMap                                                                                                              |
|   add_neighbors_to_send_list()          3         0.0008      0.000254    0.0009      0.000307    0.07     0.08     |
|   build_sparsity()                      3         0.0364      0.012124    0.0425      0.014166    3.17     3.70     |
|   create_dof_constraints()              3         0.0081      0.002695    0.0130      0.004327    0.70     1.13     |
|   distribute_dofs()                     3         0.0018      0.000606    0.0185      0.006175    0.16     1.61     |
|   dof_indices()                         24234     0.0244      0.000001    0.0244      0.000001    2.12     2.12     |
|   old_dof_indices()                     1214      0.0005      0.000000    0.0005      0.000000    0.04     0.04     |
|   prepare_send_list()                   3         0.0000      0.000015    0.0000      0.000015    0.00     0.00     |
|   reinit()                              3         0.0030      0.000985    0.0030      0.000985    0.26     0.26     |
|                                                                                                                     |
| EquationSystems                                                                                                     |
|   build_discontinuous_solution_vector() 1         0.0066      0.006601    0.0078      0.007765    0.57     0.68     |
|                                                                                                                     |
| ExodusII_IO                                                                                                         |
|   write_nodal_data_discontinuous()      1         0.0041      0.004091    0.0041      0.004091    0.36     0.36     |
|                                                                                                                     |
| FE                                                                                                                  |
|   compute_shape_functions()             6758      0.0197      0.000003    0.0197      0.000003    1.72     1.72     |
|   init_shape_functions()                5582      0.0057      0.000001    0.0057      0.000001    0.49     0.49     |
|   inverse_map()                         18340     0.0323      0.000002    0.0323      0.000002    2.81     2.81     |
|                                                                                                                     |
| FEMap                                                                                                               |
|   compute_affine_map()                  6758      0.0072      0.000001    0.0072      0.000001    0.63     0.63     |
|   compute_face_map()                    2463      0.0031      0.000001    0.0031      0.000001    0.27     0.27     |
|   init_face_shape_functions()           5         0.0000      0.000005    0.0000      0.000005    0.00     0.00     |
|   init_reference_to_physical_map()      5582      0.0305      0.000005    0.0305      0.000005    2.66     2.66     |
|                                                                                                                     |
| JumpErrorEstimator                                                                                                  |
|   estimate_error()                      2         0.0054      0.002708    0.0466      0.023276    0.47     4.05     |
|                                                                                                                     |
| LocationMap                                                                                                         |
|   find()                                22344     0.0058      0.000000    0.0058      0.000000    0.50     0.50     |
|   init()                                6         0.0008      0.000139    0.0008      0.000139    0.07     0.07     |
|                                                                                                                     |
| Mesh                                                                                                                |
|   contract()                            2         0.0002      0.000089    0.0005      0.000267    0.02     0.05     |
|   find_neighbors()                      5         0.0085      0.001699    0.0173      0.003467    0.74     1.51     |
|   renumber_nodes_and_elem()             12        0.0012      0.000102    0.0012      0.000102    0.11     0.11     |
|                                                                                                                     |
| MeshCommunication                                                                                                   |
|   compute_hilbert_indices()             5         0.0067      0.001338    0.0067      0.001338    0.58     0.58     |
|   find_global_indices()                 5         0.0007      0.000141    0.0232      0.004645    0.06     2.02     |
|   parallel_sort()                       5         0.0022      0.000443    0.0074      0.001482    0.19     0.65     |
|                                                                                                                     |
| MeshRefinement                                                                                                      |
|   _coarsen_elements()                   4         0.0002      0.000060    0.0007      0.000177    0.02     0.06     |
|   _refine_elements()                    6         0.0183      0.003055    0.0864      0.014404    1.60     7.53     |
|   add_point()                           22344     0.0454      0.000002    0.0527      0.000002    3.96     4.59     |
|   make_coarsening_compatible()          4         0.0018      0.000438    0.0018      0.000438    0.15     0.15     |
|   make_refinement_compatible()          4         0.0001      0.000021    0.0097      0.002419    0.01     0.84     |
|                                                                                                                     |
| MetisPartitioner                                                                                                    |
|   partition()                           5         0.0088      0.001760    0.0339      0.006787    0.77     2.96     |
|                                                                                                                     |
| Parallel                                                                                                            |
|   allgather()                           19        0.0023      0.000120    0.0023      0.000120    0.20     0.20     |
|   broadcast()                           20        0.0001      0.000004    0.0001      0.000003    0.01     0.00     |
|   max(bool)                             16        0.0153      0.000955    0.0153      0.000955    1.33     1.33     |
|   max(scalar)                           13        0.0319      0.002456    0.0319      0.002456    2.78     2.78     |
|   max(vector)                           5         0.0002      0.000039    0.0002      0.000039    0.02     0.02     |
|   min(bool)                             8         0.0112      0.001406    0.0112      0.001406    0.98     0.98     |
|   min(scalar)                           2         0.0001      0.000025    0.0001      0.000025    0.00     0.00     |
|   min(vector)                           5         0.0017      0.000342    0.0017      0.000342    0.15     0.15     |
|   probe()                               170       0.0241      0.000142    0.0241      0.000142    2.10     2.10     |
|   receive()                             170       0.0003      0.000002    0.0244      0.000144    0.02     2.13     |
|   send()                                170       0.0002      0.000001    0.0002      0.000001    0.01     0.01     |
|   send_receive()                        180       0.0003      0.000002    0.0250      0.000139    0.03     2.18     |
|   sum()                                 25        0.0299      0.001197    0.0299      0.001197    2.61     2.61     |
|                                                                                                                     |
| Parallel::Request                                                                                                   |
|   wait()                                170       0.0001      0.000000    0.0001      0.000000    0.01     0.01     |
|                                                                                                                     |
| Partitioner                                                                                                         |
|   set_node_processor_ids()              6         0.0010      0.000160    0.0061      0.001023    0.08     0.53     |
|   set_parent_processor_ids()            5         0.0006      0.000122    0.0006      0.000122    0.05     0.05     |
|                                                                                                                     |
| PetscLinearSolver                                                                                                   |
|   solve()                               3         0.6418      0.213946    0.6418      0.213946    55.90    55.90    |
|                                                                                                                     |
| ProjectVector                                                                                                       |
|   operator()                            2         0.0211      0.010568    0.0422      0.021110    1.84     3.68     |
|                                                                                                                     |
| System                                                                                                              |
|   assemble()                            3         0.0573      0.019087    0.1315      0.043838    4.99     11.45    |
|   project_vector()                      2         0.0183      0.009157    0.0609      0.030460    1.60     5.31     |
|                                                                                                                     |
| XdrIO                                                                                                               |
|   read()                                1         0.0001      0.000136    0.0002      0.000154    0.01     0.01     |
 ---------------------------------------------------------------------------------------------------------------------
| Totals:                                 116704    1.1482                                          100.00            |
 ---------------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 6 ./miscellaneous_ex5-opt --element_type XYZ -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
