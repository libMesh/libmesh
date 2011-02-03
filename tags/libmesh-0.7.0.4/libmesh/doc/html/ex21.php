<?php $root=""; ?>
<?php require($root."navigation.php"); ?>
<html>
<head>
  <?php load_style($root); ?>
</head>
 
<body>
 
<?php make_navigation("ex21",$root)?>
 
<div class="content">
<a name="comments"></a> 
<div class = "comment">
<h1>Example 21 - Discontinuous Galerkin in the L-Shaped Domain</h1>

<br><br>By Lorenzo Botti

<br><br>This example is based on example 14, but uses an Interior Penalty
Discontinuous Galerkin formulation.


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
        #include "boundary_info.h"
        #include "gmv_io.h"
        #include "fe_interface.h"
        #include "getpot.h"
        #include "mesh_refinement.h"
        #include "error_vector.h"
        #include "kelly_error_estimator.h"
        #include "discontinuity_measure.h"
        
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
The \p PetscMatrix::add_matrix() members do this for us.
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
right-hand-side vector.  The \p PetscMatrix::add_matrix()
and \p PetscVector::add_vector() members do this for us.
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
          GetPot input_file("ex21.in");
        
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
Refine and coaarsen the flagged elements
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
to a GMV-formatted plot file.
</div>

<div class ="fragment">
<pre>
          GMVIO (mesh).write_discontinuous_gmv("lshaped_dg.gmv",equation_system,true);
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
  #include <B><FONT COLOR="#BC8F8F">&quot;boundary_info.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;gmv_io.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;fe_interface.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;getpot.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;mesh_refinement.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;error_vector.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;kelly_error_estimator.h&quot;</FONT></B>
  #include <B><FONT COLOR="#BC8F8F">&quot;discontinuity_measure.h&quot;</FONT></B>
  
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
  
    GetPot input_file(<B><FONT COLOR="#BC8F8F">&quot;ex21.in&quot;</FONT></B>);
  
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
  
    LinearImplicitSystem&amp; ellipticdg_system = equation_system.add_system&lt;LinearImplicitSystem&gt; (<B><FONT COLOR="#BC8F8F">&quot;EllipticDG&quot;</FONT></B>);
   
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
  
    GMVIO (mesh).write_discontinuous_gmv(<B><FONT COLOR="#BC8F8F">&quot;lshaped_dg.gmv&quot;</FONT></B>,equation_system,true);
  }
</pre> 
<a name="output"></a> 
<br><br><br> <h1> The console output of the program: </h1> 
<pre>
***************************************************************
* Running Example  mpirun -np 2 ./ex21-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
***************************************************************
 
  Beginning Solve 0
Number of elements: 219
 assembling elliptic dg system... done
System has: 768 degrees of freedom.
Linear solver converged at step: 26, final residual: 7.49373e-12
L2-Error is: 0.00666744
  Beginning Solve 1
Number of elements: 827
 assembling elliptic dg system... done
System has: 2896 degrees of freedom.
Linear solver converged at step: 34, final residual: 1.5477e-11
L2-Error is: 0.00264921
  Beginning Solve 2
Number of elements: 3003
 assembling elliptic dg system... done
System has: 10512 degrees of freedom.
Linear solver converged at step: 45, final residual: 1.75113e-11
L2-Error is: 0.0016323
************************************************************************************************************************
***             WIDEN YOUR WINDOW TO 120 CHARACTERS.  Use 'enscript -r -fCourier9' to print this document            ***
************************************************************************************************************************

---------------------------------------------- PETSc Performance Summary: ----------------------------------------------

./ex21-opt on a gcc-4.5-l named daedalus with 2 processors, by roystgnr Thu Feb  3 12:11:26 2011
Using Petsc Release Version 3.1.0, Patch 5, Mon Sep 27 11:51:54 CDT 2010

                         Max       Max/Min        Avg      Total 
Time (sec):           1.699e+00      1.07132   1.643e+00
Objects:              2.120e+02      1.00000   2.120e+02
Flops:                8.836e+08      1.10309   8.423e+08  1.685e+09
Flops/sec:            5.200e+08      1.02966   5.125e+08  1.025e+09
MPI Messages:         1.720e+02      1.00000   1.720e+02  3.440e+02
MPI Message Lengths:  5.762e+05      1.00000   3.350e+03  1.152e+06
MPI Reductions:       4.230e+02      1.00000

Flop counting convention: 1 flop = 1 real number operation of type (multiply/divide/add/subtract)
                            e.g., VecAXPY() for real vectors of length N --> 2N flops
                            and VecAXPY() for complex vectors of length N --> 8N flops

Summary of Stages:   ----- Time ------  ----- Flops -----  --- Messages ---  -- Message Lengths --  -- Reductions --
                        Avg     %Total     Avg     %Total   counts   %Total     Avg         %Total   counts   %Total 
 0:      Main Stage: 1.6426e+00 100.0%  1.6846e+09 100.0%  3.440e+02 100.0%  3.350e+03      100.0%  3.690e+02  87.2% 

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

VecMDot              105 1.0 1.0803e-02 2.4 7.79e+06 1.0 0.0e+00 0.0e+00 1.0e+02  0  1  0  0 25   0  1  0  0 28  1443
VecNorm              113 1.0 2.7018e-03 1.1 6.33e+05 1.0 0.0e+00 0.0e+00 1.1e+02  0  0  0  0 27   0  0  0  0 31   469
VecScale             110 1.0 1.7715e-04 1.1 3.10e+05 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0  3495
VecCopy               17 1.0 5.1022e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecSet               131 1.0 1.7643e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecAXPY               10 1.0 3.4406e-03 1.0 5.52e+04 1.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0    32
VecMAXPY             110 1.0 2.0900e-03 1.0 8.39e+06 1.0 0.0e+00 0.0e+00 0.0e+00  0  1  0  0  0   0  1  0  0  0  8025
VecAssemblyBegin      21 1.0 2.2638e-03 3.2 0.00e+00 0.0 0.0e+00 0.0e+00 5.1e+01  0  0  0  0 12   0  0  0  0 14     0
VecAssemblyEnd        21 1.0 1.7643e-05 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
VecScatterBegin      126 1.0 6.6376e-04 1.4 0.00e+00 0.0 2.4e+02 2.8e+03 0.0e+00  0  0 69 58  0   0  0 69 58  0     0
VecScatterEnd        126 1.0 8.7073e-02 7.4 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  3  0  0  0  0   3  0  0  0  0     0
VecNormalize         110 1.0 2.9879e-03 1.1 9.29e+05 1.0 0.0e+00 0.0e+00 1.1e+02  0  0  0  0 26   0  0  0  0 30   622
MatMult              110 1.0 9.7889e-02 4.4 1.67e+07 1.0 2.2e+02 2.8e+03 0.0e+00  4  2 64 54  0   4  2 64 54  0   340
MatSolve             110 1.0 1.3883e-01 1.1 2.16e+08 1.0 0.0e+00 0.0e+00 0.0e+00  8 25  0  0  0   8 25  0  0  0  3058
MatLUFactorNum         3 1.0 3.1846e-01 1.1 6.33e+08 1.1 0.0e+00 0.0e+00 0.0e+00 18 71  0  0  0  18 71  0  0  0  3745
MatILUFactorSym        3 1.0 4.5724e-01 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 3.0e+00 27  0  0  0  1  27  0  0  0  1     0
MatAssemblyBegin       6 1.0 8.5721e-0312.2 0.00e+00 0.0 1.8e+01 8.9e+03 1.2e+01  0  0  5 14  3   0  0  5 14  3     0
MatAssemblyEnd         6 1.0 2.2421e-03 1.0 0.00e+00 0.0 1.2e+01 6.2e+02 2.4e+01  0  0  3  1  6   0  0  3  1  7     0
MatGetRowIJ            3 1.0 1.9073e-06 2.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
MatGetOrdering         3 1.0 1.2898e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 1.2e+01  0  0  0  0  3   0  0  0  0  3     0
MatZeroEntries         9 1.0 5.8413e-04 1.0 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPGMRESOrthog       105 1.0 1.2863e-02 2.0 1.56e+07 1.0 0.0e+00 0.0e+00 1.0e+02  1  2  0  0 25   1  2  0  0 28  2424
KSPSetup               6 1.0 9.0837e-05 1.1 0.00e+00 0.0 0.0e+00 0.0e+00 0.0e+00  0  0  0  0  0   0  0  0  0  0     0
KSPSolve               3 1.0 9.5206e-01 1.0 8.84e+08 1.1 2.2e+02 2.8e+03 2.3e+02 58100 64 54 55  58100 64 54 63  1769
PCSetUp                6 1.0 7.7616e-01 1.1 6.33e+08 1.1 0.0e+00 0.0e+00 1.5e+01 45 71  0  0  4  45 71  0  0  4  1536
PCSetUpOnBlocks        3 1.0 7.7592e-01 1.1 6.33e+08 1.1 0.0e+00 0.0e+00 1.5e+01 45 71  0  0  4  45 71  0  0  4  1537
PCApply              110 1.0 1.3998e-01 1.1 2.16e+08 1.0 0.0e+00 0.0e+00 0.0e+00  8 25  0  0  0   8 25  0  0  0  3033
------------------------------------------------------------------------------------------------------------------------

Memory usage is given in bytes:

Object Type          Creations   Destructions     Memory  Descendants' Mem.
Reports information only for process 0.

--- Event Stage 0: Main Stage

                 Vec   134            134      2499424     0
         Vec Scatter    15             15        13020     0
           Index Set    36             36       126648     0
   IS L to G Mapping     3              3        34604     0
              Matrix    12             12     31495620     0
       Krylov Solver     6              6        56640     0
      Preconditioner     6              6         4224     0
========================================================================================================================
Average time to get PetscTime(): 9.53674e-08
Average time for MPI_Barrier(): 1.7643e-06
Average time for zero size MPI_Send(): 5.96046e-06
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
| Time:           Thu Feb  3 12:11:26 2011                         |
| OS:             Linux                                            |
| HostName:       daedalus                                         |
| OS Release:     2.6.32-26-generic                                |
| OS Version:     #46-Ubuntu SMP Tue Oct 26 16:47:18 UTC 2010      |
| Machine:        x86_64                                           |
| Username:       roystgnr                                         |
| Configuration:  ./configure run on Tue Feb  1 12:58:27 CST 2011  |
-------------------------------------------------------------------
 ------------------------------------------------------------------------------------------------------------
| libMesh Performance: Alive time=1.70801, Active time=1.59165                                               |
 ------------------------------------------------------------------------------------------------------------
| Event                          nCalls    Total Time  Avg Time    Total Time  Avg Time    % of Active Time  |
|                                          w/o Sub     w/o Sub     With Sub    With Sub    w/o S    With S   |
|------------------------------------------------------------------------------------------------------------|
|                                                                                                            |
|                                                                                                            |
| DofMap                                                                                                     |
|   add_neighbors_to_send_list() 3         0.0012      0.000399    0.0014      0.000482    0.08     0.09     |
|   compute_sparsity()           3         0.0472      0.015721    0.0636      0.021198    2.96     4.00     |
|   create_dof_constraints()     3         0.0003      0.000090    0.0003      0.000090    0.02     0.02     |
|   distribute_dofs()            3         0.0028      0.000937    0.0071      0.002371    0.18     0.45     |
|   dof_indices()                60982     0.0217      0.000000    0.0217      0.000000    1.36     1.36     |
|   old_dof_indices()            3512      0.0011      0.000000    0.0011      0.000000    0.07     0.07     |
|   prepare_send_list()          3         0.0001      0.000031    0.0001      0.000031    0.01     0.01     |
|   reinit()                     3         0.0041      0.001375    0.0041      0.001375    0.26     0.26     |
|                                                                                                            |
| FE                                                                                                         |
|   compute_affine_map()         20176     0.0209      0.000001    0.0209      0.000001    1.31     1.31     |
|   compute_face_map()           7172      0.0074      0.000001    0.0074      0.000001    0.47     0.47     |
|   compute_shape_functions()    20176     0.0084      0.000000    0.0084      0.000000    0.53     0.53     |
|   init_face_shape_functions()  7148      0.0164      0.000002    0.0164      0.000002    1.03     1.03     |
|   init_shape_functions()       16638     0.0741      0.000004    0.0741      0.000004    4.65     4.65     |
|   inverse_map()                87552     0.2234      0.000003    0.2234      0.000003    14.04    14.04    |
|                                                                                                            |
| JumpErrorEstimator                                                                                         |
|   estimate_error()             2         0.0154      0.007697    0.0616      0.030788    0.97     3.87     |
|                                                                                                            |
| LocationMap                                                                                                |
|   find()                       22344     0.0062      0.000000    0.0062      0.000000    0.39     0.39     |
|   init()                       6         0.0012      0.000196    0.0012      0.000196    0.07     0.07     |
|                                                                                                            |
| Mesh                                                                                                       |
|   contract()                   2         0.0003      0.000140    0.0006      0.000299    0.02     0.04     |
|   find_neighbors()             5         0.0101      0.002011    0.0101      0.002016    0.63     0.63     |
|   renumber_nodes_and_elem()    12        0.0009      0.000077    0.0009      0.000077    0.06     0.06     |
|                                                                                                            |
| MeshCommunication                                                                                          |
|   compute_hilbert_indices()    5         0.0057      0.001148    0.0057      0.001148    0.36     0.36     |
|   find_global_indices()        5         0.0010      0.000194    0.0074      0.001473    0.06     0.46     |
|   parallel_sort()              5         0.0005      0.000092    0.0005      0.000103    0.03     0.03     |
|                                                                                                            |
| MeshRefinement                                                                                             |
|   _coarsen_elements()          4         0.0004      0.000094    0.0004      0.000097    0.02     0.02     |
|   _refine_elements()           6         0.0180      0.003007    0.0418      0.006972    1.13     2.63     |
|   add_point()                  22344     0.0138      0.000001    0.0214      0.000001    0.87     1.35     |
|   make_coarsening_compatible() 4         0.0024      0.000610    0.0024      0.000610    0.15     0.15     |
|   make_refinement_compatible() 4         0.0001      0.000025    0.0001      0.000028    0.01     0.01     |
|                                                                                                            |
| MetisPartitioner                                                                                           |
|   partition()                  5         0.0056      0.001116    0.0130      0.002599    0.35     0.82     |
|                                                                                                            |
| Parallel                                                                                                   |
|   allgather()                  14        0.0001      0.000009    0.0001      0.000009    0.01     0.01     |
|   broadcast()                  20        0.0000      0.000002    0.0000      0.000001    0.00     0.00     |
|   max(bool)                    16        0.0009      0.000054    0.0009      0.000054    0.05     0.05     |
|   max(scalar)                  15        0.0002      0.000013    0.0002      0.000013    0.01     0.01     |
|   max(vector)                  5         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   min(bool)                    8         0.0000      0.000004    0.0000      0.000004    0.00     0.00     |
|   min(scalar)                  2         0.0000      0.000002    0.0000      0.000002    0.00     0.00     |
|   min(vector)                  5         0.0000      0.000008    0.0000      0.000008    0.00     0.00     |
|   probe()                      34        0.0001      0.000002    0.0001      0.000002    0.00     0.00     |
|   receive()                    34        0.0001      0.000003    0.0002      0.000005    0.01     0.01     |
|   send()                       34        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|   send_receive()               44        0.0001      0.000002    0.0003      0.000007    0.00     0.02     |
|   sum()                        20        0.0024      0.000119    0.0024      0.000119    0.15     0.15     |
|   wait()                       34        0.0000      0.000001    0.0000      0.000001    0.00     0.00     |
|                                                                                                            |
| Partitioner                                                                                                |
|   set_node_processor_ids()     6         0.0016      0.000270    0.0018      0.000292    0.10     0.11     |
|   set_parent_processor_ids()   5         0.0005      0.000093    0.0005      0.000093    0.03     0.03     |
|                                                                                                            |
| PetscLinearSolver                                                                                          |
|   solve()                      3         0.9633      0.321114    0.9633      0.321114    60.52    60.52    |
|                                                                                                            |
| ProjectVector                                                                                              |
|   operator()                   2         0.0157      0.007858    0.0697      0.034875    0.99     4.38     |
|                                                                                                            |
| System                                                                                                     |
|   assemble()                   3         0.0918      0.030593    0.3045      0.101497    5.77     19.13    |
|   project_vector()             2         0.0040      0.002012    0.0744      0.037194    0.25     4.67     |
|                                                                                                            |
| XdrIO                                                                                                      |
|   read()                       1         0.0001      0.000136    0.0001      0.000147    0.01     0.01     |
 ------------------------------------------------------------------------------------------------------------
| Totals:                        268434    1.5917                                          100.00            |
 ------------------------------------------------------------------------------------------------------------

 
***************************************************************
* Done Running Example  mpirun -np 2 ./ex21-opt -pc_type bjacobi -sub_pc_type ilu -sub_pc_factor_levels 4 -sub_pc_factor_zeropivot 0 -ksp_right_pc -log_summary
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
