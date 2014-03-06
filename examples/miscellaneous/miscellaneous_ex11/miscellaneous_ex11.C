/* The Next Great Finite Element Library. */
/* Copyright (C) 2003  Benjamin S. Kirk */

/* This library is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU Lesser General Public */
/* License as published by the Free Software Foundation; either */
/* version 2.1 of the License, or (at your option) any later version. */

/* This library is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU */
/* Lesser General Public License for more details. */

/* You should have received a copy of the GNU Lesser General Public */
/* License along with this library; if not, write to the Free Software */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */



 // <h1>Miscellaneous Example 11 - Using Loop Subdivision Shell Elements</h1>
 //
 // This example demonstrates how subdivision surface shell elements
 // are used, and how boundary conditions can be applied to them.  To
 // keep it simple, we solve the static deflection of a clamped, square,
 // linearly elastic Kirchhoff-Love plate subject to a uniform
 // load distribution.  Refer to Cirak et al., Int. J. Numer. Meth.
 // Engng. 2000; 47: 2039-2072, for a detailed description of what's
 // implemented.  In fact, this example follows that paper very closely.


// C++ include files that we need
#include <iostream>

// LibMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature.h"
#include "libmesh/node.h"
#include "libmesh/elem.h"
#include "libmesh/dof_map.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_value.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/vtk_io.h"
#include "libmesh/exodusII_io.h"

// These are the include files typically needed for subdivision elements.
#include "libmesh/face_tri3_subdivision.h"
#include "libmesh/mesh_subdivision_support.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Function prototype.  This is the function that will assemble
// the stiffness matrix and the right-hand-side vector ready
// for solution.
void assemble_shell (EquationSystems& es, const std::string& system_name);

// Begin the main program.
int main (int argc, char** argv)
{
  // Initialize libMesh.
  LibMeshInit init (argc, argv);

  // Skip this 3D example if libMesh was compiled as 1D/2D-only.
  libmesh_example_assert (3 == LIBMESH_DIM, "3D support");

  // Create a 2D mesh distributed across the default MPI communicator.
  Mesh mesh (init.comm(), 2);

  // Read the coarse square mesh.
  mesh.read ("square_mesh.off");

  // Resize the square plate to edge length L.
  const Real L = 100.;
  MeshTools::Modification::scale(mesh, L, L, L);

  // Quadrisect the mesh triangles a few times to obtain a
  // finer mesh.  Subdivision surface elements require the
  // refinement data to be removed afterwards.
  MeshRefinement mesh_refinement (mesh);
  mesh_refinement.uniformly_refine (3);
  MeshTools::Modification::flatten (mesh);

  // Write the mesh before the ghost elements are added.
#if defined(LIBMESH_HAVE_VTK)
  VTKIO(mesh).write ("without_ghosts.pvtu");
#endif
#if defined(LIBMESH_HAVE_EXODUS_API)
  ExodusII_IO(mesh).write ("without_ghosts.e");
#endif

  // Print information about the triangulated mesh to the screen.
  mesh.print_info();

  // Turn the triangulated mesh into a subdivision mesh
  // and add an additional row of "ghost" elements around
  // it in order to complete the extended local support of
  // the triangles at the boundaries.  If the second
  // argument is set to true, the outermost existing
  // elements are converted into ghost elements, and the
  // actual physical mesh is thus getting smaller.
  MeshTools::Subdivision::prepare_subdivision_mesh (mesh, false);

  // Print information about the subdivision mesh to the screen.
  mesh.print_info();

  // Write the mesh with the ghost elements added.
  // Compare this to the original mesh to see the difference.
#if defined(LIBMESH_HAVE_VTK)
  VTKIO(mesh).write ("with_ghosts.pvtu");
#endif
#if defined(LIBMESH_HAVE_EXODUS_API)
  ExodusII_IO(mesh).write ("with_ghosts.e");
#endif

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Create a linear implicit system named "Shell".
  LinearImplicitSystem & system = equation_systems.add_system<LinearImplicitSystem> ("Shell");

  // Add the three translational deformation variables
  // "u", "v", "w" to "Shell".  Since subdivision shell
  // elements meet the C1-continuity requirement, no
  // rotational or other auxiliary variables are needed.
  // Loop Subdivision Elements are always interpolated
  // by quartic box splines, hence the order must always
  // be \p FOURTH.
  system.add_variable ("u", FOURTH, SUBDIVISION);
  system.add_variable ("v", FOURTH, SUBDIVISION);
  system.add_variable ("w", FOURTH, SUBDIVISION);

  // Give the system a pointer to the matrix and rhs assembly
  // function.
  system.attach_assemble_function (assemble_shell);

  // Use the parameters of the equation systems object to
  // tell the shell system about the material properties, the
  // shell thickness, and the external load.
  const Real h  = 1.;
  const Real E  = 1.e7;
  const Real nu = 0.;
  const Real q  = 1.;
  equation_systems.parameters.set<Real> ("thickness")       = h;
  equation_systems.parameters.set<Real> ("young's modulus") = E;
  equation_systems.parameters.set<Real> ("poisson ratio")   = nu;
  equation_systems.parameters.set<Real> ("uniform load")    = q;

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Print information about the system to the screen.
  equation_systems.print_info();

  // Solve the linear system.
  system.solve();

  // After solving the system, write the solution to a VTK
  // or ExodusII output file ready for import in, e.g.,
  // Paraview.
#if defined(LIBMESH_HAVE_VTK)
  VTKIO(mesh).write_equation_systems ("out.pvtu", equation_systems);
#endif
#if defined(LIBMESH_HAVE_EXODUS_API)
  ExodusII_IO(mesh).write_equation_systems ("out.e", equation_systems);
#endif

  // Find the center node to measure the maximum deformation of the plate.
  Node* center_node = 0;
  Real nearest_dist_sq = mesh.point(0).size_sq();
  for (unsigned int nid=1; nid<mesh.n_nodes(); ++nid)
  {
    const Real dist_sq = mesh.point(nid).size_sq();
    if (dist_sq < nearest_dist_sq)
    {
      nearest_dist_sq = dist_sq;
      center_node = mesh.node_ptr(nid);
    }
  }

  // Finally, we evaluate the z-displacement "w" at the center node.
  const unsigned int w_var = system.variable_number ("w");
  unsigned int w_dof = center_node->dof_number (system.number(), w_var, 0);
  std::vector<Real> soln;
  system.current_local_solution->localize_to_one(soln);
  Real w = soln[w_dof];

  // The analytic solution for the maximum displacement of
  // a clamped square plate in pure bending, from Taylor,
  // Govindjee, Commun. Numer. Meth. Eng. 20, 757-765, 2004.
  const Real D = E * h*h*h / (12*(1-nu*nu));
  const Real w_analytic = 0.001265319 * L*L*L*L * q / D;

  // Print the finite element solution and the analytic
  // prediction of the maximum displacement of the clamped
  // square plate to the screen.
  std::cout << "z-displacement of the center point: " << w << std::endl;
  std::cout << "Analytic solution for pure bending: " << w_analytic << std::endl;

  // All done.
  return 0;
}

// We now define the matrix and rhs vector assembly function
// for the shell system.  This function implements the
// linear Kirchhoff-Love theory for thin shells.  At the
// end we also take into account the boundary conditions
// here, using the penalty method.
void assemble_shell (EquationSystems& es, const std::string& system_name)
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "Shell");

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();

  // Get a reference to the shell system object.
  LinearImplicitSystem & system = es.get_system<LinearImplicitSystem> ("Shell");

  // Get the shell parameters that we need during assembly.
  const Real h  = es.parameters.get<Real> ("thickness");
  const Real E  = es.parameters.get<Real> ("young's modulus");
  const Real nu = es.parameters.get<Real> ("poisson ratio");
  const Real q  = es.parameters.get<Real> ("uniform load");

  // Compute the membrane stiffness \p K and the bending
  // rigidity \p D from these parameters.
  const Real K = E * h     /     (1-nu*nu);
  const Real D = E * h*h*h / (12*(1-nu*nu));

  // Numeric ids corresponding to each variable in the system.
  const unsigned int u_var = system.variable_number ("u");
  const unsigned int v_var = system.variable_number ("v");
  const unsigned int w_var = system.variable_number ("w");

  // Get the Finite Element type for "u".  Note this will be
  // the same as the type for "v" and "w".
  FEType fe_type = system.variable_type (u_var);

  // Build a Finite Element object of the specified type.
  AutoPtr<FEBase> fe (FEBase::build(2, fe_type));

  // A Gauss quadrature rule for numerical integration.
  // For subdivision shell elements, a single Gauss point per
  // element is sufficient, hence we use extraorder = 0.
  const int extraorder = 0;
  AutoPtr<QBase> qrule (fe_type.default_quadrature_rule (2, extraorder));

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (qrule.get());

  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<Real>& JxW = fe->get_JxW();

  // The surface tangents in both directions at the quadrature points.
  const std::vector<RealGradient>& dxyzdxi  = fe->get_dxyzdxi();
  const std::vector<RealGradient>& dxyzdeta = fe->get_dxyzdeta();

  // The second partial derivatives at the quadrature points.
  const std::vector<RealGradient>& d2xyzdxi2    = fe->get_d2xyzdxi2();
  const std::vector<RealGradient>& d2xyzdeta2   = fe->get_d2xyzdeta2();
  const std::vector<RealGradient>& d2xyzdxideta = fe->get_d2xyzdxideta();

  // The element shape function and its derivatives evaluated at the
  // quadrature points.
  const std::vector<std::vector<Real> >&          phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  const std::vector<std::vector<RealTensor> >&  d2phi = fe->get_d2phi();

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  const DofMap & dof_map = system.get_dof_map();

  // Define data structures to contain the element stiffness matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  DenseSubMatrix<Number>
    Kuu(Ke), Kuv(Ke), Kuw(Ke),
    Kvu(Ke), Kvv(Ke), Kvw(Ke),
    Kwu(Ke), Kwv(Ke), Kww(Ke);

  DenseSubVector<Number>
    Fu(Fe),
    Fv(Fe),
    Fw(Fe);

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_u;
  std::vector<unsigned int> dof_indices_v;
  std::vector<unsigned int> dof_indices_w;

  // Now we will loop over all the elements in the mesh.  We will
  // compute the element matrix and right-hand-side contribution.
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem* elem = *el;

    // The ghost elements at the boundaries need to be excluded
    // here, as they don't belong to the physical shell,
    // but serve for a proper boundary treatment only.
    libmesh_assert_equal_to (elem->type(), TRI3SUBDIVISION);
    const Tri3Subdivision* sd_elem = static_cast<const Tri3Subdivision*> (elem);
    if (sd_elem->is_ghost())
      continue;

    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices (elem, dof_indices);
    dof_map.dof_indices (elem, dof_indices_u, u_var);
    dof_map.dof_indices (elem, dof_indices_v, v_var);
    dof_map.dof_indices (elem, dof_indices_w, w_var);

    const unsigned int n_dofs   = dof_indices.size();
    const unsigned int n_u_dofs = dof_indices_u.size();
    const unsigned int n_v_dofs = dof_indices_v.size();
    const unsigned int n_w_dofs = dof_indices_w.size();

    // Compute the element-specific data for the current
    // element.  This involves computing the location of the
    // quadrature points and the shape functions
    // (phi, dphi, d2phi) for the current element.
    fe->reinit (elem);

    // Zero the element matrix and right-hand side before
    // summing them.  We use the resize member here because
    // the number of degrees of freedom might have changed from
    // the last element.
    Ke.resize (n_dofs, n_dofs);
    Fe.resize (n_dofs);

    // Reposition the submatrices...  The idea is this:
    //
    //         -           -          -  -
    //        | Kuu Kuv Kuw |        | Fu |
    //   Ke = | Kvu Kvv Kvw |;  Fe = | Fv |
    //        | Kwu Kwv Kww |        | Fw |
    //         -           -          -  -
    //
    // The \p DenseSubMatrix.repostition () member takes the
    // (row_offset, column_offset, row_size, column_size).
    //
    // Similarly, the \p DenseSubVector.reposition () member
    // takes the (row_offset, row_size)
    Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
    Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
    Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);

    Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
    Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
    Kvw.reposition (v_var*n_v_dofs, w_var*n_v_dofs, n_v_dofs, n_w_dofs);

    Kwu.reposition (w_var*n_w_dofs, u_var*n_w_dofs, n_w_dofs, n_u_dofs);
    Kwv.reposition (w_var*n_w_dofs, v_var*n_w_dofs, n_w_dofs, n_v_dofs);
    Kww.reposition (w_var*n_w_dofs, w_var*n_w_dofs, n_w_dofs, n_w_dofs);

    Fu.reposition (u_var*n_u_dofs, n_u_dofs);
    Fv.reposition (v_var*n_u_dofs, n_v_dofs);
    Fw.reposition (w_var*n_u_dofs, n_w_dofs);

    // Now we will build the element matrix and right-hand-side.
    for (unsigned int qp=0; qp<qrule->n_points(); ++qp)
    {
      // First, we compute the external force resulting
      // from a load q distributed uniformly across the plate.
      // Since the load is supposed to be transverse to the plate,
      // it affects the z-direction, i.e. the "w" variable.
      for (unsigned int i=0; i<n_u_dofs; ++i)
        Fw(i) += JxW[qp] * phi[i][qp] * q;

      // Next, we assemble the stiffness matrix.  This is only valid
      // for the linear theory, i.e., for small deformations, where
      // reference and deformed surface metrics are indistinguishable.

      // Get the three surface basis vectors.
      const RealVectorValue & a1 = dxyzdxi[qp];
      const RealVectorValue & a2 = dxyzdeta[qp];
            RealVectorValue   a3 = a1.cross(a2);
      const Real jac = a3.size(); // the surface Jacobian
      libmesh_assert_greater (jac, 0);
      a3 /= jac; // the shell director a3 is normalized to unit length

      // Get the derivatives of the surface tangents.
      const RealVectorValue & a11 = d2xyzdxi2[qp];
      const RealVectorValue & a22 = d2xyzdeta2[qp];
      const RealVectorValue & a12 = d2xyzdxideta[qp];

      // Compute the three covariant components of the first
      // fundamental form of the surface.
      const RealVectorValue a(a1*a1, a2*a2, a1*a2);

      // The elastic H matrix in Voigt's notation, computed from the
      // covariant components of the first fundamental form rather
      // than the contravariant components, exploiting that the
      // contravariant first fundamental form is the inverse of the
      // covatiant first fundamental form (hence the determinant etc.).
      RealTensorValue H;
      H(0,0)          =  a(1) * a(1);
      H(0,1) = H(1,0) =   nu  * a(1) * a(0) + (1-nu) * a(2) * a(2);
      H(0,2) = H(2,0) = -a(1) * a(2);
      H(1,1)          =  a(0) * a(0);
      H(1,2) = H(2,1) = -a(0) * a(2);
      H(2,2)          = 0.5 * ((1-nu) * a(1) * a(0) + (1+nu) * a(2) * a(2));
      const Real det = a(0) * a(1) - a(2) * a(2);
      libmesh_assert_not_equal_to (det * det, 0);
      H /= det * det;

      // Precompute come cross products for the bending part below.
      const RealVectorValue a11xa2 = a11.cross(a2);
      const RealVectorValue a22xa2 = a22.cross(a2);
      const RealVectorValue a12xa2 = a12.cross(a2);
      const RealVectorValue a1xa11 =  a1.cross(a11);
      const RealVectorValue a1xa22 =  a1.cross(a22);
      const RealVectorValue a1xa12 =  a1.cross(a12);
      const RealVectorValue a2xa3  =  a2.cross(a3);
      const RealVectorValue a3xa1  =  a3.cross(a1);

      // Loop over all pairs of nodes I,J.
      for (unsigned int i=0; i<n_u_dofs; ++i)
      {
        for (unsigned int j=0; j<n_u_dofs; ++j)
        {
          // The membrane strain matrices in Voigt's notation.
          RealTensorValue MI, MJ;
          for (unsigned int k=0; k<3; ++k)
          {
            MI(0,k) = dphi[i][qp](0) * a1(k);
            MI(1,k) = dphi[i][qp](1) * a2(k);
            MI(2,k) = dphi[i][qp](1) * a1(k)
                    + dphi[i][qp](0) * a2(k);

            MJ(0,k) = dphi[j][qp](0) * a1(k);
            MJ(1,k) = dphi[j][qp](1) * a2(k);
            MJ(2,k) = dphi[j][qp](1) * a1(k)
                    + dphi[j][qp](0) * a2(k);
          }

          // The bending strain matrices in Voigt's notation.
          RealTensorValue BI, BJ;
          for (unsigned int k=0; k<3; ++k)
          {
            const Real term_ik = dphi[i][qp](0) * a2xa3(k)
                               + dphi[i][qp](1) * a3xa1(k);
            BI(0,k) = -d2phi[i][qp](0,0) * a3(k)
                      +(dphi[i][qp](0) * a11xa2(k)
                      + dphi[i][qp](1) * a1xa11(k)
                      + (a3*a11) * term_ik) / jac;
            BI(1,k) = -d2phi[i][qp](1,1) * a3(k)
                      +(dphi[i][qp](0) * a22xa2(k)
                      + dphi[i][qp](1) * a1xa22(k)
                      + (a3*a22) * term_ik) / jac;
            BI(2,k) = 2 * (-d2phi[i][qp](0,1) * a3(k)
                           +(dphi[i][qp](0) * a12xa2(k)
                           + dphi[i][qp](1) * a1xa12(k)
                           + (a3*a12) * term_ik) / jac);

            const Real term_jk = dphi[j][qp](0) * a2xa3(k)
                               + dphi[j][qp](1) * a3xa1(k);
            BJ(0,k) = -d2phi[j][qp](0,0) * a3(k)
                      +(dphi[j][qp](0) * a11xa2(k)
                      + dphi[j][qp](1) * a1xa11(k)
                      + (a3*a11) * term_jk) / jac;
            BJ(1,k) = -d2phi[j][qp](1,1) * a3(k)
                      +(dphi[j][qp](0) * a22xa2(k)
                      + dphi[j][qp](1) * a1xa22(k)
                      + (a3*a22) * term_jk) / jac;
            BJ(2,k) = 2 * (-d2phi[j][qp](0,1) * a3(k)
                           +(dphi[j][qp](0) * a12xa2(k)
                           + dphi[j][qp](1) * a1xa12(k)
                           + (a3*a12) * term_jk) / jac);
          }

          // The total stiffness matrix coupling the nodes
          // I and J is a sum of membrane and bending
          // contributions according to the following formula.
          const RealTensorValue KIJ = JxW[qp] * K * MI.transpose() * H * MJ
                                    + JxW[qp] * D * BI.transpose() * H * BJ;

          // Insert the components of the coupling stiffness
          // matrix \p KIJ into the corresponding directional
          // submatrices.
          Kuu(i,j) += KIJ(0,0);
          Kuv(i,j) += KIJ(0,1);
          Kuw(i,j) += KIJ(0,2);

          Kvu(i,j) += KIJ(1,0);
          Kvv(i,j) += KIJ(1,1);
          Kvw(i,j) += KIJ(1,2);

          Kwu(i,j) += KIJ(2,0);
          Kwv(i,j) += KIJ(2,1);
          Kww(i,j) += KIJ(2,2);
        }
      }

    } // end of the quadrature point qp-loop

    // The element matrix and right-hand-side are now built
    // for this element.  Add them to the global matrix and
    // right-hand-side vector.  The \p NumericMatrix::add_matrix()
    // and \p NumericVector::add_vector() members do this for us.
    system.matrix->add_matrix (Ke, dof_indices);
    system.rhs->add_vector    (Fe, dof_indices);
  } // end of non-ghost element loop

  // Next, we apply the boundary conditions.  In this case,
  // all boundaries are clamped by the penalty method, using
  // the special "ghost" nodes along the boundaries.  Note
  // that there are better ways to implement boundary conditions
  // for subdivision shells.  We use the simplest way here,
  // which is known to be overly restrictive and will lead to
  // a slightly too small deformation of the plate.
  el = mesh.active_local_elements_begin();

  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem* elem = *el;

    // For the boundary conditions, we only need to loop over
    // the ghost elements.
    libmesh_assert_equal_to (elem->type(), TRI3SUBDIVISION);
    const Tri3Subdivision* gh_elem = static_cast<const Tri3Subdivision*> (elem);
    if (!gh_elem->is_ghost())
      continue;

    // Find the side which is part of the physical plate boundary,
    // that is, the boundary of the original mesh without ghosts.
    for (unsigned int s=0; s<elem->n_sides(); ++s)
    {
      const Tri3Subdivision* nb_elem = static_cast<const Tri3Subdivision*> (elem->neighbor(s));
      if (nb_elem == NULL || nb_elem->is_ghost())
        continue;

      /*
       * Determine the four nodes involved in the boundary
       * condition treatment of this side.  The \p MeshTools::Subdiv
       * namespace provides lookup tables \p next and \p prev
       * for an efficient determination of the next and previous
       * nodes of an element, respectively.
       *
       *      n4
       *     /  \
       *    / gh \
       *  n2 ---- n3
       *    \ nb /
       *     \  /
       *      n1
       */
      Node* nodes [4]; // n1, n2, n3, n4
      nodes[1] = gh_elem->get_node(s); // n2
      nodes[2] = gh_elem->get_node(MeshTools::Subdivision::next[s]); // n3
      nodes[3] = gh_elem->get_node(MeshTools::Subdivision::prev[s]); // n4

      // The node in the interior of the domain, \p n1, is the
      // hardest to find.  Walk along the edges of element \p nb until
      // we have identified it.
      unsigned int n = 0;
      nodes[0] = nb_elem->get_node(0);
      while (nodes[0]->id() == nodes[1]->id() || nodes[0]->id() == nodes[2]->id())
        nodes[0] = nb_elem->get_node(++n);

      // The penalty value.  \f$ \frac{1}{\epsilon} \f$
      const Real penalty = 1.e10;

      // With this simple method, clamped boundary conditions are
      // obtained by penalizing the displacements of all four nodes.
      // This ensures that the displacement field vanishes on the
      // boundary side \p s.
      for (unsigned int n=0; n<4; ++n)
      {
        const unsigned int u_dof = nodes[n]->dof_number (system.number(), u_var, 0);
        const unsigned int v_dof = nodes[n]->dof_number (system.number(), v_var, 0);
        const unsigned int w_dof = nodes[n]->dof_number (system.number(), w_var, 0);
        system.matrix->add (u_dof, u_dof, penalty);
        system.matrix->add (v_dof, v_dof, penalty);
        system.matrix->add (w_dof, w_dof, penalty);
      }
    }
  } // end of ghost element loop
}
