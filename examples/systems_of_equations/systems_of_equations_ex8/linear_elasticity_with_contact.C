// Local includes
#include "linear_elasticity_with_contact.h"

// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <cmath>

// Various include files needed for the mesh & solver functionality.
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_compute_data.h"
#include "libmesh/petsc_matrix.h"
#ifdef LIBMESH_HAVE_PETSC
#  include "petscmat.h"
#endif
#include LIBMESH_INCLUDE_UNORDERED_SET

// The nonlinear solver and system we will be using
#include "libmesh/nonlinear_solver.h"
#include "libmesh/nonlinear_implicit_system.h"

using namespace libMesh;

bool operator< (QuadraturePointOnSideId const& a, QuadraturePointOnSideId const& b)
{
  if( (a._element_id == b._element_id) && (a._side_index == b._side_index) )
  {
    return (a._qp < b._qp);
  }
  else if(a._element_id == b._element_id)
  {
    return (a._side_index < b._side_index);
  }
  else
  {
    return (a._element_id < b._element_id);
  }
}

LinearElasticityWithContact::LinearElasticityWithContact
  (NonlinearImplicitSystem &sys_in,
   Real contact_penalty_in,
   Real contact_proximity_tol_in) :
    _sys(sys_in),
    _contact_penalty(contact_penalty_in),
    _contact_proximity_tol(contact_proximity_tol_in),
    _contains_point_tol(TOLERANCE),
    _augment_sparsity(_sys)
{}

AugmentSparsityOnContact& LinearElasticityWithContact::get_augment_sparsity()
{
  return _augment_sparsity;
}

void LinearElasticityWithContact::set_contact_penalty(
  Real contact_penalty_in)
{
  _contact_penalty = contact_penalty_in;
}

Real LinearElasticityWithContact::get_contact_penalty() const
{
  return _contact_penalty;
}

void LinearElasticityWithContact::clear_contact_data()
{
  _contact_intersection_data.clear();
}

void LinearElasticityWithContact::set_contact_data(
  dof_id_type element_id,
  unsigned char side_index,
  unsigned int qp,
  IntersectionPointData intersection_pt_data)
{
  QuadraturePointOnSideId elem_side_qp(element_id, side_index, qp);
  _contact_intersection_data[elem_side_qp] = intersection_pt_data;
}

bool LinearElasticityWithContact::is_contact_detected(
  dof_id_type element_id,
  unsigned char side_index,
  unsigned int qp)
{
  QuadraturePointOnSideId elem_side_qp(element_id, side_index, qp);

  std::map<QuadraturePointOnSideId, IntersectionPointData>::iterator it =
    _contact_intersection_data.find(elem_side_qp);

  return (it != _contact_intersection_data.end());
}

IntersectionPointData LinearElasticityWithContact::get_contact_data(
  dof_id_type element_id,
  unsigned char side_index,
  unsigned int qp)
{
  QuadraturePointOnSideId elem_side_qp(element_id, side_index, qp);

  std::map<QuadraturePointOnSideId, IntersectionPointData>::iterator it =
    _contact_intersection_data.find(elem_side_qp);

  if(it == _contact_intersection_data.end())
  {
    libmesh_error_msg("Contact intersection point not found.");
  }

  return it->second;
}

Real LinearElasticityWithContact::kronecker_delta(
  unsigned int i,
  unsigned int j)
{
  return i == j ? 1. : 0.;
}

Real LinearElasticityWithContact::elasticity_tensor(
  Real young_modulus,
  Real poisson_ratio,
  unsigned int i,
  unsigned int j,
  unsigned int k,
  unsigned int l)
{
  // Define the Lame constants
  const Real lambda_1 = (young_modulus*poisson_ratio)/((1.+poisson_ratio)*(1.-2.*poisson_ratio));
  const Real lambda_2 = young_modulus/(2.*(1.+poisson_ratio));

  return lambda_1 * kronecker_delta(i,j) * kronecker_delta(k,l) +
         lambda_2 * (kronecker_delta(i,k) * kronecker_delta(j,l) + kronecker_delta(i,l) * kronecker_delta(j,k));
}

void LinearElasticityWithContact::move_mesh(
  MeshBase& input_mesh,
  const NumericVector<Number>& input_solution)
{
  // Maintain a set of node ids that we've encountered.
  LIBMESH_BEST_UNORDERED_SET<dof_id_type> encountered_node_ids;

  // Localize input_solution so that we have the data to move all
  // elements (not just elements local to this processor).
  UniquePtr< NumericVector<Number> > localized_input_solution =
    NumericVector<Number>::build(input_solution.comm());
  localized_input_solution->init (
    input_solution.size(), false, SERIAL);
  input_solution.localize(*localized_input_solution);

  MeshBase::const_element_iterator       el     = input_mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = input_mesh.active_elements_end();

  for ( ; el != end_el; ++el)
  {
    Elem* elem = *el;
    Elem* orig_elem = _sys.get_mesh().elem(elem->id());

    for(unsigned int node_id=0; node_id<elem->n_nodes(); node_id++)
    {
      Node* node = elem->get_node(node_id);

      if(encountered_node_ids.find(node->id()) != encountered_node_ids.end())
      {
        continue;
      }
      encountered_node_ids.insert(node->id());

      std::vector<std::string> uvw_names(3);
      uvw_names[0] = "u";
      uvw_names[1] = "v";
      uvw_names[2] = "w";

      {
        // Inverse map node to reference element
        // Get local coordinates to feed these into compute_data().
        // Note that the fe_type can safely be used from the 0-variable,
        // since the inverse mapping is the same for all FEFamilies.
        const Point reference_point (
          FEInterface::inverse_map (
            elem->dim(),
            _sys.get_dof_map().variable_type(0),
            elem,
            *node));

        Point uvw;
        for (unsigned int index=0; index<uvw_names.size(); index++)
        {
          const unsigned int var = _sys.variable_number(uvw_names[index]);
          const FEType& fe_type = _sys.get_dof_map().variable_type(var);

          FEComputeData data (_sys.get_equation_systems(), reference_point);

          FEInterface::compute_data(elem->dim(),
                                    fe_type,
                                    elem,
                                    data);

          std::vector<dof_id_type> dof_indices_var;
          _sys.get_dof_map().dof_indices (orig_elem, dof_indices_var, var);

          for (unsigned int i=0; i<dof_indices_var.size(); i++)
            {
              Number value = (*localized_input_solution)(dof_indices_var[i]) * data.shape[i];

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
              // We explicitly store the real part in uvw
              uvw(index) += value.real();
#else
              uvw(index) += value;
#endif
            }
        }

        // Update the node's location
        *node += uvw;
      }
    }
  }
}

void LinearElasticityWithContact::residual_and_jacobian (
  const NumericVector<Number>& soln,
  NumericVector<Number>* residual,
  SparseMatrix<Number>* jacobian,
  NonlinearImplicitSystem& /*sys*/)
{
  EquationSystems& es = _sys.get_equation_systems();
  const Real young_modulus = es.parameters.get<Real>("young_modulus");
  const Real poisson_ratio = es.parameters.get<Real>("poisson_ratio");
  const Real forcing_magnitude = es.parameters.get<Real>("forcing_magnitude");

  const MeshBase& mesh = _sys.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  const unsigned int u_var = _sys.variable_number ("u");

  DofMap& dof_map = _sys.get_dof_map();

  FEType fe_type = dof_map.variable_type(u_var);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
  QGauss qface (dim-1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule (&qface);

  UniquePtr<FEBase> fe_neighbor_face (FEBase::build(dim, fe_type));
  fe_neighbor_face->attach_quadrature_rule (&qface);

  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  const std::vector<Real>& JxW_face = fe_face->get_JxW();
  const std::vector<std::vector<Real> >& phi_face = fe_face->get_phi();
  const std::vector<Point>& face_normals = fe_face->get_normals();
  const std::vector<Point>& face_xyz = fe_face->get_xyz();

  const std::vector<std::vector<Real> >& phi_neighbor_face = fe_neighbor_face->get_phi();

  // 1. Move mesh_clone based on soln.
  // 2. Compute and store all contact forces.
  // 3. Augment the sparsity pattern.
  {
    UniquePtr<MeshBase> mesh_clone = mesh.clone();
    move_mesh(*mesh_clone, soln);

    _augment_sparsity.clear_contact_element_map();
    clear_contact_data();

    MeshBase::const_element_iterator       el     = mesh_clone->active_elements_begin();
    const MeshBase::const_element_iterator end_el = mesh_clone->active_elements_end();

    for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      for (unsigned int side=0; side<elem->n_sides(); side++)
      {
        if (elem->neighbor(side) == NULL)
        {
          bool on_lower_contact_surface =
            mesh_clone->get_boundary_info().has_boundary_id
                (elem, side, CONTACT_BOUNDARY_LOWER);

          bool on_upper_contact_surface =
            mesh_clone->get_boundary_info().has_boundary_id
                (elem, side, CONTACT_BOUNDARY_UPPER);

          if( on_lower_contact_surface && on_upper_contact_surface )
          {
            libmesh_error_msg("Should not be on both surfaces at the same time");
          }

          if( on_lower_contact_surface || on_upper_contact_surface )
          {
            fe_face->reinit(elem, side);

            // Let's stash xyz and normals because we reinit on other_elem below
            std::vector<Point> face_normals_stashed = fe_face->get_normals();
            std::vector<Point> face_xyz_stashed = fe_face->get_xyz();

            for (unsigned int qp=0; qp<qface.n_points(); qp++)
            {
              Point line_point = face_xyz_stashed[qp];
              Point line_direction = face_normals_stashed[qp];

              // find an element which the line intersects, based on the plane
              // defined by the normal at the centroid of the other element
              bool found_other_elem = false;

              // Note that here we loop over all elements (not just local elements)
              // to be sure we find the appropriate other element.
              MeshBase::const_element_iterator       other_el     = mesh_clone->active_elements_begin();
              const MeshBase::const_element_iterator other_end_el = mesh_clone->active_elements_end();
              for ( ; other_el != other_end_el; ++other_el)
              {
                const Elem* other_elem = *other_el;

                if( other_elem->close_to_point(line_point, _contact_proximity_tol) )
                {
                  for (unsigned int other_side=0; other_side<other_elem->n_sides(); other_side++)
                    if (other_elem->neighbor(other_side) == NULL)
                    {
                      boundary_id_type other_surface_id =
                        on_lower_contact_surface ?
                          CONTACT_BOUNDARY_UPPER : CONTACT_BOUNDARY_LOWER;

                      if( mesh_clone->get_boundary_info().has_boundary_id
                            (other_elem, other_side, other_surface_id) )
                      {
                        UniquePtr<Elem> other_side_elem = other_elem->build_side(other_side);

                        // Define a plane based on the normal at the centroid of other_side_elem
                        // and check where line_point + s * line_direction (s \in R) intersects
                        // this plane.

                        const Point reference_centroid (
                          FEInterface::inverse_map (
                            other_elem->dim(),
                            _sys.get_dof_map().variable_type(0),
                            other_elem,
                            other_side_elem->centroid()));

                        std::vector<Point> reference_centroid_vector;
                        reference_centroid_vector.push_back(reference_centroid);

                        fe_face->reinit(
                          other_elem,
                          other_side,
                          TOLERANCE,
                          &reference_centroid_vector);

                        Point plane_normal = face_normals[0];
                        Point plane_point = face_xyz[0];

                        // line_direction.dot(plane_normal) == 0.0 if the line and plane are
                        // parallel. Ignore this case since it should give zero contact force.
                        if( (line_direction * plane_normal) != 0. )
                        {
                          // The signed distance between the line and the plane
                          Real signed_distance =
                            ( (plane_point - line_point) * plane_normal) /
                              (line_direction * plane_normal);

                          Point intersection_point = line_point + signed_distance * line_direction;

                          // Note that signed_distance = (intersection_point - line_point) dot line_direction
                          // since line_direction is a unit vector.

                          if(other_side_elem->close_to_point(intersection_point, _contains_point_tol))
                          {
                            // If the signed distance is negative then we have overlapping elements
                            // i.e. we have detected contact.
                            if(signed_distance < 0.0)
                            {
                              // We need to store the intersection point, and the element/side
                              // that it belongs to. We can use this to calculate the contact
                              // force later on.

                              std::vector<Point> intersection_point_vec;
                              intersection_point_vec.push_back(intersection_point);
                              std::vector<Point> inverse_intersection_point_vec;

                              FEInterface::inverse_map(
                                other_elem->dim(),
                                fe->get_fe_type(),
                                other_elem,
                                intersection_point_vec,
                                inverse_intersection_point_vec);

                              IntersectionPointData contact_intersection_pt_data(
                                other_elem->id(),
                                other_side,
                                intersection_point,
                                inverse_intersection_point_vec[0],
                                line_point,
                                line_direction);

                              set_contact_data(
                                elem->id(),
                                side,
                                qp,
                                contact_intersection_pt_data);

                              // We also need to keep track of which elements
                              // are coupled in order to augment the sparsity pattern
                              // appropriately later on.
                              _augment_sparsity.add_contact_element(
                                elem->id(),
                                other_elem->id());
                            }

                            // If we've found an element that contains
                            // intersection_point then we're done for
                            // the current quadrature point hence break
                            // out to the next qp.
                            found_other_elem = true;
                            break;
                          }
                        }
                      }
                    }
                }

                if(found_other_elem)
                {
                  break;
                }
              } // end for other_el
            } // end for qp
          } // end if on_contact_surface
        } // end if nieghbor(side_) != NULL
      } // end for side
    } // end for el
  }

  // Clear the Jacobian matrix and reinitialize it so that
  // we get the updated sparsity pattern
  if(jacobian)
  {
    dof_map.clear_sparsity();
    dof_map.compute_sparsity(mesh);

#ifdef LIBMESH_HAVE_PETSC
    PetscMatrix<Number>* petsc_jacobian = cast_ptr<PetscMatrix<Number>*>(jacobian);
    petsc_jacobian->update_preallocation_and_zero();
#else
    libmesh_error();
#endif
  }

  if(residual)
  {
    residual->zero();
  }

  // Do jacobian and residual assembly, including contact forces
  DenseVector<Number> Re;
  DenseSubVector<Number> Re_var[3] =
    {DenseSubVector<Number>(Re), DenseSubVector<Number>(Re), DenseSubVector<Number>(Re)};

  DenseMatrix<Number> Ke;
  DenseSubMatrix<Number> Ke_var[3][3] =
    {
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}
    };

  DenseMatrix<Number> Ke_en;
  DenseSubMatrix<Number> Ke_var_en[3][3] =
    {
      {DenseSubMatrix<Number>(Ke_en), DenseSubMatrix<Number>(Ke_en), DenseSubMatrix<Number>(Ke_en)},
      {DenseSubMatrix<Number>(Ke_en), DenseSubMatrix<Number>(Ke_en), DenseSubMatrix<Number>(Ke_en)},
      {DenseSubMatrix<Number>(Ke_en), DenseSubMatrix<Number>(Ke_en), DenseSubMatrix<Number>(Ke_en)}
    };

  std::vector<dof_id_type> dof_indices;
  std::vector< std::vector<dof_id_type> > dof_indices_var(3);

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
  {
    const Elem* elem = *el;

    dof_map.dof_indices (elem, dof_indices);
    for(unsigned int var=0; var<3; var++)
    {
      dof_map.dof_indices (elem, dof_indices_var[var], var);
    }

    const unsigned int n_dofs   = dof_indices.size();
    const unsigned int n_var_dofs = dof_indices_var[0].size();

    fe->reinit (elem);

    Re.resize (n_dofs);
    for(unsigned int var=0; var<3; var++)
    {
      Re_var[var].reposition (var*n_var_dofs, n_var_dofs);
    }

    Ke.resize (n_dofs,n_dofs);
    for(unsigned int var_i=0; var_i<3; var_i++)
      for(unsigned int var_j=0; var_j<3; var_j++)
      {
        Ke_var[var_i][var_j].reposition (var_i*n_var_dofs, var_j*n_var_dofs, n_var_dofs, n_var_dofs);
      }

    for (unsigned int qp=0; qp<qrule.n_points(); qp++)
    {
      DenseMatrix<Number> grad_u(3,3);
      for(unsigned int var_i=0; var_i<3; var_i++)
      {
        for(unsigned int var_j=0; var_j<3; var_j++)
        {
          for (unsigned int j=0; j<n_var_dofs; j++)
          {
            // Row is variable u, v, or w column is x, y, or z
            grad_u(var_i,var_j) += dphi[j][qp](var_j)*soln(dof_indices_var[var_i][j]);
          }
        }
      }

      // - C_ijkl u_k,l v_i,j
      for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
      {
        for(unsigned int i=0; i<3; i++)
          for(unsigned int j=0; j<3; j++)
            for(unsigned int k=0; k<3; k++)
              for(unsigned int l=0; l<3; l++)
                {
                  Re_var[i](dof_i) -= JxW[qp] *
                    elasticity_tensor(young_modulus,poisson_ratio,i,j,k,l) *
                    grad_u(k,l) * dphi[dof_i][qp](j);
                }
      }

      if( (elem->subdomain_id() == TOP_SUBDOMAIN) )
      {
        // assemble \int_Omega f_i v_i \dx
        DenseVector<Number> f_vec(3);
        f_vec(0) =  forcing_magnitude/10.;
        f_vec(1) =  0.0;
        f_vec(2) = -forcing_magnitude;
        for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
        {
          for(unsigned int i=0; i<3; i++)
          {
            Re_var[i](dof_i) += JxW[qp] *
              ( f_vec(i) * phi[dof_i][qp] );
          }
        }
      }

      // assemble \int_Omega C_ijkl u_k,l v_i,j \dx
      for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
        for (unsigned int dof_j=0; dof_j<n_var_dofs; dof_j++)
        {
          for(unsigned int i=0; i<3; i++)
            for(unsigned int j=0; j<3; j++)
              for(unsigned int k=0; k<3; k++)
                for(unsigned int l=0; l<3; l++)
                {
                  Ke_var[i][k](dof_i,dof_j) -= JxW[qp] *
                    elasticity_tensor(young_modulus,poisson_ratio,i,j,k,l) *
                    dphi[dof_j][qp](l) *
                    dphi[dof_i][qp](j);
                }
        }
    }

    // Add contribution due to contact penalty forces
    for (unsigned int side=0; side<elem->n_sides(); side++)
      if (elem->neighbor(side) == NULL)
      {
        bool on_lower_contact_surface =
          mesh.get_boundary_info().has_boundary_id
              (elem, side, CONTACT_BOUNDARY_LOWER);

        bool on_upper_contact_surface =
          mesh.get_boundary_info().has_boundary_id
              (elem, side, CONTACT_BOUNDARY_UPPER);

        if( on_lower_contact_surface && on_upper_contact_surface )
        {
          libmesh_error_msg("Should not be on both surfaces at the same time");
        }

        if( on_lower_contact_surface || on_upper_contact_surface )
        {
          fe_face->reinit(elem, side);
          for (unsigned int qp=0; qp<qface.n_points(); qp++)
          {
            bool contact_detected =
              is_contact_detected(
                elem->id(),
                side,
                qp);

            if(contact_detected)
            {
              IntersectionPointData intersection_pt_data =
                get_contact_data(
                  elem->id(),
                  side,
                  qp);
              Point intersection_point = intersection_pt_data._intersection_point;
              Point line_point = intersection_pt_data._line_point;
              Point line_direction = intersection_pt_data._line_direction;

              // signed_distance = (intersection_point - line_point) dot line_direction,
              // hence we use this to get the contact force
              Real contact_force =
                get_contact_penalty() *
                ( (intersection_point - line_point) * line_direction );

              for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
              {
                for(unsigned int i=0; i<3; i++)
                {
                  Re_var[i](dof_i) += JxW_face[qp] *
                    ( contact_force * face_normals[qp](i) * phi_face[dof_i][qp] );
                }
              }

              // Differentiate contact_force wrt solution coefficients to get
              // the Jacobian entries.
              //
              // Note that intersection_point and line_point
              // are linear functions of the solution.
              //
              // Also, line_direction is a function of the solution, but let's
              // neglect that for now to make it easier to get the Jacobian.
              // This is probably fine anyway, since this approximation will
              // have negligible effect as we approach convergence.

              // dofs local to this element, due to differentiating line_point
              for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
                for (unsigned int dof_j=0; dof_j<n_var_dofs; dof_j++)
                {
                  for(unsigned int i=0; i<3; i++)
                    for(unsigned int j=0; j<3; j++)
                    {
                      Ke_var[i][j](dof_i,dof_j) += JxW_face[qp] *
                        ( get_contact_penalty() * (-phi_face[dof_j][qp] * line_direction(j)) *
                           face_normals[qp](i) * phi_face[dof_i][qp] );
                    }
                }

              // contribution due to dofs on the remote element, due to
              // differentiating intersection_point
              {
                dof_id_type neighbor_elem_id = intersection_pt_data._neighbor_element_id;
                const Elem* contact_neighbor = mesh.elem(neighbor_elem_id);

                unsigned char neighbor_side_index = intersection_pt_data._neighbor_side_index;

                std::vector<Point> inverse_intersection_point_vec;
                inverse_intersection_point_vec.push_back(
                  intersection_pt_data._inverse_mapped_intersection_point);

                fe_neighbor_face->reinit(
                  contact_neighbor,
                  neighbor_side_index,
                  TOLERANCE,
                  &inverse_intersection_point_vec);

                std::vector<dof_id_type> neighbor_dof_indices;
                std::vector< std::vector<unsigned int> > neighbor_dof_indices_var(3);
                dof_map.dof_indices (contact_neighbor, neighbor_dof_indices);
                for(unsigned int var=0; var<3; var++)
                {
                  dof_map.dof_indices (contact_neighbor, neighbor_dof_indices_var[var], var);
                }
                const unsigned int n_neighbor_dofs = neighbor_dof_indices.size();
                const unsigned int n_neighbor_var_dofs = neighbor_dof_indices_var[0].size();

                Ke_en.resize (n_dofs,n_neighbor_dofs);
                for(unsigned int var_i=0; var_i<3; var_i++)
                  for(unsigned int var_j=0; var_j<3; var_j++)
                  {
                    Ke_var_en[var_i][var_j].reposition(
                      var_i*n_var_dofs,
                      var_j*n_neighbor_var_dofs,
                      n_var_dofs,
                      n_neighbor_var_dofs);
                  }

                for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
                  for (unsigned int dof_j=0; dof_j<n_neighbor_var_dofs; dof_j++)
                  {
                    for(unsigned int i=0; i<3; i++)
                      for(unsigned int j=0; j<3; j++)
                      {
                        Ke_var_en[i][j](dof_i,dof_j) += JxW_face[qp] *
                          ( get_contact_penalty() * (phi_neighbor_face[dof_j][0] * line_direction(j)) *
                            face_normals[qp](i) * phi_face[dof_i][qp] );
                      }
                  }

                if(jacobian)
                {
                  jacobian->add_matrix(Ke_en,dof_indices,neighbor_dof_indices);
                }
              }

            }

          }
        }
      }

    dof_map.constrain_element_matrix_and_vector (Ke, Re, dof_indices);

    if(jacobian)
    {
      jacobian->add_matrix (Ke, dof_indices);
    }

    if(residual)
    {
      residual->add_vector (Re, dof_indices);
    }
  }

}

void LinearElasticityWithContact::compute_stresses()
{
  EquationSystems& es = _sys.get_equation_systems();
  const Real young_modulus = es.parameters.get<Real>("young_modulus");
  const Real poisson_ratio = es.parameters.get<Real>("poisson_ratio");

  const MeshBase& mesh = _sys.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  unsigned int displacement_vars[3];
  displacement_vars[0] = _sys.variable_number ("u");
  displacement_vars[1] = _sys.variable_number ("v");
  displacement_vars[2] = _sys.variable_number ("w");
  const unsigned int u_var = _sys.variable_number ("u");

  const DofMap& dof_map = _sys.get_dof_map();
  FEType fe_type = dof_map.variable_type(u_var);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();

  // Also, get a reference to the ExplicitSystem
  ExplicitSystem& stress_system = es.get_system<ExplicitSystem>("StressSystem");
  const DofMap& stress_dof_map = stress_system.get_dof_map();
  unsigned int sigma_vars[6];
  sigma_vars[0] = stress_system.variable_number ("sigma_00");
  sigma_vars[1] = stress_system.variable_number ("sigma_01");
  sigma_vars[2] = stress_system.variable_number ("sigma_02");
  sigma_vars[3] = stress_system.variable_number ("sigma_11");
  sigma_vars[4] = stress_system.variable_number ("sigma_12");
  sigma_vars[5] = stress_system.variable_number ("sigma_22");
  unsigned int vonMises_var = stress_system.variable_number ("vonMises");

  // Storage for the stress dof indices on each element
  std::vector< std::vector<dof_id_type> > dof_indices_var(_sys.n_vars());
  std::vector<dof_id_type> stress_dof_indices_var;

  // To store the stress tensor on each element
  DenseMatrix<Number> elem_avg_stress_tensor(3,3);

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      for(unsigned int var=0; var<3; var++)
        {
          dof_map.dof_indices (elem, dof_indices_var[var], displacement_vars[var]);
        }

      const unsigned int n_var_dofs = dof_indices_var[0].size();

      fe->reinit (elem);

      // clear the stress tensor
      elem_avg_stress_tensor.resize(3,3);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          DenseMatrix<Number> grad_u(3,3);
          for(unsigned int var_i=0; var_i<3; var_i++)
            for(unsigned int var_j=0; var_j<3; var_j++)
              {
                for (unsigned int j=0; j<n_var_dofs; j++)
                  {
                    // Row is variable u1, u2, or u3, column is x, y, or z
                    grad_u(var_i,var_j) += dphi[j][qp](var_j)*
                      _sys.current_solution(dof_indices_var[var_i][j]);
                  }
              }

          DenseMatrix<Number> stress_tensor(3,3);
          for(unsigned int i=0; i<3; i++)
            for(unsigned int j=0; j<3; j++)
              for(unsigned int k=0; k<3; k++)
                for(unsigned int l=0; l<3; l++)
                  {
                    stress_tensor(i,j) +=
                      elasticity_tensor(young_modulus,poisson_ratio,i,j,k,l) *
                      grad_u(k,l);
                  }

          // We want to plot the average stress on each element, hence
          // we integrate stress_tensor
          elem_avg_stress_tensor.add(JxW[qp], stress_tensor);
        }

      // Get the average stress per element by dividing by volume
      elem_avg_stress_tensor.scale(1./elem->volume());

      // load elem_sigma data into stress_system
      unsigned int stress_var_index = 0;
      for(unsigned int i=0; i<3; i++)
        for(unsigned int j=i; j<3; j++)
          {
            stress_dof_map.dof_indices (elem, stress_dof_indices_var, sigma_vars[stress_var_index]);

            // We are using CONSTANT MONOMIAL basis functions, hence we only need to get
            // one dof index per variable
            dof_id_type dof_index = stress_dof_indices_var[0];

            if( (stress_system.solution->first_local_index() <= dof_index) &&
                (dof_index < stress_system.solution->last_local_index()) )
              {
                stress_system.solution->set(dof_index, elem_avg_stress_tensor(i,j));
              }

            stress_var_index++;
          }

      // Also, the von Mises stress
      Number vonMises_value = std::sqrt( 0.5*( pow(elem_avg_stress_tensor(0,0) - elem_avg_stress_tensor(1,1),2.) +
                                               pow(elem_avg_stress_tensor(1,1) - elem_avg_stress_tensor(2,2),2.) +
                                               pow(elem_avg_stress_tensor(2,2) - elem_avg_stress_tensor(0,0),2.) +
                                               6.*(pow(elem_avg_stress_tensor(0,1),2.) +
                                                   pow(elem_avg_stress_tensor(1,2),2.) +
                                                   pow(elem_avg_stress_tensor(2,0),2.))
                                               ) );
      stress_dof_map.dof_indices (elem, stress_dof_indices_var, vonMises_var);
      dof_id_type dof_index = stress_dof_indices_var[0];
      if( (stress_system.solution->first_local_index() <= dof_index) &&
          (dof_index < stress_system.solution->last_local_index()) )
        {
          stress_system.solution->set(dof_index, vonMises_value);
        }
    }

  // Should call close and update when we set vector entries directly
  stress_system.solution->close();
  stress_system.update();
}
