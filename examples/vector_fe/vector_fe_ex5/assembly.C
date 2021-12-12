#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/mesh_base.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/fe_type.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/fe.h"
#include "libmesh/tensor_value.h"
#include "libmesh/equation_systems.h"
#include "libmesh/parameters.h"
#include "libmesh/elem.h"
#include "libmesh/fe_interface.h"
#include "libmesh/elem_side_builder.h"

#include <memory>

using namespace libMesh;

// Function prototype for the exact solution.
Real exact_solution(const int component, const Real x, const Real y, const Real z = 0.);

// Forcing function
Real forcing_function(const int component, const Real x, const Real y, const Real z = 0.);

void
compute_residual(const NumericVector<Number> & X,
                 NumericVector<Number> & R,
                 NonlinearImplicitSystem & system)
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to(system.name(), "Poisson");

  // Get the DG parameters
  auto & es = system.get_equation_systems();
  const auto sigma = es.parameters.get<Real>("sigma");
  const auto epsilon = es.parameters.get<Real>("epsilon");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = system.get_mesh();

  // The dimension that we are running
  const auto dim = mesh.mesh_dimension();

  // A reference to the  DofMap object for this system.  The  DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the  DofMap
  // in future examples.
  const DofMap & dof_map = system.get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.
  // Note that FEVectorBase is a typedef for the templated FE
  // class.
  std::unique_ptr<FEVectorBase> fe(FEVectorBase::build(dim, fe_type));

  // An automatically determined Gauss quadrature rule for numerical integration.
  QGauss qrule(dim, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule(&qrule);

  // Declare a special finite element object for face integration.
  std::unique_ptr<FEVectorBase> fe_face(FEVectorBase::build(dim, fe_type));
  // And for neighbor integration
  std::unique_ptr<FEVectorBase> fe_neighbor_face(FEVectorBase::build(dim, fe_type));

  // Boundary integration requires one quadrature rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  QGauss qface(dim - 1, fe_type.default_quadrature_order());

  // Tell the face finite element objects to use our
  // quadrature rule.
  fe_face->attach_quadrature_rule(&qface);
  fe_neighbor_face->attach_quadrature_rule(&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.
  const auto & JxW = fe->get_JxW();

  // element integration points
  const auto & xyz = fe->get_xyz();

  // The element shape function values evaluated at the quadrature points
  const auto & phi = fe->get_phi();

  // The element shape function gradients evaluated at the quadrature points.
  const auto & dphi = fe->get_dphi();

  // face integration points
  const auto & face_xyz = fe_face->get_xyz();

  // Face shape function values
  const auto & phi_face = fe_face->get_phi();

  // Face shape function gradients
  const auto & dphi_face = fe_face->get_dphi();

  // Neighbor shape function values
  const auto & phi_neighbor = fe_neighbor_face->get_phi();

  // Neighbor face shape function gradients
  const auto & dphi_neighbor = fe_neighbor_face->get_dphi();

  // face normals
  const auto & normals = fe_face->get_normals();

  // face JxW
  const auto & JxW_face = fe_face->get_JxW();

  DenseVector<Number> Fe;
  DenseVector<Number> Fn;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_neighbor;

  /// Vectors to hold the local solution degree of freedom values
  std::vector<Number> dof_u;
  std::vector<Number> dof_u_neighbor;

  /// Vector to hold the local solution
  std::vector<VectorValue<Number>> u;
  std::vector<TensorValue<Number>> grad_u;
  std::vector<VectorValue<Number>> u_neighbor;
  std::vector<TensorValue<Number>> grad_u_neighbor;

  // Now we will loop over all the elements in the mesh.
  // We will compute the element matrix and right-hand-side
  // contribution.
  //
  // Element iterators are a nice way to iterate through all the
  // elements, or all the elements that have some property.  The
  // iterator el will iterate from the first to the last element on
  // the local processor.  The iterator end_el tells us when to stop.
  // It is smart to make this one const so that we don't accidentally
  // mess it up!  In case users later modify this program to include
  // refinement, we will be safe and will only consider the active
  // elements; hence we use a variant of the active_elem_iterator.
  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices(elem, dof_indices);

    // Compute the element-specific data for the current
    // element.  This involves computing the location of the
    // quadrature points (q_point) and the shape functions
    // (phi, dphi) for the current element.
    fe->reinit(elem);

    libmesh_assert_msg(dphi.size() == dof_indices.size(),
                       "dphi size doesn't match dof_indices size");

    // DenseVector::resize() member will automatically zero out the vector.
    Fe.resize(dof_indices.size());

    // Get the local solution vector
    X.get(dof_indices, dof_u);

    // build the element solution gradient
    grad_u.resize(qrule.n_points());
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      grad_u[qp] = 0;
      for (std::size_t i = 0; i < dof_indices.size(); i++)
        grad_u[qp] += dof_u[i] * dphi[i][qp];
    }

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      auto forcing_func = VectorValue<Number>(forcing_function(0, xyz[qp](0), xyz[qp](1)),
                                              forcing_function(1, xyz[qp](0), xyz[qp](1)),
                                              0);
      for (std::size_t i = 0; i < dof_indices.size(); i++)
      {
        // diffusion elemental residual
        Fe(i) += JxW[qp] * dphi[i][qp].contract(grad_u[qp]);

        // forcing function
        Fe(i) += JxW[qp] * phi[i][qp] * forcing_func;
      }
    }

    // To avoid extraneous allocation when building element sides (to compute their volumes)
    ElemSideBuilder side_builder;

    // Now we consider residual contributions from the sides
    for (auto side : elem->side_index_range())
    {
      // We need to compute h for penalty terms
      const auto side_volume = side_builder(*elem, side).volume();
      fe_face->reinit(elem, side);
      const auto elem_b_order = static_cast<unsigned int>(fe_face->get_order());
      const auto h_elem = elem->volume() / side_volume / std::pow(elem_b_order, 2);

      // build the face solution value and gradient
      u.resize(qface.n_points());
      grad_u.resize(qface.n_points());
      for (unsigned int qp = 0; qp < qface.n_points(); qp++)
      {
        u[qp] = 0;
        grad_u[qp] = 0;
        for (std::size_t i = 0; i < dof_indices.size(); i++)
        {
          u[qp] += dof_u[i] * phi_face[i][qp];
          grad_u[qp] += dof_u[i] * dphi_face[i][qp];
        }
      }

      // No neighbor means we must be on a boundary
      if (!elem->neighbor_ptr(side))
      {
        for (unsigned int qp = 0; qp < qface.n_points(); qp++)
        {
          auto fn = VectorValue<Number>(exact_solution(0, face_xyz[qp](0), face_xyz[qp](1)),
                                        exact_solution(1, face_xyz[qp](0), face_xyz[qp](1)),
                                        0);
          for (std::size_t i = 0; i < dof_indices.size(); i++)
          {
            Fe(i) -= grad_u[qp] * normals[qp] * phi_face[i][qp] * JxW_face[qp];
            Fe(i) += epsilon * (u[qp] - fn) * dphi_face[i][qp] * normals[qp] * JxW_face[qp];
            Fe(i) += sigma / h_elem * (u[qp] - fn) * phi_face[i][qp] * JxW_face[qp];
          }
        }
      }
      else // We must be on an interior side
      {
        const Elem * neighbor = elem->neighbor_ptr(side);

        const auto elem_id = elem->id();
        const auto neighbor_id = neighbor->id();

        // We don't want to erroneously add multiple contributions from the same interior face
        if ((neighbor->active() && (neighbor->level() == elem->level()) &&
             (elem_id < neighbor_id)) ||
            (neighbor->level() < elem->level()))
        {
          dof_map.dof_indices(neighbor, dof_indices_neighbor);

          // Make sure we have the matching quadrature points on face and neighbor
          std::vector<Point> neighbor_xyz;
          FEMap::inverse_map(elem->dim(), neighbor, face_xyz,
                             neighbor_xyz);

          fe_neighbor_face->reinit(neighbor, &neighbor_xyz);

          libmesh_assert_msg(dphi_neighbor.size() == dof_indices_neighbor.size(),
                             "dphi_neighbor size doesn't match dof_indices_neighbor size");

          Fn.resize(dof_indices_neighbor.size());
          X.get(dof_indices_neighbor, dof_u_neighbor);

          // build the neighbor solution value and gradient
          u_neighbor.resize(qface.n_points());
          grad_u_neighbor.resize(qface.n_points());
          for (unsigned int qp = 0; qp < qface.n_points(); qp++)
          {
            u_neighbor[qp] = 0;
            grad_u_neighbor[qp] = 0;
            for (std::size_t i = 0; i < dof_indices_neighbor.size(); i++)
            {
              u_neighbor[qp] += dof_u_neighbor[i] * phi_neighbor[i][qp];
              grad_u_neighbor[qp] += dof_u_neighbor[i] * dphi_neighbor[i][qp];
            }
          }

          // Now add the DG contribution to the local residual
          for (unsigned int qp = 0; qp < qface.n_points(); qp++)
          {
            // element contribution
            for (std::size_t i = 0; i < dof_indices.size(); i++)
            {
              Fe(i) -= 0.5 * (grad_u[qp] * normals[qp] + grad_u_neighbor[qp] * normals[qp]) *
                       phi_face[i][qp] * JxW_face[qp];
              Fe(i) += epsilon * 0.5 * (u[qp] - u_neighbor[qp]) * dphi_face[i][qp] * normals[qp] *
                       JxW_face[qp];
              Fe(i) += sigma / h_elem * (u[qp] - u_neighbor[qp]) * phi_face[i][qp] * JxW_face[qp];
            }
            // Neighbor contribution
            for (std::size_t i = 0; i < dof_indices_neighbor.size(); i++)
            {
              Fn(i) += 0.5 * (grad_u[qp] * normals[qp] + grad_u_neighbor[qp] * normals[qp]) *
                       phi_neighbor[i][qp] * JxW_face[qp];
              Fn(i) -= epsilon * 0.5 * (u[qp] - u_neighbor[qp]) * dphi_neighbor[i][qp] *
                       normals[qp] * JxW_face[qp];
              Fn(i) -=
                  sigma / h_elem * (u[qp] - u_neighbor[qp]) * phi_neighbor[i][qp] * JxW_face[qp];
            }
          }

          R.add_vector(Fn, dof_indices_neighbor);
        } // whether we've done this internal side integral before
      }   // whether we're on an internal side
    }     // loop over sides

    R.add_vector(Fe, dof_indices);
  } // element loop
}

// We actually have no use for X here, which is an indication that this is actually a linear system,
// but it's good to have the non-linear demo
void
compute_jacobian(const NumericVector<Number> &,
                 SparseMatrix<Number> & J,
                 NonlinearImplicitSystem & system)
{
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to(system.name(), "Poisson");

  // Get the DG parameters
  auto & es = system.get_equation_systems();
  const auto sigma = es.parameters.get<Real>("sigma");
  const auto epsilon = es.parameters.get<Real>("epsilon");

  // Get a constant reference to the mesh object.
  const MeshBase & mesh = system.get_mesh();

  // The dimension that we are running
  const auto dim = mesh.mesh_dimension();

  // A reference to the  DofMap object for this system.  The  DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the  DofMap
  // in future examples.
  const DofMap & dof_map = system.get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.
  // Note that FEVectorBase is a typedef for the templated FE
  // class.
  std::unique_ptr<FEVectorBase> fe(FEVectorBase::build(dim, fe_type));

  // An automatically determined Gauss quadrature rule for numerical integration.
  QGauss qrule(dim, fe_type.default_quadrature_order());

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule(&qrule);

  // Declare a special finite element object for face integration.
  std::unique_ptr<FEVectorBase> fe_face(FEVectorBase::build(dim, fe_type));
  // And for neighbor integration
  std::unique_ptr<FEVectorBase> fe_neighbor_face(FEVectorBase::build(dim, fe_type));

  // Boundary integration requires one quadrature rule,
  // with dimensionality one less than the dimensionality
  // of the element.
  QGauss qface(dim - 1, fe_type.default_quadrature_order());

  // Tell the face finite element objects to use our
  // quadrature rule.
  fe_face->attach_quadrature_rule(&qface);
  fe_neighbor_face->attach_quadrature_rule(&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.
  const auto & JxW = fe->get_JxW();

  // The element shape function gradients evaluated at the quadrature points.
  const auto & dphi = fe->get_dphi();

  // face integration points
  const auto & face_xyz = fe_face->get_xyz();

  // Face shape function values
  const auto & phi_face = fe_face->get_phi();

  // Face shape function gradients
  const auto & dphi_face = fe_face->get_dphi();

  // Neighbor shape function values
  const auto & phi_neighbor = fe_neighbor_face->get_phi();

  // Neighbor face shape function gradients
  const auto & dphi_neighbor = fe_neighbor_face->get_dphi();

  // face normals
  const auto & normals = fe_face->get_normals();

  // face JxW
  const auto & JxW_face = fe_face->get_JxW();

  DenseMatrix<Number> Kee;
  DenseMatrix<Number> Ken;
  DenseMatrix<Number> Kne;
  DenseMatrix<Number> Knn;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_neighbor;

  // To avoid extraneous allocation when building element sides (to compute their volumes)
  ElemSideBuilder side_builder;

  // Now we will loop over all the elements in the mesh.
  // We will compute the element matrix and right-hand-side
  // contribution.
  //
  // Element iterators are a nice way to iterate through all the
  // elements, or all the elements that have some property.  The
  // iterator el will iterate from the first to the last element on
  // the local processor.  The iterator end_el tells us when to stop.
  // It is smart to make this one const so that we don't accidentally
  // mess it up!  In case users later modify this program to include
  // refinement, we will be safe and will only consider the active
  // elements; hence we use a variant of the active_elem_iterator.
  for (const auto & elem : mesh.active_local_element_ptr_range())
  {
    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices(elem, dof_indices);

    // Compute the element-specific data for the current
    // element.  This involves computing the location of the
    // quadrature points (q_point) and the shape functions
    // (phi, dphi) for the current element.
    fe->reinit(elem);

    libmesh_assert_msg(dphi.size() == dof_indices.size(),
                       "dphi size doesn't match dof_indices size");

    // DenseMatrix::resize() member will automatically zero out the matrix.
    Kee.resize(dof_indices.size(), dof_indices.size());

    // diffusion elemental jacobian
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
      for (std::size_t i = 0; i < dof_indices.size(); i++)
        for (std::size_t j = 0; j < dof_indices.size(); j++)
          Kee(i, j) += JxW[qp] * dphi[i][qp].contract(dphi[j][qp]);

    // Now we consider jacobian contributions from the sides
    for (auto side : elem->side_index_range())
    {
      // We need to compute h for penalty terms
      const auto side_volume = side_builder(*elem, side).volume();
      fe_face->reinit(elem, side);
      const auto elem_b_order = static_cast<unsigned int>(fe_face->get_order());
      const auto h_elem = elem->volume() / side_volume / std::pow(elem_b_order, 2);

      // No neighbor means we must be on a boundary
      if (!elem->neighbor_ptr(side))
      {
        for (unsigned int qp = 0; qp < qface.n_points(); qp++)
          for (std::size_t i = 0; i < dof_indices.size(); i++)
            for (std::size_t j = 0; j < dof_indices.size(); j++)

            {
              Kee(i, j) -= dphi_face[j][qp] * normals[qp] * phi_face[i][qp] * JxW_face[qp];
              Kee(i, j) +=
                  epsilon * phi_face[j][qp] * dphi_face[i][qp] * normals[qp] * JxW_face[qp];
              Kee(i, j) += sigma / h_elem * phi_face[j][qp] * phi_face[i][qp] * JxW_face[qp];
            }
      }
      else // We must be on an interior side
      {
        const Elem * neighbor = elem->neighbor_ptr(side);

        const auto elem_id = elem->id();
        const auto neighbor_id = neighbor->id();

        // We don't want to erroneously add multiple contributions from the same interior face
        if ((neighbor->active() && (neighbor->level() == elem->level()) &&
             (elem_id < neighbor_id)) ||
            (neighbor->level() < elem->level()))
        {
          dof_map.dof_indices(neighbor, dof_indices_neighbor);

          // Make sure we have the matching quadrature points on face and neighbor
          std::vector<Point> neighbor_xyz;
          FEMap::inverse_map(elem->dim(), neighbor, face_xyz,
                             neighbor_xyz);

          fe_neighbor_face->reinit(neighbor, &neighbor_xyz);

          libmesh_assert_msg(dphi_neighbor.size() == dof_indices_neighbor.size(),
                             "dphi_neighbor size doesn't match dof_indices_neighbor size");

          Ken.resize(dof_indices.size(), dof_indices_neighbor.size());
          Kne.resize(dof_indices_neighbor.size(), dof_indices.size());
          Knn.resize(dof_indices_neighbor.size(), dof_indices_neighbor.size());

          // Now add the DG contribution to the local jacobian
          for (unsigned int qp = 0; qp < qface.n_points(); qp++)
          {
            // element-element contribution
            for (std::size_t i = 0; i < dof_indices.size(); i++)
              for (std::size_t j = 0; j < dof_indices.size(); j++)
              {
                Kee(i, j) -= 0.5 * dphi_face[j][qp] * normals[qp] * phi_face[i][qp] * JxW_face[qp];
                Kee(i, j) +=
                    epsilon * 0.5 * phi_face[j][qp] * dphi_face[i][qp] * normals[qp] * JxW_face[qp];
                Kee(i, j) += sigma / h_elem * phi_face[j][qp] * phi_face[i][qp] * JxW_face[qp];
              }
            // element-neighbor contribution
            for (std::size_t i = 0; i < dof_indices.size(); i++)
              for (std::size_t j = 0; j < dof_indices_neighbor.size(); j++)
              {
                Ken(i, j) -=
                    0.5 * dphi_neighbor[j][qp] * normals[qp] * phi_face[i][qp] * JxW_face[qp];
                Ken(i, j) += epsilon * 0.5 * -phi_neighbor[j][qp] * dphi_face[i][qp] * normals[qp] *
                             JxW_face[qp];
                Ken(i, j) += sigma / h_elem * -phi_neighbor[j][qp] * phi_face[i][qp] * JxW_face[qp];
              }
            // Neighbor-element contribution
            for (std::size_t i = 0; i < dof_indices_neighbor.size(); i++)
              for (std::size_t j = 0; j < dof_indices_neighbor.size(); j++)
              {
                Kne(i, j) +=
                    0.5 * dphi_face[j][qp] * normals[qp] * phi_neighbor[i][qp] * JxW_face[qp];
                Kne(i, j) -= epsilon * 0.5 * phi_face[j][qp] * dphi_neighbor[i][qp] * normals[qp] *
                             JxW_face[qp];
                Kne(i, j) -= sigma / h_elem * phi_face[j][qp] * phi_neighbor[i][qp] * JxW_face[qp];
              }
            // Neighbor-neighbor contribution
            for (std::size_t i = 0; i < dof_indices_neighbor.size(); i++)
              for (std::size_t j = 0; j < dof_indices_neighbor.size(); j++)
              {
                Knn(i, j) +=
                    0.5 * dphi_neighbor[j][qp] * normals[qp] * phi_neighbor[i][qp] * JxW_face[qp];
                Knn(i, j) -= epsilon * 0.5 * -phi_neighbor[j][qp] * dphi_neighbor[i][qp] *
                             normals[qp] * JxW_face[qp];
                Knn(i, j) -=
                    sigma / h_elem * -phi_neighbor[j][qp] * phi_neighbor[i][qp] * JxW_face[qp];
              }
          }

          J.add_matrix(Ken, dof_indices, dof_indices_neighbor);
          J.add_matrix(Kne, dof_indices_neighbor, dof_indices);
          J.add_matrix(Knn, dof_indices_neighbor, dof_indices_neighbor);
        } // whether we've done this internal side integral before
      }   // whether we're on an internal side
    }     // loop over sides

    J.add_matrix(Kee, dof_indices, dof_indices);
  } // element loop
}
