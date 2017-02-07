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
#include "libmesh/edge_edge2.h"
#include LIBMESH_INCLUDE_UNORDERED_SET

// The nonlinear solver and system we will be using
#include "libmesh/nonlinear_solver.h"
#include "libmesh/nonlinear_implicit_system.h"

using namespace libMesh;

LinearElasticityWithContact::LinearElasticityWithContact (NonlinearImplicitSystem & sys_in,
                                                          Real contact_penalty_in) :
  _sys(sys_in),
  _contact_penalty(contact_penalty_in)
{
}

void LinearElasticityWithContact::set_contact_penalty(Real contact_penalty_in)
{
  _contact_penalty = contact_penalty_in;
}

Real LinearElasticityWithContact::get_contact_penalty() const
{
  return _contact_penalty;
}

Real LinearElasticityWithContact::kronecker_delta(unsigned int i,
                                                  unsigned int j)
{
  return i == j ? 1. : 0.;
}

Real LinearElasticityWithContact::elasticity_tensor(Real young_modulus,
                                                    Real poisson_ratio,
                                                    unsigned int i,
                                                    unsigned int j,
                                                    unsigned int k,
                                                    unsigned int l)
{
  // Define the Lame constants
  const Real lambda_1 = (young_modulus*poisson_ratio)/((1.+poisson_ratio)*(1.-2.*poisson_ratio));
  const Real lambda_2 = young_modulus/(2.*(1.+poisson_ratio));

  return lambda_1 * kronecker_delta(i, j) * kronecker_delta(k, l) +
    lambda_2 * (kronecker_delta(i, k) * kronecker_delta(j, l) + kronecker_delta(i, l) * kronecker_delta(j, k));
}

void LinearElasticityWithContact::move_mesh (MeshBase & input_mesh,
                                             const NumericVector<Number> & input_solution)
{
  // Maintain a set of node ids that we've encountered.
  LIBMESH_BEST_UNORDERED_SET<dof_id_type> encountered_node_ids;

  // Localize input_solution so that we have the data to move all
  // elements (not just elements local to this processor).
  UniquePtr< NumericVector<Number> > localized_input_solution =
    NumericVector<Number>::build(input_solution.comm());

  localized_input_solution->init (input_solution.size(), false, SERIAL);
  input_solution.localize(*localized_input_solution);

  MeshBase::const_element_iterator       el     = input_mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = input_mesh.active_elements_end();

  for ( ; el != end_el; ++el)
    {
      Elem * elem = *el;
      Elem * orig_elem = _sys.get_mesh().elem_ptr(elem->id());

      for (unsigned int node_id=0; node_id<elem->n_nodes(); node_id++)
        {
          Node & node = elem->node_ref(node_id);

          if (encountered_node_ids.find(node.id()) != encountered_node_ids.end())
            continue;

          encountered_node_ids.insert(node.id());

          std::vector<std::string> uvw_names(3);
          uvw_names[0] = "u";
          uvw_names[1] = "v";
          uvw_names[2] = "w";

          {
            const Point master_point = elem->master_point(node_id);

            Point uvw;
            for (std::size_t index=0; index<uvw_names.size(); index++)
              {
                const unsigned int var = _sys.variable_number(uvw_names[index]);
                const FEType & fe_type = _sys.get_dof_map().variable_type(var);

                FEComputeData data (_sys.get_equation_systems(), master_point);

                FEInterface::compute_data(elem->dim(),
                                          fe_type,
                                          elem,
                                          data);

                std::vector<dof_id_type> dof_indices_var;
                _sys.get_dof_map().dof_indices (orig_elem, dof_indices_var, var);

                for (std::size_t i=0; i<dof_indices_var.size(); i++)
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
            node += uvw;
          }
        }
    }
}

void LinearElasticityWithContact::initialize_contact_load_paths()
{
  const MeshBase & mesh = _sys.get_mesh();

  std::vector<dof_id_type> nodes_on_lower_surface;
  std::vector<dof_id_type> nodes_on_upper_surface;

  MeshBase::const_element_iterator       el     = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  _lambdas.clear();
  for ( ; el != end_el; ++el)
    {
      const Elem * elem = *el;

      for (unsigned int side=0; side<elem->n_sides(); side++)
        {
          if (elem->neighbor_ptr(side) == libmesh_nullptr)
            {
              bool on_lower_contact_surface =
                mesh.get_boundary_info().has_boundary_id (elem, side, CONTACT_BOUNDARY_LOWER);

              bool on_upper_contact_surface =
                mesh.get_boundary_info().has_boundary_id (elem, side, CONTACT_BOUNDARY_UPPER);

              if (on_lower_contact_surface && on_upper_contact_surface)
                libmesh_error_msg("Should not be on both surfaces at the same time");

              if (on_lower_contact_surface || on_upper_contact_surface)
                {
                  for (unsigned int node_index=0; node_index<elem->n_nodes(); node_index++)
                    if (elem->is_node_on_side(node_index, side))
                      {
                        if (on_lower_contact_surface)
                          nodes_on_lower_surface.push_back(elem->node_id(node_index));
                        else
                          {
                            _lambdas[elem->node_id(node_index)] = 0.;
                            nodes_on_upper_surface.push_back(elem->node_id(node_index));
                          }
                      }
                }

            } // end if nieghbor(side_) != libmesh_nullptr
        } // end for side
    } // end for el

  // In this example, we expect the number of upper and lower nodes to match
  libmesh_assert(nodes_on_lower_surface.size() == nodes_on_upper_surface.size());

  // Do an N^2 search to match the contact nodes
  _contact_node_map.clear();
  for (std::size_t i=0; i<nodes_on_lower_surface.size(); i++)
    {
      dof_id_type lower_node_id = nodes_on_lower_surface[i];
      Point p_lower = mesh.point(lower_node_id);

      Real min_distance = std::numeric_limits<Real>::max();

      for (std::size_t j=0; j<nodes_on_upper_surface.size(); j++)
        {
          dof_id_type upper_node_id = nodes_on_upper_surface[j];
          Point p_upper = mesh.point(upper_node_id);

          Real distance = (p_upper - p_lower).norm();

          if (distance < min_distance)
            {
              _contact_node_map[lower_node_id] = upper_node_id;
              min_distance = distance;
            }
        }
    }
}

void LinearElasticityWithContact::add_contact_edge_elements()
{
  MeshBase & mesh = _sys.get_mesh();

  std::map<dof_id_type, dof_id_type>::iterator it = _contact_node_map.begin();
  std::map<dof_id_type, dof_id_type>::iterator it_end = _contact_node_map.end();
  for( ; it != it_end ; ++it )
    {
      dof_id_type master_node_id = it->first;
      dof_id_type slave_node_id = it->second;

      Node& master_node = mesh.node(master_node_id);
      Node& slave_node = mesh.node(slave_node_id);

      Elem* connector_elem = mesh.add_elem (new Edge2);
      connector_elem->set_node(0) = &master_node;
      connector_elem->set_node(1) = &slave_node;

      connector_elem->subdomain_id() = 10;
    }

  mesh.prepare_for_use();
}

void LinearElasticityWithContact::residual_and_jacobian (const NumericVector<Number> & soln,
                                                         NumericVector<Number> * residual,
                                                         SparseMatrix<Number> * jacobian,
                                                         NonlinearImplicitSystem & /*sys*/)
{
  EquationSystems & es = _sys.get_equation_systems();
  const Real young_modulus = es.parameters.get<Real>("young_modulus");
  const Real poisson_ratio = es.parameters.get<Real>("poisson_ratio");

  const MeshBase & mesh = _sys.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  const unsigned int u_var = _sys.variable_number ("u");

  DofMap & dof_map = _sys.get_dof_map();

  FEType fe_type = dof_map.variable_type(u_var);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_type));
  QGauss qface (dim-1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule (&qface);

  UniquePtr<FEBase> fe_neighbor_face (FEBase::build(dim, fe_type));
  fe_neighbor_face->attach_quadrature_rule (&qface);

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  if (jacobian)
    jacobian->zero();

  if (residual)
    residual->zero();

  // Do jacobian and residual assembly, including contact forces
  DenseVector<Number> Re;

  DenseSubVector<Number> Re_var[3] =
    {DenseSubVector<Number>(Re),
     DenseSubVector<Number>(Re),
     DenseSubVector<Number>(Re)};

  DenseMatrix<Number> Ke;
  DenseSubMatrix<Number> Ke_var[3][3] =
    {
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}
    };

  std::vector<dof_id_type> dof_indices;
  std::vector< std::vector<dof_id_type> > dof_indices_var(3);

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem * elem = *el;

      if( elem->type() == EDGE2 )
        {
          // We do not do any assembly on the contact connector elements.
          // The contact force assembly is handled in a separate loop.
          continue;
        }

      dof_map.dof_indices (elem, dof_indices);
      for (unsigned int var=0; var<3; var++)
        dof_map.dof_indices (elem, dof_indices_var[var], var);

      const unsigned int n_dofs   = dof_indices.size();
      const unsigned int n_var_dofs = dof_indices_var[0].size();

      fe->reinit (elem);

      Re.resize (n_dofs);
      for (unsigned int var=0; var<3; var++)
        Re_var[var].reposition (var*n_var_dofs, n_var_dofs);

      Ke.resize (n_dofs, n_dofs);
      for (unsigned int var_i=0; var_i<3; var_i++)
        for (unsigned int var_j=0; var_j<3; var_j++)
          Ke_var[var_i][var_j].reposition (var_i*n_var_dofs, var_j*n_var_dofs, n_var_dofs, n_var_dofs);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Row is variable u, v, or w column is x, y, or z
          DenseMatrix<Number> grad_u(3, 3);
          for (unsigned int var_i=0; var_i<3; var_i++)
            for (unsigned int var_j=0; var_j<3; var_j++)
              for (unsigned int j=0; j<n_var_dofs; j++)
                grad_u(var_i, var_j) += dphi[j][qp](var_j)*soln(dof_indices_var[var_i][j]);

          // - C_ijkl u_k,l v_i,j
          for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
            for (unsigned int i=0; i<3; i++)
              for (unsigned int j=0; j<3; j++)
                for (unsigned int k=0; k<3; k++)
                  for (unsigned int l=0; l<3; l++)
                    Re_var[i](dof_i) -= JxW[qp] *
                      elasticity_tensor(young_modulus, poisson_ratio, i, j, k, l) * grad_u(k,l) * dphi[dof_i][qp](j);

          // assemble \int_Omega C_ijkl u_k,l v_i,j \dx
          for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
            for (unsigned int dof_j=0; dof_j<n_var_dofs; dof_j++)
              for (unsigned int i=0; i<3; i++)
                for (unsigned int j=0; j<3; j++)
                  for (unsigned int k=0; k<3; k++)
                    for (unsigned int l=0; l<3; l++)
                      Ke_var[i][k](dof_i, dof_j) -= JxW[qp] *
                        elasticity_tensor(young_modulus, poisson_ratio, i, j, k, l) * dphi[dof_j][qp](l) * dphi[dof_i][qp](j);
        }

      dof_map.constrain_element_matrix_and_vector (Ke, Re, dof_indices);

      if (jacobian)
        jacobian->add_matrix (Ke, dof_indices);

      if (residual)
        residual->add_vector (Re, dof_indices);
    }

  // Move a copy of the mesh based on the solution. This is used to get
  // the contact forces. We could compute the contact forces based on the
  // current displacement solution, which would not require a clone of the
  // mesh. Avoiding the mesh clone would be important for production-scale
  // contact solves, but for the sake of this example, using the clone is
  // simple and fast enough.
  UniquePtr<MeshBase> mesh_clone = mesh.clone();
  move_mesh(*mesh_clone, soln);

  // Add contributions due to contact penalty forces. Only need to do this on
  // one processor.
  _lambda_plus_penalty_values.clear();

  std::map<dof_id_type, dof_id_type>::iterator it = _contact_node_map.begin();
  std::map<dof_id_type, dof_id_type>::iterator it_end = _contact_node_map.end();

  for ( ; it != it_end; ++it)
    {
      dof_id_type lower_point_id = it->first;
      dof_id_type upper_point_id = it->second;

      Point upper_to_lower;
      {
        Point lower_node_moved = mesh_clone->point(lower_point_id);
        Point upper_node_moved = mesh_clone->point(upper_point_id);

        upper_to_lower = (lower_node_moved - upper_node_moved);
      }

      // Set the contact force direction. Usually this would be calculated
      // separately on each master node, but here we just set it to (0, 0, 1)
      // everywhere.
      Point contact_force_direction(0., 0., 1.);

      // gap_function > 0. means that contact has been detected
      // gap_function < 0. means that we have a gap
      // This sign convention matches Simo & Laursen (1992).
      Real gap_function = upper_to_lower * contact_force_direction;

      // We use the sign of lambda_plus_penalty to determine whether or
      // not we need to impose a contact force.
      Real lambda_plus_penalty =
        (_lambdas[upper_point_id] + gap_function * _contact_penalty);

      if (lambda_plus_penalty < 0.)
        lambda_plus_penalty = 0.;

      // Store lambda_plus_penalty, we'll need to use it later to update _lambdas
      _lambda_plus_penalty_values[upper_point_id] = lambda_plus_penalty;

      const Node & lower_node = mesh.node_ref(lower_point_id);
      const Node & upper_node = mesh.node_ref(upper_point_id);

      std::vector<dof_id_type> dof_indices_on_lower_node(3);
      std::vector<dof_id_type> dof_indices_on_upper_node(3);
      DenseVector<Number> lower_contact_force_vec(3);
      DenseVector<Number> upper_contact_force_vec(3);

      for (unsigned int var=0; var<3; var++)
        {
          dof_indices_on_lower_node[var] = lower_node.dof_number(_sys.number(), var, 0);
          lower_contact_force_vec(var) = -lambda_plus_penalty * contact_force_direction(var);

          dof_indices_on_upper_node[var] = upper_node.dof_number(_sys.number(), var, 0);
          upper_contact_force_vec(var) = lambda_plus_penalty * contact_force_direction(var);
        }

      if (lambda_plus_penalty > 0.)
        {
          if (residual && (_sys.comm().rank() == 0))
            {
              residual->add_vector (lower_contact_force_vec, dof_indices_on_lower_node);
              residual->add_vector (upper_contact_force_vec, dof_indices_on_upper_node);
            }

          // Add the Jacobian terms due to the contact forces. The lambda contribution
          // is not relevant here because it doesn't depend on the solution.
          if (jacobian && (_sys.comm().rank() == 0))
            for (unsigned int var=0; var<3; var++)
              for (unsigned int j=0; j<3; j++)
                {
                  jacobian->add(dof_indices_on_lower_node[var],
                                dof_indices_on_upper_node[j],
                                _contact_penalty * contact_force_direction(j) * contact_force_direction(var));

                  jacobian->add(dof_indices_on_lower_node[var],
                                dof_indices_on_lower_node[j],
                                -_contact_penalty * contact_force_direction(j) * contact_force_direction(var));

                  jacobian->add(dof_indices_on_upper_node[var],
                                dof_indices_on_lower_node[j],
                                _contact_penalty * contact_force_direction(j) * contact_force_direction(var));

                  jacobian->add(dof_indices_on_upper_node[var],
                                dof_indices_on_upper_node[j],
                                -_contact_penalty * contact_force_direction(j) * contact_force_direction(var));
                }
        }
      else
        {
          // We add zeros to the matrix even when lambda_plus_penalty = 0.
          // We do this because some linear algebra libraries (e.g. PETSc)
          // will condense out any unused entries from the sparsity pattern,
          // so adding these zeros in ensures that these entries are not
          // condensed out.
          if (jacobian && (_sys.comm().rank() == 0))
            for (unsigned int var=0; var<3; var++)
              for (unsigned int j=0; j<3; j++)
                {
                  jacobian->add(dof_indices_on_lower_node[var],
                                dof_indices_on_upper_node[j],
                                0.);

                  jacobian->add(dof_indices_on_lower_node[var],
                                dof_indices_on_lower_node[j],
                                0.);

                  jacobian->add(dof_indices_on_upper_node[var],
                                dof_indices_on_lower_node[j],
                                0.);

                  jacobian->add(dof_indices_on_upper_node[var],
                                dof_indices_on_upper_node[j],
                                0.);
                }
        }
    }
}

void LinearElasticityWithContact::compute_stresses()
{
  EquationSystems & es = _sys.get_equation_systems();
  const Real young_modulus = es.parameters.get<Real>("young_modulus");
  const Real poisson_ratio = es.parameters.get<Real>("poisson_ratio");

  const MeshBase & mesh = _sys.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  unsigned int displacement_vars[3];
  displacement_vars[0] = _sys.variable_number ("u");
  displacement_vars[1] = _sys.variable_number ("v");
  displacement_vars[2] = _sys.variable_number ("w");
  const unsigned int u_var = _sys.variable_number ("u");

  const DofMap & dof_map = _sys.get_dof_map();
  FEType fe_type = dof_map.variable_type(u_var);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  // Also, get a reference to the ExplicitSystem
  ExplicitSystem & stress_system = es.get_system<ExplicitSystem>("StressSystem");
  const DofMap & stress_dof_map = stress_system.get_dof_map();
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
  DenseMatrix<Number> elem_avg_stress_tensor(3, 3);

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem * elem = *el;

      if( elem->type() == EDGE2 )
        {
          // We do not compute stress on the contact connector elements.
          continue;
        }

      for (unsigned int var=0; var<3; var++)
        dof_map.dof_indices (elem, dof_indices_var[var], displacement_vars[var]);

      const unsigned int n_var_dofs = dof_indices_var[0].size();

      fe->reinit (elem);

      // clear the stress tensor
      elem_avg_stress_tensor.resize(3, 3);

      for (unsigned int qp=0; qp<qrule.n_points(); qp++)
        {
          // Row is variable u1, u2, or u3, column is x, y, or z
          DenseMatrix<Number> grad_u(3, 3);
          for (unsigned int var_i=0; var_i<3; var_i++)
            for (unsigned int var_j=0; var_j<3; var_j++)
              for (unsigned int j=0; j<n_var_dofs; j++)
                grad_u(var_i, var_j) +=
                  dphi[j][qp](var_j) * _sys.current_solution(dof_indices_var[var_i][j]);

          DenseMatrix<Number> stress_tensor(3, 3);
          for (unsigned int i=0; i<3; i++)
            for (unsigned int j=0; j<3; j++)
              for (unsigned int k=0; k<3; k++)
                for (unsigned int l=0; l<3; l++)
                  stress_tensor(i,j) +=
                    elasticity_tensor(young_modulus, poisson_ratio, i, j, k, l) * grad_u(k,l);

          // We want to plot the average stress on each element, hence
          // we integrate stress_tensor
          elem_avg_stress_tensor.add(JxW[qp], stress_tensor);
        }

      // Get the average stress per element by dividing by volume
      elem_avg_stress_tensor.scale(1./elem->volume());

      // load elem_sigma data into stress_system
      unsigned int stress_var_index = 0;
      for (unsigned int i=0; i<3; i++)
        for (unsigned int j=i; j<3; j++)
          {
            stress_dof_map.dof_indices (elem, stress_dof_indices_var, sigma_vars[stress_var_index]);

            // We are using CONSTANT MONOMIAL basis functions, hence we only need to get
            // one dof index per variable
            dof_id_type dof_index = stress_dof_indices_var[0];

            if ((stress_system.solution->first_local_index() <= dof_index) &&
                (dof_index < stress_system.solution->last_local_index()))
              stress_system.solution->set(dof_index, elem_avg_stress_tensor(i,j));

            stress_var_index++;
          }

      // Also, the von Mises stress
      Number vonMises_value = std::sqrt(0.5*(pow(elem_avg_stress_tensor(0,0) - elem_avg_stress_tensor(1,1), 2.) +
                                             pow(elem_avg_stress_tensor(1,1) - elem_avg_stress_tensor(2,2), 2.) +
                                             pow(elem_avg_stress_tensor(2,2) - elem_avg_stress_tensor(0,0), 2.) +
                                             6.*(pow(elem_avg_stress_tensor(0,1), 2.) +
                                                 pow(elem_avg_stress_tensor(1,2), 2.) +
                                                 pow(elem_avg_stress_tensor(2,0), 2.))));

      stress_dof_map.dof_indices (elem, stress_dof_indices_var, vonMises_var);
      dof_id_type dof_index = stress_dof_indices_var[0];

      if ((stress_system.solution->first_local_index() <= dof_index) &&
          (dof_index < stress_system.solution->last_local_index()))
        stress_system.solution->set(dof_index, vonMises_value);
    }

  // Should call close and update when we set vector entries directly
  stress_system.solution->close();
  stress_system.update();
}

std::pair<Real, Real> LinearElasticityWithContact::update_lambdas()
{
  Real max_delta_lambda = 0.;
  Real max_new_lambda = 0.;

  std::map<dof_id_type, Real>::iterator it = _lambdas.begin();
  std::map<dof_id_type, Real>::iterator it_end = _lambdas.end();
  for ( ; it != it_end; ++it)
    {
      dof_id_type upper_node_id = it->first;

      std::map<dof_id_type, Real>::iterator new_lambda_it = _lambda_plus_penalty_values.find(upper_node_id);
      if (new_lambda_it == _lambda_plus_penalty_values.end())
        libmesh_error_msg("New lambda value not found");

      Real new_lambda = new_lambda_it->second;
      Real old_lambda = it->second;

      it->second = new_lambda;

      Real delta_lambda = std::abs(new_lambda-old_lambda);
      if (delta_lambda > max_delta_lambda)
        max_delta_lambda = delta_lambda;

      if (std::abs(new_lambda) > max_new_lambda)
        max_new_lambda = std::abs(new_lambda);
    }

  return std::make_pair(max_delta_lambda, max_new_lambda);
}

std::pair<Real, Real> LinearElasticityWithContact::get_least_and_max_gap_function()
{
  UniquePtr<MeshBase> mesh_clone = _sys.get_mesh().clone();
  move_mesh(*mesh_clone, *_sys.solution);

  Real least_value = std::numeric_limits<Real>::max();
  Real max_value = std::numeric_limits<Real>::min();

  std::map<dof_id_type, dof_id_type>::iterator it = _contact_node_map.begin();
  std::map<dof_id_type, dof_id_type>::iterator it_end = _contact_node_map.end();
  for ( ; it != it_end; ++it)
    {
      dof_id_type lower_point_id = it->first;
      dof_id_type upper_point_id = it->second;

      Point upper_to_lower;
      {
        Point lower_node_moved = mesh_clone->point(lower_point_id);
        Point upper_node_moved = mesh_clone->point(upper_point_id);

        upper_to_lower = (lower_node_moved - upper_node_moved);
      }

      // Set the contact force direction. Usually this would be calculated
      // separately on each master node, but here we just set it to (0, 0, 1)
      // everywhere.
      Point contact_force_direction(0., 0., 1.);

      // gap_function > 0. means that contact has been detected
      // gap_function < 0. means that we have a gap
      // This sign convention matches Simo & Laursen (1992).
      Real gap_function = upper_to_lower * contact_force_direction;

      if (gap_function < least_value)
        least_value = gap_function;

      if (gap_function > max_value)
        max_value = gap_function;
    }

  return std::make_pair(least_value, max_value);
}
