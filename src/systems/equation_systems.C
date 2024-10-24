// The libMesh Finite Element Library.
// Copyright (C) 2002-2024 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


// Local Includes
#include "libmesh/default_coupling.h" // For downconversion
#include "libmesh/dof_map.h"
#include "libmesh/eigen_system.h"
#include "libmesh/elem.h"
#include "libmesh/explicit_system.h"
#include "libmesh/fe_interface.h"
#include "libmesh/frequency_system.h"
#include "libmesh/int_range.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/newmark_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/parallel.h"
#include "libmesh/rb_construction.h"
#include "libmesh/remote_elem.h"
#include "libmesh/transient_rb_construction.h"
#include "libmesh/transient_system.h"

// System includes
#include <functional> // std::plus
#include <numeric> // std::iota
#include <sstream>

// Include the systems before this one to avoid
// overlapping forward declarations.
#include "libmesh/equation_systems.h"

namespace libMesh
{

EquationSystems::EquationSystems (MeshBase & m) :
  ParallelObject (m),
  _mesh          (m),
  _refine_in_reinit(true),
  _enable_default_ghosting(true)
{
  // Set default parameters
  this->parameters.set<Real>        ("linear solver tolerance") = TOLERANCE * TOLERANCE;
  this->parameters.set<unsigned int>("linear solver maximum iterations") = 5000;
}



EquationSystems::~EquationSystems () = default;



void EquationSystems::clear ()
{
  // Clear any additional parameters
  parameters.clear ();

  // Clear the systems.
  _systems.clear();
}



void EquationSystems::init ()
{
  libmesh_assert(std::all_of(make_range(this->n_systems()).begin(),
                             make_range(this->n_systems()).end(),
                             [&](unsigned i) { return !this->get_system(i).is_initialized(); }));

  this->reinit_mesh();
}



void EquationSystems::reinit ()
{
  const bool mesh_changed = this->reinit_solutions();

  // If the mesh has changed, systems will need to reinitialize their
  // own data on the new mesh.
  if (mesh_changed)
    this->reinit_systems();
}

void EquationSystems::reinit_mesh ()
{
  const unsigned int n_sys = this->n_systems();

  libmesh_assert_not_equal_to (n_sys, 0);

  // Tell all the \p DofObject entities how many systems
  // there are.
  for (auto & node : _mesh.node_ptr_range())
    node->set_n_systems(n_sys);

  for (auto & elem : _mesh.element_ptr_range())
    elem->set_n_systems(n_sys);

  //for (auto i : make_range(this->n_systems()))
    //this->get_system(i).init();

#ifdef LIBMESH_ENABLE_AMR
  MeshRefinement mesh_refine(_mesh);
  mesh_refine.clean_refinement_flags();
#endif

 // Now loop over all the systems belonging to this ES
 // and call reinit_mesh for each system
 for (auto i : make_range(this->n_systems()))
    this->get_system(i).reinit_mesh();

}

bool EquationSystems::reinit_solutions ()
{
  parallel_object_only();

  const unsigned int n_sys = this->n_systems();
  libmesh_assert_not_equal_to (n_sys, 0);

  // And any new systems will need initialization
  for (unsigned int i=0; i != n_sys; ++i)
    if (!this->get_system(i).is_initialized())
      this->get_system(i).init();

  // We used to assert that all nodes and elements *already* had
  // n_systems() properly set; however this is false in the case where
  // user code has manually added nodes and/or elements to an
  // already-initialized system.

  // Make sure all the \p DofObject entities know how many systems
  // there are.
  {
    // All the nodes
    for (auto & node : _mesh.node_ptr_range())
      node->set_n_systems(n_sys);

    // All the elements
    for (auto & elem : _mesh.element_ptr_range())
      elem->set_n_systems(n_sys);
  }

  // Localize each system's vectors
  for (unsigned int i=0; i != n_sys; ++i)
    this->get_system(i).re_update();

#ifdef LIBMESH_ENABLE_AMR

  bool mesh_changed = false;

  // FIXME: For backwards compatibility, assume
  // refine_and_coarsen_elements or refine_uniformly have already
  // been called
  {
    for (unsigned int i=0; i != n_sys; ++i)
      {
        System & sys = this->get_system(i);

        // Even if the system doesn't have any variables in it we want
        // consistent behavior; e.g. distribute_dofs should have the
        // opportunity to count up zero dofs on each processor.
        //
        // Who's been adding zero-var systems anyway, outside of my
        // unit tests? - RHS
        // if (!sys.n_vars())
        // continue;

        sys.get_dof_map().distribute_dofs(_mesh);

        // Recreate any user or internal constraints
        sys.reinit_constraints();

        sys.get_dof_map().prepare_send_list();

        sys.prolong_vectors();
      }
    mesh_changed = true;
  }

  if (this->_refine_in_reinit)
    {
      // Don't override any user refinement settings
      MeshRefinement mesh_refine(_mesh);
      mesh_refine.face_level_mismatch_limit() = 0; // unlimited
      mesh_refine.overrefined_boundary_limit() = -1; // unlimited
      mesh_refine.underrefined_boundary_limit() = -1; // unlimited

      // Try to coarsen the mesh, then restrict each system's vectors
      // if necessary
      if (mesh_refine.coarsen_elements())
        {
          for (auto i : make_range(this->n_systems()))
            {
              System & sys = this->get_system(i);
              sys.get_dof_map().distribute_dofs(_mesh);
              sys.reinit_constraints();
              sys.get_dof_map().prepare_send_list();
              sys.restrict_vectors();
            }
          mesh_changed = true;
        }

      // Once vectors are all restricted, we can delete
      // children of coarsened elements
      if (mesh_changed)
        this->get_mesh().contract();

      // Try to refine the mesh, then prolong each system's vectors
      // if necessary
      if (mesh_refine.refine_elements())
        {
          for (auto i : make_range(this->n_systems()))
            {
              System & sys = this->get_system(i);
              sys.get_dof_map().distribute_dofs(_mesh);
              sys.reinit_constraints();
              sys.get_dof_map().prepare_send_list();
              sys.prolong_vectors();
            }
          mesh_changed = true;
        }
    }

  return mesh_changed;

#endif // #ifdef LIBMESH_ENABLE_AMR

  return false;
}



void EquationSystems::reinit_systems()
{
  for (auto i : make_range(this->n_systems()))
    this->get_system(i).reinit();
}



void EquationSystems::allgather ()
{
  // A serial mesh means nothing needs to be done
  if (_mesh.is_serial())
    return;

  const unsigned int n_sys = this->n_systems();

  libmesh_assert_not_equal_to (n_sys, 0);

  // Gather the mesh
  _mesh.allgather();

  // Tell all the \p DofObject entities how many systems
  // there are.
  for (auto & node : _mesh.node_ptr_range())
    node->set_n_systems(n_sys);

  for (auto & elem : _mesh.element_ptr_range())
    elem->set_n_systems(n_sys);

  // And distribute each system's dofs
  for (auto i : make_range(this->n_systems()))
    {
      System & sys = this->get_system(i);
      DofMap & dof_map = sys.get_dof_map();
      dof_map.distribute_dofs(_mesh);

      // The user probably won't need constraint equations or the
      // send_list after an allgather, but let's keep it in consistent
      // shape just in case.
      sys.reinit_constraints();
      dof_map.prepare_send_list();
    }
}



void EquationSystems::enable_default_ghosting (bool enable)
{
  _enable_default_ghosting = enable;
  MeshBase &mesh = this->get_mesh();

  if (enable)
    mesh.add_ghosting_functor(mesh.default_ghosting());
  else
    mesh.remove_ghosting_functor(mesh.default_ghosting());

  for (auto i : make_range(this->n_systems()))
    {
      DofMap & dof_map = this->get_system(i).get_dof_map();
      if (enable)
        dof_map.add_default_ghosting();
      else
        dof_map.remove_default_ghosting();
    }
}



void EquationSystems::update ()
{
  LOG_SCOPE("update()", "EquationSystems");

  // Localize each system's vectors
  for (auto i : make_range(this->n_systems()))
    this->get_system(i).update();
}



System & EquationSystems::add_system (std::string_view sys_type,
                                      std::string_view name)
{
  // If the user already built a system with this name, we'll
  // trust them and we'll use it.  That way they can pre-add
  // non-standard derived system classes, and if their restart file
  // has some non-standard sys_type we won't throw an error.
  if (_systems.count(name))
    {
      return this->get_system(name);
    }
  // Build a basic System
  else if (sys_type == "Basic")
    this->add_system<System> (name);

  // Build a Newmark system
  else if (sys_type == "Newmark")
    this->add_system<NewmarkSystem> (name);

  // Build an Explicit system
  else if ((sys_type == "Explicit"))
    this->add_system<ExplicitSystem> (name);

  // Build an Implicit system
  else if ((sys_type == "Implicit") ||
           (sys_type == "Steady"  ))
    this->add_system<ImplicitSystem> (name);

  // build a transient implicit linear system
  else if ((sys_type == "Transient") ||
           (sys_type == "TransientImplicit") ||
           (sys_type == "TransientLinearImplicit"))
    this->add_system<TransientLinearImplicitSystem> (name);

  // build a transient implicit nonlinear system
  else if (sys_type == "TransientNonlinearImplicit")
    this->add_system<TransientNonlinearImplicitSystem> (name);

  // build a transient explicit system
  else if (sys_type == "TransientExplicit")
    this->add_system<TransientExplicitSystem> (name);

  // build a linear implicit system
  else if (sys_type == "LinearImplicit")
    this->add_system<LinearImplicitSystem> (name);

  // build a nonlinear implicit system
  else if (sys_type == "NonlinearImplicit")
    this->add_system<NonlinearImplicitSystem> (name);

  // build a Reduced Basis Construction system
  else if (sys_type == "RBConstruction")
    this->add_system<RBConstruction> (name);

  // build a transient Reduced Basis Construction system
  else if (sys_type == "TransientRBConstruction")
    this->add_system<TransientRBConstruction> (name);

#ifdef LIBMESH_HAVE_SLEPC
  // build an eigen system
  else if (sys_type == "Eigen")
    this->add_system<EigenSystem> (name);
  else if (sys_type == "TransientEigenSystem")
    this->add_system<TransientEigenSystem> (name);
#endif

#if defined(LIBMESH_USE_COMPLEX_NUMBERS)
  // build a frequency system
  else if (sys_type == "Frequency")
    this->add_system<FrequencySystem> (name);
#endif

  else
    libmesh_error_msg("ERROR: Unknown system type: " << sys_type);

  // Return a reference to the new system
  //return (*this)(name);
  return this->get_system(name);
}



void EquationSystems::solve ()
{
  libmesh_assert (this->n_systems());

  for (auto i : make_range(this->n_systems()))
    this->get_system(i).solve();
}



void EquationSystems::sensitivity_solve (const ParameterVector & parameters_in)
{
  libmesh_assert (this->n_systems());

  for (auto i : make_range(this->n_systems()))
    this->get_system(i).sensitivity_solve(parameters_in);
}



void EquationSystems::adjoint_solve (const QoISet & qoi_indices)
{
  libmesh_assert (this->n_systems());

  for (unsigned int i=this->n_systems(); i != 0; --i)
    this->get_system(i-1).adjoint_solve(qoi_indices);
}



void EquationSystems::build_variable_names (std::vector<std::string> & var_names,
                                            const FEType * type,
                                            const std::set<std::string> * system_names) const
{
  // start indexing at end of possibly non-empty vector of variable names to avoid overwriting them
  unsigned int var_num = var_names.size();

  // We'll want to double-check that we don't have any naming
  // conflicts; this API causes problems down the line if so.
  std::unordered_multiset<std::string> seen_names;

  // Need to size var_names by scalar variables plus all the
  // vector components for all the vector variables
  //Could this be replaced by a/some convenience methods?[PB]
  {
    unsigned int n_scalar_vars = 0;
    unsigned int n_vector_vars = 0;

    for (const auto & [sys_name, sys_ptr] : _systems)
      {
        // Check current system is listed in system_names, and skip pos if not
        bool use_current_system = (system_names == nullptr);
        if (!use_current_system)
          use_current_system = system_names->count(sys_name);
        if (!use_current_system || sys_ptr->hide_output())
          {
            for (auto vn : make_range(sys_ptr->n_vars()))
              seen_names.insert(sys_ptr->variable_name(vn));
            continue;
          }

        for (auto vn : make_range(sys_ptr->n_vars()))
          {
            seen_names.insert(sys_ptr->variable_name(vn));
            if (FEInterface::field_type(sys_ptr->variable_type(vn)) == TYPE_VECTOR)
              n_vector_vars++;
            else
              n_scalar_vars++;
          }
      }

    // Here, we're assuming the number of vector components is the same
    // as the mesh dimension. Will break for mixed dimension meshes.
    unsigned int dim = this->get_mesh().mesh_dimension();
    unsigned int nv = n_scalar_vars + dim*n_vector_vars;

    // We'd better not have more than dim*his->n_vars() (all vector variables)
    // Treat the NodeElem-only mesh case as dim=1
    libmesh_assert_less_equal ( nv, (dim > 0 ? dim : 1)*this->n_vars() );

    // 'nv' represents the max possible number of output variables, so allocate enough memory for
    // all variables in the system to be populated here. When this is called more than once on a
    // single 'var_names' vector, different filters should be used such that duplicates don't occur.
    var_names.resize( nv );
  }

  for (const auto & [sys_name, sys_ptr] : _systems)
    {
      // Check current system is listed in system_names, and skip pos if not
      bool use_current_system = (system_names == nullptr);
      if (!use_current_system)
        use_current_system = system_names->count(sys_name);
      if (!use_current_system || sys_ptr->hide_output())
        continue;

      for (auto vn : make_range(sys_ptr->n_vars()))
        {
          const std::string & var_name = sys_ptr->variable_name(vn);
          const FEType & fe_type = sys_ptr->variable_type(vn);

          unsigned int n_vec_dim = FEInterface::n_vec_dim( sys_ptr->get_mesh(), fe_type);

          // Filter on the type if requested
          if (type == nullptr || (type && *type == fe_type))
            {
              if (FEInterface::field_type(fe_type) == TYPE_VECTOR)
                {
                  switch(n_vec_dim)
                    {
                    case 0:
                    case 1:
                      var_names[var_num++] = var_name;
                      libmesh_error_msg_if(seen_names.count(var_name) > 1,
                                           "Duplicate variable name "+var_name);
                      break;
                    case 2:
                      var_names[var_num++] = var_name+"_x";
                      var_names[var_num++] = var_name+"_y";
                      libmesh_error_msg_if(seen_names.count(var_name+"_x"),
                                           "Duplicate variable name "+var_name+"_x");
                      libmesh_error_msg_if(seen_names.count(var_name+"_y"),
                                           "Duplicate variable name "+var_name+"_y");
                      break;
                    case 3:
                      var_names[var_num++] = var_name+"_x";
                      var_names[var_num++] = var_name+"_y";
                      var_names[var_num++] = var_name+"_z";
                      libmesh_error_msg_if(seen_names.count(var_name+"_x"),
                                           "Duplicate variable name "+var_name+"_x");
                      libmesh_error_msg_if(seen_names.count(var_name+"_y"),
                                           "Duplicate variable name "+var_name+"_y");
                      libmesh_error_msg_if(seen_names.count(var_name+"_z"),
                                           "Duplicate variable name "+var_name+"_z");
                      break;
                    default:
                      libmesh_error_msg("Invalid dim in build_variable_names");
                    }
                }
              else
                var_names[var_num++] = var_name;
            }
        }
    }
  // Now resize again in case we filtered any names
  var_names.resize(var_num);
}



void EquationSystems::build_solution_vector (std::vector<Number> &,
                                             std::string_view,
                                             std::string_view) const
{
  // TODO:[BSK] re-implement this from the method below
  libmesh_not_implemented();
}




std::unique_ptr<NumericVector<Number>>
EquationSystems::build_parallel_solution_vector(const std::set<std::string> * system_names,
                                                bool add_sides) const
{
  LOG_SCOPE("build_parallel_solution_vector()", "EquationSystems");

  // This function must be run on all processors at once
  parallel_object_only();

  const unsigned int dim = _mesh.mesh_dimension();
  const dof_id_type max_nn   = _mesh.max_node_id();

  // allocate vector storage to hold
  // (max_node_id)*(number_of_variables) entries.
  //
  // If node renumbering is disabled and adaptive coarsening has
  // created gaps between node numbers, then this vector will be
  // sparse.
  //
  // We have to differentiate between between scalar and vector
  // variables. We intercept vector variables and treat each
  // component as a scalar variable (consistently with build_solution_names).

  unsigned int nv = 0;

  //Could this be replaced by a/some convenience methods?[PB]
  {
    unsigned int n_scalar_vars = 0;
    unsigned int n_vector_vars = 0;
    for (const auto & [sys_name, sys_ptr] : _systems)
      {
        // Check current system is listed in system_names, and skip pos if not
        bool use_current_system = (system_names == nullptr);
        if (!use_current_system)
          use_current_system = system_names->count(sys_name);
        if (!use_current_system || sys_ptr->hide_output())
          continue;

        for (auto vn : make_range(sys_ptr->n_vars()))
          {
            if (FEInterface::field_type(sys_ptr->variable_type(vn)) == TYPE_VECTOR)
              n_vector_vars++;
            else
              n_scalar_vars++;
          }
      }
    // Here, we're assuming the number of vector components is the same
    // as the mesh dimension. Will break for mixed dimension meshes.
    nv = n_scalar_vars + dim*n_vector_vars;
  }

  // Get the number of nodes to store locally.
  dof_id_type n_local_nodes = cast_int<dof_id_type>
    (std::distance(_mesh.local_nodes_begin(),
                   _mesh.local_nodes_end()));

  // If node renumbering has been disabled, nodes may not be numbered
  // contiguously, and the number of nodes might not match the
  // max_node_id.  In this case we just do our best.
  dof_id_type n_total_nodes = n_local_nodes;
  _mesh.comm().sum(n_total_nodes);

  const processor_id_type n_proc = _mesh.comm().size();
  const processor_id_type my_pid = _mesh.comm().rank();
  const dof_id_type n_gaps = max_nn - n_total_nodes;
  const dof_id_type gaps_per_processor = n_gaps / n_proc;
  const dof_id_type remainder_gaps = n_gaps % n_proc;

  n_local_nodes = n_local_nodes +      // Actual nodes
                  gaps_per_processor + // Our even share of gaps
                  (my_pid < remainder_gaps); // Leftovers

  // If we've been asked to build added sides' data, we need space to
  // add it.  Keep track of how much space.
  dof_id_type local_added_side_nodes = 0,
              added_side_nodes = 0;

  // others_added_side_nodes[p]: local_added_side_nodes on rank p
  std::vector<dof_id_type> others_added_side_nodes;

  // A map of (element_id, side, side_node) pairs to the corresponding
  // added side node index.
  std::map<std::tuple<dof_id_type, unsigned short, unsigned short>,
           dof_id_type> discontinuous_node_indices;

  // If we don't have any added side nodes, we'll have no offsets from
  // them, and we won't care about which offsets apply to which node
  // ids either.

  // Number of true nodes on processors [0,p]
  std::vector<dof_id_type> true_node_offsets;
  // Number of added (fake) nodes on processors [0,p)
  std::vector<dof_id_type> added_node_offsets;

  auto node_id_to_vec_id =
    [&true_node_offsets, &added_node_offsets]
    (const dof_id_type node_id)
    {
      if (true_node_offsets.empty())
        return node_id; // O(1) in the common !add_sides case

      // Find the processor id that has node_id in the parallel vec
      const auto lb = std::upper_bound(true_node_offsets.begin(),
                                       true_node_offsets.end(), node_id);
      libmesh_assert(lb != true_node_offsets.end());
      const processor_id_type p = lb - true_node_offsets.begin();

      return node_id + added_node_offsets[p];
    };

  if (add_sides)
    {
      true_node_offsets.resize(n_proc);
      added_node_offsets.resize(n_proc);

      // One loop to count everyone's new side nodes
      for (const auto & elem : _mesh.active_element_ptr_range())
        {
          for (auto s : elem->side_index_range())
            {
              if (redundant_added_side(*elem,s))
                continue;

              const std::vector<unsigned int> side_nodes =
                elem->nodes_on_side(s);

              if (elem->processor_id() == this->processor_id())
                local_added_side_nodes += side_nodes.size();
            }
        }

      others_added_side_nodes.resize(n_proc);
      _mesh.comm().allgather(local_added_side_nodes,
                             others_added_side_nodes);

      added_side_nodes = std::accumulate(others_added_side_nodes.begin(),
                                         others_added_side_nodes.end(), 0,
                                         std::plus<>());

      _mesh.comm().allgather(n_local_nodes, true_node_offsets);
      for (auto p : make_range(n_proc-1))
        true_node_offsets[p+1] += true_node_offsets[p];
      libmesh_assert_equal_to(true_node_offsets[n_proc-1], _mesh.max_node_id());

      // For nodes that exist in the mesh, we just need an offset to
      // tell where to put their solutions.
      added_node_offsets[0] = 0;
      for (auto p : make_range(n_proc-1))
        added_node_offsets[p+1] =
          added_node_offsets[p] + others_added_side_nodes[p];

      // For added side nodes, we need to fill a map.  Start after all
      // the true node for our pid plus all the side nodes for
      // previous pids
      dof_id_type node_counter = true_node_offsets[my_pid];
      for (auto p : make_range(my_pid))
        node_counter += others_added_side_nodes[p];

      // One loop to figure out whose added side nodes get which index
      for (const auto & elem : _mesh.active_local_element_ptr_range())
        {
          for (auto s : elem->side_index_range())
            {
              if (redundant_added_side(*elem,s))
                continue;

              const std::vector<unsigned int> side_nodes =
                elem->nodes_on_side(s);

              for (auto n : index_range(side_nodes))
                discontinuous_node_indices
                  [std::make_tuple(elem->id(),s,n)] = node_counter++;
            }
        }
    }

  const dof_id_type
    n_global_vals = (max_nn + added_side_nodes) * nv,
    n_local_vals = (n_local_nodes + local_added_side_nodes) * nv;

  // Create a NumericVector to hold the parallel solution
  std::unique_ptr<NumericVector<Number>> parallel_soln_ptr = NumericVector<Number>::build(_communicator);
  NumericVector<Number> & parallel_soln = *parallel_soln_ptr;
  parallel_soln.init(n_global_vals, n_local_vals, false, PARALLEL);

  // Create a NumericVector to hold the "repeat_count" for each node - this is essentially
  // the number of elements contributing to that node's value
  std::unique_ptr<NumericVector<Number>> repeat_count_ptr = NumericVector<Number>::build(_communicator);
  NumericVector<Number> & repeat_count = *repeat_count_ptr;
  repeat_count.init(n_global_vals, n_local_vals, false, PARALLEL);

  repeat_count.close();

  unsigned int var_num=0;

  // For each system in this EquationSystems object,
  // update the global solution and if we are on processor 0,
  // loop over the elements and build the nodal solution
  // from the element solution.  Then insert this nodal solution
  // into the vector passed to build_solution_vector.
  for (const auto & [sys_name, sys_ptr] : _systems)
    {
      // Check current system is listed in system_names, and skip pos if not
      bool use_current_system = (system_names == nullptr);
      if (!use_current_system)
        use_current_system = system_names->count(sys_name);
      if (!use_current_system || sys_ptr->hide_output())
        continue;

      const System & system  = *sys_ptr;
      const unsigned int nv_sys = system.n_vars();
      const unsigned int sys_num = system.number();

      //Could this be replaced by a/some convenience methods?[PB]
      unsigned int n_scalar_vars = 0;
      unsigned int n_vector_vars = 0;
      for (auto vn : make_range(sys_ptr->n_vars()))
        {
          if (FEInterface::field_type(sys_ptr->variable_type(vn)) == TYPE_VECTOR)
            n_vector_vars++;
          else
            n_scalar_vars++;
        }

      // Here, we're assuming the number of vector components is the same
      // as the mesh dimension. Will break for mixed dimension meshes.
      unsigned int nv_sys_split = n_scalar_vars + dim*n_vector_vars;

      // Update the current_local_solution
      {
        System & non_const_sys = const_cast<System &>(system);
        // We used to simply call non_const_sys.solution->close()
        // here, but that is not allowed when the solution vector is
        // locked read-only, for example when printing the solution
        // during the middle of a solve...  So try to be a bit
        // more careful about calling close() unnecessarily.
        libmesh_assert(this->comm().verify(non_const_sys.solution->closed()));
        if (!non_const_sys.solution->closed())
          non_const_sys.solution->close();
        non_const_sys.update();
      }

      NumericVector<Number> & sys_soln(*system.current_local_solution);

      const DofMap & dof_map = system.get_dof_map();

      std::vector<Number>      elem_soln;   // The finite element solution
      std::vector<Number>      nodal_soln;  // The FE solution interpolated to the nodes
      std::vector<dof_id_type> dof_indices; // The DOF indices for the finite element

      unsigned var_inc = 0;
      for (unsigned int var=0; var<nv_sys; var++)
        {
          const FEType & fe_type           = system.variable_type(var);
          const Variable & var_description = system.variable(var);
          unsigned int n_vec_dim = FEInterface::n_vec_dim( sys_ptr->get_mesh(), fe_type );
          const auto vg = dof_map.var_group_from_var_number(var);
          const bool add_p_level = dof_map.should_p_refine(vg);

          for (const auto & elem : _mesh.active_local_element_ptr_range())
            {
              if (var_description.active_on_subdomain(elem->subdomain_id()))
                {
                  dof_map.dof_indices (elem, dof_indices, var);
                  sys_soln.get(dof_indices, elem_soln);

                  FEInterface::nodal_soln (elem->dim(),
                                           fe_type,
                                           elem,
                                           elem_soln,
                                           nodal_soln,
                                           add_p_level);

                  // infinite elements should be skipped...
                  if (!elem->infinite())
                    {
                      libmesh_assert_equal_to (nodal_soln.size(), n_vec_dim*elem->n_nodes());

                      for (auto n : elem->node_index_range())
                        {
                          const Node & node = elem->node_ref(n);

                          const dof_id_type node_idx =
                            nv * node_id_to_vec_id(node.id());

                          for (unsigned int d=0; d < n_vec_dim; d++)
                            {
                              // For vector-valued elements, all components are in nodal_soln. For each
                              // node, the components are stored in order, i.e. node_0 -> s0_x, s0_y, s0_z
                              parallel_soln.add(node_idx + (var_inc+d + var_num), nodal_soln[n_vec_dim*n+d]);

                              // Increment the repeat count for this position
                              repeat_count.add(node_idx + (var_inc+d + var_num), 1);
                            }
                        }

                      if (add_sides)
                        {
                          for (auto s : elem->side_index_range())
                            {
                              if (redundant_added_side(*elem,s))
                                continue;

                              // Compute the FE solution at all the
                              // side nodes
                              FEInterface::side_nodal_soln
                                (fe_type, elem, s, elem_soln,
                                 nodal_soln, add_p_level);

#ifdef DEBUG
                              const std::vector<unsigned int> side_nodes =
                                elem->nodes_on_side(s);

                              libmesh_assert_equal_to
                                  (nodal_soln.size(),
                                   side_nodes.size());
#endif

                              for (auto n : index_range(nodal_soln))
                                {
                                  // Retrieve index into global solution vector.
                                  std::size_t node_index =
                                    nv * libmesh_map_find(discontinuous_node_indices,
                                                          std::make_tuple(elem->id(), s, n));

                                  for (unsigned int d=0; d < n_vec_dim; d++)
                                    {
                                      parallel_soln.add(node_index + (var_inc+d + var_num), nodal_soln[n_vec_dim*n+d]);
                                      repeat_count.add(node_index + (var_inc+d + var_num), 1);
                                    }
                                }
                            }
                        }
                    }
                }
              else // If this variable doesn't exist on this subdomain we have to still increment repeat_count so that we won't divide by 0 later:
                for (auto n : elem->node_index_range())
                  {
                    const Node & node = elem->node_ref(n);
                    // Only do this if this variable has NO DoFs at
                    // this node... it might have some from an
                    // adjoining element...
                    if (!node.n_dofs(sys_num, var))
                      {
                        const dof_id_type node_idx =
                          nv * node_id_to_vec_id(node.id());

                        for (unsigned int d=0; d < n_vec_dim; d++)
                          repeat_count.add(node_idx + (var_inc+d + var_num), 1);
                      }
                  }

            } // end loop over elements
          var_inc += n_vec_dim;
        } // end loop on variables in this system

      var_num += nv_sys_split;
    } // end loop over systems

  // Sum the nodal solution values and repeat counts.
  parallel_soln.close();
  repeat_count.close();

  // If there were gaps in the node numbering, there will be
  // corresponding zeros in the parallel_soln and repeat_count
  // vectors.  We need to set those repeat_count entries to 1
  // in order to avoid dividing by zero.
  if (n_gaps)
    {
      for (numeric_index_type i=repeat_count.first_local_index();
           i<repeat_count.last_local_index(); ++i)
        {
          // repeat_count entries are integral values but let's avoid a
          // direct floating point comparison with 0 just in case some
          // roundoff noise crept in during vector assembly?
          if (std::abs(repeat_count(i)) < TOLERANCE)
            repeat_count.set(i, 1.);
        }

      // Make sure the repeat_count vector is up-to-date on all
      // processors.
      repeat_count.close();
    }

  // Divide to get the average value at the nodes
  parallel_soln /= repeat_count;

  return parallel_soln_ptr;
}



void EquationSystems::build_solution_vector (std::vector<Number> & soln,
                                             const std::set<std::string> * system_names,
                                             bool add_sides) const
{
  LOG_SCOPE("build_solution_vector()", "EquationSystems");

  // Call the parallel implementation
  std::unique_ptr<NumericVector<Number>> parallel_soln =
    this->build_parallel_solution_vector(system_names, add_sides);

  // Localize the NumericVector into the provided std::vector.
  parallel_soln->localize_to_one(soln);
}



void EquationSystems::get_vars_active_subdomains(const std::vector<std::string> & names,
                                                 std::vector<std::set<subdomain_id_type>> & vars_active_subdomains) const
{
  vars_active_subdomains.clear();
  vars_active_subdomains.resize(names.size());

  for (const auto & pr : _systems)
    {
      const auto & sys_ptr = pr.second;
      for (auto vn : make_range(sys_ptr->n_vars()))
        {
          const std::string & var_name = sys_ptr->variable_name(vn);

          auto names_it = std::find(names.begin(), names.end(), var_name);
          if(names_it != names.end())
            {
              const Variable & variable = sys_ptr->variable(vn);
              const std::set<subdomain_id_type> & active_subdomains = variable.active_subdomains();
              vars_active_subdomains[std::distance(names.begin(), names_it)] = active_subdomains;
            }
        }
    }
}



void EquationSystems::get_solution (std::vector<Number> & soln,
                                    std::vector<std::string> & names) const
{
  libmesh_deprecated();
  this->build_elemental_solution_vector(soln, names);
}



void
EquationSystems::build_elemental_solution_vector (std::vector<Number> & soln,
                                                  std::vector<std::string> & names) const
{
  // Call the parallel version of this function
  std::unique_ptr<NumericVector<Number>> parallel_soln =
    this->build_parallel_elemental_solution_vector(names);

  // Localize into 'soln', provided that parallel_soln is not empty.
  // Note: parallel_soln will be empty in the event that none of the
  // input names were CONSTANT, MONOMIAL nor components of CONSTANT,
  // MONOMIAL_VEC variables, or there were simply none of these in
  // the EquationSystems object.
  soln.clear();
  if (parallel_soln)
    parallel_soln->localize_to_one(soln);
}

std::vector<std::pair<unsigned int, unsigned int>>
EquationSystems::find_variable_numbers
  (std::vector<std::string> & names, const FEType * type, const std::vector<FEType> * types) const
{
  // This function must be run on all processors at once
  parallel_object_only();

  libmesh_assert (this->n_systems());

  // Resolve class of type input and assert that at least one of them is null
  libmesh_assert_msg(!type || !types,
                     "Input 'type', 'types', or neither in find_variable_numbers, but not both.");

  std::vector<FEType> type_filter;
  if (type)
    type_filter.push_back(*type);
  else if (types)
    type_filter = *types;

  // Store a copy of the valid variable names, if any. The names vector will be repopulated with any
  // valid names (or all if 'is_names_empty') in the system that passes through the type filter. If
  // the variable is a vector, its name will be decomposed into its separate components in
  // accordance with build_variable_names().
  std::vector<std::string> name_filter = names;
  bool is_names_empty = name_filter.empty();
  names.clear();

  // initialize convenience variables
  FEType var_type;
  std::string name;

  const std::vector<std::string> component_suffix = {"_x", "_y", "_z"};
  unsigned int dim = _mesh.mesh_dimension();
  libmesh_error_msg_if(dim > 3, "Invalid dim in find_variable_numbers");

  // Now filter through the variables in each system and store the system index and their index
  // within that system. This way, we know where to find their data even after we sort them.
  std::vector<std::pair<unsigned int, unsigned int>> var_nums;

  for (const auto & pr : _systems)
    {
      const System & system = *(pr.second);

      for (auto var : make_range(system.n_vars()))
        {
          // apply the type filter
          var_type = system.variable_type(var);
          if (type_filter.size() &&
              std::find(type_filter.begin(), type_filter.end(), var_type) == type_filter.end())
            continue;

          // apply the name filter (note that all variables pass if it is empty)
          if (FEInterface::field_type(var_type) == TYPE_VECTOR)
            {
              std::vector<std::string> component_names;
              for (unsigned int comp = 0; comp < dim; ++comp)
                {
                  name = system.variable_name(var) + component_suffix[comp];
                  if (is_names_empty ||
                      (std::find(name_filter.begin(), name_filter.end(), name) != name_filter.end()))
                    component_names.push_back(name);
                }

              if (! component_names.empty())
                names.insert(names.end(), component_names.begin(), component_names.end());
              else
                continue;
            }
          else /*scalar-valued variable*/
            {
              name = system.variable_name(var);
              if (is_names_empty ||
                  (std::find(name_filter.begin(), name_filter.end(), name) != name_filter.end()))
                names.push_back(name);
              else
                continue;
            }

          // if the variable made it through both filters get its system indices
          var_nums.emplace_back(system.number(), var);
        }
    }

  // Sort the var_nums vector pairs alphabetically based on the variable name
  std::vector<unsigned int> sort_index(var_nums.size());
  std::iota(sort_index.begin(), sort_index.end(), 0);
  std::sort(sort_index.begin(), sort_index.end(),
            [&](const unsigned int & lhs, const unsigned int & rhs)
            {return this->get_system(var_nums[lhs].first).variable_name(var_nums[lhs].second) <
                    this->get_system(var_nums[rhs].first).variable_name(var_nums[rhs].second);});

  std::vector<std::pair<unsigned int, unsigned int>> var_nums_sorted(var_nums.size());
  for (auto i : index_range(var_nums_sorted))
    {
      var_nums_sorted[i].first = var_nums[sort_index[i]].first;
      var_nums_sorted[i].second = var_nums[sort_index[i]].second;
    }

  // Also sort the names vector
  std::sort(names.begin(), names.end());

  // Return the sorted vector pairs
  return var_nums_sorted;
}


std::unique_ptr<NumericVector<Number>>
EquationSystems::build_parallel_elemental_solution_vector (std::vector<std::string> & names) const
{
  // Filter any names that aren't elemental variables and get the system indices for those that are.
  // Note that it's probably fine if the names vector is empty since we'll still at least filter
  // out all non-monomials. If there are no monomials, then nothing is output here.
  std::vector<FEType> type = {FEType(CONSTANT, MONOMIAL), FEType(CONSTANT, MONOMIAL_VEC)};
  std::vector<std::pair<unsigned int, unsigned int>> var_nums =
    this->find_variable_numbers(names, /*type=*/nullptr, &type);

  const std::size_t nv = names.size(); /*total number of vars including vector components*/
  const dof_id_type ne = _mesh.n_elem();
  libmesh_assert_equal_to (ne, _mesh.max_elem_id());

  // If there are no variables to write out don't do anything...
  if (!nv)
    return std::unique_ptr<NumericVector<Number>>(nullptr);

  // We can handle the case where there are nullptrs in the Elem vector
  // by just having extra zeros in the solution vector.
  numeric_index_type parallel_soln_global_size = ne*nv;

  numeric_index_type div = parallel_soln_global_size / this->n_processors();
  numeric_index_type mod = parallel_soln_global_size % this->n_processors();

  // Initialize all processors to the average size.
  numeric_index_type parallel_soln_local_size = div;

  // The first "mod" processors get an extra entry.
  if (this->processor_id() < mod)
    parallel_soln_local_size = div+1;

  // Create a NumericVector to hold the parallel solution
  std::unique_ptr<NumericVector<Number>> parallel_soln_ptr = NumericVector<Number>::build(_communicator);
  NumericVector<Number> & parallel_soln = *parallel_soln_ptr;
  parallel_soln.init(parallel_soln_global_size,
                     parallel_soln_local_size,
                     /*fast=*/false,
                     /*ParallelType=*/PARALLEL);

  unsigned int sys_ctr = 0;
  unsigned int var_ctr = 0;
  for (auto i : index_range(var_nums))
    {
      std::pair<unsigned int, unsigned int> var_num = var_nums[i];
      const System & system = this->get_system(var_num.first);

      // Update the current_local_solution if necessary
      if (sys_ctr != var_num.first)
        {
          System & non_const_sys = const_cast<System &>(system);
          // We used to simply call non_const_sys.solution->close()
          // here, but that is not allowed when the solution vector is
          // locked read-only, for example when printing the solution
          // during during the middle of a solve...  So try to be a bit
          // more careful about calling close() unnecessarily.
          libmesh_assert(this->comm().verify(non_const_sys.solution->closed()));
          if (!non_const_sys.solution->closed())
            non_const_sys.solution->close();
          non_const_sys.update();
          sys_ctr = var_num.first;
        }

      NumericVector<Number> & sys_soln(*system.current_local_solution);

      // The DOF indices for the finite element
      std::vector<dof_id_type> dof_indices;

      const unsigned int var = var_num.second;

      const Variable & variable = system.variable(var);
      const DofMap & dof_map = system.get_dof_map();

      // We need to check if the constant monomial is a scalar or a vector and set the number of
      // components as the mesh dimension for the latter case as per 'find_variable_numbers()'.
      // Even for the case where a variable is not active on any subdomain belonging to the
      // processor, we still need to know this number to update 'var_ctr'.
      const unsigned int n_comps =
        (system.variable_type(var) == type[1]) ? _mesh.mesh_dimension() : 1;

      // Loop over all elements in the mesh and index all components of the variable if it's active
      for (const auto & elem : _mesh.active_local_element_ptr_range())
        if (variable.active_on_subdomain(elem->subdomain_id()))
          {
            dof_map.dof_indices(elem, dof_indices, var);

            // The number of DOF components needs to be equal to the expected number so that we know
            // where to store data to correctly correspond to variable names.
            libmesh_assert_equal_to(dof_indices.size(), n_comps);

            for (unsigned int comp = 0; comp < n_comps; comp++)
              parallel_soln.set(ne * (var_ctr + comp) + elem->id(), sys_soln(dof_indices[comp]));
          }

      var_ctr += n_comps;
    } // end loop over var_nums

  // NOTE: number of output names might not be equal to the number passed to this function. Any that
  // aren't CONSTANT MONOMIALS or components of CONSTANT MONOMIAL_VECS have been filtered out (see
  // EquationSystems::find_variable_numbers).
  //
  // But, if everything is accounted for properly, then names.size() == var_ctr
  libmesh_assert_equal_to(names.size(), var_ctr);

  parallel_soln.close();
  return parallel_soln_ptr;
}



void
EquationSystems::build_discontinuous_solution_vector
(std::vector<Number> & soln,
 const std::set<std::string> * system_names,
 const std::vector<std::string> * var_names,
 bool vertices_only,
 bool add_sides) const
{
  LOG_SCOPE("build_discontinuous_solution_vector()", "EquationSystems");

  libmesh_assert (this->n_systems());

  // Get the number of variables (nv) by counting the number of variables
  // in each system listed in system_names
  unsigned int nv = 0;

  for (const auto & [sys_name, sys_ptr] : _systems)
    {
      // Check current system is listed in system_names, and skip pos if not
      bool use_current_system = (system_names == nullptr);
      if (!use_current_system)
        use_current_system = system_names->count(sys_name);
      if (!use_current_system || sys_ptr->hide_output())
        continue;

      // Loop over all variables in this System and check whether we
      // are supposed to use each one.
      for (auto var_id : make_range(sys_ptr->n_vars()))
        {
          bool use_current_var = (var_names == nullptr);
          if (!use_current_var)
            use_current_var = std::count(var_names->begin(),
                                         var_names->end(),
                                         sys_ptr->variable_name(var_id));

          // Only increment the total number of vars if we are
          // supposed to use this one.
          if (use_current_var)
            nv++;
        }
    }

  // get the total "weight" - the number of nodal values to write for
  // each variable.
  unsigned int tw=0;
  for (const auto & elem : _mesh.active_element_ptr_range())
    {
      tw += vertices_only ? elem->n_vertices() : elem->n_nodes();

      if (add_sides)
        {
          for (auto s : elem->side_index_range())
            {
              if (redundant_added_side(*elem,s))
                continue;

              const std::vector<unsigned int> side_nodes =
                elem->nodes_on_side(s);

              if (!vertices_only)
                tw += side_nodes.size();
              else
                for (auto n : index_range(side_nodes))
                  if (elem->is_vertex(side_nodes[n]))
                    ++tw;
            }
        }
    }

  // Only if we are on processor zero, allocate the storage
  // to hold (number_of_nodes)*(number_of_variables) entries.
  if (_mesh.processor_id() == 0)
    soln.resize(tw*nv);

  std::vector<Number> sys_soln;

  // Keep track of the variable "offset". This is used for indexing
  // into the global solution vector.
  unsigned int var_offset = 0;

  // For each system in this EquationSystems object,
  // update the global solution and if we are on processor 0,
  // loop over the elements and build the nodal solution
  // from the element solution.  Then insert this nodal solution
  // into the vector passed to build_solution_vector.
  for (const auto & [sys_name, system] : _systems)
    {
      // Check current system is listed in system_names, and skip pos if not
      bool use_current_system = (system_names == nullptr);
      if (!use_current_system)
        use_current_system = system_names->count(sys_name);
      if (!use_current_system || system->hide_output())
        continue;

      const unsigned int nv_sys = system->n_vars();
      const auto & dof_map = system->get_dof_map();

      system->update_global_solution (sys_soln, 0);

      // Keep track of the number of vars actually written.
      unsigned int n_vars_written_current_system = 0;

      if (_mesh.processor_id() == 0)
        {
          std::vector<Number>       soln_coeffs; // The finite element solution coeffs
          std::vector<Number>       nodal_soln;  // The FE solution interpolated to the nodes
          std::vector<dof_id_type>  dof_indices; // The DOF indices for the finite element

          // For each variable, determine if we are supposed to
          // write it, then loop over the active elements, compute
          // the nodal_soln and store it to the "soln" vector. We
          // store zeros for subdomain-restricted variables on
          // elements where they are not active.
          for (unsigned int var=0; var<nv_sys; var++)
            {
              bool use_current_var = (var_names == nullptr);
              if (!use_current_var)
                use_current_var = std::count(var_names->begin(),
                                             var_names->end(),
                                             system->variable_name(var));

              // If we aren't supposed to write this var, go to the
              // next loop iteration.
              if (!use_current_var)
                continue;

              const FEType & fe_type = system->variable_type(var);
              const Variable & var_description = system->variable(var);
              const auto vg = dof_map.var_group_from_var_number(var);
              const bool add_p_level = dof_map.should_p_refine(vg);

              unsigned int nn=0;

              for (auto & elem : _mesh.active_element_ptr_range())
                {
                  if (var_description.active_on_subdomain(elem->subdomain_id()))
                    {
                      system->get_dof_map().dof_indices (elem, dof_indices, var);

                      soln_coeffs.resize(dof_indices.size());

                      for (auto i : index_range(dof_indices))
                        soln_coeffs[i] = sys_soln[dof_indices[i]];

                      // Compute the FE solution at all the nodes, but
                      // only use the first n_vertices() entries if
                      // vertices_only == true.
                      FEInterface::nodal_soln (elem->dim(),
                                               fe_type,
                                               elem,
                                               soln_coeffs,
                                               nodal_soln,
                                               add_p_level);

                      // infinite elements should be skipped...
                      if (!elem->infinite())
                        {
                          libmesh_assert_equal_to (nodal_soln.size(), elem->n_nodes());

                          const unsigned int n_vals =
                            vertices_only ? elem->n_vertices() : elem->n_nodes();

                          for (unsigned int n=0; n<n_vals; n++)
                            {
                              // Compute index into global solution vector.
                              std::size_t index =
                                nv * (nn++) + (n_vars_written_current_system + var_offset);

                              soln[index] += nodal_soln[n];
                            }
                        }
                    }
                  else
                    nn += vertices_only ? elem->n_vertices() : elem->n_nodes();
                } // end loop over active elements writing interiors

              // Loop writing "fake" sides, if requested
              if (add_sides)
                {
                  // We don't build discontinuous solution vectors in
                  // parallel yet, but we'll do ordering of fake side
                  // values as if we did, for consistency with the
                  // parallel continuous ordering and for future
                  // compatibility.
                  std::vector<std::vector<const Elem *>>
                    elems_by_pid(_mesh.n_processors());

                  for (const auto & elem : _mesh.active_element_ptr_range())
                    elems_by_pid[elem->processor_id()].push_back(elem);

                  for (auto p : index_range(elems_by_pid))
                    for (const Elem * elem : elems_by_pid[p])
                      {
                        if (var_description.active_on_subdomain(elem->subdomain_id()))
                          {
                            system->get_dof_map().dof_indices (elem, dof_indices, var);

                            soln_coeffs.resize(dof_indices.size());

                            for (auto i : index_range(dof_indices))
                              soln_coeffs[i] = sys_soln[dof_indices[i]];

                            for (auto s : elem->side_index_range())
                              {
                                if (redundant_added_side(*elem,s))
                                  continue;

                                const std::vector<unsigned int> side_nodes =
                                  elem->nodes_on_side(s);

                                // Compute the FE solution at all the
                                // side nodes, but only use those for
                                // which is_vertex() == true if
                                // vertices_only == true.
                                FEInterface::side_nodal_soln
                                  (fe_type, elem, s, soln_coeffs,
                                   nodal_soln, add_p_level);

                                libmesh_assert_equal_to
                                    (nodal_soln.size(),
                                     side_nodes.size());

                                // If we don't have a continuous FE
                                // then we want to average between
                                // sides, at least in the equal-level
                                // case where it's easy.  This is
                                // analogous to our repeat_count
                                // behavior elsewhere.
                                const FEContinuity cont =
                                  FEInterface::get_continuity(fe_type);
                                const Elem * const neigh = elem->neighbor_ptr(s);

                                if ((cont == DISCONTINUOUS || cont == H_CURL) &&
                                    neigh &&
                                    neigh->level() == elem->level() &&
                                    var_description.active_on_subdomain(neigh->subdomain_id()))
                                  {
                                    std::vector<dof_id_type> neigh_indices;
                                    system->get_dof_map().dof_indices (neigh, neigh_indices, var);
                                    std::vector<Number> neigh_coeffs(neigh_indices.size());

                                    for (auto i : index_range(neigh_indices))
                                      neigh_coeffs[i] = sys_soln[neigh_indices[i]];

                                    const unsigned int s_neigh =
                                      neigh->which_neighbor_am_i(elem);
                                    std::vector<Number> neigh_soln;
                                    FEInterface::side_nodal_soln
                                      (fe_type, neigh, s_neigh,
                                       neigh_coeffs, neigh_soln, add_p_level);

                                    const std::vector<unsigned int> neigh_nodes =
                                      neigh->nodes_on_side(s_neigh);
                                    for (auto n : index_range(side_nodes))
                                      for (auto neigh_n : index_range(neigh_nodes))
                                        if (neigh->node_ptr(neigh_nodes[neigh_n])
                                            == elem->node_ptr(side_nodes[n]))
                                          {
                                            nodal_soln[n] += neigh_soln[neigh_n];
                                            nodal_soln[n] /= 2;
                                          }
                                  }

                                for (auto n : index_range(side_nodes))
                                  {
                                    if (vertices_only &&
                                        !elem->is_vertex(n))
                                      continue;

                                    // Compute index into global solution vector.
                                    std::size_t index =
                                      nv * (nn++) + (n_vars_written_current_system + var_offset);

                                    soln[index] += nodal_soln[n];
                                  }
                              }
                          }
                        else
                          {
                            nn += vertices_only ? elem->n_vertices() : elem->n_nodes();

                            for (auto s : elem->side_index_range())
                              {
                                if (redundant_added_side(*elem,s))
                                  continue;

                                const std::vector<unsigned int> side_nodes =
                                  elem->nodes_on_side(s);

                                for (auto n : index_range(side_nodes))
                                  {
                                    if (vertices_only &&
                                        !elem->is_vertex(n))
                                      continue;
                                    nn++;
                                  }
                              }
                          }
                      } // end loop over active elements, writing "fake" sides
                }
              // If we made it here, we actually wrote a variable, so increment
              // the number of variables actually written for the current system.
              n_vars_written_current_system++;

            } // end loop over vars
        } // end if proc 0

      // Update offset for next loop iteration.
      var_offset += n_vars_written_current_system;
    } // end loop over systems
}



bool EquationSystems::redundant_added_side(const Elem & elem, unsigned int side)
{
  libmesh_assert(elem.active());

  const Elem * neigh = elem.neighbor_ptr(side);

  // Write boundary sides.
  if (!neigh)
    return false;

  // Write ghost sides in Nemesis
  if (neigh == remote_elem)
    return false;

  // Don't write a coarser side if a finer side exists
  if (!neigh->active())
    return true;

  // Don't write a side redundantly from both of the
  // elements sharing it.  We'll disambiguate with id().
  return (neigh->id() < elem.id());
}



bool EquationSystems::compare (const EquationSystems & other_es,
                               const Real threshold,
                               const bool verbose) const
{
  // safety check, whether we handle at least the same number
  // of systems
  std::vector<bool> os_result;

  if (this->n_systems() != other_es.n_systems())
    {
      if (verbose)
        {
          libMesh::out << "  Fatal difference. This system handles "
                       << this->n_systems() << " systems," << std::endl
                       << "  while the other system handles "
                       << other_es.n_systems()
                       << " systems." << std::endl
                       << "  Aborting comparison." << std::endl;
        }
      return false;
    }
  else
    {
      // start comparing each system
      for (const auto & [sys_name, sys_ptr] : _systems)
        {
          // get the other system
          const System & other_system   = other_es.get_system (sys_name);

          os_result.push_back (sys_ptr->compare (other_system, threshold, verbose));

        }

    }


  // sum up the results
  if (os_result.size()==0)
    return true;
  else
    {
      bool os_identical;
      unsigned int n = 0;
      do
        {
          os_identical = os_result[n];
          n++;
        }
      while (os_identical && n<os_result.size());
      return os_identical;
    }
}



std::string EquationSystems::get_info () const
{
  std::ostringstream oss;

  oss << " EquationSystems\n"
      << "  n_systems()=" << this->n_systems() << '\n';

  // Print the info for the individual systems
  for (const auto & pr : _systems)
    oss << pr.second->get_info();


  //   // Possibly print the parameters
  //   if (!this->parameters.empty())
  //     {
  //       oss << "  n_parameters()=" << this->n_parameters() << '\n';
  //       oss << "   Parameters:\n";

  //       for (const auto & [key, val] : _parameters)
  //         oss << "    "
  //             << "\""
  //             << key
  //             << "\""
  //             << "="
  //             << val
  //             << '\n';
  //     }

  return oss.str();
}



void EquationSystems::print_info (std::ostream & os) const
{
  os << this->get_info()
     << std::endl;
}



std::ostream & operator << (std::ostream & os,
                            const EquationSystems & es)
{
  es.print_info(os);
  return os;
}



unsigned int EquationSystems::n_vars () const
{
  unsigned int tot=0;

  for (const auto & pr : _systems)
    tot += pr.second->n_vars();

  return tot;
}



std::size_t EquationSystems::n_dofs () const
{
  std::size_t tot=0;

  for (const auto & pr : _systems)
    tot += pr.second->n_dofs();

  return tot;
}




std::size_t EquationSystems::n_active_dofs () const
{
  std::size_t tot=0;

  for (const auto & pr : _systems)
    tot += pr.second->n_active_dofs();

  return tot;
}


void EquationSystems::_add_system_to_nodes_and_elems()
{
  // All the nodes
  for (auto & node : _mesh.node_ptr_range())
    node->add_system();

  // All the elements
  for (auto & elem : _mesh.element_ptr_range())
    elem->add_system();
}

void EquationSystems::_remove_default_ghosting(unsigned int sys_num)
{
  this->get_system(sys_num).get_dof_map().remove_default_ghosting();
}

} // namespace libMesh
