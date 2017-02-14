// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// System includes
#include <sstream>

// Local Includes
#include "libmesh/explicit_system.h"
#include "libmesh/fe_interface.h"
#include "libmesh/frequency_system.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/newmark_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/rb_construction.h"
#include "libmesh/transient_rb_construction.h"
#include "libmesh/eigen_system.h"
#include "libmesh/parallel.h"
#include "libmesh/transient_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/mesh_base.h"
#include "libmesh/elem.h"
#include "libmesh/libmesh_logging.h"

// Include the systems before this one to avoid
// overlapping forward declarations.
#include "libmesh/equation_systems.h"

namespace libMesh
{

// Forward Declarations




// ------------------------------------------------------------
// EquationSystems class implementation
EquationSystems::EquationSystems (MeshBase & m, MeshData * mesh_data) :
  ParallelObject (m),
  _mesh          (m),
  _mesh_data     (mesh_data)
{
  // Set default parameters
  this->parameters.set<Real>        ("linear solver tolerance") = TOLERANCE * TOLERANCE;
  this->parameters.set<unsigned int>("linear solver maximum iterations") = 5000;
  this->_refine_in_reinit = true; // default value
}



EquationSystems::~EquationSystems ()
{
  this->clear ();
}



void EquationSystems::clear ()
{
  // Clear any additional parameters
  parameters.clear ();

  // clear the systems.  We must delete them
  // since we newed them!
  while (!_systems.empty())
    {
      system_iterator pos = _systems.begin();

      System * sys = pos->second;
      delete sys;
      sys = libmesh_nullptr;

      _systems.erase (pos);
    }
}



void EquationSystems::init ()
{
  const unsigned int n_sys = this->n_systems();

  libmesh_assert_not_equal_to (n_sys, 0);

  // Distribute the mesh if possible
  if (this->n_processors() > 1)
    _mesh.delete_remote_elements();

  // Tell all the \p DofObject entities how many systems
  // there are.
  {
    MeshBase::node_iterator       node_it  = _mesh.nodes_begin();
    const MeshBase::node_iterator node_end = _mesh.nodes_end();

    for ( ; node_it != node_end; ++node_it)
      (*node_it)->set_n_systems(n_sys);

    MeshBase::element_iterator       elem_it  = _mesh.elements_begin();
    const MeshBase::element_iterator elem_end = _mesh.elements_end();

    for ( ; elem_it != elem_end; ++elem_it)
      (*elem_it)->set_n_systems(n_sys);
  }

  for (unsigned int i=0; i != this->n_systems(); ++i)
    this->get_system(i).init();

#ifdef LIBMESH_ENABLE_AMR
  MeshRefinement mesh_refine(_mesh);
  mesh_refine.clean_refinement_flags();
#endif
}



void EquationSystems::reinit ()
{
  parallel_object_only();

  const unsigned int n_sys = this->n_systems();
  libmesh_assert_not_equal_to (n_sys, 0);

  // We may have added new systems since our last
  // EquationSystems::(re)init call
  bool _added_new_systems = false;
  for (unsigned int i=0; i != n_sys; ++i)
    if (!this->get_system(i).is_initialized())
      _added_new_systems = true;

  if (_added_new_systems)
    {
      // Our DofObjects will need space for the additional systems
      MeshBase::node_iterator       node_it  = _mesh.nodes_begin();
      const MeshBase::node_iterator node_end = _mesh.nodes_end();

      for ( ; node_it != node_end; ++node_it)
        (*node_it)->set_n_systems(n_sys);

      MeshBase::element_iterator       elem_it  = _mesh.elements_begin();
      const MeshBase::element_iterator elem_end = _mesh.elements_end();

      for ( ; elem_it != elem_end; ++elem_it)
        (*elem_it)->set_n_systems(n_sys);

      // And any new systems will need initialization
      for (unsigned int i=0; i != n_sys; ++i)
        if (!this->get_system(i).is_initialized())
          this->get_system(i).init();
    }


  // We used to assert that all nodes and elements *already* had
  // n_systems() properly set; however this is false in the case where
  // user code has manually added nodes and/or elements to an
  // already-initialized system.

  // Make sure all the \p DofObject entities know how many systems
  // there are.
  {
    // All the nodes
    MeshBase::node_iterator       node_it  = _mesh.nodes_begin();
    const MeshBase::node_iterator node_end = _mesh.nodes_end();

    for ( ; node_it != node_end; ++node_it)
      {
        Node * node = *node_it;
        node->set_n_systems(this->n_systems());
      }

    // All the elements
    MeshBase::element_iterator       elem_it  = _mesh.elements_begin();
    const MeshBase::element_iterator elem_end = _mesh.elements_end();

    for ( ; elem_it != elem_end; ++elem_it)
      {
        Elem * elem = *elem_it;
        elem->set_n_systems(this->n_systems());
      }
  }

  // Localize each system's vectors
  for (unsigned int i=0; i != this->n_systems(); ++i)
    this->get_system(i).re_update();

#ifdef LIBMESH_ENABLE_AMR

  bool dof_constraints_created = false;
  bool mesh_changed = false;

  // FIXME: For backwards compatibility, assume
  // refine_and_coarsen_elements or refine_uniformly have already
  // been called
  {
    for (unsigned int i=0; i != this->n_systems(); ++i)
      {
        System & sys = this->get_system(i);

        // Even if the system doesn't have any variables in it we want
        // consistent behavior; e.g. distribute_dofs should have the
        // opportunity to count up zero dofs on each processor.
        //
        // Who's been adding zero-var systems anyway, outside of my
        // unit tests? - RHS
        // if(!sys.n_vars())
        // continue;

        sys.get_dof_map().distribute_dofs(_mesh);

        // Recreate any user or internal constraints
        sys.reinit_constraints();

        sys.prolong_vectors();
      }
    mesh_changed = true;
    dof_constraints_created = true;
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
          for (unsigned int i=0; i != this->n_systems(); ++i)
            {
              System & sys = this->get_system(i);
              if (!dof_constraints_created)
                {
                  sys.get_dof_map().distribute_dofs(_mesh);
                  sys.reinit_constraints();
                }
              sys.restrict_vectors();
            }
          mesh_changed = true;
          dof_constraints_created = true;
        }

      // Once vectors are all restricted, we can delete
      // children of coarsened elements
      if (mesh_changed)
        this->get_mesh().contract();

      // Try to refine the mesh, then prolong each system's vectors
      // if necessary
      if (mesh_refine.refine_elements())
        {
          for (unsigned int i=0; i != this->n_systems(); ++i)
            {
              System & sys = this->get_system(i);
              if (!dof_constraints_created)
                {
                  sys.get_dof_map().distribute_dofs(_mesh);
                  sys.reinit_constraints();
                }
              sys.prolong_vectors();
            }
          mesh_changed = true;
          // dof_constraints_created = true;
        }
    }

  // If the mesh has changed, systems will need to create new dof
  // constraints and update their global solution vectors
  if (mesh_changed)
    {
      for (unsigned int i=0; i != this->n_systems(); ++i)
        this->get_system(i).reinit();
    }
#endif // #ifdef LIBMESH_ENABLE_AMR
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
  {
    MeshBase::node_iterator       node_it  = _mesh.nodes_begin();
    const MeshBase::node_iterator node_end = _mesh.nodes_end();

    for ( ; node_it != node_end; ++node_it)
      (*node_it)->set_n_systems(n_sys);

    MeshBase::element_iterator       elem_it  = _mesh.elements_begin();
    const MeshBase::element_iterator elem_end = _mesh.elements_end();

    for ( ; elem_it != elem_end; ++elem_it)
      (*elem_it)->set_n_systems(n_sys);
  }

  // And distribute each system's dofs
  for (unsigned int i=0; i != this->n_systems(); ++i)
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




void EquationSystems::update ()
{
  LOG_SCOPE("update()", "EquationSystems");

  // Localize each system's vectors
  for (unsigned int i=0; i != this->n_systems(); ++i)
    this->get_system(i).update();
}



System & EquationSystems::add_system (const std::string & sys_type,
                                      const std::string & name)
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






void EquationSystems::delete_system (const std::string & name)
{
  libmesh_deprecated();

  if (!_systems.count(name))
    libmesh_error_msg("ERROR: no system named " << name);

  delete _systems[name];

  _systems.erase (name);
}



void EquationSystems::solve ()
{
  libmesh_assert (this->n_systems());

  for (unsigned int i=0; i != this->n_systems(); ++i)
    this->get_system(i).solve();
}



void EquationSystems::sensitivity_solve (const ParameterVector & parameters_in)
{
  libmesh_assert (this->n_systems());

  for (unsigned int i=0; i != this->n_systems(); ++i)
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
  unsigned int var_num=0;

  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  // Need to size var_names by scalar variables plus all the
  // vector components for all the vector variables
  //Could this be replaced by a/some convenience methods?[PB]
  {
    unsigned int n_scalar_vars = 0;
    unsigned int n_vector_vars = 0;

    for (; pos != end; ++pos)
      {
        // Check current system is listed in system_names, and skip pos if not
        bool use_current_system = (system_names == libmesh_nullptr);
        if (!use_current_system)
          use_current_system = system_names->count(pos->first);
        if (!use_current_system)
          continue;

        for (unsigned int vn=0; vn<pos->second->n_vars(); vn++)
          {
            if( FEInterface::field_type(pos->second->variable_type(vn)) ==
                TYPE_VECTOR )
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
    libmesh_assert_less_equal ( nv, dim*this->n_vars() );

    // Here, we're assuming the number of vector components is the same
    // as the mesh dimension. Will break for mixed dimension meshes.

    var_names.resize( nv );
  }

  // reset
  pos = _systems.begin();

  for (; pos != end; ++pos)
    {
      // Check current system is listed in system_names, and skip pos if not
      bool use_current_system = (system_names == libmesh_nullptr);
      if (!use_current_system)
        use_current_system = system_names->count(pos->first);
      if (!use_current_system)
        continue;

      for (unsigned int vn=0; vn<pos->second->n_vars(); vn++)
        {
          std::string var_name = pos->second->variable_name(vn);
          FEType fe_type = pos->second->variable_type(vn);

          unsigned int n_vec_dim = FEInterface::n_vec_dim( pos->second->get_mesh(), fe_type);

          // Filter on the type if requested
          if (type == libmesh_nullptr || (type && *type == fe_type))
            {
              if( FEInterface::field_type(fe_type) == TYPE_VECTOR )
                {
                  switch(n_vec_dim)
                    {
                    case 0:
                    case 1:
                      var_names[var_num++] = var_name;
                      break;
                    case 2:
                      var_names[var_num++] = var_name+"_x";
                      var_names[var_num++] = var_name+"_y";
                      break;
                    case 3:
                      var_names[var_num++] = var_name+"_x";
                      var_names[var_num++] = var_name+"_y";
                      var_names[var_num++] = var_name+"_z";
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
                                             const std::string &,
                                             const std::string &) const
{
  //TODO:[BSK] re-implement this from the method below
  libmesh_not_implemented();

  //   // Get a reference to the named system
  //   const System & system = this->get_system(system_name);

  //   // Get the number associated with the variable_name we are passed
  //   const unsigned short int variable_num = system.variable_number(variable_name);

  //   // Get the dimension of the current mesh
  //   const unsigned int dim = _mesh.mesh_dimension();

  //   // If we're on processor 0, allocate enough memory to hold the solution.
  //   // Since we're only looking at one variable, there will be one solution value
  //   // for each node in the mesh.
  //   if (_mesh.processor_id() == 0)
  //     soln.resize(_mesh.n_nodes());

  //   // Vector to hold the global solution from all processors
  //   std::vector<Number> sys_soln;

  //   // Update the global solution from all processors
  //   system.update_global_solution (sys_soln, 0);

  //   // Temporary vector to store the solution on an individual element.
  //   std::vector<Number>       elem_soln;

  //   // The FE solution interpolated to the nodes
  //   std::vector<Number>       nodal_soln;

  //   // The DOF indices for the element
  //   std::vector<dof_id_type> dof_indices;

  //   // Determine the finite/infinite element type used in this system
  //   const FEType & fe_type    = system.variable_type(variable_num);

  //   // Define iterators to iterate over all the elements of the mesh
  //   const_active_elem_iterator       it (_mesh.elements_begin());
  //   const const_active_elem_iterator end(_mesh.elements_end());

  //   // Loop over elements
  //   for ( ; it != end; ++it)
  //     {
  //       // Convenient shortcut to the element pointer
  //       const Elem * elem = *it;

  //       // Fill the dof_indices vector for this variable
  //       system.get_dof_map().dof_indices(elem,
  //        dof_indices,
  //        variable_num);

  //       // Resize the element solution vector to fit the
  //       // dof_indices for this element.
  //       elem_soln.resize(dof_indices.size());

  //       // Transfer the system solution to the element
  //       // solution by mapping it through the dof_indices vector.
  //       for (std::size_t i=0; i<dof_indices.size(); i++)
  // elem_soln[i] = sys_soln[dof_indices[i]];

  //       // Using the FE interface, compute the nodal_soln
  //       // for the current elemnt type given the elem_soln
  //       FEInterface::nodal_soln (dim,
  //        fe_type,
  //        elem,
  //        elem_soln,
  //        nodal_soln);

  //       // Sanity check -- make sure that there are the same number
  //       // of entries in the nodal_soln as there are nodes in the
  //       // element!
  //       libmesh_assert_equal_to (nodal_soln.size(), elem->n_nodes());

  //       // Copy the nodal solution over into the correct place in
  //       // the global soln vector which will be returned to the user.
  //       for (unsigned int n=0; n<elem->n_nodes(); n++)
  // soln[elem->node_id(n)] = nodal_soln[n];
  //     }
}




UniquePtr<NumericVector<Number> >
EquationSystems::build_parallel_solution_vector(const std::set<std::string> * system_names) const
{
  LOG_SCOPE("build_parallel_solution_vector()", "EquationSystems");

  // This function must be run on all processors at once
  parallel_object_only();

  const unsigned int dim = _mesh.mesh_dimension();
  const dof_id_type nn   = _mesh.n_nodes();

  // We'd better have a contiguous node numbering
  libmesh_assert_equal_to (nn, _mesh.max_node_id());

  // allocate storage to hold
  // (number_of_nodes)*(number_of_variables) entries.
  // We have to differentiate between between scalar and vector
  // variables. We intercept vector variables and treat each
  // component as a scalar variable (consistently with build_solution_names).

  unsigned int nv = 0;

  //Could this be replaced by a/some convenience methods?[PB]
  {
    unsigned int n_scalar_vars = 0;
    unsigned int n_vector_vars = 0;
    const_system_iterator       pos = _systems.begin();
    const const_system_iterator end = _systems.end();

    for (; pos != end; ++pos)
      {
        // Check current system is listed in system_names, and skip pos if not
        bool use_current_system = (system_names == libmesh_nullptr);
        if (!use_current_system)
          use_current_system = system_names->count(pos->first);
        if (!use_current_system)
          continue;

        for (unsigned int vn=0; vn<pos->second->n_vars(); vn++)
          {
            if( FEInterface::field_type(pos->second->variable_type(vn)) ==
                TYPE_VECTOR )
              n_vector_vars++;
            else
              n_scalar_vars++;
          }
      }
    // Here, we're assuming the number of vector components is the same
    // as the mesh dimension. Will break for mixed dimension meshes.
    nv = n_scalar_vars + dim*n_vector_vars;
  }

  // We handle NULL entries in the Mesh's _nodes vector by just
  // carrying around extra zeroes in the solution vector.
  numeric_index_type parallel_soln_global_size = nn*nv;

  numeric_index_type div = parallel_soln_global_size / this->n_processors();
  numeric_index_type mod = parallel_soln_global_size % this->n_processors();

  // Initialize all processors to the average size.
  numeric_index_type parallel_soln_local_size = div;

  // The first "mod" processors get an extra entry.
  if (this->processor_id() < mod)
    parallel_soln_local_size = div+1;

  // Create a NumericVector to hold the parallel solution
  UniquePtr<NumericVector<Number> > parallel_soln_ptr = NumericVector<Number>::build(_communicator);
  NumericVector<Number> & parallel_soln = *parallel_soln_ptr;
  parallel_soln.init(parallel_soln_global_size, parallel_soln_local_size, /*fast=*/false, PARALLEL);

  // Create a NumericVector to hold the "repeat_count" for each node - this is essentially
  // the number of elements contributing to that node's value
  UniquePtr<NumericVector<Number> > repeat_count_ptr = NumericVector<Number>::build(_communicator);
  NumericVector<Number> & repeat_count = *repeat_count_ptr;
  repeat_count.init(parallel_soln_global_size, parallel_soln_local_size, /*fast=*/false, PARALLEL);

  repeat_count.close();

  unsigned int var_num=0;

  // For each system in this EquationSystems object,
  // update the global solution and if we are on processor 0,
  // loop over the elements and build the nodal solution
  // from the element solution.  Then insert this nodal solution
  // into the vector passed to build_solution_vector.
  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    {
      // Check current system is listed in system_names, and skip pos if not
      bool use_current_system = (system_names == libmesh_nullptr);
      if (!use_current_system)
        use_current_system = system_names->count(pos->first);
      if (!use_current_system)
        continue;

      const System & system  = *(pos->second);
      const unsigned int nv_sys = system.n_vars();
      const unsigned int sys_num = system.number();

      //Could this be replaced by a/some convenience methods?[PB]
      unsigned int n_scalar_vars = 0;
      unsigned int n_vector_vars = 0;
      for (unsigned int vn=0; vn<pos->second->n_vars(); vn++)
        {
          if( FEInterface::field_type(pos->second->variable_type(vn)) ==
              TYPE_VECTOR )
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
        // during during the middle of a solve...  So try to be a bit
        // more careful about calling close() unnecessarily.
        libmesh_assert(this->comm().verify(non_const_sys.solution->closed()));
        if (!non_const_sys.solution->closed())
          non_const_sys.solution->close();
        non_const_sys.update();
      }

      NumericVector<Number> & sys_soln(*system.current_local_solution);

      std::vector<Number>      elem_soln;   // The finite element solution
      std::vector<Number>      nodal_soln;  // The FE solution interpolated to the nodes
      std::vector<dof_id_type> dof_indices; // The DOF indices for the finite element

      for (unsigned int var=0; var<nv_sys; var++)
        {
          const FEType & fe_type           = system.variable_type(var);
          const Variable & var_description = system.variable(var);
          const DofMap & dof_map           = system.get_dof_map();

          unsigned int n_vec_dim = FEInterface::n_vec_dim( pos->second->get_mesh(), fe_type );

          MeshBase::element_iterator       it       = _mesh.active_local_elements_begin();
          const MeshBase::element_iterator end_elem = _mesh.active_local_elements_end();

          for ( ; it != end_elem; ++it)
            {
              const Elem * elem = *it;

              if (var_description.active_on_subdomain((*it)->subdomain_id()))
                {
                  dof_map.dof_indices (elem, dof_indices, var);

                  elem_soln.resize(dof_indices.size());

                  for (std::size_t i=0; i<dof_indices.size(); i++)
                    elem_soln[i] = sys_soln(dof_indices[i]);

                  FEInterface::nodal_soln (dim,
                                           fe_type,
                                           elem,
                                           elem_soln,
                                           nodal_soln);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
                  // infinite elements should be skipped...
                  if (!elem->infinite())
#endif
                    {
                      libmesh_assert_equal_to (nodal_soln.size(), n_vec_dim*elem->n_nodes());

                      for (unsigned int n=0; n<elem->n_nodes(); n++)
                        {
                          for( unsigned int d=0; d < n_vec_dim; d++ )
                            {
                              // For vector-valued elements, all components are in nodal_soln. For each
                              // node, the components are stored in order, i.e. node_0 -> s0_x, s0_y, s0_z
                              parallel_soln.add(nv*(elem->node_id(n)) + (var+d + var_num), nodal_soln[n_vec_dim*n+d]);

                              // Increment the repeat count for this position
                              repeat_count.add(nv*(elem->node_id(n)) + (var+d + var_num), 1);
                            }
                        }
                    }
                }
              else // If this variable doesn't exist on this subdomain we have to still increment repeat_count so that we won't divide by 0 later:
                for (unsigned int n=0; n<elem->n_nodes(); n++)
                  // Only do this if this variable has NO DoFs at this node... it might have some from an ajoining element...
                  if(!elem->node_ptr(n)->n_dofs(sys_num, var))
                    for( unsigned int d=0; d < n_vec_dim; d++ )
                      repeat_count.add(nv*(elem->node_id(n)) + (var+d + var_num), 1);

            } // end loop over elements
        } // end loop on variables in this system

      var_num += nv_sys_split;
    } // end loop over systems

  // Communicate the nodal solution values and repeat counts.
  parallel_soln.close();
  repeat_count.close();

  // If there were gaps in the node numbering, there will be
  // corresponding zeros in the parallel_soln and repeat_count
  // vectors.  We need to set those repeat_count entries to 1
  // in order to avoid dividing by zero.
  for (numeric_index_type i=repeat_count.first_local_index(); i<repeat_count.last_local_index(); ++i)
    {
      // repeat_count entries are integral values but let's avoid a
      // direct floating point comparison with 0 just in case some
      // roundoff noise crept in during vector assembly?
      if (std::abs(repeat_count(i)) < TOLERANCE)
        repeat_count.set(i, 1.);
    }

  // Make sure the repeat_count vector is up-to-date on all processors.
  repeat_count.close();

  // Divide to get the average value at the nodes
  parallel_soln /= repeat_count;

  return UniquePtr<NumericVector<Number> >(parallel_soln_ptr.release());
}



void EquationSystems::build_solution_vector (std::vector<Number> & soln,
                                             const std::set<std::string> * system_names) const
{
  LOG_SCOPE("build_solution_vector()", "EquationSystems");

  // Call the parallel implementation
  UniquePtr<NumericVector<Number> > parallel_soln =
    this->build_parallel_solution_vector(system_names);

  // Localize the NumericVector into the provided std::vector.
  parallel_soln->localize_to_one(soln);
}



void EquationSystems::get_solution (std::vector<Number> & soln,
                                    std::vector<std::string> & names) const
{
  // This function must be run on all processors at once
  parallel_object_only();

  libmesh_assert (this->n_systems());

  const dof_id_type ne = _mesh.n_elem();

  libmesh_assert_equal_to (ne, _mesh.max_elem_id());

  // If the names vector has entries, we will only populate the soln vector
  // with names included in that list.  Note: The names vector may be
  // reordered upon exiting this function
  std::vector<std::string> filter_names = names;
  bool is_filter_names = !filter_names.empty();

  soln.clear();
  names.clear();

  const FEType type(CONSTANT, MONOMIAL);

  dof_id_type nv = 0;

  // Find the total number of variables to output
  std::vector<std::vector<unsigned> > do_output(_systems.size());
  {
    const_system_iterator       pos = _systems.begin();
    const const_system_iterator end = _systems.end();
    unsigned sys_ctr = 0;

    for (; pos != end; ++pos, ++sys_ctr)
      {
        const System & system = *(pos->second);
        const unsigned int nv_sys = system.n_vars();

        do_output[sys_ctr].resize(nv_sys);

        for (unsigned int var=0; var < nv_sys; ++var)
          {
            if (system.variable_type(var) != type ||
                (is_filter_names && std::find(filter_names.begin(), filter_names.end(), system.variable_name(var)) == filter_names.end()))
              continue;

            // Otherwise, this variable should be output
            nv++;
            do_output[sys_ctr][var] = 1;
          }
      }
  }

  // If there are no variables to write out don't do anything...
  if (!nv)
    return;

  // We can handle the case where there are NULLs in the Elem vector
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
  UniquePtr<NumericVector<Number> > parallel_soln_ptr = NumericVector<Number>::build(_communicator);
  NumericVector<Number> & parallel_soln = *parallel_soln_ptr;
  parallel_soln.init(parallel_soln_global_size,
                     parallel_soln_local_size,
                     /*fast=*/false,
                     /*ParallelType=*/PARALLEL);

  dof_id_type var_num = 0;

  // For each system in this EquationSystems object,
  // update the global solution and collect the
  // CONSTANT MONOMIALs.  The entries are in variable-major
  // format.
  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();
  unsigned sys_ctr = 0;

  for (; pos != end; ++pos, ++sys_ctr)
    {
      const System & system  = *(pos->second);
      const unsigned int nv_sys = system.n_vars();

      // Update the current_local_solution
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
      }

      NumericVector<Number> & sys_soln(*system.current_local_solution);

      // The DOF indices for the finite element
      std::vector<dof_id_type> dof_indices;

      // Loop over the variable names and load them in order
      for (unsigned int var=0; var < nv_sys; ++var)
        {
          // Skip this variable if we are not outputting it.
          if (!do_output[sys_ctr][var])
            continue;

          names.push_back(system.variable_name(var));

          const Variable & variable = system.variable(var);
          const DofMap & dof_map = system.get_dof_map();

          MeshBase::element_iterator       it       = _mesh.active_local_elements_begin();
          const MeshBase::element_iterator end_elem = _mesh.active_local_elements_end();

          for ( ; it != end_elem; ++it)
            {
              const Elem * elem = *it;

              if (variable.active_on_subdomain(elem->subdomain_id()))
                {
                  dof_map.dof_indices (elem, dof_indices, var);

                  libmesh_assert_equal_to (1, dof_indices.size());

                  parallel_soln.set((ne*var_num)+elem->id(), sys_soln(dof_indices[0]));
                }
            }

          var_num++;
        } // end loop on variables in this system
    } // end loop over systems

  parallel_soln.close();

  parallel_soln.localize_to_one(soln);
}



void EquationSystems::build_discontinuous_solution_vector (std::vector<Number> & soln,
                                                           const std::set<std::string> * system_names) const
{
  LOG_SCOPE("build_discontinuous_solution_vector()", "EquationSystems");

  libmesh_assert (this->n_systems());

  const unsigned int dim = _mesh.mesh_dimension();

  // Get the number of variables (nv) by counting the number of variables
  // in each system listed in system_names
  unsigned int nv = 0;

  {
    const_system_iterator       pos = _systems.begin();
    const const_system_iterator end = _systems.end();

    for (; pos != end; ++pos)
      {
        // Check current system is listed in system_names, and skip pos if not
        bool use_current_system = (system_names == libmesh_nullptr);
        if (!use_current_system)
          use_current_system = system_names->count(pos->first);
        if (!use_current_system)
          continue;

        const System & system  = *(pos->second);
        nv += system.n_vars();
      }
  }

  unsigned int tw=0;

  // get the total weight
  {
    MeshBase::element_iterator       it  = _mesh.active_elements_begin();
    const MeshBase::element_iterator end = _mesh.active_elements_end();

    for ( ; it != end; ++it)
      tw += (*it)->n_nodes();
  }


  // Only if we are on processor zero, allocate the storage
  // to hold (number_of_nodes)*(number_of_variables) entries.
  if (_mesh.processor_id() == 0)
    soln.resize(tw*nv);

  std::vector<Number> sys_soln;


  unsigned int var_num=0;

  // For each system in this EquationSystems object,
  // update the global solution and if we are on processor 0,
  // loop over the elements and build the nodal solution
  // from the element solution.  Then insert this nodal solution
  // into the vector passed to build_solution_vector.
  {
    const_system_iterator       pos = _systems.begin();
    const const_system_iterator end = _systems.end();

    for (; pos != end; ++pos)
      {
        // Check current system is listed in system_names, and skip pos if not
        bool use_current_system = (system_names == libmesh_nullptr);
        if (!use_current_system)
          use_current_system = system_names->count(pos->first);
        if (!use_current_system)
          continue;

        const System & system  = *(pos->second);
        const unsigned int nv_sys = system.n_vars();

        system.update_global_solution (sys_soln, 0);

        if (_mesh.processor_id() == 0)
          {
            std::vector<Number>       elem_soln;   // The finite element solution
            std::vector<Number>       nodal_soln;  // The FE solution interpolated to the nodes
            std::vector<dof_id_type>  dof_indices; // The DOF indices for the finite element

            for (unsigned int var=0; var<nv_sys; var++)
              {
                const FEType & fe_type    = system.variable_type(var);

                MeshBase::element_iterator       it       = _mesh.active_elements_begin();
                const MeshBase::element_iterator end_elem = _mesh.active_elements_end();

                unsigned int nn=0;

                for ( ; it != end_elem; ++it)
                  {
                    const Elem * elem = *it;
                    system.get_dof_map().dof_indices (elem, dof_indices, var);

                    elem_soln.resize(dof_indices.size());

                    for (std::size_t i=0; i<dof_indices.size(); i++)
                      elem_soln[i] = sys_soln[dof_indices[i]];

                    FEInterface::nodal_soln (dim,
                                             fe_type,
                                             elem,
                                             elem_soln,
                                             nodal_soln);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
                    // infinite elements should be skipped...
                    if (!elem->infinite())
#endif
                      {
                        libmesh_assert_equal_to (nodal_soln.size(), elem->n_nodes());

                        for (unsigned int n=0; n<elem->n_nodes(); n++)
                          {
                            soln[nv*(nn++) + (var + var_num)] +=
                              nodal_soln[n];
                          }
                      }
                  }
              }
          }

        var_num += nv_sys;
      }
  }
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
      const_system_iterator       pos = _systems.begin();
      const const_system_iterator end = _systems.end();

      for (; pos != end; ++pos)
        {
          const std::string & sys_name = pos->first;
          const System &  system        = *(pos->second);

          // get the other system
          const System & other_system   = other_es.get_system (sys_name);

          os_result.push_back (system.compare (other_system, threshold, verbose));

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
  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    oss << pos->second->get_info();


  //   // Possibly print the parameters
  //   if (!this->parameters.empty())
  //     {
  //       oss << "  n_parameters()=" << this->n_parameters() << '\n';
  //       oss << "   Parameters:\n";

  //       for (std::map<std::string, Real>::const_iterator
  //      param = _parameters.begin(); param != _parameters.end();
  //    ++param)
  // oss << "    "
  //     << "\""
  //     << param->first
  //     << "\""
  //     << "="
  //     << param->second
  //     << '\n';
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

  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    tot += pos->second->n_vars();

  return tot;
}



std::size_t EquationSystems::n_dofs () const
{
  std::size_t tot=0;

  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    tot += pos->second->n_dofs();

  return tot;
}




std::size_t EquationSystems::n_active_dofs () const
{
  std::size_t tot=0;

  const_system_iterator       pos = _systems.begin();
  const const_system_iterator end = _systems.end();

  for (; pos != end; ++pos)
    tot += pos->second->n_active_dofs();

  return tot;
}


void EquationSystems::_add_system_to_nodes_and_elems()
{
  // All the nodes
  MeshBase::node_iterator       node_it  = _mesh.nodes_begin();
  const MeshBase::node_iterator node_end = _mesh.nodes_end();

  for ( ; node_it != node_end; ++node_it)
    (*node_it)->add_system();

  // All the elements
  MeshBase::element_iterator       elem_it  = _mesh.elements_begin();
  const MeshBase::element_iterator elem_end = _mesh.elements_end();

  for ( ; elem_it != elem_end; ++elem_it)
    (*elem_it)->add_system();
}

} // namespace libMesh
