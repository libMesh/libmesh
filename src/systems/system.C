// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// Local includes
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/int_range.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/point.h"              // For point_value
#include "libmesh/point_locator_base.h" // For point_value
#include "libmesh/qoi_set.h"
#include "libmesh/enum_to_string.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/system.h"
#include "libmesh/system_norm.h"
#include "libmesh/utility.h"
#include "libmesh/elem.h"
#include "libmesh/fe_type.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_compute_data.h"
#include "libmesh/static_condensation.h"

// includes for calculate_norm, point_*
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/quadrature.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_tools.h"
#include "libmesh/enum_norm_type.h"
#include "libmesh/enum_fe_family.h"

// C++ includes
#include <sstream>   // for std::ostringstream

namespace libMesh
{


// ------------------------------------------------------------
// System implementation
System::System (EquationSystems & es,
                const std::string & name_in,
                const unsigned int number_in) :

  ParallelObject                    (es),
  assemble_before_solve             (true),
  use_fixed_solution                (false),
  extra_quadrature_order            (0),
  solution                          (NumericVector<Number>::build(this->comm())),
  current_local_solution            (NumericVector<Number>::build(this->comm())),
  time                              (0.),
  _init_system_function             (nullptr),
  _init_system_object               (nullptr),
  _assemble_system_function         (nullptr),
  _assemble_system_object           (nullptr),
  _constrain_system_function        (nullptr),
  _constrain_system_object          (nullptr),
  _qoi_evaluate_function            (nullptr),
  _qoi_evaluate_object              (nullptr),
  _qoi_evaluate_derivative_function (nullptr),
  _qoi_evaluate_derivative_object   (nullptr),
  _dof_map                          (std::make_unique<DofMap>(number_in, es.get_mesh())),
  _equation_systems                 (es),
  _mesh                             (es.get_mesh()),
  _sys_name                         (name_in),
  _sys_number                       (number_in),
  _active                           (true),
  _matrices_initialized             (false),
  _solution_projection              (true),
  _basic_system_only                (false),
  _is_initialized                   (false),
  _additional_data_written          (false),
  adjoint_already_solved            (false),
  _hide_output                      (false),
  project_with_constraints          (true),
  _prefer_hash_table_matrix_assembly(false),
  _require_sparsity_pattern         (false),
  _prefix_with_name                 (false)
{
  if (libMesh::on_command_line("--solver-system-names"))
    this->prefix_with_name(true);
  if (libMesh::on_command_line("--" + name_in + "-static-condensation"))
    this->create_static_condensation();
}



System::~System ()
{
  libmesh_exceptionless_assert (!libMesh::closed());
}



dof_id_type System::n_dofs() const
{
  return _dof_map->n_dofs();
}



dof_id_type System::n_constrained_dofs() const
{
#ifdef LIBMESH_ENABLE_CONSTRAINTS

  return _dof_map->n_constrained_dofs();

#else

  return 0;

#endif
}



dof_id_type System::n_local_constrained_dofs() const
{
#ifdef LIBMESH_ENABLE_CONSTRAINTS

  return _dof_map->n_local_constrained_dofs();

#else

  return 0;

#endif
}



dof_id_type System::n_local_dofs() const
{
  return _dof_map->n_local_dofs();
}



Number System::current_solution (const dof_id_type global_dof_number) const
{
  // Check the sizes
  libmesh_assert_less (global_dof_number, _dof_map->n_dofs());
  libmesh_assert_less (global_dof_number, current_local_solution->size());

  return (*current_local_solution)(global_dof_number);
}



void System::clear ()
{
  _dof_map->clear ();
  solution->clear ();
  current_local_solution->clear ();

  // clear any user-added vectors
  _vectors.clear();
  _vector_projections.clear();
  _vector_is_adjoint.clear();
  _is_initialized = false;

  // clear any user-added matrices
  _matrices.clear();
  _matrices_initialized = false;
}



void System::init ()
{
  // Calling init() twice on the same system currently works evil
  // magic, whether done directly or via EquationSystems::read()
  libmesh_assert(!this->is_initialized());

  this->reinit_mesh();
}



void System::init_data ()
{
  parallel_object_only();

  MeshBase & mesh = this->get_mesh();

  // Distribute the degrees of freedom on the mesh
  auto total_dofs = _dof_map->distribute_dofs (mesh);

  // Throw an error if the total number of DOFs is not capable of
  // being indexed by our solution vector.
  auto max_allowed_id = solution->max_allowed_id();
  libmesh_error_msg_if(total_dofs > max_allowed_id,
                       "Cannot allocate a NumericVector with " << total_dofs << " degrees of freedom. "
                       "The vector can only index up to " << max_allowed_id << " entries.");

  // Recreate any user or internal constraints
  this->reinit_constraints();

  // Even if there weren't any constraint changes,
  // reinit_constraints() did prepare_send_list() for us.

  // Now finally after dof distribution and construction of any
  // possible constraints, we may init any static condensation
  // data
  _dof_map->reinit_static_condensation();

  // Resize the solution conformal to the current mesh
  solution->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

  // Resize the current_local_solution for the current mesh
#ifdef LIBMESH_ENABLE_GHOSTED
  current_local_solution->init (this->n_dofs(), this->n_local_dofs(),
                                _dof_map->get_send_list(), /*fast=*/false,
                                GHOSTED);
#else
  current_local_solution->init (this->n_dofs(), false, SERIAL);
#endif

  // from now on, adding additional vectors or variables can't be done
  // without immediately initializing them
  _is_initialized = true;

  // initialize & zero other vectors, if necessary
  for (auto & [vec_name, vec] : _vectors)
    {
      libmesh_ignore(vec_name); // spurious warning from old gcc
      const ParallelType type = vec->type();

      if (type == GHOSTED)
        {
#ifdef LIBMESH_ENABLE_GHOSTED
          vec->init (this->n_dofs(), this->n_local_dofs(),
                           _dof_map->get_send_list(), /*fast=*/false,
                           GHOSTED);
#else
          libmesh_error_msg("Cannot initialize ghosted vectors when they are not enabled.");
#endif
        }
      else if (type == SERIAL)
        {
          vec->init (this->n_dofs(), false, type);
        }
      else
        {
          libmesh_assert_equal_to(type, PARALLEL);
          vec->init (this->n_dofs(), this->n_local_dofs(), false, type);
        }
    }

  // Add matrices
  this->add_matrices();

  // Clear any existing matrices
  for (auto & pr : _matrices)
    pr.second->clear();

  // Initialize the matrices for the system
  if (!_basic_system_only)
   this->init_matrices();
}

void System::reinit_mesh ()
{
  parallel_object_only();

  // First initialize any required data:
  // either only the basic System data
  if (_basic_system_only)
    System::init_data();
  // or all the derived class' data too
  else
    this->init_data();

  // If no variables have been added to this system
  // don't do anything
  if (!this->n_vars())
    return;

  // Then call the user-provided initialization function
  this->user_initialization();

}

void System::init_matrices ()
{
  parallel_object_only();

  // No matrices to init
  if (_matrices.empty())
    {
      // any future matrices to be added will need their own
      // initialization
      _matrices_initialized = true;

      return;
    }

  // Check for quick return in case the first matrix
  // (and by extension all the matrices) has already
  // been initialized
  if (_matrices.begin()->second->initialized())
    {
      libmesh_assert(_matrices_initialized);
      return;
    }

  _matrices_initialized = true;

  // Tell the matrices about the dof map, and vice versa
  for (auto & pr : _matrices)
    {
      SparseMatrix<Number> & m = *(pr.second);
      libmesh_assert (!m.initialized());

      // We want to allow repeated init() on systems, but we don't
      // want to attach the same matrix to the DofMap twice
      if (!this->get_dof_map().is_attached(m))
        this->get_dof_map().attach_matrix(m);

      // If the user has already explicitly requested that this matrix use a hash table, then we
      // always honor that
      const bool use_hash =
          pr.second->use_hash_table() ||
          (this->_prefer_hash_table_matrix_assembly && pr.second->supports_hash_table());
      pr.second->use_hash_table(use_hash);
      // Make this call after we've determined whether the matrix is using a hash table
      if (pr.second->require_sparsity_pattern())
        this->_require_sparsity_pattern = true;
    }

  // Compute the sparsity pattern for the current
  // mesh and DOF distribution.  This also updates
  // additional matrices, \p DofMap now knows them
  if (this->_require_sparsity_pattern)
    this->get_dof_map().compute_sparsity(this->get_mesh());

  // Initialize matrices and set to zero
  for (auto & [name, mat] : _matrices)
    {
      mat->init(_matrix_types[name]);
      mat->zero();
    }
}



void System::restrict_vectors ()
{
  parallel_object_only();

#ifdef LIBMESH_ENABLE_AMR
  // Restrict the _vectors on the coarsened cells
  for (auto & [vec_name, vec] : _vectors)
    {
      NumericVector<Number> * v = vec.get();

      if (_vector_projections[vec_name])
        {
          this->project_vector (*v, this->vector_is_adjoint(vec_name));
        }
      else
        {
          const ParallelType type = vec->type();

          if (type == GHOSTED)
            {
#ifdef LIBMESH_ENABLE_GHOSTED
              vec->init (this->n_dofs(), this->n_local_dofs(),
                               _dof_map->get_send_list(), /*fast=*/false,
                               GHOSTED);
#else
              libmesh_error_msg("Cannot initialize ghosted vectors when they are not enabled.");
#endif
            }
          else
            vec->init (this->n_dofs(), this->n_local_dofs(), false, type);
        }
    }

  const std::vector<dof_id_type> & send_list = _dof_map->get_send_list ();

  // Restrict the solution on the coarsened cells
  if (_solution_projection)
    this->project_vector (*solution);
  // Or at least make sure the solution vector is the correct size
  else
    solution->init (this->n_dofs(), this->n_local_dofs(), true, PARALLEL);

#ifdef LIBMESH_ENABLE_GHOSTED
  current_local_solution->init(this->n_dofs(),
                               this->n_local_dofs(), send_list,
                               false, GHOSTED);
#else
  current_local_solution->init(this->n_dofs());
#endif

  if (_solution_projection)
    solution->localize (*current_local_solution, send_list);

#endif // LIBMESH_ENABLE_AMR
}



void System::prolong_vectors ()
{
#ifdef LIBMESH_ENABLE_AMR
  // Currently project_vector handles both restriction and prolongation
  this->restrict_vectors();
#endif
}



void System::reinit ()
{
  parallel_object_only();

  // project_vector handles vector initialization now
  libmesh_assert_equal_to (solution->size(), current_local_solution->size());

  // Make sure our static condensation dof map is up-to-date before we init any
  // static condensation matrices
  this->get_dof_map().reinit_static_condensation();

  if (!_matrices.empty() && !_basic_system_only)
    {
      // Clear the matrices
      for (auto & pr : _matrices)
        {
          pr.second->clear();
          pr.second->attach_dof_map(this->get_dof_map());
        }

      if (this->_require_sparsity_pattern)
        {
          // Clear the sparsity pattern
          this->get_dof_map().clear_sparsity();

          // Compute the sparsity pattern for the current
          // mesh and DOF distribution.  This also updates
          // additional matrices, \p DofMap now knows them
          this->get_dof_map().compute_sparsity (this->get_mesh());
        }

      // Initialize matrices and set to zero
      for (auto & pr : _matrices)
        {
          pr.second->init();
          pr.second->zero();
        }
    }
}


void System::reinit_constraints()
{
  parallel_object_only();

#ifdef LIBMESH_ENABLE_CONSTRAINTS
  get_dof_map().create_dof_constraints(_mesh, this->time);
  user_constrain();
  get_dof_map().process_constraints(_mesh);
  if (libMesh::on_command_line ("--print-constraints"))
    get_dof_map().print_dof_constraints(libMesh::out);
#endif
  get_dof_map().prepare_send_list();
}


void System::update ()
{
  parallel_object_only();

  libmesh_assert(solution->closed());

  const std::vector<dof_id_type> & send_list = _dof_map->get_send_list ();

  // Check sizes
  libmesh_assert_equal_to (current_local_solution->size(), solution->size());
  // More processors than elements => empty send_list
  //  libmesh_assert (!send_list.empty());
  libmesh_assert_less_equal (send_list.size(), solution->size());

  // Create current_local_solution from solution.  This will
  // put a local copy of solution into current_local_solution.
  // Only the necessary values (specified by the send_list)
  // are copied to minimize communication
  solution->localize (*current_local_solution, send_list);
}



void System::re_update ()
{
  parallel_object_only();

  // If this system is empty... don't do anything!
  if (!this->n_vars())
    return;

  const std::vector<dof_id_type> & send_list = this->get_dof_map().get_send_list ();

  // Check sizes
  libmesh_assert_equal_to (current_local_solution->size(), solution->size());
  // Not true with ghosted vectors
  // libmesh_assert_equal_to (current_local_solution->local_size(), solution->size());
  // libmesh_assert (!send_list.empty());
  libmesh_assert_less_equal (send_list.size(), solution->size());

  // Create current_local_solution from solution.  This will
  // put a local copy of solution into current_local_solution.
  solution->localize (*current_local_solution, send_list);
}



void System::restrict_solve_to (const SystemSubset * subset,
                                const SubsetSolveMode /*subset_solve_mode*/)
{
  if (subset != nullptr)
    libmesh_not_implemented();
}



void System::assemble ()
{
  // Log how long the user's assembly code takes
  LOG_SCOPE("assemble()", "System");

  libmesh_assert(this->get_mesh().is_prepared());
#ifdef DEBUG
  MeshTools::libmesh_assert_valid_is_prepared(this->get_mesh());
#endif

  // Call the user-specified assembly function
  this->user_assembly();
}



void System::assemble_qoi (const QoISet & qoi_indices)
{
  // Log how long the user's assembly code takes
  LOG_SCOPE("assemble_qoi()", "System");

  libmesh_assert(this->get_mesh().is_prepared());
#ifdef DEBUG
  MeshTools::libmesh_assert_valid_is_prepared(this->get_mesh());
#endif

  // Call the user-specified quantity of interest function
  this->user_QOI(qoi_indices);
}



void System::assemble_qoi_derivative(const QoISet & qoi_indices,
                                     bool include_liftfunc,
                                     bool apply_constraints)
{
  // Log how long the user's assembly code takes
  LOG_SCOPE("assemble_qoi_derivative()", "System");

  libmesh_assert(this->get_mesh().is_prepared());
#ifdef DEBUG
  MeshTools::libmesh_assert_valid_is_prepared(this->get_mesh());
#endif

  // Call the user-specified quantity of interest function
  this->user_QOI_derivative(qoi_indices, include_liftfunc,
                            apply_constraints);
}



void System::qoi_parameter_sensitivity (const QoISet & qoi_indices,
                                        const ParameterVector & parameters_vec,
                                        SensitivityData & sensitivities)
{
  // Forward sensitivities are more efficient for Nq > Np
  if (qoi_indices.size(*this) > parameters_vec.size())
    forward_qoi_parameter_sensitivity(qoi_indices, parameters_vec, sensitivities);
  // Adjoint sensitivities are more efficient for Np > Nq,
  // and an adjoint may be more reusable than a forward
  // solution sensitivity in the Np == Nq case.
  else
    adjoint_qoi_parameter_sensitivity(qoi_indices, parameters_vec, sensitivities);
}



bool System::compare (const System & other_system,
                      const Real threshold,
                      const bool verbose) const
{
  // we do not care for matrices, but for vectors
  libmesh_assert (_is_initialized);
  libmesh_assert (other_system._is_initialized);

  if (verbose)
    {
      libMesh::out << "  Systems \"" << _sys_name << "\"" << std::endl;
      libMesh::out << "   comparing matrices not supported." << std::endl;
      libMesh::out << "   comparing names...";
    }

  // compare the name: 0 means identical
  const int name_result = _sys_name.compare(other_system.name());
  if (verbose)
    {
      if (name_result == 0)
        libMesh::out << " identical." << std::endl;
      else
        libMesh::out << "  names not identical." << std::endl;
      libMesh::out << "   comparing solution vector...";
    }


  // compare the solution: -1 means identical
  const int solu_result = solution->compare (*other_system.solution.get(),
                                             threshold);

  if (verbose)
    {
      if (solu_result == -1)
        libMesh::out << " identical up to threshold." << std::endl;
      else
        libMesh::out << "  first difference occurred at index = "
                     << solu_result << "." << std::endl;
    }


  // safety check, whether we handle at least the same number
  // of vectors
  std::vector<int> ov_result;

  if (this->n_vectors() != other_system.n_vectors())
    {
      if (verbose)
        {
          libMesh::out << "   Fatal difference. This system handles "
                       << this->n_vectors() << " add'l vectors," << std::endl
                       << "   while the other system handles "
                       << other_system.n_vectors()
                       << " add'l vectors." << std::endl
                       << "   Aborting comparison." << std::endl;
        }
      return false;
    }
  else if (this->n_vectors() == 0)
    {
      // there are no additional vectors...
      ov_result.clear ();
    }
  else
    {
      // compare other vectors
      for (auto & [vec_name, vec] : _vectors)
        {
          if (verbose)
            libMesh::out << "   comparing vector \""
                         << vec_name << "\" ...";

          // assume they have the same name
          const NumericVector<Number> & other_system_vector =
            other_system.get_vector(vec_name);

          ov_result.push_back(vec->compare(other_system_vector, threshold));

          if (verbose)
            {
              if (ov_result[ov_result.size()-1] == -1)
                libMesh::out << " identical up to threshold." << std::endl;
              else
                libMesh::out << " first difference occurred at" << std::endl
                             << "   index = " << ov_result[ov_result.size()-1] << "." << std::endl;
            }
        }
    } // finished comparing additional vectors


  bool overall_result;

  // sum up the results
  if ((name_result==0) && (solu_result==-1))
    {
      if (ov_result.size()==0)
        overall_result = true;
      else
        {
          bool ov_identical;
          unsigned int n    = 0;
          do
            {
              ov_identical = (ov_result[n]==-1);
              n++;
            }
          while (ov_identical && n<ov_result.size());
          overall_result = ov_identical;
        }
    }
  else
    overall_result = false;

  if (verbose)
    {
      libMesh::out << "   finished comparisons, ";
      if (overall_result)
        libMesh::out << "found no differences." << std::endl << std::endl;
      else
        libMesh::out << "found differences." << std::endl << std::endl;
    }

  return overall_result;
}



void System::update_global_solution (std::vector<Number> & global_soln) const
{
  parallel_object_only();

  global_soln.resize (solution->size());

  solution->localize (global_soln);
}



void System::update_global_solution (std::vector<Number> & global_soln,
                                     const processor_id_type dest_proc) const
{
  parallel_object_only();

  global_soln.resize        (solution->size());

  solution->localize_to_one (global_soln, dest_proc);
}



NumericVector<Number> & System::add_vector (std::string_view vec_name,
                                            const bool projections,
                                            const ParallelType type)
{
  parallel_object_only();

  libmesh_assert(this->comm().verify(std::string(vec_name)));
  libmesh_assert(this->comm().verify(int(type)));
  libmesh_assert(this->comm().verify(projections));

  // Return the vector if it is already there.
  if (auto it = this->_vectors.find(vec_name);
      it != this->_vectors.end())
    {
      // If the projection setting has *upgraded*, change it.
      if (projections) // only do expensive lookup if needed
        libmesh_map_find(_vector_projections, vec_name) = projections;

      NumericVector<Number> & vec = *it->second;

      // If we're in serial, our vectors are effectively SERIAL, so
      // we'll ignore any type setting.  If we're in parallel, we
      // might have a type change to deal with.

      if (this->n_processors() > 1)
        {
          // If the type setting has changed in a way we can't
          // perceive as an upgrade or a downgrade, scream.
          libmesh_assert_equal_to(type == SERIAL,
                                  vec.type() == SERIAL);

          // If the type setting has *upgraded*, change it.
          if (type == GHOSTED && vec.type() == PARALLEL)
            {
              // A *really* late upgrade is expensive, but better not
              // to risk zeroing data.
              if (vec.initialized())
                {
                  if (!vec.closed())
                    vec.close();

                  // Ideally we'd move parallel coefficients and then
                  // add ghosted coefficients, but copy and swap is
                  // simpler.  If anyone actually ever uses this case
                  // for real we can look into optimizing it.
                  auto new_vec = NumericVector<Number>::build(this->comm());
#ifdef LIBMESH_ENABLE_GHOSTED
                  new_vec->init (this->n_dofs(), this->n_local_dofs(),
                                 _dof_map->get_send_list(), /*fast=*/false,
                                 GHOSTED);
#else
                  libmesh_error_msg("Cannot initialize ghosted vectors when they are not enabled.");
#endif

                  *new_vec = vec;
                  vec.swap(*new_vec);
                }
              else
                // The PARALLEL vec is not yet initialized, so we can
                // just "upgrade" it to GHOSTED.
                vec.set_type(type);
            }
        }

      // Any upgrades are done; we're happy here.
      return vec;
    }

  // Otherwise, build the vector. The following emplace() is
  // guaranteed to succeed because, if we made it here, we don't
  // already have a vector named "vec_name". We pass the user's
  // requested ParallelType directly to NumericVector::build() so
  // that, even if the vector is not initialized now, it will get the
  // right type when it is initialized later.
  auto pr =
    _vectors.emplace(vec_name,
                     NumericVector<Number>::build(this->comm(),
                                                  libMesh::default_solver_package(),
                                                  type));
  auto buf = pr.first->second.get();
  _vector_projections.emplace(vec_name, projections);

  // Vectors are primal by default
  _vector_is_adjoint.emplace(vec_name, -1);

  // Initialize it if necessary
  if (_is_initialized)
    {
      if (type == GHOSTED)
        {
#ifdef LIBMESH_ENABLE_GHOSTED
          buf->init (this->n_dofs(), this->n_local_dofs(),
                     _dof_map->get_send_list(), /*fast=*/false,
                     GHOSTED);
#else
          libmesh_error_msg("Cannot initialize ghosted vectors when they are not enabled.");
#endif
        }
      else
        buf->init (this->n_dofs(), this->n_local_dofs(), false, type);
    }

  return *buf;
}

void System::remove_vector (std::string_view vec_name)
{
  parallel_object_only();  // Not strictly needed, but the only safe way to keep in sync

  if (const auto pos = _vectors.find(vec_name);
      pos != _vectors.end())
    {
      _vectors.erase(pos);
      auto proj_it = _vector_projections.find(vec_name);
      libmesh_assert(proj_it != _vector_projections.end());
      _vector_projections.erase(proj_it);

      auto adj_it = _vector_is_adjoint.find(vec_name);
      libmesh_assert(adj_it != _vector_is_adjoint.end());
      _vector_is_adjoint.erase(adj_it);
    }
}

const NumericVector<Number> * System::request_vector (std::string_view vec_name) const
{
  if (const auto pos = _vectors.find(vec_name);
      pos != _vectors.end())
    return pos->second.get();

  // Otherwise, vec_name was not found
  return nullptr;
}



NumericVector<Number> * System::request_vector (std::string_view vec_name)
{
  if (auto pos = _vectors.find(vec_name);
      pos != _vectors.end())
    return pos->second.get();

  // Otherwise, vec_name was not found
  return nullptr;
}



const NumericVector<Number> * System::request_vector (const unsigned int vec_num) const
{
  // If we don't have that many vectors, return nullptr
  if (vec_num >= _vectors.size())
    return nullptr;

  // Otherwise return a pointer to the vec_num'th vector
  auto it = vectors_begin();
  std::advance(it, vec_num);
  return it->second.get();
}



NumericVector<Number> * System::request_vector (const unsigned int vec_num)
{
  // If we don't have that many vectors, return nullptr
  if (vec_num >= _vectors.size())
    return nullptr;

  // Otherwise return a pointer to the vec_num'th vector
  auto it = vectors_begin();
  std::advance(it, vec_num);
  return it->second.get();
}



const NumericVector<Number> & System::get_vector (std::string_view vec_name) const
{
  return *(libmesh_map_find(_vectors, vec_name));
}



NumericVector<Number> & System::get_vector (std::string_view vec_name)
{
  return *(libmesh_map_find(_vectors, vec_name));
}



const NumericVector<Number> & System::get_vector (const unsigned int vec_num) const
{
  // If we don't have that many vectors, throw an error
  libmesh_assert_less(vec_num, _vectors.size());

  // Otherwise return a reference to the vec_num'th vector
  auto it = vectors_begin();
  std::advance(it, vec_num);
  return *(it->second);
}



NumericVector<Number> & System::get_vector (const unsigned int vec_num)
{
  // If we don't have that many vectors, throw an error
  libmesh_assert_less(vec_num, _vectors.size());

  // Otherwise return a reference to the vec_num'th vector
  auto it = vectors_begin();
  std::advance(it, vec_num);
  return *(it->second);
}



const std::string & System::vector_name (const unsigned int vec_num) const
{
  // If we don't have that many vectors, throw an error
  libmesh_assert_less(vec_num, _vectors.size());

  // Otherwise return a reference to the vec_num'th vector name
  auto it = vectors_begin();
  std::advance(it, vec_num);
  return it->first;
}

const std::string & System::vector_name (const NumericVector<Number> & vec_reference) const
{
  // Linear search for a vector whose pointer matches vec_reference
  auto it = std::find_if(vectors_begin(), vectors_end(),
                         [&vec_reference](const decltype(_vectors)::value_type & pr)
                         { return &vec_reference == pr.second.get(); });

  // Before returning, make sure we didn't loop till the end and not find any match
  libmesh_assert (it != vectors_end());

  // Return the string associated with the current vector
  return it->first;
}



SparseMatrix<Number> & System::add_matrix (std::string_view mat_name,
                                           const ParallelType type,
                                           const MatrixBuildType mat_build_type)
{
  parallel_object_only();

  libmesh_assert(this->comm().verify(std::string(mat_name)));
  libmesh_assert(this->comm().verify(int(type)));
  libmesh_assert(this->comm().verify(int(mat_build_type)));

  // Return the matrix if it is already there.
  if (auto it = this->_matrices.find(mat_name);
      it != this->_matrices.end())
    return *it->second;

  // Otherwise build the matrix to return.
  std::unique_ptr<SparseMatrix<Number>> matrix;
  if (this->has_static_condensation())
    {
      if (mat_build_type == MatrixBuildType::DIAGONAL)
        libmesh_error_msg(
            "We do not currently support static condensation of the diagonal matrix type");
      matrix = std::make_unique<StaticCondensation>(this->get_mesh(),
                                                    *this,
                                                    this->get_dof_map(),
                                                    this->get_dof_map().get_static_condensation());
    }
  else
    matrix = SparseMatrix<Number>::build(this->comm(), libMesh::default_solver_package());
  auto & mat = *matrix;

  _matrices.emplace(mat_name, std::move(matrix));

  _matrix_types.emplace(mat_name, type);

  // Initialize it first if we've already initialized the others.
  this->late_matrix_init(mat, type);

  return mat;
}



SparseMatrix<Number> & System::add_matrix (std::string_view mat_name,
                                           std::unique_ptr<SparseMatrix<Number>> matrix,
                                           const ParallelType type)
{
  parallel_object_only();

  libmesh_assert(this->comm().verify(std::string(mat_name)));
  libmesh_assert(this->comm().verify(int(type)));

  auto [it, inserted] = _matrices.emplace(mat_name, std::move(matrix));
  libmesh_error_msg_if(!inserted,
                       "Tried to add '" << mat_name << "' but the matrix already exists");

  _matrix_types.emplace(mat_name, type);

  SparseMatrix<Number> & mat = *(it->second);

  // Initialize it first if we've already initialized the others.
  this->late_matrix_init(mat, type);

  return mat;
}

void System::late_matrix_init(SparseMatrix<Number> & mat,
                              ParallelType type)
{
  if (_matrices_initialized)
    {
      this->get_dof_map().attach_matrix(mat);
      mat.init(type);
    }
}




void System::remove_matrix (std::string_view mat_name)
{
  parallel_object_only();  // Not strictly needed, but the only safe way to keep in sync

  if (const auto pos = _matrices.find(mat_name);
      pos != _matrices.end())
    _matrices.erase(pos); // erase()'d entries are destroyed
}



const SparseMatrix<Number> * System::request_matrix (std::string_view mat_name) const
{
  if (const auto pos = _matrices.find(mat_name);
      pos != _matrices.end())
    return pos->second.get();

  // Otherwise, mat_name does not exist
  return nullptr;
}



SparseMatrix<Number> * System::request_matrix (std::string_view mat_name)
{
  if (auto pos = _matrices.find(mat_name);
      pos != _matrices.end())
    return pos->second.get();

  // Otherwise, mat_name does not exist
  return nullptr;
}



const SparseMatrix<Number> & System::get_matrix (std::string_view mat_name) const
{
  return *libmesh_map_find(_matrices, mat_name);
}



SparseMatrix<Number> & System::get_matrix (std::string_view mat_name)
{
  return *libmesh_map_find(_matrices, mat_name);
}



void System::set_vector_preservation (const std::string & vec_name,
                                      bool preserve)
{
  parallel_object_only();  // Not strictly needed, but the only safe way to keep in sync

  _vector_projections[vec_name] = preserve;
}



bool System::vector_preservation (std::string_view vec_name) const
{
  if (auto it = _vector_projections.find(vec_name);
      it != _vector_projections.end())
    return it->second;

  // vec_name was not in the map, return false
  return false;
}



void System::set_vector_as_adjoint (const std::string & vec_name,
                                    int qoi_num)
{
  parallel_object_only();  // Not strictly needed, but the only safe way to keep in sync

  // We reserve -1 for vectors which get primal constraints, -2 for
  // vectors which get no constraints
  libmesh_assert_greater_equal(qoi_num, -2);
  _vector_is_adjoint[vec_name] = qoi_num;
}



int System::vector_is_adjoint (std::string_view vec_name) const
{
  const auto it = _vector_is_adjoint.find(vec_name);
  libmesh_assert(it != _vector_is_adjoint.end());
  return it->second;
}



NumericVector<Number> & System::add_sensitivity_solution (unsigned int i)
{
  std::ostringstream sensitivity_name;
  sensitivity_name << "sensitivity_solution" << i;

  return this->add_vector(sensitivity_name.str());
}



NumericVector<Number> & System::get_sensitivity_solution (unsigned int i)
{
  std::ostringstream sensitivity_name;
  sensitivity_name << "sensitivity_solution" << i;

  return this->get_vector(sensitivity_name.str());
}



const NumericVector<Number> & System::get_sensitivity_solution (unsigned int i) const
{
  std::ostringstream sensitivity_name;
  sensitivity_name << "sensitivity_solution" << i;

  return this->get_vector(sensitivity_name.str());
}



NumericVector<Number> & System::add_weighted_sensitivity_solution ()
{
  return this->add_vector("weighted_sensitivity_solution");
}



NumericVector<Number> & System::get_weighted_sensitivity_solution ()
{
  return this->get_vector("weighted_sensitivity_solution");
}



const NumericVector<Number> & System::get_weighted_sensitivity_solution () const
{
  return this->get_vector("weighted_sensitivity_solution");
}



NumericVector<Number> & System::add_adjoint_solution (unsigned int i)
{
  std::ostringstream adjoint_name;
  adjoint_name << "adjoint_solution" << i;

  NumericVector<Number> & returnval = this->add_vector(adjoint_name.str());
  this->set_vector_as_adjoint(adjoint_name.str(), i);
  return returnval;
}



NumericVector<Number> & System::get_adjoint_solution (unsigned int i)
{
  std::ostringstream adjoint_name;
  adjoint_name << "adjoint_solution" << i;

  return this->get_vector(adjoint_name.str());
}



const NumericVector<Number> & System::get_adjoint_solution (unsigned int i) const
{
  std::ostringstream adjoint_name;
  adjoint_name << "adjoint_solution" << i;

  return this->get_vector(adjoint_name.str());
}



NumericVector<Number> & System::add_weighted_sensitivity_adjoint_solution (unsigned int i)
{
  std::ostringstream adjoint_name;
  adjoint_name << "weighted_sensitivity_adjoint_solution" << i;

  NumericVector<Number> & returnval = this->add_vector(adjoint_name.str());
  this->set_vector_as_adjoint(adjoint_name.str(), i);
  return returnval;
}



NumericVector<Number> & System::get_weighted_sensitivity_adjoint_solution (unsigned int i)
{
  std::ostringstream adjoint_name;
  adjoint_name << "weighted_sensitivity_adjoint_solution" << i;

  return this->get_vector(adjoint_name.str());
}



const NumericVector<Number> & System::get_weighted_sensitivity_adjoint_solution (unsigned int i) const
{
  std::ostringstream adjoint_name;
  adjoint_name << "weighted_sensitivity_adjoint_solution" << i;

  return this->get_vector(adjoint_name.str());
}



NumericVector<Number> & System::add_adjoint_rhs (unsigned int i)
{
  std::ostringstream adjoint_rhs_name;
  adjoint_rhs_name << "adjoint_rhs" << i;

  return this->add_vector(adjoint_rhs_name.str(), false);
}



NumericVector<Number> & System::get_adjoint_rhs (unsigned int i)
{
  std::ostringstream adjoint_rhs_name;
  adjoint_rhs_name << "adjoint_rhs" << i;

  return this->get_vector(adjoint_rhs_name.str());
}



const NumericVector<Number> & System::get_adjoint_rhs (unsigned int i) const
{
  std::ostringstream adjoint_rhs_name;
  adjoint_rhs_name << "adjoint_rhs" << i;

  return this->get_vector(adjoint_rhs_name.str());
}



NumericVector<Number> & System::add_sensitivity_rhs (unsigned int i)
{
  std::ostringstream sensitivity_rhs_name;
  sensitivity_rhs_name << "sensitivity_rhs" << i;

  return this->add_vector(sensitivity_rhs_name.str(), false);
}



NumericVector<Number> & System::get_sensitivity_rhs (unsigned int i)
{
  std::ostringstream sensitivity_rhs_name;
  sensitivity_rhs_name << "sensitivity_rhs" << i;

  return this->get_vector(sensitivity_rhs_name.str());
}



const NumericVector<Number> & System::get_sensitivity_rhs (unsigned int i) const
{
  std::ostringstream sensitivity_rhs_name;
  sensitivity_rhs_name << "sensitivity_rhs" << i;

  return this->get_vector(sensitivity_rhs_name.str());
}



unsigned int System::add_variable (std::string_view var,
                                   const FEType & type,
                                   const std::set<subdomain_id_type> * const active_subdomains)
{
  return this->get_dof_map().add_variable(*this, var, type, active_subdomains);
}



unsigned int System::add_variable (std::string_view var,
                                   const Order order,
                                   const FEFamily family,
                                   const std::set<subdomain_id_type> * const active_subdomains)
{
  return this->add_variable(var,
                            FEType(order, family),
                            active_subdomains);
}



unsigned int System::add_variables (const std::vector<std::string> & vars,
                                    const FEType & type,
                                    const std::set<subdomain_id_type> * const active_subdomains)
{
  return this->get_dof_map().add_variables(*this, vars, type, active_subdomains);
}



unsigned int System::add_variables (const std::vector<std::string> & vars,
                                    const Order order,
                                    const FEFamily family,
                                    const std::set<subdomain_id_type> * const active_subdomains)
{
  return this->add_variables(vars,
                             FEType(order, family),
                             active_subdomains);
}

unsigned int System::add_variable_array (const std::vector<std::string> & vars,
                                         const FEType & type,
                                         const std::set<subdomain_id_type> * const active_subdomains)
{
  return this->get_dof_map().add_variable_array(*this, vars, type, active_subdomains);
}

bool System::has_variable (std::string_view var) const
{
  return this->get_dof_map().has_variable(var);
}

unsigned int System::variable_number (std::string_view var) const
{
  return this->get_dof_map().variable_number(var);
}

void System::get_all_variable_numbers(std::vector<unsigned int> & all_variable_numbers) const
{
  this->get_dof_map().get_all_variable_numbers(all_variable_numbers);
}


void System::local_dof_indices(const unsigned int var,
                               std::set<dof_id_type> & var_indices) const
{
  // Make sure the set is clear
  var_indices.clear();

  std::vector<dof_id_type> dof_indices;

  const dof_id_type
    first_local = this->get_dof_map().first_dof(),
    end_local   = this->get_dof_map().end_dof();

  // Begin the loop over the elements
  for (const auto & elem : this->get_mesh().active_local_element_ptr_range())
    {
      this->get_dof_map().dof_indices (elem, dof_indices, var);

      for (dof_id_type dof : dof_indices)
        //If the dof is owned by the local processor
        if (first_local <= dof && dof < end_local)
          var_indices.insert(dof);
    }

  // we may have missed assigning DOFs to nodes that we own
  // but to which we have no connected elements matching our
  // variable restriction criterion.  this will happen, for example,
  // if variable V is restricted to subdomain S.  We may not own
  // any elements which live in S, but we may own nodes which are
  // *connected* to elements which do.
  for (const auto & node : this->get_mesh().local_node_ptr_range())
    {
      libmesh_assert(node);
      this->get_dof_map().dof_indices (node, dof_indices, var);
      for (auto dof : dof_indices)
        if (first_local <= dof && dof < end_local)
          var_indices.insert(dof);
    }
}



void System::zero_variable (NumericVector<Number> & v,
                            unsigned int var_num) const
{
  /* Make sure the call makes sense.  */
  libmesh_assert_less (var_num, this->n_vars());

  /* Get a reference to the mesh.  */
  const MeshBase & mesh = this->get_mesh();

  /* Check which system we are.  */
  const unsigned int sys_num = this->number();

  // Loop over nodes.
  for (const auto & node : mesh.local_node_ptr_range())
    {
      unsigned int n_comp = node->n_comp(sys_num,var_num);
      for (unsigned int i=0; i<n_comp; i++)
        {
          const dof_id_type index = node->dof_number(sys_num,var_num,i);
          v.set(index,0.0);
        }
    }

  // Loop over elements.
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      unsigned int n_comp = elem->n_comp(sys_num,var_num);
      for (unsigned int i=0; i<n_comp; i++)
        {
          const dof_id_type index = elem->dof_number(sys_num,var_num,i);
          v.set(index,0.0);
        }
    }
}



Real System::discrete_var_norm(const NumericVector<Number> & v,
                               unsigned int var,
                               FEMNormType norm_type) const
{
  std::set<dof_id_type> var_indices;
  local_dof_indices(var, var_indices);

  if (norm_type == DISCRETE_L1)
    return v.subset_l1_norm(var_indices);
  if (norm_type == DISCRETE_L2)
    return v.subset_l2_norm(var_indices);
  if (norm_type == DISCRETE_L_INF)
    return v.subset_linfty_norm(var_indices);
  else
    libmesh_error_msg("Invalid norm_type = " << Utility::enum_to_string(norm_type));
}



Real System::calculate_norm(const NumericVector<Number> & v,
                            unsigned int var,
                            FEMNormType norm_type,
                            std::set<unsigned int> * skip_dimensions) const
{
  //short circuit to save time
  if (norm_type == DISCRETE_L1 ||
      norm_type == DISCRETE_L2 ||
      norm_type == DISCRETE_L_INF)
    return discrete_var_norm(v,var,norm_type);

  // Not a discrete norm
  std::vector<FEMNormType> norms(this->n_vars(), L2);
  std::vector<Real> weights(this->n_vars(), 0.0);
  norms[var] = norm_type;
  weights[var] = 1.0;
  Real val = this->calculate_norm(v, SystemNorm(norms, weights), skip_dimensions);
  return val;
}



Real System::calculate_norm(const NumericVector<Number> & v,
                            const SystemNorm & norm,
                            std::set<unsigned int> * skip_dimensions) const
{
  // This function must be run on all processors at once
  parallel_object_only();

  LOG_SCOPE ("calculate_norm()", "System");

  // Zero the norm before summation
  Real v_norm = 0.;

  if (norm.is_discrete())
    {
      //Check to see if all weights are 1.0 and all types are equal
      FEMNormType norm_type0 = norm.type(0);
      unsigned int check_var = 0, check_end = this->n_vars();
      for (; check_var != check_end; ++check_var)
        if ((norm.weight(check_var) != 1.0) || (norm.type(check_var) != norm_type0))
          break;

      //All weights were 1.0 so just do the full vector discrete norm
      if (check_var == this->n_vars())
        {
          if (norm_type0 == DISCRETE_L1)
            return v.l1_norm();
          if (norm_type0 == DISCRETE_L2)
            return v.l2_norm();
          if (norm_type0 == DISCRETE_L_INF)
            return v.linfty_norm();
          else
            libmesh_error_msg("Invalid norm_type0 = " << Utility::enum_to_string(norm_type0));
        }

      for (auto var : make_range(this->n_vars()))
        {
          // Skip any variables we don't need to integrate
          if (norm.weight(var) == 0.0)
            continue;

          v_norm += norm.weight(var) * discrete_var_norm(v, var, norm.type(var));
        }

      return v_norm;
    }

  // Localize the potentially parallel vector
  std::unique_ptr<NumericVector<Number>> local_v = NumericVector<Number>::build(this->comm());
  local_v->init(v.size(), v.local_size(), _dof_map->get_send_list(),
                true, GHOSTED);
  v.localize (*local_v, _dof_map->get_send_list());

  // I'm not sure how best to mix Hilbert norms on some variables (for
  // which we'll want to square then sum then square root) with norms
  // like L_inf (for which we'll just want to take an absolute value
  // and then sum).
  bool using_hilbert_norm = true,
    using_nonhilbert_norm = true;

  // Loop over all variables
  for (auto var : make_range(this->n_vars()))
    {
      // Skip any variables we don't need to integrate
      Real norm_weight_sq = norm.weight_sq(var);
      if (norm_weight_sq == 0.0)
        continue;
      Real norm_weight = norm.weight(var);

      // Check for unimplemented norms (rather than just returning 0).
      FEMNormType norm_type = norm.type(var);
      if ((norm_type==H1) ||
          (norm_type==H2) ||
          (norm_type==L2) ||
          (norm_type==H1_SEMINORM) ||
          (norm_type==H2_SEMINORM))
        {
          if (!using_hilbert_norm)
            libmesh_not_implemented();
          using_nonhilbert_norm = false;
        }
      else if ((norm_type==L1) ||
               (norm_type==L_INF) ||
               (norm_type==W1_INF_SEMINORM) ||
               (norm_type==W2_INF_SEMINORM))
        {
          if (!using_nonhilbert_norm)
            libmesh_not_implemented();
          using_hilbert_norm = false;
        }
      else
        libmesh_not_implemented();

      const FEType & fe_type = this->get_dof_map().variable_type(var);

      // Allow space for dims 0-3, and for both scalar and vector
      // elements, even if we don't use them all
      std::vector<std::unique_ptr<FEBase>> fe_ptrs(4);
      std::vector<std::unique_ptr<FEVectorBase>> vec_fe_ptrs(4);
      std::vector<std::unique_ptr<QBase>> q_rules(4);

      const std::set<unsigned char> & elem_dims = _mesh.elem_dimensions();

      // Prepare finite elements for each dimension present in the mesh
      for (const auto & dim : elem_dims)
        {
          if (skip_dimensions && skip_dimensions->find(dim) != skip_dimensions->end())
            continue;

          // Construct quadrature and finite element objects
          q_rules[dim] = fe_type.default_quadrature_rule (dim);

          const FEFieldType field_type = FEInterface::field_type(fe_type);
          if (field_type == TYPE_SCALAR)
            {
              fe_ptrs[dim] = FEBase::build(dim, fe_type);
              fe_ptrs[dim]->attach_quadrature_rule (q_rules[dim].get());
            }
          else
            {
              vec_fe_ptrs[dim] = FEVectorBase::build(dim, fe_type);
              vec_fe_ptrs[dim]->attach_quadrature_rule (q_rules[dim].get());
              libmesh_assert_equal_to(field_type, TYPE_VECTOR);
            }

        }

      std::vector<dof_id_type> dof_indices;

      // Begin the loop over the elements
      for (const auto & elem : this->get_mesh().active_local_element_ptr_range())
        {
          const unsigned int dim = elem->dim();

          // One way for implementing this would be to exchange the fe with the FEInterface- class.
          // However, it needs to be discussed whether integral-norms make sense for infinite elements.
          // or in which sense they could make sense.
          if (elem->infinite() )
            libmesh_not_implemented();

          if (skip_dimensions && skip_dimensions->find(dim) != skip_dimensions->end())
            continue;

          QBase * qrule = q_rules[dim].get();
          libmesh_assert(qrule);

          this->get_dof_map().dof_indices (elem, dof_indices, var);

          auto element_calculation = [&dof_indices, &elem,
               norm_type, norm_weight, norm_weight_sq, &qrule,
               &local_v, &v_norm](auto & fe) {
          typedef typename std::remove_reference<decltype(fe)>::type::OutputShape OutputShape;
          typedef typename TensorTools::MakeNumber<OutputShape>::type OutputNumberShape;
          typedef typename std::remove_reference<decltype(fe)>::type::OutputGradient OutputGradient;
          typedef typename TensorTools::MakeNumber<OutputGradient>::type OutputNumberGradient;

          const std::vector<Real> &                     JxW = fe.get_JxW();
          const std::vector<std::vector<OutputShape>> * phi = nullptr;
          if (norm_type == H1 ||
              norm_type == H2 ||
              norm_type == L2 ||
              norm_type == L1 ||
              norm_type == L_INF)
            phi = &(fe.get_phi());

          const std::vector<std::vector<OutputGradient>> * dphi = nullptr;
          if (norm_type == H1 ||
              norm_type == H2 ||
              norm_type == H1_SEMINORM ||
              norm_type == W1_INF_SEMINORM)
            dphi = &(fe.get_dphi());

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          typedef typename std::remove_reference<decltype(fe)>::type::OutputTensor OutputTensor;

          const std::vector<std::vector<OutputTensor>> *  d2phi = nullptr;
          if (norm_type == H2 ||
              norm_type == H2_SEMINORM ||
              norm_type == W2_INF_SEMINORM)
            d2phi = &(fe.get_d2phi());
#endif

          fe.reinit (elem);

          const unsigned int n_qp = qrule->n_points();

          const unsigned int n_sf = cast_int<unsigned int>
            (dof_indices.size());

          // Begin the loop over the Quadrature points.
          for (unsigned int qp=0; qp<n_qp; qp++)
            {
              if (norm_type == L1)
                {
                  OutputNumberShape u_h = 0.;
                  for (unsigned int i=0; i != n_sf; ++i)
                    u_h += (*phi)[i][qp] * (*local_v)(dof_indices[i]);
                  v_norm += norm_weight *
                    JxW[qp] * TensorTools::norm(u_h);
                }

              if (norm_type == L_INF)
                {
                  OutputNumberShape u_h = 0.;
                  for (unsigned int i=0; i != n_sf; ++i)
                    u_h += (*phi)[i][qp] * (*local_v)(dof_indices[i]);
                  v_norm = std::max(v_norm, norm_weight * TensorTools::norm(u_h));
                }

              if (norm_type == H1 ||
                  norm_type == H2 ||
                  norm_type == L2)
                {
                  OutputNumberShape u_h = 0.;
                  for (unsigned int i=0; i != n_sf; ++i)
                    u_h += (*phi)[i][qp] * (*local_v)(dof_indices[i]);
                  v_norm += norm_weight_sq *
                    JxW[qp] * TensorTools::norm_sq(u_h);
                }

              if (norm_type == H1 ||
                  norm_type == H2 ||
                  norm_type == H1_SEMINORM)
                {
                  OutputNumberGradient grad_u_h;
                  for (unsigned int i=0; i != n_sf; ++i)
                    grad_u_h.add_scaled((*dphi)[i][qp], (*local_v)(dof_indices[i]));
                  v_norm += norm_weight_sq *
                    JxW[qp] * grad_u_h.norm_sq();
                }

              if (norm_type == W1_INF_SEMINORM)
                {
                  OutputNumberGradient grad_u_h;
                  for (unsigned int i=0; i != n_sf; ++i)
                    grad_u_h.add_scaled((*dphi)[i][qp], (*local_v)(dof_indices[i]));
                  v_norm = std::max(v_norm, norm_weight * grad_u_h.norm());
                }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
              typedef typename TensorTools::MakeNumber<OutputTensor>::type OutputNumberTensor;

              if (norm_type == H2 ||
                  norm_type == H2_SEMINORM)
                {
                  OutputNumberTensor hess_u_h;
                  for (unsigned int i=0; i != n_sf; ++i)
                    hess_u_h.add_scaled((*d2phi)[i][qp], (*local_v)(dof_indices[i]));
                  v_norm += norm_weight_sq *
                    JxW[qp] * hess_u_h.norm_sq();
                }

              if (norm_type == W2_INF_SEMINORM)
                {
                  OutputNumberTensor hess_u_h;
                  for (unsigned int i=0; i != n_sf; ++i)
                    hess_u_h.add_scaled((*d2phi)[i][qp], (*local_v)(dof_indices[i]));
                  v_norm = std::max(v_norm, norm_weight * hess_u_h.norm());
                }
#endif
            }
          };

          FEBase * scalar_fe = fe_ptrs[dim].get();
          FEVectorBase * vec_fe = vec_fe_ptrs[dim].get();

          if (scalar_fe)
            {
              libmesh_assert(!vec_fe);
              element_calculation(*scalar_fe);
            }

          if (vec_fe)
            {
              libmesh_assert(!scalar_fe);
              element_calculation(*vec_fe);
            }
        }
    }

  if (using_hilbert_norm)
    {
      this->comm().sum(v_norm);
      v_norm = std::sqrt(v_norm);
    }
  else
    {
      this->comm().max(v_norm);
    }

  return v_norm;
}



std::string System::get_info() const
{
  std::ostringstream oss;


  const std::string & sys_name = this->name();

  oss << "   System #"  << this->number() << ", \"" << sys_name << "\"\n"
      << "    Type \""  << this->system_type() << "\"\n"
      << "    Variables=";

  for (auto vg : make_range(this->n_variable_groups()))
    {
      const VariableGroup & vg_description (this->variable_group(vg));

      if (vg_description.n_variables() > 1) oss << "{ ";
      for (auto vn : make_range(vg_description.n_variables()))
        oss << "\"" << vg_description.name(vn) << "\" ";
      if (vg_description.n_variables() > 1) oss << "} ";
    }

  oss << '\n';

  oss << "    Finite Element Types=";
#ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS
  for (auto vg : make_range(this->n_variable_groups()))
    oss << "\""
        << Utility::enum_to_string<FEFamily>(this->get_dof_map().variable_group(vg).type().family)
        << "\" ";
#else
  for (auto vg : make_range(this->n_variable_groups()))
    {
      oss << "\""
          << Utility::enum_to_string<FEFamily>(this->get_dof_map().variable_group(vg).type().family)
          << "\", \""
          << Utility::enum_to_string<FEFamily>(this->get_dof_map().variable_group(vg).type().radial_family)
          << "\" ";
    }

  oss << '\n' << "    Infinite Element Mapping=";
  for (auto vg : make_range(this->n_variable_groups()))
    oss << "\""
        << Utility::enum_to_string<InfMapType>(this->get_dof_map().variable_group(vg).type().inf_map)
        << "\" ";
#endif

  oss << '\n';

  oss << "    Approximation Orders=";
  for (auto vg : make_range(this->n_variable_groups()))
    {
#ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS
      oss << "\""
          << Utility::enum_to_string<Order>(this->get_dof_map().variable_group(vg).type().order)
          << "\" ";
#else
      oss << "\""
          << Utility::enum_to_string<Order>(this->get_dof_map().variable_group(vg).type().order)
          << "\", \""
          << Utility::enum_to_string<Order>(this->get_dof_map().variable_group(vg).type().radial_order)
          << "\" ";
#endif
    }

  oss << '\n';

  oss << "    n_dofs()="             << this->n_dofs()             << '\n';
  dof_id_type local_dofs = this->n_local_dofs();
  oss << "    n_local_dofs()="       << local_dofs                 << '\n';
  this->comm().max(local_dofs);
  oss << "    max(n_local_dofs())="       << local_dofs                 << '\n';
#ifdef LIBMESH_ENABLE_CONSTRAINTS
  oss << "    n_constrained_dofs()=" << this->n_constrained_dofs() << '\n';
  oss << "    n_local_constrained_dofs()=" << this->n_local_constrained_dofs() << '\n';
  dof_id_type local_unconstrained_dofs = this->n_local_dofs() - this->n_local_constrained_dofs();
  this->comm().max(local_unconstrained_dofs);
  oss << "    max(local unconstrained dofs)=" << local_unconstrained_dofs << '\n';
#endif

  oss << "    " << "n_vectors()="  << this->n_vectors()  << '\n';
  oss << "    " << "n_matrices()="  << this->n_matrices()  << '\n';
  //   oss << "    " << "n_additional_matrices()=" << this->n_additional_matrices() << '\n';

  oss << this->get_dof_map().get_info();

  return oss.str();
}



void System::attach_init_function (void fptr(EquationSystems & es,
                                             const std::string & name))
{
  libmesh_assert(fptr);

  if (_init_system_object != nullptr)
    {
      libmesh_warning("WARNING:  Cannot specify both initialization function and object!");

      _init_system_object = nullptr;
    }

  _init_system_function = fptr;
}



void System::attach_init_object (System::Initialization & init_in)
{
  if (_init_system_function != nullptr)
    {
      libmesh_warning("WARNING:  Cannot specify both initialization object and function!");

      _init_system_function = nullptr;
    }

  _init_system_object = &init_in;
}



void System::attach_assemble_function (void fptr(EquationSystems & es,
                                                 const std::string & name))
{
  libmesh_assert(fptr);

  if (_assemble_system_object != nullptr)
    {
      libmesh_warning("WARNING:  Cannot specify both assembly function and object!");

      _assemble_system_object = nullptr;
    }

  _assemble_system_function = fptr;
}



void System::attach_assemble_object (System::Assembly & assemble_in)
{
  if (_assemble_system_function != nullptr)
    {
      libmesh_warning("WARNING:  Cannot specify both assembly object and function!");

      _assemble_system_function = nullptr;
    }

  _assemble_system_object = &assemble_in;
}



void System::attach_constraint_function(void fptr(EquationSystems & es,
                                                  const std::string & name))
{
  libmesh_assert(fptr);

  if (_constrain_system_object != nullptr)
    {
      libmesh_warning("WARNING:  Cannot specify both constraint function and object!");

      _constrain_system_object = nullptr;
    }

  _constrain_system_function = fptr;
}



void System::attach_constraint_object (System::Constraint & constrain)
{
  if (_constrain_system_function != nullptr)
    {
      libmesh_warning("WARNING:  Cannot specify both constraint object and function!");

      _constrain_system_function = nullptr;
    }

  _constrain_system_object = &constrain;
}

bool System::has_constraint_object () const
{
  return _constrain_system_object != nullptr;
}

System::Constraint& System::get_constraint_object ()
{
  libmesh_assert_msg(_constrain_system_object,"No constraint object available.");
  return *_constrain_system_object;
}



void System::attach_QOI_function(void fptr(EquationSystems &,
                                           const std::string &,
                                           const QoISet &))
{
  libmesh_assert(fptr);

  if (_qoi_evaluate_object != nullptr)
    {
      libmesh_warning("WARNING:  Cannot specify both QOI function and object!");

      _qoi_evaluate_object = nullptr;
    }

  _qoi_evaluate_function = fptr;
}



void System::attach_QOI_object (QOI & qoi_in)
{
  if (_qoi_evaluate_function != nullptr)
    {
      libmesh_warning("WARNING:  Cannot specify both QOI object and function!");

      _qoi_evaluate_function = nullptr;
    }

  _qoi_evaluate_object = &qoi_in;
}



void System::attach_QOI_derivative(void fptr(EquationSystems &, const std::string &,
                                             const QoISet &, bool, bool))
{
  libmesh_assert(fptr);

  if (_qoi_evaluate_derivative_object != nullptr)
    {
      libmesh_warning("WARNING:  Cannot specify both QOI derivative function and object!");

      _qoi_evaluate_derivative_object = nullptr;
    }

  _qoi_evaluate_derivative_function = fptr;
}



void System::attach_QOI_derivative_object (QOIDerivative & qoi_derivative)
{
  if (_qoi_evaluate_derivative_function != nullptr)
    {
      libmesh_warning("WARNING:  Cannot specify both QOI derivative object and function!");

      _qoi_evaluate_derivative_function = nullptr;
    }

  _qoi_evaluate_derivative_object = &qoi_derivative;
}



void System::user_initialization ()
{
  // Call the user-provided initialization function,
  // if it was provided
  if (_init_system_function != nullptr)
    this->_init_system_function (_equation_systems, this->name());

  // ...or the user-provided initialization object.
  else if (_init_system_object != nullptr)
    this->_init_system_object->initialize();
}



void System::user_assembly ()
{
  // Call the user-provided assembly function,
  // if it was provided
  if (_assemble_system_function != nullptr)
    this->_assemble_system_function (_equation_systems, this->name());

  // ...or the user-provided assembly object.
  else if (_assemble_system_object != nullptr)
    this->_assemble_system_object->assemble();
}



void System::user_constrain ()
{
  // Call the user-provided constraint function,
  // if it was provided
  if (_constrain_system_function!= nullptr)
    this->_constrain_system_function(_equation_systems, this->name());

  // ...or the user-provided constraint object.
  else if (_constrain_system_object != nullptr)
    this->_constrain_system_object->constrain();
}



void System::user_QOI (const QoISet & qoi_indices)
{
  // Call the user-provided quantity of interest function,
  // if it was provided
  if (_qoi_evaluate_function != nullptr)
    this->_qoi_evaluate_function(_equation_systems, this->name(), qoi_indices);

  // ...or the user-provided QOI function object.
  else if (_qoi_evaluate_object != nullptr)
    this->_qoi_evaluate_object->qoi(qoi_indices);
}



void System::user_QOI_derivative(const QoISet & qoi_indices,
                                 bool include_liftfunc,
                                 bool apply_constraints)
{
  // Call the user-provided quantity of interest derivative,
  // if it was provided
  if (_qoi_evaluate_derivative_function != nullptr)
    this->_qoi_evaluate_derivative_function
      (_equation_systems, this->name(), qoi_indices, include_liftfunc,
       apply_constraints);

  // ...or the user-provided QOI derivative function object.
  else if (_qoi_evaluate_derivative_object != nullptr)
    this->_qoi_evaluate_derivative_object->qoi_derivative
      (qoi_indices, include_liftfunc, apply_constraints);
}


void System::init_qois(unsigned int n_qois)
{
  _qoi.resize(n_qois);
  _qoi_error_estimates.resize(n_qois);
}


void System::set_qoi(unsigned int qoi_index, Number qoi_value)
{
  libmesh_assert(qoi_index < _qoi.size());

  _qoi[qoi_index] = qoi_value;
}


Number System::get_qoi_value(unsigned int qoi_index) const
{
  libmesh_assert(qoi_index < _qoi.size());
  return _qoi[qoi_index];
}


std::vector<Number> System::get_qoi_values() const
{
  return this->_qoi;
}


void System::set_qoi(std::vector<Number> new_qoi)
{
  libmesh_assert_equal_to(this->_qoi.size(), new_qoi.size());
  this->_qoi = std::move(new_qoi);
}


void System::set_qoi_error_estimate(unsigned int qoi_index, Number qoi_error_estimate)
{
  libmesh_assert(qoi_index < _qoi_error_estimates.size());

  _qoi_error_estimates[qoi_index] = qoi_error_estimate;
}

Number System::get_qoi_error_estimate_value(unsigned int qoi_index) const
{
  libmesh_assert(qoi_index < _qoi_error_estimates.size());
  return _qoi_error_estimates[qoi_index];
}



Number System::point_value(unsigned int var,
                           const Point & p,
                           const bool insist_on_success,
                           const NumericVector<Number> *sol) const
{
  // This function must be called on every processor; there's no
  // telling where in the partition p falls.
  parallel_object_only();

  // And every processor had better agree about which point we're
  // looking for
#ifndef NDEBUG
  libmesh_assert(this->comm().verify(p(0)));
#if LIBMESH_DIM > 1
  libmesh_assert(this->comm().verify(p(1)));
#endif
#if LIBMESH_DIM > 2
  libmesh_assert(this->comm().verify(p(2)));
#endif
#endif // NDEBUG

  // Get a reference to the mesh object associated with the system object that calls this function
  const MeshBase & mesh = this->get_mesh();

  // Use an existing PointLocator or create a new one
  std::unique_ptr<PointLocatorBase> locator_ptr = mesh.sub_point_locator();
  PointLocatorBase & locator = *locator_ptr;

  if (!insist_on_success || !mesh.is_serial())
    locator.enable_out_of_mesh_mode();

  // Get a pointer to an element that contains p and allows us to
  // evaluate var
  const std::set<subdomain_id_type> & raw_subdomains =
    this->variable(var).active_subdomains();
  const std::set<subdomain_id_type> * implicit_subdomains =
    raw_subdomains.empty() ? nullptr : &raw_subdomains;
  const Elem * e = locator(p, implicit_subdomains);

  Number u = 0;

  if (e && this->get_dof_map().is_evaluable(*e, var))
    u = point_value(var, p, *e, sol);

  // If I have an element containing p, then let's let everyone know
  processor_id_type lowest_owner =
    (e && (e->processor_id() == this->processor_id())) ?
    this->processor_id() : this->n_processors();
  this->comm().min(lowest_owner);

  // Everybody should get their value from a processor that was able
  // to compute it.
  // If nobody admits owning the point, we have a problem.
  if (lowest_owner != this->n_processors())
    this->comm().broadcast(u, lowest_owner);
  else
    libmesh_assert(!insist_on_success);

  return u;
}

Number System::point_value(unsigned int var,
                           const Point & p,
                           const Elem & e,
                           const NumericVector<Number> *sol) const
{
  // Ensuring that the given point is really in the element is an
  // expensive assert, but as long as debugging is turned on we might
  // as well try to catch a particularly nasty potential error
  libmesh_assert (e.contains_point(p));

  if (!sol)
    sol = this->current_local_solution.get();

  // Get the dof map to get the proper indices for our computation
  const DofMap & dof_map = this->get_dof_map();

  // Make sure we can evaluate on this element.
  libmesh_assert (dof_map.is_evaluable(e, var));

  // Need dof_indices for phi[i][j]
  std::vector<dof_id_type> dof_indices;

  // Fill in the dof_indices for our element
  dof_map.dof_indices (&e, dof_indices, var);

  // Get the no of dofs associated with this point
  const unsigned int num_dofs = cast_int<unsigned int>
    (dof_indices.size());

  FEType fe_type = dof_map.variable_type(var);

  // Map the physical co-ordinates to the master co-ordinates
  Point coor = FEMap::inverse_map(e.dim(), &e, p);

  // get the shape function value via the FEInterface to also handle the case
  // of infinite elements correctly, the shape function is not fe->phi().
  FEComputeData fe_data(this->get_equation_systems(), coor);
  FEInterface::compute_data(e.dim(), fe_type, &e, fe_data);

  // Get ready to accumulate a value
  Number u = 0;

  for (unsigned int l=0; l<num_dofs; l++)
    {
      u += fe_data.shape[l] * (*sol)(dof_indices[l]);
    }

  return u;
}



Number System::point_value(unsigned int var, const Point & p, const Elem * e) const
{
  libmesh_assert(e);
  return this->point_value(var, p, *e);
}



Number System::point_value(unsigned int var, const Point & p, const NumericVector<Number> * sol) const
{
  return this->point_value(var, p, true, sol);
}




Gradient System::point_gradient(unsigned int var,
                                const Point & p,
                                const bool insist_on_success,
                                const NumericVector<Number> *sol) const
{
  // This function must be called on every processor; there's no
  // telling where in the partition p falls.
  parallel_object_only();

  // And every processor had better agree about which point we're
  // looking for
#ifndef NDEBUG
  libmesh_assert(this->comm().verify(p(0)));
#if LIBMESH_DIM > 1
  libmesh_assert(this->comm().verify(p(1)));
#endif
#if LIBMESH_DIM > 2
  libmesh_assert(this->comm().verify(p(2)));
#endif
#endif // NDEBUG

  // Get a reference to the mesh object associated with the system object that calls this function
  const MeshBase & mesh = this->get_mesh();

  // Use an existing PointLocator or create a new one
  std::unique_ptr<PointLocatorBase> locator_ptr = mesh.sub_point_locator();
  PointLocatorBase & locator = *locator_ptr;

  if (!insist_on_success || !mesh.is_serial())
    locator.enable_out_of_mesh_mode();

  // Get a pointer to an element that contains p and allows us to
  // evaluate var
  const std::set<subdomain_id_type> & raw_subdomains =
    this->variable(var).active_subdomains();
  const std::set<subdomain_id_type> * implicit_subdomains =
    raw_subdomains.empty() ? nullptr : &raw_subdomains;
  const Elem * e = locator(p, implicit_subdomains);

  Gradient grad_u;

  if (e && this->get_dof_map().is_evaluable(*e, var))
    grad_u = point_gradient(var, p, *e, sol);

  // If I have an element containing p, then let's let everyone know
  processor_id_type lowest_owner =
    (e && (e->processor_id() == this->processor_id())) ?
    this->processor_id() : this->n_processors();
  this->comm().min(lowest_owner);

  // Everybody should get their value from a processor that was able
  // to compute it.
  // If nobody admits owning the point, we may have a problem.
  if (lowest_owner != this->n_processors())
    this->comm().broadcast(grad_u, lowest_owner);
  else
    libmesh_assert(!insist_on_success);

  return grad_u;
}


Gradient System::point_gradient(unsigned int var,
                                const Point & p,
                                const Elem & e,
                                const NumericVector<Number> *sol) const
{
  // Ensuring that the given point is really in the element is an
  // expensive assert, but as long as debugging is turned on we might
  // as well try to catch a particularly nasty potential error
  libmesh_assert (e.contains_point(p));

  if (!sol)
    sol = this->current_local_solution.get();

  // Get the dof map to get the proper indices for our computation
  const DofMap & dof_map = this->get_dof_map();

  // write the element dimension into a separate variable.
  const unsigned int dim = e.dim();

  // Make sure we can evaluate on this element.
  libmesh_assert (dof_map.is_evaluable(e, var));

  // Need dof_indices for phi[i][j]
  std::vector<dof_id_type> dof_indices;

  // Fill in the dof_indices for our element
  dof_map.dof_indices (&e, dof_indices, var);

  // Get the no of dofs associated with this point
  const unsigned int num_dofs = cast_int<unsigned int>
    (dof_indices.size());

  FEType fe_type = dof_map.variable_type(var);

  // Map the physical co-ordinates to the master co-ordinates
  Point coor = FEMap::inverse_map(dim, &e, p);

  // get the shape function value via the FEInterface to also handle the case
  // of infinite elements correctly, the shape function is not fe->phi().
  FEComputeData fe_data(this->get_equation_systems(), coor);
  fe_data.enable_derivative();
  FEInterface::compute_data(dim, fe_type, &e, fe_data);

  // Get ready to accumulate a gradient
  Gradient grad_u;

  for (unsigned int l=0; l<num_dofs; l++)
    {
      // Chartesian coordinates have always LIBMESH_DIM entries,
      // local coordinates have as many coordinates as the element has.
      for (std::size_t v=0; v<dim; v++)
        for (std::size_t xyz=0; xyz<LIBMESH_DIM; xyz++)
          {
            // FIXME: this needs better syntax: It is matrix-vector multiplication.
            grad_u(xyz) += fe_data.local_transform[v][xyz]
              * fe_data.dshape[l](v)
              * (*sol)(dof_indices[l]);
          }
    }

  return grad_u;
}



Gradient System::point_gradient(unsigned int var, const Point & p, const Elem * e) const
{
  libmesh_assert(e);
  return this->point_gradient(var, p, *e);
}



Gradient System::point_gradient(unsigned int var, const Point & p, const NumericVector<Number> * sol) const
{
  return this->point_gradient(var, p, true, sol);
}



// We can only accumulate a hessian with --enable-second
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor System::point_hessian(unsigned int var,
                             const Point & p,
                             const bool insist_on_success,
                             const NumericVector<Number> *sol) const
{
  // This function must be called on every processor; there's no
  // telling where in the partition p falls.
  parallel_object_only();

  // And every processor had better agree about which point we're
  // looking for
#ifndef NDEBUG
  libmesh_assert(this->comm().verify(p(0)));
#if LIBMESH_DIM > 1
  libmesh_assert(this->comm().verify(p(1)));
#endif
#if LIBMESH_DIM > 2
  libmesh_assert(this->comm().verify(p(2)));
#endif
#endif // NDEBUG

  // Get a reference to the mesh object associated with the system object that calls this function
  const MeshBase & mesh = this->get_mesh();

  // Use an existing PointLocator or create a new one
  std::unique_ptr<PointLocatorBase> locator_ptr = mesh.sub_point_locator();
  PointLocatorBase & locator = *locator_ptr;

  if (!insist_on_success || !mesh.is_serial())
    locator.enable_out_of_mesh_mode();

  // Get a pointer to an element that contains p and allows us to
  // evaluate var
  const std::set<subdomain_id_type> & raw_subdomains =
    this->variable(var).active_subdomains();
  const std::set<subdomain_id_type> * implicit_subdomains =
    raw_subdomains.empty() ? nullptr : &raw_subdomains;
  const Elem * e = locator(p, implicit_subdomains);

  Tensor hess_u;

  if (e && this->get_dof_map().is_evaluable(*e, var))
    hess_u = point_hessian(var, p, *e, sol);

  // If I have an element containing p, then let's let everyone know
  processor_id_type lowest_owner =
    (e && (e->processor_id() == this->processor_id())) ?
    this->processor_id() : this->n_processors();
  this->comm().min(lowest_owner);

  // Everybody should get their value from a processor that was able
  // to compute it.
  // If nobody admits owning the point, we may have a problem.
  if (lowest_owner != this->n_processors())
    this->comm().broadcast(hess_u, lowest_owner);
  else
    libmesh_assert(!insist_on_success);

  return hess_u;
}

Tensor System::point_hessian(unsigned int var,
                             const Point & p,
                             const Elem & e,
                             const NumericVector<Number> *sol) const
{
  // Ensuring that the given point is really in the element is an
  // expensive assert, but as long as debugging is turned on we might
  // as well try to catch a particularly nasty potential error
  libmesh_assert (e.contains_point(p));

  if (!sol)
    sol = this->current_local_solution.get();

  if (e.infinite())
    libmesh_not_implemented();

  // Get the dof map to get the proper indices for our computation
  const DofMap & dof_map = this->get_dof_map();

  // Make sure we can evaluate on this element.
  libmesh_assert (dof_map.is_evaluable(e, var));

  // Need dof_indices for phi[i][j]
  std::vector<dof_id_type> dof_indices;

  // Fill in the dof_indices for our element
  dof_map.dof_indices (&e, dof_indices, var);

  // Get the no of dofs associated with this point
  const unsigned int num_dofs = cast_int<unsigned int>
    (dof_indices.size());

  FEType fe_type = dof_map.variable_type(var);

  // Build a FE again so we can calculate u(p)
  std::unique_ptr<FEBase> fe (FEBase::build(e.dim(), fe_type));

  // Map the physical co-ordinates to the master co-ordinates
  // Build a vector of point co-ordinates to send to reinit
  std::vector<Point> coor(1, FEMap::inverse_map(e.dim(), &e, p));

  // Get the values of the shape function derivatives
  const std::vector<std::vector<RealTensor>> &  d2phi = fe->get_d2phi();

  // Reinitialize the element and compute the shape function values at coor
  fe->reinit (&e, &coor);

  // Get ready to accumulate a hessian
  Tensor hess_u;

  for (unsigned int l=0; l<num_dofs; l++)
    {
      hess_u.add_scaled (d2phi[l][0], (*sol)(dof_indices[l]));
    }

  return hess_u;
}



Tensor System::point_hessian(unsigned int var, const Point & p, const Elem * e) const
{
  libmesh_assert(e);
  return this->point_hessian(var, p, *e);
}



Tensor System::point_hessian(unsigned int var, const Point & p, const NumericVector<Number> * sol) const
{
  return this->point_hessian(var, p, true, sol);
}

#else

Tensor System::point_hessian(unsigned int, const Point &, const bool,
                             const NumericVector<Number> *) const
{
  libmesh_error_msg("We can only accumulate a hessian with --enable-second");

  // Avoid compiler warnings
  return Tensor();
}

Tensor System::point_hessian(unsigned int, const Point &, const Elem &,
                             const NumericVector<Number> *) const
{
  libmesh_error_msg("We can only accumulate a hessian with --enable-second");

  // Avoid compiler warnings
  return Tensor();
}

Tensor System::point_hessian(unsigned int, const Point &, const Elem *) const
{
  libmesh_error_msg("We can only accumulate a hessian with --enable-second");

  // Avoid compiler warnings
  return Tensor();
}

Tensor System::point_hessian(unsigned int, const Point &, const NumericVector<Number> *) const
{
  libmesh_error_msg("We can only accumulate a hessian with --enable-second");

  // Avoid compiler warnings
  return Tensor();
}

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

void System::create_static_condensation()
{
  this->get_dof_map().create_static_condensation(this->get_mesh(), *this);
}

bool System::has_static_condensation() const
{
  return this->get_dof_map().has_static_condensation();
}

unsigned int System::n_vars() const
{
  return this->get_dof_map().n_vars();
}

const std::string & System::variable_name (const unsigned int i) const
{
  return this->get_dof_map().variable_name(i);
}

bool System::identify_variable_groups () const
{
  return this->get_dof_map().identify_variable_groups();
}

void System::identify_variable_groups (const bool ivg)
{
  this->get_dof_map().identify_variable_groups(ivg);
}

unsigned int System::n_components() const
{
  return this->get_dof_map().n_components(this->get_mesh());
}

unsigned int System::n_variable_groups() const
{
  return this->get_dof_map().n_variable_groups();
}

const Variable & System::variable (const unsigned int i) const
{
  return this->get_dof_map().variable(i);
}

const VariableGroup & System::variable_group (const unsigned int vg) const
{
  return this->get_dof_map().variable_group(vg);
}

unsigned int
System::variable_scalar_number (unsigned int var_num,
                                unsigned int component) const
{
  return this->get_dof_map().variable_scalar_number(var_num, component);
}

const FEType & System::variable_type (const unsigned int i) const
{
  return this->get_dof_map().variable_type(i);
}

const FEType & System::variable_type (std::string_view var) const
{
  return this->get_dof_map().variable_type(var);
}

} // namespace libMesh
