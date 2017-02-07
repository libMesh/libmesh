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



// C++ includes
#include <sstream>   // for std::ostringstream


// Local includes
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parameter_vector.h"
#include "libmesh/point.h"              // For point_value
#include "libmesh/point_locator_base.h" // For point_value
#include "libmesh/qoi_set.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/system.h"
#include "libmesh/system_norm.h"
#include "libmesh/utility.h"
#include "libmesh/elem.h"
#include "libmesh/fe_type.h"

// includes for calculate_norm, point_*
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/quadrature.h"
#include "libmesh/tensor_value.h"
#include "libmesh/vector_value.h"
#include "libmesh/tensor_tools.h"

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
  qoi                               (0),
  _init_system_function             (libmesh_nullptr),
  _init_system_object               (libmesh_nullptr),
  _assemble_system_function         (libmesh_nullptr),
  _assemble_system_object           (libmesh_nullptr),
  _constrain_system_function        (libmesh_nullptr),
  _constrain_system_object          (libmesh_nullptr),
  _qoi_evaluate_function            (libmesh_nullptr),
  _qoi_evaluate_object              (libmesh_nullptr),
  _qoi_evaluate_derivative_function (libmesh_nullptr),
  _qoi_evaluate_derivative_object   (libmesh_nullptr),
  _dof_map                          (new DofMap(number_in, es.get_mesh())),
  _equation_systems                 (es),
  _mesh                             (es.get_mesh()),
  _sys_name                         (name_in),
  _sys_number                       (number_in),
  _active                           (true),
  _solution_projection              (true),
  _basic_system_only                (false),
  _is_initialized                   (false),
  _identify_variable_groups         (true),
  _additional_data_written          (false),
  adjoint_already_solved            (false),
  _hide_output                      (false)
{
}



// No copy construction of System objects!
System::System (const System & other) :
  ReferenceCountedObject<System>(),
  ParallelObject(other),
  _equation_systems(other._equation_systems),
  _mesh(other._mesh),
  _sys_number(other._sys_number)
{
  libmesh_not_implemented();
}



System & System::operator= (const System &)
{
  libmesh_not_implemented();
}


System::~System ()
{
  // Null-out the function pointers.  Since this
  // class is getting destructed it is pointless,
  // but a good habit.
  _init_system_function =
    _assemble_system_function =
    _constrain_system_function = libmesh_nullptr;

  _qoi_evaluate_function = libmesh_nullptr;
  _qoi_evaluate_derivative_function =  libmesh_nullptr;

  // libmesh_nullptr-out user-provided objects.
  _init_system_object             = libmesh_nullptr;
  _assemble_system_object         = libmesh_nullptr;
  _constrain_system_object        = libmesh_nullptr;
  _qoi_evaluate_object            = libmesh_nullptr;
  _qoi_evaluate_derivative_object = libmesh_nullptr;

  // Clear data
  this->clear ();

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
  return _dof_map->n_dofs_on_processor (this->processor_id());
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
  _variables.clear();

  _variable_numbers.clear();

  _dof_map->clear ();

  solution->clear ();

  current_local_solution->clear ();

  // clear any user-added vectors
  {
    for (vectors_iterator pos = _vectors.begin(); pos != _vectors.end(); ++pos)
      {
        pos->second->clear ();
        delete pos->second;
        pos->second = libmesh_nullptr;
      }

    _vectors.clear();
    _vector_projections.clear();
    _vector_is_adjoint.clear();
    _vector_types.clear();
    _is_initialized = false;
  }

}



void System::init ()
{
  // Calling init() twice on the same system currently works evil
  // magic, whether done directly or via EquationSystems::read()
  libmesh_assert(!this->is_initialized());

  // First initialize any required data:
  // either only the basic System data
  if (_basic_system_only)
    System::init_data();
  // or all the derived class' data too
  else
    this->init_data();

  // If no variables have been added to this system
  // don't do anything
  if(!this->n_vars())
    return;

  // Then call the user-provided intialization function
  this->user_initialization();
}



void System::init_data ()
{
  MeshBase & mesh = this->get_mesh();

  // Add all variable groups to our underlying DofMap
  for (unsigned int vg=0; vg<this->n_variable_groups(); vg++)
    _dof_map->add_variable_group(this->variable_group(vg));

  // Distribute the degrees of freedom on the mesh
  _dof_map->distribute_dofs (mesh);

  // Recreate any user or internal constraints
  this->reinit_constraints();

  // And clean up the send_list before we first use it
  _dof_map->prepare_send_list();

  // Resize the solution conformal to the current mesh
  solution->init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);

  // Resize the current_local_solution for the current mesh
#ifdef LIBMESH_ENABLE_GHOSTED
  current_local_solution->init (this->n_dofs(), this->n_local_dofs(),
                                _dof_map->get_send_list(), false,
                                GHOSTED);
#else
  current_local_solution->init (this->n_dofs(), false, SERIAL);
#endif

  // from now on, adding additional vectors or variables can't be done
  // without immediately initializing them
  _is_initialized = true;

  // initialize & zero other vectors, if necessary
  for (vectors_iterator pos = _vectors.begin(); pos != _vectors.end(); ++pos)
    {
      ParallelType type = _vector_types[pos->first];

      if (type == GHOSTED)
        {
#ifdef LIBMESH_ENABLE_GHOSTED
          pos->second->init (this->n_dofs(), this->n_local_dofs(),
                             _dof_map->get_send_list(), false,
                             GHOSTED);
#else
          libmesh_error_msg("Cannot initialize ghosted vectors when they are not enabled.");
#endif
        }
      else if (type == SERIAL)
        {
          pos->second->init (this->n_dofs(), false, type);
        }
      else
        {
          libmesh_assert_equal_to(type, PARALLEL);
          pos->second->init (this->n_dofs(), this->n_local_dofs(), false, type);
        }
    }
}



void System::restrict_vectors ()
{
#ifdef LIBMESH_ENABLE_AMR
  // Restrict the _vectors on the coarsened cells
  for (vectors_iterator pos = _vectors.begin(); pos != _vectors.end(); ++pos)
    {
      NumericVector<Number> * v = pos->second;

      if (_vector_projections[pos->first])
        {
          this->project_vector (*v, this->vector_is_adjoint(pos->first));
        }
      else
        {
          ParallelType type = _vector_types[pos->first];

          if(type == GHOSTED)
            {
#ifdef LIBMESH_ENABLE_GHOSTED
              pos->second->init (this->n_dofs(), this->n_local_dofs(),
                                 _dof_map->get_send_list(), false,
                                 GHOSTED);
#else
              libmesh_error_msg("Cannot initialize ghosted vectors when they are not enabled.");
#endif
            }
          else
            pos->second->init (this->n_dofs(), this->n_local_dofs(), false, type);
        }
    }

  const std::vector<dof_id_type> & send_list = _dof_map->get_send_list ();

  // Restrict the solution on the coarsened cells
  if (_solution_projection)
    this->project_vector (*solution);

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
  //If no variables have been added to this system
  //don't do anything
  if(!this->n_vars())
    return;

  // Constraints get handled in EquationSystems::reinit now
  //  _dof_map->create_dof_constraints(this->get_mesh());

  // Update the solution based on the projected
  // current_local_solution.
  solution->init (this->n_dofs(), this->n_local_dofs(), true, PARALLEL);

  libmesh_assert_equal_to (solution->size(), current_local_solution->size());
  // Not true with ghosted vectors
  // libmesh_assert_equal_to (solution->size(), current_local_solution->local_size());

  const dof_id_type first_local_dof = solution->first_local_index();
  const dof_id_type local_size      = solution->local_size();

  for (dof_id_type i=0; i<local_size; i++)
    solution->set(i+first_local_dof,
                  (*current_local_solution)(i+first_local_dof));

  solution->close();
}


void System::reinit_constraints()
{
#ifdef LIBMESH_ENABLE_CONSTRAINTS
  get_dof_map().create_dof_constraints(_mesh, this->time);
  user_constrain();
  get_dof_map().process_constraints(_mesh);
#endif
  get_dof_map().prepare_send_list();
}


void System::update ()
{
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
  if(!this->n_vars())
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
  if (subset != libmesh_nullptr)
    libmesh_not_implemented();
}



void System::assemble ()
{
  // Log how long the user's assembly code takes
  LOG_SCOPE("assemble()", "System");

  // Call the user-specified assembly function
  this->user_assembly();
}



void System::assemble_qoi (const QoISet & qoi_indices)
{
  // Log how long the user's assembly code takes
  LOG_SCOPE("assemble_qoi()", "System");

  // Call the user-specified quantity of interest function
  this->user_QOI(qoi_indices);
}



void System::assemble_qoi_derivative(const QoISet & qoi_indices,
                                     bool include_liftfunc,
                                     bool apply_constraints)
{
  // Log how long the user's assembly code takes
  LOG_SCOPE("assemble_qoi_derivative()", "System");

  // Call the user-specified quantity of interest function
  this->user_QOI_derivative(qoi_indices, include_liftfunc,
                            apply_constraints);
}



void System::qoi_parameter_sensitivity (const QoISet & qoi_indices,
                                        const ParameterVector & parameters,
                                        SensitivityData & sensitivities)
{
  // Forward sensitivities are more efficient for Nq > Np
  if (qoi_indices.size(*this) > parameters.size())
    forward_qoi_parameter_sensitivity(qoi_indices, parameters, sensitivities);
  // Adjoint sensitivities are more efficient for Np > Nq,
  // and an adjoint may be more reusable than a forward
  // solution sensitivity in the Np == Nq case.
  else
    adjoint_qoi_parameter_sensitivity(qoi_indices, parameters, sensitivities);
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
        libMesh::out << "  first difference occured at index = "
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
      for (const_vectors_iterator pos = _vectors.begin();
           pos != _vectors.end(); ++pos)
        {
          if (verbose)
            libMesh::out << "   comparing vector \""
                         << pos->first << "\" ...";

          // assume they have the same name
          const NumericVector<Number> & other_system_vector =
            other_system.get_vector(pos->first);

          ov_result.push_back(pos->second->compare (other_system_vector,
                                                    threshold));

          if (verbose)
            {
              if (ov_result[ov_result.size()-1] == -1)
                libMesh::out << " identical up to threshold." << std::endl;
              else
                libMesh::out << " first difference occured at" << std::endl
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
  global_soln.resize (solution->size());

  solution->localize (global_soln);
}



void System::update_global_solution (std::vector<Number> & global_soln,
                                     const processor_id_type dest_proc) const
{
  global_soln.resize        (solution->size());

  solution->localize_to_one (global_soln, dest_proc);
}



NumericVector<Number> & System::add_vector (const std::string & vec_name,
                                            const bool projections,
                                            const ParallelType type)
{
  // Return the vector if it is already there.
  if (this->have_vector(vec_name))
    return *(_vectors[vec_name]);

  // Otherwise build the vector
  NumericVector<Number> * buf = NumericVector<Number>::build(this->comm()).release();
  _vectors.insert (std::make_pair (vec_name, buf));
  _vector_projections.insert (std::make_pair (vec_name, projections));

  _vector_types.insert (std::make_pair (vec_name, type));

  // Vectors are primal by default
  _vector_is_adjoint.insert (std::make_pair (vec_name, -1));

  // Initialize it if necessary
  if (_is_initialized)
    {
      if(type == GHOSTED)
        {
#ifdef LIBMESH_ENABLE_GHOSTED
          buf->init (this->n_dofs(), this->n_local_dofs(),
                     _dof_map->get_send_list(), false,
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

void System::remove_vector (const std::string & vec_name)
{
  //Return if the vector does not exist
  if ( !(this->have_vector(vec_name)) )
    return;

  _vectors[vec_name]->clear();
  delete _vectors[vec_name];
  _vectors[vec_name] = libmesh_nullptr;

  _vectors.erase(vec_name);
  _vector_projections.erase(vec_name);
  _vector_is_adjoint.erase(vec_name);
  _vector_types.erase(vec_name);
}

const NumericVector<Number> * System::request_vector (const std::string & vec_name) const
{
  const_vectors_iterator pos = _vectors.find(vec_name);

  if (pos == _vectors.end())
    return libmesh_nullptr;

  return pos->second;
}



NumericVector<Number> * System::request_vector (const std::string & vec_name)
{
  vectors_iterator pos = _vectors.find(vec_name);

  if (pos == _vectors.end())
    return libmesh_nullptr;

  return pos->second;
}



const NumericVector<Number> * System::request_vector (const unsigned int vec_num) const
{
  const_vectors_iterator v = vectors_begin();
  const_vectors_iterator v_end = vectors_end();
  unsigned int num = 0;
  while((num<vec_num) && (v!=v_end))
    {
      num++;
      ++v;
    }
  if (v==v_end)
    return libmesh_nullptr;
  return v->second;
}



NumericVector<Number> * System::request_vector (const unsigned int vec_num)
{
  vectors_iterator v = vectors_begin();
  vectors_iterator v_end = vectors_end();
  unsigned int num = 0;
  while((num<vec_num) && (v!=v_end))
    {
      num++;
      ++v;
    }
  if (v==v_end)
    return libmesh_nullptr;
  return v->second;
}



const NumericVector<Number> & System::get_vector (const std::string & vec_name) const
{
  // Make sure the vector exists
  const_vectors_iterator pos = _vectors.find(vec_name);

  if (pos == _vectors.end())
    libmesh_error_msg("ERROR: vector " << vec_name << " does not exist in this system!");

  return *(pos->second);
}



NumericVector<Number> & System::get_vector (const std::string & vec_name)
{
  // Make sure the vector exists
  vectors_iterator pos = _vectors.find(vec_name);

  if (pos == _vectors.end())
    libmesh_error_msg("ERROR: vector " << vec_name << " does not exist in this system!");

  return *(pos->second);
}



const NumericVector<Number> & System::get_vector (const unsigned int vec_num) const
{
  const_vectors_iterator v = vectors_begin();
  const_vectors_iterator v_end = vectors_end();
  unsigned int num = 0;
  while((num<vec_num) && (v!=v_end))
    {
      num++;
      ++v;
    }
  libmesh_assert (v != v_end);
  return *(v->second);
}



NumericVector<Number> & System::get_vector (const unsigned int vec_num)
{
  vectors_iterator v = vectors_begin();
  vectors_iterator v_end = vectors_end();
  unsigned int num = 0;
  while((num<vec_num) && (v!=v_end))
    {
      num++;
      ++v;
    }
  libmesh_assert (v != v_end);
  return *(v->second);
}



const std::string & System::vector_name (const unsigned int vec_num) const
{
  const_vectors_iterator v = vectors_begin();
  const_vectors_iterator v_end = vectors_end();
  unsigned int num = 0;
  while((num<vec_num) && (v!=v_end))
    {
      num++;
      ++v;
    }
  libmesh_assert (v != v_end);
  return v->first;
}

const std::string & System::vector_name (const NumericVector<Number> & vec_reference) const
{
  const_vectors_iterator v = vectors_begin();
  const_vectors_iterator v_end = vectors_end();

  for(; v != v_end; ++v)
    {
      // Check if the current vector is the one whose name we want
      if(&vec_reference == v->second)
        break; // exit loop if it is
    }

  // Before returning, make sure we didnt loop till the end and not find any match
  libmesh_assert (v != v_end);

  // Return the string associated with the current vector
  return v->first;
}



void System::set_vector_preservation (const std::string & vec_name,
                                      bool preserve)
{
  _vector_projections[vec_name] = preserve;
}



bool System::vector_preservation (const std::string & vec_name) const
{
  if (_vector_projections.find(vec_name) == _vector_projections.end())
    return false;

  return _vector_projections.find(vec_name)->second;
}



void System::set_vector_as_adjoint (const std::string & vec_name,
                                    int qoi_num)
{
  // We reserve -1 for vectors which get primal constraints, -2 for
  // vectors which get no constraints
  libmesh_assert_greater_equal(qoi_num, -2);
  _vector_is_adjoint[vec_name] = qoi_num;
}



int System::vector_is_adjoint (const std::string & vec_name) const
{
  libmesh_assert(_vector_is_adjoint.find(vec_name) !=
                 _vector_is_adjoint.end());

  return _vector_is_adjoint.find(vec_name)->second;
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



unsigned int System::add_variable (const std::string & var,
                                   const FEType & type,
                                   const std::set<subdomain_id_type> * const active_subdomains)
{
  libmesh_assert(!this->is_initialized());

  // Make sure the variable isn't there already
  // or if it is, that it's the type we want
  for (unsigned int v=0; v<this->n_vars(); v++)
    if (this->variable_name(v) == var)
      {
        if (this->variable_type(v) == type)
          return _variables[v].number();

        libmesh_error_msg("ERROR: incompatible variable " << var << " has already been added for this system!");
      }

  // Optimize for VariableGroups here - if the user is adding multiple
  // variables of the same FEType and subdomain restriction, catch
  // that here and add them as members of the same VariableGroup.
  //
  // start by setting this flag to whatever the user has requested
  // and then consider the conditions which should negate it.
  bool should_be_in_vg = this->identify_variable_groups();

  // No variable groups, nothing to add to
  if (!this->n_variable_groups())
    should_be_in_vg = false;

  else
    {
      VariableGroup & vg(_variable_groups.back());

      // get a pointer to their subdomain restriction, if any.
      const std::set<subdomain_id_type> * const
        their_active_subdomains (vg.implicitly_active() ?
                                 libmesh_nullptr : &vg.active_subdomains());

      // Different types?
      if (vg.type() != type)
        should_be_in_vg = false;

      // they are restricted, we aren't?
      if (their_active_subdomains && !active_subdomains)
        should_be_in_vg = false;

      // they aren't restriced, we are?
      if (!their_active_subdomains && active_subdomains)
        should_be_in_vg = false;

      if (their_active_subdomains && active_subdomains)
        // restricted to different sets?
        if (*their_active_subdomains != *active_subdomains)
          should_be_in_vg = false;

      // OK, after all that, append the variable to the vg if none of the conditions
      // were violated
      if (should_be_in_vg)
        {
          const unsigned short curr_n_vars = cast_int<unsigned short>
            (this->n_vars());

          vg.append (var);

          _variables.push_back(vg(vg.n_variables()-1));
          _variable_numbers[var] = curr_n_vars;
          return curr_n_vars;
        }
    }

  // otherwise, fall back to adding a single variable group
  return this->add_variables (std::vector<std::string>(1, var),
                              type,
                              active_subdomains);
}



unsigned int System::add_variable (const std::string & var,
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
  libmesh_assert(!this->is_initialized());

  // Make sure the variable isn't there already
  // or if it is, that it's the type we want
  for (std::size_t ov=0; ov<vars.size(); ov++)
    for (unsigned int v=0; v<this->n_vars(); v++)
      if (this->variable_name(v) == vars[ov])
        {
          if (this->variable_type(v) == type)
            return _variables[v].number();

          libmesh_error_msg("ERROR: incompatible variable " << vars[ov] << " has already been added for this system!");
        }

  const unsigned short curr_n_vars = cast_int<unsigned short>
    (this->n_vars());

  const unsigned int next_first_component = this->n_components();

  // Add the variable group to the list
  _variable_groups.push_back((active_subdomains == libmesh_nullptr) ?
                             VariableGroup(this, vars, curr_n_vars,
                                           next_first_component, type) :
                             VariableGroup(this, vars, curr_n_vars,
                                           next_first_component, type, *active_subdomains));

  const VariableGroup & vg (_variable_groups.back());

  // Add each component of the group individually
  for (std::size_t v=0; v<vars.size(); v++)
    {
      _variables.push_back (vg(v));
      _variable_numbers[vars[v]] = cast_int<unsigned short>
        (curr_n_vars+v);
    }

  libmesh_assert_equal_to ((curr_n_vars+vars.size()), this->n_vars());

  // BSK - Defer this now to System::init_data() so we can detect
  // VariableGroups 12/28/2012
  // // Add the variable group to the _dof_map
  // _dof_map->add_variable_group (vg);

  // Return the number of the new variable
  return cast_int<unsigned int>(curr_n_vars+vars.size()-1);
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



bool System::has_variable (const std::string & var) const
{
  return _variable_numbers.count(var);
}



unsigned short int System::variable_number (const std::string & var) const
{
  // Make sure the variable exists
  std::map<std::string, unsigned short int>::const_iterator
    pos = _variable_numbers.find(var);

  if (pos == _variable_numbers.end())
    libmesh_error_msg("ERROR: variable " << var << " does not exist in this system!");

  libmesh_assert_equal_to (_variables[pos->second].name(), var);

  return pos->second;
}


void System::get_all_variable_numbers(std::vector<unsigned int> & all_variable_numbers) const
{
  all_variable_numbers.resize(n_vars());

  // Make sure the variable exists
  std::map<std::string, unsigned short int>::const_iterator
    it = _variable_numbers.begin();
  std::map<std::string, unsigned short int>::const_iterator
    it_end = _variable_numbers.end();

  unsigned int count = 0;
  for( ; it != it_end; ++it)
    {
      all_variable_numbers[count] = it->second;
      count++;
    }
}


void System::local_dof_indices(const unsigned int var,
                               std::set<dof_id_type> & var_indices) const
{
  // Make sure the set is clear
  var_indices.clear();

  std::vector<dof_id_type> dof_indices;

  // Begin the loop over the elements
  MeshBase::const_element_iterator       el     =
    this->get_mesh().active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    this->get_mesh().active_local_elements_end();

  const dof_id_type
    first_local = this->get_dof_map().first_dof(),
    end_local   = this->get_dof_map().end_dof();

  for ( ; el != end_el; ++el)
    {
      const Elem * elem = *el;
      this->get_dof_map().dof_indices (elem, dof_indices, var);

      for (std::size_t i=0; i<dof_indices.size(); i++)
        {
          dof_id_type dof = dof_indices[i];

          //If the dof is owned by the local processor
          if(first_local <= dof && dof < end_local)
            var_indices.insert(dof_indices[i]);
        }
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

  /* Loop over nodes.  */
  {
    MeshBase::const_node_iterator it = mesh.local_nodes_begin();
    const MeshBase::const_node_iterator end_it = mesh.local_nodes_end();
    for ( ; it != end_it; ++it)
      {
        const Node * node = *it;
        unsigned int n_comp = node->n_comp(sys_num,var_num);
        for(unsigned int i=0; i<n_comp; i++)
          {
            const dof_id_type index = node->dof_number(sys_num,var_num,i);
            v.set(index,0.0);
          }
      }
  }

  /* Loop over elements.  */
  {
    MeshBase::const_element_iterator it = mesh.active_local_elements_begin();
    const MeshBase::const_element_iterator end_it = mesh.active_local_elements_end();
    for ( ; it != end_it; ++it)
      {
        const Elem * elem = *it;
        unsigned int n_comp = elem->n_comp(sys_num,var_num);
        for(unsigned int i=0; i<n_comp; i++)
          {
            const dof_id_type index = elem->dof_number(sys_num,var_num,i);
            v.set(index,0.0);
          }
      }
  }
}



Real System::discrete_var_norm(const NumericVector<Number> & v,
                               unsigned int var,
                               FEMNormType norm_type) const
{
  std::set<dof_id_type> var_indices;
  local_dof_indices(var, var_indices);

  if(norm_type == DISCRETE_L1)
    return v.subset_l1_norm(var_indices);
  if(norm_type == DISCRETE_L2)
    return v.subset_l2_norm(var_indices);
  if(norm_type == DISCRETE_L_INF)
    return v.subset_linfty_norm(var_indices);
  else
    libmesh_error_msg("Invalid norm_type = " << norm_type);
}



Real System::calculate_norm(const NumericVector<Number> & v,
                            unsigned int var,
                            FEMNormType norm_type,
                            std::set<unsigned int> * skip_dimensions) const
{
  //short circuit to save time
  if(norm_type == DISCRETE_L1 ||
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
      unsigned int check_var = 0;
      for (; check_var != this->n_vars(); ++check_var)
        if((norm.weight(check_var) != 1.0) || (norm.type(check_var) != norm_type0))
          break;

      //All weights were 1.0 so just do the full vector discrete norm
      if(check_var == this->n_vars())
        {
          if(norm_type0 == DISCRETE_L1)
            return v.l1_norm();
          if(norm_type0 == DISCRETE_L2)
            return v.l2_norm();
          if(norm_type0 == DISCRETE_L_INF)
            return v.linfty_norm();
          else
            libmesh_error_msg("Invalid norm_type0 = " << norm_type0);
        }

      for (unsigned int var=0; var != this->n_vars(); ++var)
        {
          // Skip any variables we don't need to integrate
          if (norm.weight(var) == 0.0)
            continue;

          v_norm += norm.weight(var) * discrete_var_norm(v, var, norm.type(var));
        }

      return v_norm;
    }

  // Localize the potentially parallel vector
  UniquePtr<NumericVector<Number> > local_v = NumericVector<Number>::build(this->comm());
  local_v->init(v.size(), true, SERIAL);
  v.localize (*local_v, _dof_map->get_send_list());

  // I'm not sure how best to mix Hilbert norms on some variables (for
  // which we'll want to square then sum then square root) with norms
  // like L_inf (for which we'll just want to take an absolute value
  // and then sum).
  bool using_hilbert_norm = true,
    using_nonhilbert_norm = true;

  // Loop over all variables
  for (unsigned int var=0; var != this->n_vars(); ++var)
    {
      // Skip any variables we don't need to integrate
      Real norm_weight_sq = norm.weight_sq(var);
      if (norm_weight_sq == 0.0)
        continue;
      Real norm_weight = norm.weight(var);

      // Check for unimplemented norms (rather than just returning 0).
      FEMNormType norm_type = norm.type(var);
      if((norm_type==H1) ||
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

      // Allow space for dims 0-3, even if we don't use them all
      std::vector<FEBase *> fe_ptrs(4,libmesh_nullptr);
      std::vector<QBase *> q_rules(4,libmesh_nullptr);

      const std::set<unsigned char> & elem_dims = _mesh.elem_dimensions();

      // Prepare finite elements for each dimension present in the mesh
      for (std::set<unsigned char>::const_iterator d_it = elem_dims.begin();
           d_it != elem_dims.end(); ++d_it)
        {
          if (skip_dimensions && skip_dimensions->find(*d_it) != skip_dimensions->end())
            continue;

          q_rules[*d_it] =
            fe_type.default_quadrature_rule (*d_it).release();

          // Construct finite element object

          fe_ptrs[*d_it] = FEBase::build(*d_it, fe_type).release();

          // Attach quadrature rule to FE object
          fe_ptrs[*d_it]->attach_quadrature_rule (q_rules[*d_it]);
        }

      std::vector<dof_id_type> dof_indices;

      // Begin the loop over the elements
      MeshBase::const_element_iterator       el     =
        this->get_mesh().active_local_elements_begin();
      const MeshBase::const_element_iterator end_el =
        this->get_mesh().active_local_elements_end();

      for ( ; el != end_el; ++el)
        {
          const Elem * elem = *el;
          const unsigned int dim = elem->dim();

          if (skip_dimensions && skip_dimensions->find(dim) != skip_dimensions->end())
            continue;

          FEBase * fe = fe_ptrs[dim];
          QBase * qrule = q_rules[dim];
          libmesh_assert(fe);
          libmesh_assert(qrule);

          const std::vector<Real> &               JxW = fe->get_JxW();
          const std::vector<std::vector<Real> > * phi = libmesh_nullptr;
          if (norm_type == H1 ||
              norm_type == H2 ||
              norm_type == L2 ||
              norm_type == L1 ||
              norm_type == L_INF)
            phi = &(fe->get_phi());

          const std::vector<std::vector<RealGradient> > * dphi = libmesh_nullptr;
          if (norm_type == H1 ||
              norm_type == H2 ||
              norm_type == H1_SEMINORM ||
              norm_type == W1_INF_SEMINORM)
            dphi = &(fe->get_dphi());
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
          const std::vector<std::vector<RealTensor> > *   d2phi = libmesh_nullptr;
          if (norm_type == H2 ||
              norm_type == H2_SEMINORM ||
              norm_type == W2_INF_SEMINORM)
            d2phi = &(fe->get_d2phi());
#endif

          fe->reinit (elem);

          this->get_dof_map().dof_indices (elem, dof_indices, var);

          const unsigned int n_qp = qrule->n_points();

          const unsigned int n_sf = cast_int<unsigned int>
            (dof_indices.size());

          // Begin the loop over the Quadrature points.
          for (unsigned int qp=0; qp<n_qp; qp++)
            {
              if (norm_type == L1)
                {
                  Number u_h = 0.;
                  for (unsigned int i=0; i != n_sf; ++i)
                    u_h += (*phi)[i][qp] * (*local_v)(dof_indices[i]);
                  v_norm += norm_weight *
                    JxW[qp] * std::abs(u_h);
                }

              if (norm_type == L_INF)
                {
                  Number u_h = 0.;
                  for (unsigned int i=0; i != n_sf; ++i)
                    u_h += (*phi)[i][qp] * (*local_v)(dof_indices[i]);
                  v_norm = std::max(v_norm, norm_weight * std::abs(u_h));
                }

              if (norm_type == H1 ||
                  norm_type == H2 ||
                  norm_type == L2)
                {
                  Number u_h = 0.;
                  for (unsigned int i=0; i != n_sf; ++i)
                    u_h += (*phi)[i][qp] * (*local_v)(dof_indices[i]);
                  v_norm += norm_weight_sq *
                    JxW[qp] * TensorTools::norm_sq(u_h);
                }

              if (norm_type == H1 ||
                  norm_type == H2 ||
                  norm_type == H1_SEMINORM)
                {
                  Gradient grad_u_h;
                  for (unsigned int i=0; i != n_sf; ++i)
                    grad_u_h.add_scaled((*dphi)[i][qp], (*local_v)(dof_indices[i]));
                  v_norm += norm_weight_sq *
                    JxW[qp] * grad_u_h.norm_sq();
                }

              if (norm_type == W1_INF_SEMINORM)
                {
                  Gradient grad_u_h;
                  for (unsigned int i=0; i != n_sf; ++i)
                    grad_u_h.add_scaled((*dphi)[i][qp], (*local_v)(dof_indices[i]));
                  v_norm = std::max(v_norm, norm_weight * grad_u_h.norm());
                }

#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
              if (norm_type == H2 ||
                  norm_type == H2_SEMINORM)
                {
                  Tensor hess_u_h;
                  for (unsigned int i=0; i != n_sf; ++i)
                    hess_u_h.add_scaled((*d2phi)[i][qp], (*local_v)(dof_indices[i]));
                  v_norm += norm_weight_sq *
                    JxW[qp] * hess_u_h.norm_sq();
                }

              if (norm_type == W2_INF_SEMINORM)
                {
                  Tensor hess_u_h;
                  for (unsigned int i=0; i != n_sf; ++i)
                    hess_u_h.add_scaled((*d2phi)[i][qp], (*local_v)(dof_indices[i]));
                  v_norm = std::max(v_norm, norm_weight * hess_u_h.norm());
                }
#endif
            }
        }

      // Need to delete the FE and quadrature objects to prevent a memory leak
      for (std::size_t i=0; i<fe_ptrs.size(); i++)
        if (fe_ptrs[i])
          delete fe_ptrs[i];

      for (std::size_t i=0; i<q_rules.size(); i++)
        if (q_rules[i])
          delete q_rules[i];
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

  for (unsigned int vg=0; vg<this->n_variable_groups(); vg++)
    {
      const VariableGroup & vg_description (this->variable_group(vg));

      if (vg_description.n_variables() > 1) oss << "{ ";
      for (unsigned int vn=0; vn<vg_description.n_variables(); vn++)
        oss << "\"" << vg_description.name(vn) << "\" ";
      if (vg_description.n_variables() > 1) oss << "} ";
    }

  oss << '\n';

  oss << "    Finite Element Types=";
#ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS
  for (unsigned int vg=0; vg<this->n_variable_groups(); vg++)
    oss << "\""
        << Utility::enum_to_string<FEFamily>(this->get_dof_map().variable_group(vg).type().family)
        << "\" ";
#else
  for (unsigned int vg=0; vg<this->n_variable_groups(); vg++)
    {
      oss << "\""
          << Utility::enum_to_string<FEFamily>(this->get_dof_map().variable_group(vg).type().family)
          << "\", \""
          << Utility::enum_to_string<FEFamily>(this->get_dof_map().variable_group(vg).type().radial_family)
          << "\" ";
    }

  oss << '\n' << "    Infinite Element Mapping=";
  for (unsigned int vg=0; vg<this->n_variable_groups(); vg++)
    oss << "\""
        << Utility::enum_to_string<InfMapType>(this->get_dof_map().variable_group(vg).type().inf_map)
        << "\" ";
#endif

  oss << '\n';

  oss << "    Approximation Orders=";
  for (unsigned int vg=0; vg<this->n_variable_groups(); vg++)
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
  oss << "    n_local_dofs()="       << this->n_local_dofs()       << '\n';
#ifdef LIBMESH_ENABLE_CONSTRAINTS
  oss << "    n_constrained_dofs()=" << this->n_constrained_dofs() << '\n';
  oss << "    n_local_constrained_dofs()=" << this->n_local_constrained_dofs() << '\n';
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

  if (_init_system_object != libmesh_nullptr)
    {
      libmesh_here();
      libMesh::out << "WARNING:  Cannot specify both initialization function and object!"
                   << std::endl;

      _init_system_object = libmesh_nullptr;
    }

  _init_system_function = fptr;
}



void System::attach_init_object (System::Initialization & init_in)
{
  if (_init_system_function != libmesh_nullptr)
    {
      libmesh_here();
      libMesh::out << "WARNING:  Cannot specify both initialization object and function!"
                   << std::endl;

      _init_system_function = libmesh_nullptr;
    }

  _init_system_object = &init_in;
}



void System::attach_assemble_function (void fptr(EquationSystems & es,
                                                 const std::string & name))
{
  libmesh_assert(fptr);

  if (_assemble_system_object != libmesh_nullptr)
    {
      libmesh_here();
      libMesh::out << "WARNING:  Cannot specify both assembly function and object!"
                   << std::endl;

      _assemble_system_object = libmesh_nullptr;
    }

  _assemble_system_function = fptr;
}



void System::attach_assemble_object (System::Assembly & assemble_in)
{
  if (_assemble_system_function != libmesh_nullptr)
    {
      libmesh_here();
      libMesh::out << "WARNING:  Cannot specify both assembly object and function!"
                   << std::endl;

      _assemble_system_function = libmesh_nullptr;
    }

  _assemble_system_object = &assemble_in;
}



void System::attach_constraint_function(void fptr(EquationSystems & es,
                                                  const std::string & name))
{
  libmesh_assert(fptr);

  if (_constrain_system_object != libmesh_nullptr)
    {
      libmesh_here();
      libMesh::out << "WARNING:  Cannot specify both constraint function and object!"
                   << std::endl;

      _constrain_system_object = libmesh_nullptr;
    }

  _constrain_system_function = fptr;
}



void System::attach_constraint_object (System::Constraint & constrain)
{
  if (_constrain_system_function != libmesh_nullptr)
    {
      libmesh_here();
      libMesh::out << "WARNING:  Cannot specify both constraint object and function!"
                   << std::endl;

      _constrain_system_function = libmesh_nullptr;
    }

  _constrain_system_object = &constrain;
}



void System::attach_QOI_function(void fptr(EquationSystems &,
                                           const std::string &,
                                           const QoISet &))
{
  libmesh_assert(fptr);

  if (_qoi_evaluate_object != libmesh_nullptr)
    {
      libmesh_here();
      libMesh::out << "WARNING:  Cannot specify both QOI function and object!"
                   << std::endl;

      _qoi_evaluate_object = libmesh_nullptr;
    }

  _qoi_evaluate_function = fptr;
}



void System::attach_QOI_object (QOI & qoi_in)
{
  if (_qoi_evaluate_function != libmesh_nullptr)
    {
      libmesh_here();
      libMesh::out << "WARNING:  Cannot specify both QOI object and function!"
                   << std::endl;

      _qoi_evaluate_function = libmesh_nullptr;
    }

  _qoi_evaluate_object = &qoi_in;
}



void System::attach_QOI_derivative(void fptr(EquationSystems &, const std::string &,
                                             const QoISet &, bool, bool))
{
  libmesh_assert(fptr);

  if (_qoi_evaluate_derivative_object != libmesh_nullptr)
    {
      libmesh_here();
      libMesh::out << "WARNING:  Cannot specify both QOI derivative function and object!"
                   << std::endl;

      _qoi_evaluate_derivative_object = libmesh_nullptr;
    }

  _qoi_evaluate_derivative_function = fptr;
}



void System::attach_QOI_derivative_object (QOIDerivative & qoi_derivative)
{
  if (_qoi_evaluate_derivative_function != libmesh_nullptr)
    {
      libmesh_here();
      libMesh::out << "WARNING:  Cannot specify both QOI derivative object and function!"
                   << std::endl;

      _qoi_evaluate_derivative_function = libmesh_nullptr;
    }

  _qoi_evaluate_derivative_object = &qoi_derivative;
}



void System::user_initialization ()
{
  // Call the user-provided intialization function,
  // if it was provided
  if (_init_system_function != libmesh_nullptr)
    this->_init_system_function (_equation_systems, this->name());

  // ...or the user-provided initialization object.
  else if (_init_system_object != libmesh_nullptr)
    this->_init_system_object->initialize();
}



void System::user_assembly ()
{
  // Call the user-provided assembly function,
  // if it was provided
  if (_assemble_system_function != libmesh_nullptr)
    this->_assemble_system_function (_equation_systems, this->name());

  // ...or the user-provided assembly object.
  else if (_assemble_system_object != libmesh_nullptr)
    this->_assemble_system_object->assemble();
}



void System::user_constrain ()
{
  // Call the user-provided constraint function,
  // if it was provided
  if (_constrain_system_function!= libmesh_nullptr)
    this->_constrain_system_function(_equation_systems, this->name());

  // ...or the user-provided constraint object.
  else if (_constrain_system_object != libmesh_nullptr)
    this->_constrain_system_object->constrain();
}



void System::user_QOI (const QoISet & qoi_indices)
{
  // Call the user-provided quantity of interest function,
  // if it was provided
  if (_qoi_evaluate_function != libmesh_nullptr)
    this->_qoi_evaluate_function(_equation_systems, this->name(), qoi_indices);

  // ...or the user-provided QOI function object.
  else if (_qoi_evaluate_object != libmesh_nullptr)
    this->_qoi_evaluate_object->qoi(qoi_indices);
}



void System::user_QOI_derivative(const QoISet & qoi_indices,
                                 bool include_liftfunc,
                                 bool apply_constraints)
{
  // Call the user-provided quantity of interest derivative,
  // if it was provided
  if (_qoi_evaluate_derivative_function != libmesh_nullptr)
    this->_qoi_evaluate_derivative_function
      (_equation_systems, this->name(), qoi_indices, include_liftfunc,
       apply_constraints);

  // ...or the user-provided QOI derivative function object.
  else if (_qoi_evaluate_derivative_object != libmesh_nullptr)
    this->_qoi_evaluate_derivative_object->qoi_derivative
      (qoi_indices, include_liftfunc, apply_constraints);
}



Number System::point_value(unsigned int var, const Point & p, const bool insist_on_success) const
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
  UniquePtr<PointLocatorBase> locator_ptr = mesh.sub_point_locator();
  PointLocatorBase & locator = *locator_ptr;

  if (!insist_on_success)
    locator.enable_out_of_mesh_mode();

  // Get a pointer to the element that contains P
  const Elem * e = locator(p);

  Number u = 0;

  if (e && this->get_dof_map().is_evaluable(*e, var))
    u = point_value(var, p, *e);

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

Number System::point_value(unsigned int var, const Point & p, const Elem & e) const
{
  // Ensuring that the given point is really in the element is an
  // expensive assert, but as long as debugging is turned on we might
  // as well try to catch a particularly nasty potential error
  libmesh_assert (e.contains_point(p));

  // Get the dof map to get the proper indices for our computation
  const DofMap & dof_map = this->get_dof_map();

  // Make sure we can evaluate on this element.
  libmesh_assert (dof_map.is_evaluable(e, var));

  // Need dof_indices for phi[i][j]
  std::vector<dof_id_type> dof_indices;

  // Fill in the dof_indices for our element
  dof_map.dof_indices (&e, dof_indices, var);

  // Get the no of dofs assciated with this point
  const unsigned int num_dofs = cast_int<unsigned int>
    (dof_indices.size());

  FEType fe_type = dof_map.variable_type(var);

  // Build a FE so we can calculate u(p)
  UniquePtr<FEBase> fe (FEBase::build(e.dim(), fe_type));

  // Map the physical co-ordinates to the master co-ordinates using the inverse_map from fe_interface.h
  // Build a vector of point co-ordinates to send to reinit
  std::vector<Point> coor(1, FEInterface::inverse_map(e.dim(), fe_type, &e, p));

  // Get the shape function values
  const std::vector<std::vector<Real> > & phi = fe->get_phi();

  // Reinitialize the element and compute the shape function values at coor
  fe->reinit (&e, &coor);

  // Get ready to accumulate a value
  Number u = 0;

  for (unsigned int l=0; l<num_dofs; l++)
    {
      u += phi[l][0]*this->current_solution (dof_indices[l]);
    }

  return u;
}



Number System::point_value(unsigned int var, const Point & p, const Elem * e) const
{
  libmesh_assert(e);
  return this->point_value(var, p, *e);
}



Gradient System::point_gradient(unsigned int var, const Point & p, const bool insist_on_success) const
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
  UniquePtr<PointLocatorBase> locator_ptr = mesh.sub_point_locator();
  PointLocatorBase & locator = *locator_ptr;

  if (!insist_on_success)
    locator.enable_out_of_mesh_mode();

  // Get a pointer to the element that contains P
  const Elem * e = locator(p);

  Gradient grad_u;

  if (e && this->get_dof_map().is_evaluable(*e, var))
    grad_u = point_gradient(var, p, *e);

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


Gradient System::point_gradient(unsigned int var, const Point & p, const Elem & e) const
{
  // Ensuring that the given point is really in the element is an
  // expensive assert, but as long as debugging is turned on we might
  // as well try to catch a particularly nasty potential error
  libmesh_assert (e.contains_point(p));

  // Get the dof map to get the proper indices for our computation
  const DofMap & dof_map = this->get_dof_map();

  // Make sure we can evaluate on this element.
  libmesh_assert (dof_map.is_evaluable(e, var));

  // Need dof_indices for phi[i][j]
  std::vector<dof_id_type> dof_indices;

  // Fill in the dof_indices for our element
  dof_map.dof_indices (&e, dof_indices, var);

  // Get the no of dofs assciated with this point
  const unsigned int num_dofs = cast_int<unsigned int>
    (dof_indices.size());

  FEType fe_type = dof_map.variable_type(var);

  // Build a FE again so we can calculate u(p)
  UniquePtr<FEBase> fe (FEBase::build(e.dim(), fe_type));

  // Map the physical co-ordinates to the master co-ordinates using the inverse_map from fe_interface.h
  // Build a vector of point co-ordinates to send to reinit
  std::vector<Point> coor(1, FEInterface::inverse_map(e.dim(), fe_type, &e, p));

  // Get the values of the shape function derivatives
  const std::vector<std::vector<RealGradient> > &  dphi = fe->get_dphi();

  // Reinitialize the element and compute the shape function values at coor
  fe->reinit (&e, &coor);

  // Get ready to accumulate a gradient
  Gradient grad_u;

  for (unsigned int l=0; l<num_dofs; l++)
    {
      grad_u.add_scaled (dphi[l][0], this->current_solution (dof_indices[l]));
    }

  return grad_u;
}



Gradient System::point_gradient(unsigned int var, const Point & p, const Elem * e) const
{
  libmesh_assert(e);
  return this->point_gradient(var, p, *e);
}



// We can only accumulate a hessian with --enable-second
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
Tensor System::point_hessian(unsigned int var, const Point & p, const bool insist_on_success) const
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
  UniquePtr<PointLocatorBase> locator_ptr = mesh.sub_point_locator();
  PointLocatorBase & locator = *locator_ptr;

  if (!insist_on_success)
    locator.enable_out_of_mesh_mode();

  // Get a pointer to the element that contains P
  const Elem * e = locator(p);

  Tensor hess_u;

  if (e && this->get_dof_map().is_evaluable(*e, var))
    hess_u = point_hessian(var, p, *e);

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

Tensor System::point_hessian(unsigned int var, const Point & p, const Elem & e) const
{
  // Ensuring that the given point is really in the element is an
  // expensive assert, but as long as debugging is turned on we might
  // as well try to catch a particularly nasty potential error
  libmesh_assert (e.contains_point(p));

  // Get the dof map to get the proper indices for our computation
  const DofMap & dof_map = this->get_dof_map();

  // Make sure we can evaluate on this element.
  libmesh_assert (dof_map.is_evaluable(e, var));

  // Need dof_indices for phi[i][j]
  std::vector<dof_id_type> dof_indices;

  // Fill in the dof_indices for our element
  dof_map.dof_indices (&e, dof_indices, var);

  // Get the no of dofs assciated with this point
  const unsigned int num_dofs = cast_int<unsigned int>
    (dof_indices.size());

  FEType fe_type = dof_map.variable_type(var);

  // Build a FE again so we can calculate u(p)
  UniquePtr<FEBase> fe (FEBase::build(e.dim(), fe_type));

  // Map the physical co-ordinates to the master co-ordinates using the inverse_map from fe_interface.h
  // Build a vector of point co-ordinates to send to reinit
  std::vector<Point> coor(1, FEInterface::inverse_map(e.dim(), fe_type, &e, p));

  // Get the values of the shape function derivatives
  const std::vector<std::vector<RealTensor> > &  d2phi = fe->get_d2phi();

  // Reinitialize the element and compute the shape function values at coor
  fe->reinit (&e, &coor);

  // Get ready to accumulate a hessian
  Tensor hess_u;

  for (unsigned int l=0; l<num_dofs; l++)
    {
      hess_u.add_scaled (d2phi[l][0], this->current_solution (dof_indices[l]));
    }

  return hess_u;
}



Tensor System::point_hessian(unsigned int var, const Point & p, const Elem * e) const
{
  libmesh_assert(e);
  return this->point_hessian(var, p, *e);
}



#else
Tensor System::point_hessian(unsigned int, const Point &, const bool) const
{
  libmesh_error_msg("We can only accumulate a hessian with --enable-second");

  // Avoid compiler warnings
  return Tensor();
}

Tensor System::point_hessian(unsigned int, const Point &, const Elem &) const
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

#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES

} // namespace libMesh
