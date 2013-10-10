// The libMesh Finite Element Library.
// Copyright (C) 2002-2012 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_base.h"
#include "libmesh/fem_context.h"
#include "libmesh/fem_system.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel.h"
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel_ghost_sync.h"
#include "libmesh/quadrature.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/time_solver.h"
#include "libmesh/unsteady_solver.h" // For eulerian_residual
#include "libmesh/fe_interface.h"

namespace {
  using namespace libMesh;

  // give this guy some scope since there
  // is underlying vector allocation upon
  // creation/deletion
  ConstElemRange elem_range;

  typedef Threads::spin_mutex femsystem_mutex;
  femsystem_mutex assembly_mutex;

  void assemble_unconstrained_element_system
    (const FEMSystem& _sys,
     const bool _get_jacobian,
     FEMContext &_femcontext)
    {
      bool jacobian_computed =
	_sys.time_solver->element_residual(_get_jacobian, _femcontext);

      // Compute a numeric jacobian if we have to
      if (_get_jacobian && !jacobian_computed)
        {
          // Make sure we didn't compute a jacobian and lie about it
          libmesh_assert_equal_to (_femcontext.get_elem_jacobian().l1_norm(), 0.0);
          // Logging of numerical jacobians is done separately
          _sys.numerical_elem_jacobian(_femcontext);
        }

      // Compute a numeric jacobian if we're asked to verify the
      // analytic jacobian we got
      if (_get_jacobian && jacobian_computed &&
          _sys.verify_analytic_jacobians != 0.0)
        {
          DenseMatrix<Number> analytic_jacobian(_femcontext.get_elem_jacobian());

          _femcontext.get_elem_jacobian().zero();
          // Logging of numerical jacobians is done separately
          _sys.numerical_elem_jacobian(_femcontext);

          Real analytic_norm = analytic_jacobian.l1_norm();
          Real numerical_norm = _femcontext.get_elem_jacobian().l1_norm();

          // If we can continue, we'll probably prefer the analytic jacobian
          analytic_jacobian.swap(_femcontext.get_elem_jacobian());

          // The matrix "analytic_jacobian" will now hold the error matrix
          analytic_jacobian.add(-1.0, _femcontext.get_elem_jacobian());
          Real error_norm = analytic_jacobian.l1_norm();

          Real relative_error = error_norm /
                                std::max(analytic_norm, numerical_norm);

          if (relative_error > _sys.verify_analytic_jacobians)
            {
              libMesh::err << "Relative error " << relative_error
                           << " detected in analytic jacobian on element "
                           << _femcontext.get_elem().id() << '!' << std::endl;

	      std::streamsize old_precision = libMesh::out.precision();
              libMesh::out.precision(16);
	      libMesh::out << "J_analytic " << _femcontext.get_elem().id() << " = "
                           << _femcontext.get_elem_jacobian() << std::endl;
              analytic_jacobian.add(1.0, _femcontext.get_elem_jacobian());
	      libMesh::out << "J_numeric " << _femcontext.get_elem().id() << " = "
                           << analytic_jacobian << std::endl;

              libMesh::out.precision(old_precision);

              libmesh_error();
            }
        }

      for (_femcontext.side = 0;
           _femcontext.side != _femcontext.get_elem().n_sides();
           ++_femcontext.side)
        {
          // Don't compute on non-boundary sides unless requested
          if (!_sys.get_physics()->compute_internal_sides &&
              _femcontext.get_elem().neighbor(_femcontext.side) != NULL)
            continue;

          // Any mesh movement has already been done (and restored,
          // if the TimeSolver isn't broken), but
          // reinitializing the side FE objects is still necessary
          _femcontext.side_fe_reinit();

          DenseMatrix<Number> old_jacobian;
          // If we're in DEBUG mode, we should always verify that the
	  // user's side_residual function doesn't alter our existing
          // jacobian and then lie about it
#ifndef DEBUG
	  // Even if we're not in DEBUG mode, when we're verifying
	  // analytic jacobians we'll want to verify each side's
          // jacobian contribution separately.
	  /* PB: We also need to account for the case when the user wants to
		 use numerical Jacobians and not analytic Jacobians */
          if ( (_sys.verify_analytic_jacobians != 0.0 && _get_jacobian) ||
	       (!jacobian_computed && _get_jacobian) )
#endif // ifndef DEBUG
            {
              old_jacobian = _femcontext.get_elem_jacobian();
              _femcontext.get_elem_jacobian().zero();
            }
	  jacobian_computed =
            _sys.time_solver->side_residual(_get_jacobian, _femcontext);

          // Compute a numeric jacobian if we have to
          if (_get_jacobian && !jacobian_computed)
            {
	      // In DEBUG mode, we've already set elem_jacobian == 0,
	      // so we can make sure side_residual didn't compute a
              // jacobian and lie about it
#ifdef DEBUG
              libmesh_assert_equal_to (_femcontext.get_elem_jacobian().l1_norm(), 0.0);
#endif
              // Logging of numerical jacobians is done separately
              _sys.numerical_side_jacobian(_femcontext);

	      // Add back in element interior numerical Jacobian
	      _femcontext.get_elem_jacobian() += old_jacobian;
            }

          // Compute a numeric jacobian if we're asked to verify the
          // analytic jacobian we got
	  if (_get_jacobian && jacobian_computed &&
              _sys.verify_analytic_jacobians != 0.0)
            {
	      DenseMatrix<Number> analytic_jacobian(_femcontext.get_elem_jacobian());

              _femcontext.get_elem_jacobian().zero();
              // Logging of numerical jacobians is done separately
              _sys.numerical_side_jacobian(_femcontext);

              Real analytic_norm = analytic_jacobian.l1_norm();
              Real numerical_norm = _femcontext.get_elem_jacobian().l1_norm();

              // If we can continue, we'll probably prefer the analytic jacobian
              analytic_jacobian.swap(_femcontext.get_elem_jacobian());

              // The matrix "analytic_jacobian" will now hold the error matrix
              analytic_jacobian.add(-1.0, _femcontext.get_elem_jacobian());
              Real error_norm = analytic_jacobian.l1_norm();

              Real relative_error = error_norm /
                                    std::max(analytic_norm, numerical_norm);

              if (relative_error > _sys.verify_analytic_jacobians)
                {
                  libMesh::err << "Relative error " << relative_error
                               << " detected in analytic jacobian on element "
                               << _femcontext.get_elem().id()
			       << ", side "
                               << static_cast<unsigned int>(_femcontext.side) << '!' << std::endl;

		  std::streamsize old_precision = libMesh::out.precision();
                  libMesh::out.precision(16);
	          libMesh::out << "J_analytic " << _femcontext.get_elem().id() << " = "
                               << _femcontext.get_elem_jacobian() << std::endl;
                  analytic_jacobian.add(1.0, _femcontext.get_elem_jacobian());
	          libMesh::out << "J_numeric " << _femcontext.get_elem().id() << " = "
                               << analytic_jacobian << std::endl;
                  libMesh::out.precision(old_precision);

                  libmesh_error();
                }
              // Once we've verified a side, we'll want to add back the
              // rest of the accumulated jacobian
              _femcontext.get_elem_jacobian() += old_jacobian;
            }
	  // In DEBUG mode, we've set elem_jacobian == 0, and we
          // may still need to add the old jacobian back
#ifdef DEBUG
	  if (_get_jacobian && jacobian_computed &&
              _sys.verify_analytic_jacobians == 0.0)
            {
              _femcontext.get_elem_jacobian() += old_jacobian;
            }
#endif // ifdef DEBUG
        }
    }

  class AssemblyContributions
  {
  public:
    /**
     * constructor to set context
     */
    AssemblyContributions(FEMSystem &sys,
                          bool get_residual,
                          bool get_jacobian) :
      _sys(sys),
      _get_residual(get_residual),
      _get_jacobian(get_jacobian) {}

    /**
     * operator() for use with Threads::parallel_for().
     */
    void operator()(const ConstElemRange &range) const
    {
      AutoPtr<DiffContext> con = _sys.build_context();
      FEMContext &_femcontext = libmesh_cast_ref<FEMContext&>(*con);
      _sys.init_context(_femcontext);

      for (ConstElemRange::const_iterator elem_it = range.begin();
           elem_it != range.end(); ++elem_it)
        {
          Elem *el = const_cast<Elem *>(*elem_it);

          _femcontext.pre_fe_reinit(_sys, el);
          _femcontext.elem_fe_reinit();

          assemble_unconstrained_element_system
            (_sys, _get_jacobian, _femcontext);

#ifdef LIBMESH_ENABLE_CONSTRAINTS
          // We turn off the asymmetric constraint application;
          // enforce_constraints_exactly() should be called in the solver
          if (_get_residual && _get_jacobian)
            _sys.get_dof_map().constrain_element_matrix_and_vector
              (_femcontext.get_elem_jacobian(), _femcontext.get_elem_residual(),
               _femcontext.get_dof_indices(), false);
          else if (_get_residual)
            _sys.get_dof_map().constrain_element_vector
              (_femcontext.get_elem_residual(), _femcontext.get_dof_indices(), false);
          else if (_get_jacobian)
            _sys.get_dof_map().constrain_element_matrix
              (_femcontext.get_elem_jacobian(), _femcontext.get_dof_indices(), false);
#endif // #ifdef LIBMESH_ENABLE_CONSTRAINTS

          if (_get_jacobian && _sys.print_element_jacobians)
            {
	      std::streamsize old_precision = libMesh::out.precision();
              libMesh::out.precision(16);
	      libMesh::out << "J_elem " << _femcontext.get_elem().id() << " = "
                        << _femcontext.get_elem_jacobian() << std::endl;
              libMesh::out.precision(old_precision);
            }

          { // A lock is necessary around access to the global system
            femsystem_mutex::scoped_lock lock(assembly_mutex);

            if (_get_jacobian)
              _sys.matrix->add_matrix (_femcontext.get_elem_jacobian(),
                                       _femcontext.get_dof_indices());
            if (_get_residual)
              _sys.rhs->add_vector (_femcontext.get_elem_residual(),
                                    _femcontext.get_dof_indices());
          } // Scope for assembly mutex

        }
    }

  private:

    FEMSystem& _sys;

    const bool _get_residual, _get_jacobian;
  };

  class PostprocessContributions
  {
  public:
    /**
     * constructor to set context
     */
    explicit
    PostprocessContributions(FEMSystem &sys) : _sys(sys) {}

    /**
     * operator() for use with Threads::parallel_for().
     */
    void operator()(const ConstElemRange &range) const
    {
      AutoPtr<DiffContext> con = _sys.build_context();
      FEMContext &_femcontext = libmesh_cast_ref<FEMContext&>(*con);
      _sys.init_context(_femcontext);

      for (ConstElemRange::const_iterator elem_it = range.begin();
           elem_it != range.end(); ++elem_it)
        {
          Elem *el = const_cast<Elem *>(*elem_it);
          _femcontext.pre_fe_reinit(_sys, el);

          // Optionally initialize all the interior FE objects on elem.
          if (_sys.fe_reinit_during_postprocess)
            _femcontext.elem_fe_reinit();

          _sys.element_postprocess(_femcontext);

          for (_femcontext.side = 0;
               _femcontext.side != _femcontext.get_elem().n_sides();
               ++_femcontext.side)
            {
              // Don't compute on non-boundary sides unless requested
              if (!_sys.postprocess_sides ||
                  (!_sys.get_physics()->compute_internal_sides &&
                   _femcontext.get_elem().neighbor(_femcontext.side) != NULL))
                continue;

              // Optionally initialize all the FE objects on this side.
              if (_sys.fe_reinit_during_postprocess)
                _femcontext.side_fe_reinit();

              _sys.side_postprocess(_femcontext);
            }
        }
    }

  private:

    FEMSystem& _sys;
  };

  class QoIContributions
  {
  public:
    /**
     * constructor to set context
     */
    explicit
    QoIContributions(FEMSystem &sys, DifferentiableQoI &diff_qoi,
		     const QoISet &qoi_indices) :
      qoi(sys.qoi.size(), 0.), _sys(sys), _diff_qoi(diff_qoi),_qoi_indices(qoi_indices) {}

    /**
     * splitting constructor
     */
    QoIContributions(const QoIContributions &other,
                     Threads::split) :
      qoi(other._sys.qoi.size(), 0.), _sys(other._sys), _diff_qoi(other._diff_qoi) {}

    /**
     * operator() for use with Threads::parallel_reduce().
     */
    void operator()(const ConstElemRange &range)
    {
      AutoPtr<DiffContext> con = _sys.build_context();
      FEMContext &_femcontext = libmesh_cast_ref<FEMContext&>(*con);
      _diff_qoi.init_context(_femcontext);

      for (ConstElemRange::const_iterator elem_it = range.begin();
           elem_it != range.end(); ++elem_it)
        {
          Elem *el = const_cast<Elem *>(*elem_it);

          _femcontext.pre_fe_reinit(_sys, el);

          if (_diff_qoi.assemble_qoi_elements)
            {
              _femcontext.elem_fe_reinit();

              _diff_qoi.element_qoi(_femcontext, _qoi_indices);
            }

          for (_femcontext.side = 0;
               _femcontext.side != _femcontext.get_elem().n_sides();
               ++_femcontext.side)
            {
              // Don't compute on non-boundary sides unless requested
              if (!_diff_qoi.assemble_qoi_sides ||
                  (!_diff_qoi.assemble_qoi_internal_sides &&
                   _femcontext.get_elem().neighbor(_femcontext.side) != NULL))
                continue;

              _femcontext.side_fe_reinit();

              _diff_qoi.side_qoi(_femcontext, _qoi_indices);
            }
        }

      this->_diff_qoi.thread_join( this->qoi, _femcontext.get_qois(), _qoi_indices );
    }

    void join (const QoIContributions& other)
    {
      libmesh_assert_equal_to (this->qoi.size(), other.qoi.size());
      this->_diff_qoi.thread_join( this->qoi, other.qoi, _qoi_indices );
    }

    std::vector<Number> qoi;

  private:

    FEMSystem& _sys;
    DifferentiableQoI& _diff_qoi;

    const QoISet _qoi_indices;
  };

  class QoIDerivativeContributions
  {
  public:
    /**
     * constructor to set context
     */
    QoIDerivativeContributions(FEMSystem &sys, const QoISet& qoi_indices,
			       DifferentiableQoI &qoi ) :
      _sys(sys), _qoi_indices(qoi_indices), _qoi(qoi) {}

    /**
     * operator() for use with Threads::parallel_for().
     */
    void operator()(const ConstElemRange &range) const
    {
      AutoPtr<DiffContext> con = _sys.build_context();
      FEMContext &_femcontext = libmesh_cast_ref<FEMContext&>(*con);
      _qoi.init_context(_femcontext);

      for (ConstElemRange::const_iterator elem_it = range.begin();
           elem_it != range.end(); ++elem_it)
        {
          Elem *el = const_cast<Elem *>(*elem_it);

          _femcontext.pre_fe_reinit(_sys, el);

          if (_qoi.assemble_qoi_elements)
            {
              _femcontext.elem_fe_reinit();

              _qoi.element_qoi_derivative(_femcontext, _qoi_indices);
            }

          for (_femcontext.side = 0;
               _femcontext.side != _femcontext.get_elem().n_sides();
               ++_femcontext.side)
            {
              // Don't compute on non-boundary sides unless requested
              if (!_qoi.assemble_qoi_sides ||
                  (!_qoi.assemble_qoi_internal_sides &&
                   _femcontext.get_elem().neighbor(_femcontext.side) != NULL))
                continue;

              _femcontext.side_fe_reinit();

              _qoi.side_qoi_derivative(_femcontext, _qoi_indices);
            }

          // We need some unmodified indices to use for constraining
          // multiple vector
          // FIXME - there should be a DofMap::constrain_element_vectors
          // to do this more efficiently
          std::vector<dof_id_type> original_dofs = _femcontext.get_dof_indices();

          { // A lock is necessary around access to the global system
            femsystem_mutex::scoped_lock lock(assembly_mutex);

#ifdef LIBMESH_ENABLE_CONSTRAINTS
            // We'll need to see if any heterogenous constraints apply
            // to the QoI dofs on this element *or* to any of the dofs
            // they depend on, so let's get those dependencies
            _sys.get_dof_map().constrain_nothing(_femcontext.get_dof_indices());
#endif

            for (unsigned int i=0; i != _sys.qoi.size(); ++i)
              if (_qoi_indices.has_index(i))
                {
#ifdef LIBMESH_ENABLE_CONSTRAINTS
                  bool has_heterogenous_constraint = false;
                  for (unsigned int d=0; 
                       d != _femcontext.get_dof_indices().size(); ++d)
                    if (_sys.get_dof_map().has_heterogenous_adjoint_constraint
                        (i, _femcontext.get_dof_indices()[d]))
                      {
                        has_heterogenous_constraint = true;
                        break;
                      }

                  _femcontext.get_dof_indices() = original_dofs;

                  // If we're going to need K to impose a heterogenous
                  // constraint then we either already have it or we
                  // need to compute it
                  if (has_heterogenous_constraint)
                    {
                      assemble_unconstrained_element_system
                        (_sys, true, _femcontext);

                      _sys.get_dof_map().heterogenously_constrain_element_vector
                        (_femcontext.get_elem_jacobian(),
                         _femcontext.get_qoi_derivatives()[i],
                         _femcontext.get_dof_indices(), false, i);
                    }
                  else
                    {
                      _sys.get_dof_map().constrain_element_vector
                        (_femcontext.get_qoi_derivatives()[i],
                         _femcontext.get_dof_indices(), false);
                    }
#endif

                  _sys.get_adjoint_rhs(i).add_vector
                    (_femcontext.get_qoi_derivatives()[i], _femcontext.get_dof_indices());
                }
          }
        }
    }

  private:

    FEMSystem& _sys;
    const QoISet& _qoi_indices;
    DifferentiableQoI& _qoi;
  };


}


namespace libMesh
{





FEMSystem::FEMSystem (EquationSystems& es,
                      const std::string& name_in,
                      const unsigned int number_in)
  : Parent(es, name_in, number_in),
    fe_reinit_during_postprocess(true),
    numerical_jacobian_h(TOLERANCE),
    verify_analytic_jacobians(0.0)
{
}


FEMSystem::~FEMSystem ()
{
  this->clear();
}



void FEMSystem::clear()
{
  Parent::clear();
}



void FEMSystem::init_data ()
{
  // First initialize LinearImplicitSystem data
  Parent::init_data();
}


void FEMSystem::assembly (bool get_residual, bool get_jacobian)
{
  libmesh_assert(get_residual || get_jacobian);
  std::string log_name;
  if (get_residual && get_jacobian)
    log_name = "assembly()";
  else if (get_residual)
    log_name = "assembly(get_residual)";
  else
    log_name = "assembly(get_jacobian)";

  START_LOG(log_name, "FEMSystem");

  const MeshBase& mesh = this->get_mesh();

//  this->get_vector("_nonlinear_solution").localize
//    (*current_local_nonlinear_solution,
//     dof_map.get_send_list());
  this->update();

  if (print_solution_norms)
    {
//      this->get_vector("_nonlinear_solution").close();
      this->solution->close();

      std::streamsize old_precision = libMesh::out.precision();
      libMesh::out.precision(16);
      libMesh::out << "|U| = "
//                    << this->get_vector("_nonlinear_solution").l1_norm()
                    << this->solution->l1_norm()
                    << std::endl;
      libMesh::out.precision(old_precision);
    }
  if (print_solutions)
    {
      std::streamsize old_precision = libMesh::out.precision();
      libMesh::out.precision(16);
//      libMesh::out << "U = [" << this->get_vector("_nonlinear_solution")
      libMesh::out << "U = [" << *(this->solution)
                    << "];" << std::endl;
      libMesh::out.precision(old_precision);
    }

  // Is this definitely necessary? [RHS]
  // Yes. [RHS 2012]
  if (get_jacobian)
    matrix->zero();
  if (get_residual)
    rhs->zero();

  // Stupid C++ lets you set *Real* verify_analytic_jacobians = true!
  if (verify_analytic_jacobians > 0.5)
    {
      libMesh::err << "WARNING!  verify_analytic_jacobians was set "
                    << "to absurdly large value of "
                    << verify_analytic_jacobians << std::endl;
      libMesh::err << "Resetting to 1e-6!" << std::endl;
      verify_analytic_jacobians = 1e-6;
    }

  // In time-dependent problems, the nonlinear function we're trying
  // to solve at each timestep may depend on the particular solver
  // we're using
  libmesh_assert(time_solver.get());

  // Build the residual and jacobian contributions on every active
  // mesh element on this processor
  Threads::parallel_for(elem_range.reset(mesh.active_local_elements_begin(),
                                         mesh.active_local_elements_end()),
                        AssemblyContributions(*this, get_residual, get_jacobian));


  if (get_residual && (print_residual_norms || print_residuals))
    this->rhs->close();
  if (get_residual && print_residual_norms)
    {
      std::streamsize old_precision = libMesh::out.precision();
      libMesh::out.precision(16);
      libMesh::out << "|F| = " << this->rhs->l1_norm() << std::endl;
      libMesh::out.precision(old_precision);
    }
  if (get_residual && print_residuals)
    {
      std::streamsize old_precision = libMesh::out.precision();
      libMesh::out.precision(16);
      libMesh::out << "F = [" << *(this->rhs) << "];" << std::endl;
      libMesh::out.precision(old_precision);
    }

  if (get_jacobian && (print_jacobian_norms || print_jacobians))
    this->matrix->close();
  if (get_jacobian && print_jacobian_norms)
    {
      std::streamsize old_precision = libMesh::out.precision();
      libMesh::out.precision(16);
      libMesh::out << "|J| = " << this->matrix->l1_norm() << std::endl;
      libMesh::out.precision(old_precision);
    }
  if (get_jacobian && print_jacobians)
    {
      std::streamsize old_precision = libMesh::out.precision();
      libMesh::out.precision(16);
      libMesh::out << "J = [" << *(this->matrix) << "];" << std::endl;
      libMesh::out.precision(old_precision);
    }
  STOP_LOG(log_name, "FEMSystem");
}



void FEMSystem::solve()
{
  // We are solving the primal problem
  Parent::solve();

  // On a moving mesh we want the mesh to reflect the new solution
  this->mesh_position_set();
}



void FEMSystem::mesh_position_set()
{
  // If we don't need to move the mesh, we're done
  if (_mesh_sys != this)
    return;

  MeshBase& mesh = this->get_mesh();

  AutoPtr<DiffContext> con = this->build_context();
  FEMContext &_femcontext = libmesh_cast_ref<FEMContext&>(*con);
  this->init_context(_femcontext);

  // Move every mesh element we can
  MeshBase::const_element_iterator el =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      // We need the algebraic data
      _femcontext.pre_fe_reinit(*this, *el);
      // And when asserts are on, we also need the FE so
      // we can assert that the mesh data is of the right type.
#ifndef NDEBUG
      _femcontext.elem_fe_reinit();
#endif

      // This code won't handle moving subactive elements
      libmesh_assert(!_femcontext.get_elem().has_children());

      _femcontext.elem_position_set(1.);
    }

  // We've now got positions set on all local nodes (and some
  // semilocal nodes); let's request positions for non-local nodes
  // from their processors.

  SyncNodalPositions sync_object(mesh);
  Parallel::sync_dofobject_data_by_id
    (this->comm(), mesh.nodes_begin(), mesh.nodes_end(), sync_object);
}



void FEMSystem::postprocess ()
{
  START_LOG("postprocess()", "FEMSystem");

  const MeshBase& mesh = this->get_mesh();

  this->update();

  // Get the time solver object associated with the system, and tell it that
  // we are not solving the adjoint problem
  this->get_time_solver().set_is_adjoint(false);

  // Loop over every active mesh element on this processor
  Threads::parallel_for(elem_range.reset(mesh.active_local_elements_begin(),
                                         mesh.active_local_elements_end()),
                        PostprocessContributions(*this));

  STOP_LOG("postprocess()", "FEMSystem");
}



void FEMSystem::assemble_qoi (const QoISet &qoi_indices)
{
  START_LOG("assemble_qoi()", "FEMSystem");

  const MeshBase& mesh = this->get_mesh();

  this->update();

  const unsigned int Nq = libmesh_cast_int<unsigned int>(qoi.size());

  // the quantity of interest is assumed to be a sum of element and
  // side terms
  for (unsigned int i=0; i != Nq; ++i)
    if (qoi_indices.has_index(i))
      qoi[i] = 0;

  // Create a non-temporary qoi_contributions object, so we can query
  // its results after the reduction
  QoIContributions qoi_contributions(*this, *(this->diff_qoi), qoi_indices);

  // Loop over every active mesh element on this processor
  Threads::parallel_reduce(elem_range.reset(mesh.active_local_elements_begin(),
                                            mesh.active_local_elements_end()),
                           qoi_contributions);

  this->diff_qoi->parallel_op( this->comm(), this->qoi, qoi_contributions.qoi, qoi_indices );

  STOP_LOG("assemble_qoi()", "FEMSystem");
}



void FEMSystem::assemble_qoi_derivative (const QoISet& qoi_indices)
{
  START_LOG("assemble_qoi_derivative()", "FEMSystem");

  const MeshBase& mesh = this->get_mesh();

  this->update();

  // The quantity of interest derivative assembly accumulates on
  // initially zero vectors
  for (unsigned int i=0; i != qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      this->add_adjoint_rhs(i).zero();

  // Loop over every active mesh element on this processor
  Threads::parallel_for(elem_range.reset(mesh.active_local_elements_begin(),
                                         mesh.active_local_elements_end()),
                        QoIDerivativeContributions(*this, qoi_indices,
						   *(this->diff_qoi)));

  STOP_LOG("assemble_qoi_derivative()", "FEMSystem");
}



void FEMSystem::numerical_jacobian (TimeSolverResPtr res,
                                    FEMContext &context) const
{
  // Logging is done by numerical_elem_jacobian
  // or numerical_side_jacobian

  DenseVector<Number> original_residual(context.get_elem_residual());
  DenseVector<Number> backwards_residual(context.get_elem_residual());
  DenseMatrix<Number> numeric_jacobian(context.get_elem_jacobian());
#ifdef DEBUG
  DenseMatrix<Number> old_jacobian(context.get_elem_jacobian());
#endif

  Real numerical_point_h = 0.;
  if (_mesh_sys == this)
    numerical_point_h = numerical_jacobian_h * context.get_elem().hmin();

  for (unsigned int j = 0; j != context.get_dof_indices().size(); ++j)
    {
      // Take the "minus" side of a central differenced first derivative
      Number original_solution = context.get_elem_solution()(j);
      context.get_elem_solution()(j) -= numerical_jacobian_h;

      // Make sure to catch any moving mesh terms
      // FIXME - this could be less ugly
      Real *coord = NULL;
      if (_mesh_sys == this)
        {
          if (_mesh_x_var != libMesh::invalid_uint)
            for (unsigned int k = 0;
                 k != context.get_dof_indices( _mesh_x_var ).size(); ++k)
              if (context.get_dof_indices( _mesh_x_var )[k] ==
                  context.get_dof_indices()[j])
                coord = &(context.get_elem().point(k)(0));
          if (_mesh_y_var != libMesh::invalid_uint)
            for (unsigned int k = 0;
                 k != context.get_dof_indices( _mesh_y_var ).size(); ++k)
              if (context.get_dof_indices( _mesh_y_var )[k] ==
                  context.get_dof_indices()[j])
                coord = &(context.get_elem().point(k)(1));
          if (_mesh_z_var != libMesh::invalid_uint)
            for (unsigned int k = 0;
                 k != context.get_dof_indices( _mesh_z_var ).size(); ++k)
	      if (context.get_dof_indices( _mesh_z_var )[k] ==
                  context.get_dof_indices()[j])
                coord = &(context.get_elem().point(k)(2));
        }
      if (coord)
        {
          // We have enough information to scale the perturbations
          // here appropriately
          context.get_elem_solution()(j) = original_solution - numerical_point_h;
          *coord = libmesh_real(context.get_elem_solution()(j));
        }

      context.get_elem_residual().zero();
      ((*time_solver).*(res))(false, context);
#ifdef DEBUG
      libmesh_assert_equal_to (old_jacobian, context.get_elem_jacobian());
#endif
      backwards_residual = context.get_elem_residual();

      // Take the "plus" side of a central differenced first derivative
      context.get_elem_solution()(j) = original_solution + numerical_jacobian_h;
      if (coord)
        {
          context.get_elem_solution()(j) = original_solution + numerical_point_h;
          *coord = libmesh_real(context.get_elem_solution()(j));
        }
      context.get_elem_residual().zero();
      ((*time_solver).*(res))(false, context);
#ifdef DEBUG
      libmesh_assert_equal_to (old_jacobian, context.get_elem_jacobian());
#endif

      context.get_elem_solution()(j) = original_solution;
      if (coord)
        {
          *coord = libmesh_real(context.get_elem_solution()(j));
	  for (unsigned int i = 0; i != context.get_dof_indices().size(); ++i)
            {
              numeric_jacobian(i,j) =
                (context.get_elem_residual()(i) - backwards_residual(i)) /
                2. / numerical_point_h;
            }
        }
      else
        {
          for (unsigned int i = 0; i != context.get_dof_indices().size(); ++i)
            {
              numeric_jacobian(i,j) =
                (context.get_elem_residual()(i) - backwards_residual(i)) /
                2. / numerical_jacobian_h;
            }
        }
    }

  context.get_elem_residual() = original_residual;
  context.get_elem_jacobian() = numeric_jacobian;
}



void FEMSystem::numerical_elem_jacobian (FEMContext &context) const
{
  START_LOG("numerical_elem_jacobian()", "FEMSystem");
  this->numerical_jacobian(&TimeSolver::element_residual, context);
  STOP_LOG("numerical_elem_jacobian()", "FEMSystem");
}



void FEMSystem::numerical_side_jacobian (FEMContext &context) const
{
  START_LOG("numerical_side_jacobian()", "FEMSystem");
  this->numerical_jacobian(&TimeSolver::side_residual, context);
  STOP_LOG("numerical_side_jacobian()", "FEMSystem");
}



AutoPtr<DiffContext> FEMSystem::build_context ()
{
  FEMContext* fc = new FEMContext(*this);

  AutoPtr<DiffContext> ap(fc);

  DifferentiablePhysics* phys = this->get_physics();

  libmesh_assert (phys);

  // If we are solving a moving mesh problem, tell that to the Context
  fc->set_mesh_system(phys->get_mesh_system());
  fc->set_mesh_x_var(phys->get_mesh_x_var());
  fc->set_mesh_y_var(phys->get_mesh_y_var());
  fc->set_mesh_z_var(phys->get_mesh_z_var());

  ap->set_deltat_pointer( &deltat );

  // If we are solving the adjoint problem, tell that to the Context
  ap->is_adjoint() = this->get_time_solver().is_adjoint();

  return ap;
}



void FEMSystem::init_context(DiffContext &c)
{
  // Parent::init_context(c);  // may be a good idea in derived classes

  // Although we do this in DiffSystem::build_context() and
  // FEMSystem::build_context() as well, we do it here just to be
  // extra sure that the deltat pointer gets set.  Since the
  // intended behavior is for classes derived from FEMSystem to
  // call Parent::init_context() in their own init_context()
  // overloads, we can ensure that those classes get the correct
  // deltat pointers even if they have different build_context()
  // overloads.
  c.set_deltat_pointer ( &deltat );

  FEMContext &context = libmesh_cast_ref<FEMContext&>(c);

  // Make sure we're prepared to do mass integration
  for (unsigned int var = 0; var != this->n_vars(); ++var)
    if (this->get_physics()->is_time_evolving(var))
      {
        // Request shape functions based on FEType
        switch( FEInterface::field_type( this->variable_type(var) ) )
          {
          case( TYPE_SCALAR ):
            {
              FEBase* elem_fe = NULL;
              context.get_element_fe(var, elem_fe);
              elem_fe->get_JxW();
              elem_fe->get_phi();
            }
          break;
          case( TYPE_VECTOR ):
            {
              FEGenericBase<RealGradient>* elem_fe = NULL;
              context.get_element_fe(var, elem_fe);
              elem_fe->get_JxW();
              elem_fe->get_phi();
            }
          break;
          default:
            {
              libmesh_error();
            }
          }
      }
}



void FEMSystem::mesh_position_get()
{
  // This function makes no sense unless we've already picked out some
  // variable(s) to reflect mesh position coordinates
  if (!_mesh_sys)
    libmesh_error();

  // We currently assume mesh variables are in our own system
  if (_mesh_sys != this)
    libmesh_not_implemented();

  // Loop over every active mesh element on this processor
  const MeshBase& mesh = this->get_mesh();

  MeshBase::const_element_iterator el =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  AutoPtr<DiffContext> con = this->build_context();
  FEMContext &_femcontext = libmesh_cast_ref<FEMContext&>(*con);
  this->init_context(_femcontext);

  // Get the solution's mesh variables from every element
  for ( ; el != end_el; ++el)
    {
      _femcontext.pre_fe_reinit(*this, *el);

      _femcontext.elem_position_get();

      if (_mesh_x_var != libMesh::invalid_uint)
        this->solution->insert(_femcontext.get_elem_solution(_mesh_x_var),
                               _femcontext.get_dof_indices(_mesh_x_var) );
      if (_mesh_y_var != libMesh::invalid_uint)
        this->solution->insert(_femcontext.get_elem_solution(_mesh_y_var),
                               _femcontext.get_dof_indices(_mesh_y_var));
      if (_mesh_z_var != libMesh::invalid_uint)
        this->solution->insert(_femcontext.get_elem_solution(_mesh_z_var),
                               _femcontext.get_dof_indices(_mesh_z_var));
    }

  this->solution->close();

  // And make sure the current_local_solution is up to date too
  this->System::update();
}

} // namespace libMesh
