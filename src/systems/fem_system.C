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

void assemble_unconstrained_element_system(const FEMSystem & _sys,
                                           const bool _get_jacobian,
                                           const bool _constrain_heterogeneously,
                                           FEMContext & _femcontext)
{
  if (_sys.print_element_solutions)
    {
      std::streamsize old_precision = libMesh::out.precision();
      libMesh::out.precision(16);
      if (_femcontext.has_elem())
        libMesh::out << "U_elem " << _femcontext.get_elem().id();
      else
        libMesh::out << "U_scalar ";
      libMesh::out << " = " << _femcontext.get_elem_solution() << std::endl;

      if (_sys.use_fixed_solution)
        {
          if (_femcontext.has_elem())
            libMesh::out << "Ufixed_elem " << _femcontext.get_elem().id();
          else
            libMesh::out << "Ufixed_scalar ";
          libMesh::out << " = " << _femcontext.get_elem_fixed_solution() << std::endl;
          libMesh::out.precision(old_precision);
        }
    }

  // We need jacobians to do heterogeneous residual constraints
  const bool need_jacobian =
    _get_jacobian || _constrain_heterogeneously;

  bool jacobian_computed =
    _sys.time_solver->element_residual(need_jacobian, _femcontext);

  // Compute a numeric jacobian if we have to
  if (need_jacobian && !jacobian_computed)
    {
      // Make sure we didn't compute a jacobian and lie about it
      libmesh_assert_equal_to (_femcontext.get_elem_jacobian().l1_norm(), 0.0);
      // Logging of numerical jacobians is done separately
      _sys.numerical_elem_jacobian(_femcontext);
    }

  // Compute a numeric jacobian if we're asked to verify the
  // analytic jacobian we got
  if (need_jacobian && jacobian_computed &&
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

          libmesh_error_msg("Relative error too large, exiting!");
        }
    }

  for (_femcontext.side = 0;
       _femcontext.side != _femcontext.get_elem().n_sides();
       ++_femcontext.side)
    {
      // Don't compute on non-boundary sides unless requested
      if (!_sys.get_physics()->compute_internal_sides &&
          _femcontext.get_elem().neighbor_ptr(_femcontext.side) != libmesh_nullptr)
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
      if (_sys.verify_analytic_jacobians != 0.0 && need_jacobian)
#endif // ifndef DEBUG
        {
          old_jacobian = _femcontext.get_elem_jacobian();
          _femcontext.get_elem_jacobian().zero();
        }

      jacobian_computed =
        _sys.time_solver->side_residual(need_jacobian, _femcontext);

      // Compute a numeric jacobian if we have to
      if (need_jacobian && !jacobian_computed)
        {
          // If we have already backed up old_jacobian,
          // we can make sure side_residual didn't compute a
          // jacobian and lie about it.
          //
          // If we haven't, then we need to, to let
          // numerical_side_jacobian work.
          if (old_jacobian.m())
            libmesh_assert_equal_to (_femcontext.get_elem_jacobian().l1_norm(), 0.0);
          else
            {
              old_jacobian = _femcontext.get_elem_jacobian();
              _femcontext.get_elem_jacobian().zero();
            }

          // Logging of numerical jacobians is done separately
          _sys.numerical_side_jacobian(_femcontext);

          // Add back in element interior numerical Jacobian
          _femcontext.get_elem_jacobian() += old_jacobian;
        }

      // Compute a numeric jacobian if we're asked to verify the
      // analytic jacobian we got
      else if (need_jacobian && jacobian_computed &&
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

              libmesh_error_msg("Relative error too large, exiting!");
            }
          // Once we've verified a side, we'll want to add back the
          // rest of the accumulated jacobian
          _femcontext.get_elem_jacobian() += old_jacobian;
        }

      // In DEBUG mode, we've set elem_jacobian == 0, and we
      // may have yet to add the old jacobian back
#ifdef DEBUG
      else
        {
          _femcontext.get_elem_jacobian() += old_jacobian;
        }
#endif // ifdef DEBUG
    }
}

void add_element_system(const FEMSystem & _sys,
                        const bool _get_residual,
                        const bool _get_jacobian,
                        const bool _constrain_heterogeneously,
                        FEMContext & _femcontext)
{
#ifdef LIBMESH_ENABLE_CONSTRAINTS
  if (_get_residual && _sys.print_element_residuals)
    {
      std::streamsize old_precision = libMesh::out.precision();
      libMesh::out.precision(16);
      if (_femcontext.has_elem())
        libMesh::out << "Rraw_elem " << _femcontext.get_elem().id();
      else
        libMesh::out << "Rraw_scalar ";
      libMesh::out << " = " << _femcontext.get_elem_residual() << std::endl;
      libMesh::out.precision(old_precision);
    }

  // We turn off the asymmetric constraint application;
  // enforce_constraints_exactly() should be called in the solver
  if (_get_residual && _get_jacobian)
    {
      if (_constrain_heterogeneously)
        _sys.get_dof_map().heterogenously_constrain_element_matrix_and_vector
          (_femcontext.get_elem_jacobian(),
           _femcontext.get_elem_residual(),
           _femcontext.get_dof_indices(), false);
      else
        _sys.get_dof_map().constrain_element_matrix_and_vector
          (_femcontext.get_elem_jacobian(),
           _femcontext.get_elem_residual(),
           _femcontext.get_dof_indices(), false);
    }
  else if (_get_residual)
    {
      if (_constrain_heterogeneously)
        _sys.get_dof_map().heterogenously_constrain_element_vector
          (_femcontext.get_elem_jacobian(),
           _femcontext.get_elem_residual(),
           _femcontext.get_dof_indices(), false);
      else
        _sys.get_dof_map().constrain_element_vector
          (_femcontext.get_elem_residual(), _femcontext.get_dof_indices(), false);
    }
  else if (_get_jacobian)
    {
      // Heterogeneous and homogeneous constraints are the same on the
      // matrix
      _sys.get_dof_map().constrain_element_matrix
        (_femcontext.get_elem_jacobian(), _femcontext.get_dof_indices(), false);
    }
#endif // #ifdef LIBMESH_ENABLE_CONSTRAINTS

  if (_get_residual && _sys.print_element_residuals)
    {
      std::streamsize old_precision = libMesh::out.precision();
      libMesh::out.precision(16);
      if (_femcontext.has_elem())
        libMesh::out << "R_elem " << _femcontext.get_elem().id();
      else
        libMesh::out << "R_scalar ";
      libMesh::out << " = " << _femcontext.get_elem_residual() << std::endl;
      libMesh::out.precision(old_precision);
    }

  if (_get_jacobian && _sys.print_element_jacobians)
    {
      std::streamsize old_precision = libMesh::out.precision();
      libMesh::out.precision(16);
      if (_femcontext.has_elem())
        libMesh::out << "J_elem " << _femcontext.get_elem().id();
      else
        libMesh::out << "J_scalar ";
      libMesh::out << " = " << _femcontext.get_elem_jacobian() << std::endl;
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



class AssemblyContributions
{
public:
  /**
   * constructor to set context
   */
  AssemblyContributions(FEMSystem & sys,
                        bool get_residual,
                        bool get_jacobian,
                        bool constrain_heterogeneously) :
    _sys(sys),
    _get_residual(get_residual),
    _get_jacobian(get_jacobian),
    _constrain_heterogeneously(constrain_heterogeneously) {}

  /**
   * operator() for use with Threads::parallel_for().
   */
  void operator()(const ConstElemRange & range) const
  {
    UniquePtr<DiffContext> con = _sys.build_context();
    FEMContext & _femcontext = cast_ref<FEMContext &>(*con);
    _sys.init_context(_femcontext);

    for (ConstElemRange::const_iterator elem_it = range.begin();
         elem_it != range.end(); ++elem_it)
      {
        Elem * el = const_cast<Elem *>(*elem_it);

        _femcontext.pre_fe_reinit(_sys, el);
        _femcontext.elem_fe_reinit();

        assemble_unconstrained_element_system
          (_sys, _get_jacobian, _constrain_heterogeneously,
           _femcontext);

        add_element_system
          (_sys, _get_residual, _get_jacobian,
           _constrain_heterogeneously, _femcontext);
      }
  }

private:

  FEMSystem & _sys;

  const bool _get_residual, _get_jacobian, _constrain_heterogeneously;
};

class PostprocessContributions
{
public:
  /**
   * constructor to set context
   */
  explicit
  PostprocessContributions(FEMSystem & sys) : _sys(sys) {}

  /**
   * operator() for use with Threads::parallel_for().
   */
  void operator()(const ConstElemRange & range) const
  {
    UniquePtr<DiffContext> con = _sys.build_context();
    FEMContext & _femcontext = cast_ref<FEMContext &>(*con);
    _sys.init_context(_femcontext);

    for (ConstElemRange::const_iterator elem_it = range.begin();
         elem_it != range.end(); ++elem_it)
      {
        Elem * el = const_cast<Elem *>(*elem_it);
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
                 _femcontext.get_elem().neighbor_ptr(_femcontext.side) != libmesh_nullptr))
              continue;

            // Optionally initialize all the FE objects on this side.
            if (_sys.fe_reinit_during_postprocess)
              _femcontext.side_fe_reinit();

            _sys.side_postprocess(_femcontext);
          }
      }
  }

private:

  FEMSystem & _sys;
};

class QoIContributions
{
public:
  /**
   * constructor to set context
   */
  explicit
  QoIContributions(FEMSystem & sys,
                   DifferentiableQoI & diff_qoi,
                   const QoISet & qoi_indices) :
    qoi(sys.qoi.size(), 0.), _sys(sys), _diff_qoi(diff_qoi),_qoi_indices(qoi_indices) {}

  /**
   * splitting constructor
   */
  QoIContributions(const QoIContributions & other,
                   Threads::split) :
    qoi(other._sys.qoi.size(), 0.), _sys(other._sys), _diff_qoi(other._diff_qoi) {}

  /**
   * operator() for use with Threads::parallel_reduce().
   */
  void operator()(const ConstElemRange & range)
  {
    UniquePtr<DiffContext> con = _sys.build_context();
    FEMContext & _femcontext = cast_ref<FEMContext &>(*con);
    _diff_qoi.init_context(_femcontext);

    bool have_some_heterogenous_qoi_bc = false;
#ifdef LIBMESH_ENABLE_CONSTRAINTS
    std::vector<bool> have_heterogenous_qoi_bc(_sys.qoi.size(), false);
    for (std::size_t q=0; q != _sys.qoi.size(); ++q)
      if (_qoi_indices.has_index(q) &&
          _sys.get_dof_map().has_heterogenous_adjoint_constraints(q))
        {
          have_heterogenous_qoi_bc[q] = true;
          have_some_heterogenous_qoi_bc = true;
        }
#endif

    if (have_some_heterogenous_qoi_bc)
      _sys.init_context(_femcontext);

    for (ConstElemRange::const_iterator elem_it = range.begin();
         elem_it != range.end(); ++elem_it)
      {
        Elem * el = const_cast<Elem *>(*elem_it);

        _femcontext.pre_fe_reinit(_sys, el);

        // We might have some heterogenous dofs here; let's see for
        // certain
#ifdef LIBMESH_ENABLE_CONSTRAINTS
        bool elem_has_some_heterogenous_qoi_bc = false;
        std::vector<bool> elem_has_heterogenous_qoi_bc(_sys.qoi.size(), false);
        if (have_some_heterogenous_qoi_bc)
          {
            for (std::size_t q=0; q != _sys.qoi.size(); ++q)
              {
                if (have_heterogenous_qoi_bc[q])
                  {
                    for (std::size_t d=0; d != _femcontext.get_dof_indices().size(); ++d)
                      if (_sys.get_dof_map().has_heterogenous_adjoint_constraint
                          (q, _femcontext.get_dof_indices()[d]) != Number(0))
                        {
                          elem_has_some_heterogenous_qoi_bc = true;
                          elem_has_heterogenous_qoi_bc[q] = true;
                          break;
                        }
                  }
              }
          }
#endif

        if (_diff_qoi.assemble_qoi_elements ||
            elem_has_some_heterogenous_qoi_bc)
          _femcontext.elem_fe_reinit();

        if (_diff_qoi.assemble_qoi_elements)
          _diff_qoi.element_qoi(_femcontext, _qoi_indices);

        // If we have some heterogenous dofs here, those are
        // themselves part of a regularized flux QoI which the library
        // handles integrating
#ifdef LIBMESH_ENABLE_CONSTRAINTS
        if (elem_has_some_heterogenous_qoi_bc)
          {
            _sys.time_solver->element_residual(false, _femcontext);

            for (std::size_t q=0; q != _sys.qoi.size(); ++q)
              {
                if (elem_has_heterogenous_qoi_bc[q])
                  {
                    for (std::size_t d=0; d != _femcontext.get_dof_indices().size(); ++d)
                      this->qoi[q] -= _femcontext.get_elem_residual()(d) *
                        _sys.get_dof_map().has_heterogenous_adjoint_constraint(q, _femcontext.get_dof_indices()[d]);

                  }
              }
          }
#endif

        for (_femcontext.side = 0;
             _femcontext.side != _femcontext.get_elem().n_sides();
             ++_femcontext.side)
          {
            // Don't compute on non-boundary sides unless requested
            if (!_diff_qoi.assemble_qoi_sides ||
                (!_diff_qoi.assemble_qoi_internal_sides &&
                 _femcontext.get_elem().neighbor_ptr(_femcontext.side) != libmesh_nullptr))
              continue;

            _femcontext.side_fe_reinit();

            _diff_qoi.side_qoi(_femcontext, _qoi_indices);
          }
      }

    this->_diff_qoi.thread_join( this->qoi, _femcontext.get_qois(), _qoi_indices );
  }

  void join (const QoIContributions & other)
  {
    libmesh_assert_equal_to (this->qoi.size(), other.qoi.size());
    this->_diff_qoi.thread_join( this->qoi, other.qoi, _qoi_indices );
  }

  std::vector<Number> qoi;

private:

  FEMSystem & _sys;
  DifferentiableQoI & _diff_qoi;

  const QoISet _qoi_indices;
};

class QoIDerivativeContributions
{
public:
  /**
   * constructor to set context
   */
  QoIDerivativeContributions(FEMSystem & sys,
                             const QoISet & qoi_indices,
                             DifferentiableQoI & qoi,
                             bool include_liftfunc,
                             bool apply_constraints) :
    _sys(sys),
    _qoi_indices(qoi_indices),
    _qoi(qoi),
    _include_liftfunc(include_liftfunc),
    _apply_constraints(apply_constraints) {}

  /**
   * operator() for use with Threads::parallel_for().
   */
  void operator()(const ConstElemRange & range) const
  {
    UniquePtr<DiffContext> con = _sys.build_context();
    FEMContext & _femcontext = cast_ref<FEMContext &>(*con);
    _qoi.init_context(_femcontext);

    bool have_some_heterogenous_qoi_bc = false;
#ifdef LIBMESH_ENABLE_CONSTRAINTS
    std::vector<bool> have_heterogenous_qoi_bc(_sys.qoi.size(), false);
    if (_include_liftfunc || _apply_constraints)
      for (std::size_t q=0; q != _sys.qoi.size(); ++q)
        if (_qoi_indices.has_index(q) &&
            _sys.get_dof_map().has_heterogenous_adjoint_constraints(q))
          {
            have_heterogenous_qoi_bc[q] = true;
            have_some_heterogenous_qoi_bc = true;
          }
#endif

    if (have_some_heterogenous_qoi_bc)
      _sys.init_context(_femcontext);

    for (ConstElemRange::const_iterator elem_it = range.begin();
         elem_it != range.end(); ++elem_it)
      {
        Elem * el = const_cast<Elem *>(*elem_it);

        _femcontext.pre_fe_reinit(_sys, el);

        // We might have some heterogenous dofs here; let's see for
        // certain
#ifdef LIBMESH_ENABLE_CONSTRAINTS
        bool elem_has_some_heterogenous_qoi_bc = false;
        std::vector<bool> elem_has_heterogenous_qoi_bc(_sys.qoi.size(), false);
        if (have_some_heterogenous_qoi_bc)
          {
            for (std::size_t q=0; q != _sys.qoi.size(); ++q)
              {
                if (have_heterogenous_qoi_bc[q])
                  {
                    for (std::size_t d=0; d != _femcontext.get_dof_indices().size(); ++d)
                      if (_sys.get_dof_map().has_heterogenous_adjoint_constraint
                          (q, _femcontext.get_dof_indices()[d]) != Number(0))
                        {
                          elem_has_some_heterogenous_qoi_bc = true;
                          elem_has_heterogenous_qoi_bc[q] = true;
                          break;
                        }
                  }
              }
          }
#endif

        // If we're going to call a user integral, then we need FE
        // information to call element_qoi.
        // If we're going to evaluate lift-function-based components
        // of a QoI, then we need FE information to assemble the
        // element residual.
        if (_qoi.assemble_qoi_elements ||
            ((_include_liftfunc || _apply_constraints) &&
             elem_has_some_heterogenous_qoi_bc))
          _femcontext.elem_fe_reinit();

        if (_qoi.assemble_qoi_elements)
          _qoi.element_qoi_derivative(_femcontext, _qoi_indices);

#ifdef LIBMESH_ENABLE_CONSTRAINTS
        // If we need to use heterogenous dofs here, we need the
        // Jacobian either for the regularized flux QoI integration
        // and/or for constraint application.
        if ((_include_liftfunc || _apply_constraints) &&
            elem_has_some_heterogenous_qoi_bc)
          {
            bool jacobian_computed = _sys.time_solver->element_residual(true, _femcontext);

            // If we're using numerical jacobians, above wont compute them
            if (!jacobian_computed)
              {
                // Make sure we didn't compute a jacobian and lie about it
                libmesh_assert_equal_to (_femcontext.get_elem_jacobian().l1_norm(), 0.0);
                // Logging of numerical jacobians is done separately
                _sys.numerical_elem_jacobian(_femcontext);
              }
          }

        // If we have some heterogenous dofs here, those are
        // themselves part of a regularized flux QoI which the library
        // may handle integrating
        if (_include_liftfunc && elem_has_some_heterogenous_qoi_bc)
          {
            for (std::size_t q=0; q != _sys.qoi.size(); ++q)
              {
                if (elem_has_heterogenous_qoi_bc[q])
                  {
                    for (std::size_t i=0; i != _femcontext.get_dof_indices().size(); ++i)
                      {
                        Number liftfunc_val =
                          _sys.get_dof_map().has_heterogenous_adjoint_constraint(q, _femcontext.get_dof_indices()[i]);

                        if (liftfunc_val != Number(0))
                          {
                            for (std::size_t j=0; j != _femcontext.get_dof_indices().size(); ++j)
                              _femcontext.get_qoi_derivatives()[q](j) -=
                                _femcontext.get_elem_jacobian()(i,j) *
                                liftfunc_val;
                          }
                      }
                  }
              }
          }
#endif


        for (_femcontext.side = 0;
             _femcontext.side != _femcontext.get_elem().n_sides();
             ++_femcontext.side)
          {
            // Don't compute on non-boundary sides unless requested
            if (!_qoi.assemble_qoi_sides ||
                (!_qoi.assemble_qoi_internal_sides &&
                 _femcontext.get_elem().neighbor_ptr(_femcontext.side) != libmesh_nullptr))
              continue;

            _femcontext.side_fe_reinit();

            _qoi.side_qoi_derivative(_femcontext, _qoi_indices);
          }

        // We need some unmodified indices to use for constraining
        // multiple vector
        // FIXME - there should be a DofMap::constrain_element_vectors
        // to do this more efficiently
#ifdef LIBMESH_ENABLE_CONSTRAINTS
        std::vector<dof_id_type> original_dofs = _femcontext.get_dof_indices();
#endif

        { // A lock is necessary around access to the global system
          femsystem_mutex::scoped_lock lock(assembly_mutex);

#ifdef LIBMESH_ENABLE_CONSTRAINTS
          // We'll need to see if any heterogenous constraints apply
          // to the QoI dofs on this element *or* to any of the dofs
          // they depend on, so let's get those dependencies
          if (_apply_constraints)
            _sys.get_dof_map().constrain_nothing(_femcontext.get_dof_indices());
#endif

          for (std::size_t i=0; i != _sys.qoi.size(); ++i)
            if (_qoi_indices.has_index(i))
              {
#ifdef LIBMESH_ENABLE_CONSTRAINTS
                if (_apply_constraints)
                  {
#ifndef NDEBUG
                    bool has_heterogenous_constraint = false;
                    for (std::size_t d=0; d != _femcontext.get_dof_indices().size(); ++d)
                      if (_sys.get_dof_map().has_heterogenous_adjoint_constraint
                          (i, _femcontext.get_dof_indices()[d]) != Number(0))
                        {
                          has_heterogenous_constraint = true;
                          libmesh_assert(elem_has_heterogenous_qoi_bc[i]);
                          libmesh_assert(elem_has_some_heterogenous_qoi_bc);
                          break;
                        }
#else
                    bool has_heterogenous_constraint =
                      elem_has_heterogenous_qoi_bc[i];
#endif

                    _femcontext.get_dof_indices() = original_dofs;

                    if (has_heterogenous_constraint)
                      {
                        // Q_u gets used for *adjoint* solves, so we
                        // need K^T here.
                        DenseMatrix<Number> elem_jacobian_transpose;
                        _femcontext.get_elem_jacobian().get_transpose
                          (elem_jacobian_transpose);

                        _sys.get_dof_map().heterogenously_constrain_element_vector
                          (elem_jacobian_transpose,
                           _femcontext.get_qoi_derivatives()[i],
                           _femcontext.get_dof_indices(), false, i);
                      }
                    else
                      {
                        _sys.get_dof_map().constrain_element_vector
                          (_femcontext.get_qoi_derivatives()[i],
                           _femcontext.get_dof_indices(), false);
                      }
                  }
#endif

                _sys.get_adjoint_rhs(i).add_vector
                  (_femcontext.get_qoi_derivatives()[i], _femcontext.get_dof_indices());
              }
        }
      }
  }

private:

  FEMSystem & _sys;
  const QoISet & _qoi_indices;
  DifferentiableQoI & _qoi;
  bool _include_liftfunc, _apply_constraints;
};


}


namespace libMesh
{





FEMSystem::FEMSystem (EquationSystems & es,
                      const std::string & name_in,
                      const unsigned int number_in)
  : Parent(es, name_in, number_in),
    fe_reinit_during_postprocess(true),
    numerical_jacobian_h(TOLERANCE),
    verify_analytic_jacobians(0.0)
{
}


FEMSystem::~FEMSystem ()
{
}



void FEMSystem::init_data ()
{
  // First initialize LinearImplicitSystem data
  Parent::init_data();
}


void FEMSystem::assembly (bool get_residual, bool get_jacobian,
                          bool apply_heterogeneous_constraints)
{
  libmesh_assert(get_residual || get_jacobian);
  std::string log_name;
  if (get_residual && get_jacobian)
    log_name = "assembly()";
  else if (get_residual)
    log_name = "assembly(get_residual)";
  else
    log_name = "assembly(get_jacobian)";

  LOG_SCOPE(log_name.c_str(), "FEMSystem");

  const MeshBase & mesh = this->get_mesh();

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
  Threads::parallel_for
    (elem_range.reset(mesh.active_local_elements_begin(),
                      mesh.active_local_elements_end()),
     AssemblyContributions(*this, get_residual, get_jacobian,
                           apply_heterogeneous_constraints));

  // Check and see if we have SCALAR variables
  bool have_scalar = false;
  for(unsigned int i=0; i != this->n_variable_groups(); ++i)
    {
      if( this->variable_group(i).type().family == SCALAR )
        {
          have_scalar = true;
          break;
        }
    }

  // SCALAR dofs are stored on the last processor, so we'll evaluate
  // their equation terms there and only if we have a SCALAR variable
  if ( this->processor_id() == (this->n_processors()-1) && have_scalar )
    {
      UniquePtr<DiffContext> con = this->build_context();
      FEMContext & _femcontext = cast_ref<FEMContext &>(*con);
      this->init_context(_femcontext);
      _femcontext.pre_fe_reinit(*this, libmesh_nullptr);

      bool jacobian_computed =
        this->time_solver->nonlocal_residual(get_jacobian, _femcontext);

      // Nonlocal residuals are likely to be length 0, in which case we
      // don't need to do any more.  And we shouldn't try to do any
      // more; lots of DenseVector/DenseMatrix code assumes rank>0.
      if (_femcontext.get_elem_residual().size())
        {
          // Compute a numeric jacobian if we have to
          if (get_jacobian && !jacobian_computed)
            {
              // Make sure we didn't compute a jacobian and lie about it
              libmesh_assert_equal_to (_femcontext.get_elem_jacobian().l1_norm(), 0.0);
              // Logging of numerical jacobians is done separately
              this->numerical_nonlocal_jacobian(_femcontext);
            }

          // Compute a numeric jacobian if we're asked to verify the
          // analytic jacobian we got
          if (get_jacobian && jacobian_computed &&
              this->verify_analytic_jacobians != 0.0)
            {
              DenseMatrix<Number> analytic_jacobian(_femcontext.get_elem_jacobian());

              _femcontext.get_elem_jacobian().zero();
              // Logging of numerical jacobians is done separately
              this->numerical_nonlocal_jacobian(_femcontext);

              Real analytic_norm = analytic_jacobian.l1_norm();
              Real numerical_norm = _femcontext.get_elem_jacobian().l1_norm();

              // If we can continue, we'll probably prefer the analytic jacobian
              analytic_jacobian.swap(_femcontext.get_elem_jacobian());

              // The matrix "analytic_jacobian" will now hold the error matrix
              analytic_jacobian.add(-1.0, _femcontext.get_elem_jacobian());
              Real error_norm = analytic_jacobian.l1_norm();

              Real relative_error = error_norm /
                std::max(analytic_norm, numerical_norm);

              if (relative_error > this->verify_analytic_jacobians)
                {
                  libMesh::err << "Relative error " << relative_error
                               << " detected in analytic jacobian on nonlocal dofs!"
                               << std::endl;

                  std::streamsize old_precision = libMesh::out.precision();
                  libMesh::out.precision(16);
                  libMesh::out << "J_analytic nonlocal = "
                               << _femcontext.get_elem_jacobian() << std::endl;
                  analytic_jacobian.add(1.0, _femcontext.get_elem_jacobian());
                  libMesh::out << "J_numeric nonlocal = "
                               << analytic_jacobian << std::endl;

                  libMesh::out.precision(old_precision);

                  libmesh_error_msg("Relative error too large, exiting!");
                }
            }

          add_element_system
            (*this, get_residual, get_jacobian,
             apply_heterogeneous_constraints, _femcontext);
        }
    }

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

  MeshBase & mesh = this->get_mesh();

  UniquePtr<DiffContext> con = this->build_context();
  FEMContext & _femcontext = cast_ref<FEMContext &>(*con);
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
  LOG_SCOPE("postprocess()", "FEMSystem");

  const MeshBase & mesh = this->get_mesh();

  this->update();

  // Get the time solver object associated with the system, and tell it that
  // we are not solving the adjoint problem
  this->get_time_solver().set_is_adjoint(false);

  // Loop over every active mesh element on this processor
  Threads::parallel_for(elem_range.reset(mesh.active_local_elements_begin(),
                                         mesh.active_local_elements_end()),
                        PostprocessContributions(*this));
}



void FEMSystem::assemble_qoi (const QoISet & qoi_indices)
{
  LOG_SCOPE("assemble_qoi()", "FEMSystem");

  const MeshBase & mesh = this->get_mesh();

  this->update();

  const unsigned int Nq = cast_int<unsigned int>(qoi.size());

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
}



void FEMSystem::assemble_qoi_derivative (const QoISet & qoi_indices,
                                         bool include_liftfunc,
                                         bool apply_constraints)
{
  LOG_SCOPE("assemble_qoi_derivative()", "FEMSystem");

  const MeshBase & mesh = this->get_mesh();

  this->update();

  // The quantity of interest derivative assembly accumulates on
  // initially zero vectors
  for (std::size_t i=0; i != qoi.size(); ++i)
    if (qoi_indices.has_index(i))
      this->add_adjoint_rhs(i).zero();

  // Loop over every active mesh element on this processor
  Threads::parallel_for(elem_range.reset(mesh.active_local_elements_begin(),
                                         mesh.active_local_elements_end()),
                        QoIDerivativeContributions(*this, qoi_indices,
                                                   *(this->diff_qoi),
                                                   include_liftfunc,
                                                   apply_constraints));
}



void FEMSystem::numerical_jacobian (TimeSolverResPtr res,
                                    FEMContext & context) const
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

  for (unsigned int v = 0; v != context.n_vars(); ++v)
    {
      const Real my_h = this->numerical_jacobian_h_for_var(v);

      unsigned int j_offset = libMesh::invalid_uint;

      if (!context.get_dof_indices(v).empty())
        {
          for (std::size_t i = 0; i != context.get_dof_indices().size(); ++i)
            if (context.get_dof_indices()[i] ==
                context.get_dof_indices(v)[0])
              j_offset = i;

          libmesh_assert_not_equal_to(j_offset, libMesh::invalid_uint);
        }

      for (std::size_t j = 0; j != context.get_dof_indices(v).size(); ++j)
        {
          const unsigned int total_j = j + j_offset;

          // Take the "minus" side of a central differenced first derivative
          Number original_solution = context.get_elem_solution(v)(j);
          context.get_elem_solution(v)(j) -= my_h;

          // Make sure to catch any moving mesh terms
          Real * coord = libmesh_nullptr;
          if (_mesh_sys == this)
            {
              if (_mesh_x_var == v)
                coord = &(context.get_elem().point(j)(0));
              else if (_mesh_y_var == v)
                coord = &(context.get_elem().point(j)(1));
              else if (_mesh_z_var == v)
                coord = &(context.get_elem().point(j)(2));
            }
          if (coord)
            {
              // We have enough information to scale the perturbations
              // here appropriately
              context.get_elem_solution(v)(j) = original_solution - numerical_point_h;
              *coord = libmesh_real(context.get_elem_solution(v)(j));
            }

          context.get_elem_residual().zero();
          ((*time_solver).*(res))(false, context);
#ifdef DEBUG
          libmesh_assert_equal_to (old_jacobian, context.get_elem_jacobian());
#endif
          backwards_residual = context.get_elem_residual();

          // Take the "plus" side of a central differenced first derivative
          context.get_elem_solution(v)(j) = original_solution + my_h;
          if (coord)
            {
              context.get_elem_solution()(j) = original_solution + numerical_point_h;
              *coord = libmesh_real(context.get_elem_solution(v)(j));
            }
          context.get_elem_residual().zero();
          ((*time_solver).*(res))(false, context);
#ifdef DEBUG
          libmesh_assert_equal_to (old_jacobian, context.get_elem_jacobian());
#endif

          context.get_elem_solution(v)(j) = original_solution;
          if (coord)
            {
              *coord = libmesh_real(context.get_elem_solution(v)(j));
              for (std::size_t i = 0; i != context.get_dof_indices().size(); ++i)
                {
                  numeric_jacobian(i,total_j) =
                    (context.get_elem_residual()(i) - backwards_residual(i)) /
                    2. / numerical_point_h;
                }
            }
          else
            {
              for (std::size_t i = 0; i != context.get_dof_indices().size(); ++i)
                {
                  numeric_jacobian(i,total_j) =
                    (context.get_elem_residual()(i) - backwards_residual(i)) /
                    2. / my_h;
                }
            }
        }
    }

  context.get_elem_residual() = original_residual;
  context.get_elem_jacobian() = numeric_jacobian;
}



void FEMSystem::numerical_elem_jacobian (FEMContext & context) const
{
  LOG_SCOPE("numerical_elem_jacobian()", "FEMSystem");
  this->numerical_jacobian(&TimeSolver::element_residual, context);
}



void FEMSystem::numerical_side_jacobian (FEMContext & context) const
{
  LOG_SCOPE("numerical_side_jacobian()", "FEMSystem");
  this->numerical_jacobian(&TimeSolver::side_residual, context);
}



void FEMSystem::numerical_nonlocal_jacobian (FEMContext & context) const
{
  LOG_SCOPE("numerical_nonlocal_jacobian()", "FEMSystem");
  this->numerical_jacobian(&TimeSolver::nonlocal_residual, context);
}



UniquePtr<DiffContext> FEMSystem::build_context ()
{
  FEMContext * fc = new FEMContext(*this);

  DifferentiablePhysics * phys = this->get_physics();

  libmesh_assert (phys);

  // If we are solving a moving mesh problem, tell that to the Context
  fc->set_mesh_system(phys->get_mesh_system());
  fc->set_mesh_x_var(phys->get_mesh_x_var());
  fc->set_mesh_y_var(phys->get_mesh_y_var());
  fc->set_mesh_z_var(phys->get_mesh_z_var());

  fc->set_deltat_pointer( &deltat );

  // If we are solving the adjoint problem, tell that to the Context
  fc->is_adjoint() = this->get_time_solver().is_adjoint();

  return UniquePtr<DiffContext>(fc);
}



void FEMSystem::init_context(DiffContext & c)
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

  FEMContext & context = cast_ref<FEMContext &>(c);

  // Make sure we're prepared to do mass integration
  for (unsigned int var = 0; var != this->n_vars(); ++var)
    if (this->get_physics()->is_time_evolving(var))
      {
        // Request shape functions based on FEType
        switch( FEInterface::field_type( this->variable_type(var) ) )
          {
          case( TYPE_SCALAR ):
            {
              FEBase * elem_fe = libmesh_nullptr;
              context.get_element_fe(var, elem_fe);
              elem_fe->get_JxW();
              elem_fe->get_phi();
            }
            break;
          case( TYPE_VECTOR ):
            {
              FEGenericBase<RealGradient> * elem_fe = libmesh_nullptr;
              context.get_element_fe(var, elem_fe);
              elem_fe->get_JxW();
              elem_fe->get_phi();
            }
            break;
          default:
            libmesh_error_msg("Unrecognized field type!");
          }
      }
}



void FEMSystem::mesh_position_get()
{
  // This function makes no sense unless we've already picked out some
  // variable(s) to reflect mesh position coordinates
  if (!_mesh_sys)
    libmesh_error_msg("_mesh_sys was NULL!");

  // We currently assume mesh variables are in our own system
  if (_mesh_sys != this)
    libmesh_not_implemented();

  // Loop over every active mesh element on this processor
  const MeshBase & mesh = this->get_mesh();

  MeshBase::const_element_iterator el =
    mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
    mesh.active_local_elements_end();

  UniquePtr<DiffContext> con = this->build_context();
  FEMContext & _femcontext = cast_ref<FEMContext &>(*con);
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
