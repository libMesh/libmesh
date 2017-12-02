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



#ifndef GENERIC_PROJECTOR_H
#define GENERIC_PROJECTOR_H

// C++ includes
#include <vector>

// Local includes
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature.h"
#include "libmesh/system.h"
#include "libmesh/threads.h"

namespace libMesh
{

/**
 * The GenericProjector class implements the core of other projection
 * operations, using two input functors to read values to be projected
 * and an output functor to set degrees of freedom in the result.
 *
 * This may be executed in parallel on multiple threads.
 *
 * \author Roy H. Stogner
 * \date 2016
 */
template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
class GenericProjector
{
private:
  const System & system;

  // For TBB compatibility and thread safety we'll copy these in
  // operator()
  const FFunctor & master_f;
  const GFunctor * master_g;  // Needed for C1 type elements only
  bool g_was_copied;
  const ProjectionAction & master_action;
  const std::vector<unsigned int> & variables;

public:
  GenericProjector (const System & system_in,
                    const FFunctor & f_in,
                    const GFunctor * g_in,
                    const ProjectionAction & act_in,
                    const std::vector<unsigned int> & variables_in) :
    system(system_in),
    master_f(f_in),
    master_g(g_in),
    g_was_copied(false),
    master_action(act_in),
    variables(variables_in)
  {}

  GenericProjector (const GenericProjector & in) :
    system(in.system),
    master_f(in.master_f),
    master_g(in.master_g ? new GFunctor(*in.master_g) : libmesh_nullptr),
    g_was_copied(in.master_g),
    master_action(in.master_action),
    variables(in.variables)
  {}

  ~GenericProjector()
  {
    if (g_was_copied)
      delete master_g;
  }

  void operator() (const ConstElemRange & range) const;
};


/**
 * The VectorSetAction output functor class can be used with a
 * GenericProjector to set projection values (which must be of type
 * Val) as coefficients of the given NumericVector.
 *
 * \author Roy H. Stogner
 * \date 2016
 */

template <typename Val>
class VectorSetAction
{
private:
  NumericVector<Val> & target_vector;

public:
  VectorSetAction(NumericVector<Val> & target_vec) :
    target_vector(target_vec) {}

  void insert(const FEMContext & c,
              unsigned int var_num,
              const DenseVector<Val> & Ue)
  {
    const numeric_index_type
      first = target_vector.first_local_index(),
      last  = target_vector.last_local_index();

    const std::vector<dof_id_type> & dof_indices =
      c.get_dof_indices(var_num);

    unsigned int size = Ue.size();

    libmesh_assert_equal_to(size, dof_indices.size());

    // Lock the new vector since it is shared among threads.
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

      for (unsigned int i = 0; i != size; ++i)
        if ((dof_indices[i] >= first) && (dof_indices[i] <  last))
          target_vector.set(dof_indices[i], Ue(i));
    }
  }
};


/**
 * The FEMFunctionWrapper input functor class can be used with a
 * GenericProjector to read values from an FEMFunction.
 *
 * \author Roy H. Stogner
 * \date 2016
 */

template <typename Output>
class FEMFunctionWrapper
{
public:
  FEMFunctionWrapper(const FEMFunctionBase<Output> & f) : _f(f.clone()) {}

  FEMFunctionWrapper(const FEMFunctionWrapper<Output> & fw) :
    _f(fw._f->clone()) {}

  void init_context (FEMContext & c) { _f->init_context(c); }

  Output eval_at_node (const FEMContext & c,
                       unsigned int i,
                       unsigned int /*elem_dim*/,
                       const Node & n,
                       const Real time)
  { return _f->component(c, i, n, time); }

  Output eval_at_point (const FEMContext & c,
                        unsigned int i,
                        const Point & n,
                        const Real time)
  { return _f->component(c, i, n, time); }

  bool is_grid_projection() { return false; }

  void eval_old_dofs (const FEMContext & /* c */,
                      unsigned int /* var_component */,
                      std::vector<Output> /* values */)
  { libmesh_error(); }

private:
  std::unique_ptr<FEMFunctionBase<Output>> _f;
};


/**
 * The OldSolutionBase input functor abstract base class is the root
 * of the OldSolutionValue and OldSolutionCoefs classes which allow a
 * GenericProjector to read old solution values or solution
 * interpolation coefficients for a just-refined-and-coarsened mesh.
 *
 * \author Roy H. Stogner
 * \date 2016
 */

#ifdef LIBMESH_ENABLE_AMR
template <typename Output,
          void (FEMContext::*point_output) (unsigned int,
                                            const Point &,
                                            Output &,
                                            const Real) const>
class OldSolutionBase
{
public:
  OldSolutionBase(const libMesh::System & sys_in) :
    last_elem(libmesh_nullptr),
    sys(sys_in),
    old_context(sys_in)
  {
  }

  OldSolutionBase(const OldSolutionBase & in) :
    last_elem(libmesh_nullptr),
    sys(in.sys),
    old_context(sys)
  {
  }

  static void get_shape_outputs(FEBase & fe);

  // Integrating on new mesh elements, we won't yet have an up to date
  // current_local_solution.
  void init_context (FEMContext & c)
  {
    c.set_algebraic_type(FEMContext::DOFS_ONLY);

    // Loop over variables, to prerequest
    for (unsigned int var=0; var!=sys.n_vars(); ++var)
      {
        FEBase * fe = libmesh_nullptr;
        const std::set<unsigned char> & elem_dims =
          old_context.elem_dimensions();

        for (std::set<unsigned char>::const_iterator dim_it =
               elem_dims.begin(); dim_it != elem_dims.end(); ++dim_it)
          {
            const unsigned char dim = *dim_it;
            old_context.get_element_fe( var, fe, dim );
            get_shape_outputs(*fe);
          }
      }
  }

  bool is_grid_projection() { return true; }

protected:
  void check_old_context (const FEMContext & c)
  {
    LOG_SCOPE ("check_old_context(c)", "OldSolutionBase");
    const Elem & elem = c.get_elem();
    if (last_elem != &elem)
      {
        if (elem.refinement_flag() == Elem::JUST_REFINED)
          {
            old_context.pre_fe_reinit(sys, elem.parent());
          }
        else if (elem.refinement_flag() == Elem::JUST_COARSENED)
          {
            libmesh_error();
          }
        else
          {
            if (!elem.old_dof_object)
              {
                libmesh_error();
              }

            old_context.pre_fe_reinit(sys, &elem);
          }

        last_elem = &elem;
      }
    else
      {
        libmesh_assert(old_context.has_elem());
      }
  }


  bool check_old_context (const FEMContext & c, const Point & p)
  {
    LOG_SCOPE ("check_old_context(c,p)", "OldSolutionBase");
    const Elem & elem = c.get_elem();
    if (last_elem != &elem)
      {
        if (elem.refinement_flag() == Elem::JUST_REFINED)
          {
            old_context.pre_fe_reinit(sys, elem.parent());
          }
        else if (elem.refinement_flag() == Elem::JUST_COARSENED)
          {
            // Find the child with this point.  Use out_of_elem_tol
            // (in physical space, which may correspond to a large
            // tolerance in master space!) to allow for out-of-element
            // finite differencing of mixed gradient terms.  Pray we
            // have no quadrature locations which are within 1e-5 of
            // the element subdivision boundary but are not exactly on
            // that boundary.
            const Real master_tol = out_of_elem_tol / elem.hmax() * 2;

            for (auto & child : elem.child_ref_range())
              if (child.close_to_point(p, master_tol))
                {
                  old_context.pre_fe_reinit(sys, &child);
                  break;
                }

            libmesh_assert
              (old_context.get_elem().close_to_point(p, master_tol));
          }
        else
          {
            if (!elem.old_dof_object)
              return false;

            old_context.pre_fe_reinit(sys, &elem);
          }

        last_elem = &elem;
      }
    else
      {
        libmesh_assert(old_context.has_elem());

        const Real master_tol = out_of_elem_tol / elem.hmax() * 2;

        if (!old_context.get_elem().close_to_point(p, master_tol))
          {
            libmesh_assert_equal_to
              (elem.refinement_flag(), Elem::JUST_COARSENED);

            for (auto & child : elem.child_ref_range())
              if (child.close_to_point(p, master_tol))
                {
                  old_context.pre_fe_reinit(sys, &child);
                  break;
                }

            libmesh_assert
              (old_context.get_elem().close_to_point(p, master_tol));
          }
      }

    return true;
  }

protected:
  const Elem * last_elem;
  const System & sys;
  FEMContext old_context;

  static const Real out_of_elem_tol;
};


/**
 * The OldSolutionValue input functor class can be used with
 * GenericProjector to read values from a solution on a
 * just-refined-and-coarsened mesh.
 *
 * \author Roy H. Stogner
 * \date 2016
 */
template <typename Output,
          void (FEMContext::*point_output) (unsigned int,
                                            const Point &,
                                            Output &,
                                            const Real) const>
class OldSolutionValue : public OldSolutionBase<Output, point_output>
{
public:
  OldSolutionValue(const libMesh::System & sys_in,
                   const NumericVector<Number> & old_sol) :
    OldSolutionBase<Output, point_output>(sys_in),
    old_solution(old_sol)
  {
    this->old_context.set_algebraic_type(FEMContext::OLD);
    this->old_context.set_custom_solution(&old_solution);
  }

  OldSolutionValue(const OldSolutionValue & in) :
    OldSolutionBase<Output, point_output>(in.sys),
    old_solution(in.old_solution)
  {
    this->old_context.set_algebraic_type(FEMContext::OLD);
    this->old_context.set_custom_solution(&old_solution);
  }


  Output eval_at_node (const FEMContext & c,
                       unsigned int i,
                       unsigned int elem_dim,
                       const Node & n,
                       Real /* time */ =0.);

  Output eval_at_point(const FEMContext & c,
                       unsigned int i,
                       const Point & p,
                       Real /* time */ =0.)
  {
    LOG_SCOPE ("eval_at_point()", "OldSolutionValue");

    if (!this->check_old_context(c, p))
      return 0;

    Output n;
    (this->old_context.*point_output)(i, p, n, this->out_of_elem_tol);
    return n;
  }

  void eval_old_dofs (const FEMContext & c,
                      unsigned int var,
                      std::vector<Output> & values)
  {
    LOG_SCOPE ("eval_old_dofs()", "OldSolutionValue");

    this->check_old_context(c);

    const std::vector<dof_id_type> & old_dof_indices =
      this->old_context.get_dof_indices(var);

    libmesh_assert_equal_to (old_dof_indices.size(), values.size());

    old_solution.get(old_dof_indices, values);
  }

private:
  const NumericVector<Number> & old_solution;
};


template<>
inline void
OldSolutionBase<Number, &FEMContext::point_value>::get_shape_outputs(FEBase & fe)
{
  fe.get_phi();
}


template<>
inline void
OldSolutionBase<Gradient, &FEMContext::point_gradient>::get_shape_outputs(FEBase & fe)
{
  fe.get_dphi();
}


#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template<>
inline void
OldSolutionBase<Real, &FEMContext::point_value>::get_shape_outputs(FEBase & fe)
{
  fe.get_phi();
}


template<>
inline void
OldSolutionBase<RealGradient, &FEMContext::point_gradient>::get_shape_outputs(FEBase & fe)
{
  fe.get_dphi();
}
#endif // LIBMESH_USE_COMPLEX_NUMBERS


template<>
inline
Number
OldSolutionValue<Number, &FEMContext::point_value>::
eval_at_node(const FEMContext & c,
             unsigned int i,
             unsigned int /* elem_dim */,
             const Node & n,
             Real /* time */)
{
  LOG_SCOPE ("Number eval_at_node()", "OldSolutionValue");

  // Optimize for the common case, where this node was part of the
  // old solution.
  //
  // Be sure to handle cases where the variable wasn't defined on
  // this node (due to changing subdomain support) or where the
  // variable has no components on this node (due to Elem order
  // exceeding FE order)
  if (n.old_dof_object &&
      n.old_dof_object->n_vars(sys.number()) &&
      n.old_dof_object->n_comp(sys.number(), i))
    {
      const dof_id_type old_id =
        n.old_dof_object->dof_number(sys.number(), i, 0);
      return old_solution(old_id);
    }

  return this->eval_at_point(c, i, n, 0);
}



template<>
inline
Gradient
OldSolutionValue<Gradient, &FEMContext::point_gradient>::
eval_at_node(const FEMContext & c,
             unsigned int i,
             unsigned int elem_dim,
             const Node & n,
             Real /* time */)
{
  LOG_SCOPE ("Gradient eval_at_node()", "OldSolutionValue");

  // Optimize for the common case, where this node was part of the
  // old solution.
  //
  // Be sure to handle cases where the variable wasn't defined on
  // this node (due to changing subdomain support) or where the
  // variable has no components on this node (due to Elem order
  // exceeding FE order)
  if (n.old_dof_object &&
      n.old_dof_object->n_vars(sys.number()) &&
      n.old_dof_object->n_comp(sys.number(), i))
    {
      Gradient g;
      for (unsigned int d = 0; d != elem_dim; ++d)
        {
          const dof_id_type old_id =
            n.old_dof_object->dof_number(sys.number(), i, d+1);
          g(d) = old_solution(old_id);
        }
      return g;
    }

  return this->eval_at_point(c, i, n, 0);
}





template <>
const Real OldSolutionBase <Number, &FEMContext::point_value>::out_of_elem_tol = 10*TOLERANCE;

template <>
const Real OldSolutionBase <Gradient, &FEMContext::point_gradient>::out_of_elem_tol = 10*TOLERANCE;

#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template <>
const Real OldSolutionBase <Real, &FEMContext::point_value>::out_of_elem_tol = 10*TOLERANCE;

template <>
const Real OldSolutionBase <RealGradient, &FEMContext::point_gradient>::out_of_elem_tol = 10*TOLERANCE;
#endif // LIBMESH_USE_COMPLEX_NUMBERS

#endif // LIBMESH_ENABLE_AMR

/**
 * Function definitions
 */


template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
void GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::operator()
  (const ConstElemRange & range) const
{
  LOG_SCOPE ("operator()", "GenericProjector");

  ProjectionAction action(master_action);
  FFunctor f(master_f);
  std::unique_ptr<GFunctor> g;
  if (master_g)
    g.reset(new GFunctor(*master_g));

  // The DofMap for this system
  const DofMap & dof_map = system.get_dof_map();

  // The element matrix and RHS for projections.
  // Note that Ke is always real-valued, whereas
  // Fe may be complex valued if complex number
  // support is enabled
  DenseMatrix<Real> Ke;
  DenseVector<FValue> Fe;
  // The new element degree of freedom coefficients
  DenseVector<FValue> Ue;

  // Context objects to contain all our required FE objects
  FEMContext context( system );

  // Loop over all the variables we've been requested to project, to
  // pre-request
  for (std::size_t v=0; v!=variables.size(); v++)
    {
      const unsigned int var = variables[v];

      // FIXME: Need to generalize this to vector-valued elements. [PB]
      FEBase * fe = libmesh_nullptr;
      FEBase * side_fe = libmesh_nullptr;
      FEBase * edge_fe = libmesh_nullptr;

      const std::set<unsigned char> & elem_dims =
        context.elem_dimensions();

      for (std::set<unsigned char>::const_iterator dim_it =
             elem_dims.begin(); dim_it != elem_dims.end(); ++dim_it)
        {
          const unsigned char dim = *dim_it;

          context.get_element_fe( var, fe, dim );
          if (fe->get_fe_type().family == SCALAR)
            continue;
          if (dim > 1)
            context.get_side_fe( var, side_fe, dim );
          if (dim > 2)
            context.get_edge_fe( var, edge_fe );

          fe->get_xyz();

          fe->get_phi();
          if (dim > 1)
            side_fe->get_phi();
          if (dim > 2)
            edge_fe->get_phi();

          const FEContinuity cont = fe->get_continuity();
          if (cont == C_ONE)
            {
              // Our C1 elements need gradient information
              libmesh_assert(g);

              fe->get_dphi();
              if (dim > 1)
                side_fe->get_dphi();
              if (dim > 2)
                edge_fe->get_dphi();
            }
        }
    }

  // Now initialize any user requested shape functions, xyz vals, etc
  f.init_context(context);
  if (g.get())
    g->init_context(context);

  // this->init_context(context);

  // Iterate over all the elements in the range
  for (ConstElemRange::const_iterator elem_it=range.begin(); elem_it != range.end();
       ++elem_it)
    {
      const Elem * elem = *elem_it;

      unsigned int dim = elem->dim();

      context.pre_fe_reinit(system, elem);

      // If we're doing AMR, this might be a grid projection with a cheap
      // early exit.
#ifdef LIBMESH_ENABLE_AMR
      // If this element doesn't have an old_dof_object, but it
      // wasn't just refined or just coarsened into activity, then
      // it must be newly added, so the user is responsible for
      // setting the new dofs on it during a grid projection.
      if (!elem->old_dof_object &&
          elem->refinement_flag() != Elem::JUST_REFINED &&
          elem->refinement_flag() != Elem::JUST_COARSENED &&
          f.is_grid_projection())
        continue;
#endif // LIBMESH_ENABLE_AMR

      // Loop over all the variables we've been requested to project, to
      // do the projection
      for (std::size_t v=0; v!=variables.size(); v++)
        {
          const unsigned int var = variables[v];

          const Variable & variable = dof_map.variable(var);

          const FEType & base_fe_type = variable.type();

          FEType fe_type = base_fe_type;

          // This may be a p refined element
          fe_type.order =
            libMesh::Order (fe_type.order + elem->p_level());

          if (fe_type.family == SCALAR)
            continue;

          // Per-subdomain variables don't need to be projected on
          // elements where they're not active
          if (!variable.active_on_subdomain(elem->subdomain_id()))
            continue;

          const std::vector<dof_id_type> & dof_indices =
            context.get_dof_indices(var);

          // The number of DOFs on the element
          const unsigned int n_dofs =
            cast_int<unsigned int>(dof_indices.size());

          const unsigned int var_component =
            system.variable_scalar_number(var, 0);

          // Zero the interpolated values
          Ue.resize (n_dofs); Ue.zero();

          // If we're projecting from an old grid
#ifdef LIBMESH_ENABLE_AMR
          if (f.is_grid_projection() &&
              // And either this is an unchanged element
              ((elem->refinement_flag() != Elem::JUST_REFINED &&
                elem->refinement_flag() != Elem::JUST_COARSENED &&
                elem->p_refinement_flag() != Elem::JUST_REFINED &&
                elem->p_refinement_flag() != Elem::JUST_COARSENED) ||
               // Or this is a low order monomial element which has merely
               // been h refined
               (fe_type.family == MONOMIAL &&
                fe_type.order == CONSTANT &&
                elem->p_level() == 0 &&
                elem->refinement_flag() != Elem::JUST_COARSENED &&
                elem->p_refinement_flag() != Elem::JUST_COARSENED))
              )
            // then we can simply copy its old dof
            // values to new indices.
            {
              LOG_SCOPE ("copy_dofs", "GenericProjector");

              f.eval_old_dofs(context, var_component, Ue.get_values());

              action.insert(context, var, Ue);

              continue;
            }
#endif // LIBMESH_ENABLE_AMR

          FEBase * fe = libmesh_nullptr;
          FEBase * side_fe = libmesh_nullptr;
          FEBase * edge_fe = libmesh_nullptr;

          context.get_element_fe( var, fe, dim );
          if (fe->get_fe_type().family == SCALAR)
            continue;
          if (dim > 1)
            context.get_side_fe( var, side_fe, dim );
          if (dim > 2)
            context.get_edge_fe( var, edge_fe );

          const FEContinuity cont = fe->get_continuity();

          std::vector<unsigned int> side_dofs;

          // Fixed vs. free DoFs on edge/face projections
          std::vector<char> dof_is_fixed(n_dofs, false); // bools
          std::vector<int> free_dof(n_dofs, 0);

          // The element type
          const ElemType elem_type = elem->type();

          // The number of nodes on the new element
          const unsigned int n_nodes = elem->n_nodes();

          START_LOG ("project_nodes","GenericProjector");

          // When interpolating C1 elements we need to know which
          // vertices were also parent vertices; we'll cache an
          // intermediate fact outside the nodes loop.
          int i_am_child = -1;
#ifdef LIBMESH_ENABLE_AMR
          if (elem->parent())
            i_am_child = elem->parent()->which_child_am_i(elem);
#endif
          // In general, we need a series of
          // projections to ensure a unique and continuous
          // solution.  We start by interpolating nodes, then
          // hold those fixed and project edges, then
          // hold those fixed and project faces, then
          // hold those fixed and project interiors
          //
          // In the LAGRANGE case, we will save a lot of solution
          // evaluations (at a slight cost in accuracy) by simply
          // interpolating all nodes rather than projecting.

          // Interpolate vertex (or for LAGRANGE, all node) values first.
          unsigned int current_dof = 0;
          for (unsigned int n=0; n!= n_nodes; ++n)
            {
              // FIXME: this should go through the DofMap,
              // not duplicate dof_indices code badly!
              const unsigned int nc =
                FEInterface::n_dofs_at_node (dim, fe_type, elem_type, n);

              if (!elem->is_vertex(n) &&
                  fe_type.family != LAGRANGE)
                {
                  current_dof += nc;
                  continue;
                }

              if (cont == DISCONTINUOUS)
                {
                  libmesh_assert_equal_to (nc, 0);
                }
              else if (!nc)
                {
                  // This should only occur for first-order LAGRANGE
                  // FE on non-vertices of higher-order elements
                  libmesh_assert (!elem->is_vertex(n));
                  libmesh_assert_equal_to(fe_type.family, LAGRANGE);
                }
              // Assume that C_ZERO elements have a single nodal
              // value shape function at vertices
              else if (cont == C_ZERO)
                {
                  Ue(current_dof) = f.eval_at_node(context,
                                                   var_component,
                                                   dim,
                                                   elem->node_ref(n),
                                                   system.time);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                }
              else if (cont == C_ONE)
                {
                  const bool is_parent_vertex = (i_am_child == -1) ||
                    elem->parent()->is_vertex_on_parent(i_am_child, n);

                  // The hermite element vertex shape functions are weird
                  if (fe_type.family == HERMITE)
                    {
                      Ue(current_dof) =
                        f.eval_at_node(context,
                                       var_component,
                                       dim,
                                       elem->node_ref(n),
                                       system.time);
                      dof_is_fixed[current_dof] = true;
                      current_dof++;
                      VectorValue<FValue> grad =
                        is_parent_vertex ?
                        g->eval_at_node(context,
                                        var_component,
                                        dim,
                                        elem->node_ref(n),
                                        system.time) :
                        g->eval_at_point(context,
                                         var_component,
                                         elem->node_ref(n),
                                         system.time);
                      // x derivative
                      Ue(current_dof) = grad(0);
                      dof_is_fixed[current_dof] = true;
                      current_dof++;
                      if (dim > 1)
                        {
                          // We'll finite difference mixed derivatives
                          Real delta_x = TOLERANCE * elem->hmin();

                          Point nxminus = elem->point(n),
                            nxplus = elem->point(n);
                          nxminus(0) -= delta_x;
                          nxplus(0) += delta_x;
                          VectorValue<FValue> gxminus =
                            g->eval_at_point(context,
                                             var_component,
                                             nxminus,
                                             system.time);
                          VectorValue<FValue> gxplus =
                            g->eval_at_point(context,
                                             var_component,
                                             nxplus,
                                             system.time);
                          // y derivative
                          Ue(current_dof) = grad(1);
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                          // xy derivative
                          Ue(current_dof) = (gxplus(1) - gxminus(1))
                            / 2. / delta_x;
                          dof_is_fixed[current_dof] = true;
                          current_dof++;

                          if (dim > 2)
                            {
                              // z derivative
                              Ue(current_dof) = grad(2);
                              dof_is_fixed[current_dof] = true;
                              current_dof++;
                              // xz derivative
                              Ue(current_dof) = (gxplus(2) - gxminus(2))
                                / 2. / delta_x;
                              dof_is_fixed[current_dof] = true;
                              current_dof++;
                              // We need new points for yz
                              Point nyminus = elem->point(n),
                                nyplus = elem->point(n);
                              nyminus(1) -= delta_x;
                              nyplus(1) += delta_x;
                              VectorValue<FValue> gyminus =
                                g->eval_at_point(context,
                                                 var_component,
                                                 nyminus,
                                                 system.time);
                              VectorValue<FValue> gyplus =
                                g->eval_at_point(context,
                                                 var_component,
                                                 nyplus,
                                                 system.time);
                              // xz derivative
                              Ue(current_dof) = (gyplus(2) - gyminus(2))
                                / 2. / delta_x;
                              dof_is_fixed[current_dof] = true;
                              current_dof++;
                              // Getting a 2nd order xyz is more tedious
                              Point nxmym = elem->point(n),
                                nxmyp = elem->point(n),
                                nxpym = elem->point(n),
                                nxpyp = elem->point(n);
                              nxmym(0) -= delta_x;
                              nxmym(1) -= delta_x;
                              nxmyp(0) -= delta_x;
                              nxmyp(1) += delta_x;
                              nxpym(0) += delta_x;
                              nxpym(1) -= delta_x;
                              nxpyp(0) += delta_x;
                              nxpyp(1) += delta_x;
                              VectorValue<FValue> gxmym =
                                g->eval_at_point(context,
                                                 var_component,
                                                 nxmym,
                                                 system.time);
                              VectorValue<FValue> gxmyp =
                                g->eval_at_point(context,
                                                 var_component,
                                                 nxmyp,
                                                 system.time);
                              VectorValue<FValue> gxpym =
                                g->eval_at_point(context,
                                                 var_component,
                                                 nxpym,
                                                 system.time);
                              VectorValue<FValue> gxpyp =
                                g->eval_at_point(context,
                                                 var_component,
                                                 nxpyp,
                                                 system.time);
                              FValue gxzplus = (gxpyp(2) - gxmyp(2))
                                / 2. / delta_x;
                              FValue gxzminus = (gxpym(2) - gxmym(2))
                                / 2. / delta_x;
                              // xyz derivative
                              Ue(current_dof) = (gxzplus - gxzminus)
                                / 2. / delta_x;
                              dof_is_fixed[current_dof] = true;
                              current_dof++;
                            }
                        }
                    }
                  else
                    {
                      // Assume that other C_ONE elements have a single nodal
                      // value shape function and nodal gradient component
                      // shape functions
                      libmesh_assert_equal_to (nc, 1 + dim);
                      Ue(current_dof) = f.eval_at_node(context,
                                                       var_component,
                                                       dim,
                                                       elem->node_ref(n),
                                                       system.time);
                      dof_is_fixed[current_dof] = true;
                      current_dof++;
                      VectorValue<FValue> grad =
                        is_parent_vertex ?
                        g->eval_at_node(context,
                                        var_component,
                                        dim,
                                        elem->node_ref(n),
                                        system.time) :
                        g->eval_at_point(context,
                                         var_component,
                                         elem->node_ref(n),
                                         system.time);
                      for (unsigned int i=0; i!= dim; ++i)
                        {
                          Ue(current_dof) = grad(i);
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                        }
                    }
                }
              else
                libmesh_error_msg("Unknown continuity " << cont);
            }

          STOP_LOG ("project_nodes","GenericProjector");

          START_LOG ("project_edges","GenericProjector");

          // In 3D with non-LAGRANGE, project any edge values next
          if (dim > 2 &&
              cont != DISCONTINUOUS &&
              (fe_type.family != LAGRANGE
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
               || elem->infinite()
#endif
               ))
            {
              // If we're JUST_COARSENED we'll need a custom
              // evaluation, not just the standard edge FE
              const std::vector<Point> & xyz_values =
#ifdef LIBMESH_ENABLE_AMR
                (elem->refinement_flag() == Elem::JUST_COARSENED) ?
                fe->get_xyz() :
#endif
                edge_fe->get_xyz();
              const std::vector<Real> & JxW =
#ifdef LIBMESH_ENABLE_AMR
                (elem->refinement_flag() == Elem::JUST_COARSENED) ?
                fe->get_JxW() :
#endif
                edge_fe->get_JxW();

              const std::vector<std::vector<Real>> & phi =
#ifdef LIBMESH_ENABLE_AMR
                (elem->refinement_flag() == Elem::JUST_COARSENED) ?
                fe->get_phi() :
#endif
                edge_fe->get_phi();
              const std::vector<std::vector<RealGradient>> * dphi = libmesh_nullptr;
              if (cont == C_ONE)
                dphi =
#ifdef LIBMESH_ENABLE_AMR
                  (elem->refinement_flag() == Elem::JUST_COARSENED) ?
                  &(fe->get_dphi()) :
#endif
                  &(edge_fe->get_dphi());

              for (auto e : elem->edge_index_range())
                {
                  context.edge = e;

#ifdef LIBMESH_ENABLE_AMR
                  if (elem->refinement_flag() == Elem::JUST_COARSENED)
                    {
                      std::vector<Point> fine_points;

                      std::unique_ptr<FEBase> fine_fe
                        (FEBase::build (dim, base_fe_type));

                      std::unique_ptr<QBase> qrule
                        (base_fe_type.default_quadrature_rule(1));
                      fine_fe->attach_quadrature_rule(qrule.get());

                      const std::vector<Point> & child_xyz =
                        fine_fe->get_xyz();

                      for (unsigned int c = 0, nc = elem->n_children();
                           c != nc; ++c)
                        {
                          if (!elem->is_child_on_edge(c, e))
                            continue;

                          fine_fe->edge_reinit(elem->child_ptr(c), e);
                          fine_points.insert(fine_points.end(),
                                             child_xyz.begin(),
                                             child_xyz.end());
                        }

                      std::vector<Point> fine_qp;
                      FEInterface::inverse_map (dim, base_fe_type, elem,
                                                fine_points, fine_qp);

                      context.elem_fe_reinit(&fine_qp);
                    }
                  else
#endif // LIBMESH_ENABLE_AMR
                    context.edge_fe_reinit();

                  const unsigned int n_qp = xyz_values.size();

                  FEInterface::dofs_on_edge(elem, dim, base_fe_type,
                                            e, side_dofs);

                  // Some edge dofs are on nodes and already
                  // fixed, others are free to calculate
                  unsigned int free_dofs = 0;
                  for (std::size_t i=0; i != side_dofs.size(); ++i)
                    if (!dof_is_fixed[side_dofs[i]])
                      free_dof[free_dofs++] = i;

                  // There may be nothing to project
                  if (!free_dofs)
                    continue;

                  Ke.resize (free_dofs, free_dofs); Ke.zero();
                  Fe.resize (free_dofs); Fe.zero();
                  // The new edge coefficients
                  DenseVector<FValue> Uedge(free_dofs);

                  // Loop over the quadrature points
                  for (unsigned int qp=0; qp<n_qp; qp++)
                    {
                      // solution at the quadrature point
                      FValue fineval = f.eval_at_point(context,
                                                       var_component,
                                                       xyz_values[qp],
                                                       system.time);
                      // solution grad at the quadrature point
                      VectorValue<FValue> finegrad;
                      if (cont == C_ONE)
                        finegrad = g->eval_at_point(context,
                                                    var_component,
                                                    xyz_values[qp],
                                                    system.time);

                      // Form edge projection matrix
                      for (std::size_t sidei=0, freei=0;
                           sidei != side_dofs.size(); ++sidei)
                        {
                          unsigned int i = side_dofs[sidei];
                          // fixed DoFs aren't test functions
                          if (dof_is_fixed[i])
                            continue;
                          for (std::size_t sidej=0, freej=0;
                               sidej != side_dofs.size(); ++sidej)
                            {
                              unsigned int j = side_dofs[sidej];
                              if (dof_is_fixed[j])
                                Fe(freei) -= phi[i][qp] * phi[j][qp] *
                                  JxW[qp] * Ue(j);
                              else
                                Ke(freei,freej) += phi[i][qp] *
                                  phi[j][qp] * JxW[qp];
                              if (cont == C_ONE)
                                {
                                  if (dof_is_fixed[j])
                                    Fe(freei) -= ( (*dphi)[i][qp] *
                                                   (*dphi)[j][qp] ) *
                                      JxW[qp] * Ue(j);
                                  else
                                    Ke(freei,freej) += ( (*dphi)[i][qp] *
                                                         (*dphi)[j][qp] )
                                      * JxW[qp];
                                }
                              if (!dof_is_fixed[j])
                                freej++;
                            }
                          Fe(freei) += phi[i][qp] * fineval * JxW[qp];
                          if (cont == C_ONE)
                            Fe(freei) += (finegrad * (*dphi)[i][qp] ) *
                              JxW[qp];
                          freei++;
                        }
                    }

                  Ke.cholesky_solve(Fe, Uedge);

                  // Transfer new edge solutions to element
                  for (unsigned int i=0; i != free_dofs; ++i)
                    {
                      FValue & ui = Ue(side_dofs[free_dof[i]]);
                      libmesh_assert(bool(std::abs(ui) < TOLERANCE) ||
                                     bool(std::abs(ui - Uedge(i)) < TOLERANCE));
                      ui = Uedge(i);
                      dof_is_fixed[side_dofs[free_dof[i]]] = true;
                    }
                }
            } // end if (dim > 2, !DISCONTINUOUS, !LAGRANGE)

          STOP_LOG ("project_edges","GenericProjector");

          START_LOG ("project_sides","GenericProjector");

          // With non-LAGRANGE, project any side values (edges in 2D,
          // faces in 3D) next.
          if (dim > 1 &&
              cont != DISCONTINUOUS &&
              (fe_type.family != LAGRANGE
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
               || elem->infinite()
#endif
               ))
            {
              // If we're JUST_COARSENED we'll need a custom
              // evaluation, not just the standard side FE
              const std::vector<Point> & xyz_values =
#ifdef LIBMESH_ENABLE_AMR
                (elem->refinement_flag() == Elem::JUST_COARSENED) ?
                fe->get_xyz() :
#endif // LIBMESH_ENABLE_AMR
                side_fe->get_xyz();
              const std::vector<Real> & JxW =
#ifdef LIBMESH_ENABLE_AMR
                (elem->refinement_flag() == Elem::JUST_COARSENED) ?
                fe->get_JxW() :
#endif // LIBMESH_ENABLE_AMR
                side_fe->get_JxW();
              const std::vector<std::vector<Real>> & phi =
#ifdef LIBMESH_ENABLE_AMR
                (elem->refinement_flag() == Elem::JUST_COARSENED) ?
                fe->get_phi() :
#endif // LIBMESH_ENABLE_AMR
                side_fe->get_phi();
              const std::vector<std::vector<RealGradient>> * dphi = libmesh_nullptr;
              if (cont == C_ONE)
                dphi =
#ifdef LIBMESH_ENABLE_AMR
                  (elem->refinement_flag() == Elem::JUST_COARSENED) ? &(fe->get_dphi()) :
#endif // LIBMESH_ENABLE_AMR
                  &(side_fe->get_dphi());

              for (auto s : elem->side_index_range())
                {
                  FEInterface::dofs_on_side(elem, dim, base_fe_type,
                                            s, side_dofs);

                  // Some side dofs are on nodes/edges and already
                  // fixed, others are free to calculate
                  unsigned int free_dofs = 0;
                  for (std::size_t i=0; i != side_dofs.size(); ++i)
                    if (!dof_is_fixed[side_dofs[i]])
                      free_dof[free_dofs++] = i;

                  // There may be nothing to project
                  if (!free_dofs)
                    continue;

                  Ke.resize (free_dofs, free_dofs); Ke.zero();
                  Fe.resize (free_dofs); Fe.zero();
                  // The new side coefficients
                  DenseVector<FValue> Uside(free_dofs);

                  context.side = s;

#ifdef LIBMESH_ENABLE_AMR
                  if (elem->refinement_flag() == Elem::JUST_COARSENED)
                    {
                      std::vector<Point> fine_points;

                      std::unique_ptr<FEBase> fine_fe
                        (FEBase::build (dim, base_fe_type));

                      std::unique_ptr<QBase> qrule
                        (base_fe_type.default_quadrature_rule(dim-1));
                      fine_fe->attach_quadrature_rule(qrule.get());

                      const std::vector<Point> & child_xyz =
                        fine_fe->get_xyz();

                      for (unsigned int c = 0, nc = elem->n_children();
                           c != nc; ++c)
                        {
                          if (!elem->is_child_on_side(c, s))
                            continue;

                          fine_fe->reinit(elem->child_ptr(c), s);
                          fine_points.insert(fine_points.end(),
                                             child_xyz.begin(),
                                             child_xyz.end());
                        }

                      std::vector<Point> fine_qp;
                      FEInterface::inverse_map (dim, base_fe_type, elem,
                                                fine_points, fine_qp);

                      context.elem_fe_reinit(&fine_qp);
                    }
                  else
#endif // LIBMESH_ENABLE_AMR
                    context.side_fe_reinit();

                  const unsigned int n_qp = xyz_values.size();

                  // Loop over the quadrature points
                  for (unsigned int qp=0; qp<n_qp; qp++)
                    {
                      // solution at the quadrature point
                      FValue fineval = f.eval_at_point(context,
                                                       var_component,
                                                       xyz_values[qp],
                                                       system.time);
                      // solution grad at the quadrature point
                      VectorValue<FValue> finegrad;
                      if (cont == C_ONE)
                        finegrad = g->eval_at_point(context,
                                                    var_component,
                                                    xyz_values[qp],
                                                    system.time);

                      // Form side projection matrix
                      for (std::size_t sidei=0, freei=0;
                           sidei != side_dofs.size(); ++sidei)
                        {
                          unsigned int i = side_dofs[sidei];
                          // fixed DoFs aren't test functions
                          if (dof_is_fixed[i])
                            continue;
                          for (std::size_t sidej=0, freej=0;
                               sidej != side_dofs.size(); ++sidej)
                            {
                              unsigned int j = side_dofs[sidej];
                              if (dof_is_fixed[j])
                                Fe(freei) -= phi[i][qp] * phi[j][qp] *
                                  JxW[qp] * Ue(j);
                              else
                                Ke(freei,freej) += phi[i][qp] *
                                  phi[j][qp] * JxW[qp];
                              if (cont == C_ONE)
                                {
                                  if (dof_is_fixed[j])
                                    Fe(freei) -= ( (*dphi)[i][qp] *
                                                   (*dphi)[j][qp] ) *
                                      JxW[qp] * Ue(j);
                                  else
                                    Ke(freei,freej) += ( (*dphi)[i][qp] *
                                                         (*dphi)[j][qp] )
                                      * JxW[qp];
                                }
                              if (!dof_is_fixed[j])
                                freej++;
                            }
                          Fe(freei) += (fineval * phi[i][qp]) * JxW[qp];
                          if (cont == C_ONE)
                            Fe(freei) += (finegrad * (*dphi)[i][qp] ) *
                              JxW[qp];
                          freei++;
                        }
                    }

                  Ke.cholesky_solve(Fe, Uside);

                  // Transfer new side solutions to element
                  for (unsigned int i=0; i != free_dofs; ++i)
                    {
                      FValue & ui = Ue(side_dofs[free_dof[i]]);
                      libmesh_assert(bool(std::abs(ui) < TOLERANCE) ||
                                     bool(std::abs(ui - Uside(i)) < TOLERANCE));
                      ui = Uside(i);
                      dof_is_fixed[side_dofs[free_dof[i]]] = true;
                    }
                }
            } // end if (dim > 1, !DISCONTINUOUS, !LAGRANGE)

          STOP_LOG ("project_sides","GenericProjector");

          START_LOG ("project_interior","GenericProjector");

          // Project the interior values, finally

          // Some interior dofs are on nodes/edges/sides and
          // already fixed, others are free to calculate
          unsigned int free_dofs = 0;
          for (unsigned int i=0; i != n_dofs; ++i)
            if (!dof_is_fixed[i])
              free_dof[free_dofs++] = i;

          // Project any remaining (interior) dofs in the non-LAGRANGE
          // case.
          if (free_dofs &&
              (fe_type.family != LAGRANGE
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
               || elem->infinite()
#endif
               ))
            {
              const std::vector<Point> & xyz_values = fe->get_xyz();
              const std::vector<Real> & JxW = fe->get_JxW();

              const std::vector<std::vector<Real>> & phi = fe->get_phi();
              const std::vector<std::vector<RealGradient>> * dphi = libmesh_nullptr;
              if (cont == C_ONE)
                dphi = &(fe->get_dphi());

#ifdef LIBMESH_ENABLE_AMR
              if (elem->refinement_flag() == Elem::JUST_COARSENED)
                {
                  std::vector<Point> fine_points;

                  std::unique_ptr<FEBase> fine_fe
                    (FEBase::build (dim, base_fe_type));

                  std::unique_ptr<QBase> qrule
                    (base_fe_type.default_quadrature_rule(dim));
                  fine_fe->attach_quadrature_rule(qrule.get());

                  const std::vector<Point> & child_xyz =
                    fine_fe->get_xyz();

                  for (auto & child : elem->child_ref_range())
                    {
                      fine_fe->reinit(&child);
                      fine_points.insert(fine_points.end(),
                                         child_xyz.begin(),
                                         child_xyz.end());
                    }

                  std::vector<Point> fine_qp;
                  FEInterface::inverse_map (dim, base_fe_type, elem,
                                            fine_points, fine_qp);

                  context.elem_fe_reinit(&fine_qp);
                }
              else
#endif // LIBMESH_ENABLE_AMR
                context.elem_fe_reinit();

              const unsigned int n_qp = xyz_values.size();

              Ke.resize (free_dofs, free_dofs); Ke.zero();
              Fe.resize (free_dofs); Fe.zero();
              // The new interior coefficients
              DenseVector<FValue> Uint(free_dofs);

              // Loop over the quadrature points
              for (unsigned int qp=0; qp<n_qp; qp++)
                {
                  // solution at the quadrature point
                  FValue fineval = f.eval_at_point(context,
                                                   var_component,
                                                   xyz_values[qp],
                                                   system.time);
                  // solution grad at the quadrature point
                  VectorValue<FValue> finegrad;
                  if (cont == C_ONE)
                    finegrad = g->eval_at_point(context,
                                                var_component,
                                                xyz_values[qp],
                                                system.time);

                  // Form interior projection matrix
                  for (unsigned int i=0, freei=0; i != n_dofs; ++i)
                    {
                      // fixed DoFs aren't test functions
                      if (dof_is_fixed[i])
                        continue;
                      for (unsigned int j=0, freej=0; j != n_dofs; ++j)
                        {
                          if (dof_is_fixed[j])
                            Fe(freei) -= phi[i][qp] * phi[j][qp] * JxW[qp]
                              * Ue(j);
                          else
                            Ke(freei,freej) += phi[i][qp] * phi[j][qp] *
                              JxW[qp];
                          if (cont == C_ONE)
                            {
                              if (dof_is_fixed[j])
                                Fe(freei) -= ( (*dphi)[i][qp] *
                                               (*dphi)[j][qp] ) * JxW[qp] *
                                  Ue(j);
                              else
                                Ke(freei,freej) += ( (*dphi)[i][qp] *
                                                     (*dphi)[j][qp] ) *
                                  JxW[qp];
                            }
                          if (!dof_is_fixed[j])
                            freej++;
                        }
                      Fe(freei) += phi[i][qp] * fineval * JxW[qp];
                      if (cont == C_ONE)
                        Fe(freei) += (finegrad * (*dphi)[i][qp] ) * JxW[qp];
                      freei++;
                    }
                }
              Ke.cholesky_solve(Fe, Uint);

              // Transfer new interior solutions to element
              for (unsigned int i=0; i != free_dofs; ++i)
                {
                  FValue & ui = Ue(free_dof[i]);
                  libmesh_assert(bool(std::abs(ui) < TOLERANCE) ||
                                 bool(std::abs(ui - Uint(i)) < TOLERANCE));
                  ui = Uint(i);
                  dof_is_fixed[free_dof[i]] = true;
                }

            } // if there are free interior dofs

          STOP_LOG ("project_interior","GenericProjector");

          // Make sure every DoF got reached!
          for (unsigned int i=0; i != n_dofs; ++i)
            libmesh_assert(dof_is_fixed[i]);

          action.insert(context, var, Ue);
        } // end variables loop
    } // end elements loop
}

} // namespace libMesh

#endif // GENERIC_PROJECTOR_H
