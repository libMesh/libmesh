// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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
#include "libmesh/mesh_tools.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parallel_sync.h"
#include "libmesh/quadrature.h"
#include "libmesh/system.h"
#include "libmesh/threads.h"

namespace libMesh
{

/**
 * For ease of communication, we allow users to translate their own
 * value types to a more easily computable (typically a vector of some
 * fixed-size type) output, by specializing these calls using
 * different types.
 */
template <typename T>
struct TypeToSend {
  typedef T type;
};

template <typename T>
const typename TypeToSend<T>::type convert_to_send (const T& in)
{ return in; }

template <typename SendT, typename T>
void convert_from_receive (SendT & received, T & converted)
{ converted = received; }

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
  bool g_was_copied, map_was_created;
  const ProjectionAction & master_action;
  const std::vector<unsigned int> & variables;
  std::unordered_map<dof_id_type, std::vector<dof_id_type>> * nodes_to_elem;

public:
  GenericProjector (const System & system_in,
                    const FFunctor & f_in,
                    const GFunctor * g_in,
                    const ProjectionAction & act_in,
                    const std::vector<unsigned int> & variables_in,
                    std::unordered_map<dof_id_type, std::vector<dof_id_type>> *
                      nodes_to_elem_in = nullptr) :
    system(system_in),
    master_f(f_in),
    master_g(g_in),
    g_was_copied(false),
    map_was_created(!nodes_to_elem_in),
    master_action(act_in),
    variables(variables_in),
    nodes_to_elem(nodes_to_elem_in)
  {
    if (map_was_created) // past tense misnomer here
      {
        nodes_to_elem = new
          std::unordered_map<dof_id_type, std::vector<dof_id_type>>;
        MeshTools::build_nodes_to_elem_map (system.get_mesh(), *nodes_to_elem);
      }
  }

  GenericProjector (const GenericProjector & in) :
    system(in.system),
    master_f(in.master_f),
    master_g(in.master_g ? new GFunctor(*in.master_g) : nullptr),
    g_was_copied(in.master_g),
    master_action(in.master_action),
    variables(in.variables),
    nodes_to_elem(in.nodes_to_elem)
  {}

  ~GenericProjector()
  {
    if (g_was_copied)
      delete master_g;
    if (map_was_created)
      delete nodes_to_elem;
  }

  void operator() (const ConstElemRange & range) const;

  void find_dofs_to_send
    (const Node & node,
     std::unordered_map<dof_id_type, std::pair<FValue, processor_id_type>> & ids_to_push) const;

  void send_and_insert_dof_values
    (std::unordered_map<dof_id_type, std::pair<FValue, processor_id_type>> & ids_to_push,
     ProjectionAction & action) const;

  template <typename InsertId>
  void construct_projection
    (const std::vector<dof_id_type> & dof_indices_var,
     const std::vector<unsigned int> & involved_dofs,
     unsigned int var_component,
     const FEMContext & context,
     const Node * node,
     InsertId & insert_id,
     const FEBase & fe,
     FFunctor & f,
     GFunctor * g,
     const std::unordered_map<dof_id_type, FValue> & ids_to_save
     ) const;
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

  void insert(dof_id_type id,
              Val val)
  {
    // Lock the new vector since it is shared among threads.
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
      target_vector.set(id, val);
    }
  }


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
    last_elem(nullptr),
    sys(sys_in),
    old_context(sys_in)
  {
    // We'll be queried for components but we'll typically be looking
    // up data by variables, and those indices don't always match
    for (auto v : IntRange<unsigned int>(0, sys.n_vars()))
      {
        const unsigned int vcomp = sys.variable_scalar_number(v,0);
        if (vcomp >= component_to_var.size())
          component_to_var.resize(vcomp+1, static_cast<unsigned int>(-1));
        component_to_var[vcomp] = v;
      }
  }

  OldSolutionBase(const OldSolutionBase & in) :
    last_elem(nullptr),
    sys(in.sys),
    old_context(sys),
    component_to_var(in.component_to_var)
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
        FEBase * fe = nullptr;
        const std::set<unsigned char> & elem_dims =
          old_context.elem_dimensions();

        for (const auto & dim : elem_dims)
          {
            old_context.get_element_fe(var, fe, dim);
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
  std::vector<unsigned int> component_to_var;

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

    // Handle offset from non-scalar components in previous variables
    libmesh_assert_less(i, this->component_to_var.size());
    unsigned int var = this->component_to_var[i];

    Output n;
    (this->old_context.*point_output)(var, p, n, this->out_of_elem_tol);
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

  // Handle offset from non-scalar components in previous variables
  libmesh_assert_less(i, this->component_to_var.size());
  unsigned int var = this->component_to_var[i];

  // Optimize for the common case, where this node was part of the
  // old solution.
  //
  // Be sure to handle cases where the variable wasn't defined on
  // this node (due to changing subdomain support) or where the
  // variable has no components on this node (due to Elem order
  // exceeding FE order)
  if (n.old_dof_object &&
      n.old_dof_object->n_vars(sys.number()) &&
      n.old_dof_object->n_comp(sys.number(), var))
    {
      const dof_id_type old_id =
        n.old_dof_object->dof_number(sys.number(), var, 0);
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

  // Handle offset from non-scalar components in previous variables
  libmesh_assert_less(i, this->component_to_var.size());
  unsigned int var = this->component_to_var[i];

  // Optimize for the common case, where this node was part of the
  // old solution.
  //
  // Be sure to handle cases where the variable wasn't defined on
  // this node (due to changing subdomain support) or where the
  // variable has no components on this node (due to Elem order
  // exceeding FE order)
  if (n.old_dof_object &&
      n.old_dof_object->n_vars(sys.number()) &&
      n.old_dof_object->n_comp(sys.number(), var))
    {
      Gradient g;
      for (unsigned int d = 0; d != elem_dim; ++d)
        {
          const dof_id_type old_id =
            n.old_dof_object->dof_number(sys.number(), var, d+1);
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

  // Context objects to contain all our required FE objects
  FEMContext context( system );

  std::vector<FEContinuity> conts(system.n_vars());

  // Loop over all the variables we've been requested to project, to
  // pre-request
  for (const auto & var : variables)
    {
      // FIXME: Need to generalize this to vector-valued elements. [PB]
      FEBase * fe = nullptr;
      FEBase * side_fe = nullptr;
      FEBase * edge_fe = nullptr;

      const std::set<unsigned char> & elem_dims =
        context.elem_dimensions();

      for (const auto & dim : elem_dims)
        {
          context.get_element_fe( var, fe, dim );
          if (fe->get_fe_type().family == SCALAR)
            continue;
          if (dim > 1)
            context.get_side_fe( var, side_fe, dim );
          if (dim > 2)
            context.get_edge_fe( var, edge_fe );

          fe->get_JxW();
          fe->get_xyz();
          fe->get_JxW();

          fe->get_phi();
          if (dim > 1)
            {
              side_fe->get_JxW();
              side_fe->get_xyz();
              side_fe->get_phi();
            }
          if (dim > 2)
            {
              edge_fe->get_JxW();
              edge_fe->get_xyz();
              edge_fe->get_phi();
            }

          const FEContinuity cont = fe->get_continuity();
          conts[var] = cont;
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

  // Look at all the elements in the range.  Determine sets of
  // vertices, edge nodes, and side nodes to project.
  //
  // As per our other weird nomenclature, "sides" means faces in 3D
  // and edges in 2D, and "edges" gets skipped in 2D
  //
  // This gets tricky in the case of subdomain-restricted
  // variables, for which we might need to do the same projection
  // from different elements when evaluating different variables.
  // We'll keep track of which variables can be projected from which
  // elements.

  typedef std::unordered_set<unsigned int> var_set;

  std::unordered_multimap<const Node *, std::pair<const Elem *, var_set>> vertices;
  std::unordered_multimap<const Node *, std::tuple<const Elem *, unsigned short, var_set>> edges;
  std::unordered_multimap<const Node *, std::tuple<const Elem *, unsigned short, var_set>> sides;
  std::vector<const Elem *> interiors;

  for (const auto & elem : range)
    {
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

      const int dim = elem->dim();

      const unsigned int n_vertices = elem->n_vertices();
      const unsigned int n_edges = elem->n_edges();
      const unsigned int n_nodes = elem->n_nodes();

      // In 1-D we already handle our sides as vertices
      const unsigned int n_sides = (dim > 1) * elem->n_sides();

      // What variables are supported on each kind of node on elem?
      var_set vertex_vars, edge_vars, side_vars;

      // If we have non-vertex nodes, the first is an edge node, but
      // if we're in 2D we'll call that a side node
      const bool has_edge_nodes = (n_nodes > n_vertices && dim > 2);

      // If we have even more nodes, the next is a side node.
      const bool has_side_nodes =
        (n_nodes > n_vertices + ((dim > 2) * n_edges));

      // We may be out of nodes at this point or we have interior
      // nodes which may have DoFs to project too
      const bool has_interior_nodes =
        (n_nodes > n_vertices + ((dim > 2) * n_edges) + n_sides);

      for (auto v_num : variables)
        {
          const Variable & var = system.variable(v_num);
          if (!var.active_on_subdomain(elem->subdomain_id()))
            continue;
          FEType fe_type = var.type();
          fe_type.order =
            libMesh::Order (fe_type.order + elem->p_level());
          const ElemType elem_type = elem->type();

          if (FEInterface::n_dofs_at_node(dim, fe_type, elem_type, 0))
            vertex_vars.insert(vertex_vars.end(), v_num);

          // The first non-vertex node is always an edge node if those
          // exist.  All edge nodes have the same number of DoFs
          if (has_edge_nodes)
            if (FEInterface::n_dofs_at_node(dim, fe_type, elem_type, n_vertices))
              edge_vars.insert(edge_vars.end(), v_num);

          if (has_side_nodes)
            {
              if (dim != 3)
                {
                  if (FEInterface::n_dofs_at_node(dim, fe_type, elem_type, n_vertices))
                    side_vars.insert(side_vars.end(), v_num);
                }
              else
                // In 3D, not all face nodes always have the same number of
                // DoFs!  We'll loop over all sides to be safe.
                for (unsigned int n = 0; n != n_nodes; ++n)
                  if (elem->is_face(n))
                    if (FEInterface::n_dofs_at_node(dim, fe_type,
                                                    elem_type, n))
                      {
                        side_vars.insert(side_vars.end(), v_num);
                        break;
                      }
            }

          if (FEInterface::n_dofs_per_elem(dim, fe_type, elem_type) ||
              (has_interior_nodes &&
               FEInterface::n_dofs_at_node(dim, fe_type, elem_type, n_nodes-1)))
            if (interiors.empty() || interiors.back() != elem)
              interiors.push_back(elem);
        }

      // We'll use a greedy algorithm in most cases: if another
      // element has already claimed some of our DoFs, we'll let it do
      // the work.
      for (unsigned int v=0; v != n_vertices; ++v)
        {
          const Node * node = elem->node_ptr(v);

          auto remaining_vars = vertex_vars;

          auto set_range = vertices.equal_range(node);
          for (const auto & v_ent : as_range(set_range))
            for (const unsigned int var_covered :
                 v_ent.second.second)
              remaining_vars.erase(var_covered);

          if (!remaining_vars.empty())
            vertices.emplace
              (node, std::make_pair(elem, std::move(remaining_vars)));
        }

      if (has_edge_nodes)
        {
          for (unsigned int e=0; e != n_edges; ++e)
            {
              const Node * node = elem->node_ptr(n_vertices+e);

              auto remaining_vars = edge_vars;

              auto set_range = edges.equal_range(node);
              for (const auto & v_ent : as_range(set_range))
                for (const unsigned int var_covered :
                     std::get<2>(v_ent.second))
                  remaining_vars.erase(var_covered);

              if (!remaining_vars.empty())
                edges.emplace
                  (node, std::make_tuple(elem, e, std::move(remaining_vars)));
            }
        }

      if (has_side_nodes)
        {
          for (unsigned int side=0; side != n_sides; ++side)
            {
              const Node * node = nullptr;
              if (dim != 3)
                node = elem->node_ptr(n_vertices+(dim>2)*n_edges+side);
              else
                {
                  // In 3D only some sides may have nodes
                  for (unsigned int n = 0; n != n_nodes; ++n)
                    {
                      if (!elem->is_face(n))
                        continue;

                      if (elem->is_node_on_side(n, side))
                        {
                          node = elem->node_ptr(n);
                          break;
                        }
                    }
                }

              if (!node)
                continue;

              auto remaining_vars = side_vars;

              auto set_range = sides.equal_range(node);
              for (const auto & v_ent : as_range(set_range))
                for (const unsigned int var_covered :
                     std::get<2>(v_ent.second))
                  remaining_vars.erase(var_covered);

              if (!remaining_vars.empty())
                sides.emplace
                  (node, std::make_tuple(elem, side, std::move(remaining_vars)));
            }
        }
    }

  // While we're looping over nodes, also figure out which ghosted
  // nodes will have data we might need to send to their owners
  // instead of being acted on by ourselves.
  //
  // We keep track of which dof ids we might need to send, and what
  // values those ids should get (along with a pprocessor_id to leave
  // invalid in case *we* can't compute those values either).
  std::unordered_map<dof_id_type, std::pair<FValue, processor_id_type>> ids_to_push;

  // We generally need to hang on to every value we've calculated
  // until we're all done, because later projection calculations
  // depend on boundary data from earlier calculations.
  std::unordered_map<dof_id_type, FValue> ids_to_save;
  bool done_saving_ids = edges.empty() && sides.empty() && interiors.empty();

  // When we have new data to act on, we may also need to save it
  // and get ready to push it.
  auto insert_id = [&ids_to_push, &ids_to_save, &action, &done_saving_ids]
    (dof_id_type id, const FValue & val, processor_id_type pid)
    {
      auto iter = ids_to_push.find(id);
      if (iter == ids_to_push.end())
        action.insert(id, val);
      else
        {
          libmesh_assert(pid != DofObject::invalid_processor_id);
          iter->second = std::make_pair(val, pid);
        }
      if (!done_saving_ids)
        {
          libmesh_assert(!ids_to_save.count(id));
          ids_to_save[id] = val;
        }
    };

  START_LOG ("project_vertices","GenericProjector");
  for (const auto & v_pair : vertices)
    {
      const Node & vertex = *v_pair.first;
      const Elem & elem = *v_pair.second.first;
      context.pre_fe_reinit(system, &elem);

      this->find_dofs_to_send(vertex, ids_to_push);

      // Look at all the variables we're supposed to interpolate from
      // this element on this vertex
      for (const auto & var : v_pair.second.second)
        {
          const Variable & variable = system.variable(var);
          const FEType & base_fe_type = variable.type();
          const unsigned int var_component =
            system.variable_scalar_number(var, 0);

          if (base_fe_type.family == SCALAR)
            continue;

          const FEContinuity & cont = conts[var];
          if (cont == DISCONTINUOUS)
            {
              libmesh_assert_equal_to(vertex.n_comp(system.number(), var), 0);
            }
          else if (cont == C_ZERO)
            {
              libmesh_assert(vertex.n_comp(system.number(), var));
              const dof_id_type id = vertex.dof_number(system.number(), var, 0);
              // C_ZERO elements have a single nodal value DoF at vertices
              const FValue val = f.eval_at_node
                (context, var_component, /*dim=*/ 0, // Don't care w/C0
                 vertex, system.time);
              insert_id(id, val, vertex.processor_id());
            }
          else if (cont == C_ONE)
            {
              libmesh_assert(vertex.n_comp(system.number(), var));
              const dof_id_type first_id = vertex.dof_number(system.number(), var, 0);

              // C_ONE elements have a single nodal value and dim
              // gradient values at vertices, as well as cross
              // gradients for HERMITE.  We need to have an element in
              // hand to figure out dim and to have in case this
              // vertex is a new vertex.
              const int dim = elem.dim();
#ifndef NDEBUG
              // For now all C1 elements at a vertex had better have
              // the same dimension.  If anyone hits these asserts let
              // me know; we could probably support a mixed-dimension
              // mesh IFF the 2D elements were all parallel to xy and
              // the 1D elements all parallel to x.
              for (const auto e_id : (*nodes_to_elem)[vertex.id()])
                {
                  const Elem & e = system.get_mesh().elem_ref(e_id);
                  libmesh_assert_equal_to(dim, e.dim());
                }
#endif
#ifdef LIBMESH_ENABLE_AMR
              bool is_parent_vertex = false;
              if (elem.parent())
                {
                  const int i_am_child =
                    elem.parent()->which_child_am_i(&elem);
                  const unsigned int n = elem.get_node_index(&vertex);
                  is_parent_vertex =
                    elem.parent()->is_vertex_on_parent(i_am_child, n);
                }
#else
              const bool is_parent_vertex = false;
#endif

              // The hermite element vertex shape functions are weird
              if (base_fe_type.family == HERMITE)
                {
                  const unsigned int n = elem.get_node_index(&vertex);
                  const FValue val =
                    f.eval_at_node(context,
                                   var_component,
                                   dim,
                                   vertex,
                                   system.time);
                  insert_id(first_id, val, vertex.processor_id());

                  VectorValue<FValue> grad =
                    is_parent_vertex ?
                    g->eval_at_node(context,
                                    var_component,
                                    dim,
                                    vertex,
                                    system.time) :
                    g->eval_at_point(context,
                                     var_component,
                                     vertex,
                                     system.time);
                  // x derivative
                  insert_id(first_id+1, grad(0),
                            vertex.processor_id());
                  if (dim > 1)
                    {
                      // We'll finite difference mixed derivatives
                      Real delta_x = TOLERANCE * elem.hmin();

                      Point nxminus = elem.point(n),
                            nxplus = elem.point(n);
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
                      insert_id(first_id+2, grad(1),
                                vertex.processor_id());
                      // xy derivative
                      insert_id(first_id+3,
                        (gxplus(1) - gxminus(1)) / 2. / delta_x,
                        vertex.processor_id());

                      if (dim > 2)
                        {
                          // z derivative
                          insert_id(first_id+4, grad(2),
                                    vertex.processor_id());
                          // xz derivative
                          insert_id(first_id+5,
                            (gxplus(2) - gxminus(2)) / 2. / delta_x,
                            vertex.processor_id());

                          // We need new points for yz
                          Point nyminus = elem.point(n),
                            nyplus = elem.point(n);
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
                          // yz derivative
                          insert_id(first_id+6,
                            (gyplus(2) - gyminus(2)) / 2. / delta_x,
                            vertex.processor_id());
                          // Getting a 2nd order xyz is more tedious
                          Point nxmym = elem.point(n),
                            nxmyp = elem.point(n),
                            nxpym = elem.point(n),
                            nxpyp = elem.point(n);
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
                          insert_id(first_id+7,
                            (gxzplus - gxzminus) / 2. / delta_x,
                            vertex.processor_id());
                        }
                    }
                }
              else
                {
                  // Currently other C_ONE elements have a single nodal
                  // value shape function and nodal gradient component
                  // shape functions
                  libmesh_assert_equal_to
                    (FEInterface::n_dofs_at_node
                      (dim, base_fe_type, elem.type(),
                       elem.get_node_index(&vertex)),
                    (unsigned int)(1 + dim));
                  const FValue val =
                    f.eval_at_node(context, var_component, dim,
                                   vertex, system.time);
                  insert_id(first_id, val, vertex.processor_id());
                  VectorValue<FValue> grad =
                    is_parent_vertex ?
                    g->eval_at_node(context, var_component, dim,
                                    vertex, system.time) :
                    g->eval_at_point(context, var_component, vertex,
                                     system.time);
                  for (int i=0; i!= dim; ++i)
                    insert_id(first_id + i + 1, grad(i),
                              vertex.processor_id());
                }
            }
          else
            libmesh_error_msg("Unknown continuity " << cont);
        }
    }
  STOP_LOG ("project_vertices","GenericProjector");

  done_saving_ids = sides.empty() && interiors.empty();

  this->send_and_insert_dof_values(ids_to_push, action);
  ids_to_push.clear();

  START_LOG ("project_edges","GenericProjector");
  for (const auto & e_pair : edges)
    {
      const Elem & elem = *std::get<0>(e_pair.second);

      // If this is an unchanged element then we already copied all
      // its dofs
#ifdef LIBMESH_ENABLE_AMR
      if (f.is_grid_projection() &&
          (elem.refinement_flag() != Elem::JUST_REFINED &&
           elem.refinement_flag() != Elem::JUST_COARSENED &&
           elem.p_refinement_flag() != Elem::JUST_REFINED &&
           elem.p_refinement_flag() != Elem::JUST_COARSENED))
        continue;
#endif // LIBMESH_ENABLE_AMR

      const Node & edge_node = *e_pair.first;
      const int dim = elem.dim();

      context.edge = std::get<1>(e_pair.second);
      context.pre_fe_reinit(system, &elem);

      this->find_dofs_to_send(edge_node, ids_to_push);

      // Look at all the variables we're supposed to interpolate from
      // this element on this edge
      for (const auto & var : std::get<2>(e_pair.second))
        {
          const Variable & variable = system.variable(var);
          const FEType & base_fe_type = variable.type();
          const unsigned int var_component =
            system.variable_scalar_number(var, 0);

          if (base_fe_type.family == SCALAR)
            continue;

          FEType fe_type = base_fe_type;

          // This may be a p refined element
          fe_type.order =
            libMesh::Order (fe_type.order + elem.p_level());

          // If this is a Lagrange element with DoFs on edges then by
          // convention we interpolate at the node rather than project
          // along the edge.
          if (fe_type.family == LAGRANGE)
            {
              if (fe_type.order > 1)
                {
                  const dof_id_type dof_id =
                    edge_node.dof_number(system.number(), var, 0);
                  const processor_id_type pid =
                    edge_node.processor_id();
                  FValue fval = f.eval_at_point
                    (context, var_component, edge_node, system.time);
                  insert_id(dof_id, fval, pid);
                }
              continue;
            }

          // If this is a low order monomial element which has merely
          // been h refined then we already copied all its dofs
          if (fe_type.family == MONOMIAL &&
              fe_type.order == CONSTANT &&
              elem.refinement_flag() != Elem::JUST_COARSENED &&
              elem.p_refinement_flag() != Elem::JUST_COARSENED)
            continue;

          // FIXME: Need to generalize this to vector-valued elements. [PB]
          FEBase * fe = nullptr;
          context.get_element_fe( var, fe, dim );
          FEBase * edge_fe = nullptr;
          context.get_edge_fe( var, edge_fe );

          // If we're JUST_COARSENED we'll need a custom
          // evaluation, not just the standard edge FE
          const FEBase & proj_fe =
#ifdef LIBMESH_ENABLE_AMR
            (elem.refinement_flag() == Elem::JUST_COARSENED) ?
            *fe :
#endif
            *edge_fe;

#ifdef LIBMESH_ENABLE_AMR
          if (elem.refinement_flag() == Elem::JUST_COARSENED)
            {
              std::vector<Point> fine_points;

              std::unique_ptr<FEBase> fine_fe
                (FEBase::build (dim, base_fe_type));

              std::unique_ptr<QBase> qrule
                (base_fe_type.default_quadrature_rule(1));
              fine_fe->attach_quadrature_rule(qrule.get());

              const std::vector<Point> & child_xyz =
                        fine_fe->get_xyz();

              for (unsigned int c = 0, nc = elem.n_children();
                   c != nc; ++c)
                {
                  if (!elem.is_child_on_edge(c, context.edge))
                    continue;

                  fine_fe->edge_reinit(elem.child_ptr(c), context.edge);
                  fine_points.insert(fine_points.end(),
                                     child_xyz.begin(),
                                     child_xyz.end());
                }

              std::vector<Point> fine_qp;
              FEInterface::inverse_map (dim, base_fe_type, &elem,
                                        fine_points, fine_qp);

              context.elem_fe_reinit(&fine_qp);
            }
          else
#endif // LIBMESH_ENABLE_AMR
            context.edge_fe_reinit();

          const std::vector<dof_id_type> & dof_indices =
            context.get_dof_indices(var);

          std::vector<unsigned int> edge_dofs;
          FEInterface::dofs_on_edge(&elem, dim, base_fe_type,
                                    context.edge, edge_dofs);

          this->construct_projection
            (dof_indices, edge_dofs, var_component, context,
             &edge_node, insert_id, proj_fe, f, g.get(), ids_to_save);
        }
    }
  STOP_LOG ("project_edges","GenericProjector");

  done_saving_ids = !interiors.empty();

  this->send_and_insert_dof_values(ids_to_push, action);
  ids_to_push.clear();

  START_LOG ("project_sides","GenericProjector");
  for (const auto & s_pair : sides)
    {
      const Elem & elem = *std::get<0>(s_pair.second);

      // If this is an unchanged element then we already copied all
      // its dofs
#ifdef LIBMESH_ENABLE_AMR
      if (f.is_grid_projection() &&
          (elem.refinement_flag() != Elem::JUST_REFINED &&
           elem.refinement_flag() != Elem::JUST_COARSENED &&
           elem.p_refinement_flag() != Elem::JUST_REFINED &&
           elem.p_refinement_flag() != Elem::JUST_COARSENED))
        continue;
#endif // LIBMESH_ENABLE_AMR

      const Node & side_node = *s_pair.first;
      const int dim = elem.dim();

      context.side = std::get<1>(s_pair.second);
      context.pre_fe_reinit(system, &elem);

      this->find_dofs_to_send(side_node, ids_to_push);

      // Look at all the variables we're supposed to interpolate from
      // this element on this side
      for (const auto & var : std::get<2>(s_pair.second))
        {
          const Variable & variable = system.variable(var);
          const FEType & base_fe_type = variable.type();
          const unsigned int var_component =
            system.variable_scalar_number(var, 0);

          if (base_fe_type.family == SCALAR)
            continue;

          FEType fe_type = base_fe_type;

          // This may be a p refined element
          fe_type.order =
            libMesh::Order (fe_type.order + elem.p_level());

          // If this is a Lagrange element with DoFs on sides then by
          // convention we interpolate at the node rather than project
          // along the side.
          if (fe_type.family == LAGRANGE)
            {
              if (fe_type.order > 1)
                {
                  const dof_id_type dof_id =
                    side_node.dof_number(system.number(), var, 0);
                  const processor_id_type pid =
                    side_node.processor_id();
                  FValue fval = f.eval_at_point
                    (context, var_component, side_node, system.time);
                  insert_id(dof_id, fval, pid);
                }
              continue;
            }

          // If this is a low order monomial element which has merely
          // been h refined then we already copied all its dofs
          if (fe_type.family == MONOMIAL &&
              fe_type.order == CONSTANT &&
              elem.refinement_flag() != Elem::JUST_COARSENED &&
              elem.p_refinement_flag() != Elem::JUST_COARSENED)
            continue;

          // FIXME: Need to generalize this to vector-valued elements. [PB]
          FEBase * fe = nullptr;
          context.get_element_fe( var, fe, dim );
          FEBase * side_fe = nullptr;
          context.get_side_fe( var, side_fe );

          // If we're JUST_COARSENED we'll need a custom
          // evaluation, not just the standard side FE
          const FEBase & proj_fe =
#ifdef LIBMESH_ENABLE_AMR
            (elem.refinement_flag() == Elem::JUST_COARSENED) ?
            *fe :
#endif
            *side_fe;

#ifdef LIBMESH_ENABLE_AMR
          if (elem.refinement_flag() == Elem::JUST_COARSENED)
            {
              std::vector<Point> fine_points;

              std::unique_ptr<FEBase> fine_fe
                (FEBase::build (dim, base_fe_type));

              std::unique_ptr<QBase> qrule
                (base_fe_type.default_quadrature_rule(1));
              fine_fe->attach_quadrature_rule(qrule.get());

              const std::vector<Point> & child_xyz =
                        fine_fe->get_xyz();

              for (unsigned int c = 0, nc = elem.n_children();
                   c != nc; ++c)
                {
                  if (!elem.is_child_on_side(c, context.side))
                    continue;

                  fine_fe->reinit(elem.child_ptr(c), context.side);
                  fine_points.insert(fine_points.end(),
                                     child_xyz.begin(),
                                     child_xyz.end());
                }

              std::vector<Point> fine_qp;
              FEInterface::inverse_map (dim, base_fe_type, &elem,
                                        fine_points, fine_qp);

              context.elem_fe_reinit(&fine_qp);
            }
          else
#endif // LIBMESH_ENABLE_AMR
            context.side_fe_reinit();

          const std::vector<dof_id_type> & dof_indices =
            context.get_dof_indices(var);

          std::vector<unsigned int> side_dofs;
          FEInterface::dofs_on_side(&elem, dim, base_fe_type,
                                    context.side, side_dofs);

          this->construct_projection
            (dof_indices, side_dofs, var_component, context,
             &side_node, insert_id, proj_fe, f, g.get(), ids_to_save);
        }
    }
  STOP_LOG ("project_sides","GenericProjector");

  done_saving_ids = true;

  this->send_and_insert_dof_values(ids_to_push, action);
  ids_to_push.clear();

  START_LOG ("project_interiors","GenericProjector");
  // Iterate over all dof-bearing element interiors in the range
  for (const auto & elem : interiors)
    {
      unsigned char dim = cast_int<unsigned char>(elem->dim());

      context.pre_fe_reinit(system, elem);

      // Loop over all the variables we've been requested to project, to
      // do the projection
      for (const auto & var : variables)
        {
          const Variable & variable = system.variable(var);

          if (!variable.active_on_subdomain(elem->subdomain_id()))
            continue;

          const FEType & base_fe_type = variable.type();

          if (base_fe_type.family == SCALAR)
            continue;

          FEBase * fe = nullptr;
          context.get_element_fe( var, fe, dim );

          FEType fe_type = base_fe_type;

          // This may be a p refined element
          fe_type.order =
            libMesh::Order (fe_type.order + elem->p_level());

          const unsigned int var_component =
            system.variable_scalar_number(var, 0);

          // If this is a Lagrange element with interior DoFs then by
          // convention we interpolate at the interior node rather
          // than project along the interior.
          if (fe_type.family == LAGRANGE)
            {
              if (fe_type.order > 1)
                {
                  const unsigned int first_interior_node =
                    (elem->n_vertices() +
                     ((elem->dim() > 2) * elem->n_edges()) +
                     ((elem->dim() > 1) * elem->n_sides()));
                  const unsigned int n_nodes = elem->n_nodes();

                  // < vs != is important here for HEX20, QUAD8!
                  for (unsigned int n = first_interior_node; n < n_nodes; ++n)
                    {
                      const Node & interior_node = elem->node_ref(n);
                      const dof_id_type dof_id =
                        interior_node.dof_number(system.number(), var, 0);
                      const processor_id_type pid =
                        interior_node.processor_id();
                      FValue fval = f.eval_at_point
                        (context, var_component, interior_node, system.time);
                      insert_id(dof_id, fval, pid);
                    }
                }
              continue;
            }

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

          const std::vector<dof_id_type> & dof_indices =
            context.get_dof_indices(var);

          const unsigned int n_dofs =
            cast_int<unsigned int>(dof_indices.size());

          std::vector<unsigned int> all_dofs(n_dofs);
          std::iota(all_dofs.begin(), all_dofs.end(), 0);

          this->construct_projection
            (dof_indices, all_dofs, var_component, context, nullptr,
             insert_id, *fe, f, g.get(), ids_to_save);
        } // end variables loop
    } // end elements loop
  STOP_LOG ("project_interiors","GenericProjector");
}



template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
void
GenericProjector<FFunctor, GFunctor, FValue,
                 ProjectionAction>::find_dofs_to_send
  (const Node & node,
   std::unordered_map<dof_id_type, std::pair<FValue, processor_id_type>> & ids_to_push) const
{
  // Any ghosted node in our range might have an owner who needs our
  // data
  const processor_id_type owner = node.processor_id();
  if (owner != system.processor_id())
    {
      const MeshBase & mesh = system.get_mesh();
      const DofMap & dof_map = system.get_dof_map();

      // But let's check and see if we can be certain the owner can
      // compute any or all of its own dof coefficients on that node
      std::vector<dof_id_type> node_dof_ids, elem_dof_ids;
      dof_map.dof_indices(&node, node_dof_ids);
      libmesh_assert(std::is_sorted(node_dof_ids.begin(),
                                    node_dof_ids.end()));
      const std::vector<dof_id_type> & patch =
        (*nodes_to_elem)[node.id()];
      for (const auto & elem_id : patch)
        {
          const Elem * elem = mesh.query_elem_ptr(elem_id);
          if (!elem->active())
            continue;
          dof_map.dof_indices(elem, elem_dof_ids);
          std::sort(elem_dof_ids.begin(), elem_dof_ids.end());

          std::vector<dof_id_type> diff_ids(node_dof_ids.size());
          auto it = std::set_difference(node_dof_ids.begin(), node_dof_ids.end(),
                                        elem_dof_ids.begin(), elem_dof_ids.end(), diff_ids.begin());
          diff_ids.resize(it-diff_ids.begin());
          node_dof_ids.swap(diff_ids);
          if (node_dof_ids.empty())
            break;
        }

      // Give ids_to_push default invalid pid: not yet computed
      for (auto id : node_dof_ids)
        ids_to_push[id].second = DofObject::invalid_processor_id;
    }
}



template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
void
GenericProjector<FFunctor, GFunctor, FValue,
                 ProjectionAction>::send_and_insert_dof_values
  (std::unordered_map<dof_id_type, std::pair<FValue, processor_id_type>> & ids_to_push,
   ProjectionAction & action) const
{
  // See if we calculated any ids that need to be pushed; get them
  // ready to push.
  std::unordered_map<processor_id_type, std::vector<dof_id_type>>
    pushed_dof_ids, received_dof_ids;
  std::unordered_map<processor_id_type, std::vector<typename TypeToSend<FValue>::type>>
    pushed_dof_values, received_dof_values;
  for (auto & id_val_pid : ids_to_push)
    {
      processor_id_type pid = id_val_pid.second.second;
      if (pid != DofObject::invalid_processor_id)
        {
          pushed_dof_ids[pid].push_back(id_val_pid.first);
          pushed_dof_values[pid].push_back(convert_to_send(id_val_pid.second.first));
        }
    }

  // If and when we get ids pushed to us, act on them if we have
  // corresponding values or save them if not
  auto ids_action_functor =
    [&action, &received_dof_ids, &received_dof_values]
    (processor_id_type pid,
     const std::vector<dof_id_type> & data)
    {
      auto iter = received_dof_values.find(pid);
      if (iter == received_dof_values.end())
        {
          libmesh_assert(received_dof_ids.find(pid) ==
                         received_dof_ids.end());
          received_dof_ids[pid] = data;
        }
      else
        {
          auto & vals = iter->second;
          std::size_t vals_size = vals.size();
          libmesh_assert_equal_to(vals_size, data.size());
          for (std::size_t i=0; i != vals_size; ++i)
            {
              FValue val;
              convert_from_receive(vals[i], val);
              action.insert(data[i], val);
            }
          received_dof_values.erase(iter);
        }
    };

  // If and when we get values pushed to us, act on them if we have
  // corresponding ids or save them if not
  auto values_action_functor =
    [&action, &received_dof_ids, &received_dof_values]
    (processor_id_type pid,
     const std::vector<typename TypeToSend<FValue>::type> & data)
    {
      auto iter = received_dof_ids.find(pid);
      if (iter == received_dof_ids.end())
        {
          // We get no more than 1 values vector from anywhere
          libmesh_assert(received_dof_values.find(pid) ==
                         received_dof_values.end());
          received_dof_values[pid] = data;
        }
      else
        {
          auto & ids = iter->second;
          std::size_t ids_size = ids.size();
          libmesh_assert_equal_to(ids_size, data.size());
          for (std::size_t i=0; i != ids_size; ++i)
            {
              FValue val;
              convert_from_receive(data[i], val);
              action.insert(ids[i], val);
            }
          received_dof_ids.erase(iter);
        }
    };

  Parallel::push_parallel_vector_data
    (system.comm(), pushed_dof_ids, ids_action_functor);

  Parallel::push_parallel_vector_data
    (system.comm(), pushed_dof_values, values_action_functor);

  // At this point we shouldn't have any unprocessed data left
  libmesh_assert(received_dof_ids.empty());
  libmesh_assert(received_dof_values.empty());

}



template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
template <typename InsertId>
void
GenericProjector<FFunctor, GFunctor, FValue,
                 ProjectionAction>::construct_projection
  (const std::vector<dof_id_type> & dof_indices_var,
   const std::vector<unsigned int> & involved_dofs,
   unsigned int var_component,
   const FEMContext & context,
   const Node * node,
   InsertId & insert_id,
   const FEBase & fe,
   FFunctor & f,
   GFunctor * g,
   const std::unordered_map<dof_id_type, FValue> & ids_to_save
   ) const
{
  const std::vector<Real> & JxW = fe.get_JxW();
  const std::vector<std::vector<Real>> & phi = fe.get_phi();
  const std::vector<std::vector<RealGradient>> * dphi = nullptr;
  const std::vector<Point> & xyz_values = fe.get_xyz();
  const FEContinuity cont = fe.get_continuity();

  if (cont == C_ONE)
    dphi = &(fe.get_dphi());

  const unsigned int n_involved_dofs =
    cast_int<unsigned int>(involved_dofs.size());

  std::vector<dof_id_type> free_dof_ids;
  DenseVector<FValue> Uinvolved(n_involved_dofs);
  std::vector<char> dof_is_fixed(n_involved_dofs, false); // bools

  for (auto i : IntRange<unsigned int>(0, n_involved_dofs))
    {
      const dof_id_type id = dof_indices_var[involved_dofs[i]];
      auto iter = ids_to_save.find(id);
      if (iter == ids_to_save.end())
        free_dof_ids.push_back(id);
      else
        {
          dof_is_fixed[i] = true;
          Uinvolved(i) = iter->second;
        }
    }

  const unsigned int free_dofs = free_dof_ids.size();

  // There may be nothing to project
  if (!free_dofs)
    return;

  // The element matrix and RHS for projections.
  // Note that Ke is always real-valued, whereas
  // Fe may be complex valued if complex number
  // support is enabled
  DenseMatrix<Real> Ke(free_dofs, free_dofs);
  DenseVector<FValue> Fe(free_dofs);
  // The new degree of freedom coefficients to solve for
  DenseVector<FValue> Ufree(free_dofs);

  const unsigned int n_qp =
    cast_int<unsigned int>(xyz_values.size());

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
      for (unsigned int sidei=0, freei=0;
           sidei != n_involved_dofs; ++sidei)
        {
          unsigned int i = involved_dofs[sidei];
          // fixed DoFs aren't test functions
          if (dof_is_fixed[sidei])
            continue;
          for (unsigned int sidej=0, freej=0;
               sidej != n_involved_dofs; ++sidej)
            {
              unsigned int j = involved_dofs[sidej];
              if (dof_is_fixed[sidej])
                Fe(freei) -= phi[i][qp] * phi[j][qp] *
                  JxW[qp] * Uinvolved(sidej);
              else
                Ke(freei,freej) += phi[i][qp] *
                  phi[j][qp] * JxW[qp];
              if (cont == C_ONE)
                {
                  if (dof_is_fixed[sidej])
                    Fe(freei) -= ( (*dphi)[i][qp] *
                                   (*dphi)[j][qp] ) *
                      JxW[qp] * Uinvolved(sidej);
                  else
                    Ke(freei,freej) += ( (*dphi)[i][qp] *
                                         (*dphi)[j][qp] )
                      * JxW[qp];
                }
              if (!dof_is_fixed[sidej])
                freej++;
            }
          Fe(freei) += phi[i][qp] * fineval * JxW[qp];
          if (cont == C_ONE)
            Fe(freei) += (finegrad * (*dphi)[i][qp] ) *
              JxW[qp];
          freei++;
        }
    }

  Ke.cholesky_solve(Fe, Ufree);

  // Transfer new edge solutions to element
  const processor_id_type pid = node ?
    node->processor_id() : DofObject::invalid_processor_id;
  for (unsigned int i=0; i != free_dofs; ++i)
    insert_id(free_dof_ids[i], Ufree(i), pid);
}


} // namespace libMesh

#endif // GENERIC_PROJECTOR_H
