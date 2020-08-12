// The libMesh Finite Element Library.
// Copyright (C) 2002-2020 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// libMesh includes
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/int_range.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature.h"
#include "libmesh/system.h"
#include "libmesh/threads.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique
#include "libmesh/tensor_tools.h"
#include "libmesh/libmesh_common.h"

// TIMPI includes
#include "timpi/parallel_sync.h"


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
  bool done_saving_ids;

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

  // The single "project" function handles all the threaded projection
  // calculations and intervening MPI communications
  void project(const ConstElemRange & range);

  // We generally need to hang on to every value we've calculated
  // until we're all done, because later projection calculations
  // depend on boundary data from earlier calculations.
  std::unordered_map<dof_id_type, typename FFunctor::ValuePushType> ids_to_save;

  // This needs to be sorted so we can get sorted dof indices cheaply
  // later
  typedef std::set<unsigned int> var_set;

  // We now need to divide our threadable work into stages, to allow
  // for the potential necessity of MPI communication in between them.
  struct SubFunctor {
    GenericProjector & projector;

    SubFunctor (GenericProjector & p);

    // While we're working on nodes, also figure out which ghosted
    // nodes will have data we might need to send to their owners
    // instead of being acted on by ourselves.
    //
    // We keep track of which dof ids we might need to send, and what
    // values those ids should get (along with a pprocessor_id to
    // leave invalid in case *we* can't compute those values either).
    std::unordered_map<dof_id_type,
                       std::pair<typename FFunctor::ValuePushType,
                                 processor_id_type>> new_ids_to_push;

    // We'll hang on to any new ids to save on a per-thread basis
    // because we won't need them until subsequent phases
    std::unordered_map<dof_id_type, typename FFunctor::ValuePushType> new_ids_to_save;

    // Helper function for filling new_ids_to_push
    void find_dofs_to_send (const Node & node,
                            const Elem & elem,
                            unsigned short node_num,
                            const var_set & vars);

    // When we have new data to act on, we may also need to save it
    // and get ready to push it.
    template <typename InsertInput,
              typename std::enable_if<
                  std::is_same<typename ProjectionAction::InsertInput, InsertInput>::value,
                  int>::type = 0>
    void insert_id(dof_id_type id, const InsertInput & val, processor_id_type pid);

    template <typename InsertInput,
              typename std::enable_if<
                  !std::is_same<typename ProjectionAction::InsertInput, InsertInput>::value,
                  int>::type = 0>
    void insert_id(dof_id_type id, const InsertInput & val, processor_id_type pid);

    template <typename InsertInput,
              typename std::enable_if<
                  std::is_same<typename ProjectionAction::InsertInput, InsertInput>::value,
                  int>::type = 0>
    void insert_ids(const std::vector<dof_id_type> & ids,
                    const std::vector<InsertInput> & vals,
                    processor_id_type pid);

    template <typename InsertInput,
              typename std::enable_if<
                  !std::is_same<typename ProjectionAction::InsertInput, InsertInput>::value,
                  int>::type = 0>
    void insert_ids(const std::vector<dof_id_type> & ids,
                    const std::vector<InsertInput> & vals,
                    processor_id_type pid);

    // Our subclasses will need to merge saved ids
    void join (const SubFunctor & other);

  protected:
    // For thread safety and efficiency we'll want thread-local copies
    // of various functors
    ProjectionAction action;
    FFunctor f;

    // Each stage of the work we do will require us to prepare
    // FEMContext objects to assist in the work.
    FEMContext context;

    // These get used often enough to cache them
    std::vector<FEContinuity> conts;
    std::vector<FEFieldType> field_types;

    const System & system;
  };


  typedef std::pair<const Node *, std::tuple<const Elem *, unsigned short, var_set>>
    node_projection;

  typedef StoredRange<std::vector<node_projection>::const_iterator,
                      node_projection> node_range;


  struct SubProjector : public SubFunctor {
    SubProjector (GenericProjector & p);

    using SubFunctor::action;
    using SubFunctor::f;
    using SubFunctor::system;
    using SubFunctor::context;
    using SubFunctor::conts;
    using SubFunctor::field_types;
    using SubFunctor::insert_id;
    using SubFunctor::insert_ids;

  protected:
    // Projections of C1 elements require a gradient as well
    std::unique_ptr<GFunctor> g;

    void construct_projection
      (const std::vector<dof_id_type> & dof_indices_var,
       const std::vector<unsigned int> & involved_dofs,
       unsigned int var_component,
       const Node * node,
       const FEGenericBase<typename FFunctor::RealType> & fe);
  };


  // First we'll copy DoFs where we can, while sorting remaining DoFs
  // for interpolation and projection later.
  struct SortAndCopy : public SubFunctor {
    SortAndCopy (GenericProjector & p) : SubFunctor(p) {}

    SortAndCopy (SortAndCopy & other, Threads::split) : SubFunctor(other.projector) {}

    using SubFunctor::action;
    using SubFunctor::f;
    using SubFunctor::system;
    using SubFunctor::context;
    using SubFunctor::insert_ids;

    void operator() (const ConstElemRange & range);

    // We need to merge saved multimaps when working from multiple
    // threads
    void join (const SortAndCopy & other);

    // For vertices, edges and sides, we need to know what element,
    // what local vertex/edge/side number, and what set of variables
    // to evaluate.
    typedef std::unordered_multimap<const Node *, std::tuple<const Elem *, unsigned short, var_set>> ves_multimap;

    ves_multimap vertices, edges, sides;
    std::vector<const Elem *> interiors;
  };


  struct ProjectVertices : public SubProjector {
    ProjectVertices (GenericProjector & p) : SubProjector(p) {}

    ProjectVertices (ProjectVertices & p_v, Threads::split) : SubProjector(p_v.projector) {}

    using SubProjector::action;
    using SubProjector::f;
    using SubProjector::g;
    using SubProjector::system;
    using SubProjector::context;
    using SubProjector::conts;
    using SubFunctor::field_types;

    using SubProjector::insert_id;
    using SubProjector::insert_ids;

    void operator() (const node_range & range);
  };


  struct ProjectEdges : public SubProjector {
    ProjectEdges (GenericProjector & p) : SubProjector(p) {}

    ProjectEdges (ProjectEdges & p_e, Threads::split) : SubProjector(p_e.projector) {}

    using SubProjector::action;
    using SubProjector::f;
    using SubProjector::g;
    using SubProjector::system;
    using SubProjector::context;
    using SubProjector::conts;
    using SubFunctor::field_types;

    using SubProjector::insert_id;
    using SubProjector::insert_ids;

    void operator() (const node_range & range);
  };

  struct ProjectSides : public SubProjector {
    ProjectSides (GenericProjector & p) : SubProjector(p) {}

    ProjectSides (ProjectSides & p_s, Threads::split) : SubProjector(p_s.projector) {}

    using SubProjector::action;
    using SubProjector::f;
    using SubProjector::g;
    using SubProjector::system;
    using SubProjector::context;
    using SubProjector::conts;
    using SubFunctor::field_types;

    using SubProjector::insert_id;
    using SubProjector::insert_ids;

    void operator() (const node_range & range);
  };


  typedef StoredRange<std::vector<const Elem *>::const_iterator,
                      const Elem *> interior_range;

  struct ProjectInteriors : public SubProjector {
    ProjectInteriors (GenericProjector & p) : SubProjector(p) {}

    ProjectInteriors (ProjectInteriors & p_i, Threads::split) : SubProjector(p_i.projector) {}

    using SubProjector::action;
    using SubProjector::f;
    using SubProjector::g;
    using SubProjector::system;
    using SubProjector::context;
    using SubProjector::conts;
    using SubFunctor::field_types;

    using SubProjector::insert_id;
    using SubProjector::insert_ids;

    void operator() (const interior_range & range);
  };

  template <typename Value>
  void send_and_insert_dof_values
    (std::unordered_map<dof_id_type, std::pair<Value, processor_id_type>> & ids_to_push,
     ProjectionAction & action) const;
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
  typedef Val InsertInput;

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


  void insert(const std::vector<dof_id_type> & dof_indices,
              const DenseVector<Val> & Ue)
  {
    const numeric_index_type
      first = target_vector.first_local_index(),
      last  = target_vector.last_local_index();

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
protected:
  typedef typename TensorTools::MakeBaseNumber<Output>::type DofValueType;

public:
  typedef typename TensorTools::MakeReal<Output>::type RealType;
  typedef DofValueType ValuePushType;
  typedef Output FunctorValue;

  FEMFunctionWrapper(const FEMFunctionBase<Output> & f) : _f(f.clone()) {}

  FEMFunctionWrapper(const FEMFunctionWrapper<Output> & fw) :
    _f(fw._f->clone()) {}

  void init_context (FEMContext & c) { _f->init_context(c); }

  Output eval_at_node (const FEMContext & c,
                       unsigned int i,
                       unsigned int /*elem_dim*/,
                       const Node & n,
                       bool /*extra_hanging_dofs*/,
                       const Real time)
  { return _f->component(c, i, n, time); }

  Output eval_at_point (const FEMContext & c,
                        unsigned int i,
                        const Point & n,
                        const Real time,
                        bool /*skip_context_check*/)
  { return _f->component(c, i, n, time); }

  void eval_mixed_derivatives (const FEMContext & /*c*/,
                               unsigned int /*i*/,
                               unsigned int /*dim*/,
                               const Node & /*n*/,
                               std::vector<Output> & /*derivs*/)
  { libmesh_error(); } // this is only for grid projections

  bool is_grid_projection() { return false; }

  void eval_old_dofs (const Elem &,
                      unsigned int,
                      unsigned int,
                      std::vector<dof_id_type> &,
                      std::vector<Output> &)
  { libmesh_error(); }

  void eval_old_dofs (const Elem &,
                      const FEType &,
                      unsigned int,
                      unsigned int,
                      std::vector<dof_id_type> &,
                      std::vector<Output> &)
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
protected:
  typedef typename TensorTools::MakeBaseNumber<Output>::type DofValueType;
public:
  typedef typename TensorTools::MakeReal<Output>::type RealType;

  OldSolutionBase(const libMesh::System & sys_in) :
    last_elem(nullptr),
    sys(sys_in),
    old_context(sys_in)
  {
    // We'll be queried for components but we'll typically be looking
    // up data by variables, and those indices don't always match
    for (auto v : make_range(sys.n_vars()))
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

  static void get_shape_outputs(FEAbstract & fe);

  // Integrating on new mesh elements, we won't yet have an up to date
  // current_local_solution.
  void init_context (FEMContext & c)
  {
    c.set_algebraic_type(FEMContext::DOFS_ONLY);

    // Loop over variables, to prerequest
    for (auto var : make_range(sys.n_vars()))
      {
        FEAbstract * fe = nullptr;
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
  using typename OldSolutionBase<Output, point_output>::RealType;
  using typename OldSolutionBase<Output, point_output>::DofValueType;
  typedef Output FunctorValue;
  typedef DofValueType ValuePushType;

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
                       bool /* extra_hanging_dofs */,
                       Real /* time */ =0.);

  Output eval_at_point(const FEMContext & c,
                       unsigned int i,
                       const Point & p,
                       Real /* time */,
                       bool skip_context_check)
  {
    LOG_SCOPE ("eval_at_point()", "OldSolutionValue");

    if (!skip_context_check)
      if (!this->check_old_context(c, p))
        return Output(0);

    // Handle offset from non-scalar components in previous variables
    libmesh_assert_less(i, this->component_to_var.size());
    unsigned int var = this->component_to_var[i];

    Output n;
    (this->old_context.*point_output)(var, p, n, this->out_of_elem_tol);
    return n;
  }

  template <typename T = Output,
            typename std::enable_if<std::is_same<T, Number>::value, int>::type = 0>
  void eval_mixed_derivatives(const FEMContext & libmesh_dbg_var(c),
                              unsigned int i,
                              unsigned int dim,
                              const Node & n,
                              std::vector<Output> & derivs)
  {
    LOG_SCOPE ("eval_mixed_derivatives", "OldSolutionValue");

    // This should only be called on vertices
    libmesh_assert_less(c.get_elem().get_node_index(&n),
                        c.get_elem().n_vertices());

    // Handle offset from non-scalar components in previous variables
    libmesh_assert_less(i, this->component_to_var.size());
    unsigned int var = this->component_to_var[i];

    // We have 1 mixed derivative in 2D, 4 in 3D
    const unsigned int n_mixed = (dim-1) * (dim-1);
    derivs.resize(n_mixed);

    // Be sure to handle cases where the variable wasn't defined on
    // this node (e.g. due to changing subdomain support)
    if (n.old_dof_object &&
        n.old_dof_object->n_vars(this->sys.number()) &&
        n.old_dof_object->n_comp(this->sys.number(), var))
      {
        const dof_id_type first_old_id =
          n.old_dof_object->dof_number(this->sys.number(), var, dim);
        std::vector<dof_id_type> old_ids(n_mixed);
        std::iota(old_ids.begin(), old_ids.end(), first_old_id);
        old_solution.get(old_ids, derivs);
      }
    else
      {
        std::fill(derivs.begin(), derivs.end(), 0);
      }
  }

  template <typename T = Output,
            typename std::enable_if<!std::is_same<T, Number>::value, int>::type = 0>
  void eval_mixed_derivatives(
      const FEMContext &, unsigned int, unsigned int, const Node &, std::vector<Output> &)
  {
    libmesh_error_msg("eval_mixed_derivatives should only be applicable for Hermite finite element "
                      "types. I don't know how you got here");
  }

  void eval_old_dofs (const Elem & elem,
                      unsigned int node_num,
                      unsigned int var_num,
                      std::vector<dof_id_type> & indices,
                      std::vector<DofValueType> & values)
  {
    LOG_SCOPE ("eval_old_dofs(node)", "OldSolutionValue");

    this->sys.get_dof_map().dof_indices(elem, node_num, indices, var_num);

    std::vector<dof_id_type> old_indices;

    this->sys.get_dof_map().old_dof_indices(elem, node_num, old_indices, var_num);

    libmesh_assert_equal_to (old_indices.size(), indices.size());

    // We may have invalid_id in cases where no old DoF existed, e.g.
    // due to expansion of a subdomain-restricted variable's subdomain
    bool invalid_old_index = false;
    for (const auto & di : old_indices)
      if (di == DofObject::invalid_id)
        invalid_old_index = true;

    values.resize(old_indices.size());
    if (invalid_old_index)
      {
        for (auto i : index_range(old_indices))
          {
            const dof_id_type di = old_indices[i];
            if (di == DofObject::invalid_id)
              values[i] = 0;
            else
              values[i] = old_solution(di);
          }
      }
    else
      old_solution.get(old_indices, values);
  }


  void eval_old_dofs (const Elem & elem,
                      const FEType & fe_type,
                      unsigned int sys_num,
                      unsigned int var_num,
                      std::vector<dof_id_type> & indices,
                      std::vector<DofValueType> & values)
  {
    LOG_SCOPE ("eval_old_dofs(elem)", "OldSolutionValue");

    // We're only to be asked for old dofs on elements that can copy
    // them through DO_NOTHING or through refinement.
    const Elem & old_elem =
      (elem.refinement_flag() == Elem::JUST_REFINED) ?
      *elem.parent() : elem;

    // If there are any element-based DOF numbers, get them
    const unsigned int nc = FEInterface::n_dofs_per_elem(fe_type, &elem);

    std::vector<dof_id_type> old_dof_indices(nc);
    indices.resize(nc);

    // We should never have fewer dofs than necessary on an
    // element unless we're getting indices on a parent element,
    // and we should never need those indices
    if (nc != 0)
      {
        libmesh_assert(old_elem.old_dof_object);

        const std::pair<unsigned int, unsigned int>
          vg_and_offset = elem.var_to_vg_and_offset(sys_num,var_num);
        const unsigned int vg = vg_and_offset.first;
        const unsigned int vig = vg_and_offset.second;

        const unsigned int n_comp = elem.n_comp_group(sys_num,vg);
        libmesh_assert_greater(elem.n_systems(), sys_num);
        libmesh_assert_greater_equal(n_comp, nc);

        for (unsigned int i=0; i<nc; i++)
          {
            const dof_id_type d_old =
              old_elem.old_dof_object->dof_number(sys_num, vg, vig, i, n_comp);
            const dof_id_type d_new =
              elem.dof_number(sys_num, vg, vig, i, n_comp);
            libmesh_assert_not_equal_to (d_old, DofObject::invalid_id);
            libmesh_assert_not_equal_to (d_new, DofObject::invalid_id);

            old_dof_indices[i] = d_old;
            indices[i] = d_new;
          }
      }

    values.resize(nc);

    old_solution.get(old_dof_indices, values);
  }

private:
  const NumericVector<Number> & old_solution;
};

template<>
inline void
OldSolutionBase<Number, &FEMContext::point_value>::get_shape_outputs(FEAbstract & fe)
{
  fe.request_phi();
}


template<>
inline void
OldSolutionBase<Gradient, &FEMContext::point_gradient>::get_shape_outputs(FEAbstract & fe)
{
  fe.request_dphi();
}

template<>
inline void
OldSolutionBase<Gradient, &FEMContext::point_value>::get_shape_outputs(FEAbstract & fe)
{
  fe.request_phi();
}


template<>
inline void
OldSolutionBase<Tensor, &FEMContext::point_gradient>::get_shape_outputs(FEAbstract & fe)
{
  fe.request_dphi();
}


#ifdef LIBMESH_USE_COMPLEX_NUMBERS
template<>
inline void
OldSolutionBase<Real, &FEMContext::point_value>::get_shape_outputs(FEAbstract & fe)
{
  fe.request_phi();
}


template<>
inline void
OldSolutionBase<RealGradient, &FEMContext::point_gradient>::get_shape_outputs(FEAbstract & fe)
{
  fe.request_dphi();
}

template<>
inline void
OldSolutionBase<RealGradient, &FEMContext::point_value>::get_shape_outputs(FEAbstract & fe)
{
  fe.request_phi();
}


template<>
inline void
OldSolutionBase<RealTensor, &FEMContext::point_gradient>::get_shape_outputs(FEAbstract & fe)
{
  fe.request_dphi();
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
             bool extra_hanging_dofs,
             Real /* time */)
{
  LOG_SCOPE ("Number eval_at_node()", "OldSolutionValue");

  // This should only be called on vertices
  libmesh_assert_less(c.get_elem().get_node_index(&n),
                      c.get_elem().n_vertices());

  // Handle offset from non-scalar components in previous variables
  libmesh_assert_less(i, this->component_to_var.size());
  unsigned int var = this->component_to_var[i];

  // Optimize for the common case, where this node was part of the
  // old solution.
  //
  // Be sure to handle cases where the variable wasn't defined on
  // this node (due to changing subdomain support) or where the
  // variable has no components on this node (due to Elem order
  // exceeding FE order) or where the old_dof_object dofs might
  // correspond to non-vertex dofs (due to extra_hanging_dofs and
  // refinement)

  const Elem::RefinementState flag = c.get_elem().refinement_flag();

  if (n.old_dof_object &&
      (!extra_hanging_dofs ||
       flag == Elem::JUST_COARSENED ||
       flag == Elem::DO_NOTHING) &&
      n.old_dof_object->n_vars(sys.number()) &&
      n.old_dof_object->n_comp(sys.number(), var))
    {
      const dof_id_type old_id =
        n.old_dof_object->dof_number(sys.number(), var, 0);
      return old_solution(old_id);
    }

  return this->eval_at_point(c, i, n, 0, false);
}

template <>
inline
Gradient
OldSolutionValue<Gradient, &FEMContext::point_value>::
eval_at_node(const FEMContext & c,
             unsigned int i,
             unsigned int /* elem_dim */,
             const Node & n,
             bool extra_hanging_dofs,
             Real /* time */)
{
  LOG_SCOPE ("Number eval_at_node()", "OldSolutionValue");

  // This should only be called on vertices
  libmesh_assert_less(c.get_elem().get_node_index(&n),
                      c.get_elem().n_vertices());

  // Handle offset from non-scalar components in previous variables
  libmesh_assert_less(i, this->component_to_var.size());
  unsigned int var = this->component_to_var[i];

  // Optimize for the common case, where this node was part of the
  // old solution.
  //
  // Be sure to handle cases where the variable wasn't defined on
  // this node (due to changing subdomain support) or where the
  // variable has no components on this node (due to Elem order
  // exceeding FE order) or where the old_dof_object dofs might
  // correspond to non-vertex dofs (due to extra_hanging_dofs and
  // refinement)

  const auto & elem = c.get_elem();

  const Elem::RefinementState flag = elem.refinement_flag();

  if (n.old_dof_object &&
      (!extra_hanging_dofs ||
       flag == Elem::JUST_COARSENED ||
       flag == Elem::DO_NOTHING) &&
      n.old_dof_object->n_vars(sys.number()) &&
      n.old_dof_object->n_comp(sys.number(), var))
    {
      Gradient return_val;

      for (unsigned int dim = 0; dim < elem.dim(); ++dim)
      {
        const dof_id_type old_id =
          n.old_dof_object->dof_number(sys.number(), var, dim);
        return_val(dim) = old_solution(old_id);
      }

      return return_val;
    }

  return this->eval_at_point(c, i, n, 0, false);
}



template<>
inline
Gradient
OldSolutionValue<Gradient, &FEMContext::point_gradient>::
eval_at_node(const FEMContext & c,
             unsigned int i,
             unsigned int elem_dim,
             const Node & n,
             bool extra_hanging_dofs,
             Real /* time */)
{
  LOG_SCOPE ("Gradient eval_at_node()", "OldSolutionValue");

  // This should only be called on vertices
  libmesh_assert_less(c.get_elem().get_node_index(&n),
                      c.get_elem().n_vertices());

  // Handle offset from non-scalar components in previous variables
  libmesh_assert_less(i, this->component_to_var.size());
  unsigned int var = this->component_to_var[i];

  // Optimize for the common case, where this node was part of the
  // old solution.
  //
  // Be sure to handle cases where the variable wasn't defined on
  // this node (due to changing subdomain support) or where the
  // variable has no components on this node (due to Elem order
  // exceeding FE order) or where the old_dof_object dofs might
  // correspond to non-vertex dofs (due to extra_hanging_dofs and
  // refinement)

  const Elem::RefinementState flag = c.get_elem().refinement_flag();

  if (n.old_dof_object &&
      (!extra_hanging_dofs ||
       flag == Elem::JUST_COARSENED ||
       flag == Elem::DO_NOTHING) &&
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

  return this->eval_at_point(c, i, n, 0, false);
}

template<>
inline
Tensor
OldSolutionValue<Tensor, &FEMContext::point_gradient>::
eval_at_node(const FEMContext &,
             unsigned int,
             unsigned int,
             const Node &,
             bool,
             Real)
{
  libmesh_error_msg("You shouldn't need to call eval_at_node for the gradient "
                    "functor for a vector-valued finite element type");
  return Tensor();
}

template <typename Output,
          void (FEMContext::*point_output) (unsigned int,
                                            const Point &,
                                            Output &,
                                            const Real) const>
const Real OldSolutionBase<Output, point_output>::out_of_elem_tol = 10 * TOLERANCE;

#endif // LIBMESH_ENABLE_AMR

/**
 * Function definitions
 */

template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
void GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::project
  (const ConstElemRange & range)
{
  LOG_SCOPE ("project", "GenericProjector");

  // Unless we split sort and copy into two passes we can't know for
  // sure ahead of time whether we need to save the copied ids
  done_saving_ids = false;

  SortAndCopy sort_work(*this);
  Threads::parallel_reduce (range, sort_work);
  ProjectionAction action(master_action);

  // Keep track of dof ids and values to send to other ranks
  std::unordered_map<dof_id_type, std::pair<typename FFunctor::ValuePushType, processor_id_type>>
      ids_to_push;

  ids_to_push.insert(sort_work.new_ids_to_push.begin(),
                     sort_work.new_ids_to_push.end());
  ids_to_save.insert(sort_work.new_ids_to_save.begin(),
                     sort_work.new_ids_to_save.end());

  std::vector<node_projection> vertices(sort_work.vertices.begin(),
                                        sort_work.vertices.end());

  done_saving_ids = sort_work.edges.empty() &&
    sort_work.sides.empty() && sort_work.interiors.empty();
  system.comm().max(done_saving_ids);

  {
    ProjectVertices project_vertices(*this);
    Threads::parallel_reduce (node_range(&vertices), project_vertices);
    ids_to_push.insert(project_vertices.new_ids_to_push.begin(),
                       project_vertices.new_ids_to_push.end());
    ids_to_save.insert(project_vertices.new_ids_to_save.begin(),
                       project_vertices.new_ids_to_save.end());
  }

  done_saving_ids = sort_work.sides.empty() && sort_work.interiors.empty();
  system.comm().max(done_saving_ids);

  this->send_and_insert_dof_values(ids_to_push, action);

  {
    std::vector<node_projection> edges(sort_work.edges.begin(), sort_work.edges.end());
    ProjectEdges project_edges(*this);
    Threads::parallel_reduce (node_range(&edges), project_edges);
    ids_to_push.insert(project_edges.new_ids_to_push.begin(),
                       project_edges.new_ids_to_push.end());
    ids_to_save.insert(project_edges.new_ids_to_save.begin(),
                       project_edges.new_ids_to_save.end());
  }

  done_saving_ids = sort_work.interiors.empty();
  system.comm().max(done_saving_ids);

  this->send_and_insert_dof_values(ids_to_push, action);

  {
    std::vector<node_projection> sides(sort_work.sides.begin(), sort_work.sides.end());
    ProjectSides project_sides(*this);
    Threads::parallel_reduce (node_range(&sides), project_sides);
    ids_to_push.insert(project_sides.new_ids_to_push.begin(),
                       project_sides.new_ids_to_push.end());
    ids_to_save.insert(project_sides.new_ids_to_save.begin(),
                       project_sides.new_ids_to_save.end());
  }

  done_saving_ids = true;
  this->send_and_insert_dof_values(ids_to_push, action);

  // No ids to save or push this time, but we still use a reduce since
  // nominally ProjectInteriors still has non-const operator()
  ProjectInteriors project_interiors(*this);
  Threads::parallel_reduce (interior_range(&sort_work.interiors),
                            project_interiors);
}


template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::SubFunctor::SubFunctor
  (GenericProjector & p) :
  projector(p), action(p.master_action), f(p.master_f),
  context(p.system), conts(p.system.n_vars()), field_types(p.system.n_vars()), system(p.system)
{
  // Loop over all the variables we've been requested to project, to
  // pre-request
  for (const auto & var : this->projector.variables)
    {
      // FIXME: Need to generalize this to vector-valued elements. [PB]
      FEAbstract * fe = nullptr;
      FEAbstract * side_fe = nullptr;
      FEAbstract * edge_fe = nullptr;

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

          fe->request_phi();
          if (dim > 1)
            {
              side_fe->get_JxW();
              side_fe->get_xyz();
              side_fe->request_phi();
            }
          if (dim > 2)
            {
              edge_fe->get_JxW();
              edge_fe->get_xyz();
              edge_fe->request_phi();
            }

          const FEContinuity cont = fe->get_continuity();
          this->conts[var] = cont;
          if (cont == C_ONE)
            {
              fe->request_dphi();
              if (dim > 1)
                side_fe->request_dphi();
              if (dim > 2)
                edge_fe->request_dphi();
            }

          this->field_types[var] = FEInterface::field_type(fe->get_fe_type());
        }
    }

  // Now initialize any user requested shape functions, xyz vals, etc
  f.init_context(context);
}


template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::SubProjector::SubProjector
  (GenericProjector & p) :
  SubFunctor(p)
{
  if (p.master_g)
    g = libmesh_make_unique<GFunctor>(*p.master_g);

#ifndef NDEBUG
  // Our C1 elements need gradient information
  for (const auto & var : this->projector.variables)
    if (this->conts[var] == C_ONE)
      libmesh_assert(g);
#endif

  if (g)
    g->init_context(context);
}

template <typename FFunctor, typename GFunctor, typename FValue, typename ProjectionAction>
template <typename InsertInput,
          typename std::enable_if<
              !std::is_same<typename ProjectionAction::InsertInput, InsertInput>::value,
              int>::type>
void
GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::SubFunctor::insert_id(
    dof_id_type, const InsertInput & , processor_id_type)
{
  libmesh_error_msg("ID insertion should only occur when the provided input matches that "
                    "expected by the ProjectionAction");
}

template <typename FFunctor, typename GFunctor, typename FValue, typename ProjectionAction>
template <typename InsertInput,
          typename std::enable_if<
              std::is_same<typename ProjectionAction::InsertInput, InsertInput>::value,
              int>::type>
void
GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::SubFunctor::insert_id(
    dof_id_type id, const InsertInput & val, processor_id_type pid)
{
  auto iter = new_ids_to_push.find(id);
  if (iter == new_ids_to_push.end())
    action.insert(id, val);
  else
    {
      libmesh_assert(pid != DofObject::invalid_processor_id);
      iter->second = std::make_pair(val, pid);
    }
  if (!this->projector.done_saving_ids)
    {
      libmesh_assert(!new_ids_to_save.count(id));
      new_ids_to_save[id] = val;
    }
}

template <typename FFunctor, typename GFunctor, typename FValue, typename ProjectionAction>
template <typename InsertInput,
          typename std::enable_if<
              !std::is_same<typename ProjectionAction::InsertInput, InsertInput>::value,
              int>::type>
void
GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::SubFunctor::insert_ids(
    const std::vector<dof_id_type> &,
    const std::vector<InsertInput> &,
    processor_id_type)
{
  libmesh_error_msg("ID insertion should only occur when the provided input matches that "
                    "expected by the ProjectionAction");
}

template <typename FFunctor, typename GFunctor, typename FValue, typename ProjectionAction>
template <typename InsertInput,
          typename std::enable_if<
              std::is_same<typename ProjectionAction::InsertInput, InsertInput>::value,
              int>::type>
void
GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::SubFunctor::insert_ids(
    const std::vector<dof_id_type> & ids,
    const std::vector<InsertInput> & vals,
    processor_id_type pid)
{
  libmesh_assert_equal_to(ids.size(), vals.size());
  for (auto i : index_range(ids))
    {
      const dof_id_type id = ids[i];
      const InsertInput & val = vals[i];

      auto iter = new_ids_to_push.find(id);
      if (iter == new_ids_to_push.end())
        action.insert(id, val);
      else
        {
          libmesh_assert(pid != DofObject::invalid_processor_id);
          iter->second = std::make_pair(val, pid);
        }
      if (!this->projector.done_saving_ids)
        {
          libmesh_assert(!new_ids_to_save.count(id));
          new_ids_to_save[id] = val;
        }
    }
}

template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
void GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::SubFunctor::join
  (const GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::SubFunctor & other)
{
  new_ids_to_push.insert(other.new_ids_to_push.begin(), other.new_ids_to_push.end());
  new_ids_to_save.insert(other.new_ids_to_save.begin(), other.new_ids_to_save.end());
}


template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
void GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::SortAndCopy::operator()
  (const ConstElemRange & range)
{
  // Look at all the elements in the range.  Directly copy values from
  // unchanged elements.  For other elements, determine sets of
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

  // If we copy DoFs on a node, keep track of it so we don't bother
  // with any redundant interpolations or projections later.
  //
  // We still need to keep track of *which* variables get copied,
  // since we may be able to copy elements which lack some
  // subdomain-restricted variables.
  //
  // For extra_hanging_dofs FE types, we'll keep track of which
  // variables were copied from vertices, and which from edges/sides,
  // because those types need copies made from *both* on hanging
  // nodes.
  std::unordered_map<const Node *, std::pair<var_set, var_set>> copied_nodes;

  const unsigned int sys_num = system.number();

  // At hanging nodes for variables with extra hanging dofs we'll need
  // to do projections *separately* from vertex elements and side/edge
  // elements.
  std::vector<unsigned short> extra_hanging_dofs;
  bool all_extra_hanging_dofs = true;
  for (auto v_num : this->projector.variables)
    {
      if (extra_hanging_dofs.size() <= v_num)
        extra_hanging_dofs.resize(v_num+1, false);
      extra_hanging_dofs[v_num] =
        FEInterface::extra_hanging_dofs(system.variable(v_num).type());

      if (!extra_hanging_dofs[v_num])
        all_extra_hanging_dofs = false;
    }

  for (const auto & elem : range)
    {
      // If we're doing AMR, we might be able to copy more DoFs than
      // we interpolate or project.
      bool copy_this_elem = false;

#ifdef LIBMESH_ENABLE_AMR
      // If we're projecting from an old grid
      if (f.is_grid_projection())
        {
          // If this element doesn't have an old_dof_object, but it
          // wasn't just refined or just coarsened into activity, then
          // it must be newly added, so the user is responsible for
          // setting the new dofs on it during a grid projection.
          if (!elem->old_dof_object &&
              elem->refinement_flag() != Elem::JUST_REFINED &&
              elem->refinement_flag() != Elem::JUST_COARSENED)
            continue;

          // If this is an unchanged element, just copy everything
          if ((elem->refinement_flag() != Elem::JUST_REFINED &&
              elem->refinement_flag() != Elem::JUST_COARSENED &&
              elem->p_refinement_flag() != Elem::JUST_REFINED &&
              elem->p_refinement_flag() != Elem::JUST_COARSENED))
            copy_this_elem = true;
          else
            {
              bool reinitted = false;

              // If this element has a low order monomial which has
              // merely been h refined, copy it.
              for (auto v_num : this->projector.variables)
                {
                  const Variable & var = system.variable(v_num);
                  if (!var.active_on_subdomain(elem->subdomain_id()))
                    continue;
                  FEType fe_type = var.type();

                  if (fe_type.family == MONOMIAL &&
                      fe_type.order == CONSTANT &&
                      elem->p_level() == 0 &&
                      elem->refinement_flag() != Elem::JUST_COARSENED &&
                      elem->p_refinement_flag() != Elem::JUST_COARSENED)
                    {
                      if (!reinitted)
                        {
                          reinitted = true;
                          context.pre_fe_reinit(system, elem);
                        }

                      std::vector<typename FFunctor::ValuePushType> Ue(1);
                      std::vector<dof_id_type> elem_dof_ids(1);

                      f.eval_old_dofs(*elem, fe_type, sys_num, v_num,
                                      elem_dof_ids, Ue);

                      action.insert(elem_dof_ids[0], Ue[0]);
                    }
                }
            }
        }
#endif // LIBMESH_ENABLE_AMR

      const int dim = elem->dim();

      const unsigned int n_vertices = elem->n_vertices();
      const unsigned int n_edges = elem->n_edges();
      const unsigned int n_nodes = elem->n_nodes();

      // In 1-D we already handle our sides as vertices
      const unsigned int n_sides = (dim > 1) * elem->n_sides();

      // What variables are supported on each kind of node on this elem?
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

      for (auto v_num : this->projector.variables)
        {
          const Variable & var = this->projector.system.variable(v_num);
          if (!var.active_on_subdomain(elem->subdomain_id()))
            continue;
          FEType fe_type = var.type();

          if (FEInterface::n_dofs_at_node(fe_type, elem, 0))
            vertex_vars.insert(vertex_vars.end(), v_num);

          // The first non-vertex node is always an edge node if those
          // exist.  All edge nodes have the same number of DoFs
          if (has_edge_nodes)
            if (FEInterface::n_dofs_at_node(fe_type, elem, n_vertices))
              edge_vars.insert(edge_vars.end(), v_num);

          if (has_side_nodes)
            {
              if (dim != 3)
                {
                  if (FEInterface::n_dofs_at_node(fe_type, elem, n_vertices))
                    side_vars.insert(side_vars.end(), v_num);
                }
              else
                // In 3D, not all face nodes always have the same number of
                // DoFs!  We'll loop over all sides to be safe.
                for (unsigned int n = 0; n != n_nodes; ++n)
                  if (elem->is_face(n))
                    if (FEInterface::n_dofs_at_node(fe_type, elem, n))
                      {
                        side_vars.insert(side_vars.end(), v_num);
                        break;
                      }
            }

          if (FEInterface::n_dofs_per_elem(fe_type, elem) ||
              (has_interior_nodes &&
               FEInterface::n_dofs_at_node(fe_type, elem, n_nodes-1)))
            {
#ifdef LIBMESH_ENABLE_AMR
              // We may have already just copied constant monomials,
              // or we may be about to copy the whole element
              if ((f.is_grid_projection() &&
                   fe_type.family == MONOMIAL &&
                   fe_type.order == CONSTANT &&
                   elem->p_level() == 0 &&
                   elem->refinement_flag() != Elem::JUST_COARSENED &&
                   elem->p_refinement_flag() != Elem::JUST_COARSENED)
                  || copy_this_elem
                  )
                continue;
#endif // LIBMESH_ENABLE_AMR

              // We need to project any other variables
              if (interiors.empty() || interiors.back() != elem)
                interiors.push_back(elem);
            }
        }

      // We'll use a greedy algorithm in most cases: if another
      // element has already claimed some of our DoFs, we'll let it do
      // the work.
      auto erase_covered_vars = []
        (const Node * node, var_set & remaining, const ves_multimap & covered)
        {
          auto covered_range = covered.equal_range(node);
          for (const auto & v_ent : as_range(covered_range))
            for (const unsigned int var_covered :
                 std::get<2>(v_ent.second))
              remaining.erase(var_covered);
        };

      auto erase_nonhanging_vars = [&extra_hanging_dofs]
        (const Node * node, var_set & remaining, const ves_multimap & covered)
        {
          auto covered_range = covered.equal_range(node);
          for (const auto & v_ent : as_range(covered_range))
            for (const unsigned int var_covered :
                 std::get<2>(v_ent.second))
              if (!extra_hanging_dofs[var_covered])
                remaining.erase(var_covered);
        };

      auto erase_copied_vars = [&copied_nodes, &extra_hanging_dofs]
        (const Node * node, bool is_vertex, var_set & remaining)
        {
          auto copying_range = copied_nodes.equal_range(node);
          for (const auto & v_ent : as_range(copying_range))
            {
              for (const unsigned int var_covered :
                   v_ent.second.first)
                if (is_vertex || !extra_hanging_dofs[var_covered])
                  remaining.erase(var_covered);

              for (const unsigned int var_covered :
                   v_ent.second.second)
                if (!is_vertex || !extra_hanging_dofs[var_covered])
                  remaining.erase(var_covered);
            }
        };

      for (unsigned int v=0; v != n_vertices; ++v)
        {
          const Node * node = elem->node_ptr(v);

          auto remaining_vars = vertex_vars;

          erase_covered_vars(node, remaining_vars, vertices);

          if (remaining_vars.empty())
            continue;

          if (!all_extra_hanging_dofs)
            {
              erase_nonhanging_vars(node, remaining_vars, edges);
              if (remaining_vars.empty())
                continue;

              erase_nonhanging_vars(node, remaining_vars, sides);
              if (remaining_vars.empty())
                continue;
            }

          erase_copied_vars(node, true, remaining_vars);
          if (remaining_vars.empty())
            continue;

          if (copy_this_elem)
            {
              for (auto var : remaining_vars)
                {
                  std::vector<dof_id_type> node_dof_ids;
                  std::vector<typename FFunctor::ValuePushType> values;

                  f.eval_old_dofs(*elem, v, var, node_dof_ids, values);

                  insert_ids(node_dof_ids, values, node->processor_id());
                }
              copied_nodes[node].first.insert(remaining_vars.begin(),
                                              remaining_vars.end());
              this->find_dofs_to_send(*node, *elem, v, remaining_vars);
            }
          else
            vertices.emplace
              (node, std::make_tuple(elem, v, std::move(remaining_vars)));
        }

      if (has_edge_nodes)
        {
          for (unsigned int e=0; e != n_edges; ++e)
            {
              const Node * node = elem->node_ptr(n_vertices+e);

              auto remaining_vars = edge_vars;

              erase_covered_vars(node, remaining_vars, edges);
              if (remaining_vars.empty())
                continue;

              erase_covered_vars(node, remaining_vars, sides);
              if (remaining_vars.empty())
                continue;

              if (!all_extra_hanging_dofs)
                {
                  erase_nonhanging_vars(node, remaining_vars, vertices);
                  if (remaining_vars.empty())
                    continue;
                }

              erase_copied_vars(node, false, remaining_vars);
              if (remaining_vars.empty())
                continue;

              if (copy_this_elem)
                {
                  for (auto var : remaining_vars)
                    {
                      std::vector<dof_id_type> edge_dof_ids;
                      std::vector<typename FFunctor::ValuePushType> values;

                      f.eval_old_dofs(*elem, n_vertices+e, var, edge_dof_ids, values);

                      insert_ids(edge_dof_ids, values, node->processor_id());
                    }
                  copied_nodes[node].second.insert(remaining_vars.begin(),
                                                   remaining_vars.end());
                  this->find_dofs_to_send(*node, *elem, n_vertices+e, remaining_vars);
                }
              else
                edges.emplace
                  (node, std::make_tuple(elem, e, std::move(remaining_vars)));
            }
        }

      if (has_side_nodes)
        {
          for (unsigned int side=0; side != n_sides; ++side)
            {
              const Node * node = nullptr;
              unsigned short node_num = n_vertices+(dim>2)*n_edges+side;
              if (dim != 3)
                node = elem->node_ptr(node_num);
              else
                {
                  // In 3D only some sides may have nodes
                  for (unsigned int n = 0; n != n_nodes; ++n)
                    {
                      if (!elem->is_face(n))
                        continue;

                      if (elem->is_node_on_side(n, side))
                        {
                          node_num = n;
                          node = elem->node_ptr(node_num);
                          break;
                        }
                    }
                }

              if (!node)
                continue;

              auto remaining_vars = side_vars;

              erase_covered_vars(node, remaining_vars, edges);
              if (remaining_vars.empty())
                continue;

              erase_covered_vars(node, remaining_vars, sides);
              if (remaining_vars.empty())
                continue;

              if (!all_extra_hanging_dofs)
                {
                  erase_nonhanging_vars(node, remaining_vars, vertices);
                  if (remaining_vars.empty())
                    continue;
                }

              erase_copied_vars(node, false, remaining_vars);
              if (remaining_vars.empty())
                continue;

              if (copy_this_elem)
                {
                  for (auto var : remaining_vars)
                    {
                      std::vector<dof_id_type> side_dof_ids;
                      std::vector<typename FFunctor::ValuePushType> values;

                      f.eval_old_dofs(*elem, node_num, var, side_dof_ids, values);

                      insert_ids(side_dof_ids, values, node->processor_id());
                    }
                  copied_nodes[node].second.insert(remaining_vars.begin(),
                                                   remaining_vars.end());
                  this->find_dofs_to_send(*node, *elem, node_num, remaining_vars);
                }
              else
                sides.emplace
                  (node, std::make_tuple(elem, side, std::move(remaining_vars)));
            }
        }

      // Elements with elemental dofs might need those copied too.
      if (copy_this_elem)
        {
          for (auto v_num : this->projector.variables)
            {
              const Variable & var = system.variable(v_num);
              if (!var.active_on_subdomain(elem->subdomain_id()))
                continue;
              FEType fe_type = var.type();

              std::vector<typename FFunctor::ValuePushType> Ue;
              std::vector<dof_id_type> elem_dof_ids;
              f.eval_old_dofs(*elem, fe_type, sys_num, v_num,
                              elem_dof_ids, Ue);
              action.insert(elem_dof_ids, Ue);

              if (has_interior_nodes)
                {
                  std::vector<typename FFunctor::ValuePushType> Un;
                  std::vector<dof_id_type> node_dof_ids;

                  f.eval_old_dofs(*elem, n_nodes-1, v_num, node_dof_ids, Un);
                  action.insert(node_dof_ids, Un);
                }
            }
        }
    }
}


template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
void GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::SortAndCopy::join
  (const GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::SortAndCopy & other)
{
  auto merge_multimaps = [](ves_multimap & self, const ves_multimap & other_mm)
    {
      for (const auto & pair : other_mm)
        {
          const Node * node = pair.first;
          auto remaining_vars = std::get<2>(pair.second);

          auto my_range = self.equal_range(node);
          for (const auto & v_ent : as_range(my_range))
            for (const unsigned int var_covered :
                 std::get<2>(v_ent.second))
              remaining_vars.erase(var_covered);

          if (!remaining_vars.empty())
            self.emplace
              (node, std::make_tuple(std::get<0>(pair.second),
                                     std::get<1>(pair.second),
                                     std::move(remaining_vars)));
        }
    };

  merge_multimaps(vertices, other.vertices);
  merge_multimaps(edges, other.edges);
  merge_multimaps(sides, other.sides);

  std::sort(interiors.begin(), interiors.end());
  std::vector<const Elem *> other_interiors = other.interiors;
  std::sort(other_interiors.begin(), other_interiors.end());
  std::vector<const Elem *> merged_interiors;
  std::set_union(interiors.begin(), interiors.end(),
                 other_interiors.begin(), other_interiors.end(),
                 std::back_inserter(merged_interiors));
  interiors.swap(merged_interiors);

  SubFunctor::join(other);
}

namespace
{
template <typename Output, typename Input>
typename std::enable_if<ScalarTraits<Input>::value, Output>::type
raw_value(const Input & input, unsigned int /*index*/)
{
  return input;
}

template <typename Output, template <typename> class WrapperClass, typename T>
typename std::enable_if<ScalarTraits<T>::value &&
                            TensorTools::MathWrapperTraits<WrapperClass<T>>::value,
                        Output>::type
raw_value(const WrapperClass<T> & input, unsigned int index)
{
  return input.slice(index);
}

template <typename T>
typename T::value_type
grad_component(const T &, unsigned int);

template <typename T>
typename VectorValue<T>::value_type
grad_component(const VectorValue<T> & grad, unsigned int component)
{
  return grad(component);
}

template <typename T>
typename TensorValue<T>::value_type
grad_component(const TensorValue<T> & grad, unsigned int component)
{
  libmesh_error_msg("This function should only be called for Hermites. "
                    "I don't know how you got here");
  return grad(component, component);
}


}

template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
void GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::ProjectVertices::operator()
  (const node_range & range)
{
  LOG_SCOPE ("project_vertices","GenericProjector");

  const unsigned int sys_num = system.number();

  // Variables with extra hanging dofs can't safely use eval_at_node
  // in as many places as variables without can.
  std::vector<unsigned short> extra_hanging_dofs;
  for (auto v_num : this->projector.variables)
    {
      if (extra_hanging_dofs.size() <= v_num)
        extra_hanging_dofs.resize(v_num+1, false);
      extra_hanging_dofs[v_num] =
        FEInterface::extra_hanging_dofs(system.variable(v_num).type());
    }

  for (const auto & v_pair : range)
    {
      const Node & vertex = *v_pair.first;
      const Elem & elem = *std::get<0>(v_pair.second);
      const unsigned int n = std::get<1>(v_pair.second);
      const var_set & vertex_vars = std::get<2>(v_pair.second);

      context.pre_fe_reinit(system, &elem);

      this->find_dofs_to_send(vertex, elem, n, vertex_vars);

      // Look at all the variables we're supposed to interpolate from
      // this element on this vertex
      for (const auto & var : vertex_vars)
        {
          const Variable & variable = system.variable(var);
          const FEType & base_fe_type = variable.type();
          const unsigned int var_component =
            system.variable_scalar_number(var, 0);

          if (base_fe_type.family == SCALAR)
            continue;

          const FEContinuity cont = this->conts[var];
          const FEFieldType field_type = this->field_types[var];

          if (cont == DISCONTINUOUS)
            {
              libmesh_assert_equal_to(vertex.n_comp(sys_num, var), 0);
            }
          else if (cont == C_ZERO)
            {
              libmesh_assert(vertex.n_comp(sys_num, var));
              const FValue val = f.eval_at_node
                (context, var_component, /*dim=*/ 0, // Don't care w/C0
                 vertex, extra_hanging_dofs[var], system.time);

              if (field_type == TYPE_VECTOR)
              {
                // We will have a number of nodal value DoFs equal to the elem dim
                for (unsigned int i = 0; i < elem.dim(); ++i)
                {
                  const dof_id_type id = vertex.dof_number(sys_num, var, i);

                  // Need this conversion so that this method
                  // will compile for TYPE_SCALAR instantiations
                  const auto insert_val =
                    raw_value<typename ProjectionAction::InsertInput>(val, i);

                  insert_id(id, insert_val, vertex.processor_id());
                }
              }
              else
              {
                // C_ZERO elements have a single nodal value DoF at vertices
                const dof_id_type id = vertex.dof_number(sys_num, var, 0);
                insert_id(id, val, vertex.processor_id());
              }
            }
          else if (cont == C_ONE)
            {
              libmesh_assert(vertex.n_comp(sys_num, var));
              const dof_id_type first_id = vertex.dof_number(sys_num, var, 0);

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
              for (const auto e_id : (*this->projector.nodes_to_elem)[vertex.id()])
                {
                  const Elem & e = system.get_mesh().elem_ref(e_id);
                  libmesh_assert_equal_to(dim, e.dim());
                }
#endif
#ifdef LIBMESH_ENABLE_AMR
              bool is_old_vertex = true;
              if (elem.refinement_flag() == Elem::JUST_REFINED)
                {
                  const int i_am_child =
                    elem.parent()->which_child_am_i(&elem);
                  is_old_vertex =
                    elem.parent()->is_vertex_on_parent(i_am_child, n);
                }
#else
              const bool is_old_vertex = false;
#endif

              // The hermite element vertex shape functions are weird
              if (base_fe_type.family == HERMITE)
                {
                  const FValue val =
                    f.eval_at_node(context,
                                   var_component,
                                   dim,
                                   vertex,
                                   extra_hanging_dofs[var],
                                   system.time);
                  insert_id(first_id, val, vertex.processor_id());

                  typename GFunctor::FunctorValue grad =
                    is_old_vertex ?
                    g->eval_at_node(context,
                                    var_component,
                                    dim,
                                    vertex,
                                    extra_hanging_dofs[var],
                                    system.time) :
                    g->eval_at_point(context,
                                     var_component,
                                     vertex,
                                     system.time,
                                     false);
                  // x derivative. Use slice because grad may be a tensor type
                  insert_id(first_id+1, grad.slice(0),
                            vertex.processor_id());
#if LIBMESH_DIM > 1
                  if (dim > 1 && is_old_vertex && f.is_grid_projection())
                    {
                      for (int i = 1; i < dim; ++i)
                        insert_id(first_id+i+1, grad.slice(i),
                                  vertex.processor_id());

                      // We can directly copy everything else too
                      std::vector<FValue> derivs;
                      f.eval_mixed_derivatives
                        (context, var_component, dim, vertex, derivs);
                      for (auto i : index_range(derivs))
                        insert_id(first_id+dim+1+i, derivs[i],
                                  vertex.processor_id());
                    }
                  else if (dim > 1)
                    {
                      // We'll finite difference mixed derivatives.
                      // This delta_x used to be TOLERANCE*hmin, but
                      // the factor of 10 improved the accuracy in
                      // some unit test projections
                      Real delta_x = TOLERANCE * 10 * elem.hmin();

                      Point nxminus = elem.point(n),
                            nxplus = elem.point(n);
                      nxminus(0) -= delta_x;
                      nxplus(0) += delta_x;
                      typename GFunctor::FunctorValue gxminus =
                        g->eval_at_point(context,
                                         var_component,
                                         nxminus,
                                         system.time,
                                         true);
                      typename GFunctor::FunctorValue gxplus =
                        g->eval_at_point(context,
                                         var_component,
                                         nxplus,
                                         system.time,
                                         true);
                      // y derivative
                      insert_id(first_id+2, grad.slice(1),
                                vertex.processor_id());
                      // xy derivative
                      insert_id(first_id+3,
                        (grad_component(gxplus, 1) - grad_component(gxminus, 1)) / 2. / delta_x,
                        vertex.processor_id());

#if LIBMESH_DIM > 2
                      if (dim > 2)
                        {
                          // z derivative
                          insert_id(first_id+4, grad.slice(2),
                                    vertex.processor_id());
                          // xz derivative
                          insert_id(first_id+5,
                            (grad_component(gxplus, 2) - grad_component(gxminus, 2)) / 2. / delta_x,
                            vertex.processor_id());

                          // We need new points for yz
                          Point nyminus = elem.point(n),
                            nyplus = elem.point(n);
                          nyminus(1) -= delta_x;
                          nyplus(1) += delta_x;
                          typename GFunctor::FunctorValue gyminus =
                            g->eval_at_point(context,
                                             var_component,
                                             nyminus,
                                             system.time,
                                             true);
                          typename GFunctor::FunctorValue gyplus =
                            g->eval_at_point(context,
                                             var_component,
                                             nyplus,
                                             system.time,
                                             true);
                          // yz derivative
                          insert_id(first_id+6,
                            (grad_component(gyplus, 2) - grad_component(gyminus, 2)) / 2. / delta_x,
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
                          typename GFunctor::FunctorValue gxmym =
                            g->eval_at_point(context,
                                             var_component,
                                             nxmym,
                                             system.time,
                                             true);
                          typename GFunctor::FunctorValue gxmyp =
                            g->eval_at_point(context,
                                             var_component,
                                             nxmyp,
                                             system.time,
                                             true);
                          typename GFunctor::FunctorValue gxpym =
                            g->eval_at_point(context,
                                             var_component,
                                             nxpym,
                                             system.time,
                                             true);
                          typename GFunctor::FunctorValue gxpyp =
                            g->eval_at_point(context,
                                             var_component,
                                             nxpyp,
                                             system.time,
                                             true);
                          FValue gxzplus = (grad_component(gxpyp, 2) - grad_component(gxmyp, 2))
                            / 2. / delta_x;
                          FValue gxzminus = (grad_component(gxpym, 2) - grad_component(gxmym, 2))
                            / 2. / delta_x;
                          // xyz derivative
                          insert_id(first_id+7,
                            (gxzplus - gxzminus) / 2. / delta_x,
                            vertex.processor_id());
                        }
#endif // LIBMESH_DIM > 2
                    }
#endif // LIBMESH_DIM > 1
                }
              else
                {
                  // Currently other C_ONE elements have a single nodal
                  // value shape function and nodal gradient component
                  // shape functions
                  libmesh_assert_equal_to
                    (FEInterface::n_dofs_at_node
                      (base_fe_type, &elem,
                       elem.get_node_index(&vertex)),
                    (unsigned int)(1 + dim));
                  const FValue val =
                    f.eval_at_node(context, var_component, dim,
                                   vertex, extra_hanging_dofs[var],
                                   system.time);
                  insert_id(first_id, val, vertex.processor_id());
                  typename GFunctor::FunctorValue grad =
                    is_old_vertex ?
                    g->eval_at_node(context, var_component, dim,
                                    vertex, extra_hanging_dofs[var],
                                    system.time) :
                    g->eval_at_point(context, var_component, vertex,
                                     system.time, false);
                  for (int i=0; i!= dim; ++i)
                    insert_id(first_id + i + 1, grad.slice(i),
                              vertex.processor_id());
                }
            }
          else
            libmesh_error_msg("Unknown continuity " << cont);
        }
    }
}


template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
void GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::ProjectEdges::operator()
  (const node_range & range)
{
  LOG_SCOPE ("project_edges","GenericProjector");

  const unsigned int sys_num = system.number();

  for (const auto & e_pair : range)
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
      const var_set & edge_vars = std::get<2>(e_pair.second);

      const unsigned short edge_num = std::get<1>(e_pair.second);
      const unsigned short node_num = elem.n_vertices() + edge_num;
      context.edge = edge_num;
      context.pre_fe_reinit(system, &elem);

      this->find_dofs_to_send(edge_node, elem, node_num, edge_vars);

      // Look at all the variables we're supposed to interpolate from
      // this element on this edge
      for (const auto & var : edge_vars)
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
          if (fe_type.family == LAGRANGE || fe_type.family == LAGRANGE_VEC)
            {
              if (fe_type.order > 1)
                {
                  const processor_id_type pid =
                    edge_node.processor_id();
                  FValue fval = f.eval_at_point
                    (context, var_component, edge_node, system.time,
                     false);
                  if (fe_type.family == LAGRANGE_VEC)
                  {
                    // We will have a number of nodal value DoFs equal to the elem dim
                    for (unsigned int i = 0; i < elem.dim(); ++i)
                    {
                      const dof_id_type dof_id =
                        edge_node.dof_number(sys_num, var, i);

                      // Need this conversion so that this method
                      // will compile for TYPE_SCALAR instantiations
                      const auto insert_val =
                        raw_value<typename ProjectionAction::InsertInput>(fval, i);

                      insert_id(dof_id, insert_val, pid);
                    }
                  }
                  else // We are LAGRANGE
                  {
                    const dof_id_type dof_id =
                      edge_node.dof_number(sys_num, var, 0);
                    insert_id(dof_id, fval, pid);
                  }
                }
              continue;
            }

#ifdef LIBMESH_ENABLE_AMR
          // If this is a low order monomial element which has merely
          // been h refined then we already copied all its dofs
          if (fe_type.family == MONOMIAL &&
              fe_type.order == CONSTANT &&
              elem.refinement_flag() != Elem::JUST_COARSENED &&
              elem.p_refinement_flag() != Elem::JUST_COARSENED)
            continue;
#endif // LIBMESH_ENABLE_AMR

          // FIXME: Need to generalize this to vector-valued elements. [PB]
          FEGenericBase<typename FFunctor::RealType> * fe = nullptr;
          context.get_element_fe( var, fe, dim );
          FEGenericBase<typename FFunctor::RealType> * edge_fe = nullptr;
          context.get_edge_fe( var, edge_fe );

          // If we're JUST_COARSENED we'll need a custom
          // evaluation, not just the standard edge FE
          const FEGenericBase<typename FFunctor::RealType> & proj_fe =
#ifdef LIBMESH_ENABLE_AMR
            (elem.refinement_flag() == Elem::JUST_COARSENED) ?
            *fe :
#endif
            *edge_fe;

#ifdef LIBMESH_ENABLE_AMR
          if (elem.refinement_flag() == Elem::JUST_COARSENED)
            {
              std::vector<Point> fine_points;

              std::unique_ptr<FEGenericBase<typename FFunctor::RealType>> fine_fe
                (FEGenericBase<typename FFunctor::RealType>::build (dim, base_fe_type));

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
              FEMap::inverse_map (dim, &elem, fine_points, fine_qp);

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
            (dof_indices, edge_dofs, var_component,
             &edge_node, proj_fe);
        }
    }
}


template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
void GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::ProjectSides::operator()
  (const node_range & range)
{
  LOG_SCOPE ("project_sides","GenericProjector");

  const unsigned int sys_num = system.number();

  for (const auto & s_pair : range)
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
      const var_set & side_vars = std::get<2>(s_pair.second);

      const unsigned int side_num = std::get<1>(s_pair.second);
      unsigned short node_num = elem.n_vertices()+side_num;
      // In 3D only some sides may have nodes
      if (dim == 3)
        for (auto n : make_range(elem.n_nodes()))
          {
            if (!elem.is_face(n))
              continue;

            if (elem.is_node_on_side(n, side_num))
              {
                node_num = n;
                break;
              }
          }

      context.side = side_num;
      context.pre_fe_reinit(system, &elem);

      this->find_dofs_to_send(side_node, elem, node_num, side_vars);

      // Look at all the variables we're supposed to interpolate from
      // this element on this side
      for (const auto & var : side_vars)
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
          if (fe_type.family == LAGRANGE || fe_type.family == LAGRANGE_VEC)
            {
              if (fe_type.order > 1)
                {
                  const processor_id_type pid =
                    side_node.processor_id();
                  FValue fval = f.eval_at_point
                    (context, var_component, side_node, system.time,
                     false);

                  if (fe_type.family == LAGRANGE_VEC)
                  {
                    // We will have a number of nodal value DoFs equal to the elem dim
                    for (unsigned int i = 0; i < elem.dim(); ++i)
                    {
                      const dof_id_type dof_id = side_node.dof_number(sys_num, var, i);

                      // Need this conversion so that this method
                      // will compile for TYPE_SCALAR instantiations
                      const auto insert_val =
                        raw_value<typename ProjectionAction::InsertInput>(fval, i);

                      insert_id(dof_id, insert_val, pid);
                    }
                  }
                  else // We are LAGRANGE
                  {
                    const dof_id_type dof_id =
                      side_node.dof_number(sys_num, var, 0);
                    insert_id(dof_id, fval, pid);
                  }
                }
              continue;
            }

#ifdef LIBMESH_ENABLE_AMR
          // If this is a low order monomial element which has merely
          // been h refined then we already copied all its dofs
          if (fe_type.family == MONOMIAL &&
              fe_type.order == CONSTANT &&
              elem.refinement_flag() != Elem::JUST_COARSENED &&
              elem.p_refinement_flag() != Elem::JUST_COARSENED)
            continue;
#endif // LIBMESH_ENABLE_AMR

          // FIXME: Need to generalize this to vector-valued elements. [PB]
          FEGenericBase<typename FFunctor::RealType> * fe = nullptr;
          context.get_element_fe( var, fe, dim );
          FEGenericBase<typename FFunctor::RealType> * side_fe = nullptr;
          context.get_side_fe( var, side_fe );

          // If we're JUST_COARSENED we'll need a custom
          // evaluation, not just the standard side FE
          const FEGenericBase<typename FFunctor::RealType> & proj_fe =
#ifdef LIBMESH_ENABLE_AMR
            (elem.refinement_flag() == Elem::JUST_COARSENED) ?
            *fe :
#endif
            *side_fe;

#ifdef LIBMESH_ENABLE_AMR
          if (elem.refinement_flag() == Elem::JUST_COARSENED)
            {
              std::vector<Point> fine_points;

              std::unique_ptr<FEGenericBase<typename FFunctor::RealType>> fine_fe
                (FEGenericBase<typename FFunctor::RealType>::build (dim, base_fe_type));

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
              FEMap::inverse_map (dim, &elem, fine_points, fine_qp);

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
            (dof_indices, side_dofs, var_component,
             &side_node, proj_fe);
        }
    }
}


template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
void GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::ProjectInteriors::operator()
  (const interior_range & range)
{
  LOG_SCOPE ("project_interiors","GenericProjector");

  const unsigned int sys_num = system.number();

  // Iterate over all dof-bearing element interiors in the range
  for (const auto & elem : range)
    {
      unsigned char dim = cast_int<unsigned char>(elem->dim());

      context.pre_fe_reinit(system, elem);

      // Loop over all the variables we've been requested to project, to
      // do the projection
      for (const auto & var : this->projector.variables)
        {
          const Variable & variable = system.variable(var);

          if (!variable.active_on_subdomain(elem->subdomain_id()))
            continue;

          const FEType & base_fe_type = variable.type();

          if (base_fe_type.family == SCALAR)
            continue;

          FEGenericBase<typename FFunctor::RealType> * fe = nullptr;
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
          if (fe_type.family == LAGRANGE || fe_type.family == LAGRANGE_VEC)
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

                      FValue fval = f.eval_at_point
                        (context, var_component, interior_node,
                         system.time, false);
                      const processor_id_type pid =
                        interior_node.processor_id();

                      if (fe_type.family == LAGRANGE_VEC)
                      {
                        // We will have a number of nodal value DoFs equal to the elem dim
                        for (unsigned int i = 0; i < elem->dim(); ++i)
                        {
                          const dof_id_type dof_id = interior_node.dof_number(sys_num, var, i);

                          // Need this conversion so that this method
                          // will compile for TYPE_SCALAR instantiations
                          const auto insert_val =
                            raw_value<typename ProjectionAction::InsertInput>(fval, i);

                          insert_id(dof_id, insert_val, pid);
                        }
                      }
                      else // We are LAGRANGE
                      {
                        const dof_id_type dof_id =
                          interior_node.dof_number(sys_num, var, 0);
                        insert_id(dof_id, fval, pid);
                      }
                    }
                }
              continue;
            }

#ifdef LIBMESH_ENABLE_AMR
          if (elem->refinement_flag() == Elem::JUST_COARSENED)
            {
              std::vector<Point> fine_points;

              std::unique_ptr<FEGenericBase<typename FFunctor::RealType>> fine_fe
                (FEGenericBase<typename FFunctor::RealType>::build (dim, base_fe_type));

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
              FEMap::inverse_map (dim, elem, fine_points, fine_qp);

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
            (dof_indices, all_dofs, var_component,
             nullptr, *fe);
        } // end variables loop
    } // end elements loop
}



template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
void
GenericProjector<FFunctor, GFunctor, FValue,
                 ProjectionAction>::SubFunctor::find_dofs_to_send
  (const Node & node, const Elem & elem, unsigned short node_num, const var_set & vars)
{
  libmesh_assert (&node == elem.node_ptr(node_num));

  // Any ghosted node in our range might have an owner who needs our
  // data
  const processor_id_type owner = node.processor_id();
  if (owner != system.processor_id())
    {
      const MeshBase & mesh = system.get_mesh();
      const DofMap & dof_map = system.get_dof_map();

      // But let's check and see if we can be certain the owner can
      // compute any or all of its own dof coefficients on that node.
      std::vector<dof_id_type> node_dof_ids, patch_dof_ids;
      for (const auto & var : vars)
        {
          const Variable & variable = system.variable(var);

          if (!variable.active_on_subdomain(elem.subdomain_id()))
            continue;

          dof_map.dof_indices(elem, node_num, node_dof_ids, var);
        }
      libmesh_assert(std::is_sorted(node_dof_ids.begin(),
                                    node_dof_ids.end()));
      const std::vector<dof_id_type> & patch =
        (*this->projector.nodes_to_elem)[node.id()];
      for (const auto & elem_id : patch)
        {
          const Elem & patch_elem = mesh.elem_ref(elem_id);
          if (!patch_elem.active() || owner != patch_elem.processor_id())
            continue;
          dof_map.dof_indices(&patch_elem, patch_dof_ids);
          std::sort(patch_dof_ids.begin(), patch_dof_ids.end());

          // Remove any node_dof_ids that we see can be calculated on
          // this element
          std::vector<dof_id_type> diff_ids(node_dof_ids.size());
          auto it = std::set_difference(node_dof_ids.begin(), node_dof_ids.end(),
                                        patch_dof_ids.begin(), patch_dof_ids.end(), diff_ids.begin());
          diff_ids.resize(it-diff_ids.begin());
          node_dof_ids.swap(diff_ids);
          if (node_dof_ids.empty())
            break;
        }

      // Give ids_to_push default invalid pid: not yet computed
      for (auto id : node_dof_ids)
        new_ids_to_push[id].second = DofObject::invalid_processor_id;
    }
}



template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
template <typename Value>
void
GenericProjector<FFunctor, GFunctor, FValue,
                 ProjectionAction>::send_and_insert_dof_values
  (std::unordered_map<dof_id_type, std::pair<Value, processor_id_type>> & ids_to_push,
   ProjectionAction & action) const
{
  // See if we calculated any ids that need to be pushed; get them
  // ready to push.
  std::unordered_map<processor_id_type, std::vector<dof_id_type>>
    pushed_dof_ids, received_dof_ids;
  std::unordered_map<processor_id_type, std::vector<typename TypeToSend<Value>::type>>
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
              Value val;
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
     const std::vector<typename TypeToSend<Value>::type> & data)
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
              Value val;
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
void
GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::SubProjector::construct_projection
  (const std::vector<dof_id_type> & dof_indices_var,
   const std::vector<unsigned int> & involved_dofs,
   unsigned int var_component,
   const Node * node,
   const FEGenericBase<typename FFunctor::RealType> & fe)
{
  const auto & JxW = fe.get_JxW();
  const auto & phi = fe.get_phi();
  const std::vector<std::vector<typename FEGenericBase<typename FFunctor::RealType>::OutputGradient>> * dphi = nullptr;
  const std::vector<Point> & xyz_values = fe.get_xyz();
  const FEContinuity cont = fe.get_continuity();
  const std::unordered_map<dof_id_type, typename FFunctor::ValuePushType> & ids_to_save =
    this->projector.ids_to_save;

  if (cont == C_ONE)
    dphi = &(fe.get_dphi());

  const unsigned int n_involved_dofs =
    cast_int<unsigned int>(involved_dofs.size());

  std::vector<dof_id_type> free_dof_ids;
  DenseVector<typename FFunctor::ValuePushType> Uinvolved(n_involved_dofs);
  std::vector<char> dof_is_fixed(n_involved_dofs, false); // bools

  for (auto i : make_range(n_involved_dofs))
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
  DenseVector<typename FFunctor::ValuePushType> Fe(free_dofs);
  // The new degree of freedom coefficients to solve for
  DenseVector<typename FFunctor::ValuePushType> Ufree(free_dofs);

  const unsigned int n_qp =
    cast_int<unsigned int>(xyz_values.size());

  // Loop over the quadrature points
  for (unsigned int qp=0; qp<n_qp; qp++)
    {
      // solution at the quadrature point
      FValue fineval = f.eval_at_point(context,
                                       var_component,
                                       xyz_values[qp],
                                       system.time,
                                       false);
      // solution grad at the quadrature point
      typename GFunctor::FunctorValue finegrad;
      if (cont == C_ONE)
        finegrad = g->eval_at_point(context,
                                    var_component,
                                    xyz_values[qp],
                                    system.time,
                                    false);

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
                    Fe(freei) -= ( TensorTools::inner_product((*dphi)[i][qp],
                                                              (*dphi)[j][qp]) ) *
                      JxW[qp] * Uinvolved(sidej);
                  else
                    Ke(freei,freej) += ( TensorTools::inner_product((*dphi)[i][qp],
                                                                    (*dphi)[j][qp]) )
                      * JxW[qp];
                }
              if (!dof_is_fixed[sidej])
                freej++;
            }
          Fe(freei) += phi[i][qp] * fineval * JxW[qp];
          if (cont == C_ONE)
            Fe(freei) += (TensorTools::inner_product(finegrad,
                                                     (*dphi)[i][qp]) ) *
              JxW[qp];
          freei++;
        }
    }

  Ke.cholesky_solve(Fe, Ufree);

  // Transfer new edge solutions to element
  const processor_id_type pid = node ?
    node->processor_id() : DofObject::invalid_processor_id;
  insert_ids(free_dof_ids, Ufree.get_values(), pid);
}


} // namespace libMesh

#endif // GENERIC_PROJECTOR_H
