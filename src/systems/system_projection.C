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
#include <vector>

// Local includes
#include "libmesh/boundary_info.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/system.h"
#include "libmesh/threads.h"
#include "libmesh/wrapped_function.h"
#include "libmesh/wrapped_functor.h"

namespace libMesh
{

// ------------------------------------------------------------
// Helper class definitions

/**
 * This class implements the loops to other projection operations.
 * This may be executed in parallel on multiple threads.
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
 * This action class can be used with a GenericProjector to set
 * projection values (which must be of type Val) as coefficients of
 * the given NumericVector.
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
  UniquePtr<FEMFunctionBase<Output> > _f;
};


#ifdef LIBMESH_ENABLE_AMR
template <typename Output,
          void (FEMContext::*point_output) (unsigned int,
                                            const Point &,
                                            Output &,
                                            const Real) const>
class OldSolutionValue
{
public:
  OldSolutionValue(const libMesh::System & sys_in,
                   const NumericVector<Number> & old_sol) :
    last_elem(libmesh_nullptr),
    sys(sys_in),
    old_context(sys_in),
    old_solution(old_sol)
  {
    old_context.set_algebraic_type(FEMContext::OLD);
    old_context.set_custom_solution(&old_solution);
  }

  OldSolutionValue(const OldSolutionValue & in) :
    last_elem(libmesh_nullptr),
    sys(in.sys),
    old_context(sys),
    old_solution(in.old_solution)
  {
    old_context.set_algebraic_type(FEMContext::OLD);
    old_context.set_custom_solution(&old_solution);
  }

  static void get_shape_outputs(FEBase& fe);

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
    (old_context.*point_output)(i, p, n, out_of_elem_tol);
    return n;
  }

  bool is_grid_projection() { return true; }

  void eval_old_dofs (const FEMContext & c,
                      unsigned int var,
                      std::vector<Output> & values)
  {
    LOG_SCOPE ("eval_old_dofs()", "OldSolutionValue");

    this->check_old_context(c);

    const std::vector<dof_id_type> & old_dof_indices =
      old_context.get_dof_indices(var);

    libmesh_assert_equal_to (old_dof_indices.size(), values.size());

    old_solution.get(old_dof_indices, values);
  }

protected:
  void check_old_context (const FEMContext & c)
  {
    LOG_SCOPE ("check_old_context(c)", "OldSolutionValue");
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
    LOG_SCOPE ("check_old_context(c,p)", "OldSolutionValue");
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

            for (unsigned int c=0; c != elem.n_children(); ++c)
              if (elem.child_ptr(c)->close_to_point(p, master_tol))
                {
                  old_context.pre_fe_reinit(sys, elem.child_ptr(c));
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

            for (unsigned int c=0; c != elem.n_children(); ++c)
              if (elem.child_ptr(c)->close_to_point(p, master_tol))
                {
                  old_context.pre_fe_reinit(sys, elem.child_ptr(c));
                  break;
                }

            libmesh_assert
              (old_context.get_elem().close_to_point(p, master_tol));
          }
      }

    return true;
  }

private:
  const Elem * last_elem;
  const System & sys;
  FEMContext old_context;
  const NumericVector<Number> & old_solution;

  static const Real out_of_elem_tol;
};


template<>
inline void
OldSolutionValue<Number, &FEMContext::point_value>::get_shape_outputs(FEBase& fe)
{
  fe.get_phi();
}


template<>
inline void
OldSolutionValue<Gradient, &FEMContext::point_gradient>::get_shape_outputs(FEBase& fe)
{
  fe.get_dphi();
}


template<>
inline
Number
OldSolutionValue<Number, &FEMContext::point_value>::eval_at_node
  (const FEMContext & c,
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
OldSolutionValue<Gradient, &FEMContext::point_gradient>::eval_at_node
  (const FEMContext & c,
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
const Real OldSolutionValue<Number, &FEMContext::point_value>::out_of_elem_tol = 10*TOLERANCE;

template <>
const Real OldSolutionValue <Gradient, &FEMContext::point_gradient>::out_of_elem_tol = 10*TOLERANCE;


/**
 * This class builds the send_list of old dof indices
 * whose coefficients are needed to perform a projection.
 * This may be executed in parallel on multiple threads.
 * The end result is a \p send_list vector which is
 * unsorted and may contain duplicate elements.
 * The \p unique() method can be used to sort and
 * create a unique list.
 */
class BuildProjectionList
{
private:
  const System & system;

public:
  BuildProjectionList (const System & system_in) :
    system(system_in),
    send_list()
  {}

  BuildProjectionList (BuildProjectionList & other, Threads::split) :
    system(other.system),
    send_list()
  {}

  void unique();
  void operator()(const ConstElemRange & range);
  void join (const BuildProjectionList & other);
  std::vector<dof_id_type> send_list;
};

#endif // LIBMESH_ENABLE_AMR


/**
 * This class implements projecting an arbitrary
 * boundary function to the current mesh.  This
 * may be exectued in parallel on multiple threads.
 */
class BoundaryProjectSolution
{
private:
  const std::set<boundary_id_type> & b;
  const std::vector<unsigned int>  & variables;
  const System                     & system;
  UniquePtr<FunctionBase<Number> >   f;
  UniquePtr<FunctionBase<Gradient> > g;
  const Parameters                 & parameters;
  NumericVector<Number>            & new_vector;

public:
  BoundaryProjectSolution (const std::set<boundary_id_type> & b_in,
                           const std::vector<unsigned int> & variables_in,
                           const System & system_in,
                           FunctionBase<Number> * f_in,
                           FunctionBase<Gradient> * g_in,
                           const Parameters & parameters_in,
                           NumericVector<Number> & new_v_in) :
    b(b_in),
    variables(variables_in),
    system(system_in),
    f(f_in ? f_in->clone() : UniquePtr<FunctionBase<Number> >()),
    g(g_in ? g_in->clone() : UniquePtr<FunctionBase<Gradient> >()),
    parameters(parameters_in),
    new_vector(new_v_in)
  {
    libmesh_assert(f.get());
    f->init();
    if (g.get())
      g->init();
  }

  BoundaryProjectSolution (const BoundaryProjectSolution & in) :
    b(in.b),
    variables(in.variables),
    system(in.system),
    f(in.f.get() ? in.f->clone() : UniquePtr<FunctionBase<Number> >()),
    g(in.g.get() ? in.g->clone() : UniquePtr<FunctionBase<Gradient> >()),
    parameters(in.parameters),
    new_vector(in.new_vector)
  {
    libmesh_assert(f.get());
    f->init();
    if (g.get())
      g->init();
  }

  void operator()(const ConstElemRange & range) const;
};



// ------------------------------------------------------------
// System implementation
void System::project_vector (NumericVector<Number> & vector,
                             int is_adjoint) const
{
  // Create a copy of the vector, which currently
  // contains the old data.
  UniquePtr<NumericVector<Number> >
    old_vector (vector.clone());

  // Project the old vector to the new vector
  this->project_vector (*old_vector, vector, is_adjoint);
}


/**
 * This method projects the vector
 * via L2 projections or nodal
 * interpolations on each element.
 */
void System::project_vector (const NumericVector<Number> & old_v,
                             NumericVector<Number> & new_v,
                             int
#ifdef LIBMESH_ENABLE_AMR
                             is_adjoint
#endif
                             ) const
{
  LOG_SCOPE ("project_vector(old,new)", "System");

  /**
   * This method projects a solution from an old mesh to a current, refined
   * mesh.  The input vector \p old_v gives the solution on the
   * old mesh, while the \p new_v gives the solution (to be computed)
   * on the new mesh.
   */
  new_v.clear();

#ifdef LIBMESH_ENABLE_AMR

  // Resize the new vector and get a serial version.
  NumericVector<Number> * new_vector_ptr = libmesh_nullptr;
  UniquePtr<NumericVector<Number> > new_vector_built;
  NumericVector<Number> * local_old_vector;
  UniquePtr<NumericVector<Number> > local_old_vector_built;
  const NumericVector<Number> * old_vector_ptr = libmesh_nullptr;

  ConstElemRange active_local_elem_range
    (this->get_mesh().active_local_elements_begin(),
     this->get_mesh().active_local_elements_end());

  // If the old vector was uniprocessor, make the new
  // vector uniprocessor
  if (old_v.type() == SERIAL)
    {
      new_v.init (this->n_dofs(), false, SERIAL);
      new_vector_ptr = &new_v;
      old_vector_ptr = &old_v;
    }

  // Otherwise it is a parallel, distributed vector, which
  // we need to localize.
  else if (old_v.type() == PARALLEL)
    {
      // Build a send list for efficient localization
      BuildProjectionList projection_list(*this);
      Threads::parallel_reduce (active_local_elem_range,
                                projection_list);

      // Create a sorted, unique send_list
      projection_list.unique();

      new_v.init (this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
      new_vector_built = NumericVector<Number>::build(this->comm());
      local_old_vector_built = NumericVector<Number>::build(this->comm());
      new_vector_ptr = new_vector_built.get();
      local_old_vector = local_old_vector_built.get();
      new_vector_ptr->init(this->n_dofs(), false, SERIAL);
      local_old_vector->init(old_v.size(), false, SERIAL);
      old_v.localize(*local_old_vector, projection_list.send_list);
      local_old_vector->close();
      old_vector_ptr = local_old_vector;
    }
  else if (old_v.type() == GHOSTED)
    {
      // Build a send list for efficient localization
      BuildProjectionList projection_list(*this);
      Threads::parallel_reduce (active_local_elem_range,
                                projection_list);

      // Create a sorted, unique send_list
      projection_list.unique();

      new_v.init (this->n_dofs(), this->n_local_dofs(),
                  this->get_dof_map().get_send_list(), false, GHOSTED);

      local_old_vector_built = NumericVector<Number>::build(this->comm());
      new_vector_ptr = &new_v;
      local_old_vector = local_old_vector_built.get();
      local_old_vector->init(old_v.size(), old_v.local_size(),
                             projection_list.send_list, false, GHOSTED);
      old_v.localize(*local_old_vector, projection_list.send_list);
      local_old_vector->close();
      old_vector_ptr = local_old_vector;
    }
  else // unknown old_v.type()
    libmesh_error_msg("ERROR: Unknown old_v.type() == " << old_v.type());

  // Note that the above will have zeroed the new_vector.
  // Just to be sure, assert that new_vector_ptr and old_vector_ptr
  // were successfully set before trying to deref them.
  libmesh_assert(new_vector_ptr);
  libmesh_assert(old_vector_ptr);

  NumericVector<Number> & new_vector = *new_vector_ptr;
  const NumericVector<Number> & old_vector = *old_vector_ptr;

  const unsigned int n_variables = this->n_vars();

  if (n_variables)
    {
      std::vector<unsigned int> vars(n_variables);
      for (unsigned int i=0; i != n_variables; ++i)
        vars[i] = i;

      // Use a typedef to make the calling sequence for parallel_for() a bit more readable
      typedef
        GenericProjector<OldSolutionValue<Number,   &FEMContext::point_value>,
                         OldSolutionValue<Gradient, &FEMContext::point_gradient>,
                         Number, VectorSetAction<Number> > FEMProjector;

      OldSolutionValue<Number,   &FEMContext::point_value>    f(*this, old_vector);
      OldSolutionValue<Gradient, &FEMContext::point_gradient> g(*this, old_vector);
      VectorSetAction<Number> setter(new_vector);

      Threads::parallel_for (active_local_elem_range,
                             FEMProjector(*this, f, &g, setter, vars));

      // Copy the SCALAR dofs from old_vector to new_vector
      // Note: We assume that all SCALAR dofs are on the
      // processor with highest ID
      if(this->processor_id() == (this->n_processors()-1))
        {
          const DofMap & dof_map = this->get_dof_map();
          for (unsigned int var=0; var<this->n_vars(); var++)
            if(this->variable(var).type().family == SCALAR)
              {
                // We can just map SCALAR dofs directly across
                std::vector<dof_id_type> new_SCALAR_indices, old_SCALAR_indices;
                dof_map.SCALAR_dof_indices (new_SCALAR_indices, var, false);
                dof_map.SCALAR_dof_indices (old_SCALAR_indices, var, true);
                const unsigned int new_n_dofs =
                  cast_int<unsigned int>(new_SCALAR_indices.size());

                for (unsigned int i=0; i<new_n_dofs; i++)
                  {
                    new_vector.set( new_SCALAR_indices[i], old_vector(old_SCALAR_indices[i]) );
                  }
              }
        }
    }

  new_vector.close();

  // If the old vector was serial, we probably need to send our values
  // to other processors
  //
  // FIXME: I'm not sure how to make a NumericVector do that without
  // creating a temporary parallel vector to use localize! - RHS
  if (old_v.type() == SERIAL)
    {
      UniquePtr<NumericVector<Number> > dist_v = NumericVector<Number>::build(this->comm());
      dist_v->init(this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
      dist_v->close();

      for (dof_id_type i=0; i!=dist_v->size(); i++)
        if (new_vector(i) != 0.0)
          dist_v->set(i, new_vector(i));

      dist_v->close();

      dist_v->localize (new_v, this->get_dof_map().get_send_list());
      new_v.close();
    }
  // If the old vector was parallel, we need to update it
  // and free the localized copies
  else if (old_v.type() == PARALLEL)
    {
      // We may have to set dof values that this processor doesn't
      // own in certain special cases, like LAGRANGE FIRST or
      // HERMITE THIRD elements on second-order meshes
      for (dof_id_type i=0; i!=new_v.size(); i++)
        if (new_vector(i) != 0.0)
          new_v.set(i, new_vector(i));
      new_v.close();
    }

  if (is_adjoint == -1)
    this->get_dof_map().enforce_constraints_exactly(*this, &new_v);
  else if (is_adjoint >= 0)
    this->get_dof_map().enforce_adjoint_constraints_exactly(new_v,
                                                            is_adjoint);

#else

  // AMR is disabled: simply copy the vector
  new_v = old_v;

#endif // #ifdef LIBMESH_ENABLE_AMR
}



/**
 * This method projects an arbitrary function onto the solution via L2
 * projections and nodal interpolations on each element.
 */
void System::project_solution (Number fptr(const Point & p,
                                           const Parameters & parameters,
                                           const std::string & sys_name,
                                           const std::string & unknown_name),
                               Gradient gptr(const Point & p,
                                             const Parameters & parameters,
                                             const std::string & sys_name,
                                             const std::string & unknown_name),
                               const Parameters & parameters) const
{
  WrappedFunction<Number> f(*this, fptr, &parameters);
  WrappedFunction<Gradient> g(*this, gptr, &parameters);
  this->project_solution(&f, &g);
}


/**
 * This method projects an arbitrary function onto the solution via L2
 * projections and nodal interpolations on each element.
 */
void System::project_solution (FunctionBase<Number> * f,
                               FunctionBase<Gradient> * g) const
{
  this->project_vector(*solution, f, g);

  solution->localize(*current_local_solution, _dof_map->get_send_list());
}


/**
 * This method projects an arbitrary function onto the solution via L2
 * projections and nodal interpolations on each element.
 */
void System::project_solution (FEMFunctionBase<Number> * f,
                               FEMFunctionBase<Gradient> * g) const
{
  this->project_vector(*solution, f, g);

  solution->localize(*current_local_solution, _dof_map->get_send_list());
}


/**
 * This method projects an arbitrary function via L2 projections and
 * nodal interpolations on each element.
 */
void System::project_vector (Number fptr(const Point & p,
                                         const Parameters & parameters,
                                         const std::string & sys_name,
                                         const std::string & unknown_name),
                             Gradient gptr(const Point & p,
                                           const Parameters & parameters,
                                           const std::string & sys_name,
                                           const std::string & unknown_name),
                             const Parameters & parameters,
                             NumericVector<Number> & new_vector,
                             int is_adjoint) const
{
  WrappedFunction<Number> f(*this, fptr, &parameters);
  WrappedFunction<Gradient> g(*this, gptr, &parameters);
  this->project_vector(new_vector, &f, &g, is_adjoint);
}

/**
 * This method projects an arbitrary function via L2 projections and
 * nodal interpolations on each element.
 */
void System::project_vector (NumericVector<Number> & new_vector,
                             FunctionBase<Number> * f,
                             FunctionBase<Gradient> * g,
                             int is_adjoint) const
{
  LOG_SCOPE ("project_vector(FunctionBase)", "System");

  libmesh_assert(f);

  WrappedFunctor<Number> f_fem(*f);

  if (g)
    {
      WrappedFunctor<Gradient> g_fem(*g);

      this->project_vector(new_vector, &f_fem, &g_fem, is_adjoint);
    }
  else
    this->project_vector(new_vector, &f_fem, libmesh_nullptr, is_adjoint);
}


/**
 * This method projects an arbitrary function via L2 projections and
 * nodal interpolations on each element.
 */
void System::project_vector (NumericVector<Number> & new_vector,
                             FEMFunctionBase<Number> * f,
                             FEMFunctionBase<Gradient> * g,
                             int is_adjoint) const
{
  LOG_SCOPE ("project_fem_vector()", "System");

  libmesh_assert (f);

  ConstElemRange active_local_range
    (this->get_mesh().active_local_elements_begin(),
     this->get_mesh().active_local_elements_end() );

  VectorSetAction<Number> setter(new_vector);

  const unsigned int n_variables = this->n_vars();

  std::vector<unsigned int> vars(n_variables);
  for (unsigned int i=0; i != n_variables; ++i)
    vars[i] = i;

  // Use a typedef to make the calling sequence for parallel_for() a bit more readable
  typedef
    GenericProjector<FEMFunctionWrapper<Number>, FEMFunctionWrapper<Gradient>,
                     Number, VectorSetAction<Number> > FEMProjector;

  FEMFunctionWrapper<Number> fw(*f);

  if (g)
    {
      FEMFunctionWrapper<Gradient> gw(*g);

      Threads::parallel_for
        (active_local_range,
         FEMProjector(*this, fw, &gw, setter, vars));
    }
  else
    Threads::parallel_for
      (active_local_range,
       FEMProjector(*this, fw, libmesh_nullptr, setter, vars));

  // Also, load values into the SCALAR dofs
  // Note: We assume that all SCALAR dofs are on the
  // processor with highest ID
  if(this->processor_id() == (this->n_processors()-1))
    {
      // FIXME: Do we want to first check for SCALAR vars before building this? [PB]
      FEMContext context( *this );

      const DofMap & dof_map = this->get_dof_map();
      for (unsigned int var=0; var<this->n_vars(); var++)
        if(this->variable(var).type().family == SCALAR)
          {
            // FIXME: We reinit with an arbitrary element in case the user
            //        doesn't override FEMFunctionBase::component. Is there
            //        any use case we're missing? [PB]
            Elem * el = const_cast<Elem *>(*(this->get_mesh().active_local_elements_begin()));
            context.pre_fe_reinit(*this, el);

            std::vector<dof_id_type> SCALAR_indices;
            dof_map.SCALAR_dof_indices (SCALAR_indices, var);
            const unsigned int n_SCALAR_dofs =
              cast_int<unsigned int>(SCALAR_indices.size());

            for (unsigned int i=0; i<n_SCALAR_dofs; i++)
              {
                const dof_id_type global_index = SCALAR_indices[i];
                const unsigned int component_index =
                  this->variable_scalar_number(var,i);

                new_vector.set(global_index, f->component(context, component_index, Point(), this->time));
              }
          }
    }

  new_vector.close();

#ifdef LIBMESH_ENABLE_CONSTRAINTS
  if (is_adjoint == -1)
    this->get_dof_map().enforce_constraints_exactly(*this, &new_vector);
  else if (is_adjoint >= 0)
    this->get_dof_map().enforce_adjoint_constraints_exactly(new_vector,
                                                            is_adjoint);
#endif
}


/**
 * This method projects components of an arbitrary boundary function
 * onto the solution via L2 projections and nodal interpolations on
 * each element.
 */
void System::boundary_project_solution (const std::set<boundary_id_type> & b,
                                        const std::vector<unsigned int> & variables,
                                        Number fptr(const Point & p,
                                                    const Parameters & parameters,
                                                    const std::string & sys_name,
                                                    const std::string & unknown_name),
                                        Gradient gptr(const Point & p,
                                                      const Parameters & parameters,
                                                      const std::string & sys_name,
                                                      const std::string & unknown_name),
                                        const Parameters & parameters)
{
  WrappedFunction<Number> f(*this, fptr, &parameters);
  WrappedFunction<Gradient> g(*this, gptr, &parameters);
  this->boundary_project_solution(b, variables, &f, &g);
}


/**
 * This method projects an arbitary boundary function onto the
 * solution via L2 projections and nodal interpolations on each
 * element.
 */
void System::boundary_project_solution (const std::set<boundary_id_type> & b,
                                        const std::vector<unsigned int> & variables,
                                        FunctionBase<Number> * f,
                                        FunctionBase<Gradient> * g)
{
  this->boundary_project_vector(b, variables, *solution, f, g);

  solution->localize(*current_local_solution);
}





/**
 * This method projects an arbitrary boundary function via L2
 * projections and nodal interpolations on each element.
 */
void System::boundary_project_vector (const std::set<boundary_id_type> & b,
                                      const std::vector<unsigned int> & variables,
                                      Number fptr(const Point & p,
                                                  const Parameters & parameters,
                                                  const std::string & sys_name,
                                                  const std::string & unknown_name),
                                      Gradient gptr(const Point & p,
                                                    const Parameters & parameters,
                                                    const std::string & sys_name,
                                                    const std::string & unknown_name),
                                      const Parameters & parameters,
                                      NumericVector<Number> & new_vector,
                                      int is_adjoint) const
{
  WrappedFunction<Number> f(*this, fptr, &parameters);
  WrappedFunction<Gradient> g(*this, gptr, &parameters);
  this->boundary_project_vector(b, variables, new_vector, &f, &g,
                                is_adjoint);
}

/**
 * This method projects an arbitrary function via L2 projections and
 * nodal interpolations on each element.
 */
void System::boundary_project_vector (const std::set<boundary_id_type> & b,
                                      const std::vector<unsigned int> & variables,
                                      NumericVector<Number> & new_vector,
                                      FunctionBase<Number> * f,
                                      FunctionBase<Gradient> * g,
                                      int is_adjoint) const
{
  LOG_SCOPE ("boundary_project_vector()", "System");

  Threads::parallel_for
    (ConstElemRange (this->get_mesh().active_local_elements_begin(),
                     this->get_mesh().active_local_elements_end() ),
     BoundaryProjectSolution(b, variables, *this, f, g,
                             this->get_equation_systems().parameters,
                             new_vector)
     );

  // We don't do SCALAR dofs when just projecting the boundary, so
  // we're done here.

  new_vector.close();

#ifdef LIBMESH_ENABLE_CONSTRAINTS
  if (is_adjoint == -1)
    this->get_dof_map().enforce_constraints_exactly(*this, &new_vector);
  else if (is_adjoint >= 0)
    this->get_dof_map().enforce_adjoint_constraints_exactly(new_vector,
                                                            is_adjoint);
#endif
}



template <typename FFunctor, typename GFunctor,
          typename FValue, typename ProjectionAction>
void GenericProjector<FFunctor, GFunctor, FValue, ProjectionAction>::operator()
  (const ConstElemRange & range) const
{
  LOG_SCOPE ("operator()", "GenericProjector");

  ProjectionAction action(master_action);
  FFunctor f(master_f);
  UniquePtr<GFunctor> g;
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
          if( dim > 1 )
            context.get_side_fe( var, side_fe, dim );
          if( dim > 2 )
            context.get_edge_fe( var, edge_fe );

          fe->get_xyz();

          fe->get_phi();
          if( dim > 1 )
            side_fe->get_phi();
          if( dim > 2 )
            edge_fe->get_phi();

          const FEContinuity cont = fe->get_continuity();
          if(cont == C_ONE)
            {
              // Our C1 elements need gradient information
              libmesh_assert(g);

              fe->get_dphi();
              if( dim > 1 )
                side_fe->get_dphi();
              if( dim > 2 )
                edge_fe->get_dphi();
            }
        }
    }

  // Now initialize any user requested shape functions, xyz vals, etc
  f.init_context(context);
  if(g.get())
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
          if( dim > 1 )
            context.get_side_fe( var, side_fe, dim );
          if( dim > 2 )
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
                    elem->parent()->is_vertex_on_parent
                      (i_am_child, n);

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

              const std::vector<std::vector<Real> > & phi =
#ifdef LIBMESH_ENABLE_AMR
                (elem->refinement_flag() == Elem::JUST_COARSENED) ?
                fe->get_phi() :
#endif
                edge_fe->get_phi();
              const std::vector<std::vector<RealGradient> > * dphi = libmesh_nullptr;
              if (cont == C_ONE)
                dphi =
#ifdef LIBMESH_ENABLE_AMR
                  (elem->refinement_flag() == Elem::JUST_COARSENED) ?
                  &(fe->get_dphi()) :
#endif
                  &(edge_fe->get_dphi());

              for (unsigned char e=0; e != elem->n_edges(); ++e)
                {
                  context.edge = e;

#ifdef LIBMESH_ENABLE_AMR
                  if (elem->refinement_flag() == Elem::JUST_COARSENED)
                    {
                      std::vector<Point> fine_points;

                      UniquePtr<FEBase> fine_fe
                        (FEBase::build (dim, base_fe_type));

                      UniquePtr<QBase> qrule
                        (base_fe_type.default_quadrature_rule(1));
                      fine_fe->attach_quadrature_rule(qrule.get());

                      const std::vector<Point> & child_xyz =
                        fine_fe->get_xyz();

                      for (unsigned int c = 0;
                           c != elem->n_children(); ++c)
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
                      libmesh_assert(std::abs(ui) < TOLERANCE ||
                                     std::abs(ui - Uedge(i)) < TOLERANCE);
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
              const std::vector<std::vector<Real> > & phi =
#ifdef LIBMESH_ENABLE_AMR
                (elem->refinement_flag() == Elem::JUST_COARSENED) ?
                fe->get_phi() :
#endif // LIBMESH_ENABLE_AMR
                side_fe->get_phi();
              const std::vector<std::vector<RealGradient> > * dphi = libmesh_nullptr;
              if (cont == C_ONE)
                dphi =
#ifdef LIBMESH_ENABLE_AMR
                  (elem->refinement_flag() == Elem::JUST_COARSENED) ?
                    &(fe->get_dphi()) :
#endif // LIBMESH_ENABLE_AMR
                    &(side_fe->get_dphi());

              for (unsigned char s=0; s != elem->n_sides(); ++s)
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

                      UniquePtr<FEBase> fine_fe
                        (FEBase::build (dim, base_fe_type));

                      UniquePtr<QBase> qrule
                        (base_fe_type.default_quadrature_rule(dim-1));
                      fine_fe->attach_quadrature_rule(qrule.get());

                      const std::vector<Point> & child_xyz =
                        fine_fe->get_xyz();

                      for (unsigned int c = 0;
                           c != elem->n_children(); ++c)
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
                      libmesh_assert(std::abs(ui) < TOLERANCE ||
                                     std::abs(ui - Uside(i)) < TOLERANCE);
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

              const std::vector<std::vector<Real> > & phi = fe->get_phi();
              const std::vector<std::vector<RealGradient> > * dphi = libmesh_nullptr;
              if (cont == C_ONE)
                dphi = &(fe->get_dphi());

#ifdef LIBMESH_ENABLE_AMR
              if (elem->refinement_flag() == Elem::JUST_COARSENED)
                {
                  std::vector<Point> fine_points;

                  UniquePtr<FEBase> fine_fe
                    (FEBase::build (dim, base_fe_type));

                  UniquePtr<QBase> qrule
                    (base_fe_type.default_quadrature_rule(dim));
                  fine_fe->attach_quadrature_rule(qrule.get());

                  const std::vector<Point> & child_xyz =
                    fine_fe->get_xyz();

                  for (unsigned int c = 0;
                       c != elem->n_children(); ++c)
                    {
                      fine_fe->reinit(elem->child_ptr(c));
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
                  libmesh_assert(std::abs(ui) < TOLERANCE ||
                                 std::abs(ui - Uint(i)) < TOLERANCE);
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



#ifdef LIBMESH_ENABLE_AMR
void BuildProjectionList::unique()
{
  // Sort the send list.  After this duplicated
  // elements will be adjacent in the vector
  std::sort(this->send_list.begin(),
            this->send_list.end());

  // Now use std::unique to remove duplicate entries
  std::vector<dof_id_type>::iterator new_end =
    std::unique (this->send_list.begin(),
                 this->send_list.end());

  // Remove the end of the send_list.  Use the "swap trick"
  // from Effective STL
  std::vector<dof_id_type>
    (this->send_list.begin(), new_end).swap (this->send_list);
}



void BuildProjectionList::operator()(const ConstElemRange & range)
{
  // The DofMap for this system
  const DofMap & dof_map = system.get_dof_map();

  const dof_id_type first_old_dof = dof_map.first_old_dof();
  const dof_id_type end_old_dof   = dof_map.end_old_dof();

  // We can handle all the variables at once.
  // The old global DOF indices
  std::vector<dof_id_type> di;

  // Iterate over the elements in the range
  for (ConstElemRange::const_iterator elem_it=range.begin(); elem_it != range.end(); ++elem_it)
    {
      const Elem * elem = *elem_it;
      // If this element doesn't have an old_dof_object with dofs for the
      // current system, then it must be newly added, so the user
      // is responsible for setting the new dofs.

      // ... but we need a better way to test for that; the code
      // below breaks on any FE type for which the elem stores no
      // dofs.
      // if (!elem->old_dof_object || !elem->old_dof_object->has_dofs(system.number()))
      //  continue;

      // Examining refinement flags instead should distinguish
      // between refinement-added and user-added elements lacking
      // old_dof_object
      if (!elem->old_dof_object &&
          elem->refinement_flag() != Elem::JUST_REFINED &&
          elem->refinement_flag() != Elem::JUST_COARSENED)
        continue;

      const Elem * parent = elem->parent();

      if (elem->refinement_flag() == Elem::JUST_REFINED)
        {
          libmesh_assert(parent);

          // We used to hack_p_level here, but that wasn't thread-safe
          // so now we take p refinement flags into account in
          // old_dof_indices

          dof_map.old_dof_indices (parent, di);

          for (unsigned int n=0; n != elem->n_nodes(); ++n)
            {
              const Node * node = elem->node_ptr(n);
              const DofObject * old_dofs = node->old_dof_object;

              if (old_dofs)
                {
                  const unsigned int sysnum = system.number();
                  const unsigned int nv = old_dofs->n_vars(sysnum);
                  for (unsigned int v=0; v != nv; ++v)
                    {
                      const unsigned int nc =
                        old_dofs->n_comp(sysnum, v);
                      for (unsigned int c=0; c != nc; ++c)
                        di.push_back
                          (old_dofs->dof_number(sysnum, v, c));
                    }
                }
            }

          std::sort(di.begin(), di.end());
          std::vector<dof_id_type>::iterator new_end =
            std::unique(di.begin(), di.end());
          std::vector<dof_id_type>(di.begin(), new_end).swap(di);
        }
      else if (elem->refinement_flag() == Elem::JUST_COARSENED)
        {
          std::vector<dof_id_type> di_child;
          di.clear();
          for (unsigned int c=0; c != elem->n_children(); ++c)
            {
              dof_map.old_dof_indices (elem->child_ptr(c), di_child);
              di.insert(di.end(), di_child.begin(), di_child.end());
            }
        }
      else
        dof_map.old_dof_indices (elem, di);

      for (std::size_t i=0; i != di.size(); ++i)
        if (di[i] < first_old_dof || di[i] >= end_old_dof)
          this->send_list.push_back(di[i]);
    }  // end elem loop
}



void BuildProjectionList::join(const BuildProjectionList & other)
{
  // Joining simply requires I add the dof indices from the other object
  this->send_list.insert(this->send_list.end(),
                         other.send_list.begin(),
                         other.send_list.end());
}
#endif // LIBMESH_ENABLE_AMR



void BoundaryProjectSolution::operator()(const ConstElemRange & range) const
{
  // We need data to project
  libmesh_assert(f.get());

  /**
   * This method projects an arbitrary boundary solution to the current
   * mesh.  The input function \p f gives the arbitrary solution,
   * while the \p new_vector (which should already be correctly sized)
   * gives the solution (to be computed) on the current mesh.
   */

  // The dimensionality of the current mesh
  const unsigned int dim = system.get_mesh().mesh_dimension();

  // The DofMap for this system
  const DofMap & dof_map = system.get_dof_map();

  // Boundary info for the current mesh
  const BoundaryInfo & boundary_info =
    system.get_mesh().get_boundary_info();

  // The element matrix and RHS for projections.
  // Note that Ke is always real-valued, whereas
  // Fe may be complex valued if complex number
  // support is enabled
  DenseMatrix<Real> Ke;
  DenseVector<Number> Fe;
  // The new element coefficients
  DenseVector<Number> Ue;


  // Loop over all the variables we've been requested to project
  for (std::size_t v=0; v!=variables.size(); v++)
    {
      const unsigned int var = variables[v];

      const Variable & variable = dof_map.variable(var);

      const FEType & fe_type = variable.type();

      if (fe_type.family == SCALAR)
        continue;

      const unsigned int var_component =
        system.variable_scalar_number(var, 0);

      // Get FE objects of the appropriate type
      UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));

      // Prepare variables for projection
      UniquePtr<QBase> qedgerule (fe_type.default_quadrature_rule(1));
      UniquePtr<QBase> qsiderule (fe_type.default_quadrature_rule(dim-1));

      // The values of the shape functions at the quadrature
      // points
      const std::vector<std::vector<Real> > & phi = fe->get_phi();

      // The gradients of the shape functions at the quadrature
      // points on the child element.
      const std::vector<std::vector<RealGradient> > * dphi = libmesh_nullptr;

      const FEContinuity cont = fe->get_continuity();

      if (cont == C_ONE)
        {
          // We'll need gradient data for a C1 projection
          libmesh_assert(g.get());

          const std::vector<std::vector<RealGradient> > &
            ref_dphi = fe->get_dphi();
          dphi = &ref_dphi;
        }

      // The Jacobian * quadrature weight at the quadrature points
      const std::vector<Real> & JxW =
        fe->get_JxW();

      // The XYZ locations of the quadrature points
      const std::vector<Point> & xyz_values =
        fe->get_xyz();

      // The global DOF indices
      std::vector<dof_id_type> dof_indices;
      // Side/edge DOF indices
      std::vector<unsigned int> side_dofs;

      // Container to catch IDs passed back from BoundaryInfo.
      std::vector<boundary_id_type> bc_ids;

      // Iterate over all the elements in the range
      for (ConstElemRange::const_iterator elem_it=range.begin(); elem_it != range.end(); ++elem_it)
        {
          const Elem * elem = *elem_it;

          // Per-subdomain variables don't need to be projected on
          // elements where they're not active
          if (!variable.active_on_subdomain(elem->subdomain_id()))
            continue;

          // Find out which nodes, edges and sides are on a requested
          // boundary:
          std::vector<bool> is_boundary_node(elem->n_nodes(), false),
            is_boundary_edge(elem->n_edges(), false),
            is_boundary_side(elem->n_sides(), false);
          for (unsigned char s=0; s != elem->n_sides(); ++s)
            {
              // First see if this side has been requested
              boundary_info.boundary_ids (elem, s, bc_ids);
              bool do_this_side = false;
              for (std::vector<boundary_id_type>::iterator i=bc_ids.begin();
                   i!=bc_ids.end(); ++i)
                if (b.count(*i))
                  {
                    do_this_side = true;
                    break;
                  }
              if (!do_this_side)
                continue;

              is_boundary_side[s] = true;

              // Then see what nodes and what edges are on it
              for (unsigned int n=0; n != elem->n_nodes(); ++n)
                if (elem->is_node_on_side(n,s))
                  is_boundary_node[n] = true;
              for (unsigned int e=0; e != elem->n_edges(); ++e)
                if (elem->is_edge_on_side(e,s))
                  is_boundary_edge[e] = true;
            }

          // Update the DOF indices for this element based on
          // the current mesh
          dof_map.dof_indices (elem, dof_indices, var);

          // The number of DOFs on the element
          const unsigned int n_dofs =
            cast_int<unsigned int>(dof_indices.size());

          // Fixed vs. free DoFs on edge/face projections
          std::vector<char> dof_is_fixed(n_dofs, false); // bools
          std::vector<int> free_dof(n_dofs, 0);

          // The element type
          const ElemType elem_type = elem->type();

          // The number of nodes on the new element
          const unsigned int n_nodes = elem->n_nodes();

          // Zero the interpolated values
          Ue.resize (n_dofs); Ue.zero();

          // In general, we need a series of
          // projections to ensure a unique and continuous
          // solution.  We start by interpolating boundary nodes, then
          // hold those fixed and project boundary edges, then hold
          // those fixed and project boundary faces,

          // Interpolate node values first
          unsigned int current_dof = 0;
          for (unsigned int n=0; n!= n_nodes; ++n)
            {
              // FIXME: this should go through the DofMap,
              // not duplicate dof_indices code badly!
              const unsigned int nc =
                FEInterface::n_dofs_at_node (dim, fe_type, elem_type,
                                             n);
              if (!elem->is_vertex(n) || !is_boundary_node[n])
                {
                  current_dof += nc;
                  continue;
                }
              if (cont == DISCONTINUOUS)
                {
                  libmesh_assert_equal_to (nc, 0);
                }
              // Assume that C_ZERO elements have a single nodal
              // value shape function
              else if (cont == C_ZERO)
                {
                  libmesh_assert_equal_to (nc, 1);
                  Ue(current_dof) = f->component(var_component,
                                                 elem->point(n),
                                                 system.time);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                }
              // The hermite element vertex shape functions are weird
              else if (fe_type.family == HERMITE)
                {
                  Ue(current_dof) = f->component(var_component,
                                                 elem->point(n),
                                                 system.time);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                  Gradient grad = g->component(var_component,
                                               elem->point(n),
                                               system.time);
                  // x derivative
                  Ue(current_dof) = grad(0);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                  if (dim > 1)
                    {
                      // We'll finite difference mixed derivatives
                      Point nxminus = elem->point(n),
                        nxplus = elem->point(n);
                      nxminus(0) -= TOLERANCE;
                      nxplus(0) += TOLERANCE;
                      Gradient gxminus = g->component(var_component,
                                                      nxminus,
                                                      system.time);
                      Gradient gxplus = g->component(var_component,
                                                     nxplus,
                                                     system.time);
                      // y derivative
                      Ue(current_dof) = grad(1);
                      dof_is_fixed[current_dof] = true;
                      current_dof++;
                      // xy derivative
                      Ue(current_dof) = (gxplus(1) - gxminus(1))
                        / 2. / TOLERANCE;
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
                            / 2. / TOLERANCE;
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                          // We need new points for yz
                          Point nyminus = elem->point(n),
                            nyplus = elem->point(n);
                          nyminus(1) -= TOLERANCE;
                          nyplus(1) += TOLERANCE;
                          Gradient gyminus = g->component(var_component,
                                                          nyminus,
                                                          system.time);
                          Gradient gyplus = g->component(var_component,
                                                         nyplus,
                                                         system.time);
                          // xz derivative
                          Ue(current_dof) = (gyplus(2) - gyminus(2))
                            / 2. / TOLERANCE;
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                          // Getting a 2nd order xyz is more tedious
                          Point nxmym = elem->point(n),
                            nxmyp = elem->point(n),
                            nxpym = elem->point(n),
                            nxpyp = elem->point(n);
                          nxmym(0) -= TOLERANCE;
                          nxmym(1) -= TOLERANCE;
                          nxmyp(0) -= TOLERANCE;
                          nxmyp(1) += TOLERANCE;
                          nxpym(0) += TOLERANCE;
                          nxpym(1) -= TOLERANCE;
                          nxpyp(0) += TOLERANCE;
                          nxpyp(1) += TOLERANCE;
                          Gradient gxmym = g->component(var_component,
                                                        nxmym,
                                                        system.time);
                          Gradient gxmyp = g->component(var_component,
                                                        nxmyp,
                                                        system.time);
                          Gradient gxpym = g->component(var_component,
                                                        nxpym,
                                                        system.time);
                          Gradient gxpyp = g->component(var_component,
                                                        nxpyp,
                                                        system.time);
                          Number gxzplus = (gxpyp(2) - gxmyp(2))
                            / 2. / TOLERANCE;
                          Number gxzminus = (gxpym(2) - gxmym(2))
                            / 2. / TOLERANCE;
                          // xyz derivative
                          Ue(current_dof) = (gxzplus - gxzminus)
                            / 2. / TOLERANCE;
                          dof_is_fixed[current_dof] = true;
                          current_dof++;
                        }
                    }
                }
              // Assume that other C_ONE elements have a single nodal
              // value shape function and nodal gradient component
              // shape functions
              else if (cont == C_ONE)
                {
                  libmesh_assert_equal_to (nc, 1 + dim);
                  Ue(current_dof) = f->component(var_component,
                                                 elem->point(n),
                                                 system.time);
                  dof_is_fixed[current_dof] = true;
                  current_dof++;
                  Gradient grad = g->component(var_component,
                                               elem->point(n),
                                               system.time);
                  for (unsigned int i=0; i!= dim; ++i)
                    {
                      Ue(current_dof) = grad(i);
                      dof_is_fixed[current_dof] = true;
                      current_dof++;
                    }
                }
              else
                libmesh_error_msg("Unknown continuity " << cont);
            }

          // In 3D, project any edge values next
          if (dim > 2 && cont != DISCONTINUOUS)
            for (unsigned int e=0; e != elem->n_edges(); ++e)
              {
                if (!is_boundary_edge[e])
                  continue;

                FEInterface::dofs_on_edge(elem, dim, fe_type, e,
                                          side_dofs);

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
                DenseVector<Number> Uedge(free_dofs);

                // Initialize FE data on the edge
                fe->attach_quadrature_rule (qedgerule.get());
                fe->edge_reinit (elem, e);
                const unsigned int n_qp = qedgerule->n_points();

                // Loop over the quadrature points
                for (unsigned int qp=0; qp<n_qp; qp++)
                  {
                    // solution at the quadrature point
                    Number fineval = f->component(var_component,
                                                  xyz_values[qp],
                                                  system.time);
                    // solution grad at the quadrature point
                    Gradient finegrad;
                    if (cont == C_ONE)
                      finegrad = g->component(var_component,
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
                                  Fe(freei) -= ((*dphi)[i][qp] *
                                                (*dphi)[j][qp]) *
                                    JxW[qp] * Ue(j);
                                else
                                  Ke(freei,freej) += ((*dphi)[i][qp] *
                                                      (*dphi)[j][qp])
                                    * JxW[qp];
                              }
                            if (!dof_is_fixed[j])
                              freej++;
                          }
                        Fe(freei) += phi[i][qp] * fineval * JxW[qp];
                        if (cont == C_ONE)
                          Fe(freei) += (finegrad * (*dphi)[i][qp]) *
                            JxW[qp];
                        freei++;
                      }
                  }

                Ke.cholesky_solve(Fe, Uedge);

                // Transfer new edge solutions to element
                for (unsigned int i=0; i != free_dofs; ++i)
                  {
                    Number & ui = Ue(side_dofs[free_dof[i]]);
                    libmesh_assert(std::abs(ui) < TOLERANCE ||
                                   std::abs(ui - Uedge(i)) < TOLERANCE);
                    ui = Uedge(i);
                    dof_is_fixed[side_dofs[free_dof[i]]] = true;
                  }
              }

          // Project any side values (edges in 2D, faces in 3D)
          if (dim > 1 && cont != DISCONTINUOUS)
            for (unsigned int s=0; s != elem->n_sides(); ++s)
              {
                if (!is_boundary_side[s])
                  continue;

                FEInterface::dofs_on_side(elem, dim, fe_type, s,
                                          side_dofs);

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
                DenseVector<Number> Uside(free_dofs);

                // Initialize FE data on the side
                fe->attach_quadrature_rule (qsiderule.get());
                fe->reinit (elem, s);
                const unsigned int n_qp = qsiderule->n_points();

                // Loop over the quadrature points
                for (unsigned int qp=0; qp<n_qp; qp++)
                  {
                    // solution at the quadrature point
                    Number fineval = f->component(var_component,
                                                  xyz_values[qp],
                                                  system.time);
                    // solution grad at the quadrature point
                    Gradient finegrad;
                    if (cont == C_ONE)
                      finegrad = g->component(var_component,
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
                                  Fe(freei) -= ((*dphi)[i][qp] *
                                                (*dphi)[j][qp]) *
                                    JxW[qp] * Ue(j);
                                else
                                  Ke(freei,freej) += ((*dphi)[i][qp] *
                                                      (*dphi)[j][qp])
                                    * JxW[qp];
                              }
                            if (!dof_is_fixed[j])
                              freej++;
                          }
                        Fe(freei) += (fineval * phi[i][qp]) * JxW[qp];
                        if (cont == C_ONE)
                          Fe(freei) += (finegrad * (*dphi)[i][qp]) *
                            JxW[qp];
                        freei++;
                      }
                  }

                Ke.cholesky_solve(Fe, Uside);

                // Transfer new side solutions to element
                for (unsigned int i=0; i != free_dofs; ++i)
                  {
                    Number & ui = Ue(side_dofs[free_dof[i]]);
                    libmesh_assert(std::abs(ui) < TOLERANCE ||
                                   std::abs(ui - Uside(i)) < TOLERANCE);
                    ui = Uside(i);
                    dof_is_fixed[side_dofs[free_dof[i]]] = true;
                  }
              }

          const dof_id_type
            first = new_vector.first_local_index(),
            last  = new_vector.last_local_index();

          // Lock the new_vector since it is shared among threads.
          {
            Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

            for (unsigned int i = 0; i < n_dofs; i++)
              if (dof_is_fixed[i] &&
                  (dof_indices[i] >= first) &&
                  (dof_indices[i] <  last))
                {
                  new_vector.set(dof_indices[i], Ue(i));
                }
          }
        }  // end elem loop
    } // end variables loop
}


} // namespace libMesh
