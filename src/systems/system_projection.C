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



// C++ includes
#include <vector>
#include <numeric> // std::iota

// Local includes
#include "libmesh/libmesh_config.h"

#ifdef LIBMESH_HAVE_METAPHYSICL

// With quad precision we need the shim function declarations to
// precede the MetaPhysicL use of them
#include "libmesh/libmesh_common.h"

// FIXME - having to do this with MetaPhysicL brings me shame - RHS
#include "libmesh/ignore_warnings.h"
// Template specialization declarations in here need to *precede* code
// using them.
#include "metaphysicl/dynamicsparsenumberarray_decl.h"
#include "libmesh/restore_warnings.h"

#include "libmesh/compare_types.h"
#include "libmesh/int_range.h"

using MetaPhysicL::DynamicSparseNumberArray;

namespace libMesh
{
// From the perspective of libMesh gradient vectors, a DSNA is a
// scalar component
template <typename T, typename IndexType>
struct ScalarTraits<MetaPhysicL::DynamicSparseNumberArray<T,IndexType> >
{
  static const bool value = true;
};

// And although MetaPhysicL knows how to combine DSNA with something
// else, we need to teach libMesh too.
template <typename T, typename IndexType, typename T2>
struct CompareTypes<MetaPhysicL::DynamicSparseNumberArray<T,IndexType>, T2>
{
  typedef typename
  MetaPhysicL::DynamicSparseNumberArray
  <typename CompareTypes<T,T2>::supertype,IndexType> supertype;
};

template <typename T> struct TypeToSend;

template <typename T, typename IndexType>
struct TypeToSend<MetaPhysicL::DynamicSparseNumberArray<T,IndexType>> {
  typedef std::vector<std::pair<IndexType,T>> type;
};

template <typename T, typename IndexType>
const std::vector<std::pair<IndexType,T>>
convert_to_send(MetaPhysicL::DynamicSparseNumberArray<T,IndexType> & in)
{
  const std::size_t in_size = in.size();
  std::vector<std::pair<IndexType,T>> returnval(in_size);

  for (std::size_t i=0; i != in_size; ++i)
    {
      returnval[i].first = in.raw_index(i);
      returnval[i].second = in.raw_at(i);
    }
  return returnval;
}

template <typename SendT, typename T, typename IndexType>
void convert_from_receive (SendT & received,
                           MetaPhysicL::DynamicSparseNumberArray<T,IndexType> & converted)
{
  const std::size_t received_size = received.size();
  converted.resize(received_size);
  for (std::size_t i=0; i != received_size; ++i)
    {
      converted.raw_index(i) = received[i].first;
      converted.raw_at(i) = received[i].second;
    }
}

}


#endif

#include "libmesh/boundary_info.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/generic_projector.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/system.h"
#include "libmesh/threads.h"
#include "libmesh/wrapped_function.h"
#include "libmesh/wrapped_functor.h"
#include "libmesh/fe_interface.h"



#ifdef LIBMESH_HAVE_METAPHYSICL
// FIXME - wrapping MetaPhysicL is shameful - RHS
#include "libmesh/ignore_warnings.h"
// Include MetaPhysicL definitions finally
#include "metaphysicl/dynamicsparsenumberarray.h"
#include "libmesh/restore_warnings.h"

// And make sure we instantiate the methods we'll need to use on them.
#include "libmesh/dense_matrix_impl.h"

namespace libMesh {
typedef DynamicSparseNumberArray<Real, dof_id_type> DSNAN;

template void
DenseMatrix<Real>::cholesky_solve(const DenseVector<DSNAN> &,
                                  DenseVector<DSNAN> &);
template void
DenseMatrix<Real>::_cholesky_back_substitute(const DenseVector<DSNAN> &,
                                             DenseVector<DSNAN> &) const;
}
#endif



namespace libMesh
{

// ------------------------------------------------------------
// Helper class definitions

#ifdef LIBMESH_ENABLE_AMR

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
 * may be executed in parallel on multiple threads.
 */
class BoundaryProjectSolution
{
private:
  const std::set<boundary_id_type> & b;
  const std::vector<unsigned int>  & variables;
  const System                     & system;
  std::unique_ptr<FunctionBase<Number>>   f;
  std::unique_ptr<FunctionBase<Gradient>> g;
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
    f(f_in ? f_in->clone() : std::unique_ptr<FunctionBase<Number>>()),
    g(g_in ? g_in->clone() : std::unique_ptr<FunctionBase<Gradient>>()),
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
    f(in.f.get() ? in.f->clone() : std::unique_ptr<FunctionBase<Number>>()),
    g(in.g.get() ? in.g->clone() : std::unique_ptr<FunctionBase<Gradient>>()),
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
  std::unique_ptr<NumericVector<Number>>
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
                             int is_adjoint) const
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
  NumericVector<Number> * new_vector_ptr = nullptr;
  std::unique_ptr<NumericVector<Number>> new_vector_built;
  NumericVector<Number> * local_old_vector;
  std::unique_ptr<NumericVector<Number>> local_old_vector_built;
  const NumericVector<Number> * old_vector_ptr = nullptr;

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
      new_vector_ptr->init(this->n_dofs(), this->n_local_dofs(),
                           this->get_dof_map().get_send_list(), false,
                           GHOSTED);
      local_old_vector->init(old_v.size(), old_v.local_size(),
                             projection_list.send_list, false, GHOSTED);
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
      std::iota(vars.begin(), vars.end(), 0);
      std::vector<unsigned int> regular_vars, vector_vars;
      for (auto var : vars)
      {
        if (FEInterface::field_type(this->variable_type(var)) == TYPE_SCALAR)
          regular_vars.push_back(var);
        else
          vector_vars.push_back(var);
      }

      // Use a typedef to make the calling sequence for parallel_for() a bit more readable
      typedef
        GenericProjector<OldSolutionValue<Number,   &FEMContext::point_value>,
                         OldSolutionValue<Gradient, &FEMContext::point_gradient>,
                         Number, VectorSetAction<Number>> FEMProjector;

      OldSolutionValue<Number,   &FEMContext::point_value>    f(*this, old_vector);
      OldSolutionValue<Gradient, &FEMContext::point_gradient> g(*this, old_vector);
      VectorSetAction<Number> setter(new_vector);

      FEMProjector projector(*this, f, &g, setter, regular_vars);
      projector.project(active_local_elem_range);

      typedef
        GenericProjector<OldSolutionValue<Gradient,   &FEMContext::point_value>,
                         OldSolutionValue<Tensor, &FEMContext::point_gradient>,
                         Gradient, VectorSetAction<Number>> FEMVectorProjector;

      OldSolutionValue<Gradient, &FEMContext::point_value> f_vector(*this, old_vector);
      OldSolutionValue<Tensor, &FEMContext::point_gradient> g_vector(*this, old_vector);

      FEMVectorProjector vector_projector(*this, f_vector, &g_vector, setter, vector_vars);
      vector_projector.project(active_local_elem_range);

      // Copy the SCALAR dofs from old_vector to new_vector
      // Note: We assume that all SCALAR dofs are on the
      // processor with highest ID
      if (this->processor_id() == (this->n_processors()-1))
        {
          const DofMap & dof_map = this->get_dof_map();
          for (auto var : make_range(this->n_vars()))
            if (this->variable(var).type().family == SCALAR)
              {
                // We can just map SCALAR dofs directly across
                std::vector<dof_id_type> new_SCALAR_indices, old_SCALAR_indices;
                dof_map.SCALAR_dof_indices (new_SCALAR_indices, var, false);
                dof_map.SCALAR_dof_indices (old_SCALAR_indices, var, true);
                for (auto i : index_range(new_SCALAR_indices))
                  new_vector.set(new_SCALAR_indices[i], old_vector(old_SCALAR_indices[i]));
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
      std::unique_ptr<NumericVector<Number>> dist_v = NumericVector<Number>::build(this->comm());
      dist_v->init(this->n_dofs(), this->n_local_dofs(), false, PARALLEL);
      dist_v->close();

      for (auto i : make_range(dist_v->size()))
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
      // HERMITE THIRD elements on second-order meshes?
      new_v = new_vector;
      new_v.close();
    }


  // Apply constraints only if we we are asked to
  if(this->project_with_constraints)
  {
    if (is_adjoint == -1)
    {
      this->get_dof_map().enforce_constraints_exactly(*this, &new_v);
    }
    else if (is_adjoint >= 0)
    {
      this->get_dof_map().enforce_adjoint_constraints_exactly(new_v,
                                                            is_adjoint);
    }
  }
#else

  // AMR is disabled: simply copy the vector
  new_v = old_v;

  libmesh_ignore(is_adjoint);

#endif // #ifdef LIBMESH_ENABLE_AMR
}


#ifdef LIBMESH_ENABLE_AMR
#ifdef LIBMESH_HAVE_METAPHYSICL

template <typename Output>
class DSNAOutput
{
public:
  typedef DynamicSparseNumberArray<Output, dof_id_type> type;
};

template <typename InnerOutput>
class DSNAOutput<VectorValue<InnerOutput> >
{
public:
  typedef VectorValue<DynamicSparseNumberArray<InnerOutput, dof_id_type> > type;
};

/**
 * The OldSolutionCoefs input functor class can be used with
 * GenericProjector to read solution transfer coefficients on a
 * just-refined-and-coarsened mesh.
 *
 * \author Roy H. Stogner
 * \date 2017
 */

template <typename Output,
          void (FEMContext::*point_output) (unsigned int,
                                            const Point &,
                                            Output &,
                                            const Real) const>
class OldSolutionCoefs : public OldSolutionBase<Output, point_output>
{
public:
  typedef typename DSNAOutput<Output>::type DSNA;
  typedef DSNA ValuePushType;
  typedef DSNA FunctorValue;

  OldSolutionCoefs(const libMesh::System & sys_in) :
    OldSolutionBase<Output, point_output>(sys_in)
  {
    this->old_context.set_algebraic_type(FEMContext::OLD_DOFS_ONLY);
  }

  OldSolutionCoefs(const OldSolutionCoefs & in) :
    OldSolutionBase<Output, point_output>(in.sys)
  {
    this->old_context.set_algebraic_type(FEMContext::OLD_DOFS_ONLY);
  }

  DSNA eval_at_node (const FEMContext & c,
                     unsigned int i,
                     unsigned int elem_dim,
                     const Node & n,
                     bool extra_hanging_dofs,
                     Real /* time */ = 0.);

  DSNA eval_at_point(const FEMContext & c,
                     unsigned int i,
                     const Point & p,
                     Real time,
                     bool skip_context_check);

  void eval_mixed_derivatives (const FEMContext & libmesh_dbg_var(c),
                               unsigned int i,
                               unsigned int dim,
                               const Node & n,
                               std::vector<DSNA> & derivs)
  {
    LOG_SCOPE ("eval_mixed_derivatives", "OldSolutionCoefs");

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

        for (auto d_i : index_range(derivs))
          {
            derivs[d_i].resize(1);
            derivs[d_i].raw_at(0) = 1;
            derivs[d_i].raw_index(0) = old_ids[d_i];
          }
      }
    else
      {
        std::fill(derivs.begin(), derivs.end(), 0);
      }
  }


  void eval_old_dofs (const Elem & elem,
                      unsigned int node_num,
                      unsigned int var_num,
                      std::vector<dof_id_type> & indices,
                      std::vector<DSNA> & values)
  {
    LOG_SCOPE ("eval_old_dofs(node)", "OldSolutionCoefs");

    this->sys.get_dof_map().dof_indices(elem, node_num, indices, var_num);

    std::vector<dof_id_type> old_indices;

    this->sys.get_dof_map().old_dof_indices(elem, node_num, old_indices, var_num);

    libmesh_assert_equal_to (old_indices.size(), indices.size());

    values.resize(old_indices.size());

    for (auto i : index_range(values))
      {
        values[i].resize(1);
        values[i].raw_at(0) = 1;
        values[i].raw_index(0) = old_indices[i];
      }
  }


  void eval_old_dofs (const Elem & elem,
                      const FEType & fe_type,
                      unsigned int sys_num,
                      unsigned int var_num,
                      std::vector<dof_id_type> & indices,
                      std::vector<DSNA> & values)
  {
    LOG_SCOPE ("eval_old_dofs(elem)", "OldSolutionCoefs");

    // We're only to be asked for old dofs on elements that can copy
    // them through DO_NOTHING or through refinement.
    const Elem & old_elem =
      (elem.refinement_flag() == Elem::JUST_REFINED) ?
      *elem.parent() : elem;

    // If there are any element-based DOF numbers, get them
    const unsigned int nc =
      FEInterface::n_dofs_per_elem(fe_type, &elem);

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

    values.resize(old_dof_indices.size());

    for (auto i : index_range(values))
      {
        values[i].resize(1);
        values[i].raw_at(0) = 1;
        values[i].raw_index(0) = old_dof_indices[i];
      }
  }
};



template<>
inline
DynamicSparseNumberArray<Real, dof_id_type>
OldSolutionCoefs<Real, &FEMContext::point_value>::
eval_at_point(const FEMContext & c,
              unsigned int i,
              const Point & p,
              Real /* time */,
              bool skip_context_check)
{
  LOG_SCOPE ("eval_at_point()", "OldSolutionCoefs");

  if (!skip_context_check)
    if (!this->check_old_context(c, p))
      return 0;

  // Get finite element object
  FEGenericBase<Real> * fe = nullptr;
  this->old_context.get_element_fe<Real>
    (i, fe, this->old_context.get_elem_dim());

  // Build a FE for calculating phi(p)
  FEGenericBase<Real> * fe_new =
    this->old_context.build_new_fe(fe, p);

  // Get the values and global indices of the shape functions
  const std::vector<std::vector<Real> > & phi = fe_new->get_phi();
  const std::vector<dof_id_type> & dof_indices =
    this->old_context.get_dof_indices(i);

  const std::size_t n_dofs = phi.size();
  libmesh_assert_equal_to(n_dofs, dof_indices.size());

  DynamicSparseNumberArray<Real, dof_id_type> returnval;
  returnval.resize(n_dofs);

  for (auto j : index_range(phi))
    {
      returnval.raw_at(j) = phi[j][0];
      returnval.raw_index(j) = dof_indices[j];
    }

  return returnval;
}



template<>
inline
VectorValue<DynamicSparseNumberArray<Real, dof_id_type> >
OldSolutionCoefs<RealGradient, &FEMContext::point_gradient>::
eval_at_point(const FEMContext & c,
              unsigned int i,
              const Point & p,
              Real /* time */,
              bool skip_context_check)
{
  LOG_SCOPE ("eval_at_point()", "OldSolutionCoefs");

  if (!skip_context_check)
    if (!this->check_old_context(c, p))
      return 0;

  // Get finite element object
  FEGenericBase<Real> * fe = nullptr;
  this->old_context.get_element_fe<Real>
    (i, fe, this->old_context.get_elem_dim());

  // Build a FE for calculating phi(p)
  FEGenericBase<Real> * fe_new =
    this->old_context.build_new_fe(fe, p);

  // Get the values and global indices of the shape functions
  const std::vector<std::vector<RealGradient> > & dphi = fe_new->get_dphi();
  const std::vector<dof_id_type> & dof_indices =
    this->old_context.get_dof_indices(i);

  const std::size_t n_dofs = dphi.size();
  libmesh_assert_equal_to(n_dofs, dof_indices.size());

  VectorValue<DynamicSparseNumberArray<Real, dof_id_type> > returnval;

  for (unsigned int d = 0; d != LIBMESH_DIM; ++d)
    returnval(d).resize(n_dofs);

  for (auto j : index_range(dphi))
    for (int d = 0; d != LIBMESH_DIM; ++d)
      {
        returnval(d).raw_at(j) = dphi[j][0](d);
        returnval(d).raw_index(j) = dof_indices[j];
      }

  return returnval;
}


template<>
inline
DynamicSparseNumberArray<Real, dof_id_type>
OldSolutionCoefs<Real, &FEMContext::point_value>::
eval_at_node(const FEMContext & c,
             unsigned int i,
             unsigned int /* elem_dim */,
             const Node & n,
             bool extra_hanging_dofs,
             Real /* time */)
{
  LOG_SCOPE ("Real eval_at_node()", "OldSolutionCoefs");

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
      n.old_dof_object->n_comp(sys.number(), i))
    {
      DynamicSparseNumberArray<Real, dof_id_type> returnval;
      const dof_id_type old_id =
        n.old_dof_object->dof_number(sys.number(), i, 0);
      returnval.resize(1);
      returnval.raw_at(0) = 1;
      returnval.raw_index(0) = old_id;
      return returnval;
    }

  return this->eval_at_point(c, i, n, 0, false);
}



template<>
inline
VectorValue<DynamicSparseNumberArray<Real, dof_id_type> >
OldSolutionCoefs<RealGradient, &FEMContext::point_gradient>::
eval_at_node(const FEMContext & c,
             unsigned int i,
             unsigned int elem_dim,
             const Node & n,
             bool extra_hanging_dofs,
             Real /* time */)
{
  LOG_SCOPE ("RealGradient eval_at_node()", "OldSolutionCoefs");

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
      n.old_dof_object->n_comp(sys.number(), i))
    {
      VectorValue<DynamicSparseNumberArray<Real, dof_id_type> > g;
      for (unsigned int d = 0; d != elem_dim; ++d)
        {
          const dof_id_type old_id =
            n.old_dof_object->dof_number(sys.number(), i, d+1);
          g(d).resize(1);
          g(d).raw_at(0) = 1;
          g(d).raw_index(0) = old_id;
        }
      return g;
    }

  return this->eval_at_point(c, i, n, 0, false);
}



/**
 * The MatrixFillAction output functor class can be used with
 * GenericProjector to write solution transfer coefficients into a
 * sparse matrix.
 *
 * \author Roy H. Stogner
 * \date 2017
 */
template <typename ValIn, typename ValOut>
class MatrixFillAction
{
public:
  typedef DynamicSparseNumberArray<ValIn, dof_id_type> InsertInput;
private:
  SparseMatrix<ValOut> & target_matrix;

public:
  MatrixFillAction(SparseMatrix<ValOut> & target_mat) :
    target_matrix(target_mat) {}

  void insert(dof_id_type id,
              const DynamicSparseNumberArray<ValIn, dof_id_type> & val)
  {
    // Lock the target matrix since it is shared among threads.
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

      const std::size_t dnsa_size = val.size();
      for (unsigned int j = 0; j != dnsa_size; ++j)
        {
          const dof_id_type dof_j = val.raw_index(j);
          const ValIn dof_val = val.raw_at(j);
          target_matrix.set(id, dof_j, dof_val);
        }
    }
  }


  void insert(const std::vector<dof_id_type> & dof_indices,
              const std::vector<DynamicSparseNumberArray<ValIn, dof_id_type> > & Ue)
  {
    const numeric_index_type
      begin_dof = target_matrix.row_start(),
      end_dof = target_matrix.row_stop();

    unsigned int size = Ue.size();

    libmesh_assert_equal_to(size, dof_indices.size());

    // Lock the target matrix since it is shared among threads.
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

      for (unsigned int i = 0; i != size; ++i)
        {
          const dof_id_type dof_i = dof_indices[i];
          if ((dof_i >= begin_dof) && (dof_i < end_dof))
            {
              const DynamicSparseNumberArray<ValIn,dof_id_type> & dnsa = Ue[i];
              const std::size_t dnsa_size = dnsa.size();
              for (unsigned int j = 0; j != dnsa_size; ++j)
                {
                  const dof_id_type dof_j = dnsa.raw_index(j);
                  const ValIn dof_val = dnsa.raw_at(j);
                  target_matrix.set(dof_i, dof_j, dof_val);
                }
            }
        }
    }
  }
};



/**
 * This method creates a projection matrix which corresponds to the
 * operation of project_vector between old and new solution spaces.
 */
void System::projection_matrix (SparseMatrix<Number> & proj_mat) const
{
  LOG_SCOPE ("projection_matrix()", "System");

  const unsigned int n_variables = this->n_vars();

  if (n_variables)
    {
      ConstElemRange active_local_elem_range
        (this->get_mesh().active_local_elements_begin(),
         this->get_mesh().active_local_elements_end());

      std::vector<unsigned int> vars(n_variables);
      std::iota(vars.begin(), vars.end(), 0);

      // Use a typedef to make the calling sequence for parallel_for() a bit more readable
      typedef OldSolutionCoefs<Real, &FEMContext::point_value> OldSolutionValueCoefs;
      typedef OldSolutionCoefs<RealGradient, &FEMContext::point_gradient> OldSolutionGradientCoefs;

      typedef
        GenericProjector<OldSolutionValueCoefs,
                         OldSolutionGradientCoefs,
                         DynamicSparseNumberArray<Real,dof_id_type>,
                         MatrixFillAction<Real, Number> > ProjMatFiller;

      OldSolutionValueCoefs    f(*this);
      OldSolutionGradientCoefs g(*this);
      MatrixFillAction<Real, Number> setter(proj_mat);

      ProjMatFiller mat_filler(*this, f, &g, setter, vars);
      mat_filler.project(active_local_elem_range);

      // Set the SCALAR dof transfer entries too.
      // Note: We assume that all SCALAR dofs are on the
      // processor with highest ID
      if (this->processor_id() == (this->n_processors()-1))
        {
          const DofMap & dof_map = this->get_dof_map();
          for (auto var : make_range(this->n_vars()))
            if (this->variable(var).type().family == SCALAR)
              {
                // We can just map SCALAR dofs directly across
                std::vector<dof_id_type> new_SCALAR_indices, old_SCALAR_indices;
                dof_map.SCALAR_dof_indices (new_SCALAR_indices, var, false);
                dof_map.SCALAR_dof_indices (old_SCALAR_indices, var, true);
                const unsigned int new_n_dofs =
                  cast_int<unsigned int>(new_SCALAR_indices.size());

                for (unsigned int i=0; i<new_n_dofs; i++)
                  {
                    proj_mat.set( new_SCALAR_indices[i],
                                  old_SCALAR_indices[i], 1);
                  }
              }
        }
    }
}
#endif // LIBMESH_HAVE_METAPHYSICL
#endif // LIBMESH_ENABLE_AMR



/**
 * This method projects an arbitrary function onto the solution via L2
 * projections and nodal interpolations on each element.
 */
void System::project_solution (ValueFunctionPointer fptr,
                               GradientFunctionPointer gptr,
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
void System::project_vector (ValueFunctionPointer fptr,
                             GradientFunctionPointer gptr,
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
    this->project_vector(new_vector, &f_fem, nullptr, is_adjoint);
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
  std::iota(vars.begin(), vars.end(), 0);

  // Use a typedef to make the calling sequence for parallel_for() a bit more readable
  typedef
    GenericProjector<FEMFunctionWrapper<Number>, FEMFunctionWrapper<Gradient>,
                     Number, VectorSetAction<Number>> FEMProjector;

  FEMFunctionWrapper<Number> fw(*f);

  if (g)
    {
      FEMFunctionWrapper<Gradient> gw(*g);

      FEMProjector projector(*this, fw, &gw, setter, vars);
      projector.project(active_local_range);
    }
  else
    {
      FEMProjector projector(*this, fw, nullptr, setter, vars);
      projector.project(active_local_range);
    }

  // Also, load values into the SCALAR dofs
  // Note: We assume that all SCALAR dofs are on the
  // processor with highest ID
  if (this->processor_id() == (this->n_processors()-1))
    {
      // FIXME: Do we want to first check for SCALAR vars before building this? [PB]
      FEMContext context( *this );

      const DofMap & dof_map = this->get_dof_map();
      for (auto var : make_range(this->n_vars()))
        if (this->variable(var).type().family == SCALAR)
          {
            // FIXME: We reinit with an arbitrary element in case the user
            //        doesn't override FEMFunctionBase::component. Is there
            //        any use case we're missing? [PB]
            context.pre_fe_reinit(*this, *(this->get_mesh().active_local_elements_begin()));

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
#else
  libmesh_ignore(is_adjoint);
#endif
}


/**
 * This method projects components of an arbitrary boundary function
 * onto the solution via L2 projections and nodal interpolations on
 * each element.
 */
void System::boundary_project_solution (const std::set<boundary_id_type> & b,
                                        const std::vector<unsigned int> & variables,
                                        ValueFunctionPointer fptr,
                                        GradientFunctionPointer gptr,
                                        const Parameters & parameters)
{
  WrappedFunction<Number> f(*this, fptr, &parameters);
  WrappedFunction<Gradient> g(*this, gptr, &parameters);
  this->boundary_project_solution(b, variables, &f, &g);
}


/**
 * This method projects an arbitrary boundary function onto the
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
                                      ValueFunctionPointer fptr,
                                      GradientFunctionPointer gptr,
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
#else
  libmesh_ignore(is_adjoint);
#endif
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
  for (const auto & elem : range)
    {
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

          for (auto & node : elem->node_ref_range())
            {
              const DofObject * old_dofs = node.old_dof_object;

              if (old_dofs)
                {
                  const unsigned int sysnum = system.number();
                  const unsigned int nvg = old_dofs->n_var_groups(sysnum);

                  for (unsigned int vg=0; vg != nvg; ++vg)
                    {
                      const unsigned int nvig =
                        old_dofs->n_vars(sysnum, vg);
                      for (unsigned int vig=0; vig != nvig; ++vig)
                        {
                          const unsigned int n_comp =
                            old_dofs->n_comp_group(sysnum, vg);
                          for (unsigned int c=0; c != n_comp; ++c)
                            {
                              const dof_id_type old_id =
                                old_dofs->dof_number(sysnum, vg, vig,
                                                     c, n_comp);

                              // We should either have no old id
                              // (e.g. on a newly expanded subdomain)
                              // or an id from the old system.
                              libmesh_assert(old_id < dof_map.n_old_dofs() ||
                                             old_id == DofObject::invalid_id);
                              di.push_back(old_id);
                            }
                        }
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
          for (auto & child : elem->child_ref_range())
            {
              dof_map.old_dof_indices (&child, di_child);
              di.insert(di.end(), di_child.begin(), di_child.end());
            }
        }
      else
        dof_map.old_dof_indices (elem, di);

      for (auto di_i : di)
        {
          // If we've just expanded a subdomain for a
          // subdomain-restricted variable, then we may have an
          // old_dof_object that doesn't have an old DoF for every
          // local index.
          if (di_i == DofObject::invalid_id)
            continue;

          libmesh_assert_less(di_i, dof_map.n_old_dofs());
          if (di_i < first_old_dof || di_i >= end_old_dof)
            this->send_list.push_back(di_i);
        }
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
  for (auto v : make_range(variables.size()))
    {
      const unsigned int var = variables[v];

      const Variable & variable = dof_map.variable(var);

      const FEType & fe_type = variable.type();

      if (fe_type.family == SCALAR)
        continue;

      const unsigned int var_component =
        system.variable_scalar_number(var, 0);

      // Get FE objects of the appropriate type
      std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));

      // Prepare variables for projection
      std::unique_ptr<QBase> qedgerule (fe_type.default_quadrature_rule(1));
      std::unique_ptr<QBase> qsiderule (fe_type.default_quadrature_rule(dim-1));

      // The values of the shape functions at the quadrature
      // points
      const std::vector<std::vector<Real>> & phi = fe->get_phi();

      // The gradients of the shape functions at the quadrature
      // points on the child element.
      const std::vector<std::vector<RealGradient>> * dphi = nullptr;

      const FEContinuity cont = fe->get_continuity();

      if (cont == C_ONE)
        {
          // We'll need gradient data for a C1 projection
          libmesh_assert(g.get());

          const std::vector<std::vector<RealGradient>> &
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
      for (const auto & elem : range)
        {
          // Per-subdomain variables don't need to be projected on
          // elements where they're not active
          if (!variable.active_on_subdomain(elem->subdomain_id()))
            continue;

          const unsigned short n_nodes = elem->n_nodes();
          const unsigned short n_edges = elem->n_edges();
          const unsigned short n_sides = elem->n_sides();

          // Find out which nodes, edges and sides are on a requested
          // boundary:
          std::vector<bool> is_boundary_node(n_nodes, false),
            is_boundary_edge(n_edges, false),
            is_boundary_side(n_sides, false);

          // We also maintain a separate list of nodeset-based boundary nodes
          std::vector<bool> is_boundary_nodeset(n_nodes, false);

          for (unsigned char s=0; s != n_sides; ++s)
            {
              // First see if this side has been requested
              boundary_info.boundary_ids (elem, s, bc_ids);
              bool do_this_side = false;
              for (const auto & bc_id : bc_ids)
                if (b.count(bc_id))
                  {
                    do_this_side = true;
                    break;
                  }
              if (!do_this_side)
                continue;

              is_boundary_side[s] = true;

              // Then see what nodes and what edges are on it
              for (unsigned int n=0; n != n_nodes; ++n)
                if (elem->is_node_on_side(n,s))
                  is_boundary_node[n] = true;
              for (unsigned int e=0; e != n_edges; ++e)
                if (elem->is_edge_on_side(e,s))
                  is_boundary_edge[e] = true;
            }

            // We can also project on nodes, so we should also independently
            // check whether the nodes have been requested
            for (unsigned int n=0; n != n_nodes; ++n)
              {
                boundary_info.boundary_ids (elem->node_ptr(n), bc_ids);

                for (const auto & bc_id : bc_ids)
                  if (b.count(bc_id))
                    {
                      is_boundary_node[n] = true;
                      is_boundary_nodeset[n] = true;
                    }
              }

            // We can also project on edges, so we should also independently
            // check whether the edges have been requested
            for (unsigned short e=0; e != n_edges; ++e)
              {
                boundary_info.edge_boundary_ids (elem, e, bc_ids);

                for (const auto & bc_id : bc_ids)
                  if (b.count(bc_id))
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

          // Zero the interpolated values
          Ue.resize (n_dofs); Ue.zero();

          // In general, we need a series of
          // projections to ensure a unique and continuous
          // solution.  We start by interpolating boundary nodes, then
          // hold those fixed and project boundary edges, then hold
          // those fixed and project boundary faces,

          // Interpolate node values first
          unsigned int current_dof = 0;
          for (unsigned short n = 0; n != n_nodes; ++n)
            {
              // FIXME: this should go through the DofMap,
              // not duplicate dof_indices code badly!

              // This call takes into account elem->p_level() internally.
              const unsigned int nc =
                FEInterface::n_dofs_at_node (fe_type, elem, n);

              if ((!elem->is_vertex(n) || !is_boundary_node[n]) &&
                  !is_boundary_nodeset[n])
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
#if LIBMESH_DIM > 1
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

#if LIBMESH_DIM > 2
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
#endif // LIBMESH_DIM > 2
                    }
#endif // LIBMESH_DIM > 1
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
            for (unsigned short e = 0; e != n_edges; ++e)
              {
                if (!is_boundary_edge[e])
                  continue;

                FEInterface::dofs_on_edge(elem, dim, fe_type, e,
                                          side_dofs);

                const unsigned int n_side_dofs =
                  cast_int<unsigned int>(side_dofs.size());

                // Some edge dofs are on nodes and already
                // fixed, others are free to calculate
                unsigned int free_dofs = 0;
                for (auto i : make_range(n_side_dofs))
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
                    for (unsigned int sidei=0, freei=0;
                         sidei != n_side_dofs; ++sidei)
                      {
                        unsigned int i = side_dofs[sidei];
                        // fixed DoFs aren't test functions
                        if (dof_is_fixed[i])
                          continue;
                        for (unsigned int sidej=0, freej=0;
                             sidej != n_side_dofs; ++sidej)
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
            for (unsigned short s = 0; s != n_sides; ++s)
              {
                if (!is_boundary_side[s])
                  continue;

                FEInterface::dofs_on_side(elem, dim, fe_type, s,
                                          side_dofs);

                // Some side dofs are on nodes/edges and already
                // fixed, others are free to calculate
                unsigned int free_dofs = 0;
                for (auto i : index_range(side_dofs))
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

                const unsigned int n_side_dofs =
                  cast_int<unsigned int>(side_dofs.size());

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
                    for (unsigned int sidei=0, freei=0;
                         sidei != n_side_dofs; ++sidei)
                      {
                        unsigned int i = side_dofs[sidei];
                        // fixed DoFs aren't test functions
                        if (dof_is_fixed[i])
                          continue;
                        for (unsigned int sidej=0, freej=0;
                             sidej != n_side_dofs; ++sidej)
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
