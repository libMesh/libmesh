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

// Local includes
#include "libmesh/dof_map.h"

// libMesh includes
#include "libmesh/boundary_info.h" // needed for dirichlet constraints
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/elem.h"
#include "libmesh/elem_range.h"
#include "libmesh/fe_base.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_type.h"
#include "libmesh/function_base.h"
#include "libmesh/int_range.h"
#include "libmesh/libmesh_logging.h"
#include "libmesh/mesh_base.h"
#include "libmesh/mesh_inserter_iterator.h"
#include "libmesh/mesh_tools.h" // for libmesh_assert_valid_boundary_ids()
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/numeric_vector.h" // for enforce_constraints_exactly()
#include "libmesh/parallel_algebra.h"
#include "libmesh/parallel_elem.h"
#include "libmesh/parallel_node.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/periodic_boundary_base.h"
#include "libmesh/point_locator_base.h"
#include "libmesh/quadrature.h" // for dirichlet constraints
#include "libmesh/raw_accessor.h"
#include "libmesh/sparse_matrix.h" // needed to constrain adjoint rhs
#include "libmesh/system.h" // needed by enforce_constraints_exactly()
#include "libmesh/tensor_tools.h"
#include "libmesh/threads.h"

// TIMPI includes
#include "timpi/parallel_implementation.h"
#include "timpi/parallel_sync.h"

// C++ Includes
#include <set>
#include <algorithm> // for std::count, std::fill
#include <sstream>
#include <cstdlib> // *must* precede <cmath> for proper std:abs() on PGI, Sun Studio CC
#include <cmath>
#include <numeric>

// Anonymous namespace to hold helper classes
namespace {

using namespace libMesh;

class ComputeConstraints
{
public:
  ComputeConstraints (DofConstraints & constraints,
                      DofMap & dof_map,
#ifdef LIBMESH_ENABLE_PERIODIC
                      PeriodicBoundaries & periodic_boundaries,
#endif
                      const MeshBase & mesh,
                      const unsigned int variable_number) :
    _constraints(constraints),
    _dof_map(dof_map),
#ifdef LIBMESH_ENABLE_PERIODIC
    _periodic_boundaries(periodic_boundaries),
#endif
    _mesh(mesh),
    _variable_number(variable_number)
  {}

  void operator()(const ConstElemRange & range) const
  {
    const Variable & var_description = _dof_map.variable(_variable_number);

#ifdef LIBMESH_ENABLE_PERIODIC
    std::unique_ptr<PointLocatorBase> point_locator;
    const bool have_periodic_boundaries =
      !_periodic_boundaries.empty();
    if (have_periodic_boundaries && !range.empty())
      point_locator = _mesh.sub_point_locator();
#endif

    for (const auto & elem : range)
      if (var_description.active_on_subdomain(elem->subdomain_id()))
        {
#ifdef LIBMESH_ENABLE_AMR
          FEInterface::compute_constraints (_constraints,
                                            _dof_map,
                                            _variable_number,
                                            elem);
#endif
#ifdef LIBMESH_ENABLE_PERIODIC
          // FIXME: periodic constraints won't work on a non-serial
          // mesh unless it's kept ghost elements from opposing
          // boundaries!
          if (have_periodic_boundaries)
            FEInterface::compute_periodic_constraints (_constraints,
                                                       _dof_map,
                                                       _periodic_boundaries,
                                                       _mesh,
                                                       point_locator.get(),
                                                       _variable_number,
                                                       elem);
#endif
        }
  }

private:
  DofConstraints & _constraints;
  DofMap & _dof_map;
#ifdef LIBMESH_ENABLE_PERIODIC
  PeriodicBoundaries & _periodic_boundaries;
#endif
  const MeshBase & _mesh;
  const unsigned int _variable_number;
};



#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
class ComputeNodeConstraints
{
public:
  ComputeNodeConstraints (NodeConstraints & node_constraints,
#ifdef LIBMESH_ENABLE_PERIODIC
                          PeriodicBoundaries & periodic_boundaries,
#endif
                          const MeshBase & mesh) :
    _node_constraints(node_constraints),
#ifdef LIBMESH_ENABLE_PERIODIC
    _periodic_boundaries(periodic_boundaries),
#endif
    _mesh(mesh)
  {}

  void operator()(const ConstElemRange & range) const
  {
#ifdef LIBMESH_ENABLE_PERIODIC
    std::unique_ptr<PointLocatorBase> point_locator;
    bool have_periodic_boundaries = !_periodic_boundaries.empty();
    if (have_periodic_boundaries && !range.empty())
      point_locator = _mesh.sub_point_locator();
#endif

    for (const auto & elem : range)
      {
#ifdef LIBMESH_ENABLE_AMR
        FEBase::compute_node_constraints (_node_constraints, elem);
#endif
#ifdef LIBMESH_ENABLE_PERIODIC
        // FIXME: periodic constraints won't work on a non-serial
        // mesh unless it's kept ghost elements from opposing
        // boundaries!
        if (have_periodic_boundaries)
          FEBase::compute_periodic_node_constraints (_node_constraints,
                                                     _periodic_boundaries,
                                                     _mesh,
                                                     point_locator.get(),
                                                     elem);
#endif
      }
  }

private:
  NodeConstraints & _node_constraints;
#ifdef LIBMESH_ENABLE_PERIODIC
  PeriodicBoundaries & _periodic_boundaries;
#endif
  const MeshBase & _mesh;
};
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS


#ifdef LIBMESH_ENABLE_DIRICHLET

/**
 * This functor class hierarchy adds a constraint row to a DofMap
 */
class AddConstraint
{
protected:
  DofMap                  & dof_map;

public:
  AddConstraint(DofMap & dof_map_in) : dof_map(dof_map_in) {}

  virtual void operator()(dof_id_type dof_number,
                          const DofConstraintRow & constraint_row,
                          const Number constraint_rhs) const = 0;
};

class AddPrimalConstraint : public AddConstraint
{
public:
  AddPrimalConstraint(DofMap & dof_map_in) : AddConstraint(dof_map_in) {}

  virtual void operator()(dof_id_type dof_number,
                          const DofConstraintRow & constraint_row,
                          const Number constraint_rhs) const
  {
    if (!dof_map.is_constrained_dof(dof_number))
      dof_map.add_constraint_row (dof_number, constraint_row,
                                  constraint_rhs, true);
  }
};

class AddAdjointConstraint : public AddConstraint
{
private:
  const unsigned int qoi_index;

public:
  AddAdjointConstraint(DofMap & dof_map_in, unsigned int qoi_index_in)
    : AddConstraint(dof_map_in), qoi_index(qoi_index_in) {}

  virtual void operator()(dof_id_type dof_number,
                          const DofConstraintRow & constraint_row,
                          const Number constraint_rhs) const
  {
    dof_map.add_adjoint_constraint_row
      (qoi_index, dof_number, constraint_row, constraint_rhs,
       true);
  }
};



/**
 * This class implements turning an arbitrary
 * boundary function into Dirichlet constraints.  It
 * may be executed in parallel on multiple threads.
 */
class ConstrainDirichlet
{
private:
  DofMap                  & dof_map;
  const MeshBase          & mesh;
  const Real               time;
  const DirichletBoundaries & dirichlets;

  const AddConstraint     & add_fn;

  static Number f_component (FunctionBase<Number> * f,
                             FEMFunctionBase<Number> * f_fem,
                             const FEMContext * c,
                             unsigned int i,
                             const Point & p,
                             Real time)
  {
    if (f_fem)
      {
        if (c)
          return f_fem->component(*c, i, p, time);
        else
          return std::numeric_limits<Real>::quiet_NaN();
      }
    return f->component(i, p, time);
  }

  static Gradient g_component (FunctionBase<Gradient> * g,
                               FEMFunctionBase<Gradient> * g_fem,
                               const FEMContext * c,
                               unsigned int i,
                               const Point & p,
                               Real time)
  {
    if (g_fem)
      {
        if (c)
          return g_fem->component(*c, i, p, time);
        else
          return std::numeric_limits<Number>::quiet_NaN();
      }
    return g->component(i, p, time);
  }



  /**
   * Handy struct to pass around BoundaryInfo for a single Elem.  Must
   * be created with a reference to a BoundaryInfo object and a map
   * from boundary_id -> set<DirichletBoundary *> objects involving
   * that id.
   */
  struct SingleElemBoundaryInfo
  {
    SingleElemBoundaryInfo(const BoundaryInfo & bi,
                           const std::map<boundary_id_type, std::set<std::pair<unsigned int, DirichletBoundary *>>> & ordered_map_in) :
      boundary_info(bi),
      boundary_id_to_ordered_dirichlet_boundaries(ordered_map_in),
      elem(nullptr),
      n_sides(0),
      n_edges(0),
      n_nodes(0)
    {}

    const BoundaryInfo & boundary_info;
    const std::map<boundary_id_type, std::set<std::pair<unsigned int, DirichletBoundary *>>> & boundary_id_to_ordered_dirichlet_boundaries;
    const Elem * elem;

    unsigned short n_sides;
    unsigned short n_edges;
    unsigned short n_nodes;

    // Mapping from DirichletBoundary objects which are active on this
    // element to sides/nodes/edges/shellfaces of this element which
    // they are active on.
    std::map<const DirichletBoundary *, std::vector<bool>> is_boundary_node_map;
    std::map<const DirichletBoundary *, std::vector<bool>> is_boundary_side_map;
    std::map<const DirichletBoundary *, std::vector<bool>> is_boundary_edge_map;
    std::map<const DirichletBoundary *, std::vector<bool>> is_boundary_shellface_map;

    std::map<const DirichletBoundary *, std::vector<bool>> is_boundary_nodeset_map;

    // The set of (dirichlet_id, DirichletBoundary) pairs which have at least one boundary
    // id related to this Elem.
    std::set<std::pair<unsigned int, DirichletBoundary *>> ordered_dbs;

    /**
     * Given a single Elem, fills the SingleElemBoundaryInfo struct with
     * required data.
     *
     * @return true if this Elem has _any_ boundary ids associated with
     * it, false otherwise.
     */
    bool reinit(const Elem * elem_in)
    {
      elem = elem_in;

      n_sides = elem->n_sides();
      n_edges = elem->n_edges();
      n_nodes = elem->n_nodes();

      // objects and node/side/edge/shellface ids.
      is_boundary_node_map.clear();
      is_boundary_side_map.clear();
      is_boundary_edge_map.clear();
      is_boundary_shellface_map.clear();
      is_boundary_nodeset_map.clear();

      // Clear any DirichletBoundaries from the previous Elem
      ordered_dbs.clear();

      // Update has_dirichlet_constraint below, and if it remains false then
      // we can skip this element since there are not constraints to impose.
      bool has_dirichlet_constraint = false;

      // Container to catch boundary ids handed back for sides,
      // nodes, and edges in the loops below.
      std::vector<boundary_id_type> ids_vec;

      for (unsigned char s = 0; s != n_sides; ++s)
        {
          // First see if this side has been requested
          boundary_info.boundary_ids (elem, s, ids_vec);

          bool do_this_side = false;
          for (const auto & bc_id : ids_vec)
            {
              auto it = boundary_id_to_ordered_dirichlet_boundaries.find(bc_id);
              if (it != boundary_id_to_ordered_dirichlet_boundaries.end())
                {
                  do_this_side = true;

                  // Associate every DirichletBoundary object that has this bc_id with the current Elem
                  ordered_dbs.insert(it->second.begin(), it->second.end());

                  // Turn on the flag for the current side for each DirichletBoundary
                  for (const auto & db_pair : it->second)
                    {
                      // Attempt to emplace an empty vector. If there
                      // is already an entry, the insertion will fail
                      // and we'll get an iterator back to the
                      // existing entry. Either way, we'll then set
                      // index s of that vector to "true".
                      auto pr = is_boundary_side_map.emplace(db_pair.second, std::vector<bool>(n_sides, false));
                      pr.first->second[s] = true;
                    }
                }
            }

          if (!do_this_side)
            continue;

          has_dirichlet_constraint = true;

          // Then determine what nodes are on this side
          for (unsigned int n = 0; n != n_nodes; ++n)
            if (elem->is_node_on_side(n,s))
              {
                // Attempt to emplace an empty vector. If there is
                // already an entry, the insertion will fail and we'll
                // get an iterator back to the existing entry. Either
                // way, we'll then set index n of that vector to
                // "true".
                for (const auto & db_pair : ordered_dbs)
                  {
                    // Only add this as a boundary node for this db if
                    // it is also a boundary side for this db.
                    auto side_it = is_boundary_side_map.find(db_pair.second);
                    if (side_it != is_boundary_side_map.end() && side_it->second[s])
                      {
                        auto pr = is_boundary_node_map.emplace(db_pair.second, std::vector<bool>(n_nodes, false));
                        pr.first->second[n] = true;
                      }
                  }
              }

          // Finally determine what edges are on this side
          for (unsigned int e = 0; e != n_edges; ++e)
            if (elem->is_edge_on_side(e,s))
              {
                // Attempt to emplace an empty vector. If there is
                // already an entry, the insertion will fail and we'll
                // get an iterator back to the existing entry. Either
                // way, we'll then set index e of that vector to
                // "true".
                for (const auto & db_pair : ordered_dbs)
                  {
                    // Only add this as a boundary edge for this db if
                    // it is also a boundary side for this db.
                    auto side_it = is_boundary_side_map.find(db_pair.second);
                    if (side_it != is_boundary_side_map.end() && side_it->second[s])
                      {
                        auto pr = is_boundary_edge_map.emplace(db_pair.second, std::vector<bool>(n_edges, false));
                        pr.first->second[e] = true;
                      }
                  }
              }
        } // for (s = 0..n_sides)

      // We can also impose Dirichlet boundary conditions on nodes, so we should
      // also independently check whether the nodes have been requested
      for (unsigned int n=0; n != n_nodes; ++n)
        {
          boundary_info.boundary_ids (elem->node_ptr(n), ids_vec);

          for (const auto & bc_id : ids_vec)
            {
              auto it = boundary_id_to_ordered_dirichlet_boundaries.find(bc_id);
              if (it != boundary_id_to_ordered_dirichlet_boundaries.end())
                {
                  // Associate every DirichletBoundary object that has this bc_id with the current Elem
                  ordered_dbs.insert(it->second.begin(), it->second.end());

                  // Turn on the flag for the current node for each DirichletBoundary
                  for (const auto & db_pair : it->second)
                    {
                      auto pr = is_boundary_node_map.emplace(db_pair.second, std::vector<bool>(n_nodes, false));
                      pr.first->second[n] = true;

                      auto pr2 = is_boundary_nodeset_map.emplace(db_pair.second, std::vector<bool>(n_nodes, false));
                      pr2.first->second[n] = true;
                    }

                  has_dirichlet_constraint = true;
                }
            }
        } // for (n = 0..n_nodes)

      // We can also impose Dirichlet boundary conditions on edges, so we should
      // also independently check whether the edges have been requested
      for (unsigned short e=0; e != n_edges; ++e)
        {
          boundary_info.edge_boundary_ids (elem, e, ids_vec);

          bool do_this_side = false;
          for (const auto & bc_id : ids_vec)
            {
              auto it = boundary_id_to_ordered_dirichlet_boundaries.find(bc_id);
              if (it != boundary_id_to_ordered_dirichlet_boundaries.end())
                {
                  do_this_side = true;

                  // We need to loop over all DirichletBoundary objects associated with bc_id
                  ordered_dbs.insert(it->second.begin(), it->second.end());

                  // Turn on the flag for the current edge for each DirichletBoundary
                  for (const auto & db_pair : it->second)
                    {
                      auto pr = is_boundary_edge_map.emplace(db_pair.second, std::vector<bool>(n_edges, false));
                      pr.first->second[e] = true;
                    }
                }
            }

          if (!do_this_side)
            continue;

          has_dirichlet_constraint = true;

          // Then determine what nodes are on this edge
          for (unsigned int n = 0; n != n_nodes; ++n)
            if (elem->is_node_on_edge(n,e))
              {
                // Attempt to emplace an empty vector. If there is
                // already an entry, the insertion will fail and we'll
                // get an iterator back to the existing entry. Either
                // way, we'll then set index n of that vector to
                // "true".
                for (const auto & db_pair : ordered_dbs)
                  {
                    // Only add this as a boundary node for this db if
                    // it is also a boundary edge for this db.
                    auto edge_it = is_boundary_edge_map.find(db_pair.second);
                    if (edge_it != is_boundary_edge_map.end() && edge_it->second[e])
                      {
                        auto pr = is_boundary_node_map.emplace(db_pair.second, std::vector<bool>(n_nodes, false));
                        pr.first->second[n] = true;
                      }
                  }
              }
        }

      // We can also impose Dirichlet boundary conditions on shellfaces, so we should
      // also independently check whether the shellfaces have been requested
      for (unsigned short shellface=0; shellface != 2; ++shellface)
        {
          boundary_info.shellface_boundary_ids (elem, shellface, ids_vec);
          bool do_this_shellface = false;

          for (const auto & bc_id : ids_vec)
            {
              auto it = boundary_id_to_ordered_dirichlet_boundaries.find(bc_id);
              if (it != boundary_id_to_ordered_dirichlet_boundaries.end())
                {
                  has_dirichlet_constraint = true;
                  do_this_shellface = true;

                  // We need to loop over all DirichletBoundary objects associated with bc_id
                  ordered_dbs.insert(it->second.begin(), it->second.end());

                  // Turn on the flag for the current shellface for each DirichletBoundary
                  for (const auto & db_pair : it->second)
                    {
                      auto pr = is_boundary_shellface_map.emplace(db_pair.second, std::vector<bool>(/*n_shellfaces=*/2, false));
                      pr.first->second[shellface] = true;
                    }
                }
            }

          if (do_this_shellface)
            {
              // Shellface BCs induce BCs on all the nodes of a shell Elem
              for (unsigned int n = 0; n != n_nodes; ++n)
                for (const auto & db_pair : ordered_dbs)
                  {
                    // Only add this as a boundary node for this db if
                    // it is also a boundary shellface for this db.
                    auto side_it = is_boundary_shellface_map.find(db_pair.second);
                    if (side_it != is_boundary_shellface_map.end() && side_it->second[shellface])
                      {
                        auto pr = is_boundary_node_map.emplace(db_pair.second, std::vector<bool>(n_nodes, false));
                        pr.first->second[n] = true;
                      }
                  }
            }
        } // for (shellface = 0..2)

      return has_dirichlet_constraint;
    } // SingleElemBoundaryInfo::reinit()

  }; // struct SingleElemBoundaryInfo



  template<typename OutputType>
  void apply_lagrange_dirichlet_impl(const SingleElemBoundaryInfo & sebi,
                            const Variable & variable,
                            const DirichletBoundary & dirichlet,
                            FEMContext & fem_context) const
  {
    // Get pointer to the Elem we are currently working on
    const Elem * elem = sebi.elem;

    // Per-subdomain variables don't need to be projected on
    // elements where they're not active
    if (!variable.active_on_subdomain(elem->subdomain_id()))
      return;

    FunctionBase<Number> * f = dirichlet.f.get();
    FEMFunctionBase<Number> * f_fem = dirichlet.f_fem.get();

    const System * f_system = dirichlet.f_system;

    // We need data to project
    libmesh_assert(f || f_fem);
    libmesh_assert(!(f && f_fem));

    // Iff our data depends on a system, we should have it.
    libmesh_assert(!(f && f_system));
    libmesh_assert(!(f_fem && !f_system));

    // The new element coefficients. For Lagrange FEs, these are the
    // nodal values.
    DenseVector<Number> Ue;

    // Get a reference to the fe_type associated with this variable
    const FEType & fe_type = variable.type();

    // Dimension of the vector-valued FE (1 for scalar-valued FEs)
    unsigned int n_vec_dim = FEInterface::n_vec_dim(mesh, fe_type);

    const unsigned int var_component =
      variable.first_scalar_number();

    // Get this Variable's number, as determined by the System.
    const unsigned int var = variable.number();

    // If our supplied functions require a FEMContext, and if we have
    // an initialized solution to use with that FEMContext, then
    // create one
    std::unique_ptr<FEMContext> context;
    if (f_fem)
      {
        libmesh_assert(f_system);
        if (f_system->current_local_solution->initialized())
          {
            context = libmesh_make_unique<FEMContext>(*f_system);
            f_fem->init_context(*context);
          }
      }

    if (f_system && context.get())
      context->pre_fe_reinit(*f_system, elem);

    // Also pre-init the fem_context() we were passed on the current Elem.
    fem_context.pre_fe_reinit(fem_context.get_system(), elem);

    // Get a reference to the DOF indices for the current element
    const std::vector<dof_id_type> & dof_indices =
      fem_context.get_dof_indices(var);

    // The number of DOFs on the element
    const unsigned int n_dofs =
      cast_int<unsigned int>(dof_indices.size());

    // Fixed vs. free DoFs on edge/face projections
    std::vector<char> dof_is_fixed(n_dofs, false); // bools

    // Zero the interpolated values
    Ue.resize (n_dofs); Ue.zero();

    // For Lagrange elements, side, edge, and shellface BCs all
    // "induce" boundary conditions on the nodes of those entities.
    // In SingleElemBoundaryInfo::reinit(), we therefore set entries
    // in the "is_boundary_node_map" container based on side and
    // shellface BCs, Then, when we actually apply constraints, we
    // only have to check whether any Nodes are in this container, and
    // compute values as necessary.
    unsigned int current_dof = 0;
    for (unsigned int n=0; n!= sebi.n_nodes; ++n)
      {
        // For Lagrange this can return 0 (in case of a lower-order FE
        // on a higher-order Elem) or 1. This function accounts for the
        // Elem::p_level() internally.
        const unsigned int nc =
          FEInterface::n_dofs_at_node (fe_type, elem, n);

        // If there are no DOFs at this node, then it doesn't matter
        // if it's technically a boundary node or not, there's nothing
        // to constrain.
        if (!nc)
          continue;

        // Check whether the current node is a boundary node
        auto is_boundary_node_it = sebi.is_boundary_node_map.find(&dirichlet);
        const bool is_boundary_node =
          (is_boundary_node_it != sebi.is_boundary_node_map.end() &&
           is_boundary_node_it->second[n]);

        // Check whether the current node is in a boundary nodeset
        auto is_boundary_nodeset_it = sebi.is_boundary_nodeset_map.find(&dirichlet);
        const bool is_boundary_nodeset =
          (is_boundary_nodeset_it != sebi.is_boundary_nodeset_map.end() &&
           is_boundary_nodeset_it->second[n]);

        // If node is neither a boundary node or from a boundary nodeset, go to the next one.
        if ( !(is_boundary_node || is_boundary_nodeset) )
          {
            current_dof += nc;
            continue;
          }

        // Compute function values, storing them in Ue
        libmesh_assert_equal_to (nc, n_vec_dim);
        for (unsigned int c = 0; c < n_vec_dim; c++)
          {
            Ue(current_dof+c) =
              f_component(f, f_fem, context.get(), var_component+c,
                          elem->point(n), time);
            dof_is_fixed[current_dof+c] = true;
          }
        current_dof += n_vec_dim;
      } // end for (n=0..n_nodes)

    // Lock the DofConstraints since it is shared among threads.
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

      for (unsigned int i = 0; i < n_dofs; i++)
        {
          DofConstraintRow empty_row;
          if (dof_is_fixed[i] && !libmesh_isnan(Ue(i)))
            add_fn (dof_indices[i], empty_row, Ue(i));
        }
    }

  } // apply_lagrange_dirichlet_impl



  template<typename OutputType>
  void apply_dirichlet_impl(const SingleElemBoundaryInfo & sebi,
                            const Variable & variable,
                            const DirichletBoundary & dirichlet,
                            FEMContext & fem_context) const
  {
    // Get pointer to the Elem we are currently working on
    const Elem * elem = sebi.elem;

    // Per-subdomain variables don't need to be projected on
    // elements where they're not active
    if (!variable.active_on_subdomain(elem->subdomain_id()))
      return;

    typedef OutputType                                                      OutputShape;
    typedef typename TensorTools::IncrementRank<OutputShape>::type          OutputGradient;
    //typedef typename TensorTools::IncrementRank<OutputGradient>::type       OutputTensor;
    typedef typename TensorTools::MakeNumber<OutputShape>::type             OutputNumber;
    typedef typename TensorTools::IncrementRank<OutputNumber>::type         OutputNumberGradient;
    //typedef typename TensorTools::IncrementRank<OutputNumberGradient>::type OutputNumberTensor;

    FunctionBase<Number> * f = dirichlet.f.get();
    FunctionBase<Gradient> * g = dirichlet.g.get();

    FEMFunctionBase<Number> * f_fem = dirichlet.f_fem.get();
    FEMFunctionBase<Gradient> * g_fem = dirichlet.g_fem.get();

    const System * f_system = dirichlet.f_system;

    // We need data to project
    libmesh_assert(f || f_fem);
    libmesh_assert(!(f && f_fem));

    // Iff our data depends on a system, we should have it.
    libmesh_assert(!(f && f_system));
    libmesh_assert(!(f_fem && !f_system));

    // The element matrix and RHS for projections.
    // Note that Ke is always real-valued, whereas
    // Fe may be complex valued if complex number
    // support is enabled
    DenseMatrix<Real> Ke;
    DenseVector<Number> Fe;
    // The new element coefficients
    DenseVector<Number> Ue;

    // The dimensionality of the current mesh
    const unsigned int dim = mesh.mesh_dimension();

    // Get a reference to the fe_type associated with this variable
    const FEType & fe_type = variable.type();

    // Dimension of the vector-valued FE (1 for scalar-valued FEs)
    unsigned int n_vec_dim = FEInterface::n_vec_dim(mesh, fe_type);

    const unsigned int var_component =
      variable.first_scalar_number();

    // Get this Variable's number, as determined by the System.
    const unsigned int var = variable.number();

    // The type of projections done depend on the FE's continuity.
    FEContinuity cont = FEInterface::get_continuity(fe_type);

    // Make sure we have the right data available for C1 projections
    if ((cont == C_ONE) && (fe_type.family != SUBDIVISION))
      {
        // We'll need gradient data for a C1 projection
        libmesh_assert(g || g_fem);

        // We currently demand that either neither nor both function
        // object depend on current FEM data.
        libmesh_assert(!(g && g_fem));
        libmesh_assert(!(f && g_fem));
        libmesh_assert(!(f_fem && g));
      }

    // If our supplied functions require a FEMContext, and if we have
    // an initialized solution to use with that FEMContext, then
    // create one
    std::unique_ptr<FEMContext> context;
    if (f_fem)
      {
        libmesh_assert(f_system);
        if (f_system->current_local_solution->initialized())
          {
            context = libmesh_make_unique<FEMContext>(*f_system);
            f_fem->init_context(*context);
            if (g_fem)
              g_fem->init_context(*context);
          }
      }

    // There's a chicken-and-egg problem with FEMFunction-based
    // Dirichlet constraints: we can't evaluate the FEMFunction
    // until we have an initialized local solution vector, we
    // can't initialize local solution vectors until we have a
    // send list, and we can't generate a send list until we know
    // all our constraints
    //
    // We don't generate constraints on uninitialized systems;
    // currently user code will have to reinit() before any
    // FEMFunction-based constraints will be correct.  This should
    // be fine, since user code would want to reinit() after
    // setting initial conditions anyway.
    if (f_system && context.get())
      context->pre_fe_reinit(*f_system, elem);

    // Also pre-init the fem_context() we were passed on the current Elem.
    fem_context.pre_fe_reinit(fem_context.get_system(), elem);

    // Get a reference to the DOF indices for the current element
    const std::vector<dof_id_type> & dof_indices =
      fem_context.get_dof_indices(var);

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

    // Interpolate node values first. Note that we have a special
    // case for nodes that have a boundary nodeset, since we do
    // need to interpolate them directly, even if they're non-vertex
    // nodes.
    unsigned int current_dof = 0;
    for (unsigned int n=0; n!= sebi.n_nodes; ++n)
      {
        // FIXME: this should go through the DofMap,
        // not duplicate dof_indices code badly!

        // Get the number of DOFs at this node, accounting for
        // Elem::p_level() internally.
        const unsigned int nc =
          FEInterface::n_dofs_at_node (fe_type, elem, n);

        // Get a reference to the "is_boundary_node" flags for the
        // current DirichletBoundary object.  In case the map does not
        // contain an entry for this DirichletBoundary object, it
        // means there are no boundary nodes active.
        auto is_boundary_node_it = sebi.is_boundary_node_map.find(&dirichlet);

        // The current n is not a boundary node if either there is no
        // boundary_node_map for this DirichletBoundary object, or if
        // there is but the entry in the corresponding vector is
        // false.
        const bool not_boundary_node =
          (is_boundary_node_it == sebi.is_boundary_node_map.end() ||
           !is_boundary_node_it->second[n]);

        // Same thing for nodesets
        auto is_boundary_nodeset_it = sebi.is_boundary_nodeset_map.find(&dirichlet);
        const bool not_boundary_nodeset =
          (is_boundary_nodeset_it == sebi.is_boundary_nodeset_map.end() ||
           !is_boundary_nodeset_it->second[n]);

        if ((!elem->is_vertex(n) || not_boundary_node) &&
            not_boundary_nodeset)
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
        else if ((cont == C_ZERO) || (fe_type.family == SUBDIVISION))
          {
            libmesh_assert_equal_to (nc, n_vec_dim);
            for (unsigned int c = 0; c < n_vec_dim; c++)
              {
                Ue(current_dof+c) =
                  f_component(f, f_fem, context.get(), var_component+c,
                              elem->point(n), time);
                dof_is_fixed[current_dof+c] = true;
              }
            current_dof += n_vec_dim;
          }
        // The hermite element vertex shape functions are weird
        else if (fe_type.family == HERMITE)
          {
            Ue(current_dof) =
              f_component(f, f_fem, context.get(), var_component,
                          elem->point(n), time);
            dof_is_fixed[current_dof] = true;
            current_dof++;
            Gradient grad =
              g_component(g, g_fem, context.get(), var_component,
                          elem->point(n), time);
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
                Gradient gxminus =
                  g_component(g, g_fem, context.get(), var_component,
                              nxminus, time);
                Gradient gxplus =
                  g_component(g, g_fem, context.get(), var_component,
                              nxplus, time);
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
                    Gradient gyminus =
                      g_component(g, g_fem, context.get(), var_component,
                                  nyminus, time);
                    Gradient gyplus =
                      g_component(g, g_fem, context.get(), var_component,
                                  nyplus, time);
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
                    Gradient gxmym =
                      g_component(g, g_fem, context.get(), var_component,
                                  nxmym, time);
                    Gradient gxmyp =
                      g_component(g, g_fem, context.get(), var_component,
                                  nxmyp, time);
                    Gradient gxpym =
                      g_component(g, g_fem, context.get(), var_component,
                                  nxpym, time);
                    Gradient gxpyp =
                      g_component(g, g_fem, context.get(), var_component,
                                  nxpyp, time);
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
            Ue(current_dof) =
              f_component(f, f_fem, context.get(), var_component,
                          elem->point(n), time);
            dof_is_fixed[current_dof] = true;
            current_dof++;
            Gradient grad =
              g_component(g, g_fem, context.get(), var_component,
                          elem->point(n), time);
            for (unsigned int i=0; i!= dim; ++i)
              {
                Ue(current_dof) = grad(i);
                dof_is_fixed[current_dof] = true;
                current_dof++;
              }
          }
        else
          libmesh_error_msg("Unknown continuity cont = " << cont);
      } // end for (n=0..n_nodes)

        // In 3D, project any edge values next
    if (dim > 2 && cont != DISCONTINUOUS)
      {
        // Get a pointer to the 1 dimensional (edge) FE for the current
        // var which is stored in the fem_context. This will only be
        // different from side_fe in 3D.
        FEGenericBase<OutputType> * edge_fe = nullptr;
        fem_context.get_edge_fe(var, edge_fe);

        // Set tolerance on underlying FEMap object. This will allow us to
        // avoid spurious negative Jacobian errors while imposing BCs by
        // simply ignoring them. This should only be required in certain
        // special cases, see the DirichletBoundaries comments on this
        // parameter for more information.
        edge_fe->get_fe_map().set_jacobian_tolerance(dirichlet.jacobian_tolerance);

        // Pre-request FE data
        const std::vector<std::vector<OutputShape>> & phi = edge_fe->get_phi();
        const std::vector<Point> & xyz_values = edge_fe->get_xyz();
        const std::vector<Real> & JxW = edge_fe->get_JxW();

        // Only pre-request gradients for C1 projections
        const std::vector<std::vector<OutputGradient>> * dphi = nullptr;
        if ((cont == C_ONE) && (fe_type.family != SUBDIVISION))
          {
            const std::vector<std::vector<OutputGradient>> & ref_dphi = edge_fe->get_dphi();
            dphi = &ref_dphi;
          }

        // Vector to hold edge local DOF indices
        std::vector<unsigned int> edge_dofs;

        // Get a reference to the "is_boundary_edge" flags for the
        // current DirichletBoundary object.  In case the map does not
        // contain an entry for this DirichletBoundary object, it
        // means there are no boundary edges active.
        auto is_boundary_edge_it = sebi.is_boundary_edge_map.find(&dirichlet);

        for (unsigned int e=0; e != sebi.n_edges; ++e)
          {
            if (is_boundary_edge_it == sebi.is_boundary_edge_map.end() ||
                !is_boundary_edge_it->second[e])
              continue;

            FEInterface::dofs_on_edge(elem, dim, fe_type, e,
                                      edge_dofs);

            const unsigned int n_edge_dofs =
              cast_int<unsigned int>(edge_dofs.size());

            // Some edge dofs are on nodes and already
            // fixed, others are free to calculate
            unsigned int free_dofs = 0;
            for (unsigned int i=0; i != n_edge_dofs; ++i)
              if (!dof_is_fixed[edge_dofs[i]])
                free_dof[free_dofs++] = i;

            // There may be nothing to project
            if (!free_dofs)
              continue;

            Ke.resize (free_dofs, free_dofs); Ke.zero();
            Fe.resize (free_dofs); Fe.zero();
            // The new edge coefficients
            DenseVector<Number> Uedge(free_dofs);

            // Initialize FE data on the edge
            edge_fe->edge_reinit(elem, e);
            const unsigned int n_qp = fem_context.get_edge_qrule().n_points();

            // Loop over the quadrature points
            for (unsigned int qp=0; qp<n_qp; qp++)
              {
                // solution at the quadrature point
                OutputNumber fineval(0);
                libMesh::RawAccessor<OutputNumber> f_accessor( fineval, dim );

                for (unsigned int c = 0; c < n_vec_dim; c++)
                  f_accessor(c) =
                    f_component(f, f_fem, context.get(), var_component+c,
                                xyz_values[qp], time);

                // solution grad at the quadrature point
                OutputNumberGradient finegrad;
                libMesh::RawAccessor<OutputNumberGradient> g_accessor( finegrad, dim );

                unsigned int g_rank;
                switch( FEInterface::field_type( fe_type ) )
                  {
                  case TYPE_SCALAR:
                    {
                      g_rank = 1;
                      break;
                    }
                  case TYPE_VECTOR:
                    {
                      g_rank = 2;
                      break;
                    }
                  default:
                    libmesh_error_msg("Unknown field type!");
                  }

                if (cont == C_ONE)
                  for (unsigned int c = 0; c < n_vec_dim; c++)
                    for (unsigned int d = 0; d < g_rank; d++)
                      g_accessor(c + d*dim ) =
                        g_component(g, g_fem, context.get(), var_component,
                                    xyz_values[qp], time)(c);

                // Form edge projection matrix
                for (unsigned int sidei=0, freei=0; sidei != n_edge_dofs; ++sidei)
                  {
                    unsigned int i = edge_dofs[sidei];
                    // fixed DoFs aren't test functions
                    if (dof_is_fixed[i])
                      continue;
                    for (unsigned int sidej=0, freej=0; sidej != n_edge_dofs; ++sidej)
                      {
                        unsigned int j = edge_dofs[sidej];
                        if (dof_is_fixed[j])
                          Fe(freei) -= phi[i][qp] * phi[j][qp] *
                            JxW[qp] * Ue(j);
                        else
                          Ke(freei,freej) += phi[i][qp] *
                            phi[j][qp] * JxW[qp];
                        if (cont == C_ONE)
                          {
                            if (dof_is_fixed[j])
                              Fe(freei) -= ((*dphi)[i][qp].contract((*dphi)[j][qp]) ) *
                                JxW[qp] * Ue(j);
                            else
                              Ke(freei,freej) += ((*dphi)[i][qp].contract((*dphi)[j][qp]))
                                * JxW[qp];
                          }
                        if (!dof_is_fixed[j])
                          freej++;
                      }
                    Fe(freei) += phi[i][qp] * fineval * JxW[qp];
                    if (cont == C_ONE)
                      Fe(freei) += (finegrad.contract( (*dphi)[i][qp]) ) *
                        JxW[qp];
                    freei++;
                  }
              }

            Ke.cholesky_solve(Fe, Uedge);

            // Transfer new edge solutions to element
            for (unsigned int i=0; i != free_dofs; ++i)
              {
                Number & ui = Ue(edge_dofs[free_dof[i]]);
                libmesh_assert(std::abs(ui) < TOLERANCE ||
                               std::abs(ui - Uedge(i)) < TOLERANCE);
                ui = Uedge(i);
                dof_is_fixed[edge_dofs[free_dof[i]]] = true;
              }
          } // end for (e = 0..n_edges)
      } // end if (dim > 2 && cont != DISCONTINUOUS)

        // Project any side values (edges in 2D, faces in 3D)
    if (dim > 1 && cont != DISCONTINUOUS)
      {
        FEGenericBase<OutputType> * side_fe = nullptr;
        fem_context.get_side_fe(var, side_fe);

        // Set tolerance on underlying FEMap object. This will allow us to
        // avoid spurious negative Jacobian errors while imposing BCs by
        // simply ignoring them. This should only be required in certain
        // special cases, see the DirichletBoundaries comments on this
        // parameter for more information.
        side_fe->get_fe_map().set_jacobian_tolerance(dirichlet.jacobian_tolerance);

        // Pre-request FE data
        const std::vector<std::vector<OutputShape>> & phi = side_fe->get_phi();
        const std::vector<Point> & xyz_values = side_fe->get_xyz();
        const std::vector<Real> & JxW = side_fe->get_JxW();

        // Only pre-request gradients for C1 projections
        const std::vector<std::vector<OutputGradient>> * dphi = nullptr;
        if ((cont == C_ONE) && (fe_type.family != SUBDIVISION))
          {
            const std::vector<std::vector<OutputGradient>> & ref_dphi = side_fe->get_dphi();
            dphi = &ref_dphi;
          }

        // Vector to hold side local DOF indices
        std::vector<unsigned int> side_dofs;

        // Get a reference to the "is_boundary_side" flags for the
        // current DirichletBoundary object.  In case the map does not
        // contain an entry for this DirichletBoundary object, it
        // means there are no boundary sides active.
        auto is_boundary_side_it = sebi.is_boundary_side_map.find(&dirichlet);

        for (unsigned int s=0; s != sebi.n_sides; ++s)
          {
            if (is_boundary_side_it == sebi.is_boundary_side_map.end() ||
                !is_boundary_side_it->second[s])
              continue;

            FEInterface::dofs_on_side(elem, dim, fe_type, s,
                                      side_dofs);

            const unsigned int n_side_dofs =
              cast_int<unsigned int>(side_dofs.size());

            // Some side dofs are on nodes/edges and already
            // fixed, others are free to calculate
            unsigned int free_dofs = 0;
            for (unsigned int i=0; i != n_side_dofs; ++i)
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
            side_fe->reinit(elem, s);
            const unsigned int n_qp = fem_context.get_side_qrule().n_points();

            // Loop over the quadrature points
            for (unsigned int qp=0; qp<n_qp; qp++)
              {
                // solution at the quadrature point
                OutputNumber fineval(0);
                libMesh::RawAccessor<OutputNumber> f_accessor( fineval, dim );

                for (unsigned int c = 0; c < n_vec_dim; c++)
                  f_accessor(c) =
                    f_component(f, f_fem, context.get(), var_component+c,
                                xyz_values[qp], time);

                // solution grad at the quadrature point
                OutputNumberGradient finegrad;
                libMesh::RawAccessor<OutputNumberGradient> g_accessor( finegrad, dim );

                unsigned int g_rank;
                switch( FEInterface::field_type( fe_type ) )
                  {
                  case TYPE_SCALAR:
                    {
                      g_rank = 1;
                      break;
                    }
                  case TYPE_VECTOR:
                    {
                      g_rank = 2;
                      break;
                    }
                  default:
                    libmesh_error_msg("Unknown field type!");
                  }

                if (cont == C_ONE)
                  for (unsigned int c = 0; c < n_vec_dim; c++)
                    for (unsigned int d = 0; d < g_rank; d++)
                      g_accessor(c + d*dim ) =
                        g_component(g, g_fem, context.get(), var_component,
                                    xyz_values[qp], time)(c);

                // Form side projection matrix
                for (unsigned int sidei=0, freei=0; sidei != n_side_dofs; ++sidei)
                  {
                    unsigned int i = side_dofs[sidei];
                    // fixed DoFs aren't test functions
                    if (dof_is_fixed[i])
                      continue;
                    for (unsigned int sidej=0, freej=0; sidej != n_side_dofs; ++sidej)
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
                              Fe(freei) -= ((*dphi)[i][qp].contract((*dphi)[j][qp])) *
                                JxW[qp] * Ue(j);
                            else
                              Ke(freei,freej) += ((*dphi)[i][qp].contract((*dphi)[j][qp]))
                                * JxW[qp];
                          }
                        if (!dof_is_fixed[j])
                          freej++;
                      }
                    Fe(freei) += (fineval * phi[i][qp]) * JxW[qp];
                    if (cont == C_ONE)
                      Fe(freei) += (finegrad.contract((*dphi)[i][qp])) *
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
          } // end for (s = 0..n_sides)
      } // end if (dim > 1 && cont != DISCONTINUOUS)

        // Project any shellface values
    if (dim == 2 && cont != DISCONTINUOUS)
      {
        FEGenericBase<OutputType> * fe = nullptr;
        fem_context.get_element_fe(var, fe, dim);

        // Set tolerance on underlying FEMap object. This will allow us to
        // avoid spurious negative Jacobian errors while imposing BCs by
        // simply ignoring them. This should only be required in certain
        // special cases, see the DirichletBoundaries comments on this
        // parameter for more information.
        fe->get_fe_map().set_jacobian_tolerance(dirichlet.jacobian_tolerance);

        // Pre-request FE data
        const std::vector<std::vector<OutputShape>> & phi = fe->get_phi();
        const std::vector<Point> & xyz_values = fe->get_xyz();
        const std::vector<Real> & JxW = fe->get_JxW();

        // Only pre-request gradients for C1 projections
        const std::vector<std::vector<OutputGradient>> * dphi = nullptr;
        if ((cont == C_ONE) && (fe_type.family != SUBDIVISION))
          {
            const std::vector<std::vector<OutputGradient>> & ref_dphi = fe->get_dphi();
            dphi = &ref_dphi;
          }

        // Get a reference to the "is_boundary_shellface" flags for the
        // current DirichletBoundary object.  In case the map does not
        // contain an entry for this DirichletBoundary object, it
        // means there are no boundary shellfaces active.
        auto is_boundary_shellface_it = sebi.is_boundary_shellface_map.find(&dirichlet);

        for (unsigned int shellface=0; shellface != 2; ++shellface)
          {
            if (is_boundary_shellface_it == sebi.is_boundary_shellface_map.end() ||
                !is_boundary_shellface_it->second[shellface])
              continue;

            // A shellface has the same dof indices as the element itself
            std::vector<unsigned int> shellface_dofs(n_dofs);
            std::iota(shellface_dofs.begin(), shellface_dofs.end(), 0);

            // Some shellface dofs are on nodes/edges and already
            // fixed, others are free to calculate
            unsigned int free_dofs = 0;
            for (unsigned int i=0; i != n_dofs; ++i)
              if (!dof_is_fixed[shellface_dofs[i]])
                free_dof[free_dofs++] = i;

            // There may be nothing to project
            if (!free_dofs)
              continue;

            Ke.resize (free_dofs, free_dofs); Ke.zero();
            Fe.resize (free_dofs); Fe.zero();
            // The new shellface coefficients
            DenseVector<Number> Ushellface(free_dofs);

            // Initialize FE data on the element
            fe->reinit (elem);
            const unsigned int n_qp = fem_context.get_element_qrule().n_points();

            // Loop over the quadrature points
            for (unsigned int qp=0; qp<n_qp; qp++)
              {
                // solution at the quadrature point
                OutputNumber fineval(0);
                libMesh::RawAccessor<OutputNumber> f_accessor( fineval, dim );

                for (unsigned int c = 0; c < n_vec_dim; c++)
                  f_accessor(c) =
                    f_component(f, f_fem, context.get(), var_component+c,
                                xyz_values[qp], time);

                // solution grad at the quadrature point
                OutputNumberGradient finegrad;
                libMesh::RawAccessor<OutputNumberGradient> g_accessor( finegrad, dim );

                unsigned int g_rank;
                switch( FEInterface::field_type( fe_type ) )
                  {
                  case TYPE_SCALAR:
                    {
                      g_rank = 1;
                      break;
                    }
                  case TYPE_VECTOR:
                    {
                      g_rank = 2;
                      break;
                    }
                  default:
                    libmesh_error_msg("Unknown field type!");
                  }

                if (cont == C_ONE)
                  for (unsigned int c = 0; c < n_vec_dim; c++)
                    for (unsigned int d = 0; d < g_rank; d++)
                      g_accessor(c + d*dim ) =
                        g_component(g, g_fem, context.get(), var_component,
                                    xyz_values[qp], time)(c);

                // Form shellface projection matrix
                for (unsigned int shellfacei=0, freei=0;
                     shellfacei != n_dofs; ++shellfacei)
                  {
                    unsigned int i = shellface_dofs[shellfacei];
                    // fixed DoFs aren't test functions
                    if (dof_is_fixed[i])
                      continue;
                    for (unsigned int shellfacej=0, freej=0;
                         shellfacej != n_dofs; ++shellfacej)
                      {
                        unsigned int j = shellface_dofs[shellfacej];
                        if (dof_is_fixed[j])
                          Fe(freei) -= phi[i][qp] * phi[j][qp] *
                            JxW[qp] * Ue(j);
                        else
                          Ke(freei,freej) += phi[i][qp] *
                            phi[j][qp] * JxW[qp];
                        if (cont == C_ONE)
                          {
                            if (dof_is_fixed[j])
                              Fe(freei) -= ((*dphi)[i][qp].contract((*dphi)[j][qp])) *
                                JxW[qp] * Ue(j);
                            else
                              Ke(freei,freej) += ((*dphi)[i][qp].contract((*dphi)[j][qp]))
                                * JxW[qp];
                          }
                        if (!dof_is_fixed[j])
                          freej++;
                      }
                    Fe(freei) += (fineval * phi[i][qp]) * JxW[qp];
                    if (cont == C_ONE)
                      Fe(freei) += (finegrad.contract((*dphi)[i][qp])) *
                        JxW[qp];
                    freei++;
                  }
              }

            Ke.cholesky_solve(Fe, Ushellface);

            // Transfer new shellface solutions to element
            for (unsigned int i=0; i != free_dofs; ++i)
              {
                Number & ui = Ue(shellface_dofs[free_dof[i]]);
                libmesh_assert(std::abs(ui) < TOLERANCE ||
                               std::abs(ui - Ushellface(i)) < TOLERANCE);
                ui = Ushellface(i);
                dof_is_fixed[shellface_dofs[free_dof[i]]] = true;
              }
          } // end for (shellface = 0..2)
      } // end if (dim == 2 && cont != DISCONTINUOUS)

    // Lock the DofConstraints since it is shared among threads.
    {
      Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

      for (unsigned int i = 0; i < n_dofs; i++)
        {
          DofConstraintRow empty_row;
          if (dof_is_fixed[i] && !libmesh_isnan(Ue(i)))
            add_fn (dof_indices[i], empty_row, Ue(i));
        }
    }
  } // apply_dirichlet_impl

public:
  ConstrainDirichlet (DofMap & dof_map_in,
                      const MeshBase & mesh_in,
                      const Real time_in,
                      const DirichletBoundaries & dirichlets_in,
                      const AddConstraint & add_in) :
    dof_map(dof_map_in),
    mesh(mesh_in),
    time(time_in),
    dirichlets(dirichlets_in),
    add_fn(add_in) { }

  // This class can be default copy/move constructed.
  ConstrainDirichlet (ConstrainDirichlet &&) = default;
  ConstrainDirichlet (const ConstrainDirichlet &) = default;

  // This class cannot be default copy/move assigned because it
  // contains reference members.
  ConstrainDirichlet & operator= (const ConstrainDirichlet &) = delete;
  ConstrainDirichlet & operator= (ConstrainDirichlet &&) = delete;

  void operator()(const ConstElemRange & range) const
  {
    /**
     * This method examines an arbitrary boundary solution to calculate
     * corresponding Dirichlet constraints on the current mesh.  The
     * input function \p f gives the arbitrary solution.
     */

    // Figure out which System the DirichletBoundary objects are
    // defined for. We break out of the loop as soon as we encounter a
    // valid System pointer, the assumption is thus that all Variables
    // are defined on the same System.
    System * system = nullptr;

    // Map from boundary_id -> set<pair<id,DirichletBoundary*>> objects which
    // are active on that boundary_id. Later we will use this to determine
    // which DirichletBoundary objects to loop over for each Elem.
    std::map<boundary_id_type, std::set<std::pair<unsigned int, DirichletBoundary *>>>
      boundary_id_to_ordered_dirichlet_boundaries;

    for (auto dirichlet_id : index_range(dirichlets))
      {
        // Pointer to the current DirichletBoundary object
        const auto & dirichlet = dirichlets[dirichlet_id];

        // Construct mapping from boundary_id -> (dirichlet_id, DirichletBoundary)
        for (const auto & b_id : dirichlet->b)
          boundary_id_to_ordered_dirichlet_boundaries[b_id].emplace(dirichlet_id, dirichlet);

        for (const auto & var : dirichlet->variables)
          {
            const Variable & variable = dof_map.variable(var);
            auto current_system = variable.system();

            if (!system)
              system = current_system;
            else
              libmesh_error_msg_if(current_system != system,
                                   "All variables should be defined on the same System");
          }
      }

    // If we found no System, it could be because none of the
    // Variables have one defined, or because there are
    // DirichletBoundary objects with no Variables defined on
    // them. These situations both indicate a likely error in the
    // setup of a problem, so let's throw an error in this case.
    libmesh_error_msg_if(!system, "Valid System not found for any Variables.");

    // Construct a FEMContext object for the System on which the
    // Variables in our DirichletBoundary objects are defined. This
    // will be used in the apply_dirichlet_impl() function.
    auto fem_context = libmesh_make_unique<FEMContext>(*system);

    // At the time we are using this FEMContext, the current_local_solution
    // vector is not initialized, but also we don't need it, so set
    // the algebraic_type flag to DOFS_ONLY.
    fem_context->set_algebraic_type(FEMContext::DOFS_ONLY);

    // Boundary info for the current mesh
    const BoundaryInfo & boundary_info = mesh.get_boundary_info();

    // This object keeps track of the BoundaryInfo for a single Elem
    SingleElemBoundaryInfo sebi(boundary_info, boundary_id_to_ordered_dirichlet_boundaries);

    // Iterate over all the elements in the range
    for (const auto & elem : range)
      {
        // We only calculate Dirichlet constraints on active
        // elements
        if (!elem->active())
          continue;

        // Reinitialize BoundaryInfo data structures for the current elem
        bool has_dirichlet_constraint = sebi.reinit(elem);

        // If this Elem has no boundary ids, go to the next one.
        if (!has_dirichlet_constraint)
          continue;

        for (const auto & db_pair : sebi.ordered_dbs)
          {
            // Get pointer to the DirichletBoundary object
            const auto & dirichlet = db_pair.second;

            // Loop over all the variables which this DirichletBoundary object is responsible for
            for (const auto & var : dirichlet->variables)
              {
                const Variable & variable = dof_map.variable(var);

                // Make sure that the Variable and the DofMap agree on
                // what number this variable is.
                libmesh_assert_equal_to(variable.number(), var);

                const FEType & fe_type = variable.type();

                if (fe_type.family == SCALAR)
                  continue;

                switch( FEInterface::field_type( fe_type ) )
                  {
                  case TYPE_SCALAR:
                    {
                      // For Lagrange FEs we don't need to do a full
                      // blown projection, we can just interpolate
                      // values directly.
                      if (fe_type.family == LAGRANGE)
                        this->apply_lagrange_dirichlet_impl<Real>( sebi, variable, *dirichlet, *fem_context );
                      else
                        this->apply_dirichlet_impl<Real>( sebi, variable, *dirichlet, *fem_context );
                      break;
                    }
                  case TYPE_VECTOR:
                    {
                      this->apply_dirichlet_impl<RealGradient>( sebi, variable, *dirichlet, *fem_context );
                      break;
                    }
                  default:
                    libmesh_error_msg("Unknown field type!");
                  }
              } // for (var : variables)
          } // for (db_pair : ordered_dbs)
      } // for (elem : range)
  } // operator()

}; // class ConstrainDirichlet


#endif // LIBMESH_ENABLE_DIRICHLET


} // anonymous namespace



namespace libMesh
{

// ------------------------------------------------------------
// DofMap member functions

#ifdef LIBMESH_ENABLE_CONSTRAINTS


dof_id_type DofMap::n_constrained_dofs() const
{
  parallel_object_only();

  dof_id_type nc_dofs = this->n_local_constrained_dofs();
  this->comm().sum(nc_dofs);
  return nc_dofs;
}


dof_id_type DofMap::n_local_constrained_dofs() const
{
  const DofConstraints::const_iterator lower =
    _dof_constraints.lower_bound(this->first_dof()),
    upper =
    _dof_constraints.upper_bound(this->end_dof()-1);

  return cast_int<dof_id_type>(std::distance(lower, upper));
}



void DofMap::create_dof_constraints(const MeshBase & mesh, Real time)
{
  parallel_object_only();

  LOG_SCOPE("create_dof_constraints()", "DofMap");

  libmesh_assert (mesh.is_prepared());

  // The user might have set boundary conditions after the mesh was
  // prepared; we should double-check that those boundary conditions
  // are still consistent.
#ifdef DEBUG
  MeshTools::libmesh_assert_valid_boundary_ids(mesh);
#endif

  // We might get constraint equations from AMR hanging nodes in 2D/3D
  // or from boundary conditions in any dimension
  const bool possible_local_constraints = false
    || !mesh.n_elem()
#ifdef LIBMESH_ENABLE_AMR
    || mesh.mesh_dimension() > 1
#endif
#ifdef LIBMESH_ENABLE_PERIODIC
    || !_periodic_boundaries->empty()
#endif
#ifdef LIBMESH_ENABLE_DIRICHLET
    || !_dirichlet_boundaries->empty()
#endif
    ;

  // Even if we don't have constraints, another processor might.
  bool possible_global_constraints = possible_local_constraints;
#if defined(LIBMESH_ENABLE_PERIODIC) || defined(LIBMESH_ENABLE_DIRICHLET) || defined(LIBMESH_ENABLE_AMR)
  libmesh_assert(this->comm().verify(mesh.is_serial()));

  this->comm().max(possible_global_constraints);
#endif

  if (!possible_global_constraints)
    {
      // Clear out any old constraints; maybe the user just deleted
      // their last remaining dirichlet/periodic/user constraint?
      // Note: any _stashed_dof_constraints are not cleared as it
      // may be the user's intention to restore them later.
#ifdef LIBMESH_ENABLE_CONSTRAINTS
      _dof_constraints.clear();
      _primal_constraint_values.clear();
      _adjoint_constraint_values.clear();
#endif
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
      _node_constraints.clear();
#endif

      return;
    }

  // Here we build the hanging node constraints.  This is done
  // by enforcing the condition u_a = u_b along hanging sides.
  // u_a = u_b is collocated at the nodes of side a, which gives
  // one row of the constraint matrix.

  // Processors only compute their local constraints
  ConstElemRange range (mesh.local_elements_begin(),
                        mesh.local_elements_end());

  // Global computation fails if we're using a FEMFunctionBase BC on a
  // ReplicatedMesh in parallel
  // ConstElemRange range (mesh.elements_begin(),
  //                       mesh.elements_end());

  // compute_periodic_constraints requires a point_locator() from our
  // Mesh, but point_locator() construction is parallel and threaded.
  // Rather than nest threads within threads we'll make sure it's
  // preconstructed.
#ifdef LIBMESH_ENABLE_PERIODIC
  bool need_point_locator = !_periodic_boundaries->empty() && !range.empty();

  this->comm().max(need_point_locator);

  if (need_point_locator)
    mesh.sub_point_locator();
#endif

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  // recalculate node constraints from scratch
  _node_constraints.clear();

  Threads::parallel_for (range,
                         ComputeNodeConstraints (_node_constraints,
#ifdef LIBMESH_ENABLE_PERIODIC
                                                 *_periodic_boundaries,
#endif
                                                 mesh));
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS


  // recalculate dof constraints from scratch
  // Note: any _stashed_dof_constraints are not cleared as it
  // may be the user's intention to restore them later.
  _dof_constraints.clear();
  _primal_constraint_values.clear();
  _adjoint_constraint_values.clear();

  // Look at all the variables in the system.  Reset the element
  // range at each iteration -- there is no need to reconstruct it.
  for (unsigned int variable_number=0; variable_number<this->n_variables();
       ++variable_number, range.reset())
    Threads::parallel_for (range,
                           ComputeConstraints (_dof_constraints,
                                               *this,
#ifdef LIBMESH_ENABLE_PERIODIC
                                               *_periodic_boundaries,
#endif
                                               mesh,
                                               variable_number));

#ifdef LIBMESH_ENABLE_DIRICHLET

  if (!_dirichlet_boundaries->empty())
    {
      // Sanity check that the boundary ids associated with the
      // DirichletBoundary objects are actually present in the
      // mesh.
      for (const auto & dirichlet : *_dirichlet_boundaries)
        this->check_dirichlet_bcid_consistency(mesh, *dirichlet);

      // Threaded loop over local over elems applying all Dirichlet BCs
      Threads::parallel_for
        (range,
         ConstrainDirichlet(*this, mesh, time, *_dirichlet_boundaries,
                            AddPrimalConstraint(*this)));

      // Threaded loop over local over elems per QOI applying all adjoint
      // Dirichlet BCs.  Note that the ConstElemRange is reset before each
      // execution of Threads::parallel_for().

      for (auto qoi_index : index_range(_adjoint_dirichlet_boundaries))
        {
          Threads::parallel_for
            (range.reset(),
             ConstrainDirichlet(*this, mesh, time, *(_adjoint_dirichlet_boundaries[qoi_index]),
                                AddAdjointConstraint(*this, qoi_index)));
        }
    }

#endif // LIBMESH_ENABLE_DIRICHLET
}



void DofMap::add_constraint_row (const dof_id_type dof_number,
                                 const DofConstraintRow & constraint_row,
                                 const Number constraint_rhs,
                                 const bool forbid_constraint_overwrite)
{
  // Optionally allow the user to overwrite constraints.  Defaults to false.
  libmesh_error_msg_if(forbid_constraint_overwrite && this->is_constrained_dof(dof_number),
                       "ERROR: DOF " << dof_number << " was already constrained!");

  libmesh_assert_less(dof_number, this->n_dofs());

  // There is an implied "1" on the diagonal of the constraint row, and the user
  // should not try to manually set _any_ value on the diagonal.
  libmesh_assert_msg(!constraint_row.count(dof_number),
                     "Error: constraint_row for dof " << dof_number <<
                     " should not contain an entry for dof " << dof_number);

#ifndef NDEBUG
  for (const auto & pr : constraint_row)
    libmesh_assert_less(pr.first, this->n_dofs());
#endif

  // We don't get insert_or_assign until C++17 so we make do.
  std::pair<DofConstraints::iterator, bool> it =
    _dof_constraints.emplace(dof_number, constraint_row);
  if (!it.second)
    it.first->second = constraint_row;

  std::pair<DofConstraintValueMap::iterator, bool> rhs_it =
    _primal_constraint_values.emplace(dof_number, constraint_rhs);
  if (!rhs_it.second)
    rhs_it.first->second = constraint_rhs;
}


void DofMap::add_adjoint_constraint_row (const unsigned int qoi_index,
                                         const dof_id_type dof_number,
                                         const DofConstraintRow & /*constraint_row*/,
                                         const Number constraint_rhs,
                                         const bool forbid_constraint_overwrite)
{
  // Optionally allow the user to overwrite constraints.  Defaults to false.
  if (forbid_constraint_overwrite)
    {
      libmesh_error_msg_if(!this->is_constrained_dof(dof_number),
                           "ERROR: DOF " << dof_number << " has no corresponding primal constraint!");
#ifndef NDEBUG
      // No way to do this without a non-normalized tolerance?

      // // If the user passed in more than just the rhs, let's check the
      // // coefficients for consistency
      // if (!constraint_row.empty())
      //   {
      //     DofConstraintRow row = _dof_constraints[dof_number];
      //     for (const auto & pr : row)
      //       libmesh_assert(constraint_row.find(pr.first)->second == pr.second);
      //   }
      //
      // if (_adjoint_constraint_values[qoi_index].find(dof_number) !=
      //     _adjoint_constraint_values[qoi_index].end())
      //   libmesh_assert_equal_to(_adjoint_constraint_values[qoi_index][dof_number],
      //                           constraint_rhs);

#endif
    }

  // Creates the map of rhs values if it doesn't already exist; then
  // adds the current value to that map

  // We don't get insert_or_assign until C++17 so we make do.
  std::pair<DofConstraintValueMap::iterator, bool> rhs_it =
    _adjoint_constraint_values[qoi_index].emplace(dof_number, constraint_rhs);
  if (!rhs_it.second)
    rhs_it.first->second = constraint_rhs;
}




void DofMap::print_dof_constraints(std::ostream & os,
                                   bool print_nonlocal) const
{
  parallel_object_only();

  std::string local_constraints =
    this->get_local_constraints(print_nonlocal);

  if (this->processor_id())
    {
      this->comm().send(0, local_constraints);
    }
  else
    {
      os << "Processor 0:\n";
      os << local_constraints;

      for (auto p : IntRange<processor_id_type>(1, this->n_processors()))
        {
          this->comm().receive(p, local_constraints);
          os << "Processor " << p << ":\n";
          os << local_constraints;
        }
    }
}

std::string DofMap::get_local_constraints(bool print_nonlocal) const
{
  std::ostringstream os;
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  if (print_nonlocal)
    os << "All ";
  else
    os << "Local ";

  os << "Node Constraints:"
     << std::endl;

  for (const auto & pr : _node_constraints)
    {
      const Node * node = pr.first;

      // Skip non-local nodes if requested
      if (!print_nonlocal &&
          node->processor_id() != this->processor_id())
        continue;

      const NodeConstraintRow & row = pr.second.first;
      const Point & offset = pr.second.second;

      os << "Constraints for Node id " << node->id()
         << ": \t";

      for (const auto & item : row)
        os << " (" << item.first->id() << "," << item.second << ")\t";

      os << "rhs: " << offset;

      os << std::endl;
    }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

  if (print_nonlocal)
    os << "All ";
  else
    os << "Local ";

  os << "DoF Constraints:"
     << std::endl;

  for (const auto & pr : _dof_constraints)
    {
      const dof_id_type i = pr.first;

      // Skip non-local dofs if requested
      if (!print_nonlocal && !this->local_index(i))
        continue;

      const DofConstraintRow & row = pr.second;
      DofConstraintValueMap::const_iterator rhsit =
        _primal_constraint_values.find(i);
      const Number rhs = (rhsit == _primal_constraint_values.end()) ?
        0 : rhsit->second;

      os << "Constraints for DoF " << i
         << ": \t";

      for (const auto & item : row)
        os << " (" << item.first << "," << item.second << ")\t";

      os << "rhs: " << rhs;
      os << std::endl;
    }

  for (unsigned int qoi_index = 0,
       n_qois = cast_int<unsigned int>(_adjoint_dirichlet_boundaries.size());
       qoi_index != n_qois; ++qoi_index)
    {
      os << "Adjoint " << qoi_index << " DoF rhs values:"
         << std::endl;

      AdjointDofConstraintValues::const_iterator adjoint_map_it =
        _adjoint_constraint_values.find(qoi_index);

      if (adjoint_map_it != _adjoint_constraint_values.end())
        for (const auto & pr : adjoint_map_it->second)
          {
            const dof_id_type i = pr.first;

            // Skip non-local dofs if requested
            if (!print_nonlocal && !this->local_index(i))
              continue;

            const Number rhs = pr.second;

            os << "RHS for DoF " << i
               << ": " << rhs;

            os << std::endl;
          }
    }

  return os.str();
}



void DofMap::constrain_element_matrix (DenseMatrix<Number> & matrix,
                                       std::vector<dof_id_type> & elem_dofs,
                                       bool asymmetric_constraint_rows) const
{
  libmesh_assert_equal_to (elem_dofs.size(), matrix.m());
  libmesh_assert_equal_to (elem_dofs.size(), matrix.n());

  // check for easy return
  if (this->_dof_constraints.empty())
    return;

  // The constrained matrix is built up as C^T K C.
  DenseMatrix<Number> C;


  this->build_constraint_matrix (C, elem_dofs);

  LOG_SCOPE("constrain_elem_matrix()", "DofMap");

  // It is possible that the matrix is not constrained at all.
  if ((C.m() == matrix.m()) &&
      (C.n() == elem_dofs.size())) // It the matrix is constrained
    {
      // Compute the matrix-matrix-matrix product C^T K C
      matrix.left_multiply_transpose  (C);
      matrix.right_multiply (C);


      libmesh_assert_equal_to (matrix.m(), matrix.n());
      libmesh_assert_equal_to (matrix.m(), elem_dofs.size());
      libmesh_assert_equal_to (matrix.n(), elem_dofs.size());


      for (unsigned int i=0,
           n_elem_dofs = cast_int<unsigned int>(elem_dofs.size());
           i != n_elem_dofs; i++)
        // If the DOF is constrained
        if (this->is_constrained_dof(elem_dofs[i]))
          {
            for (auto j : make_range(matrix.n()))
              matrix(i,j) = 0.;

            matrix(i,i) = 1.;

            if (asymmetric_constraint_rows)
              {
                DofConstraints::const_iterator
                  pos = _dof_constraints.find(elem_dofs[i]);

                libmesh_assert (pos != _dof_constraints.end());

                const DofConstraintRow & constraint_row = pos->second;

                // This is an overzealous assertion in the presence of
                // heterogenous constraints: we now can constrain "u_i = c"
                // with no other u_j terms involved.
                //
                // libmesh_assert (!constraint_row.empty());

                for (const auto & item : constraint_row)
                  for (unsigned int j=0; j != n_elem_dofs; j++)
                    if (elem_dofs[j] == item.first)
                      matrix(i,j) = -item.second;
              }
          }
    } // end if is constrained...
}



void DofMap::constrain_element_matrix_and_vector (DenseMatrix<Number> & matrix,
                                                  DenseVector<Number> & rhs,
                                                  std::vector<dof_id_type> & elem_dofs,
                                                  bool asymmetric_constraint_rows) const
{
  libmesh_assert_equal_to (elem_dofs.size(), matrix.m());
  libmesh_assert_equal_to (elem_dofs.size(), matrix.n());
  libmesh_assert_equal_to (elem_dofs.size(), rhs.size());

  // check for easy return
  if (this->_dof_constraints.empty())
    return;

  // The constrained matrix is built up as C^T K C.
  // The constrained RHS is built up as C^T F
  DenseMatrix<Number> C;

  this->build_constraint_matrix (C, elem_dofs);

  LOG_SCOPE("cnstrn_elem_mat_vec()", "DofMap");

  // It is possible that the matrix is not constrained at all.
  if ((C.m() == matrix.m()) &&
      (C.n() == elem_dofs.size())) // It the matrix is constrained
    {
      // Compute the matrix-matrix-matrix product C^T K C
      matrix.left_multiply_transpose  (C);
      matrix.right_multiply (C);


      libmesh_assert_equal_to (matrix.m(), matrix.n());
      libmesh_assert_equal_to (matrix.m(), elem_dofs.size());
      libmesh_assert_equal_to (matrix.n(), elem_dofs.size());


      for (unsigned int i=0,
           n_elem_dofs = cast_int<unsigned int>(elem_dofs.size());
           i != n_elem_dofs; i++)
        if (this->is_constrained_dof(elem_dofs[i]))
          {
            for (auto j : make_range(matrix.n()))
              matrix(i,j) = 0.;

            // If the DOF is constrained
            matrix(i,i) = 1.;

            // This will put a nonsymmetric entry in the constraint
            // row to ensure that the linear system produces the
            // correct value for the constrained DOF.
            if (asymmetric_constraint_rows)
              {
                DofConstraints::const_iterator
                  pos = _dof_constraints.find(elem_dofs[i]);

                libmesh_assert (pos != _dof_constraints.end());

                const DofConstraintRow & constraint_row = pos->second;

                // p refinement creates empty constraint rows
                //    libmesh_assert (!constraint_row.empty());

                for (const auto & item : constraint_row)
                  for (unsigned int j=0; j != n_elem_dofs; j++)
                    if (elem_dofs[j] == item.first)
                      matrix(i,j) = -item.second;
              }
          }


      // Compute the matrix-vector product C^T F
      DenseVector<Number> old_rhs(rhs);

      // compute matrix/vector product
      C.vector_mult_transpose(rhs, old_rhs);
    } // end if is constrained...
}



void DofMap::heterogenously_constrain_element_matrix_and_vector (DenseMatrix<Number> & matrix,
                                                                 DenseVector<Number> & rhs,
                                                                 std::vector<dof_id_type> & elem_dofs,
                                                                 bool asymmetric_constraint_rows,
                                                                 int qoi_index) const
{
  libmesh_assert_equal_to (elem_dofs.size(), matrix.m());
  libmesh_assert_equal_to (elem_dofs.size(), matrix.n());
  libmesh_assert_equal_to (elem_dofs.size(), rhs.size());

  // check for easy return
  if (this->_dof_constraints.empty())
    return;

  // The constrained matrix is built up as C^T K C.
  // The constrained RHS is built up as C^T (F - K H)
  DenseMatrix<Number> C;
  DenseVector<Number> H;

  this->build_constraint_matrix_and_vector (C, H, elem_dofs, qoi_index);

  LOG_SCOPE("hetero_cnstrn_elem_mat_vec()", "DofMap");

  // It is possible that the matrix is not constrained at all.
  if ((C.m() == matrix.m()) &&
      (C.n() == elem_dofs.size())) // It the matrix is constrained
    {
      // We may have rhs values to use later
      const DofConstraintValueMap * rhs_values = nullptr;
      if (qoi_index < 0)
        rhs_values = &_primal_constraint_values;
      else
        {
          const AdjointDofConstraintValues::const_iterator
            it = _adjoint_constraint_values.find(qoi_index);
          if (it != _adjoint_constraint_values.end())
            rhs_values = &it->second;
        }

      // Compute matrix/vector product K H
      DenseVector<Number> KH;
      matrix.vector_mult(KH, H);

      // Compute the matrix-vector product C^T (F - KH)
      DenseVector<Number> F_minus_KH(rhs);
      F_minus_KH -= KH;
      C.vector_mult_transpose(rhs, F_minus_KH);

      // Compute the matrix-matrix-matrix product C^T K C
      matrix.left_multiply_transpose  (C);
      matrix.right_multiply (C);

      libmesh_assert_equal_to (matrix.m(), matrix.n());
      libmesh_assert_equal_to (matrix.m(), elem_dofs.size());
      libmesh_assert_equal_to (matrix.n(), elem_dofs.size());

      for (unsigned int i=0,
           n_elem_dofs = cast_int<unsigned int>(elem_dofs.size());
           i != n_elem_dofs; i++)
        {
          const dof_id_type dof_id = elem_dofs[i];

          if (this->is_constrained_dof(dof_id))
            {
              for (auto j : make_range(matrix.n()))
                matrix(i,j) = 0.;

              // If the DOF is constrained
              matrix(i,i) = 1.;

              // This will put a nonsymmetric entry in the constraint
              // row to ensure that the linear system produces the
              // correct value for the constrained DOF.
              if (asymmetric_constraint_rows)
                {
                  DofConstraints::const_iterator
                    pos = _dof_constraints.find(dof_id);

                  libmesh_assert (pos != _dof_constraints.end());

                  const DofConstraintRow & constraint_row = pos->second;

                  for (const auto & item : constraint_row)
                    for (unsigned int j=0; j != n_elem_dofs; j++)
                      if (elem_dofs[j] == item.first)
                        matrix(i,j) = -item.second;

                  if (rhs_values)
                    {
                      const DofConstraintValueMap::const_iterator valpos =
                        rhs_values->find(dof_id);

                      rhs(i) = (valpos == rhs_values->end()) ?
                        0 : valpos->second;
                    }
                }
              else
                rhs(i) = 0.;
            }
        }

    } // end if is constrained...
}



void DofMap::heterogenously_constrain_element_vector (const DenseMatrix<Number> & matrix,
                                                      DenseVector<Number> & rhs,
                                                      std::vector<dof_id_type> & elem_dofs,
                                                      bool asymmetric_constraint_rows,
                                                      int qoi_index) const
{
  libmesh_assert_equal_to (elem_dofs.size(), matrix.m());
  libmesh_assert_equal_to (elem_dofs.size(), matrix.n());
  libmesh_assert_equal_to (elem_dofs.size(), rhs.size());

  // check for easy return
  if (this->_dof_constraints.empty())
    return;

  // The constrained matrix is built up as C^T K C.
  // The constrained RHS is built up as C^T (F - K H)
  DenseMatrix<Number> C;
  DenseVector<Number> H;

  this->build_constraint_matrix_and_vector (C, H, elem_dofs, qoi_index);

  LOG_SCOPE("hetero_cnstrn_elem_vec()", "DofMap");

  // It is possible that the matrix is not constrained at all.
  if ((C.m() == matrix.m()) &&
      (C.n() == elem_dofs.size())) // It the matrix is constrained
    {
      // We may have rhs values to use later
      const DofConstraintValueMap * rhs_values = nullptr;
      if (qoi_index < 0)
        rhs_values = &_primal_constraint_values;
      else
        {
          const AdjointDofConstraintValues::const_iterator
            it = _adjoint_constraint_values.find(qoi_index);
          if (it != _adjoint_constraint_values.end())
            rhs_values = &it->second;
        }

      // Compute matrix/vector product K H
      DenseVector<Number> KH;
      matrix.vector_mult(KH, H);

      // Compute the matrix-vector product C^T (F - KH)
      DenseVector<Number> F_minus_KH(rhs);
      F_minus_KH -= KH;
      C.vector_mult_transpose(rhs, F_minus_KH);

      for (unsigned int i=0,
           n_elem_dofs = cast_int<unsigned int>(elem_dofs.size());
           i != n_elem_dofs; i++)
        {
          const dof_id_type dof_id = elem_dofs[i];

          if (this->is_constrained_dof(dof_id))
            {
              // This will put a nonsymmetric entry in the constraint
              // row to ensure that the linear system produces the
              // correct value for the constrained DOF.
              if (asymmetric_constraint_rows && rhs_values)
                {
                  const DofConstraintValueMap::const_iterator valpos =
                    rhs_values->find(dof_id);

                  rhs(i) = (valpos == rhs_values->end()) ?
                    0 : valpos->second;
                }
              else
                rhs(i) = 0.;
            }
        }

    } // end if is constrained...
}




void DofMap::constrain_element_matrix (DenseMatrix<Number> & matrix,
                                       std::vector<dof_id_type> & row_dofs,
                                       std::vector<dof_id_type> & col_dofs,
                                       bool asymmetric_constraint_rows) const
{
  libmesh_assert_equal_to (row_dofs.size(), matrix.m());
  libmesh_assert_equal_to (col_dofs.size(), matrix.n());

  // check for easy return
  if (this->_dof_constraints.empty())
    return;

  // The constrained matrix is built up as R^T K C.
  DenseMatrix<Number> R;
  DenseMatrix<Number> C;

  // Safeguard against the user passing us the same
  // object for row_dofs and col_dofs.  If that is done
  // the calls to build_matrix would fail
  std::vector<dof_id_type> orig_row_dofs(row_dofs);
  std::vector<dof_id_type> orig_col_dofs(col_dofs);

  this->build_constraint_matrix (R, orig_row_dofs);
  this->build_constraint_matrix (C, orig_col_dofs);

  LOG_SCOPE("constrain_elem_matrix()", "DofMap");

  row_dofs = orig_row_dofs;
  col_dofs = orig_col_dofs;

  bool constraint_found = false;

  // K_constrained = R^T K C

  if ((R.m() == matrix.m()) &&
      (R.n() == row_dofs.size()))
    {
      matrix.left_multiply_transpose  (R);
      constraint_found = true;
    }

  if ((C.m() == matrix.n()) &&
      (C.n() == col_dofs.size()))
    {
      matrix.right_multiply (C);
      constraint_found = true;
    }

  // It is possible that the matrix is not constrained at all.
  if (constraint_found)
    {
      libmesh_assert_equal_to (matrix.m(), row_dofs.size());
      libmesh_assert_equal_to (matrix.n(), col_dofs.size());


      for (unsigned int i=0,
           n_row_dofs = cast_int<unsigned int>(row_dofs.size());
           i != n_row_dofs; i++)
        if (this->is_constrained_dof(row_dofs[i]))
          {
            for (auto j : make_range(matrix.n()))
              {
                if (row_dofs[i] != col_dofs[j])
                  matrix(i,j) = 0.;
                else // If the DOF is constrained
                  matrix(i,j) = 1.;
              }

            if (asymmetric_constraint_rows)
              {
                DofConstraints::const_iterator
                  pos = _dof_constraints.find(row_dofs[i]);

                libmesh_assert (pos != _dof_constraints.end());

                const DofConstraintRow & constraint_row = pos->second;

                libmesh_assert (!constraint_row.empty());

                for (const auto & item : constraint_row)
                  for (unsigned int j=0,
                       n_col_dofs = cast_int<unsigned int>(col_dofs.size());
                       j != n_col_dofs; j++)
                    if (col_dofs[j] == item.first)
                      matrix(i,j) = -item.second;
              }
          }
    } // end if is constrained...
}



void DofMap::constrain_element_vector (DenseVector<Number> & rhs,
                                       std::vector<dof_id_type> & row_dofs,
                                       bool) const
{
  libmesh_assert_equal_to (rhs.size(), row_dofs.size());

  // check for easy return
  if (this->_dof_constraints.empty())
    return;

  // The constrained RHS is built up as R^T F.
  DenseMatrix<Number> R;

  this->build_constraint_matrix (R, row_dofs);

  LOG_SCOPE("constrain_elem_vector()", "DofMap");

  // It is possible that the vector is not constrained at all.
  if ((R.m() == rhs.size()) &&
      (R.n() == row_dofs.size())) // if the RHS is constrained
    {
      // Compute the matrix-vector product
      DenseVector<Number> old_rhs(rhs);
      R.vector_mult_transpose(rhs, old_rhs);

      libmesh_assert_equal_to (row_dofs.size(), rhs.size());

      for (unsigned int i=0,
           n_row_dofs = cast_int<unsigned int>(row_dofs.size());
           i != n_row_dofs; i++)
        if (this->is_constrained_dof(row_dofs[i]))
          {
            // If the DOF is constrained
            libmesh_assert (_dof_constraints.find(row_dofs[i]) != _dof_constraints.end());

            rhs(i) = 0;
          }
    } // end if the RHS is constrained.
}



void DofMap::constrain_element_dyad_matrix (DenseVector<Number> & v,
                                            DenseVector<Number> & w,
                                            std::vector<dof_id_type> & row_dofs,
                                            bool) const
{
  libmesh_assert_equal_to (v.size(), row_dofs.size());
  libmesh_assert_equal_to (w.size(), row_dofs.size());

  // check for easy return
  if (this->_dof_constraints.empty())
    return;

  // The constrained RHS is built up as R^T F.
  DenseMatrix<Number> R;

  this->build_constraint_matrix (R, row_dofs);

  LOG_SCOPE("cnstrn_elem_dyad_mat()", "DofMap");

  // It is possible that the vector is not constrained at all.
  if ((R.m() == v.size()) &&
      (R.n() == row_dofs.size())) // if the RHS is constrained
    {
      // Compute the matrix-vector products
      DenseVector<Number> old_v(v);
      DenseVector<Number> old_w(w);

      // compute matrix/vector product
      R.vector_mult_transpose(v, old_v);
      R.vector_mult_transpose(w, old_w);

      libmesh_assert_equal_to (row_dofs.size(), v.size());
      libmesh_assert_equal_to (row_dofs.size(), w.size());

      /* Constrain only v, not w.  */

      for (unsigned int i=0,
           n_row_dofs = cast_int<unsigned int>(row_dofs.size());
           i != n_row_dofs; i++)
        if (this->is_constrained_dof(row_dofs[i]))
          {
            // If the DOF is constrained
            libmesh_assert (_dof_constraints.find(row_dofs[i]) != _dof_constraints.end());

            v(i) = 0;
          }
    } // end if the RHS is constrained.
}



void DofMap::constrain_nothing (std::vector<dof_id_type> & dofs) const
{
  // check for easy return
  if (this->_dof_constraints.empty())
    return;

  // All the work is done by \p build_constraint_matrix.  We just need
  // a dummy matrix.
  DenseMatrix<Number> R;
  this->build_constraint_matrix (R, dofs);
}



void DofMap::enforce_constraints_exactly (const System & system,
                                          NumericVector<Number> * v,
                                          bool homogeneous) const
{
  parallel_object_only();

  if (!this->n_constrained_dofs())
    return;

  LOG_SCOPE("enforce_constraints_exactly()","DofMap");

  if (!v)
    v = system.solution.get();

  NumericVector<Number> * v_local  = nullptr; // will be initialized below
  NumericVector<Number> * v_global = nullptr; // will be initialized below
  std::unique_ptr<NumericVector<Number>> v_built;
  if (v->type() == SERIAL)
    {
      v_built = NumericVector<Number>::build(this->comm());
      v_built->init(this->n_dofs(), this->n_local_dofs(), true, PARALLEL);
      v_built->close();

      for (dof_id_type i=v_built->first_local_index();
           i<v_built->last_local_index(); i++)
        v_built->set(i, (*v)(i));
      v_built->close();
      v_global = v_built.get();

      v_local = v;
      libmesh_assert (v_local->closed());
    }
  else if (v->type() == PARALLEL)
    {
      v_built = NumericVector<Number>::build(this->comm());
      v_built->init (v->size(), v->local_size(),
                     this->get_send_list(), true,
                     GHOSTED);
      v->localize(*v_built, this->get_send_list());
      v_built->close();
      v_local = v_built.get();

      v_global = v;
    }
  else if (v->type() == GHOSTED)
    {
      v_local = v;
      v_global = v;
    }
  else // unknown v->type()
    libmesh_error_msg("ERROR: Unknown v->type() == " << v->type());

  // We should never hit these asserts because we should error-out in
  // else clause above.  Just to be sure we don't try to use v_local
  // and v_global uninitialized...
  libmesh_assert(v_local);
  libmesh_assert(v_global);
  libmesh_assert_equal_to (this, &(system.get_dof_map()));

  for (const auto & pr : _dof_constraints)
    {
      dof_id_type constrained_dof = pr.first;
      if (!this->local_index(constrained_dof))
        continue;

      const DofConstraintRow constraint_row = pr.second;

      Number exact_value = 0;
      if (!homogeneous)
        {
          DofConstraintValueMap::const_iterator rhsit =
            _primal_constraint_values.find(constrained_dof);
          if (rhsit != _primal_constraint_values.end())
            exact_value = rhsit->second;
        }
      for (const auto & j : constraint_row)
        exact_value += j.second * (*v_local)(j.first);

      v_global->set(constrained_dof, exact_value);
    }

  // If the old vector was serial, we probably need to send our values
  // to other processors
  if (v->type() == SERIAL)
    {
#ifndef NDEBUG
      v_global->close();
#endif
      v_global->localize (*v);
    }
  v->close();
}

void DofMap::enforce_constraints_on_residual (const NonlinearImplicitSystem & system,
                                              NumericVector<Number> * rhs,
                                              NumericVector<Number> const * solution,
                                              bool homogeneous) const
{
  parallel_object_only();

  if (!this->n_constrained_dofs())
    return;

  if (!rhs)
    rhs = system.rhs;
  if (!solution)
    solution = system.solution.get();

  NumericVector<Number> const * solution_local  = nullptr; // will be initialized below
  std::unique_ptr<NumericVector<Number>> solution_built;
  if (solution->type() == SERIAL || solution->type() == GHOSTED)
      solution_local = solution;
  else if (solution->type() == PARALLEL)
    {
      solution_built = NumericVector<Number>::build(this->comm());
      solution_built->init (solution->size(), solution->local_size(),
                            this->get_send_list(), true, GHOSTED);
      solution->localize(*solution_built, this->get_send_list());
      solution_built->close();
      solution_local = solution_built.get();
    }
  else // unknown solution->type()
    libmesh_error_msg("ERROR: Unknown solution->type() == " << solution->type());

  // We should never hit these asserts because we should error-out in
  // else clause above.  Just to be sure we don't try to use solution_local
  libmesh_assert(solution_local);
  libmesh_assert_equal_to (this, &(system.get_dof_map()));

  for (const auto & pr : _dof_constraints)
    {
      dof_id_type constrained_dof = pr.first;
      if (!this->local_index(constrained_dof))
        continue;

      const DofConstraintRow constraint_row = pr.second;

      Number exact_value = 0;
      for (const auto & j : constraint_row)
        exact_value -= j.second * (*solution_local)(j.first);
      exact_value += (*solution_local)(constrained_dof);
      if (!homogeneous)
        {
          DofConstraintValueMap::const_iterator rhsit =
            _primal_constraint_values.find(constrained_dof);
          if (rhsit != _primal_constraint_values.end())
            exact_value += rhsit->second;
        }

      rhs->set(constrained_dof, exact_value);
    }
}

void DofMap::enforce_constraints_on_jacobian (const NonlinearImplicitSystem & system,
                                              SparseMatrix<Number> * jac) const
{
  parallel_object_only();

  if (!this->n_constrained_dofs())
    return;

  if (!jac)
    jac = system.matrix;

  libmesh_assert_equal_to (this, &(system.get_dof_map()));

  for (const auto & pr : _dof_constraints)
    {
      dof_id_type constrained_dof = pr.first;
      if (!this->local_index(constrained_dof))
        continue;

      const DofConstraintRow constraint_row = pr.second;

      for (const auto & j : constraint_row)
        jac->set(constrained_dof, j.first, -j.second);
      jac->set(constrained_dof, constrained_dof, 1);
    }
}


void DofMap::enforce_adjoint_constraints_exactly (NumericVector<Number> & v,
                                                  unsigned int q) const
{
  parallel_object_only();

  if (!this->n_constrained_dofs())
    return;

  LOG_SCOPE("enforce_adjoint_constraints_exactly()", "DofMap");

  NumericVector<Number> * v_local  = nullptr; // will be initialized below
  NumericVector<Number> * v_global = nullptr; // will be initialized below
  std::unique_ptr<NumericVector<Number>> v_built;
  if (v.type() == SERIAL)
    {
      v_built = NumericVector<Number>::build(this->comm());
      v_built->init(this->n_dofs(), this->n_local_dofs(), true, PARALLEL);
      v_built->close();

      for (dof_id_type i=v_built->first_local_index();
           i<v_built->last_local_index(); i++)
        v_built->set(i, v(i));
      v_built->close();
      v_global = v_built.get();

      v_local = &v;
      libmesh_assert (v_local->closed());
    }
  else if (v.type() == PARALLEL)
    {
      v_built = NumericVector<Number>::build(this->comm());
      v_built->init (v.size(), v.local_size(),
                     this->get_send_list(), true, GHOSTED);
      v.localize(*v_built, this->get_send_list());
      v_built->close();
      v_local = v_built.get();

      v_global = &v;
    }
  else if (v.type() == GHOSTED)
    {
      v_local = &v;
      v_global = &v;
    }
  else // unknown v.type()
    libmesh_error_msg("ERROR: Unknown v.type() == " << v.type());

  // We should never hit these asserts because we should error-out in
  // else clause above.  Just to be sure we don't try to use v_local
  // and v_global uninitialized...
  libmesh_assert(v_local);
  libmesh_assert(v_global);

  // Do we have any non_homogeneous constraints?
  const AdjointDofConstraintValues::const_iterator
    adjoint_constraint_map_it = _adjoint_constraint_values.find(q);
  const DofConstraintValueMap * constraint_map =
    (adjoint_constraint_map_it == _adjoint_constraint_values.end()) ?
    nullptr : &adjoint_constraint_map_it->second;

  for (const auto & pr : _dof_constraints)
    {
      dof_id_type constrained_dof = pr.first;
      if (!this->local_index(constrained_dof))
        continue;

      const DofConstraintRow constraint_row = pr.second;

      Number exact_value = 0;
      if (constraint_map)
        {
          const DofConstraintValueMap::const_iterator
            adjoint_constraint_it =
            constraint_map->find(constrained_dof);
          if (adjoint_constraint_it != constraint_map->end())
            exact_value = adjoint_constraint_it->second;
        }

      for (const auto & j : constraint_row)
        exact_value += j.second * (*v_local)(j.first);

      v_global->set(constrained_dof, exact_value);
    }

  // If the old vector was serial, we probably need to send our values
  // to other processors
  if (v.type() == SERIAL)
    {
#ifndef NDEBUG
      v_global->close();
#endif
      v_global->localize (v);
    }
  v.close();
}



std::pair<Real, Real>
DofMap::max_constraint_error (const System & system,
                              NumericVector<Number> * v) const
{
  if (!v)
    v = system.solution.get();
  NumericVector<Number> & vec = *v;

  // We'll assume the vector is closed
  libmesh_assert (vec.closed());

  Real max_absolute_error = 0., max_relative_error = 0.;

  const MeshBase & mesh = system.get_mesh();

  libmesh_assert_equal_to (this, &(system.get_dof_map()));

  // indices on each element
  std::vector<dof_id_type> local_dof_indices;

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      this->dof_indices(elem, local_dof_indices);
      std::vector<dof_id_type> raw_dof_indices = local_dof_indices;

      // Constraint matrix for each element
      DenseMatrix<Number> C;

      this->build_constraint_matrix (C, local_dof_indices);

      // Continue if the element is unconstrained
      if (!C.m())
        continue;

      libmesh_assert_equal_to (C.m(), raw_dof_indices.size());
      libmesh_assert_equal_to (C.n(), local_dof_indices.size());

      for (auto i : make_range(C.m()))
        {
          // Recalculate any constrained dof owned by this processor
          dof_id_type global_dof = raw_dof_indices[i];
          if (this->is_constrained_dof(global_dof) &&
              global_dof >= vec.first_local_index() &&
              global_dof < vec.last_local_index())
            {
#ifndef NDEBUG
              DofConstraints::const_iterator
                pos = _dof_constraints.find(global_dof);

              libmesh_assert (pos != _dof_constraints.end());
#endif

              Number exact_value = 0;
              DofConstraintValueMap::const_iterator rhsit =
                _primal_constraint_values.find(global_dof);
              if (rhsit != _primal_constraint_values.end())
                exact_value = rhsit->second;

              for (auto j : make_range(C.n()))
                {
                  if (local_dof_indices[j] != global_dof)
                    exact_value += C(i,j) *
                      vec(local_dof_indices[j]);
                }

              max_absolute_error = std::max(max_absolute_error,
                                            std::abs(vec(global_dof) - exact_value));
              max_relative_error = std::max(max_relative_error,
                                            std::abs(vec(global_dof) - exact_value)
                                            / std::abs(exact_value));
            }
        }
    }

  return std::pair<Real, Real>(max_absolute_error, max_relative_error);
}



void DofMap::build_constraint_matrix (DenseMatrix<Number> & C,
                                      std::vector<dof_id_type> & elem_dofs,
                                      const bool called_recursively) const
{
  LOG_SCOPE_IF("build_constraint_matrix()", "DofMap", !called_recursively);

  // Create a set containing the DOFs we already depend on
  typedef std::set<dof_id_type> RCSet;
  RCSet dof_set;

  bool we_have_constraints = false;

  // Next insert any other dofs the current dofs might be constrained
  // in terms of.  Note that in this case we may not be done: Those
  // may in turn depend on others.  So, we need to repeat this process
  // in that case until the system depends only on unconstrained
  // degrees of freedom.
  for (const auto & dof : elem_dofs)
    if (this->is_constrained_dof(dof))
      {
        we_have_constraints = true;

        // If the DOF is constrained
        DofConstraints::const_iterator
          pos = _dof_constraints.find(dof);

        libmesh_assert (pos != _dof_constraints.end());

        const DofConstraintRow & constraint_row = pos->second;

        // Constraint rows in p refinement may be empty
        //libmesh_assert (!constraint_row.empty());

        for (const auto & item : constraint_row)
          dof_set.insert (item.first);
      }

  // May be safe to return at this point
  // (but remember to stop the perflog)
  if (!we_have_constraints)
    return;

  for (const auto & dof : elem_dofs)
    dof_set.erase (dof);

  // If we added any DOFS then we need to do this recursively.
  // It is possible that we just added a DOF that is also
  // constrained!
  //
  // Also, we need to handle the special case of an element having DOFs
  // constrained in terms of other, local DOFs
  if (!dof_set.empty() ||  // case 1: constrained in terms of other DOFs
      !called_recursively) // case 2: constrained in terms of our own DOFs
    {
      const unsigned int old_size =
        cast_int<unsigned int>(elem_dofs.size());

      // Add new dependency dofs to the end of the current dof set
      elem_dofs.insert(elem_dofs.end(),
                       dof_set.begin(), dof_set.end());

      // Now we can build the constraint matrix.
      // Note that resize also zeros for a DenseMatrix<Number>.
      C.resize (old_size,
                cast_int<unsigned int>(elem_dofs.size()));

      // Create the C constraint matrix.
      for (unsigned int i=0; i != old_size; i++)
        if (this->is_constrained_dof(elem_dofs[i]))
          {
            // If the DOF is constrained
            DofConstraints::const_iterator
              pos = _dof_constraints.find(elem_dofs[i]);

            libmesh_assert (pos != _dof_constraints.end());

            const DofConstraintRow & constraint_row = pos->second;

            // p refinement creates empty constraint rows
            //    libmesh_assert (!constraint_row.empty());

            for (const auto & item : constraint_row)
              for (unsigned int j=0,
                   n_elem_dofs = cast_int<unsigned int>(elem_dofs.size());
                   j != n_elem_dofs; j++)
                if (elem_dofs[j] == item.first)
                  C(i,j) = item.second;
          }
        else
          {
            C(i,i) = 1.;
          }

      // May need to do this recursively.  It is possible
      // that we just replaced a constrained DOF with another
      // constrained DOF.
      DenseMatrix<Number> Cnew;

      this->build_constraint_matrix (Cnew, elem_dofs, true);

      if ((C.n() == Cnew.m()) &&
          (Cnew.n() == elem_dofs.size())) // If the constraint matrix
        C.right_multiply(Cnew);           // is constrained...

      libmesh_assert_equal_to (C.n(), elem_dofs.size());
    }
}



void DofMap::build_constraint_matrix_and_vector (DenseMatrix<Number> & C,
                                                 DenseVector<Number> & H,
                                                 std::vector<dof_id_type> & elem_dofs,
                                                 int qoi_index,
                                                 const bool called_recursively) const
{
  LOG_SCOPE_IF("build_constraint_matrix_and_vector()", "DofMap", !called_recursively);

  // Create a set containing the DOFs we already depend on
  typedef std::set<dof_id_type> RCSet;
  RCSet dof_set;

  bool we_have_constraints = false;

  // Next insert any other dofs the current dofs might be constrained
  // in terms of.  Note that in this case we may not be done: Those
  // may in turn depend on others.  So, we need to repeat this process
  // in that case until the system depends only on unconstrained
  // degrees of freedom.
  for (const auto & dof : elem_dofs)
    if (this->is_constrained_dof(dof))
      {
        we_have_constraints = true;

        // If the DOF is constrained
        DofConstraints::const_iterator
          pos = _dof_constraints.find(dof);

        libmesh_assert (pos != _dof_constraints.end());

        const DofConstraintRow & constraint_row = pos->second;

        // Constraint rows in p refinement may be empty
        //libmesh_assert (!constraint_row.empty());

        for (const auto & item : constraint_row)
          dof_set.insert (item.first);
      }

  // May be safe to return at this point
  // (but remember to stop the perflog)
  if (!we_have_constraints)
    return;

  for (const auto & dof : elem_dofs)
    dof_set.erase (dof);

  // If we added any DOFS then we need to do this recursively.
  // It is possible that we just added a DOF that is also
  // constrained!
  //
  // Also, we need to handle the special case of an element having DOFs
  // constrained in terms of other, local DOFs
  if (!dof_set.empty() ||  // case 1: constrained in terms of other DOFs
      !called_recursively) // case 2: constrained in terms of our own DOFs
    {
      const DofConstraintValueMap * rhs_values = nullptr;
      if (qoi_index < 0)
        rhs_values = &_primal_constraint_values;
      else
        {
          const AdjointDofConstraintValues::const_iterator
            it = _adjoint_constraint_values.find(qoi_index);
          if (it != _adjoint_constraint_values.end())
            rhs_values = &it->second;
        }

      const unsigned int old_size =
        cast_int<unsigned int>(elem_dofs.size());

      // Add new dependency dofs to the end of the current dof set
      elem_dofs.insert(elem_dofs.end(),
                       dof_set.begin(), dof_set.end());

      // Now we can build the constraint matrix and vector.
      // Note that resize also zeros for a DenseMatrix and DenseVector
      C.resize (old_size,
                cast_int<unsigned int>(elem_dofs.size()));
      H.resize (old_size);

      // Create the C constraint matrix.
      for (unsigned int i=0; i != old_size; i++)
        if (this->is_constrained_dof(elem_dofs[i]))
          {
            // If the DOF is constrained
            DofConstraints::const_iterator
              pos = _dof_constraints.find(elem_dofs[i]);

            libmesh_assert (pos != _dof_constraints.end());

            const DofConstraintRow & constraint_row = pos->second;

            // p refinement creates empty constraint rows
            //    libmesh_assert (!constraint_row.empty());

            for (const auto & item : constraint_row)
              for (unsigned int j=0,
                   n_elem_dofs = cast_int<unsigned int>(elem_dofs.size());
                   j != n_elem_dofs; j++)
                if (elem_dofs[j] == item.first)
                  C(i,j) = item.second;

            if (rhs_values)
              {
                DofConstraintValueMap::const_iterator rhsit =
                  rhs_values->find(elem_dofs[i]);
                if (rhsit != rhs_values->end())
                  H(i) = rhsit->second;
              }
          }
        else
          {
            C(i,i) = 1.;
          }

      // May need to do this recursively.  It is possible
      // that we just replaced a constrained DOF with another
      // constrained DOF.
      DenseMatrix<Number> Cnew;
      DenseVector<Number> Hnew;

      this->build_constraint_matrix_and_vector (Cnew, Hnew, elem_dofs,
                                                qoi_index, true);

      if ((C.n() == Cnew.m()) &&          // If the constraint matrix
          (Cnew.n() == elem_dofs.size())) // is constrained...
        {
          // If x = Cy + h and y = Dz + g
          // Then x = (CD)z + (Cg + h)
          C.vector_mult_add(H, 1, Hnew);

          C.right_multiply(Cnew);
        }

      libmesh_assert_equal_to (C.n(), elem_dofs.size());
    }
}


void DofMap::allgather_recursive_constraints(MeshBase & mesh)
{
  // This function must be run on all processors at once
  parallel_object_only();

  // Return immediately if there's nothing to gather
  if (this->n_processors() == 1)
    return;

  // We might get to return immediately if none of the processors
  // found any constraints
  unsigned int has_constraints = !_dof_constraints.empty()
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
    || !_node_constraints.empty()
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS
    ;
  this->comm().max(has_constraints);
  if (!has_constraints)
    return;

  // If we have heterogenous adjoint constraints we need to
  // communicate those too.
  const unsigned int max_qoi_num =
    _adjoint_constraint_values.empty() ?
    0 : _adjoint_constraint_values.rbegin()->first+1;

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  // We may need to send nodes ahead of data about them
  std::vector<Parallel::Request> packed_range_sends;

  // We may be receiving packed_range sends out of order with
  // parallel_sync tags, so make sure they're received correctly.
  Parallel::MessageTag range_tag = this->comm().get_unique_tag();

  // We only need to do these sends on a distributed mesh
  const bool dist_mesh = !mesh.is_serial();
#endif

  // We might have calculated constraints for constrained dofs
  // which have support on other processors.
  // Push these out first.
  {
    std::map<processor_id_type, std::set<dof_id_type>> pushed_ids;

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
    std::map<processor_id_type, std::set<dof_id_type>> pushed_node_ids;
#endif

    const unsigned int sys_num = this->sys_number();

    // Collect the constraints to push to each processor
    for (auto & elem : as_range(mesh.active_not_local_elements_begin(),
                                mesh.active_not_local_elements_end()))
      {
        const unsigned short n_nodes = elem->n_nodes();

        // Just checking dof_indices on the foreign element isn't
        // enough.  Consider a central hanging node between a coarse
        // Q2/Q1 element and its finer neighbors on a higher-ranked
        // processor.  The coarse element's processor will own the node,
        // and will thereby own the pressure dof on that node, despite
        // the fact that that pressure dof doesn't directly exist on the
        // coarse element!
        //
        // So, we loop through dofs manually.

        {
          const unsigned int n_vars = elem->n_vars(sys_num);
          for (unsigned int v=0; v != n_vars; ++v)
            {
              const unsigned int n_comp = elem->n_comp(sys_num,v);
              for (unsigned int c=0; c != n_comp; ++c)
                {
                  const unsigned int id =
                    elem->dof_number(sys_num,v,c);
                  if (this->is_constrained_dof(id))
                    pushed_ids[elem->processor_id()].insert(id);
                }
            }
        }

        for (unsigned short n = 0; n != n_nodes; ++n)
          {
            const Node & node = elem->node_ref(n);
            const unsigned int n_vars = node.n_vars(sys_num);
            for (unsigned int v=0; v != n_vars; ++v)
              {
                const unsigned int n_comp = node.n_comp(sys_num,v);
                for (unsigned int c=0; c != n_comp; ++c)
                  {
                    const unsigned int id =
                      node.dof_number(sys_num,v,c);
                    if (this->is_constrained_dof(id))
                      pushed_ids[elem->processor_id()].insert(id);
                  }
              }
          }

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
        for (unsigned short n = 0; n != n_nodes; ++n)
          if (this->is_constrained_node(elem->node_ptr(n)))
            pushed_node_ids[elem->processor_id()].insert(elem->node_id(n));
#endif
      }

    // Rewrite those id sets as vectors for sending and receiving,
    // then find the corresponding data for each id, then push it all.
    std::map<processor_id_type, std::vector<dof_id_type>>
      pushed_id_vecs, received_id_vecs;
    for (auto & p : pushed_ids)
      pushed_id_vecs[p.first].assign(p.second.begin(), p.second.end());

    std::map<processor_id_type, std::vector<std::vector<std::pair<dof_id_type,Real>>>>
      pushed_keys_vals, received_keys_vals;
    std::map<processor_id_type, std::vector<std::vector<Number>>> pushed_rhss, received_rhss;
    for (auto & p : pushed_id_vecs)
      {
        auto & keys_vals = pushed_keys_vals[p.first];
        keys_vals.reserve(p.second.size());

        auto & rhss = pushed_rhss[p.first];
        rhss.reserve(p.second.size());
        for (auto & pushed_id : p.second)
          {
            const DofConstraintRow & row = _dof_constraints[pushed_id];
            keys_vals.emplace_back(row.begin(), row.end());

            rhss.push_back(std::vector<Number>(max_qoi_num+1));
            std::vector<Number> & rhs = rhss.back();
            DofConstraintValueMap::const_iterator rhsit =
              _primal_constraint_values.find(pushed_id);
            rhs[max_qoi_num] =
              (rhsit == _primal_constraint_values.end()) ?
              0 : rhsit->second;
            for (unsigned int q = 0; q != max_qoi_num; ++q)
              {
                AdjointDofConstraintValues::const_iterator adjoint_map_it =
                  _adjoint_constraint_values.find(q);

                if (adjoint_map_it == _adjoint_constraint_values.end())
                  continue;

                const DofConstraintValueMap & constraint_map =
                  adjoint_map_it->second;

                DofConstraintValueMap::const_iterator adj_rhsit =
                  constraint_map.find(pushed_id);

                rhs[q] =
                  (adj_rhsit == constraint_map.end()) ?
                  0 : adj_rhsit->second;
              }
          }
      }

    auto ids_action_functor =
      [& received_id_vecs]
      (processor_id_type pid,
       const std::vector<dof_id_type> & data)
      {
        received_id_vecs[pid] = data;
      };

    Parallel::push_parallel_vector_data
      (this->comm(), pushed_id_vecs, ids_action_functor);

    auto keys_vals_action_functor =
      [& received_keys_vals]
      (processor_id_type pid,
       const std::vector<std::vector<std::pair<dof_id_type,Real>>> & data)
      {
        received_keys_vals[pid] = data;
      };

    Parallel::push_parallel_vector_data
      (this->comm(), pushed_keys_vals, keys_vals_action_functor);

    auto rhss_action_functor =
      [& received_rhss]
      (processor_id_type pid,
       const std::vector<std::vector<Number>> & data)
      {
        received_rhss[pid] = data;
      };

    Parallel::push_parallel_vector_data
      (this->comm(), pushed_rhss, rhss_action_functor);

    // Now we have all the DofConstraint rows and rhs values received
    // from others, so add the DoF constraints that we've been sent

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
    std::map<processor_id_type, std::vector<dof_id_type>>
      pushed_node_id_vecs, received_node_id_vecs;
    for (auto & p : pushed_node_ids)
      pushed_node_id_vecs[p.first].assign(p.second.begin(), p.second.end());

    std::map<processor_id_type, std::vector<std::vector<std::pair<dof_id_type,Real>>>>
      pushed_node_keys_vals, received_node_keys_vals;
    std::map<processor_id_type, std::vector<Point>> pushed_offsets, received_offsets;

    for (auto & p : pushed_node_id_vecs)
      {
        const processor_id_type pid = p.first;

        // FIXME - this could be an unordered set, given a
        // hash<pointers> specialization
        std::set<const Node *> nodes_requested;

        auto & node_keys_vals = pushed_node_keys_vals[pid];
        node_keys_vals.reserve(p.second.size());

        auto & offsets = pushed_offsets[pid];
        offsets.reserve(p.second.size());

        for (auto & pushed_node_id : p.second)
          {
            const Node * node = mesh.node_ptr(pushed_node_id);
            NodeConstraintRow & row = _node_constraints[node].first;
            const std::size_t row_size = row.size();
            node_keys_vals.push_back
              (std::vector<std::pair<dof_id_type,Real>>());
            std::vector<std::pair<dof_id_type,Real>> & this_node_kv =
              node_keys_vals.back();
            this_node_kv.reserve(row_size);
            for (const auto & j : row)
              {
                this_node_kv.emplace_back(j.first->id(), j.second);

                // If we're not sure whether our send
                // destination already has this node, let's give
                // it a copy.
                if (j.first->processor_id() != pid && dist_mesh)
                  nodes_requested.insert(j.first);
              }

            offsets.push_back(_node_constraints[node].second);

          }

        // Constraining nodes might not even exist on our
        // correspondant's subset of a distributed mesh, so let's
        // make them exist.
        if (dist_mesh)
          {
            packed_range_sends.push_back(Parallel::Request());
            this->comm().send_packed_range
              (pid, &mesh, nodes_requested.begin(), nodes_requested.end(),
               packed_range_sends.back(), range_tag);
          }
      }

    auto node_ids_action_functor =
      [& received_node_id_vecs]
      (processor_id_type pid,
       const std::vector<dof_id_type> & data)
      {
        received_node_id_vecs[pid] = data;
      };

    Parallel::push_parallel_vector_data
      (this->comm(), pushed_node_id_vecs, node_ids_action_functor);

    auto node_keys_vals_action_functor =
      [& received_node_keys_vals]
      (processor_id_type pid,
       const std::vector<std::vector<std::pair<dof_id_type,Real>>> & data)
      {
        received_node_keys_vals[pid] = data;
      };

    Parallel::push_parallel_vector_data
      (this->comm(), pushed_node_keys_vals,
       node_keys_vals_action_functor);

    auto node_offsets_action_functor =
      [& received_offsets]
      (processor_id_type pid,
       const std::vector<Point> & data)
      {
        received_offsets[pid] = data;
      };

    Parallel::push_parallel_vector_data
      (this->comm(), pushed_offsets, node_offsets_action_functor);

#endif

    // Add all the dof constraints that I've been sent
    for (auto & p : received_id_vecs)
      {
        const processor_id_type pid = p.first;
        const auto & pushed_ids_to_me = p.second;
        libmesh_assert(received_keys_vals.count(pid));
        libmesh_assert(received_rhss.count(pid));
        const auto & pushed_keys_vals_to_me = received_keys_vals.at(pid);
        const auto & pushed_rhss_to_me = received_rhss.at(pid);

        libmesh_assert_equal_to (pushed_ids_to_me.size(),
                                 pushed_keys_vals_to_me.size());
        libmesh_assert_equal_to (pushed_ids_to_me.size(),
                                 pushed_rhss_to_me.size());

        for (auto i : index_range(pushed_ids_to_me))
          {
            dof_id_type constrained = pushed_ids_to_me[i];

            // If we don't already have a constraint for this dof,
            // add the one we were sent
            if (!this->is_constrained_dof(constrained))
              {
                DofConstraintRow & row = _dof_constraints[constrained];
                for (auto & kv : pushed_keys_vals_to_me[i])
                  {
                    libmesh_assert_less(kv.first, this->n_dofs());
                    row[kv.first] = kv.second;
                  }

                const Number primal_rhs = pushed_rhss_to_me[i][max_qoi_num];

                if (libmesh_isnan(primal_rhs))
                  libmesh_assert(pushed_keys_vals_to_me[i].empty());

                if (primal_rhs != Number(0))
                  _primal_constraint_values[constrained] = primal_rhs;
                else
                  _primal_constraint_values.erase(constrained);

                for (unsigned int q = 0; q != max_qoi_num; ++q)
                  {
                    AdjointDofConstraintValues::iterator adjoint_map_it =
                      _adjoint_constraint_values.find(q);

                    const Number adj_rhs = pushed_rhss_to_me[i][q];

                    if ((adjoint_map_it == _adjoint_constraint_values.end()) &&
                        adj_rhs == Number(0))
                      continue;

                    if (adjoint_map_it == _adjoint_constraint_values.end())
                      adjoint_map_it = _adjoint_constraint_values.emplace
                        (q, DofConstraintValueMap()).first;

                    DofConstraintValueMap & constraint_map =
                      adjoint_map_it->second;

                    if (adj_rhs != Number(0))
                      constraint_map[constrained] = adj_rhs;
                    else
                      constraint_map.erase(constrained);
                  }
              }
          }
      }

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
    // Add all the node constraints that I've been sent
    for (auto & p : received_node_id_vecs)
      {
        const processor_id_type pid = p.first;

        // Before we act on any new constraint rows, we may need to
        // make sure we have all the nodes involved!
        if (dist_mesh)
          this->comm().receive_packed_range
            (pid, &mesh, mesh_inserter_iterator<Node>(mesh),
             (Node**)nullptr, range_tag);

        const auto & pushed_node_ids_to_me = p.second;
        libmesh_assert(received_node_keys_vals.count(pid));
        libmesh_assert(received_offsets.count(pid));
        const auto & pushed_node_keys_vals_to_me = received_node_keys_vals.at(pid);
        const auto & pushed_offsets_to_me = received_offsets.at(pid);

        libmesh_assert_equal_to (pushed_node_ids_to_me.size(),
                                 pushed_node_keys_vals_to_me.size());
        libmesh_assert_equal_to (pushed_node_ids_to_me.size(),
                                 pushed_offsets_to_me.size());

        for (auto i : index_range(pushed_node_ids_to_me))
          {
            dof_id_type constrained_id = pushed_node_ids_to_me[i];

            // If we don't already have a constraint for this node,
            // add the one we were sent
            const Node * constrained = mesh.node_ptr(constrained_id);
            if (!this->is_constrained_node(constrained))
              {
                NodeConstraintRow & row = _node_constraints[constrained].first;
                for (auto & kv : pushed_node_keys_vals_to_me[i])
                  {
                    const Node * key_node = mesh.node_ptr(kv.first);
                    libmesh_assert(key_node);
                    row[key_node] = kv.second;
                  }
                _node_constraints[constrained].second = pushed_offsets_to_me[i];
              }
          }
      }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS
  }

  // Now start checking for any other constraints we need
  // to know about, requesting them recursively.

  // Create sets containing the DOFs and nodes we already depend on
  typedef std::set<dof_id_type> DoF_RCSet;
  DoF_RCSet unexpanded_dofs;

  for (const auto & i : _dof_constraints)
    unexpanded_dofs.insert(i.first);

  // Gather all the dof constraints we need
  this->gather_constraints(mesh, unexpanded_dofs, false);

  // Gather all the node constraints we need
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  typedef std::set<const Node *> Node_RCSet;
  Node_RCSet unexpanded_nodes;

  for (const auto & i : _node_constraints)
    unexpanded_nodes.insert(i.first);

  // We have to keep recursing while the unexpanded set is
  // nonempty on *any* processor
  bool unexpanded_set_nonempty = !unexpanded_nodes.empty();
  this->comm().max(unexpanded_set_nonempty);

  while (unexpanded_set_nonempty)
    {
      // Let's make sure we don't lose sync in this loop.
      parallel_object_only();

      // Request sets
      Node_RCSet node_request_set;

      // Request sets to send to each processor
      std::map<processor_id_type, std::vector<dof_id_type>>
        requested_node_ids;

      // And the sizes of each
      std::map<processor_id_type, dof_id_type> node_ids_on_proc;

      // Fill (and thereby sort and uniq!) the main request sets
      for (const auto & i : unexpanded_nodes)
        {
          NodeConstraintRow & row = _node_constraints[i].first;
          for (const auto & j : row)
            {
              const Node * const node = j.first;
              libmesh_assert(node);

              // If it's non-local and we haven't already got a
              // constraint for it, we might need to ask for one
              if ((node->processor_id() != this->processor_id()) &&
                  !_node_constraints.count(node))
                node_request_set.insert(node);
            }
        }

      // Clear the unexpanded constraint sets; we're about to expand
      // them
      unexpanded_nodes.clear();

      // Count requests by processor
      for (const auto & node : node_request_set)
        {
          libmesh_assert(node);
          libmesh_assert_less (node->processor_id(), this->n_processors());
          node_ids_on_proc[node->processor_id()]++;
        }

      for (auto pair : node_ids_on_proc)
        requested_node_ids[pair.first].reserve(pair.second);

      // Prepare each processor's request set
      for (const auto & node : node_request_set)
        requested_node_ids[node->processor_id()].push_back(node->id());

      typedef std::vector<std::pair<dof_id_type, Real>> row_datum;

      auto node_row_gather_functor =
        [this,
         & mesh,
         dist_mesh,
         & packed_range_sends,
         & range_tag]
        (processor_id_type pid,
         const std::vector<dof_id_type> & ids,
         std::vector<row_datum> & data)
        {
          // FIXME - this could be an unordered set, given a
          // hash<pointers> specialization
          std::set<const Node *> nodes_requested;

          // Fill those requests
          const std::size_t query_size = ids.size();

          data.resize(query_size);
          for (std::size_t i=0; i != query_size; ++i)
            {
              dof_id_type constrained_id = ids[i];
              const Node * constrained_node = mesh.node_ptr(constrained_id);
              if (_node_constraints.count(constrained_node))
                {
                  const NodeConstraintRow & row = _node_constraints[constrained_node].first;
                  std::size_t row_size = row.size();
                  data[i].reserve(row_size);
                  for (const auto & j : row)
                    {
                      const Node * node = j.first;
                      data[i].emplace_back(node->id(), j.second);

                      // If we're not sure whether our send
                      // destination already has this node, let's give
                      // it a copy.
                      if (node->processor_id() != pid && dist_mesh)
                        nodes_requested.insert(node);

                      // We can have 0 nodal constraint
                      // coefficients, where no Lagrange constraint
                      // exists but non-Lagrange basis constraints
                      // might.
                      // libmesh_assert(j.second);
                    }
                }
              else
                {
                  // We have to distinguish "constraint with no
                  // constraining nodes" (e.g. due to user node
                  // constraint equations) from "no constraint".
                  // We'll use invalid_id for the latter.
                  data[i].emplace_back(DofObject::invalid_id, Real(0));
                }
            }

          // Constraining nodes might not even exist on our
          // correspondant's subset of a distributed mesh, so let's
          // make them exist.
          if (dist_mesh)
            {
              packed_range_sends.push_back(Parallel::Request());
              this->comm().send_packed_range
                (pid, &mesh, nodes_requested.begin(), nodes_requested.end(),
                 packed_range_sends.back(), range_tag);
            }
        };

      typedef Point node_rhs_datum;

      auto node_rhs_gather_functor =
        [this,
         & mesh]
        (processor_id_type,
         const std::vector<dof_id_type> & ids,
         std::vector<node_rhs_datum> & data)
        {
          // Fill those requests
          const std::size_t query_size = ids.size();

          data.resize(query_size);
          for (std::size_t i=0; i != query_size; ++i)
            {
              dof_id_type constrained_id = ids[i];
              const Node * constrained_node = mesh.node_ptr(constrained_id);
              if (_node_constraints.count(constrained_node))
                data[i] = _node_constraints[constrained_node].second;
              else
                data[i](0) = std::numeric_limits<Real>::quiet_NaN();
            }
        };

      auto node_row_action_functor =
        [this,
         & mesh,
         dist_mesh,
         & range_tag,
         & unexpanded_nodes]
        (processor_id_type pid,
         const std::vector<dof_id_type> & ids,
         const std::vector<row_datum> & data)
        {
          // Before we act on any new constraint rows, we may need to
          // make sure we have all the nodes involved!
          if (dist_mesh)
            this->comm().receive_packed_range
              (pid, &mesh, mesh_inserter_iterator<Node>(mesh),
               (Node**)nullptr, range_tag);

          // Add any new constraint rows we've found
          const std::size_t query_size = ids.size();

          for (std::size_t i=0; i != query_size; ++i)
            {
              const dof_id_type constrained_id = ids[i];

              // An empty row is an constraint with an empty row; for
              // no constraint we use a "no row" placeholder
              if (data[i].empty())
                {
                  const Node * constrained_node = mesh.node_ptr(constrained_id);
                  NodeConstraintRow & row = _node_constraints[constrained_node].first;
                  row.clear();
                }
              else if (data[i][0].first != DofObject::invalid_id)
                {
                  const Node * constrained_node = mesh.node_ptr(constrained_id);
                  NodeConstraintRow & row = _node_constraints[constrained_node].first;
                  row.clear();
                  for (auto & pair : data[i])
                    {
                      const Node * key_node =
                        mesh.node_ptr(pair.first);
                      libmesh_assert(key_node);
                      row[key_node] = pair.second;
                    }

                  // And prepare to check for more recursive constraints
                  unexpanded_nodes.insert(constrained_node);
                }
            }
        };

      auto node_rhs_action_functor =
        [this,
         & mesh]
        (processor_id_type,
         const std::vector<dof_id_type> & ids,
         const std::vector<node_rhs_datum> & data)
        {
          // Add rhs data for any new node constraint rows we've found
          const std::size_t query_size = ids.size();

          for (std::size_t i=0; i != query_size; ++i)
            {
              dof_id_type constrained_id = ids[i];
              const Node * constrained_node = mesh.node_ptr(constrained_id);

              if (!libmesh_isnan(data[i](0)))
                _node_constraints[constrained_node].second = data[i];
              else
                _node_constraints.erase(constrained_node);
            }
        };

      // Now request node constraint rows from other processors
      row_datum * node_row_ex = nullptr;
      Parallel::pull_parallel_vector_data
        (this->comm(), requested_node_ids, node_row_gather_functor,
         node_row_action_functor, node_row_ex);

      // And request node constraint right hand sides from other procesors
      node_rhs_datum * node_rhs_ex = nullptr;
      Parallel::pull_parallel_vector_data
        (this->comm(), requested_node_ids, node_rhs_gather_functor,
         node_rhs_action_functor, node_rhs_ex);


      // We have to keep recursing while the unexpanded set is
      // nonempty on *any* processor
      unexpanded_set_nonempty = !unexpanded_nodes.empty();
      this->comm().max(unexpanded_set_nonempty);
    }
  Parallel::wait(packed_range_sends);
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS
}



void DofMap::process_constraints (MeshBase & mesh)
{
  // We've computed our local constraints, but they may depend on
  // non-local constraints that we'll need to take into account.
  this->allgather_recursive_constraints(mesh);

  if (_error_on_constraint_loop)
  {
    // Optionally check for constraint loops and throw an error
    // if they're detected. We always do this check below in dbg/devel
    // mode but here we optionally do it in opt mode as well.
    check_for_constraint_loops();
  }

  // Adjoints will be constrained where the primal is
  // Therefore, we will expand the adjoint_constraint_values
  // map whenever the primal_constraint_values map is expanded

  // First, figure out the total number of QoIs
  const unsigned int max_qoi_num =
    _adjoint_constraint_values.empty() ?
    0 : _adjoint_constraint_values.rbegin()->first+1;

  // Create a set containing the DOFs we already depend on
  typedef std::set<dof_id_type> RCSet;
  RCSet unexpanded_set;

  for (const auto & i : _dof_constraints)
    unexpanded_set.insert(i.first);

  while (!unexpanded_set.empty())
    for (RCSet::iterator i = unexpanded_set.begin();
         i != unexpanded_set.end(); /* nothing */)
      {
        // If the DOF is constrained
        DofConstraints::iterator
          pos = _dof_constraints.find(*i);

        libmesh_assert (pos != _dof_constraints.end());

        DofConstraintRow & constraint_row = pos->second;

        DofConstraintValueMap::iterator rhsit =
          _primal_constraint_values.find(*i);
        Number constraint_rhs = (rhsit == _primal_constraint_values.end()) ?
          0 : rhsit->second;

        // A vector of DofConstraintValueMaps for each adjoint variable
        std::vector<DofConstraintValueMap::iterator> adjoint_rhs_iterators;
        adjoint_rhs_iterators.resize(max_qoi_num);

        // Another to hold the adjoint constraint rhs
        std::vector<Number> adjoint_constraint_rhs(max_qoi_num, 0.0);

        // Find and gather recursive constraints for each adjoint variable
        for (auto & adjoint_map : _adjoint_constraint_values)
          {
            const std::size_t q = adjoint_map.first;
            adjoint_rhs_iterators[q] = adjoint_map.second.find(*i);

            adjoint_constraint_rhs[q] =
              (adjoint_rhs_iterators[q] == adjoint_map.second.end()) ?
              0 : adjoint_rhs_iterators[q]->second;
          }

        std::vector<dof_id_type> constraints_to_expand;

        for (const auto & item : constraint_row)
          if (item.first != *i && this->is_constrained_dof(item.first))
            {
              unexpanded_set.insert(item.first);
              constraints_to_expand.push_back(item.first);
            }

        for (const auto & expandable : constraints_to_expand)
          {
            const Real this_coef = constraint_row[expandable];

            DofConstraints::const_iterator
              subpos = _dof_constraints.find(expandable);

            libmesh_assert (subpos != _dof_constraints.end());

            const DofConstraintRow & subconstraint_row = subpos->second;

            for (const auto & item : subconstraint_row)
              {
                // Assert that the constraint does not form a cycle.
                libmesh_assert(item.first != expandable);
                constraint_row[item.first] += item.second * this_coef;
              }

            DofConstraintValueMap::const_iterator subrhsit =
              _primal_constraint_values.find(expandable);
            if (subrhsit != _primal_constraint_values.end())
              constraint_rhs += subrhsit->second * this_coef;

            // Find and gather recursive constraints for each adjoint variable
            for (const auto & adjoint_map : _adjoint_constraint_values)
              {
                const std::size_t q = adjoint_map.first;

                DofConstraintValueMap::const_iterator adjoint_subrhsit =
                  adjoint_map.second.find(expandable);

                if (adjoint_subrhsit != adjoint_map.second.end())
                adjoint_constraint_rhs[q] += adjoint_subrhsit->second * this_coef;
              }

            constraint_row.erase(expandable);
          }

        if (rhsit == _primal_constraint_values.end())
          {
            if (constraint_rhs != Number(0))
              _primal_constraint_values[*i] = constraint_rhs;
            else
              _primal_constraint_values.erase(*i);
          }
        else
          {
            if (constraint_rhs != Number(0))
              rhsit->second = constraint_rhs;
            else
              _primal_constraint_values.erase(rhsit);
          }

        // Finally fill in the adjoint constraints for each adjoint variable if possible
        for (auto & adjoint_map : _adjoint_constraint_values)
          {
            const std::size_t q = adjoint_map.first;

            if(adjoint_rhs_iterators[q] == adjoint_map.second.end())
              {
                if (adjoint_constraint_rhs[q] != Number(0))
                   (adjoint_map.second)[*i] = adjoint_constraint_rhs[q];
                else
                  adjoint_map.second.erase(*i);
              }
            else
              {
                if (adjoint_constraint_rhs[q] != Number(0))
                  adjoint_rhs_iterators[q]->second = adjoint_constraint_rhs[q];
                else
                  adjoint_map.second.erase(adjoint_rhs_iterators[q]);
              }
          }

        if (constraints_to_expand.empty())
          i = unexpanded_set.erase(i);
        else
          ++i;
      }

  // In parallel we can't guarantee that nodes/dofs which constrain
  // others are on processors which are aware of that constraint, yet
  // we need such awareness for sparsity pattern generation.  So send
  // other processors any constraints they might need to know about.
  this->scatter_constraints(mesh);

  // Now that we have our root constraint dependencies sorted out, add
  // them to the send_list
  this->add_constraints_to_send_list();
}


#ifdef LIBMESH_ENABLE_CONSTRAINTS
void DofMap::check_for_cyclic_constraints()
{
  // Eventually make this officially libmesh_deprecated();
  check_for_constraint_loops();
}

void DofMap::check_for_constraint_loops()
{
  // Create a set containing the DOFs we already depend on
  typedef std::set<dof_id_type> RCSet;
  RCSet unexpanded_set;

  // Use dof_constraints_copy in this method so that we don't
  // mess with _dof_constraints.
  DofConstraints dof_constraints_copy = _dof_constraints;

  for (const auto & i : dof_constraints_copy)
    unexpanded_set.insert(i.first);

  while (!unexpanded_set.empty())
    for (RCSet::iterator i = unexpanded_set.begin();
         i != unexpanded_set.end(); /* nothing */)
      {
        // If the DOF is constrained
        DofConstraints::iterator
          pos = dof_constraints_copy.find(*i);

        libmesh_assert (pos != dof_constraints_copy.end());

        DofConstraintRow & constraint_row = pos->second;

        // Comment out "rhs" parts of this method copied from process_constraints
        // DofConstraintValueMap::iterator rhsit =
        //   _primal_constraint_values.find(*i);
        // Number constraint_rhs = (rhsit == _primal_constraint_values.end()) ?
        //   0 : rhsit->second;

        std::vector<dof_id_type> constraints_to_expand;

        for (const auto & item : constraint_row)
          if (item.first != *i && this->is_constrained_dof(item.first))
            {
              unexpanded_set.insert(item.first);
              constraints_to_expand.push_back(item.first);
            }

        for (const auto & expandable : constraints_to_expand)
          {
            const Real this_coef = constraint_row[expandable];

            DofConstraints::const_iterator
              subpos = dof_constraints_copy.find(expandable);

            libmesh_assert (subpos != dof_constraints_copy.end());

            const DofConstraintRow & subconstraint_row = subpos->second;

            for (const auto & item : subconstraint_row)
              {
                libmesh_error_msg_if(item.first == expandable, "Constraint loop detected");

                constraint_row[item.first] += item.second * this_coef;
              }

            // Comment out "rhs" parts of this method copied from process_constraints
            // DofConstraintValueMap::const_iterator subrhsit =
            //   _primal_constraint_values.find(expandable);
            // if (subrhsit != _primal_constraint_values.end())
            //   constraint_rhs += subrhsit->second * this_coef;

            constraint_row.erase(expandable);
          }

        // Comment out "rhs" parts of this method copied from process_constraints
        // if (rhsit == _primal_constraint_values.end())
        //   {
        //     if (constraint_rhs != Number(0))
        //       _primal_constraint_values[*i] = constraint_rhs;
        //     else
        //       _primal_constraint_values.erase(*i);
        //   }
        // else
        //   {
        //     if (constraint_rhs != Number(0))
        //       rhsit->second = constraint_rhs;
        //     else
        //       _primal_constraint_values.erase(rhsit);
        //   }

        if (constraints_to_expand.empty())
          i = unexpanded_set.erase(i);
        else
          ++i;
      }
}
#else
void DofMap::check_for_constraint_loops() {}
void DofMap::check_for_cyclic_constraints()
{
  // Do nothing
}
#endif


void DofMap::scatter_constraints(MeshBase & mesh)
{
  // At this point each processor with a constrained node knows
  // the corresponding constraint row, but we also need each processor
  // with a constrainer node to know the corresponding row(s).

  // This function must be run on all processors at once
  parallel_object_only();

  // Return immediately if there's nothing to gather
  if (this->n_processors() == 1)
    return;

  // We might get to return immediately if none of the processors
  // found any constraints
  unsigned int has_constraints = !_dof_constraints.empty()
#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
    || !_node_constraints.empty()
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS
    ;
  this->comm().max(has_constraints);
  if (!has_constraints)
    return;

  // We may be receiving packed_range sends out of order with
  // parallel_sync tags, so make sure they're received correctly.
  Parallel::MessageTag range_tag = this->comm().get_unique_tag();

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  std::map<processor_id_type, std::set<dof_id_type>> pushed_node_ids;
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

  std::map<processor_id_type, std::set<dof_id_type>> pushed_ids;

  // Collect the dof constraints I need to push to each processor
  dof_id_type constrained_proc_id = 0;
  for (auto & i : _dof_constraints)
    {
      const dof_id_type constrained = i.first;
      while (constrained >= _end_df[constrained_proc_id])
        constrained_proc_id++;

      if (constrained_proc_id != this->processor_id())
        continue;

      DofConstraintRow & row = i.second;
      for (auto & j : row)
        {
          const dof_id_type constraining = j.first;

          processor_id_type constraining_proc_id = 0;
          while (constraining >= _end_df[constraining_proc_id])
            constraining_proc_id++;

          if (constraining_proc_id != this->processor_id() &&
              constraining_proc_id != constrained_proc_id)
            pushed_ids[constraining_proc_id].insert(constrained);
        }
    }

  // Pack the dof constraint rows and rhs's to push

  std::map<processor_id_type,
          std::vector<std::vector<std::pair<dof_id_type, Real>>>>
    pushed_keys_vals, pushed_keys_vals_to_me;

  std::map<processor_id_type, std::vector<std::pair<dof_id_type, Number>>>
    pushed_ids_rhss, pushed_ids_rhss_to_me;

  auto gather_ids =
    [this,
     & pushed_ids,
     & pushed_keys_vals,
     & pushed_ids_rhss]
    ()
    {
      for (auto & pid_id_pair : pushed_ids)
        {
          const processor_id_type pid = pid_id_pair.first;
          const std::set<dof_id_type>
            & pid_ids = pid_id_pair.second;

          const std::size_t ids_size = pid_ids.size();
          std::vector<std::vector<std::pair<dof_id_type, Real>>> &
            keys_vals = pushed_keys_vals[pid];
          std::vector<std::pair<dof_id_type,Number>> &
            ids_rhss = pushed_ids_rhss[pid];
          keys_vals.resize(ids_size);
          ids_rhss.resize(ids_size);

          std::size_t push_i;
          std::set<dof_id_type>::const_iterator it;
          for (push_i = 0, it = pid_ids.begin();
               it != pid_ids.end(); ++push_i, ++it)
            {
              const dof_id_type constrained = *it;
              DofConstraintRow & row = _dof_constraints[constrained];
              keys_vals[push_i].assign(row.begin(), row.end());

              DofConstraintValueMap::const_iterator rhsit =
                _primal_constraint_values.find(constrained);
              ids_rhss[push_i].first = constrained;
              ids_rhss[push_i].second =
                (rhsit == _primal_constraint_values.end()) ?
                0 : rhsit->second;
            }
        }
    };

  gather_ids();

  auto ids_rhss_action_functor =
    [& pushed_ids_rhss_to_me]
    (processor_id_type pid,
     const std::vector<std::pair<dof_id_type, Number>> & data)
    {
      pushed_ids_rhss_to_me[pid] = data;
    };

  auto keys_vals_action_functor =
    [& pushed_keys_vals_to_me]
    (processor_id_type pid,
     const std::vector<std::vector<std::pair<dof_id_type, Real>>> & data)
    {
      pushed_keys_vals_to_me[pid] = data;
    };

  Parallel::push_parallel_vector_data
    (this->comm(), pushed_ids_rhss, ids_rhss_action_functor);
  Parallel::push_parallel_vector_data
    (this->comm(), pushed_keys_vals, keys_vals_action_functor);

  // Now work on traded dof constraint rows
  auto receive_dof_constraints =
    [this,
     & pushed_ids_rhss_to_me,
     & pushed_keys_vals_to_me]
    ()
    {
      for (auto & pid_id_pair : pushed_ids_rhss_to_me)
        {
          const processor_id_type pid = pid_id_pair.first;
          const auto & ids_rhss = pid_id_pair.second;
          const auto & keys_vals = pushed_keys_vals_to_me[pid];

          libmesh_assert_equal_to
            (ids_rhss.size(), keys_vals.size());

          // Add the dof constraints that I've been sent
          for (auto i : index_range(ids_rhss))
            {
              dof_id_type constrained = ids_rhss[i].first;

              // If we don't already have a constraint for this dof,
              // add the one we were sent
              if (!this->is_constrained_dof(constrained))
                {
                  DofConstraintRow & row = _dof_constraints[constrained];
                  for (auto & key_val : keys_vals[i])
                    {
                      libmesh_assert_less(key_val.first, this->n_dofs());
                      row[key_val.first] = key_val.second;
                    }
                  if (ids_rhss[i].second != Number(0))
                    _primal_constraint_values[constrained] =
                      ids_rhss[i].second;
                  else
                    _primal_constraint_values.erase(constrained);
                }
            }
        }
    };

  receive_dof_constraints();

#ifdef LIBMESH_ENABLE_NODE_CONSTRAINTS
  // Collect the node constraints to push to each processor
  for (auto & i : _node_constraints)
    {
      const Node * constrained = i.first;

      if (constrained->processor_id() != this->processor_id())
        continue;

      NodeConstraintRow & row = i.second.first;
      for (auto & j : row)
        {
          const Node * constraining = j.first;

          if (constraining->processor_id() != this->processor_id() &&
              constraining->processor_id() != constrained->processor_id())
            pushed_node_ids[constraining->processor_id()].insert(constrained->id());
        }
    }

  // Pack the node constraint rows and rhss to push
  std::map<processor_id_type,
          std::vector<std::vector<std::pair<dof_id_type,Real>>>>
    pushed_node_keys_vals, pushed_node_keys_vals_to_me;
  std::map<processor_id_type, std::vector<std::pair<dof_id_type, Point>>>
    pushed_node_ids_offsets, pushed_node_ids_offsets_to_me;
  std::map<processor_id_type, std::set<Node *>> pushed_node_sets;

  for (auto & pid_id_pair : pushed_node_ids)
    {
      const processor_id_type pid = pid_id_pair.first;
      const std::set<dof_id_type>
        & pid_ids = pid_id_pair.second;

      const std::size_t ids_size = pid_ids.size();
      std::vector<std::vector<std::pair<dof_id_type,Real>>> &
        keys_vals = pushed_node_keys_vals[pid];
      std::vector<std::pair<dof_id_type, Point>> &
        ids_offsets = pushed_node_ids_offsets[pid];
      keys_vals.resize(ids_size);
      ids_offsets.resize(ids_size);
      std::set<Node *> nodes;

      std::size_t push_i;
      std::set<dof_id_type>::const_iterator it;
      for (push_i = 0, it = pid_ids.begin();
           it != pid_ids.end(); ++push_i, ++it)
        {
          Node * constrained = mesh.node_ptr(*it);

          if (constrained->processor_id() != pid)
            nodes.insert(constrained);

          NodeConstraintRow & row = _node_constraints[constrained].first;
          std::size_t row_size = row.size();
          keys_vals[push_i].reserve(row_size);
          for (const auto & j : row)
            {
              Node * constraining = const_cast<Node *>(j.first);

              keys_vals[push_i].emplace_back(constraining->id(), j.second);

              if (constraining->processor_id() != pid)
                nodes.insert(constraining);
            }

          ids_offsets[push_i].first = *it;
          ids_offsets[push_i].second = _node_constraints[constrained].second;
        }

      if (!mesh.is_serial())
        {
          auto & pid_nodes = pushed_node_sets[pid];
          pid_nodes.insert(nodes.begin(), nodes.end());
        }
    }

  auto node_ids_offsets_action_functor =
    [& pushed_node_ids_offsets_to_me]
    (processor_id_type pid,
     const std::vector<std::pair<dof_id_type, Point>> & data)
    {
      pushed_node_ids_offsets_to_me[pid] = data;
    };

  auto node_keys_vals_action_functor =
    [& pushed_node_keys_vals_to_me]
    (processor_id_type pid,
     const std::vector<std::vector<std::pair<dof_id_type, Real>>> & data)
    {
      pushed_node_keys_vals_to_me[pid] = data;
    };

  // Trade pushed node constraint rows
  Parallel::push_parallel_vector_data
    (this->comm(), pushed_node_ids_offsets, node_ids_offsets_action_functor);
  Parallel::push_parallel_vector_data
    (this->comm(), pushed_node_keys_vals, node_keys_vals_action_functor);

  // Constraining nodes might not even exist on our subset of a
  // distributed mesh, so let's make them exist.
  auto insert_nodes_functor =
    [& mesh]
    (processor_id_type /* pid */,
     const std::set<Node *> & nodes)
    {
      for (Node * node : nodes)
        mesh.add_node(node);
    };

  if (!mesh.is_serial())
    Parallel::push_parallel_packed_range
      (this->comm(), pushed_node_sets, &mesh, insert_nodes_functor);

  for (auto & pid_id_pair : pushed_node_ids_offsets_to_me)
    {
      const processor_id_type pid = pid_id_pair.first;
      const auto & ids_offsets = pid_id_pair.second;
      const auto & keys_vals = pushed_node_keys_vals_to_me[pid];

      libmesh_assert_equal_to
        (ids_offsets.size(), keys_vals.size());

      // Add the node constraints that I've been sent
      for (auto i : index_range(ids_offsets))
        {
          dof_id_type constrained_id = ids_offsets[i].first;

          // If we don't already have a constraint for this node,
          // add the one we were sent
          const Node * constrained = mesh.node_ptr(constrained_id);
          if (!this->is_constrained_node(constrained))
            {
              NodeConstraintRow & row = _node_constraints[constrained].first;
              for (auto & key_val : keys_vals[i])
                {
                  const Node * key_node = mesh.node_ptr(key_val.first);
                  row[key_node] = key_val.second;
                }
              _node_constraints[constrained].second =
                ids_offsets[i].second;
            }
        }
    }
#endif // LIBMESH_ENABLE_NODE_CONSTRAINTS

  // Next we need to push constraints to processors which don't own
  // the constrained dof, don't own the constraining dof, but own an
  // element supporting the constraining dof.
  //
  // We need to be able to quickly look up constrained dof ids by what
  // constrains them, so that we can handle the case where we see a
  // foreign element containing one of our constraining DoF ids and we
  // need to push that constraint.
  //
  // Getting distributed adaptive sparsity patterns right is hard.

  typedef std::map<dof_id_type, std::set<dof_id_type>> DofConstrainsMap;
  DofConstrainsMap dof_id_constrains;

  for (auto & i : _dof_constraints)
    {
      const dof_id_type constrained = i.first;
      DofConstraintRow & row = i.second;
      for (const auto & j : row)
        {
          const dof_id_type constraining = j.first;

          dof_id_type constraining_proc_id = 0;
          while (constraining >= _end_df[constraining_proc_id])
            constraining_proc_id++;

          if (constraining_proc_id == this->processor_id())
            dof_id_constrains[constraining].insert(constrained);
        }
    }

  // Loop over all foreign elements, find any supporting our
  // constrained dof indices.
  pushed_ids.clear();

  for (const auto & elem : as_range(mesh.active_not_local_elements_begin(),
                                    mesh.active_not_local_elements_end()))
    {
      std::vector<dof_id_type> my_dof_indices;
      this->dof_indices (elem, my_dof_indices);

      for (const auto & dof : my_dof_indices)
        {
          DofConstrainsMap::const_iterator dcmi = dof_id_constrains.find(dof);
          if (dcmi != dof_id_constrains.end())
            {
              for (const auto & constrained : dcmi->second)
                {
                  dof_id_type the_constrained_proc_id = 0;
                  while (constrained >= _end_df[the_constrained_proc_id])
                    the_constrained_proc_id++;

                  const processor_id_type elemproc = elem->processor_id();
                  if (elemproc != the_constrained_proc_id)
                    pushed_ids[elemproc].insert(constrained);
                }
            }
        }
    }

  pushed_ids_rhss.clear();
  pushed_ids_rhss_to_me.clear();
  pushed_keys_vals.clear();
  pushed_keys_vals_to_me.clear();

  gather_ids();

  // Trade pushed dof constraint rows
  Parallel::push_parallel_vector_data
    (this->comm(), pushed_ids_rhss, ids_rhss_action_functor);
  Parallel::push_parallel_vector_data
    (this->comm(), pushed_keys_vals, keys_vals_action_functor);

  receive_dof_constraints();

  // Finally, we need to handle the case of remote dof coupling.  If a
  // processor's element is coupled to a ghost element, then the
  // processor needs to know about all constraints which affect the
  // dofs on that ghost element, so we'll have to query the ghost
  // element's owner.

  GhostingFunctor::map_type elements_to_couple;

  // Man, I wish we had guaranteed unique_ptr availability...
  std::set<CouplingMatrix*> temporary_coupling_matrices;

  this->merge_ghost_functor_outputs
    (elements_to_couple,
     temporary_coupling_matrices,
     this->coupling_functors_begin(),
     this->coupling_functors_end(),
     mesh.active_local_elements_begin(),
     mesh.active_local_elements_end(),
     this->processor_id());

  // Each ghost-coupled element's owner should get a request for its dofs
  std::set<dof_id_type> requested_dofs;

  for (const auto & pr : elements_to_couple)
    {
      const Elem * elem = pr.first;

      // FIXME - optimize for the non-fully-coupled case?
      std::vector<dof_id_type> element_dofs;
      this->dof_indices(elem, element_dofs);

      for (auto dof : element_dofs)
        requested_dofs.insert(dof);
    }

  this->gather_constraints(mesh, requested_dofs, false);
}


void DofMap::gather_constraints (MeshBase & /*mesh*/,
                                 std::set<dof_id_type> & unexpanded_dofs,
                                 bool /*look_for_constrainees*/)
{
  typedef std::set<dof_id_type> DoF_RCSet;

  // If we have heterogenous adjoint constraints we need to
  // communicate those too.
  const unsigned int max_qoi_num =
    _adjoint_constraint_values.empty() ?
    0 : _adjoint_constraint_values.rbegin()->first+1;

  // We have to keep recursing while the unexpanded set is
  // nonempty on *any* processor
  bool unexpanded_set_nonempty = !unexpanded_dofs.empty();
  this->comm().max(unexpanded_set_nonempty);

  while (unexpanded_set_nonempty)
    {
      // Let's make sure we don't lose sync in this loop.
      parallel_object_only();

      // Request sets
      DoF_RCSet   dof_request_set;

      // Request sets to send to each processor
      std::map<processor_id_type, std::vector<dof_id_type>>
        requested_dof_ids;

      // And the sizes of each
      std::map<processor_id_type, dof_id_type>
        dof_ids_on_proc;

      // Fill (and thereby sort and uniq!) the main request sets
      for (const auto & unexpanded_dof : unexpanded_dofs)
        {
          DofConstraints::const_iterator
            pos = _dof_constraints.find(unexpanded_dof);

          // If we were asked for a DoF and we don't already have a
          // constraint for it, then we need to check for one.
          if (pos == _dof_constraints.end())
            {
              if (!this->local_index(unexpanded_dof) &&
                  !_dof_constraints.count(unexpanded_dof) )
                dof_request_set.insert(unexpanded_dof);
            }
          // If we were asked for a DoF and we already have a
          // constraint for it, then we need to check if the
          // constraint is recursive.
          else
            {
              const DofConstraintRow & row = pos->second;
              for (const auto & j : row)
                {
                  const dof_id_type constraining_dof = j.first;

                  // If it's non-local and we haven't already got a
                  // constraint for it, we might need to ask for one
                  if (!this->local_index(constraining_dof) &&
                      !_dof_constraints.count(constraining_dof))
                    dof_request_set.insert(constraining_dof);
                }
            }
        }

      // Clear the unexpanded constraint set; we're about to expand it
      unexpanded_dofs.clear();

      // Count requests by processor
      processor_id_type proc_id = 0;
      for (const auto & i : dof_request_set)
        {
          while (i >= _end_df[proc_id])
            proc_id++;
          dof_ids_on_proc[proc_id]++;
        }

      for (auto & pair : dof_ids_on_proc)
        {
          requested_dof_ids[pair.first].reserve(pair.second);
        }

      // Prepare each processor's request set
      proc_id = 0;
      for (const auto & i : dof_request_set)
        {
          while (i >= _end_df[proc_id])
            proc_id++;
          requested_dof_ids[proc_id].push_back(i);
        }

      typedef std::vector<std::pair<dof_id_type, Real>> row_datum;

      typedef std::vector<Number> rhss_datum;

      auto row_gather_functor =
        [this]
        (processor_id_type,
         const std::vector<dof_id_type> & ids,
         std::vector<row_datum> & data)
        {
          // Fill those requests
          const std::size_t query_size = ids.size();

          data.resize(query_size);
          for (std::size_t i=0; i != query_size; ++i)
            {
              dof_id_type constrained = ids[i];
              if (_dof_constraints.count(constrained))
                {
                  DofConstraintRow & row = _dof_constraints[constrained];
                  std::size_t row_size = row.size();
                  data[i].reserve(row_size);
                  for (const auto & j : row)
                    {
                      data[i].push_back(j);

                      // We should never have an invalid constraining
                      // dof id
                      libmesh_assert(j.first != DofObject::invalid_id);

                      // We should never have a 0 constraint
                      // coefficient; that's implicit via sparse
                      // constraint storage
                      //
                      // But we can't easily control how users add
                      // constraints, so we can't safely assert that
                      // we're being efficient here.
                      //
                      // libmesh_assert(j.second);
                    }
                }
              else
                {
                  // We have to distinguish "constraint with no
                  // constraining dofs" (e.g. due to Dirichlet
                  // constraint equations) from "no constraint".
                  // We'll use invalid_id for the latter.
                  data[i].emplace_back(DofObject::invalid_id, Real(0));
                }
            }
        };

      auto rhss_gather_functor =
        [this,
         max_qoi_num]
        (processor_id_type,
         const std::vector<dof_id_type> & ids,
         std::vector<rhss_datum> & data)
        {
          // Fill those requests
          const std::size_t query_size = ids.size();

          data.resize(query_size);
          for (std::size_t i=0; i != query_size; ++i)
            {
              dof_id_type constrained = ids[i];
              data[i].clear();
              if (_dof_constraints.count(constrained))
                {
                  DofConstraintValueMap::const_iterator rhsit =
                    _primal_constraint_values.find(constrained);
                  data[i].push_back
                    ((rhsit == _primal_constraint_values.end()) ?
                     0 : rhsit->second);

                  for (unsigned int q = 0; q != max_qoi_num; ++q)
                    {
                      AdjointDofConstraintValues::const_iterator adjoint_map_it =
                        _adjoint_constraint_values.find(q);

                      if (adjoint_map_it == _adjoint_constraint_values.end())
                        {
                          data[i].push_back(0);
                          continue;
                        }

                      const DofConstraintValueMap & constraint_map =
                        adjoint_map_it->second;

                      DofConstraintValueMap::const_iterator adj_rhsit =
                        constraint_map.find(constrained);
                      data[i].push_back
                        ((adj_rhsit == constraint_map.end()) ?
                         0 : adj_rhsit->second);
                    }
                }
            }
        };

      auto row_action_functor =
        [this,
         & unexpanded_dofs]
        (processor_id_type,
         const std::vector<dof_id_type> & ids,
         const std::vector<row_datum> & data)
        {
          // Add any new constraint rows we've found
          const std::size_t query_size = ids.size();

          for (std::size_t i=0; i != query_size; ++i)
            {
              const dof_id_type constrained = ids[i];

              // An empty row is an constraint with an empty row; for
              // no constraint we use a "no row" placeholder
              if (data[i].empty())
                {
                  DofConstraintRow & row = _dof_constraints[constrained];
                  row.clear();
                }
              else if (data[i][0].first != DofObject::invalid_id)
                {
                  DofConstraintRow & row = _dof_constraints[constrained];
                  row.clear();
                  for (auto & pair : data[i])
                    {
                      libmesh_assert_less(pair.first, this->n_dofs());
                      row[pair.first] = pair.second;
                    }

                  // And prepare to check for more recursive constraints
                  unexpanded_dofs.insert(constrained);
                }
            }
        };

      auto rhss_action_functor =
        [this,
         max_qoi_num]
        (processor_id_type,
         const std::vector<dof_id_type> & ids,
         const std::vector<rhss_datum> & data)
        {
          // Add rhs data for any new constraint rows we've found
          const std::size_t query_size = ids.size();

          for (std::size_t i=0; i != query_size; ++i)
            {
              if (!data[i].empty())
                {
                  dof_id_type constrained = ids[i];
                  if (data[i][0] != Number(0))
                    _primal_constraint_values[constrained] = data[i][0];
                  else
                    _primal_constraint_values.erase(constrained);

                  for (unsigned int q = 0; q != max_qoi_num; ++q)
                    {
                      AdjointDofConstraintValues::iterator adjoint_map_it =
                        _adjoint_constraint_values.find(q);

                      if ((adjoint_map_it == _adjoint_constraint_values.end()) &&
                          data[i][q+1] == Number(0))
                        continue;

                      if (adjoint_map_it == _adjoint_constraint_values.end())
                        adjoint_map_it = _adjoint_constraint_values.emplace
                          (q, DofConstraintValueMap()).first;

                      DofConstraintValueMap & constraint_map =
                        adjoint_map_it->second;

                      if (data[i][q+1] != Number(0))
                        constraint_map[constrained] =
                          data[i][q+1];
                      else
                        constraint_map.erase(constrained);
                    }
                }
            }

        };

      // Now request constraint rows from other processors
      row_datum * row_ex = nullptr;
      Parallel::pull_parallel_vector_data
        (this->comm(), requested_dof_ids, row_gather_functor,
         row_action_functor, row_ex);

      // And request constraint right hand sides from other procesors
      rhss_datum * rhs_ex = nullptr;
      Parallel::pull_parallel_vector_data
        (this->comm(), requested_dof_ids, rhss_gather_functor,
         rhss_action_functor, rhs_ex);

      // We have to keep recursing while the unexpanded set is
      // nonempty on *any* processor
      unexpanded_set_nonempty = !unexpanded_dofs.empty();
      this->comm().max(unexpanded_set_nonempty);
    }
}

void DofMap::add_constraints_to_send_list()
{
  // This function must be run on all processors at once
  parallel_object_only();

  // Return immediately if there's nothing to gather
  if (this->n_processors() == 1)
    return;

  // We might get to return immediately if none of the processors
  // found any constraints
  unsigned int has_constraints = !_dof_constraints.empty();
  this->comm().max(has_constraints);
  if (!has_constraints)
    return;

  for (const auto & i : _dof_constraints)
    {
      dof_id_type constrained_dof = i.first;

      // We only need the dependencies of our own constrained dofs
      if (!this->local_index(constrained_dof))
        continue;

      const DofConstraintRow & constraint_row = i.second;
      for (const auto & j : constraint_row)
        {
          dof_id_type constraint_dependency = j.first;

          // No point in adding one of our own dofs to the send_list
          if (this->local_index(constraint_dependency))
            continue;

          _send_list.push_back(constraint_dependency);
        }
    }
}



#endif // LIBMESH_ENABLE_CONSTRAINTS


#ifdef LIBMESH_ENABLE_AMR

void DofMap::constrain_p_dofs (unsigned int var,
                               const Elem * elem,
                               unsigned int s,
                               unsigned int p)
{
  // We're constraining dofs on elem which correspond to p refinement
  // levels above p - this only makes sense if elem's p refinement
  // level is above p.
  libmesh_assert_greater (elem->p_level(), p);
  libmesh_assert_less (s, elem->n_sides());

  const unsigned int sys_num = this->sys_number();
  FEType fe_type = this->variable_type(var);

  const unsigned int n_nodes = elem->n_nodes();
  for (unsigned int n = 0; n != n_nodes; ++n)
    if (elem->is_node_on_side(n, s))
      {
        const Node & node = elem->node_ref(n);
        const unsigned int low_nc =
          FEInterface::n_dofs_at_node (fe_type, p, elem, n);
        const unsigned int high_nc =
          FEInterface::n_dofs_at_node (fe_type, elem, n);

        // since we may be running this method concurrently
        // on multiple threads we need to acquire a lock
        // before modifying the _dof_constraints object.
        Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);

        if (elem->is_vertex(n))
          {
            // Add "this is zero" constraint rows for high p vertex
            // dofs
            for (unsigned int i = low_nc; i != high_nc; ++i)
              {
                _dof_constraints[node.dof_number(sys_num,var,i)].clear();
                _primal_constraint_values.erase(node.dof_number(sys_num,var,i));
              }
          }
        else
          {
            const unsigned int total_dofs = node.n_comp(sys_num, var);
            libmesh_assert_greater_equal (total_dofs, high_nc);
            // Add "this is zero" constraint rows for high p
            // non-vertex dofs, which are numbered in reverse
            for (unsigned int j = low_nc; j != high_nc; ++j)
              {
                const unsigned int i = total_dofs - j - 1;
                _dof_constraints[node.dof_number(sys_num,var,i)].clear();
                _primal_constraint_values.erase(node.dof_number(sys_num,var,i));
              }
          }
      }
}

#endif // LIBMESH_ENABLE_AMR


#ifdef LIBMESH_ENABLE_DIRICHLET
void DofMap::add_dirichlet_boundary (const DirichletBoundary & dirichlet_boundary)
{
  _dirichlet_boundaries->push_back(new DirichletBoundary(dirichlet_boundary));
}


void DofMap::add_adjoint_dirichlet_boundary (const DirichletBoundary & dirichlet_boundary,
                                             unsigned int qoi_index)
{
  unsigned int old_size = cast_int<unsigned int>
    (_adjoint_dirichlet_boundaries.size());
  for (unsigned int i = old_size; i <= qoi_index; ++i)
    _adjoint_dirichlet_boundaries.push_back(new DirichletBoundaries());

  _adjoint_dirichlet_boundaries[qoi_index]->push_back
    (new DirichletBoundary(dirichlet_boundary));
}


bool DofMap::has_adjoint_dirichlet_boundaries(unsigned int q) const
{
  if (_adjoint_dirichlet_boundaries.size() > q)
    return true;

  return false;
}


const DirichletBoundaries *
DofMap::get_adjoint_dirichlet_boundaries(unsigned int q) const
{
  libmesh_assert_greater(_adjoint_dirichlet_boundaries.size(),q);
  return _adjoint_dirichlet_boundaries[q];
}


DirichletBoundaries *
DofMap::get_adjoint_dirichlet_boundaries(unsigned int q)
{
  unsigned int old_size = cast_int<unsigned int>
    (_adjoint_dirichlet_boundaries.size());
  for (unsigned int i = old_size; i <= q; ++i)
    _adjoint_dirichlet_boundaries.push_back(new DirichletBoundaries());

  return _adjoint_dirichlet_boundaries[q];
}


void DofMap::remove_dirichlet_boundary (const DirichletBoundary & boundary_to_remove)
{
  // Find a boundary condition matching the one to be removed
  auto lam = [&boundary_to_remove](const DirichletBoundary * bdy)
    {return bdy->b == boundary_to_remove.b && bdy->variables == boundary_to_remove.variables;};

  auto it = std::find_if(_dirichlet_boundaries->begin(), _dirichlet_boundaries->end(), lam);

  // Delete it and remove it
  libmesh_assert (it != _dirichlet_boundaries->end());
  delete *it;
  _dirichlet_boundaries->erase(it);
}


void DofMap::remove_adjoint_dirichlet_boundary (const DirichletBoundary & boundary_to_remove,
                                                unsigned int qoi_index)
{
  libmesh_assert_greater(_adjoint_dirichlet_boundaries.size(),
                         qoi_index);

  auto lam = [&boundary_to_remove](const DirichletBoundary * bdy)
    {return bdy->b == boundary_to_remove.b && bdy->variables == boundary_to_remove.variables;};

  auto it = std::find_if(_adjoint_dirichlet_boundaries[qoi_index]->begin(),
                         _adjoint_dirichlet_boundaries[qoi_index]->end(),
                         lam);

  // Delete it and remove it
  libmesh_assert (it != _adjoint_dirichlet_boundaries[qoi_index]->end());
  delete *it;
  _adjoint_dirichlet_boundaries[qoi_index]->erase(it);
}


DirichletBoundaries::~DirichletBoundaries()
{
  for (auto & item : *this)
    delete item;
}

void DofMap::check_dirichlet_bcid_consistency (const MeshBase & mesh,
                                               const DirichletBoundary & boundary) const
{
  const std::set<boundary_id_type>& mesh_side_bcids =
    mesh.get_boundary_info().get_boundary_ids();
  const std::set<boundary_id_type>& mesh_edge_bcids =
    mesh.get_boundary_info().get_edge_boundary_ids();
  const std::set<boundary_id_type>& mesh_node_bcids =
    mesh.get_boundary_info().get_node_boundary_ids();
  const std::set<boundary_id_type>& dbc_bcids = boundary.b;

  // DirichletBoundary id sets should be consistent across all ranks
  libmesh_assert(mesh.comm().verify(dbc_bcids.size()));

  for (const auto & bc_id : dbc_bcids)
    {
      // DirichletBoundary id sets should be consistent across all ranks
      libmesh_assert(mesh.comm().verify(bc_id));

      bool found_bcid = (mesh_side_bcids.find(bc_id) != mesh_side_bcids.end() ||
                         mesh_edge_bcids.find(bc_id) != mesh_edge_bcids.end() ||
                         mesh_node_bcids.find(bc_id) != mesh_node_bcids.end());

      // On a distributed mesh, boundary id sets may *not* be
      // consistent across all ranks, since not all ranks see all
      // boundaries
      mesh.comm().max(found_bcid);

      libmesh_error_msg_if(!found_bcid,
                           "Could not find Dirichlet boundary id " << bc_id << " in mesh!");
    }
}

#endif // LIBMESH_ENABLE_DIRICHLET


#ifdef LIBMESH_ENABLE_PERIODIC

void DofMap::add_periodic_boundary (const PeriodicBoundaryBase & periodic_boundary)
{
  // See if we already have a periodic boundary associated myboundary...
  PeriodicBoundaryBase * existing_boundary = _periodic_boundaries->boundary(periodic_boundary.myboundary);

  if (!existing_boundary)
    {
      // ...if not, clone the input (and its inverse) and add them to the PeriodicBoundaries object
      // Pass the pairedboundary of the original as the boundary id of the inverse clone.
      _periodic_boundaries->emplace(periodic_boundary.myboundary, periodic_boundary.clone());
      _periodic_boundaries->emplace(periodic_boundary.pairedboundary, periodic_boundary.clone(PeriodicBoundaryBase::INVERSE));
    }
  else
    {
      // ...otherwise, merge this object's variable IDs with the existing boundary object's.
      existing_boundary->merge(periodic_boundary);

      // Do the same merging process for the inverse boundary.  Note: the inverse better already exist!
      PeriodicBoundaryBase * inverse_boundary = _periodic_boundaries->boundary(periodic_boundary.pairedboundary);
      libmesh_assert(inverse_boundary);
      inverse_boundary->merge(periodic_boundary);
    }
}




void DofMap::add_periodic_boundary (const PeriodicBoundaryBase & boundary,
                                    const PeriodicBoundaryBase & inverse_boundary)
{
  libmesh_assert_equal_to (boundary.myboundary, inverse_boundary.pairedboundary);
  libmesh_assert_equal_to (boundary.pairedboundary, inverse_boundary.myboundary);

  // Store clones of the passed-in objects. These will be cleaned up
  // automatically in the _periodic_boundaries destructor.
  _periodic_boundaries->emplace(boundary.myboundary, boundary.clone());
  _periodic_boundaries->emplace(inverse_boundary.myboundary, inverse_boundary.clone());
}


#endif


} // namespace libMesh
