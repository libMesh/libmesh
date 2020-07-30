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



#ifndef LIBMESH_GHOSTING_FUNCTOR_H
#define LIBMESH_GHOSTING_FUNCTOR_H

// Local Includes
#include "libmesh/libmesh_common.h"
#include "libmesh/id_types.h"
#include "libmesh/mesh_base.h"
#include "libmesh/reference_counted_object.h"
#include "libmesh/auto_ptr.h" // libmesh_make_unique

// C++ Includes
#include <unordered_map>

namespace libMesh
{

// Forward Declarations
class CouplingMatrix;
class Elem;


/**
 * This abstract base class defines the interface by which library
 * code and user code can report associations between elements.  These
 * associations between elements may be used to keep copies of
 * non-local elements retained on a distributed mesh, may be used to
 * keep element degrees of freedom properly communicated on a
 * distributed vector, and/or may be used to insert coupling entries
 * between elements' degrees of freedom to a sparsity matrix.
 *
 *
 * We can think of three levels of "element dependencies". An element
 * K1 has a coupling dependency on K2 if the dofs on K1 might (modulo
 * the coupling matrix) need sparsity pattern entries for dofs on K2.
 * An element K1 has an algebraic dependency on K2 if a processor
 * which owns K1 might need to examine the solution dof values on K2.
 * An element K1 has a geometric dependency on K2 if a processor which
 * owns K1 might need to examine the geometry of K2. For any element
 * K, we could call the set of coupling-ghosted ("coupled") elements
 * C(K), call the set of algebraic-ghosted ("evaluable") elements E(K),
 * and call the set of geometry-ghosted ("ghosted") elements G(K).
 *
 * It should be safe to assume that, for any element K, C(K) implies
 * E(K) implies G(K). These will be one-way implications in some
 * problems and equality relations in others.
 *
 * We can think of these as operators on sets of elements in the
 * obvious way, e.g.: G({K}) = {G(Ki) for all Ki in {K}}.
 *
 * The user should omit from their implementations relations which we
 * already have enough information to understand implicitly. For
 * instance, K is obviously in C(K), so a GhostingFunctor should never
 * bother telling us so. We may have a PeriodicBoundary, a hanging
 * node constraint equation, or a user-defined constraint equation
 * which creates a dependency between two elements; if so then we
 * don't need the user to also tell us about that relation (although
 * a PeriodicBoundary use will cause the library to create its own
 * PeriodicGhostingFunctor for internal use).
 *
 * Users may only care about a subset of variables in distant
 * evaluable elements, so we could imagine defining E_v(K) for each
 * variable number v, in which case E(K) under our previous definition
 * is the union of E_v(K) forall v, and this is what would be included
 * in G(K) when deciding what to ghost, but our send_list would only
 * include the subset of variables we need, so communication and
 * memory use would be much reduced.  However, for efficiency and API
 * simplicity, we instead define the isomorphic operator E'(K) which
 * gives a set of ordered pairs of elements and variable-number-sets,
 * from which a consistent E(K) would be derived by ignoring the
 * variable-number-sets.
 *
 * For C(K), there are similar issues: e.g. we may want some
 * discontinuous variables to couple only within their own element but
 * other discontinuous variables to couple in a DG/FVM way. This could
 * extend to more than one variable index: i.e. a dof for variable v1
 * in element K1 would depend on a dof for variable v2 in element K2
 * iff K2 is in C_v1_v2(K1). This would induce a consistent E_v(K) =
 * union of C_w_v(K) forall variable indices w. Again, the equivalent
 * API alternative we use here is for C'(K) to return a set of ordered
 * pairs of elements and variable-number-pair-sets.
 *
 * We return variable-stuff-sets as pointers to CouplingMatrix.  That
 * way, in the common case where the user cares about all variables and
 * couple to all variables, all the functor needs to return for
 * variable-number-sets and variable-number-pair-sets is nullptr (which
 * as in other libMesh APIs will be interpreted and documented to mean
 * the full set). In the common case where the user wants coupling
 * between elements to match coupling within elements, the functor can
 * return the same pointer as DofMap::_coupling_matrix. Even in the
 * less common cases, the user can store common variable-number-sets
 * and variable-number-pair sets as CouplingMatrix members of a
 * functor object of their subclass, and setting up a few of those
 * matrices then setting lots of those pointers is cheap.
 *
 * The definition of the CouplingMatrix for a variable-dependent E'
 * should be consistent with the requirements that would have been
 * imposed had that matrix been used for a C'.  In other words, if the
 * returned CouplingMatrix CM has CM(i,j)==true for any i, then
 * variable j will be evaluable on the algebraically ghosted element.
 *
 * After a GhostingFunctor has been created, a reference to it can be
 * passed to MeshBase::add_ghosting_functor to expand geometric
 * ghosting, or to DofMap::add_algebraic_ghosting_functor to expand
 * both algebraic and geometric ghosting, or to
 * DofMap::add_coupling_functor to expand coupling along with both
 * types of ghosting.
 *
 * Note that when an element is specified in algebraic ghosting or
 * coupling queries, only degrees of freedom for the variables
 * supported on that element are thereby ghosted and/or coupled.
 * Any unsupported variable dofs associated with the element's nodes
 * (e.g. subdomain-restricted variables on a neighboring subdomain)
 * will be unaffected.
 *
 * Typical usage of the GhostingFunctor would be to add a geometric ghosting
 * functor before the mesh preparation is completed; progmatically, this would
 * be before MeshBase::prepare_for_use() is called, but many different libMesh
 * idioms internally call this function. The algebraic and coupling ghosting
 * functors normally are added before EquationSystems::init() is called.
 * However, in some circumstances, solution evaluation may be needed within the
 * GhostingFunctor in order to determine the ghosting, in which case the appropriate
 * functor would need to be added after EquationSystems::init(). In this case,
 * the user will need to reinitialize certain parts of the DofMap for
 * algebraic and coupling functors. For algebraic ghosting functors, the
 * user will need to call DofMap::reinit_send_list() and then reinitialize
 * any NumericVectors that are GHOSTED, e.g. the System::current_local_solution.
 * For coupling ghosting, the user will also need to recompute the sparsity
 * pattern via DofMap::clear_sparsity() and then DofMap::compute_sparsity() and
 * then reinitialize any SparseMatrix objects attached to the System, e.g.
 * the system.get_matrix("System Matrix").
 *
 * \author Roy H. Stogner
 * \date 2016
 */
class GhostingFunctor : public ReferenceCountedObject<GhostingFunctor>
{
public:

  /**
   * Constructor.  Empty in the base class.
   */
  GhostingFunctor() {}

  /**
   * Copy Constructor
   */
  GhostingFunctor(const GhostingFunctor & other) : ReferenceCountedObject<GhostingFunctor>(other) {}

  /**
   * Virtual destructor; this is an abstract base class.
   */
  virtual ~GhostingFunctor() {}

  /**
   * A clone() is needed because GhostingFunctor can not be shared between
   * different meshes. The operations in  GhostingFunctor are mesh dependent.
   */
  virtual std::unique_ptr<GhostingFunctor> clone () const
  { libmesh_not_implemented(); }

  /**
   * It should be called after cloning a ghosting functor.
   * Ghosting functor is mesh dependent
   */
  virtual void set_mesh(const MeshBase * /*mesh*/) { libmesh_not_implemented(); }

  /**
   * What elements do we care about and what variables do we care
   * about on each element?
   */
  typedef std::unordered_map<const Elem*, const CouplingMatrix*> map_type;

  /**
   * For the specified range of active elements, what other elements
   * currently living (whether local or ghosted) on this processor
   * need to be coupled/ghosted to accommodate them?  Don't bother to
   * return any results which already have processor_id p.
   *
   * This API is new, and we should replace "ignoring those on
   * processor p" with "ignoring those which match a predicate
   * functor" eventually.
   */
  virtual void operator() (const MeshBase::const_element_iterator & range_begin,
                           const MeshBase::const_element_iterator & range_end,
                           processor_id_type p,
                           map_type & coupled_elements) = 0;

  /**
   * GhostingFunctor subclasses which cache data will need to
   * initialize that cache.  We call mesh_reinit() whenever the
   * relevant Mesh has changed, but before remote elements on a
   * distributed mesh are deleted.
   */
  virtual void mesh_reinit () {};

  /**
   * For algebraic ghosting or coupling functors we also call
   * dofmap_reinit() later, after dofs have been distributed on the
   * new mesh but before the functors have been queried for send_list
   * or sparsity pattern calculations.
   */
  virtual void dofmap_reinit () {};

  /**
   * GhostingFunctor subclasses with relatively long-lasting caches
   * may want to redistribute those caches whenever the relevant Mesh
   * is redistributed; we will give them an opportunity when that
   * happens.  At the point in the code where this is called, element
   * processor ids have been set to their new destinations, and those
   * elements have been copied to their new destinations, but the
   * elements have not yet been deleted by the processors which
   * previously held them..
   */
  virtual void redistribute () {};

  /**
   * GhostingFunctor subclasses with relatively long-lasting caches
   * may want to delete the no-longer-relevant parts of those caches
   * after a redistribution is complete.
   */
  virtual void delete_remote_elements () {};
};

} // namespace libMesh

#endif // LIBMESH_GHOSTING_FUNCTOR_H
