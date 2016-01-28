// The libMesh Finite Element Library.
// Copyright (C) 2002-2016 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



#ifndef LIBMESH_SYSTEM_SUBSET_BY_SUBDOMAIN_H
#define LIBMESH_SYSTEM_SUBSET_BY_SUBDOMAIN_H

// Local Includes
#include "libmesh/system_subset.h"
#include "libmesh/id_types.h"
#include "libmesh/parallel_object.h"

// C++ includes
#include <cstddef>
#include <set>

namespace libMesh
{

// Forward Declarations

/**
 * This class represents a subset of the dofs of a \p System, selected
 * by the \p subdomain_id and possible the variable numbers.  The dofs
 * in the subset will be sorted.
 *
 * \author Tim Kroeger
 * \date 2010
 */
class SystemSubsetBySubdomain : public SystemSubset,
                                public ParallelObject
{
public:

  /**
   * Subclass for user-specified selection of subdomain ids to be
   * included in a \p SystemSubset.
   */
  class SubdomainSelection : public ReferenceCountedObject<SubdomainSelection>
  {
  public:

    /**
     * Constructor.
     */
    SubdomainSelection ();

    /**
     * Destructor.
     */
    virtual ~SubdomainSelection ();

    /**
     * Method that decides whether a given subdomain id is included in
     * the subset or nor.
     */
    virtual bool operator() (const subdomain_id_type & subdomain_id) const = 0;

  private:
    /**
     * This isn't a copyable object, so let's make sure nobody tries.
     *
     * We won't even bother implementing this; we'll just make sure
     * that the compiler doesn't implement a default.
     */
    SubdomainSelection (const SubdomainSelection &);

    /**
     * This isn't a copyable object, so let's make sure nobody tries.
     *
     * We won't even bother implementing this; we'll just make sure
     * that the compiler doesn't implement a default.
     */
    SubdomainSelection & operator = (const SubdomainSelection &);
  }; // subclass \p SubdomainSelection



  /**
   * Selection of subdomain ids by a list.
   */
  class SubdomainSelectionByList : public SubdomainSelection
  {
  public:
    /**
     * Constructor.  Does not take a copy of the \p list, so make sure
     * that the \p list does not go out of scope before \p *this does.
     */
    explicit
    SubdomainSelectionByList (const std::set<subdomain_id_type> & list);

    /**
     * Method that decides whether a given subdomain id is included in
     * the subset or nor.
     */
    virtual bool operator() (const subdomain_id_type & subdomain_id) const libmesh_override;

  protected:
    /**
     * The actual list.
     */
    const std::set<subdomain_id_type> & _list;
  };

  /**
   * Constructor.  The subset will consist of those dofs which are
   * associated to at least one mesh element that has a subdomain id
   * contained in the \p subdomain_selection.  If \p var_nums is not a
   * \p NULL pointer, dofs that are associated to a variable number
   * that is not contained in \p var_nums will not contain to the
   * subset, no matter what elements they belong to.
   */
  SystemSubsetBySubdomain (const System & system,
                           const SubdomainSelection & subdomain_selection,
                           const std::set<unsigned int> * const var_nums = libmesh_nullptr);

  /**
   * Constructor.  The subset will consist of those dofs which are
   * associated to at least one mesh element that has a subdomain id
   * contained in the set \p subdomain_ids.  If \p var_nums is not a
   * \p NULL pointer, dofs that are associated to a variable number
   * that is not contained in \p var_nums will not contain to the
   * subset, no matter what elements they belong to.
   */
  SystemSubsetBySubdomain (const System & system,
                           const std::set<subdomain_id_type> & subdomain_ids,
                           const std::set<unsigned int> * const var_nums = libmesh_nullptr);

  /**
   * Destructor.
   */
  virtual ~SystemSubsetBySubdomain ();

  /**
   * Method that returns the actual set of dofs that the subset
   * consists of.  The result will contain local dofs on each
   * processor only and will not contain duplictates.
   */
  virtual const std::vector<unsigned int> & dof_ids () const libmesh_override;

  /**
   * Initializes the class.  Will be called by the constructors.  Can
   * also be called manually to update the subset.  This is required
   * if (a) the subdomain ids of some elements have changed in the
   * meantime and you want these changes to take effect, or (b) you
   * want to use a different \p SubdomainSelection object now.
   */
  void init (const SubdomainSelection & subdomain_selection);

  /**
   * Initializes the class.  Will be called by the constructors.  Can
   * also be called manually to update the subset.  This is required
   * if (a) the subdomain ids of some elements have changed in the
   * meantime and you want these changes to take effect, or (b) you
   * want to use a different list of subdomain ids now.
   */
  void init (const std::set<subdomain_id_type> & subdomain_ids);

protected:

  /**
   * Sets \p _var_nums to either a copy of \p var_nums or, if that is
   * \p NULL, a set of all variable numbers that occur in the system.
   */
  void set_var_nums (const std::set<unsigned int> * const var_nums);

  /**
   * The set of all variable numbers that are contained in the subset.
   * This will be set by the constructor to either a copy of its \p
   * var_nums argument or, if that is \p NULL, a set of all variable
   * numbers that occur in the system.
   */
  std::set<unsigned int> _var_nums;

  /**
   * The actual set of the dof ids.
   */
  std::vector<unsigned int> _dof_ids;

}; // class SystemSubset

} // namespace libMesh

#endif // LIBMESH_SYSTEM_SUBSET_BY_SUBDOMAIN_H
