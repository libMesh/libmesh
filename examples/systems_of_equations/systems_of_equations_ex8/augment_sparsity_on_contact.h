#ifndef AUGMENT_SPARSITY_ON_CONTACT_H
#define AUGMENT_SPARSITY_ON_CONTACT_H

#include "libmesh/dof_map.h"
#include "libmesh/system.h"
#include LIBMESH_INCLUDE_UNORDERED_MAP

using libMesh::DofMap;
using libMesh::System;
using libMesh::dof_id_type;
using libMesh::boundary_id_type;

class AugmentSparsityOnContact : public DofMap::AugmentSparsityPattern
{
private:

  /**
   * The System object that we're using here.
   */
  System& _sys;

  /**
   * This maps from element --> elements it is in contact with.
   * We use a set as the value to enfore uniqueness of the values.
   */
  LIBMESH_BEST_UNORDERED_MAP<dof_id_type, std::set<dof_id_type> > _contact_element_map;

public:

  /**
   * Constructor.
   */
  AugmentSparsityOnContact(
    System& _sys);

  /**
   * Clear the contact element map.
   */
  void clear_contact_element_map();

  /**
   * Add a new entry to the multimap that defines
   * which elements element \p element_id is in
   * contact with.
   */
  void add_contact_element(
    dof_id_type element_id,
    dof_id_type other_element_id);

  /**
   * User-defined function to augment the sparsity pattern.
   */
  virtual void augment_sparsity_pattern (
    libMesh::SparsityPattern::Graph & ,
    std::vector<dof_id_type> & n_nz,
    std::vector<dof_id_type> & n_oz);

};

#endif
