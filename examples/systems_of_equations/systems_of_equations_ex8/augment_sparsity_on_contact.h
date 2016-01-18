#ifndef AUGMENT_SPARSITY_ON_CONTACT_H
#define AUGMENT_SPARSITY_ON_CONTACT_H

#include "libmesh/dof_map.h"
#include "libmesh/system.h"
#include LIBMESH_INCLUDE_UNORDERED_MAP

using libMesh::DofMap;
using libMesh::System;
using libMesh::Node;
using libMesh::dof_id_type;
using libMesh::boundary_id_type;

class AugmentSparsityOnContact : public DofMap::AugmentSparsityPattern
{
public:

  /**
   * Constructor.
   */
  AugmentSparsityOnContact(System & _sys);

  /**
   * Clear the contact element map.
   */
  void clear_contact_node_map();

  /**
   * Add a new entry to _contact_node_map.
   */
  void add_contact_node(dof_id_type lower_node_id,
                        dof_id_type upper_node_id);

  /**
   * User-defined function to augment the sparsity pattern.
   */
  virtual void augment_sparsity_pattern (libMesh::SparsityPattern::Graph &,
                                         std::vector<dof_id_type> & n_nz,
                                         std::vector<dof_id_type> & n_oz);

  /**
   * Helper function that actually sets the sparsity values.
   */
  void set_sparsity_values(const Node & this_node,
                           const Node & other_node,
                           std::vector<dof_id_type> & n_nz,
                           std::vector<dof_id_type> & n_oz);

  /**
   * The System object that we're using here.
   */
  System & _sys;

  /**
   * This provides a map between contact nodes.
   */
  LIBMESH_BEST_UNORDERED_MAP<dof_id_type, dof_id_type> _contact_node_map;
};

#endif
