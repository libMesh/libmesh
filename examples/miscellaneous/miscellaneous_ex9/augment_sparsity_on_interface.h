#ifndef AUGMENT_SPARSITY_ON_INTERFACE_H
#define AUGMENT_SPARSITY_ON_INTERFACE_H

#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"

using libMesh::DofMap;
using libMesh::EquationSystems;
using libMesh::dof_id_type;
using libMesh::boundary_id_type;

// Convenient typedef for a map for (element id, side id) --> element neighbor id
typedef std::map<std::pair<dof_id_type, unsigned char>, dof_id_type> ElementIdMap;

class AugmentSparsityOnInterface : public DofMap::AugmentSparsityPattern
{
private:

  /**
   * The EquationSystems object that we're using here.
   */
  EquationSystems & _es;

  /**
   * A map from (lower element ID, side ID) to matching upper element ID. Here "lower"
   * and "upper" refer to the lower and upper (wrt +z direction) sides of the crack in
   * our mesh.
   */
  ElementIdMap _lower_to_upper;

  /**
   * Boundary IDs for the lower and upper faces of the "crack" in the mesh.
   */
  boundary_id_type _crack_boundary_lower, _crack_boundary_upper;

public:

  /**
   * Constructor.
   */
  AugmentSparsityOnInterface(EquationSystems & es,
                             boundary_id_type crack_boundary_lower,
                             boundary_id_type crack_boundary_upper);

  /**
   * @return a const reference to the lower-to-upper element ID map.
   */
  const ElementIdMap & get_lower_to_upper() const;

  /**
   * User-defined function to augment the sparsity pattern.
   */
  virtual void augment_sparsity_pattern (libMesh::SparsityPattern::Graph & ,
                                         std::vector<dof_id_type> & n_nz,
                                         std::vector<dof_id_type> & n_oz);
};

#endif
