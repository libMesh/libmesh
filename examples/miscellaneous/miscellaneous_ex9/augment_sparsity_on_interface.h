#ifndef AUGMENT_SPARSITY_ON_INTERFACE_H
#define AUGMENT_SPARSITY_ON_INTERFACE_H

#include "libmesh/ghosting_functor.h"
#include "libmesh/mesh_base.h"

using libMesh::Elem;
using libMesh::GhostingFunctor;
using libMesh::MeshBase;
using libMesh::boundary_id_type;
using libMesh::processor_id_type;

// Convenient typedef for a map for (element, side id) --> element neighbor
typedef std::map<std::pair<const Elem *, unsigned char>, const Elem *> ElementSideMap;

// And the inverse map, but without sides in the key
typedef std::map<const Elem *, const Elem *> ElementMap;

class AugmentSparsityOnInterface : public GhostingFunctor
{
private:

  /**
   * The Mesh we're calculating on
   */
  MeshBase & _mesh;

  /**
   * A map from (lower element ID, side ID) to matching upper element ID. Here "lower"
   * and "upper" refer to the lower and upper (wrt +z direction) sides of the crack in
   * our mesh.
   */
  ElementSideMap _lower_to_upper;

  /**
   * The inverse (ignoring sides) of the above map.
   */
  ElementMap _upper_to_lower;

  /**
   * Boundary IDs for the lower and upper faces of the "crack" in the mesh.
   */
  boundary_id_type _crack_boundary_lower, _crack_boundary_upper;

  /**
   * Make sure we've been initialized before use
   */
  bool _initialized;

public:

  /**
   * Constructor.
   */
  AugmentSparsityOnInterface(MeshBase & mesh,
                             boundary_id_type crack_boundary_lower,
                             boundary_id_type crack_boundary_upper);


  virtual std::unique_ptr<GhostingFunctor> clone () const override
  {
    return std::make_unique<AugmentSparsityOnInterface>
      (_mesh, _crack_boundary_lower, _crack_boundary_upper);
  }

  /**
   * @return a const reference to the lower-to-upper element ID map.
   */
  const ElementSideMap & get_lower_to_upper() const;

  /**
   * User-defined function to augment the sparsity pattern.
   */
  virtual void operator() (const MeshBase::const_element_iterator & range_begin,
                           const MeshBase::const_element_iterator & range_end,
                           processor_id_type p,
                           map_type & coupled_elements) override;

  /**
   * Rebuild the cached _lower_to_upper map whenever our Mesh has
   * changed.
   */
  virtual void mesh_reinit () override;

  /**
   * Update the cached _lower_to_upper map whenever our Mesh has been
   * redistributed.  We'll be lazy and just recalculate from scratch.
   */
  virtual void redistribute () override
  { this->mesh_reinit(); }


};

#endif
