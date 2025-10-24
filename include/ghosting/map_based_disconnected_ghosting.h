// libMesh Includes
#include "libmesh/ghosting_functor.h" // base class

// C++ includes
#include <memory> // for std::unique_ptr
#include <map>    // for std::map

namespace libMesh
{

// Forward declarations
class MeshBase;

/**
 * @brief A GhostingFunctor subclass that uses an externally-provided map
 * to define non-manifold neighbor relationships and request mesh-level ghosting.
 *
 * This functor is designed to handle manually defined disconnected neighbor
 * relationships (e.g., those created by add_disconnected_boundaries()),
 * ensuring that remote neighbor elements are properly ghosted in parallel.
 */
class MapBasedDisconnectedGhosting : public GhostingFunctor
{
public:

  // Map from (local element ID, local side number) to neighbor element ID.
  using DisconnectedMap =
      std::map<dof_id_type, dof_id_type>;

  /**
   * @brief Constructor.
   * Requires the mesh and a reference to the neighbor relationship map.
   * Note: The lifetime of the map must exceed the lifetime of this functor.
   */
  MapBasedDisconnectedGhosting(MeshBase & mesh,
                               const DisconnectedMap & map);

  /**
   * @brief Virtual destructor.
   */
  virtual ~MapBasedDisconnectedGhosting() = default;

  /**
   * @brief clone() calls the copy constructor. Required by libMesh.
   */
  virtual std::unique_ptr<GhostingFunctor> clone () const override;

  /**
   * @brief Ghosting operator.
   */
  virtual void operator() (const MeshBase::const_element_iterator & range_begin,
                           const MeshBase::const_element_iterator & range_end,
                           processor_id_type p,
                           map_type & coupled_elements) override;

private:
  /**
   * @brief A reference to the external map defining neighbor relationships.
   */
  const DisconnectedMap & _map;
};

} // namespace libMesh
