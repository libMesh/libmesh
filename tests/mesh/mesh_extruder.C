#include <libmesh/libmesh.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/mesh_generation.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class MeshExtruderTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify proper operation of the Mesh Extruder
   * with the optional object callback for setting custom subdomain IDs.
   * We pass a custom object for generating subdomains based on the old element
   * ID and the current layer and assert the proper values.
   */
public:
  LIBMESH_CPPUNIT_TEST_SUITE( MeshExtruderTest );

#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testExtruder );
#endif

  CPPUNIT_TEST_SUITE_END();

protected:

  class QueryElemSubdomainID : public MeshTools::Generation::QueryElemSubdomainIDBase
  {
    /// The override from the base class for obtaining a new id based on the old (original) element and the specified layer
    virtual subdomain_id_type get_subdomain_for_layer(const Elem * old_elem, unsigned int layer)
    {
      // This method will assign an new id based on the old element and the layer
      return old_elem->subdomain_id() + layer;
    }
  };

public:
  void setUp() {}

  void tearDown() {}

  void testExtruder()
  {
    LOG_UNIT_TEST;

    ReplicatedMesh src_mesh(*TestCommWorld, /*dim=*/2);

    const unsigned int n_elems_per_side = 4;
    const unsigned int num_layers = 4;
    const unsigned int n_elems_per_layer = n_elems_per_side * n_elems_per_side;


    MeshTools::Generation::build_square(src_mesh, n_elems_per_side, n_elems_per_side);
    for (unsigned int i=0; i<n_elems_per_layer; ++i)
      {
        // Retrieve the element from the mesh by ID to guarantee proper ordering instead of with iterators
        Elem & elem = src_mesh.elem_ref(i);
        elem.subdomain_id() = i;
      }

    ReplicatedMesh dest_mesh(*TestCommWorld, /*dim=*/3);

    RealVectorValue extrusion_vector(0, 0, 1);

    QueryElemSubdomainID new_elem_subdomain_id;

    /**
     * The test mesh is designed to be square with the subdomain corresponding to the element number.
     * We will use this to assert the correct pattern from the custom extruder.
     */
    MeshTools::Generation::build_extrusion(dest_mesh, src_mesh, num_layers, extrusion_vector, &new_elem_subdomain_id);

    for (unsigned int i=0; i<n_elems_per_layer * num_layers; ++i)
      {
        // Retrieve the element from the mesh by ID to guarantee proper ordering instead of with iterators
        Elem & elem = dest_mesh.elem_ref(i);

        CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(elem.subdomain_id()), i%n_elems_per_layer + i/n_elems_per_layer /* integer division */);
      }
  }
};


CPPUNIT_TEST_SUITE_REGISTRATION( MeshExtruderTest );
