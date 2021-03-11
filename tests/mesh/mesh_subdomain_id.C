#include <libmesh/libmesh.h>
#include <libmesh/elem.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class MeshSubdomainIDTest : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify proper operation of the
   * MeshBase::subdomain_ids() method.
   */
public:
  CPPUNIT_TEST_SUITE( MeshSubdomainIDTest );

  CPPUNIT_TEST( testUnpartitioned );
  CPPUNIT_TEST( testMultiple );

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp() {}

  void tearDown() {}

  void testUnpartitioned()
  {
    std::unique_ptr<UnstructuredMesh> mesh = libmesh_make_unique<Mesh>(*TestCommWorld);

    mesh->add_point(Point(0, 0, 0), 0);
    mesh->add_point(Point(1, 0, 0), 1);

    Elem * elem = mesh.add_elem(Elem::build(EDGE2));
    elem->set_node(0) = mesh.node_ptr(0);
    elem->set_node(1) = mesh.node_ptr(1);

    std::set<subdomain_id_type> ids;
    mesh->subdomain_ids(ids);

    CPPUNIT_ASSERT_EQUAL(mesh->n_subdomains(), (subdomain_id_type)1);
    CPPUNIT_ASSERT_EQUAL(ids.size(), (std::size_t)1);
    CPPUNIT_ASSERT_EQUAL(*ids.begin(), (subdomain_id_type)0);
  }

  void testMultiple()
  {
    std::unique_ptr<UnstructuredMesh> mesh = libmesh_make_unique<Mesh>(*TestCommWorld);

    MeshTools::Generation::build_line(*mesh, 5, 0.0, 1.0, EDGE2);

    std::set<subdomain_id_type> actual_ids;
    for (auto & elem : mesh->active_element_ptr_range())
      {
        CPPUNIT_ASSERT(elem->id() < (dof_id_type)Elem::invalid_subdomain_id);
        actual_ids.insert((subdomain_id_type)elem->id());
        elem->subdomain_id() = (subdomain_id_type)elem->id();
      }

    std::set<subdomain_id_type> ids;
    mesh->subdomain_ids(ids);

    CPPUNIT_ASSERT_EQUAL(ids.size(), actual_ids.size());
    for (auto id : ids)
      CPPUNIT_ASSERT(actual_ids.count(id));
  }
};


CPPUNIT_TEST_SUITE_REGISTRATION( MeshSubdomainIDTest );
