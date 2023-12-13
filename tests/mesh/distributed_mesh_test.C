// libmesh includes
#include <libmesh/distributed_mesh.h>
#include <libmesh/elem.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/utility.h>

// cppunit includes
#include "test_comm.h"
#include "libmesh_cppunit.h"

// C++ includes
#include <regex>

using namespace libMesh;


class DistributedMeshTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( DistributedMeshTest );

  CPPUNIT_TEST( testRemoteElemError );

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

  void testRemoteElemError()
  {
    LOG_UNIT_TEST;

    DistributedMesh mesh(*TestCommWorld);

    // Build enough elements to be sure some will be remote on 2
    // processors
    MeshTools::Generation::build_line(mesh, 10, EDGE2);

#ifdef LIBMESH_ENABLE_EXCEPTIONS
    if (mesh.n_processors() > 1)
      {
        // We only expect a throw if we can get to a remote_elem to
        // throw from, which can't happen if we're on 11+ procs
        bool threw_expected_error = true;
        try
        {
          for (const auto & elem : mesh.element_ptr_range())
            {
              threw_expected_error = false;
              for (const Elem * neigh : elem->neighbor_ptr_range())
                if (neigh && neigh->n_sides() != 2) // this should throw if not 2
                  CPPUNIT_ASSERT(false);
            }
        }
        catch (libMesh::LogicError & e)
        {
          std::regex msg_regex("merely a shim");
          CPPUNIT_ASSERT(std::regex_search(e.what(), msg_regex));
          threw_expected_error = true;
        }

        CPPUNIT_ASSERT(threw_expected_error);
      }
#endif // LIBMESH_ENABLE_EXCEPTIONS
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( DistributedMeshTest );
