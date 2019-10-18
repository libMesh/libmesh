// Basic include files
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/parallel.h" // set_union
#include "libmesh/string_to_enum.h"
#include "libmesh/boundary_info.h"

#include "test_comm.h"
#include "libmesh_cppunit.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

class WriteEdgesetData : public CppUnit::TestCase
{
  /**
   * This test ensures you can write both vector and scalar variables
   */
public:
  CPPUNIT_TEST_SUITE(WriteEdgesetData);

#if LIBMESH_DIM > 1
  CPPUNIT_TEST(testWrite);
#endif

  CPPUNIT_TEST_SUITE_END();

  void testWrite()
  {
    Mesh mesh(*TestCommWorld);

    MeshTools::Generation::build_cube (mesh,
                                       5, 5, 5,
                                       0., 1.,
                                       0., 1.,
                                       0., 1.,
                                       HEX8);

    // Add named edgesets containing the "front" and "back" edges of
    // the cube.
    BoundaryInfo & bi = mesh.get_boundary_info();

    std::vector<BoundaryInfo::BCTuple> bc_tuples =
      bi.build_side_list();

    for (const auto & t : bc_tuples)
      {
        dof_id_type elem_id = std::get<0>(t);
        boundary_id_type b_id = std::get<2>(t);

        // Side 0 corresponds to the "back" sideset, with edges {0,1,2,3}
        // Side 5 corresponds to the "front" sideset, with edges {4,5,6,7}
        // If we are on either of those sides, add all the corresponding edges.
        if (b_id == 0)
          for (unsigned short edge=0; edge<4; ++edge)
            bi.add_edge(elem_id, edge, b_id);

        if (b_id == 5)
          for (unsigned short edge=8; edge<12; ++edge)
            bi.add_edge(elem_id, edge, b_id);
      }

    // Name the edgesets
    bi.edgeset_name(0) = "back_edgeset";
    bi.edgeset_name(5) = "front_edgeset";

#ifdef LIBMESH_HAVE_EXODUS_API

    // We write the file in the ExodusII format.
    {
      ExodusII_IO writer(mesh);
      writer.write("write_edgeset_data.e");
    }

    // Make sure that the writing is done before the reading starts.
    TestCommWorld->barrier();

    // Now read it back in
    Mesh read_mesh(*TestCommWorld);
    ExodusII_IO reader(read_mesh);
    reader.read("write_edgeset_data.e");

    // Assert that we got back out what we put in.
    BoundaryInfo & read_bi = read_mesh.get_boundary_info();
    CPPUNIT_ASSERT(read_bi.edgeset_name(0) == "back_edgeset");
    CPPUNIT_ASSERT(read_bi.edgeset_name(5) == "front_edgeset");
    CPPUNIT_ASSERT(read_bi.n_edge_conds() == 200);
    CPPUNIT_ASSERT(read_bi.get_edge_boundary_ids().size() == 2);

    // Make sure that we have the expected amount of edges in each edgeset.
    std::map<boundary_id_type, unsigned int> counts;
    std::vector<BoundaryInfo::BCTuple> edge_tuples = bi.build_edge_list();
    for (const auto & t : edge_tuples)
      counts[std::get<2>(t)]++;
    CPPUNIT_ASSERT(counts[0] == 100);
    CPPUNIT_ASSERT(counts[5] == 100);

#endif // #ifdef LIBMESH_HAVE_EXODUS_API
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(WriteEdgesetData);
