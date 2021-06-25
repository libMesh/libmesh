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

class WriteNodesetData : public CppUnit::TestCase
{
public:
  CPPUNIT_TEST_SUITE(WriteNodesetData);

#if LIBMESH_DIM > 1
  CPPUNIT_TEST(testWrite);
#endif

  CPPUNIT_TEST_SUITE_END();

  template <typename IOClass>
  void testWriteImpl(const std::string & filename)
  {
    Mesh mesh(*TestCommWorld);
    mesh.allow_renumbering(false);
    MeshTools::Generation::build_square(mesh,
                                        /*nx=*/5, /*ny=*/5,
                                        -1., 1.,
                                        -1., 1.,
                                        QUAD4);

    // Add an empty nodeset
    mesh.get_boundary_info().nodeset_name(4) = "empty";

    BoundaryInfo & bi = mesh.get_boundary_info();

    // Meshes created via build_square() don't have any nodesets
    // defined on them by default, so let's create some now. Note that
    // this function handles the synchronization of ghost element and
    // node ids in parallel internally.
    bi.build_node_list_from_side_list();

    // Get list of all (node, id) tuples
    std::vector<BoundaryInfo::NodeBCTuple> all_bc_tuples = bi.build_node_list();

    // Data structures to be passed to ExodusII_IO::write_nodeset_data()
    std::vector<std::string> var_names = {"var1", "var2", "var3"};
    std::vector<std::set<boundary_id_type>> node_boundary_ids =
      {
        {0, 2}, // var1 is defined on nodesets 0 and 2
        {1, 3}, // var2 is defined on nodesets 1 and 3
        {4}     // var3 is only defined on the empty nodeset 4
      };

    // Data structure mapping (node, id) tuples to Real values that
    // will be passed to Exodus.
    std::vector<std::map<BoundaryInfo::NodeBCTuple, Real>> bc_vals(var_names.size());

    // For each var_names[i], construct bc_vals[i]
    for (unsigned int i=0; i<var_names.size(); ++i)
      {
        // const auto & var_name = var_names[i];
        auto & vals = bc_vals[i];

        for (const auto & t : all_bc_tuples)
          {
            boundary_id_type b_id = std::get<1>(t);

            if (node_boundary_ids[i].count(b_id))
              {
                // Compute a value. This could in theory depend on
                // var_name, node_id, and/or b_id.
                Real val = static_cast<Real>(b_id);

                // Insert into the vals map.
                vals.emplace(t, val);
              }
          }

        // If we have a distributed mesh, the ExodusII writer
        // serializes it before writing, so we need to serialize our
        // nodeset data as well so that it can all be written from
        // processor 0.
        if (!mesh.is_serial())
          TestCommWorld->set_union(vals);

      } // done constructing bc_vals

    // We write the file in the ExodusII format.
    {
      IOClass writer(mesh);
      writer.write(filename);
      writer.write_nodeset_data (/*timestep=*/1, var_names, node_boundary_ids, bc_vals);
    }

    // Make sure that the writing is done before the reading starts.
    TestCommWorld->barrier();

    // Now read it back in
    Mesh read_mesh(*TestCommWorld);
    IOClass reader(read_mesh);
    reader.read(filename);

    std::vector<std::string> read_in_var_names;
    std::vector<std::set<boundary_id_type>> read_in_node_boundary_ids;
    std::vector<std::map<BoundaryInfo::NodeBCTuple, Real>> read_in_bc_vals;
    reader.read_nodeset_data
      (/*timestep=*/1, read_in_var_names, read_in_node_boundary_ids, read_in_bc_vals);

    // Assert that we got back out what we put in.
    CPPUNIT_ASSERT(read_in_var_names == var_names);
    CPPUNIT_ASSERT(read_in_node_boundary_ids == node_boundary_ids);
    CPPUNIT_ASSERT(read_in_bc_vals == bc_vals);
  }

  void testWrite()
  {
#ifdef LIBMESH_HAVE_EXODUS_API
    testWriteImpl<ExodusII_IO>("write_nodeset_data.e");
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(WriteNodesetData);
