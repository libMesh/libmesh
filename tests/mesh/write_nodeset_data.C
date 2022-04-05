// Basic include files
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/nemesis_io.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/parallel.h" // set_union
#include "libmesh/string_to_enum.h"
#include "libmesh/boundary_info.h"
#include "libmesh/utility.h" // libmesh_map_find

#include "test_comm.h"
#include "libmesh_cppunit.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

class WriteNodesetData : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(WriteNodesetData);

#if LIBMESH_DIM > 1
#ifdef LIBMESH_HAVE_EXODUS_API
  CPPUNIT_TEST(testWriteExodus);
#endif // #ifdef LIBMESH_HAVE_EXODUS_API
#ifdef LIBMESH_HAVE_NEMESIS_API
  // CPPUNIT_TEST(testWriteNemesis); // Not yet implemented
#endif
#endif

  CPPUNIT_TEST_SUITE_END();

  template <typename IOClass>
  void testWriteImpl(const std::string & filename,
                     bool write_vars)
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

    // Only used if write_vars == true
    std::vector<std::string> var_names;
    std::vector<std::set<boundary_id_type>> node_boundary_ids;
    std::vector<std::map<BoundaryInfo::NodeBCTuple, Real>> bc_vals;

    if (write_vars)
      {
        BoundaryInfo & bi = mesh.get_boundary_info();

        // Meshes created via build_square() don't have any nodesets
        // defined on them by default, so let's create some now. Note that
        // this function handles the synchronization of ghost element and
        // node ids in parallel internally.
        bi.build_node_list_from_side_list();

        // Get list of all (node, id) tuples
        std::vector<BoundaryInfo::NodeBCTuple> all_bc_tuples = bi.build_node_list();

        // Data structures to be passed to ExodusII_IO::write_nodeset_data()
        var_names = {"var1", "var2", "var3"};
        node_boundary_ids =
          {
            {0, 2}, // var1 is defined on nodesets 0 and 2
            {1, 3}, // var2 is defined on nodesets 1 and 3
            {4}     // var3 is only defined on the empty nodeset 4
          };

        // Data structure mapping (node, id) tuples to Real values that
        // will be passed to Exodus.
        bc_vals.resize(var_names.size());

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
      }

    // Make sure that the writing is done before the reading starts.
    TestCommWorld->barrier();

    // Now read it back in
    Mesh read_mesh(*TestCommWorld);
    IOClass reader(read_mesh);
    reader.read(filename);

    if (write_vars)
      {
        std::vector<std::string> read_in_var_names;
        std::vector<std::set<boundary_id_type>> read_in_node_boundary_ids;
        std::vector<std::map<BoundaryInfo::NodeBCTuple, Real>> read_in_bc_vals;
        reader.read_nodeset_data
          (/*timestep=*/1, read_in_var_names, read_in_node_boundary_ids, read_in_bc_vals);

        // Assert that we got back out what we put in.
        CPPUNIT_ASSERT(read_in_var_names == var_names);
        CPPUNIT_ASSERT(read_in_node_boundary_ids == node_boundary_ids);
        CPPUNIT_ASSERT(read_in_bc_vals == bc_vals);
      } // if (write_vars)

    // Also check that the flat indices match those in the file
    std::map<BoundaryInfo::NodeBCTuple, unsigned int> bc_array_indices;
    reader.get_nodeset_data_indices(bc_array_indices);

    // Debugging
    // for (const auto & pr : bc_array_indices)
    //   {
    //     const auto & t = pr.first;
    //     const auto & node_id = std::get<0>(t);
    //     const auto & boundary_id = std::get<1>(t);
    //     libMesh::out << "(node, boundary_id) = "
    //                  << "(" << node_id << ", " << boundary_id << ")"
    //                  << " is at array index " << pr.second
    //                  << std::endl;
    //   }

    // For this test case, the nodeset arrays are ordered as follows:
    // ns_prop1 = 0, 1, 2, 3, 4 ;
    // ns_names = "bottom", "right", "top", "left", "empty" ;
    //
    // Exodus (1-based) node ids
    // node_ns1 = 1, 2, 3, 4, 5, 6 ;
    // node_ns2 = 6, 12, 18, 24, 30, 36 ;
    // node_ns3 = 31, 32, 33, 34, 35, 36 ;
    // node_ns4 = 1, 7, 13, 19, 25, 31 ;

    // Check the node ids in the "bottom" nodeset, which contains the
    // first six nodes in order.
    for (unsigned int i=0; i<6; ++i)
      CPPUNIT_ASSERT_EQUAL
        (static_cast<unsigned int>(i),
         libmesh_map_find(bc_array_indices, std::make_tuple(/*node_id=*/cast_int<dof_id_type>(i), /*b_id=*/0)));

    // Check the node ids in the "right" sideset, which contains every
    // 6th node starting from 5, i.e. 5, 11, 17
    for (unsigned int i=0; i<6; ++i)
      CPPUNIT_ASSERT_EQUAL
        (static_cast<unsigned int>(i),
         libmesh_map_find(bc_array_indices, std::make_tuple(/*node_id=*/cast_int<dof_id_type>(6*i + 5), /*b_id=*/1)));

    // Check the node ids in the "top" sideset, which contains the last
    // six nodes in order.
    for (unsigned int i=0; i<6; ++i)
      CPPUNIT_ASSERT_EQUAL
        (static_cast<unsigned int>(i),
         libmesh_map_find(bc_array_indices, std::make_tuple(/*node_id=*/cast_int<dof_id_type>(30 + i), /*b_id=*/2)));

    // Check the node ids in the "left" sideset, which contains every
    // 6th node starting from 0, i.e. 0, 6, 12
    for (unsigned int i=0; i<6; ++i)
      CPPUNIT_ASSERT_EQUAL
        (static_cast<unsigned int>(i),
         libmesh_map_find(bc_array_indices, std::make_tuple(/*node_id=*/cast_int<dof_id_type>(6*i), /*b_id=*/3)));
  }

  void testWriteExodus()
  {
    LOG_UNIT_TEST;

    testWriteImpl<ExodusII_IO>("write_nodeset_data.e", /*write_vars=*/true);
    testWriteImpl<ExodusII_IO>("write_nodeset_data.e", /*write_vars=*/false);
  }

  void testWriteNemesis()
  {
    LOG_UNIT_TEST;

    // FIXME: Not yet implemented
    // testWriteImpl<Nemesis_IO>("write_nodeset_data.n");
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(WriteNodesetData);
