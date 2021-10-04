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

class WriteSidesetData : public CppUnit::TestCase
{
  /**
   * This test ensures you can write both vector and scalar variables
   */
public:
  CPPUNIT_TEST_SUITE(WriteSidesetData);

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
  void testWriteImpl(const std::string & filename)
  {
    Mesh mesh(*TestCommWorld);

    // We set our initial conditions based on build_square node ids
    mesh.allow_renumbering(false);

    MeshTools::Generation::build_square(mesh,
                                        /*nx=*/5, /*ny=*/5,
                                        -1., 1.,
                                        -1., 1.,
                                        QUAD4);

    // Add an empty sideset
    mesh.get_boundary_info().sideset_name(4) = "empty";

    // Get list of all (elem, side, id) tuples
    std::vector<BoundaryInfo::BCTuple> all_bc_tuples =
      mesh.get_boundary_info().build_side_list();

    // Data structures to be passed to ExodusII_IO::write_sideset_data
    std::vector<std::string> var_names = {"var1", "var2", "var3"};
    std::vector<std::set<boundary_id_type>> side_ids =
      {
        {0, 2}, // var1 is defined on sidesets 0 and 2
        {1, 3}, // var2 is defined on sidesets 1 and 3
        {4}     // var3 is only defined on the empty sideset 4
      };

    // Data structure mapping (elem, side, id) tuples to Real values that
    // will be passed to Exodus.
    std::vector<std::map<BoundaryInfo::BCTuple, Real>> bc_vals(var_names.size());

    // For each var_names[i], construct bc_vals[i]
    for (unsigned int i=0; i<var_names.size(); ++i)
      {
        // const auto & var_name = var_names[i];
        auto & vals = bc_vals[i];

        for (const auto & t : all_bc_tuples)
          {
            // dof_id_type elem_id = std::get<0>(t);
            // unsigned int side_id = std::get<1>(t);
            boundary_id_type b_id = std::get<2>(t);

            if (side_ids[i].count(b_id))
              {
                // Compute a value. This could in theory depend on
                // var_name, elem_id, side_id, and/or b_id.
                Real val = static_cast<Real>(b_id);

                // Insert into the vals map.
                vals.emplace(t, val);
              }
          }

        // If we have a distributed mesh, write_sideset_data wants our
        // ghost data too; we'll just serialize everything here.
        if (!mesh.is_serial())
          TestCommWorld->set_union(vals);

      } // done constructing bc_vals

    // We write the file in the ExodusII format.
    {
      IOClass writer(mesh);
      writer.write(filename);
      writer.write_sideset_data (/*timestep=*/1, var_names, side_ids, bc_vals);
    }

    // Make sure that the writing is done before the reading starts.
    TestCommWorld->barrier();

    // Now read it back in
    Mesh read_mesh(*TestCommWorld);
    IOClass reader(read_mesh);
    reader.read(filename);

    std::vector<std::string> read_in_var_names;
    std::vector<std::set<boundary_id_type>> read_in_side_ids;
    std::vector<std::map<BoundaryInfo::BCTuple, Real>> read_in_bc_vals;
    reader.read_sideset_data
      (/*timestep=*/1, read_in_var_names, read_in_side_ids, read_in_bc_vals);

    // Assert that we got back out what we put in.
    CPPUNIT_ASSERT(read_in_var_names == var_names);
    CPPUNIT_ASSERT(read_in_side_ids == side_ids);
    CPPUNIT_ASSERT(read_in_bc_vals == bc_vals);

    // Also check that the flat indices match those in the file
    std::map<BoundaryInfo::BCTuple, unsigned int> bc_array_indices;
    reader.get_sideset_data_indices
      (/*timestep=*/1, bc_array_indices);

    // Debugging
    // for (const auto & pr : bc_array_indices)
    //   {
    //     const auto & t = pr.first;
    //     const auto & elem_id = std::get<0>(t);
    //     const auto & side_id = std::get<1>(t);
    //     const auto & boundary_id = std::get<2>(t);
    //     libMesh::out << "(elem, side, boundary_id) = "
    //                  << "(" << elem_id << ", " << side_id << ", " << boundary_id << ")"
    //                  << std::endl;
    //   }

    // For this test case, the sideset arrays are ordered as follows:
    // elem_ss1 = 1, 2, 3, 4, 5 ;
    // side_ss1 = 1, 1, 1, 1, 1 ;
    // elem_ss2 = 5, 10, 15, 20, 25 ;
    // side_ss2 = 2, 2, 2, 2, 2 ;
    // elem_ss3 = 21, 22, 23, 24, 25 ;
    // side_ss3 = 3, 3, 3, 3, 3 ;
    // elem_ss4 = 1, 6, 11, 16, 21 ;
    // side_ss4 = 4, 4, 4, 4, 4 ;

    // Check that the 0th side of the 0th element on boundary 0
    // corresponds to array index 0.
    // CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(0),
    //                      libmesh_map_find(bc_array_indices,
    //                                       std::make_tuple(/*elem id*/static_cast<dof_id_type>(0),
    //                                                       /*side id*/static_cast<unsigned short int>(0),
    //                                                       /*boundary id*/static_cast<boundary_id_type>(0))));

    // CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(0),
    //                      libmesh_map_find(bc_array_indices, std::make_tuple(0,0,0)));
    //
    // CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(1),
    //                      libmesh_map_find(bc_array_indices, std::make_tuple(1,0,0)));

    // Check the first five elements, which all have side 0 on
    // boundary 0, and which happen to be in order in the BC array.
    for (unsigned int i=0; i<5; ++i)
      CPPUNIT_ASSERT_EQUAL
        (static_cast<unsigned int>(i),
         libmesh_map_find(bc_array_indices,
                          std::make_tuple(/*elem_id=*/cast_int<dof_id_type>(i),
                                          /*side_id=*/0,
                                          /*b_id=*/0)));

    // Check side 1 of every fifth element starting with element 4, they are all in sideset 1.
    for (unsigned int i=0; i<5; ++i)
      CPPUNIT_ASSERT_EQUAL
        (static_cast<unsigned int>(i),
         libmesh_map_find(bc_array_indices,
                          std::make_tuple(/*elem_id=*/cast_int<dof_id_type>(5*i + 4),
                                          /*side_id*/1,
                                          /*b_id=*/1)));

    // Check side 2 of the 5 consecutive elements starting with Elem 20. They are all in sideset 2.
    for (unsigned int i=0; i<5; ++i)
      CPPUNIT_ASSERT_EQUAL
        (static_cast<unsigned int>(i),
         libmesh_map_find(bc_array_indices,
                          std::make_tuple(/*elem_id=*/cast_int<dof_id_type>(20+i),
                                          /*side_id*/2,
                                          /*b_id=*/2)));

    // Check side 3 of every fifth element, they are all in sideset 3
    for (unsigned int i=0; i<5; ++i)
      CPPUNIT_ASSERT_EQUAL
        (static_cast<unsigned int>(i),
         libmesh_map_find(bc_array_indices,
                          std::make_tuple(/*elem_id=*/cast_int<dof_id_type>(5*i),
                                          /*side_id*/3,
                                          /*b_id=*/3)));
  }

  void testWriteExodus()
  {
    testWriteImpl<ExodusII_IO>("write_sideset_data.e");
  }

  void testWriteNemesis()
  {
    // FIXME: Not yet implemented
    // testWriteImpl<Nemesis_IO>("write_sideset_data.n");
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION(WriteSidesetData);
