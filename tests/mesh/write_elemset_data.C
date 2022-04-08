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
#include "libmesh/elem.h"

#include "test_comm.h"
#include "libmesh_cppunit.h"


// Bring in everything from the libMesh namespace
using namespace libMesh;

class WriteElemsetData : public CppUnit::TestCase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE(WriteElemsetData);

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

    // Allocate space for an extra integer on each element to store a "code" which
    // determines which sets an Elem belongs to. We do this before building the Mesh.
    //
    // extra_integer val               & sets elem belongs to
    // DofObject::invalid_id (default) & elem belongs to no sets
    // 1                               & elem belongs to set A only
    // 2                               & elem belongs to set B only
    // 3                               & elem belongs to sets A and B
    unsigned int elemset_index =
      mesh.add_elem_integer("elemset_code",
                            /*allocate_data=*/true);

    // We are generating this mesh, so it should not be renumbered.
    // No harm in being explicit about it, however.
    mesh.allow_renumbering(false);

    MeshTools::Generation::build_square(mesh,
                                        /*nx=*/5, /*ny=*/5,
                                        -1., 1.,
                                        -1., 1.,
                                        QUAD4);

    // Set ids for elements in elemsets 1, 2
    std::set<dof_id_type> set1 = {3, 8, 14, 24};
    std::set<dof_id_type> set2 = {3, 9, 15, 24};

    // Loop over Elems and set "elemset_code" values
    for (const auto & elem : mesh.element_ptr_range())
      {
        bool
          in1 = set1.count(elem->id()),
          in2 = set2.count(elem->id());

        dof_id_type val = DofObject::invalid_id;
        if (in1)
          val = 1;
        if (in2)
          val = 2;
        if (in1 && in2)
          val = 3;

        elem->set_extra_integer(elemset_index, val);
      }

    // Tell the Mesh about these elemsets
    mesh.add_elemset_code(1, {1});
    mesh.add_elemset_code(2, {2});
    mesh.add_elemset_code(3, {1,2});

    // Debugging: print valid elemset_codes values
    for (const auto & elem : mesh.element_ptr_range())
      {
        dof_id_type elemset_code =
          elem->get_extra_integer(elemset_index);

        if (elemset_code != DofObject::invalid_id)
          libMesh::out << "Elem " << elem->id() << ", elemset_code = " << elemset_code << std::endl;
      }

    // Write the file in the ExodusII format, including the element set information.
    // Note: elemsets should eventually be written during ExodusII_IO::write(), this
    // would match the behavior of sidesets and nodesets.
    {
      IOClass writer(mesh);
      writer.write(filename);
      writer.write_elemsets();
    }

   // Make sure that the writing is done before the reading starts.
   TestCommWorld->barrier();

   // Now read it back in
   Mesh read_mesh(*TestCommWorld);
   IOClass reader(read_mesh);
   reader.verbose(true); // additional messages while debugging
   reader.read(filename);

//    if (write_vars)
//      {
//        std::vector<std::string> read_in_var_names;
//        std::vector<std::set<boundary_id_type>> read_in_side_ids;
//        std::vector<std::map<BoundaryInfo::BCTuple, Real>> read_in_bc_vals;
//        reader.read_sideset_data
//          (/*timestep=*/1, read_in_var_names, read_in_side_ids, read_in_bc_vals);
//
//        // Assert that we got back out what we put in.
//        CPPUNIT_ASSERT(read_in_var_names == var_names);
//        CPPUNIT_ASSERT(read_in_side_ids == side_ids);
//        CPPUNIT_ASSERT(read_in_bc_vals == bc_vals);
//      } // if (write_vars)
//
//    // Also check that the flat indices match those in the file
//    std::map<BoundaryInfo::BCTuple, unsigned int> bc_array_indices;
//    reader.get_sideset_data_indices(bc_array_indices);
//
//    // Debugging
//    // for (const auto & pr : bc_array_indices)
//    //   {
//    //     const auto & t = pr.first;
//    //     const auto & elem_id = std::get<0>(t);
//    //     const auto & side_id = std::get<1>(t);
//    //     const auto & boundary_id = std::get<2>(t);
//    //     libMesh::out << "(elem, side, boundary_id) = "
//    //                  << "(" << elem_id << ", " << side_id << ", " << boundary_id << ")"
//    //                  << std::endl;
//    //   }
//
//    // For this test case, the sideset arrays are ordered as follows:
//    // elem_ss1 = 1, 2, 3, 4, 5 ;
//    // side_ss1 = 1, 1, 1, 1, 1 ;
//    // elem_ss2 = 5, 10, 15, 20, 25 ;
//    // side_ss2 = 2, 2, 2, 2, 2 ;
//    // elem_ss3 = 21, 22, 23, 24, 25 ;
//    // side_ss3 = 3, 3, 3, 3, 3 ;
//    // elem_ss4 = 1, 6, 11, 16, 21 ;
//    // side_ss4 = 4, 4, 4, 4, 4 ;
//
//    // Check that the 0th side of the 0th element on boundary 0
//    // corresponds to array index 0.
//    // CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(0),
//    //                      libmesh_map_find(bc_array_indices,
//    //                                       std::make_tuple(/*elem id*/static_cast<dof_id_type>(0),
//    //                                                       /*side id*/static_cast<unsigned short int>(0),
//    //                                                       /*boundary id*/static_cast<boundary_id_type>(0))));
//
//    // CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(0),
//    //                      libmesh_map_find(bc_array_indices, std::make_tuple(0,0,0)));
//    //
//    // CPPUNIT_ASSERT_EQUAL(static_cast<unsigned int>(1),
//    //                      libmesh_map_find(bc_array_indices, std::make_tuple(1,0,0)));
//
//    // Check the first five elements, which all have side 0 on
//    // boundary 0, and which happen to be in order in the BC array.
//    for (unsigned int i=0; i<5; ++i)
//      CPPUNIT_ASSERT_EQUAL
//        (static_cast<unsigned int>(i),
//         libmesh_map_find(bc_array_indices,
//                          std::make_tuple(/*elem_id=*/cast_int<dof_id_type>(i),
//                                          /*side_id=*/0,
//                                          /*b_id=*/0)));
//
//    // Check side 1 of every fifth element starting with element 4, they are all in sideset 1.
//    for (unsigned int i=0; i<5; ++i)
//      CPPUNIT_ASSERT_EQUAL
//        (static_cast<unsigned int>(i),
//         libmesh_map_find(bc_array_indices,
//                          std::make_tuple(/*elem_id=*/cast_int<dof_id_type>(5*i + 4),
//                                          /*side_id*/1,
//                                          /*b_id=*/1)));
//
//    // Check side 2 of the 5 consecutive elements starting with Elem 20. They are all in sideset 2.
//    for (unsigned int i=0; i<5; ++i)
//      CPPUNIT_ASSERT_EQUAL
//        (static_cast<unsigned int>(i),
//         libmesh_map_find(bc_array_indices,
//                          std::make_tuple(/*elem_id=*/cast_int<dof_id_type>(20+i),
//                                          /*side_id*/2,
//                                          /*b_id=*/2)));
//
//    // Check side 3 of every fifth element, they are all in sideset 3
//    for (unsigned int i=0; i<5; ++i)
//      CPPUNIT_ASSERT_EQUAL
//        (static_cast<unsigned int>(i),
//         libmesh_map_find(bc_array_indices,
//                          std::make_tuple(/*elem_id=*/cast_int<dof_id_type>(5*i),
//                                          /*side_id*/3,
//                                          /*b_id=*/3)));
  }

  void testWriteExodus()
  {
    LOG_UNIT_TEST;

    testWriteImpl<ExodusII_IO>("write_elemset_data.e");
  }

  void testWriteNemesis()
  {
    // LOG_UNIT_TEST;

    // FIXME: Not yet implemented
    // testWriteImpl<Nemesis_IO>("write_elemset_data.n");
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION(WriteElemsetData);
