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
    // for (const auto & elem : mesh.element_ptr_range())
    //   {
    //     dof_id_type elemset_code =
    //       elem->get_extra_integer(elemset_index);
    //
    //     if (elemset_code != DofObject::invalid_id)
    //       libMesh::out << "Elem " << elem->id() << ", elemset_code = " << elemset_code << std::endl;
    //   }

    // Set up variables defined on these elemsets
    std::vector<std::string> var_names = {"var1", "var2", "var3"};
    std::vector<std::set<elemset_id_type>> elemset_ids =
      {
        {1},  // var1 is defined on elemset 1
        {2},  // var2 is defined on elemset 2
        {1,2} // var3 is defined on elemsets 1 and 2
      };
    std::vector<std::map<std::pair<dof_id_type, elemset_id_type>, Real>> elemset_vals(var_names.size());

    // To catch values handed back by MeshBase::get_elemsets()
    std::set<elemset_id_type> id_set_to_fill;

    for (const auto & elem : mesh.element_ptr_range())
      {
        // Get list of elemset ids to which this element belongs
        mesh.get_elemsets(elem->get_extra_integer(elemset_index), id_set_to_fill);

        bool
          in1 = id_set_to_fill.count(1),
          in2 = id_set_to_fill.count(2);

        // Set the value for var1 == 1.0 on all elements in elemset 1
        if (in1)
          elemset_vals[/*var1 index=*/0].emplace( std::make_pair(elem->id(), /*elemset_id=*/1), 1.0);

        // Set the value of var2 == 2.0 on all elements in elemset 2
        if (in2)
          elemset_vals[/*var2 index=*/1].emplace( std::make_pair(elem->id(), /*elemset_id=*/2), 2.0);

        // Set the value of var3 == 3.0 on elements in the union of sets 1 and 2
        if (in1 || in2)
          for (const auto & id : id_set_to_fill)
            elemset_vals[/*var3 index=*/2].emplace( std::make_pair(elem->id(), /*elemset_id=*/id), 3.0);
      }

    // Sanity check: we should have 8 total elements in set1 and set2 combined
    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(8), elemset_vals[/*var3 index=*/2].size());

    // Lambda to help with debugging
    // auto print_map = [](const std::vector<std::map<std::pair<dof_id_type, elemset_id_type>, Real>> & input)
    //   {
    //     for (auto i : index_range(input))
    //       {
    //         libMesh::out << "Map " << i << " = " << std::endl;
    //         for (const auto & [pr, val] : input[i])
    //           {
    //             const auto & elem_id = pr.first;
    //             const auto & elemset_id = pr.second;
    //             libMesh::out << "(" << elem_id << ", " << elemset_id << ") = " << val << std::endl;
    //           }
    //       }
    //   };

    // Debugging: print the elemset_vals struct we just built
    // print_map(elemset_vals);

    // Write the file in the ExodusII format, including the element set information.
    // Note: elemsets should eventually be written during ExodusII_IO::write(), this
    // would match the behavior of sidesets and nodesets.
    {
      IOClass writer(mesh);
      writer.write(filename);
      writer.write_elemsets();
      writer.write_elemset_data(/*timestep=*/1, var_names, elemset_ids, elemset_vals);
    }

    // Make sure that the writing is done before the reading starts.
    TestCommWorld->barrier();

    // Now read it back in
    Mesh read_mesh(*TestCommWorld);

    // Do not allow renumbering on this mesh either.
    read_mesh.allow_renumbering(false);

    IOClass reader(read_mesh);
    // reader.verbose(true); // additional messages while debugging
    reader.read(filename);

    // Check that the elements in read_mesh are in the correct elemsets.
    // The elemset_codes will not in general match because they are
    // created by a generic algorithm in the Exodus reader while above
    // they were hard-coded.

    // Make sure that the mesh actually has an extra_integer for "elemset_code"
    CPPUNIT_ASSERT(read_mesh.has_elem_integer("elemset_code"));

    // Make sure the extra integer is in the same index as before
    CPPUNIT_ASSERT(read_mesh.get_elem_integer_index("elemset_code") == elemset_index);

    // Make sure the elemset_codes match what we are expecting.
    // The Exodus reader assigns the codes based on operator<
    // for std::sets, which gives us the ordering {1}, {1,2}, {2}
    CPPUNIT_ASSERT(read_mesh.get_elemset_code({1}) == 0);
    CPPUNIT_ASSERT(read_mesh.get_elemset_code({1,2}) == 1);
    CPPUNIT_ASSERT(read_mesh.get_elemset_code({2}) == 2);

    // Assert that the elemset_codes for particular elements are set as expected

    // Elements 8, 14 are in set1 which has code 0
    CPPUNIT_ASSERT(read_mesh.elem_ptr(8)->get_extra_integer(elemset_index) == 0);
    CPPUNIT_ASSERT(read_mesh.elem_ptr(14)->get_extra_integer(elemset_index) == 0);

    // Elements 3, 24 are in both set1 and set2, which has code 1
    CPPUNIT_ASSERT(read_mesh.elem_ptr(3)->get_extra_integer(elemset_index) == 1);
    CPPUNIT_ASSERT(read_mesh.elem_ptr(24)->get_extra_integer(elemset_index) == 1);

    // Elements 9, 15 are in set2 which has code 2
    CPPUNIT_ASSERT(read_mesh.elem_ptr(9)->get_extra_integer(elemset_index) == 2);
    CPPUNIT_ASSERT(read_mesh.elem_ptr(15)->get_extra_integer(elemset_index) == 2);

    // Read in the elemset variables from file
    std::vector<std::string> read_in_var_names;
    std::vector<std::set<elemset_id_type>> read_in_elemset_ids;
    std::vector<std::map<std::pair<dof_id_type, elemset_id_type>, Real>> read_in_elemset_vals;
    reader.read_elemset_data(/*timestep=*/1, read_in_var_names, read_in_elemset_ids, read_in_elemset_vals);

    // Debugging
    // print_map(read_in_elemset_vals);

    // Assert that the data we read in matches what we wrote out
    CPPUNIT_ASSERT(read_in_var_names == var_names);
    CPPUNIT_ASSERT(read_in_elemset_ids == elemset_ids);
    CPPUNIT_ASSERT_EQUAL(static_cast<std::size_t>(8), read_in_elemset_vals[/*var3 index=*/2].size());
    CPPUNIT_ASSERT(read_in_elemset_vals == elemset_vals);
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
