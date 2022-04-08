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
