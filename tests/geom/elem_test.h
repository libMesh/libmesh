#include "test_comm.h"

#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>

// Avoiding Elem::build() when we can't pass sides or n_sides to it
#include <libmesh/cell_c0polyhedron.h>
#include <libmesh/face_c0polygon.h>

#include "libmesh_cppunit.h"

#include <memory>
#include <sstream>

using namespace libMesh;

template <ElemType elem_type>
class PerElemTest : public CppUnit::TestCase
{
protected:
  std::unique_ptr<Mesh> _mesh;
  std::string libmesh_suite_name;

public:
  void setUp()
  {
    const Real minpos = 1.5, maxpos = 5.5;
    const unsigned int N = 2;

    _mesh = std::make_unique<Mesh>(*TestCommWorld);

    std::unique_ptr<Elem> test_elem;

    if (elem_type != C0POLYHEDRON)
      test_elem = Elem::build(elem_type);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
#if LIBMESH_DIM > 1
    if (test_elem.get() && test_elem->infinite())
      {
        Elem * elem = _mesh->add_elem(std::move(test_elem));

        const auto add_point =
          [this, elem](const unsigned int i,
                       const Real x,
                       const Real y,
                       const Real
#if LIBMESH_DIM == 3
                                  z
#endif
                      )
          {
#if LIBMESH_DIM == 2
            auto node = _mesh->add_point(Point(x, y), i);
#else
            auto node = _mesh->add_point(Point(x, y, z), i);
#endif
            elem->set_node(i, node);
          };

        const Real halfpos = (minpos + maxpos) / 2.;

        if (elem_type == INFQUAD4 || elem_type == INFQUAD6 ||
            elem_type == INFHEX8 || elem_type == INFHEX16 || elem_type == INFHEX18)
          {
            const bool is_quad = (elem_type == INFQUAD4 || elem_type == INFQUAD6);

            add_point(0, minpos, minpos, minpos);
            add_point(1, maxpos, minpos, minpos);
            add_point(2+is_quad, maxpos, maxpos, minpos);
            add_point(3-is_quad, minpos, maxpos, minpos);

            if (elem_type == INFQUAD6)
              {
                add_point(4, halfpos, minpos, minpos);
                add_point(5, halfpos, maxpos, minpos);
              }
          }
        if (elem_type == INFHEX8 || elem_type == INFHEX16 || elem_type == INFHEX18)
          {
            add_point(4, minpos, minpos, maxpos);
            add_point(5, maxpos, minpos, maxpos);
            add_point(6, maxpos, maxpos, maxpos);
            add_point(7, minpos, maxpos, maxpos);

            if (elem_type == INFHEX16 || elem_type == INFHEX18)
              {
                add_point(8, halfpos, minpos, minpos);
                add_point(9, maxpos, halfpos, minpos);
                add_point(10, halfpos, maxpos, minpos);
                add_point(11, minpos, halfpos, minpos);
                add_point(12, halfpos, minpos, maxpos);
                add_point(13, maxpos, halfpos, maxpos);
                add_point(14, halfpos, maxpos, maxpos);
                add_point(15, minpos, halfpos, maxpos);
              }
            if (elem_type == INFHEX18)
              {
                add_point(16, halfpos, halfpos, minpos);
                add_point(17, halfpos, halfpos, maxpos);
              }
          }
        if (elem_type == INFPRISM6 || elem_type == INFPRISM12)
          {
            add_point(0, minpos, minpos, minpos);
            add_point(1, maxpos, minpos, minpos);
            add_point(2, halfpos, maxpos, minpos);
            add_point(3, minpos, minpos, maxpos);
            add_point(4, maxpos, minpos, maxpos);
            add_point(5, halfpos, maxpos, maxpos);

            if (elem_type == INFPRISM12)
              {
                add_point(6, halfpos, minpos, minpos);
                add_point(7, (halfpos + maxpos) / 2., halfpos, minpos);
                add_point(8, (halfpos + minpos) / 2., halfpos, minpos);
                add_point(9, halfpos, minpos, maxpos);
                add_point(10, (halfpos + maxpos) / 2., halfpos, maxpos);
                add_point(11, (halfpos + minpos) / 2., halfpos, maxpos);
              }
          }

        _mesh->prepare_for_use();
      }
    else
#endif // LIBMESH_DIM > 1
#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
    if (elem_type == C0POLYGON)
      {
        // We're not going to implement build_square for e.g.
        // pentagons any time soon.
        //
        // We should probably implement build_dual_mesh, though, and
        // run it on a perturbed or triangle input so we don't just
        // get quads in the output...

        _mesh->add_point(Point(0, 0), 0);
        _mesh->add_point(Point(1, 0), 1);
        _mesh->add_point(Point(1.5, 0.5), 2);
        _mesh->add_point(Point(1, 1), 3);
        _mesh->add_point(Point(0, 1), 4);

        std::unique_ptr<Elem> polygon = std::make_unique<C0Polygon>(5);
        for (auto i : make_range(5))
          polygon->set_node(i, _mesh->node_ptr(i));
        polygon->set_id() = 0;

        _mesh->add_elem(std::move(polygon));
        _mesh->prepare_for_use();
      }
    else if (elem_type == C0POLYHEDRON)
      {
        // There's even less point in having a build_cube for general
        // polyhedra, so again we'll hand-make one for testing and
        // we'll plan on making a build_dual_mesh for the future.

#if 0
        // Try a truncated cube for relatively simple verification.

        // We have some differing orientations on the top sides, to
        // test handling of that.
        // Lower octagon points, counterclockwise
        _mesh->add_point(Point(1/Real(3), 0, 0), 0);
        _mesh->add_point(Point(2/Real(3), 0, 0), 1);
        _mesh->add_point(Point(1, 1/Real(3), 0), 2);
        _mesh->add_point(Point(1, 2/Real(3), 0), 3);
        _mesh->add_point(Point(2/Real(3), 1, 0), 4);
        _mesh->add_point(Point(1/Real(3), 1, 0), 5);
        _mesh->add_point(Point(0, 2/Real(3), 0), 6);
        _mesh->add_point(Point(0, 1/Real(3), 0), 7);

        // Lower-middle points, counterclockwise
        _mesh->add_point(Point(0, 0, 1/Real(3)), 8);
        _mesh->add_point(Point(1, 0, 1/Real(3)), 9);
        _mesh->add_point(Point(1, 1, 1/Real(3)), 10);
        _mesh->add_point(Point(0, 1, 1/Real(3)), 11);

        // Upper-middle points, counterclockwise
        _mesh->add_point(Point(0, 0, 2/Real(3)), 12);
        _mesh->add_point(Point(1, 0, 2/Real(3)), 13);
        _mesh->add_point(Point(1, 1, 2/Real(3)), 14);
        _mesh->add_point(Point(0, 1, 2/Real(3)), 15);

        // Upper octagon points, counterclockwise
        _mesh->add_point(Point(1/Real(3), 0, 1), 16);
        _mesh->add_point(Point(2/Real(3), 0, 1), 17);
        _mesh->add_point(Point(1, 1/Real(3), 1), 18);
        _mesh->add_point(Point(1, 2/Real(3), 1), 19);
        _mesh->add_point(Point(2/Real(3), 1, 1), 20);
        _mesh->add_point(Point(1/Real(3), 1, 1), 21);
        _mesh->add_point(Point(0, 2/Real(3), 1), 22);
        _mesh->add_point(Point(0, 1/Real(3), 1), 23);

        const std::vector<std::vector<unsigned int>> nodes_on_side =
          { {0, 1, 2, 3, 4, 5, 6, 7},         // min z
            {0, 1, 9, 13, 17, 16, 12, 8},     // min y
            {2, 3, 10, 14, 19, 18, 13, 9},    // max x
            {4, 5, 11, 15, 21, 20, 14, 10},   // max y
            {6, 7, 8, 12, 23, 22, 15, 11},    // min x
            {16, 17, 18, 19, 20, 21, 22, 23}, // max z
            {7, 0, 8},                        // max nothing
            {1, 2, 9},                        // max x
            {3, 4, 10},                       // max xy
            {5, 6, 11},                       // max y
            {23, 16, 12},                     // max z
            {17, 18, 13},                     // max xz
            {19, 20, 14},                     // max xyz
            {21, 22, 15} };                   // max yz
#endif

#if 1
        // Or try a plain cube, for simpler debugging
        _mesh->add_point(Point(0, 0, 0), 0);
        _mesh->add_point(Point(1, 0, 0), 1);
        _mesh->add_point(Point(1, 1, 0), 2);
        _mesh->add_point(Point(0, 1, 0), 3);
        _mesh->add_point(Point(0, 0, 1), 4);
        _mesh->add_point(Point(1, 0, 1), 5);
        _mesh->add_point(Point(1, 1, 1), 6);
        _mesh->add_point(Point(0, 1, 1), 7);

        // With some combinations of face triangulations, even a
        // simple cube has no tetrahedralization that doesn't have
        // either interior discontinuities or a 0-volume pancake tet!
        //
        // The initial "natural" way to orient our sides is commented
        // out here; permutations that fix diagonals for us are used
        // instead.
        const std::vector<std::vector<unsigned int>> nodes_on_side =
          { {0, 1, 2, 3},   // min z
            {0, 1, 5, 4},   // min y
        //     {1, 2, 6, 5},   // max x - bad
            {2, 6, 5, 1},   // max x
            {2, 3, 7, 6},   // max y
        //     {3, 0, 4, 7},   // min x - bad
            {0, 4, 7, 3},   // min x
        //     {4, 5, 6, 7} }; // max z - bad
            {5, 6, 7, 4} }; // max z
#endif

        // Build all the sides.
        std::vector<std::shared_ptr<Polygon>> sides(nodes_on_side.size());

        for (auto s : index_range(nodes_on_side))
          {
            const auto & nodes_on_s = nodes_on_side[s];
            sides[s] = std::make_shared<C0Polygon>(nodes_on_s.size());
            for (auto i : index_range(nodes_on_s))
              sides[s]->set_node(i, _mesh->node_ptr(nodes_on_s[i]));
          }

        std::unique_ptr<Elem> polyhedron = std::make_unique<C0Polyhedron>(sides);
        _mesh->add_elem(std::move(polyhedron));
        _mesh->prepare_for_use();
      }
    else
      {
        const unsigned int dim = test_elem->dim();
        const unsigned int use_x = dim > 0;
        const unsigned int use_y = dim > 1;
        const unsigned int use_z = dim > 2;

        MeshTools::Generation::build_cube (*_mesh,
                                           N*use_x, N*use_y, N*use_z,
                                           minpos, maxpos,
                                           minpos, use_y*maxpos,
                                           minpos, use_z*maxpos,
                                           elem_type);
      }

    // Use non-default subdomain ids so we can properly test their
    // preservation.
    //
    // Use long subdomain names so we can test that we're not
    // truncating as badly as we used to in ExodusII.
    for (const auto & elem :
         this->_mesh->element_ptr_range())
    {
      const subdomain_id_type sbdid = 10 + (elem->id() % 10);
      elem->subdomain_id() = sbdid;
      std::ostringstream sbdname;
      sbdname <<
        "a_very_long_subdomain_name_for_the_subdomain_with_number_" << sbdid;
      this->_mesh->subdomain_name(sbdid) = sbdname.str();
    }

    // Make sure our mesh's cache knows about them all for later
    // test comparisons
    this->_mesh->cache_elem_data();

    // We may have updated only a portion of subdomain names if we're
    // on a distributed mesh, but every processor should know about
    // every name.
    this->_mesh->sync_subdomain_name_map();
  }
};
