#include "test_comm.h"

#include <libmesh/elem.h>
#include <libmesh/enum_elem_type.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>

#include "libmesh_cppunit.h"

#include <memory>

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

    std::unique_ptr<Elem> test_elem = Elem::build(elem_type);

#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
#if LIBMESH_DIM > 1
    if (test_elem->infinite())
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
            elem->set_node(i) = node;
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
  }
};
