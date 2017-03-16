// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/edge_edge2.h>
#include <libmesh/face_quad4.h>
#include <libmesh/face_tri3.h>
#include <libmesh/cell_hex8.h>
#include <libmesh/dof_map.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/mesh_function.h>
#include <libmesh/numeric_vector.h>

#include "test_comm.h"

// THE CPPUNIT_TEST_SUITE_END macro expands to code that involves
// std::auto_ptr, which in turn produces -Wdeprecated-declarations
// warnings.  These can be ignored in GCC as long as we wrap the
// offending code in appropriate pragmas.  We can't get away with a
// single ignore_warnings.h inclusion at the beginning of this file,
// since the libmesh headers pull in a restore_warnings.h at some
// point.  We also don't bother restoring warnings at the end of this
// file since it's not a header.
#include <libmesh/ignore_warnings.h>

using namespace libMesh;

Number position_function (const Point& p,
                          const Parameters&,
                          const std::string&,
                          const std::string&)
{
  if (p(1) > 0.0)
    return 0.0;
  else
    return 1.0;
}

Number position_function2 (const Point& p,
                           const Parameters&,
                           const std::string&,
                           const std::string&)
{
  if (p(1) > 0.0)
    return 2.0 * p(1) + 1.0;
  else
    return p(1) + 1.0;
}

class MeshfunctionDFEM : public CppUnit::TestCase
{
  /**
   * The goal of this test is to verify proper operation of mesh function
   * for discontinuous shape functions and points close to the boundary
   */
public:
  CPPUNIT_TEST_SUITE( MeshfunctionDFEM );

  CPPUNIT_TEST( test_point_locator_dfem );

  CPPUNIT_TEST( test_mesh_function_dfem );

  CPPUNIT_TEST( test_mesh_function_dfem_grad );

  CPPUNIT_TEST_SUITE_END();

protected:
  ReplicatedMesh * _mesh;
  UniquePtr<PointLocatorBase> _point_locator;

  void build_mesh()
  {
    _mesh = new ReplicatedMesh(*TestCommWorld);

    // (0,1)           (1,1)
    // x---------------x
    // |               |
    // |               |
    // |               |
    // |               |
    // |               |
    // x---------------x
    // (0,0)           (1,0)
    // |               |
    // |               |
    // |               |
    // |               |
    // x---------------x
    // (0,-1)          (1,-1)

    _mesh->set_mesh_dimension(2);

    _mesh->add_point( Point(0.0,-1.0), 4 );
    _mesh->add_point( Point(1.0,-1.0), 5 );
    _mesh->add_point( Point(1.0, 0.0), 1 );
    _mesh->add_point( Point(1.0, 1.0), 2 );
    _mesh->add_point( Point(0.0, 1.0), 3 );
    _mesh->add_point( Point(0.0, 0.0), 0 );

    {
      Elem* elem_top = _mesh->add_elem( new Quad4 );
      elem_top->set_node(0) = _mesh->node_ptr(0);
      elem_top->set_node(1) = _mesh->node_ptr(1);
      elem_top->set_node(2) = _mesh->node_ptr(2);
      elem_top->set_node(3) = _mesh->node_ptr(3);

      Elem* elem_bottom = _mesh->add_elem( new Quad4 );
      elem_bottom->set_node(0) = _mesh->node_ptr(4);
      elem_bottom->set_node(1) = _mesh->node_ptr(5);
      elem_bottom->set_node(2) = _mesh->node_ptr(1);
      elem_bottom->set_node(3) = _mesh->node_ptr(0);
    }

    // libMesh will renumber, but we numbered according to its scheme
    // anyway. We do this because when we call uniformly_refine subsequenly,
    // it's going use skip_renumber=false.
    _mesh->prepare_for_use(false /*skip_renumber*/);

    // get a point locator
    _point_locator = _mesh->sub_point_locator();
  }

public:
  void setUp()
  {
    this->build_mesh();
  }

  void tearDown()
  {
    delete _mesh;
  }

  // test that point locator works correctly
  void test_point_locator_dfem()
  {
    // this point is in bottom element only
    Point interior(0.5, -0.5, 0.0);

    // this point is on the face between bottom and top
    Point face(0.5, 0.0, 0.0);

    // test interior point
    std::set<const Elem *> int_cand;
    (*_point_locator)(interior, int_cand, libmesh_nullptr);
    CPPUNIT_ASSERT (int_cand.size() == 1);
    for (std::set<const Elem *>::iterator it = int_cand.begin(); it != int_cand.end(); ++it)
      CPPUNIT_ASSERT ((*it)->id() == 1);

    // test interior point
    std::set<const Elem *> face_cand;
    (*_point_locator)(face, face_cand, libmesh_nullptr);
    CPPUNIT_ASSERT (face_cand.size() == 2);
    int array[2] = {0, 0};
    for (std::set<const Elem *>::iterator it = face_cand.begin(); it != face_cand.end(); ++it)
      {
        // first test that array entry hasn't been set before
        CPPUNIT_ASSERT (array[(*it)->id()] == 0);
        array[(*it)->id()] = 1;
      }
    CPPUNIT_ASSERT (array[0] == 1);
    CPPUNIT_ASSERT (array[1] == 1);
  }

  // test that mesh function works correctly
  void test_mesh_function_dfem()
  {
    // this point is in top element only
    Point top(0.5, 0.5, 0.0);

    // this point is in bottom element only
    Point bottom(0.5, -0.5, 0.0);

    // this point is on the face between bottom and top
    Point face(0.5, 0.0, 0.0);

    // set up some necessary objects
    EquationSystems es(*_mesh);
    System & sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", CONSTANT, MONOMIAL);

    es.init();
    sys.project_solution(position_function, NULL, es.parameters);

    std::vector<unsigned int> variables;
    sys.get_all_variable_numbers(variables);
    UniquePtr< NumericVector<Number> > mesh_function_vector =
      NumericVector<Number>::build(es.comm());
    mesh_function_vector->init(sys.n_dofs(), false, SERIAL);
    sys.solution->localize(*mesh_function_vector);

    MeshFunction mesh_function (es, *mesh_function_vector,
                                sys.get_dof_map(), variables);
    mesh_function.init();

    // test mesh function in top
    std::map<const Elem *, Number> top_val = mesh_function.discontinuous_value(top);

    // check that there is only one value
    CPPUNIT_ASSERT (top_val.size() == 1);

    // check that this one value is the right one
    for (std::map<const Elem *, Number>::const_iterator it = top_val.begin(); it != top_val.end(); ++it)
      CPPUNIT_ASSERT (it->first->id() == 0 && std::abs(it->second) < 1.0e-10);

    // test mesh function in bottom
    std::map<const Elem *, Number> bottom_val = mesh_function.discontinuous_value(bottom);

    // check that there is only one value
    CPPUNIT_ASSERT (bottom_val.size() == 1);

    // check that this one value is the right one
    for (std::map<const Elem *, Number>::const_iterator it = bottom_val.begin(); it != bottom_val.end(); ++it)
      CPPUNIT_ASSERT (it->first->id() == 1 && std::abs(it->second - 1.0) < 1.0e-10);

    // test mesh function in face
    std::map<const Elem *, Number> face_val = mesh_function.discontinuous_value(face);

    // check that there are two values
    CPPUNIT_ASSERT (face_val.size() == 2);

    // check that the two values are attached to the right element
    for (std::map<const Elem *, Number>::const_iterator it = face_val.begin(); it != face_val.end(); ++it)
      if (it->first->id() == 0)
        CPPUNIT_ASSERT (std::abs(it->second) < 1.0e-10);
      else if (it->first->id() == 1)
        CPPUNIT_ASSERT (std::abs(it->second - 1.0) < 1.0e-10);
      else
        CPPUNIT_ASSERT (false);

  }

  // test that gradient function works correctly
  void test_mesh_function_dfem_grad()
  {
    // this point is in top element only
    Point top(0.5, 0.5, 0.0);

    // this point is in bottom element only
    Point bottom(0.5, -0.5, 0.0);

    // this point is on the face between bottom and top
    Point face(0.5, 0.0, 0.0);

    // set up some necessary objects
    EquationSystems es(*_mesh);
    System & sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", FIRST, LAGRANGE);

    es.init();
    sys.project_solution(position_function2, NULL, es.parameters);

    std::vector<unsigned int> variables;
    sys.get_all_variable_numbers(variables);
    UniquePtr< NumericVector<Number> > mesh_function_vector =
      NumericVector<Number>::build(es.comm());
    mesh_function_vector->init(sys.n_dofs(), false, SERIAL);
    sys.solution->localize(*mesh_function_vector);

    MeshFunction mesh_function (es, *mesh_function_vector,
                                sys.get_dof_map(), variables);
    mesh_function.init();

    // test mesh function in top
    std::map<const Elem *, Gradient> top_val = mesh_function.discontinuous_gradient(top);

    // check that there is only one value
    CPPUNIT_ASSERT (top_val.size() == 1);

    // check that this one value is the right one
    for (std::map<const Elem *, Gradient>::const_iterator it = top_val.begin(); it != top_val.end(); ++it)
      CPPUNIT_ASSERT (it->first->id() == 0 && std::abs(it->second(1) - 2.0) < 1.0e-10);

    // test mesh function in bottom
    std::map<const Elem *, Gradient> bottom_val = mesh_function.discontinuous_gradient(bottom);

    // check that there is only one value
    CPPUNIT_ASSERT (bottom_val.size() == 1);

    // check that this one value is the right one
    for (std::map<const Elem *, Gradient>::const_iterator it = bottom_val.begin(); it != bottom_val.end(); ++it)
      CPPUNIT_ASSERT (it->first->id() == 1 && std::abs(it->second(1) - 1.0) < 1.0e-10);

    // test mesh function in face
    std::map<const Elem *, Gradient> face_val = mesh_function.discontinuous_gradient(face);

    // check that there are two values
    CPPUNIT_ASSERT (face_val.size() == 2);

    // check that the two values are attached to the right element
    for (std::map<const Elem *, Gradient>::const_iterator it = face_val.begin(); it != face_val.end(); ++it)
      if (it->first->id() == 0)
        CPPUNIT_ASSERT (std::abs(it->second(1) - 2.0) < 1.0e-10);
      else if (it->first->id() == 1)
        CPPUNIT_ASSERT (std::abs(it->second(1) - 1.0) < 1.0e-10);
      else
        CPPUNIT_ASSERT (false);
  }
};


CPPUNIT_TEST_SUITE_REGISTRATION( MeshfunctionDFEM );
