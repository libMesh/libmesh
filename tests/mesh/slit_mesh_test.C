// Ignore unused parameter warnings coming from cppuint headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/edge_edge2.h>
#include <libmesh/face_quad4.h>
#include <libmesh/cell_hex8.h>
#include <libmesh/dof_map.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh_refinement.h>

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

class SlitFunc : public FEMFunctionBase<Number>
{
public:

  SlitFunc() {}

  ~SlitFunc () {}

  virtual void init_context (const FEMContext &) libmesh_override {}

  virtual UniquePtr<FEMFunctionBase<Number> >
  clone () const libmesh_override
  {
    return UniquePtr<FEMFunctionBase<Number> > (new SlitFunc());
  }

  virtual Number operator() (const FEMContext & c,
                             const Point & p,
                             const Real /*time*/ = 0.)
  libmesh_override
  {
    const Real & x = p(0);
    const Real & y = p(1);
    const Point centroid = c.get_elem().centroid();
    const Real sign = centroid(1)/std::abs(centroid(1));

    return (1 - std::abs(1-x)) * (1-std::abs(y)) * sign;
  }

  virtual void operator() (const FEMContext & c,
                           const Point & p,
                           const Real time,
                           DenseVector<Number> & output)
  libmesh_override
  {
    for (unsigned int i=0; i != output.size(); ++i)
      output(i) = (*this)(c, p, time);
  }
};







class SlitMeshTest : public CppUnit::TestCase {
  /**
   * The goal of this test is to ensure that a 2D mesh with nodes overlapping
   * on opposite sides of an internal, "slit" edge is useable.
   */
public:
  CPPUNIT_TEST_SUITE( SlitMeshTest );

  CPPUNIT_TEST( testMesh );

  CPPUNIT_TEST_SUITE_END();

protected:

  Mesh* _mesh;

  void build_mesh()
  {
    _mesh = new Mesh(*TestCommWorld);

    /*
      (0,1)           (1,1)           (2,1)
        x---------------x---------------x
        |               |               |
        |               |               |
        |               |               |
        |               |               |
        |               |               |
        x---------------x---------------x
       (0,0)           (1,0)          (2,0)
        |               |               |
        |               |               |
        |               |               |
        |               |               |
        x---------------x---------------x
       (0,-1)          (1,-1)         (2,-1)
     */

    _mesh->set_mesh_dimension(2);

    _mesh->add_point( Point(0.0, 0.0), 0 );
    _mesh->add_point( Point(1.0, 0.0), 1 );
    _mesh->add_point( Point(1.0, 1.0), 2 );
    _mesh->add_point( Point(0.0, 1.0), 3 );
    _mesh->add_point( Point(0.0,-1.0), 4 );
    _mesh->add_point( Point(1.0,-1.0), 5 );
    _mesh->add_point( Point(1.0, 0.0), 6 );
    _mesh->add_point( Point(2.0, 0.0), 7 );
    _mesh->add_point( Point(2.0, 1.0), 8 );
    _mesh->add_point( Point(2.0,-1.0), 9 );

    {
      Elem* elem_top_left = new Quad4;
      elem_top_left->set_node(0) = _mesh->node_ptr(0);
      elem_top_left->set_node(1) = _mesh->node_ptr(1);
      elem_top_left->set_node(2) = _mesh->node_ptr(2);
      elem_top_left->set_node(3) = _mesh->node_ptr(3);
      elem_top_left->set_id() = 0;
      _mesh->add_elem(elem_top_left);

      Elem* elem_bottom_left = new Quad4;
      elem_bottom_left->set_node(0) = _mesh->node_ptr(4);
      elem_bottom_left->set_node(1) = _mesh->node_ptr(5);
      elem_bottom_left->set_node(2) = _mesh->node_ptr(6);
      elem_bottom_left->set_node(3) = _mesh->node_ptr(0);
      elem_bottom_left->set_id() = 1;
      _mesh->add_elem(elem_bottom_left);

      Elem* elem_top_right = new Quad4;
      elem_top_right->set_node(0) = _mesh->node_ptr(1);
      elem_top_right->set_node(1) = _mesh->node_ptr(7);
      elem_top_right->set_node(2) = _mesh->node_ptr(8);
      elem_top_right->set_node(3) = _mesh->node_ptr(2);
      elem_top_right->set_id() = 2;
      _mesh->add_elem(elem_top_right);

      Elem* elem_bottom_right = new Quad4;
      elem_bottom_right->set_node(0) = _mesh->node_ptr(5);
      elem_bottom_right->set_node(1) = _mesh->node_ptr(9);
      elem_bottom_right->set_node(2) = _mesh->node_ptr(7);
      elem_bottom_right->set_node(3) = _mesh->node_ptr(6);
      elem_bottom_right->set_id() = 3;
      _mesh->add_elem(elem_bottom_right);
    }

    // libMesh shouldn't renumber, or our based-on-initial-id
    // assertions later may fail.
    _mesh->allow_renumbering(false);

    _mesh->prepare_for_use();
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

  void testMesh()
  {
    // There'd better be 4 elements
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)4, _mesh->n_elem() );

    // There'd better still be a full 10 nodes
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)10, _mesh->n_nodes() );

    /* The middle nodes should still be distinct between the top and
     * bottom elements */
    if (_mesh->query_elem_ptr(0) && _mesh->query_elem_ptr(1))
      CPPUNIT_ASSERT( _mesh->elem_ref(0).node_id(1) != _mesh->elem_ref(1).node_id(2) );
    if (_mesh->query_elem_ptr(2) && _mesh->query_elem_ptr(3))
      CPPUNIT_ASSERT( _mesh->elem_ref(2).node_id(0) != _mesh->elem_ref(3).node_id(3) );

    /* The middle nodes should still be shared between left and right
     * elements on top and bottom */
    if (_mesh->query_elem_ptr(0) && _mesh->query_elem_ptr(2))
      CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(0).node_id(1),
                            _mesh->elem_ref(2).node_id(0) );
    if (_mesh->query_elem_ptr(1) && _mesh->query_elem_ptr(3))
      CPPUNIT_ASSERT_EQUAL( _mesh->elem_ref(1).node_id(2),
                            _mesh->elem_ref(3).node_id(3) );
  }

};

class SlitMeshRefinedMeshTest : public SlitMeshTest {
  /**
   * The goal of this test is the same as the previous, but now we do a
   * uniform refinement and make sure the result mesh is consistent. i.e.
   * the new node shared between the 1D elements is the same as the node
   * shared on the underlying quads, and so on.
   */
public:
  CPPUNIT_TEST_SUITE( SlitMeshRefinedMeshTest );

  CPPUNIT_TEST( testMesh );

  CPPUNIT_TEST_SUITE_END();

  // Yes, this is necessary. Somewhere in those macros is a protected/private
public:

  void setUp()
  {
    this->build_mesh();

#ifdef LIBMESH_ENABLE_AMR
    MeshRefinement(*_mesh).uniformly_refine(1);
#endif
  }

  void testMesh()
  {
#ifdef LIBMESH_ENABLE_AMR
    // We should have 20 total and 16 active elements.
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)20, _mesh->n_elem() );
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)16, _mesh->n_active_elem() );

    // We should have 28 nodes, not 25 or 26
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)28, _mesh->n_nodes() );
#endif
  }
};

class SlitMeshRefinedSystemTest : public SlitMeshTest {
  /**
   * The goal of this test is the same as the previous, but now we
   * create a system and set dof values to make sure they are properly
   * interpolated after refinement.
   */
public:
  CPPUNIT_TEST_SUITE( SlitMeshRefinedSystemTest );

  CPPUNIT_TEST( testMesh );

  CPPUNIT_TEST( testSystem );

#ifdef LIBMESH_ENABLE_UNIQUE_ID
  CPPUNIT_TEST( testRestart );
#endif

  CPPUNIT_TEST_SUITE_END();

protected:

  System* _sys;
  EquationSystems* _es;

public:

  void setUp()
  {
    this->build_mesh();

    // libMesh *should* renumber now, or a DistributedMesh might not
    // have contiguous ids, which is a requirement to write xda files.
    _mesh->allow_renumbering(true);

    _es = new EquationSystems(*_mesh);
    _sys = &_es->add_system<System> ("SimpleSystem");
    _sys->add_variable("u", FIRST);

    _es->init();
    SlitFunc slitfunc;
    _sys->project_solution(&slitfunc);

#ifdef LIBMESH_ENABLE_AMR
    MeshRefinement(*_mesh).uniformly_refine(1);
    _es->reinit();
    MeshRefinement(*_mesh).uniformly_refine(1);
    _es->reinit();
#endif
  }

  void tearDown()
  {
    delete _es;
    // _sys is owned by _es
    delete _mesh;
  }

  void testMesh()
  {
#ifdef LIBMESH_ENABLE_AMR
    // We should have 84 total and 64 active elements.
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)(4+16+64), _mesh->n_elem() );
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)64, _mesh->n_active_elem() );

    // We should have 88 nodes
    CPPUNIT_ASSERT_EQUAL( (dof_id_type)88, _mesh->n_nodes() );
#endif
  }

  void testSystem()
  {
    SlitFunc slitfunc;

    unsigned int dim = 2;

    CPPUNIT_ASSERT_EQUAL( _sys->n_vars(), 1u );

    FEMContext context(*_sys);
    FEBase * fe = NULL;
    context.get_element_fe( 0, fe, dim );
    const std::vector<Point> & xyz = fe->get_xyz();
    fe->get_phi();

    MeshBase::const_element_iterator       el     =
      _mesh->active_local_elements_begin();
    const MeshBase::const_element_iterator end_el =
      _mesh->active_local_elements_end();

    for (; el != end_el; ++el)
      {
        const Elem * elem = *el;
        context.pre_fe_reinit(*_sys, elem);
        context.elem_fe_reinit();

        const unsigned int n_qp = xyz.size();

        for (unsigned int qp=0; qp != n_qp; ++qp)
          {
            const Number exact_val = slitfunc(context, xyz[qp]);

            const Number discrete_val = context.interior_value(0, qp);

            CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(exact_val),
                                         libmesh_real(discrete_val),
                                         TOLERANCE*TOLERANCE);
          }
      }
  }

  void testRestart()
  {
    SlitFunc slitfunc;

    _mesh->write("slit_mesh.xda");
    _es->write("slit_solution.xda",
               EquationSystems::WRITE_DATA |
               EquationSystems::WRITE_SERIAL_FILES);

    Mesh mesh2(*TestCommWorld);
    mesh2.read("slit_mesh.xda");
    EquationSystems es2(mesh2);
    es2.read("slit_solution.xda");

    System & sys2 = es2.get_system<System> ("SimpleSystem");

    unsigned int dim = 2;

    CPPUNIT_ASSERT_EQUAL( sys2.n_vars(), 1u );

    FEMContext context(sys2);
    FEBase * fe = NULL;
    context.get_element_fe( 0, fe, dim );
    const std::vector<Point> & xyz = fe->get_xyz();
    fe->get_phi();

    // While we're in the middle of a unique id based test case, let's
    // make sure our unique ids were all read in correctly too.
    UniquePtr<PointLocatorBase> locator = _mesh->sub_point_locator();

    if (!_mesh->is_serial())
      locator->enable_out_of_mesh_mode();

    MeshBase::const_element_iterator       el     =
      mesh2.active_local_elements_begin();
    const MeshBase::const_element_iterator end_el =
      mesh2.active_local_elements_end();

    for (; el != end_el; ++el)
      {
        const Elem * elem = *el;

        const Elem * mesh1_elem = (*locator)(elem->centroid());
        if (mesh1_elem)
          {
            CPPUNIT_ASSERT_EQUAL( elem->unique_id(),
                                  mesh1_elem->unique_id() );

            for (unsigned int n=0; n != elem->n_nodes(); ++n)
              {
                const Node & node       = elem->node_ref(n);
                const Node & mesh1_node = mesh1_elem->node_ref(n);
                CPPUNIT_ASSERT_EQUAL( node.unique_id(),
                                      mesh1_node.unique_id() );
              }
          }

        context.pre_fe_reinit(sys2, elem);
        context.elem_fe_reinit();

        const unsigned int n_qp = xyz.size();

        for (unsigned int qp=0; qp != n_qp; ++qp)
          {
            const Number exact_val = slitfunc(context, xyz[qp]);

            const Number discrete_val = context.interior_value(0, qp);

            CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(exact_val),
                                         libmesh_real(discrete_val),
                                         TOLERANCE*TOLERANCE);
          }
      }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( SlitMeshTest );
CPPUNIT_TEST_SUITE_REGISTRATION( SlitMeshRefinedMeshTest );
CPPUNIT_TEST_SUITE_REGISTRATION( SlitMeshRefinedSystemTest );
