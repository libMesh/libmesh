// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

#include <libmesh/equation_systems.h>
#include <libmesh/mesh.h>
#include <libmesh/node.h>
#include <libmesh/dof_map.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/mesh_function.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/sparse_matrix.h>
#include "libmesh/string_to_enum.h"
#include <libmesh/cell_tet4.h>

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

Number cubic_test (const Point& p,
                   const Parameters&,
                   const std::string&,
                   const std::string&)
{
  const Real & x = p(0);
  const Real & y = p(1);
  const Real & z = p(2);

  return x*(1-x)*(1-x) + x*x*(1-y) + x*(1-y)*(1-z) + y*(1-y)*z + z*(1-z)*(1-z);
}


class SystemsTest : public CppUnit::TestCase {
public:
  CPPUNIT_TEST_SUITE( SystemsTest );

  CPPUNIT_TEST( testProjectHierarchicEdge3 );
  CPPUNIT_TEST( testProjectHierarchicQuad9 );
  CPPUNIT_TEST( testProjectHierarchicTri6 );
  CPPUNIT_TEST( testProjectHierarchicHex27 );
  CPPUNIT_TEST( testProjectMeshFunctionHex27 );

#ifdef LIBMESH_ENABLE_AMR
#ifdef LIBMESH_HAVE_METAPHYSICL
  CPPUNIT_TEST( testProjectMatrixEdge2 );
  CPPUNIT_TEST( testProjectMatrixQuad4 );
  CPPUNIT_TEST( testProjectMatrixTri3 );
  CPPUNIT_TEST( testProjectMatrixHex8 );
  CPPUNIT_TEST( testProjectMatrixTet4 );
#endif // LIBMESH_HAVE_METAPHYSICL
#endif // LIBMESH_ENABLE_AMR

  CPPUNIT_TEST_SUITE_END();

private:

public:
  void setUp()
  {}

  void tearDown()
  {}

  void testProjectLine(const ElemType elem_type)
  {
    Mesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    System &sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", THIRD, HIERARCHIC);

    MeshTools::Generation::build_line (mesh,
                                       3,
                                       0., 1.,
                                       elem_type);

    es.init();
    sys.project_solution(cubic_test, NULL, es.parameters);

    for (Real x = 0.1; x < 1; x += 0.2)
      {
        Point p(x);
        CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys.point_value(0,p)),
                                     libmesh_real(cubic_test(p,es.parameters,"","")),
                                     TOLERANCE*TOLERANCE);
      }
  }

  void testProjectSquare(const ElemType elem_type)
  {
    Mesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    System &sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", THIRD, HIERARCHIC);

    MeshTools::Generation::build_square (mesh,
                                         3, 3,
                                         0., 1., 0., 1.,
                                         elem_type);

    es.init();
    sys.project_solution(cubic_test, NULL, es.parameters);

    for (Real x = 0.1; x < 1; x += 0.2)
      for (Real y = 0.1; y < 1; y += 0.2)
        {
          Point p(x,y);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys.point_value(0,p)),
                                       libmesh_real(cubic_test(p,es.parameters,"","")),
                                       TOLERANCE*TOLERANCE);
        }
  }

  void testProjectCube(const ElemType elem_type)
  {
    Mesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    System &sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", THIRD, HIERARCHIC);

    MeshTools::Generation::build_cube (mesh,
                                       3, 3, 3,
                                       0., 1., 0., 1., 0., 1.,
                                       elem_type);

    es.init();
    sys.project_solution(cubic_test, NULL, es.parameters);

    for (Real x = 0.1; x < 1; x += 0.2)
      for (Real y = 0.1; y < 1; y += 0.2)
        for (Real z = 0.1; z < 1; z += 0.2)
          {
            Point p(x,y,z);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(sys.point_value(0,p)),
                                         libmesh_real(cubic_test(p,es.parameters,"","")),
                                         TOLERANCE*TOLERANCE);
          }
  }

  void testProjectCubeWithMeshFunction(const ElemType elem_type)
  {
    // The source mesh needs to exist everywhere it's queried, so we
    // use a ReplicatedMesh
    ReplicatedMesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    System &sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", THIRD, HIERARCHIC);

    MeshTools::Generation::build_cube (mesh,
                                       3, 3, 3,
                                       0., 1., 0., 1., 0., 1.,
                                       elem_type);

    es.init();
    sys.project_solution(cubic_test, NULL, es.parameters);

    std::vector<unsigned int> variables;
    sys.get_all_variable_numbers(variables);
    std::sort(variables.begin(),variables.end());

    std::unique_ptr< NumericVector<Number> > mesh_function_vector =
      NumericVector<Number>::build(es.comm());
    mesh_function_vector->init(sys.n_dofs(), false, SERIAL);
    sys.solution->localize( *mesh_function_vector );

    MeshFunction mesh_function(es,
                               *mesh_function_vector,
                               sys.get_dof_map(),
                               variables);
    mesh_function.init();

    // Make a second system and project onto it using a MeshFunction
    Mesh proj_mesh(*TestCommWorld);
    EquationSystems proj_es(proj_mesh);

    System &proj_sys = proj_es.add_system<System> ("ProjectionSystem");
    proj_sys.add_variable("u", SECOND, LAGRANGE);

    MeshTools::Generation::build_cube (proj_mesh,
                                       5, 5, 5,
                                       0., 1., 0., 1., 0., 1.,
                                       elem_type);

    proj_es.init();
    proj_sys.project_solution(&mesh_function);

    for (Real x = 0.1; x < 1; x += 0.2)
      for (Real y = 0.1; y < 1; y += 0.2)
        for (Real z = 0.1; z < 1; z += 0.2)
          {
            Point p(x,y,z);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(libmesh_real(proj_sys.point_value(0,p)),
                                         libmesh_real(cubic_test(p,es.parameters,"","")),
                                         TOLERANCE*TOLERANCE);
          }
  }

#ifdef LIBMESH_ENABLE_AMR
#ifdef LIBMESH_HAVE_METAPHYSICL
  void testProjectMatrix1D(const ElemType elem_type)
  {
    // Use ReplicatedMesh to get consistent child element node
    // numbering during refinement
    ReplicatedMesh mesh(*TestCommWorld);

    // fix the node numbering to resolve dof_id numbering issues in parallel tests
    mesh.allow_renumbering(false);

    // init a simple 1d system
    EquationSystems es(mesh);
    System &sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", FIRST, LAGRANGE);

    MeshTools::Generation::build_line (mesh,
                                       4, 0., 1.,
                                       elem_type);

    es.init();

    // static set of coarse nodes / order of fine grid nodes from x=0 to x=1 going left to right
    std::set<dof_id_type> coarse_nodes({0,1,2,3,4});
    std::vector<dof_id_type> node_order_f({0,5,1,6,2,7,3,8,4});

    // stash number of dofs on coarse grid for projection sizing
    int n_old_dofs = sys.n_dofs();

    // save old coarse dof_ids in order of coarse nodes
    std::map <dof_id_type, dof_id_type> node2dof_c;
    for ( const auto & node : mesh.node_ptr_range() )
      {
        dof_id_type cdof_id = node->dof_number(0,0,0);
        node2dof_c.insert( std::pair<dof_id_type,dof_id_type>( node->id() , cdof_id) );
      }

    // refine the mesh so we can utilize old_dof_indices for projection_matrix
    MeshRefinement mr(mesh);
    mr.uniformly_refine(1);
    sys.get_dof_map().distribute_dofs(mesh);

    // fine node to dof map
    std::map <dof_id_type, dof_id_type> node2dof_f;
    for ( const auto & node : mesh.local_node_ptr_range() )
      {
        dof_id_type fdof_id = node->dof_number(0,0,0);
        node2dof_f.insert( std::pair<dof_id_type,dof_id_type>(node->id() , fdof_id) );
      }

    // local and global projection_matrix sizes infos
    int n_new_dofs = sys.n_dofs();
    int n_new_dofs_local = sys.get_dof_map().n_dofs_on_processor(sys.processor_id());
    int ndofs_old_first = sys.get_dof_map().first_old_dof(sys.processor_id());
    int ndofs_old_end   = sys.get_dof_map().end_old_dof(sys.processor_id());
    int n_old_dofs_local = ndofs_old_end - ndofs_old_first;

    // init and compute the projection matrix using GenericProjector
    std::unique_ptr<SparseMatrix<Number> > proj_mat_ptr =
      SparseMatrix<Number>::build(*TestCommWorld);
    SparseMatrix<Number> & proj_mat = *proj_mat_ptr;
    proj_mat.init(n_new_dofs, n_old_dofs, n_new_dofs_local, n_old_dofs_local);
    sys.projection_matrix(proj_mat);
    proj_mat.close();

    // init the gold standard projection matrix
    std::unique_ptr<SparseMatrix<Number> > gold_mat_ptr =
      SparseMatrix<Number>::build(*TestCommWorld);
    SparseMatrix<Number> & gold_mat = *gold_mat_ptr;
    gold_mat.init(n_new_dofs, n_old_dofs, n_new_dofs_local, n_old_dofs);

    // construct the gold projection matrix using static node numbering as reference info
    for ( const auto & node : mesh.local_node_ptr_range() )
      {
        dof_id_type node_id = node->id();
        dof_id_type fdof_id = (node2dof_f.find(node_id))->second;

        if (coarse_nodes.find(node_id) != coarse_nodes.end() )
          { //direct inject coarse nodes
            if (fdof_id >= gold_mat.row_start() && fdof_id < gold_mat.row_stop())
              {
                auto cdof_id = node2dof_c.find(node_id);
                gold_mat.set(fdof_id, cdof_id->second, 1.0);
              }
          }
        else
          { // new nodes with old_dof neighbor contributions
            if (fdof_id >= gold_mat.row_start() && fdof_id < gold_mat.row_stop())
              {
                auto node_loc = std::find(node_order_f.begin(), node_order_f.end(), node_id);
                auto node_n = *std::next(node_loc, 1);
                auto node_p = *std::prev(node_loc, 1);
                auto dof_p = node2dof_c.find(node_p);
                auto dof_n = node2dof_c.find(node_n);

                gold_mat.set(fdof_id, dof_p->second, 0.5);
                gold_mat.set(fdof_id, dof_n->second, 0.5);
              }
          }
      } // end gold mat build
    gold_mat.close();

    // calculate relative difference norm between the two projection matrices
    Real gold_norm = gold_mat.linfty_norm();
    gold_mat.add(-1.0, proj_mat);
    Real diff_norm = gold_mat.linfty_norm();
    CPPUNIT_ASSERT(diff_norm/gold_norm < TOLERANCE*TOLERANCE);
  }

  void testProjectMatrix2D(const ElemType elem_type)
  {
    // Use ReplicatedMesh to get consistent child element node
    // numbering during refinement
    ReplicatedMesh mesh(*TestCommWorld);

    // fix the node numbering to resolve dof_id numbering issues in parallel tests
    mesh.allow_renumbering(false);

    // init a simple 1d system
    EquationSystems es(mesh);
    System &sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", FIRST, LAGRANGE);

    if (elem_type == Utility::string_to_enum<ElemType>("QUAD4"))
      MeshTools::Generation::build_square (mesh,
                                           2, 2,
                                           0., 1., 0., 1.,
                                           elem_type);
    else if (elem_type == Utility::string_to_enum<ElemType>("TRI3"))
      MeshTools::Generation::build_square (mesh,
                                           1, 1,
                                           0., 1., 0., 1.,
                                           elem_type);

    es.init();

    // static sets of nodes and their neighbors
    std::set<dof_id_type> coarse_nodes;
    std::map<dof_id_type, std::vector<dof_id_type>> side_nbr_nodes;
    std::map<dof_id_type, std::vector<dof_id_type>> int_nbr_nodes;

    // fill neighbor maps based on static node numbering
    if (elem_type == Utility::string_to_enum<ElemType>("QUAD4"))
      {
        coarse_nodes.insert({0,1,2,3,4,5,6,7,8});

        side_nbr_nodes.insert({9, {0,1}});
        side_nbr_nodes.insert({14, {1,2}});
        side_nbr_nodes.insert({11, {0,3}});
        side_nbr_nodes.insert({12, {1,4}});
        side_nbr_nodes.insert({16, {2,5}});
        side_nbr_nodes.insert({13, {3,4}});
        side_nbr_nodes.insert({17, {4,5}});
        side_nbr_nodes.insert({19, {3,6}});
        side_nbr_nodes.insert({20, {4,7}});
        side_nbr_nodes.insert({23, {5,8}});
        side_nbr_nodes.insert({21, {6,7}});
        side_nbr_nodes.insert({24, {7,8}});

        int_nbr_nodes.insert({10, {0,1,3,4}});
        int_nbr_nodes.insert({15, {1,2,4,5}});
        int_nbr_nodes.insert({18, {3,4,6,7}});
        int_nbr_nodes.insert({22, {4,5,7,8}});
      }
    else if (elem_type == Utility::string_to_enum<ElemType>("TRI3"))
      {
        coarse_nodes.insert({0,1,2,3});

        side_nbr_nodes.insert({4, {0,1}});
        side_nbr_nodes.insert({5, {0,3}});
        side_nbr_nodes.insert({6, {1,3}});
        side_nbr_nodes.insert({7, {0,2}});
        side_nbr_nodes.insert({8, {2,3}});
      }

    // stash number of dofs on coarse grid for projection sizing
    int n_old_dofs = sys.n_dofs();

    // save old coarse dof_ids in order of coarse nodes
    std::map <dof_id_type, dof_id_type> node2dof_c;
    for ( const auto & node : mesh.node_ptr_range() )
      {
        dof_id_type cdof_id = node->dof_number(0,0,0);
        node2dof_c.insert( std::pair<dof_id_type,dof_id_type>( node->id() , cdof_id) );
      }

    // refine the mesh so we can utilize old_dof_indices for projection_matrix
    MeshRefinement mr(mesh);
    mr.uniformly_refine(1);
    sys.get_dof_map().distribute_dofs(mesh);

    // fine node to dof map
    std::map <dof_id_type, dof_id_type> node2dof_f;
    for ( const auto & node : mesh.local_node_ptr_range() )
      {
        dof_id_type fdof_id = node->dof_number(0,0,0);
        node2dof_f.insert( std::pair<dof_id_type,dof_id_type>(node->id() , fdof_id) );
      }

    // local and global projection_matrix sizes infos
    int n_new_dofs = sys.n_dofs();
    int n_new_dofs_local = sys.get_dof_map().n_dofs_on_processor(sys.processor_id());
    int ndofs_old_first = sys.get_dof_map().first_old_dof(sys.processor_id());
    int ndofs_old_end   = sys.get_dof_map().end_old_dof(sys.processor_id());
    int n_old_dofs_local = ndofs_old_end - ndofs_old_first;

    // init and compute the projection matrix using GenericProjector
    std::unique_ptr<SparseMatrix<Number> > proj_mat_ptr =
      SparseMatrix<Number>::build(*TestCommWorld);
    SparseMatrix<Number> & proj_mat = *proj_mat_ptr;
    proj_mat.init(n_new_dofs, n_old_dofs, n_new_dofs_local, n_old_dofs_local);
    sys.projection_matrix(proj_mat);
    proj_mat.close();

    // init the gold standard projection matrix
    std::unique_ptr<SparseMatrix<Number> > gold_mat_ptr =
      SparseMatrix<Number>::build(*TestCommWorld);
    SparseMatrix<Number> & gold_mat = *gold_mat_ptr;
    gold_mat.init(n_new_dofs, n_old_dofs, n_new_dofs_local, n_old_dofs);

    // construct the gold projection matrix using static node numbering as reference info
    for ( const auto & node : mesh.local_node_ptr_range() )
      {
        dof_id_type node_id = node->id();
        dof_id_type fdof_id = (node2dof_f.find(node_id))->second;

        if (coarse_nodes.find(node_id) != coarse_nodes.end() )
          { // direct inject coarse nodes
            if (fdof_id >= gold_mat.row_start() && fdof_id < gold_mat.row_stop())
              {
                auto cdof_id = node2dof_c.find(node_id);
                gold_mat.set(fdof_id, cdof_id->second, 1.0);
              }
          }
        else if ( side_nbr_nodes.find(node_id) != side_nbr_nodes.end() )
          { // new side nodes with old_dof neighbor contributions
            if (fdof_id >= gold_mat.row_start() && fdof_id < gold_mat.row_stop())
              {
                auto node_nbrs = side_nbr_nodes.find(node_id);
                for (auto nbr : node_nbrs->second)
                  {
                    auto nbr_dof = node2dof_c.find(nbr);
                    gold_mat.set(fdof_id, nbr_dof->second, 0.5);
                  }
              }
          }
        else
          { // new interior nodes with old_dof neighbor contributions
            if (fdof_id >= gold_mat.row_start() && fdof_id < gold_mat.row_stop())
              {
                auto node_nbrs = int_nbr_nodes.find(node_id);
                for (auto nbr : node_nbrs->second)
                  {
                    auto nbr_dof = node2dof_c.find(nbr);
                    gold_mat.set(fdof_id, nbr_dof->second, 0.25);
                  }
              }
          }
      } // end gold mat build
    gold_mat.close();

    // calculate relative difference norm between the two projection matrices
    Real gold_norm = gold_mat.linfty_norm();
    proj_mat.add(-1.0, gold_mat);
    Real diff_norm = proj_mat.linfty_norm();
    CPPUNIT_ASSERT(diff_norm/gold_norm < TOLERANCE*TOLERANCE);
  }

  void testProjectMatrix3D(const ElemType elem_type)
  {
    // Use ReplicatedMesh to get consistent child element node
    // numbering during refinement
    ReplicatedMesh mesh(*TestCommWorld);

    // fix the node numbering to resolve dof_id numbering issues in parallel tests
    mesh.allow_renumbering(false);

    // init a simple 1d system
    EquationSystems es(mesh);
    System &sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", FIRST, LAGRANGE);

    if (elem_type == Utility::string_to_enum<ElemType>("HEX8"))
      MeshTools::Generation::build_cube (mesh,
                                         1, 1, 1,
                                         0., 1., 0., 1., 0., 1.,
                                         elem_type);
    else if (elem_type == Utility::string_to_enum<ElemType>("TET4"))
      {
        // manually build a Tet4 element
        mesh.add_point( Point(0,0,0), 0 );
        mesh.add_point( Point(1,0,0), 1 );
        mesh.add_point( Point(0,1,0), 2 );
        mesh.add_point( Point(1./3.,1./3.,1), 3 );

        Elem * elem = new Tet4();
        elem->set_id(0);
        elem = mesh.add_elem(elem);
        elem->set_node(0) = mesh.node_ptr(0);
        elem->set_node(1) = mesh.node_ptr(1);
        elem->set_node(2) = mesh.node_ptr(2);
        elem->set_node(3) = mesh.node_ptr(3);

        mesh.prepare_for_use();
      }
    es.init();

    // static sets of nodes and their neighbors
    std::set<dof_id_type> coarse_nodes;
    std::map<dof_id_type, std::vector<dof_id_type>> side_nbr_nodes;
    std::map<dof_id_type, std::vector<dof_id_type>> face_nbr_nodes;
    std::map<dof_id_type, std::vector<dof_id_type>> int_nbr_nodes;

    if (elem_type == Utility::string_to_enum<ElemType>("HEX8"))
      {
        coarse_nodes.insert({0,1,2,3,4,5,6,7});

        // fill neighbor maps based on static node numbering
        side_nbr_nodes.insert({8, {0,1}});
        side_nbr_nodes.insert({10, {0,2}});
        side_nbr_nodes.insert({15, {1,3}});
        side_nbr_nodes.insert({18, {2,3}});
        side_nbr_nodes.insert({11, {0,4}});
        side_nbr_nodes.insert({16, {1,5}});
        side_nbr_nodes.insert({21, {3,7}});
        side_nbr_nodes.insert({20, {2,6}});
        side_nbr_nodes.insert({22, {4,5}});
        side_nbr_nodes.insert({24, {4,6}});
        side_nbr_nodes.insert({25, {5,7}});
        side_nbr_nodes.insert({26, {6,7}});

        face_nbr_nodes.insert({12, {0,1,4,5}});
        face_nbr_nodes.insert({9 , {0,1,2,3}});
        face_nbr_nodes.insert({14, {0,2,4,6}});
        face_nbr_nodes.insert({17, {1,3,5,7}});
        face_nbr_nodes.insert({19, {2,3,6,7}});
        face_nbr_nodes.insert({23, {4,5,6,7}});

        int_nbr_nodes.insert({13, {0,1,2,3,4,5,6,7}});
      }
    else if (elem_type == Utility::string_to_enum<ElemType>("TET4"))
      {
        coarse_nodes.insert({0,1,2,3});

        // fill neighbor maps based on static node numbering
        side_nbr_nodes.insert({4, {0,1}});
        side_nbr_nodes.insert({5, {0,2}});
        side_nbr_nodes.insert({6, {0,3}});
        side_nbr_nodes.insert({7, {1,2}});
        side_nbr_nodes.insert({8, {1,3}});
        side_nbr_nodes.insert({9, {2,3}});
      }

    // stash number of dofs on coarse grid for projection sizing
    int n_old_dofs = sys.n_dofs();

    // save old coarse dof_ids in order of coarse nodes
    std::map <dof_id_type, dof_id_type> node2dof_c;
    for ( const auto & node : mesh.node_ptr_range() )
      {
        dof_id_type cdof_id = node->dof_number(0,0,0);
        node2dof_c.insert( std::pair<dof_id_type,dof_id_type>( node->id() , cdof_id) );
      }

    // refine the mesh so we can utilize old_dof_indices for projection_matrix
    MeshRefinement mr(mesh);
    mr.uniformly_refine(1);
    sys.get_dof_map().distribute_dofs(mesh);

    // fine node to dof map
    std::map <dof_id_type, dof_id_type> node2dof_f;
    for ( const auto & node : mesh.local_node_ptr_range() )
      {
        dof_id_type fdof_id = node->dof_number(0,0,0);
        node2dof_f.insert( std::pair<dof_id_type,dof_id_type>(node->id() , fdof_id) );
      }

    // local and global projection_matrix sizes infos
    int n_new_dofs = sys.n_dofs();
    int n_new_dofs_local = sys.get_dof_map().n_dofs_on_processor(sys.processor_id());
    int ndofs_old_first = sys.get_dof_map().first_old_dof(sys.processor_id());
    int ndofs_old_end   = sys.get_dof_map().end_old_dof(sys.processor_id());
    int n_old_dofs_local = ndofs_old_end - ndofs_old_first;

    // init and compute the projection matrix using GenericProjector
    std::unique_ptr<SparseMatrix<Number> > proj_mat_ptr =
      SparseMatrix<Number>::build(*TestCommWorld);
    SparseMatrix<Number> & proj_mat = *proj_mat_ptr;
    proj_mat.init(n_new_dofs, n_old_dofs, n_new_dofs_local, n_old_dofs_local);
    sys.projection_matrix(proj_mat);
    proj_mat.close();

    // init the gold standard projection matrix
    std::unique_ptr<SparseMatrix<Number> > gold_mat_ptr =
      SparseMatrix<Number>::build(*TestCommWorld);
    SparseMatrix<Number> & gold_mat = *gold_mat_ptr;
    gold_mat.init(n_new_dofs, n_old_dofs, n_new_dofs_local, n_old_dofs);

    // construct the gold projection matrix using static node numbering as reference info
    for ( const auto & node : mesh.local_node_ptr_range() )
      {
        dof_id_type node_id = node->id();
        dof_id_type fdof_id = (node2dof_f.find(node_id))->second;

        if (coarse_nodes.find(node_id) != coarse_nodes.end() )
          { // direct inject coarse nodes
            if (fdof_id >= gold_mat.row_start() && fdof_id < gold_mat.row_stop())
              {
                auto cdof_id = node2dof_c.find(node_id);
                gold_mat.set(fdof_id, cdof_id->second, 1.0);
              }
          }
        else if ( side_nbr_nodes.find(node_id) != side_nbr_nodes.end() )
          { // new side nodes with old_dof neighbor contributions
            if (fdof_id >= gold_mat.row_start() && fdof_id < gold_mat.row_stop())
              {
                auto node_nbrs = side_nbr_nodes.find(node_id);
                for (auto nbr : node_nbrs->second)
                  {
                    auto nbr_dof = node2dof_c.find(nbr);
                    gold_mat.set(fdof_id, nbr_dof->second, 0.5);
                  }
              }
          }
        else if ( face_nbr_nodes.find(node_id) != face_nbr_nodes.end() )
          { // new face nodes with old_dof neighbor contributions
            if (fdof_id >= gold_mat.row_start() && fdof_id < gold_mat.row_stop())
              {
                auto node_nbrs = face_nbr_nodes.find(node_id);
                for (auto nbr : node_nbrs->second)
                  {
                    auto nbr_dof = node2dof_c.find(nbr);
                    gold_mat.set(fdof_id, nbr_dof->second, 0.25);
                  }
              }
          }
        else
          { // new interior nodes with old_dof neighbor contributions
            if (fdof_id >= gold_mat.row_start() && fdof_id < gold_mat.row_stop())
              {
                auto node_nbrs = int_nbr_nodes.find(node_id);
                for (auto nbr : node_nbrs->second)
                  {
                    auto nbr_dof = node2dof_c.find(nbr);
                    gold_mat.set(fdof_id, nbr_dof->second, 0.125);
                  }
              }
          }
      } // end gold mat build
    gold_mat.close();

    // calculate relative difference norm between the two projection matrices
    Real gold_norm = gold_mat.linfty_norm();
    proj_mat.add(-1.0, gold_mat);
    Real diff_norm = proj_mat.linfty_norm();
    CPPUNIT_ASSERT(diff_norm/gold_norm < TOLERANCE*TOLERANCE);
  }
#endif // LIBMESH_HAVE_METAPHYSICL
#endif // LIBMESH_ENABLE_AMR


  void testProjectHierarchicEdge3() { testProjectLine(EDGE3); }
  void testProjectHierarchicQuad9() { testProjectSquare(QUAD9); }
  void testProjectHierarchicTri6()  { testProjectSquare(TRI6); }
  void testProjectHierarchicHex27() { testProjectCube(HEX27); }
  void testProjectMeshFunctionHex27() { testProjectCubeWithMeshFunction(HEX27); }

#ifdef LIBMESH_ENABLE_AMR
#ifdef LIBMESH_HAVE_METAPHYSICL
  // projection matrix tests
  void testProjectMatrixEdge2() { testProjectMatrix1D(EDGE2); }
  void testProjectMatrixQuad4() { testProjectMatrix2D(QUAD4); }
  void testProjectMatrixTri3() { testProjectMatrix2D(TRI3); }
  void testProjectMatrixHex8() { testProjectMatrix3D(HEX8); }
  void testProjectMatrixTet4() { testProjectMatrix3D(TET4); }
#endif // LIBMESH_HAVE_METAPHYSICL
#endif // LIBMESH_ENABLE_AMR

};

CPPUNIT_TEST_SUITE_REGISTRATION( SystemsTest );
