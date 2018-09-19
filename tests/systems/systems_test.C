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
#include <libmesh/zero_function.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/node_elem.h>
#include "libmesh/edge_edge2.h"
#include "libmesh/edge_edge3.h"

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

// Sparsity pattern augmentation class used in testDofCouplingWithVarGroups
class AugmentSparsityOnNodes : public GhostingFunctor
{
private:

  /**
   * The Mesh we're calculating on
   */
  MeshBase & _mesh;

public:

  /**
   * Constructor.
   */
  AugmentSparsityOnNodes(MeshBase & mesh)
  :
  _mesh(mesh)
  {}

  /**
   * User-defined function to augment the sparsity pattern.
   */
  virtual void operator() (const MeshBase::const_element_iterator & range_begin,
                           const MeshBase::const_element_iterator & range_end,
                           processor_id_type p,
                           map_type & coupled_elements) override
  {
    dof_id_type node_elem_id_1 = 2;
    dof_id_type node_elem_id_2 = 3;

    const CouplingMatrix * const null_mat = nullptr;

    for (const auto & elem : as_range(range_begin, range_end))
      {
        if (elem->id() == node_elem_id_1)
          {
            if (elem->processor_id() != p)
              {
                coupled_elements.insert (std::make_pair(elem, null_mat));

                const Elem * neighbor = _mesh.elem_ptr(node_elem_id_2);
                if (neighbor->processor_id() != p)
                  coupled_elements.insert (std::make_pair(neighbor, null_mat));
              }
          }
        if (elem->id() == node_elem_id_2)
          {
            if (elem->processor_id() != p)
              {
                coupled_elements.insert (std::make_pair(elem, null_mat));

                const Elem * neighbor = _mesh.elem_ptr(node_elem_id_1);
                if (neighbor->processor_id() != p)
                  coupled_elements.insert (std::make_pair(neighbor, null_mat));
              }
          }
      }
  }

  /**
   * Rebuild the cached _lower_to_upper map whenever our Mesh has
   * changed.
   */
  virtual void mesh_reinit () override
  {
  }

  /**
   * Update the cached _lower_to_upper map whenever our Mesh has been
   * redistributed.  We'll be lazy and just recalculate from scratch.
   */
  virtual void redistribute () override
  { this->mesh_reinit(); }

};

// Assembly function used in testDofCouplingWithVarGroups
void assemble_matrix_and_rhs(EquationSystems& es,
                             const std::string& system_name)
{
  const MeshBase& mesh = es.get_mesh();
  LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("test");
  const DofMap& dof_map = system.get_dof_map();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  std::vector<dof_id_type> dof_indices;

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      if(elem->type() == NODEELEM)
      {
        continue;
      }

      dof_map.dof_indices (elem, dof_indices);
      const unsigned int n_dofs = dof_indices.size();

      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      for(unsigned int i=0; i<n_dofs; i++)
      {
        Ke(i,i) = 1.;
        Fe(i) = 1.;
      }

      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }

  // Add matrix for extra coupled dofs
  {
    const Node & node_1 = mesh.node_ref(1);
    const Node & node_2 = mesh.node_ref(2);
    dof_indices.resize(6);
    dof_indices[0] =
      node_1.dof_number(system.number(), system.variable_number("u"), 0);
    dof_indices[1] =
      node_1.dof_number(system.number(), system.variable_number("v"), 0);
    dof_indices[2] =
      node_1.dof_number(system.number(), system.variable_number("w"), 0);

    dof_indices[3] =
      node_2.dof_number(system.number(), system.variable_number("u"), 0);
    dof_indices[4] =
      node_2.dof_number(system.number(), system.variable_number("v"), 0);
    dof_indices[5] =
      node_2.dof_number(system.number(), system.variable_number("w"), 0);

    const unsigned int n_dofs = dof_indices.size();
    Ke.resize (n_dofs, n_dofs);
    Fe.resize (n_dofs);

    for(unsigned int i=0; i<n_dofs; i++)
    {
      Ke(i,i) = 1.;
      Fe(i) = 1.;
    }

    system.matrix->add_matrix (Ke, dof_indices);
    system.rhs->add_vector    (Fe, dof_indices);
  }

  system.rhs->close();
  system.matrix->close();
}

// Assembly function used in testSubdomainVariableOrder
void assemble_matrix_and_rhs_variable_order(EquationSystems& es,
                                            const std::string& system_name)
{
  const MeshBase& mesh = es.get_mesh();
  LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("test");
  const DofMap& dof_map = system.get_dof_map();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  std::vector<dof_id_type> dof_indices;

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      if(elem->type() == NODEELEM)
      {
        continue;
      }

      dof_map.dof_indices (elem, dof_indices);
      const unsigned int n_dofs = dof_indices.size();

      Ke.resize (n_dofs, n_dofs);
      Fe.resize (n_dofs);

      for(unsigned int i=0; i<n_dofs; i++)
      {
        Ke(i,i) = 1.;
        Fe(i) = 1.;
      }

      system.matrix->add_matrix (Ke, dof_indices);
      system.rhs->add_vector    (Fe, dof_indices);
    }

  system.rhs->close();
  system.matrix->close();
}

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
  CPPUNIT_TEST( testBoundaryProjectCube );
  CPPUNIT_TEST( testDofCouplingWithVarGroups );
  CPPUNIT_TEST( testSubdomainVariableOrder );

#ifdef LIBMESH_ENABLE_AMR
#ifdef LIBMESH_HAVE_METAPHYSICL
#ifdef LIBMESH_HAVE_PETSC
  CPPUNIT_TEST( testProjectMatrixEdge2 );
  CPPUNIT_TEST( testProjectMatrixQuad4 );
  CPPUNIT_TEST( testProjectMatrixTri3 );
  CPPUNIT_TEST( testProjectMatrixHex8 );
  CPPUNIT_TEST( testProjectMatrixTet4 );
#endif // LIBMESH_HAVE_PETSC
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
    sys.project_solution(cubic_test, nullptr, es.parameters);

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
    sys.project_solution(cubic_test, nullptr, es.parameters);

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
    sys.project_solution(cubic_test, nullptr, es.parameters);

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
    sys.project_solution(cubic_test, nullptr, es.parameters);

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

  void testBoundaryProjectCube()
  {
    Mesh mesh(*TestCommWorld);

    // const boundary_id_type BOUNDARY_ID_MIN_Z  = 0;
    const boundary_id_type BOUNDARY_ID_MIN_Y  = 1;
    const boundary_id_type BOUNDARY_ID_MAX_X  = 2;
    const boundary_id_type BOUNDARY_ID_MAX_Y  = 3;
    const boundary_id_type BOUNDARY_ID_MIN_X  = 4;
    const boundary_id_type BOUNDARY_ID_MAX_Z  = 5;
    const boundary_id_type NODE_BOUNDARY_ID  = 10;
    const boundary_id_type EDGE_BOUNDARY_ID  = 20;
    const boundary_id_type SIDE_BOUNDARY_ID  = BOUNDARY_ID_MIN_X;

    EquationSystems es(mesh);
    System &sys = es.add_system<System> ("SimpleSystem");
    unsigned int u_var = sys.add_variable("u", FIRST, LAGRANGE);

    MeshTools::Generation::build_cube (mesh,
                                       3, 3, 3,
                                       0., 1., 0., 1., 0., 1.,
                                       HEX8);

    // Count the number of nodes on SIDE_BOUNDARY_ID
    std::set<dof_id_type> projected_nodes_set;

    for (const auto & elem : mesh.element_ptr_range())
      {
        for (auto side : elem->side_index_range())
          {
            std::vector<boundary_id_type> vec_to_fill;
            mesh.get_boundary_info().boundary_ids(elem, side, vec_to_fill);

            auto vec_it = std::find(vec_to_fill.begin(), vec_to_fill.end(), SIDE_BOUNDARY_ID);
            if (vec_it != vec_to_fill.end())
              {
                for (unsigned int node_index=0; node_index<elem->n_nodes(); node_index++)
                  {
                    if( elem->is_node_on_side(node_index, side))
                      {
                        projected_nodes_set.insert(elem->node_id(node_index));
                      }
                  }
              }
          }
      }

  // Also add some edge and node boundary IDs
  for (const auto & elem : mesh.element_ptr_range())
    {
      unsigned int
        side_max_x = 0, side_min_y = 0,
        side_max_y = 0, side_max_z = 0;

      bool
        found_side_max_x = false, found_side_max_y = false,
        found_side_min_y = false, found_side_max_z = false;

      for (auto side : elem->side_index_range())
        {
          if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MAX_X))
            {
              side_max_x = side;
              found_side_max_x = true;
            }

          if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MIN_Y))
            {
              side_min_y = side;
              found_side_min_y = true;
            }

          if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MAX_Y))
            {
              side_max_y = side;
              found_side_max_y = true;
            }

          if (mesh.get_boundary_info().has_boundary_id(elem, side, BOUNDARY_ID_MAX_Z))
            {
              side_max_z = side;
              found_side_max_z = true;
            }
        }

      // If elem has sides on boundaries
      // BOUNDARY_ID_MAX_X, BOUNDARY_ID_MAX_Y, BOUNDARY_ID_MAX_Z
      // then let's set a node boundary condition
      if (found_side_max_x && found_side_max_y && found_side_max_z)
        for (auto n : elem->node_index_range())
          if (elem->is_node_on_side(n, side_max_x) &&
              elem->is_node_on_side(n, side_max_y) &&
              elem->is_node_on_side(n, side_max_z))
            {
              projected_nodes_set.insert(elem->node_id(n));
              mesh.get_boundary_info().add_node(elem->node_ptr(n), NODE_BOUNDARY_ID);
            }

      // If elem has sides on boundaries
      // BOUNDARY_ID_MAX_X and BOUNDARY_ID_MIN_Y
      // then let's set an edge boundary condition
      if (found_side_max_x && found_side_min_y)
        for (auto e : elem->edge_index_range())
          if (elem->is_edge_on_side(e, side_max_x) &&
              elem->is_edge_on_side(e, side_min_y))
            {
              mesh.get_boundary_info().add_edge(elem, e, EDGE_BOUNDARY_ID);

              for (unsigned int node_index=0; node_index<elem->n_nodes(); node_index++)
                {
                  if (elem->is_node_on_edge(node_index, e))
                    {
                      projected_nodes_set.insert(elem->node_id(node_index));
                    }
                }
            }
    }

    es.init();

    sys.solution->add(1.0);

    std::set<boundary_id_type> boundary_ids;
    boundary_ids.insert(NODE_BOUNDARY_ID);
    boundary_ids.insert(EDGE_BOUNDARY_ID);
    boundary_ids.insert(SIDE_BOUNDARY_ID);
    std::vector<unsigned int> variables;
    variables.push_back(u_var);
    ZeroFunction<> zf;
    sys.boundary_project_solution(boundary_ids, variables, &zf);

    // On a distributed mesh we may not have every node on every
    // processor
    TestCommWorld->set_union(projected_nodes_set);

    // We set the solution to be 1 everywhere and then zero on specific boundary
    // nodes, so the final l1 norm of the solution is the difference between the
    // number of nodes in the mesh and the number of nodes we zero on.
    Real ref_l1_norm = static_cast<Real>(mesh.n_nodes()) - static_cast<Real>(projected_nodes_set.size());

    CPPUNIT_ASSERT_DOUBLES_EQUAL(sys.solution->l1_norm(), ref_l1_norm, TOLERANCE*TOLERANCE);
  }

  void testDofCouplingWithVarGroups()
  {
    ReplicatedMesh mesh(*TestCommWorld);

    MeshTools::Generation::build_cube (mesh,
                                      1,
                                      0,
                                      0,
                                      0., 1.,
                                      0., 0.,
                                      0., 0.,
                                      EDGE2);

    Point new_point_a(2.);
    Point new_point_b(3.);
    Node* new_node_a = mesh.add_point( new_point_a );
    Node* new_node_b = mesh.add_point( new_point_b );
    Elem* new_edge_elem = mesh.add_elem (new Edge2);
    new_edge_elem->set_node(0) = new_node_a;
    new_edge_elem->set_node(1) = new_node_b;

    mesh.elem_ref(0).subdomain_id() = 10;
    mesh.elem_ref(1).subdomain_id() = 10;

    // Add NodeElems for coupling purposes
    Elem* node_elem_1 = mesh.add_elem (new NodeElem);
    node_elem_1->set_node(0) = mesh.elem_ref(0).node_ptr(1);
    Elem* node_elem_2 = mesh.add_elem (new NodeElem);
    node_elem_2->set_node(0) = new_node_a;

    mesh.prepare_for_use();

    // Create an equation systems object.
    EquationSystems equation_systems (mesh);
    ExplicitSystem& system =
      equation_systems.add_system<LinearImplicitSystem> ("test");

    system.add_variable ("u", libMesh::FIRST);
    system.add_variable ("v", libMesh::FIRST);
    system.add_variable ("w", libMesh::FIRST);

    std::set<subdomain_id_type> theta_subdomains;
    theta_subdomains.insert(10);
    system.add_variable ("theta_x", libMesh::FIRST, &theta_subdomains);
    system.add_variable ("theta_y", libMesh::FIRST, &theta_subdomains);
    system.add_variable ("theta_z", libMesh::FIRST, &theta_subdomains);

    system.attach_assemble_function (assemble_matrix_and_rhs);

    AugmentSparsityOnNodes augment_sparsity(mesh);
    system.get_dof_map().add_coupling_functor(augment_sparsity);

    equation_systems.init ();

    system.solve();

    // We set the solution to be 1 everywhere, so the final l1 norm of the
    // solution is the product of the number of variables and number of nodes.
    Real ref_l1_norm = static_cast<Real>(mesh.n_nodes() * system.n_vars());

    CPPUNIT_ASSERT_DOUBLES_EQUAL(system.solution->l1_norm(), ref_l1_norm, TOLERANCE*TOLERANCE);
  }

  void testSubdomainVariableOrder()
  {
    ReplicatedMesh mesh(*TestCommWorld);

    MeshTools::Generation::build_cube (mesh,
                                       3,
                                       0,
                                       0,
                                       0., 1.,
                                       0., 0.,
                                       0., 0.,
                                       EDGE3);

    mesh.elem_ref(0).subdomain_id() = 10;
    mesh.elem_ref(1).subdomain_id() = 10;
    mesh.elem_ref(2).subdomain_id() = 20;

    // Create an equation systems object.
    EquationSystems equation_systems (mesh);
    ExplicitSystem& system =
      equation_systems.add_system<LinearImplicitSystem> ("test");

    system.add_variable ("u", libMesh::FIRST);
    system.add_variable ("v", libMesh::FIRST);
    system.add_variable ("w", libMesh::FIRST);

    std::map<subdomain_id_type, unsigned char> theta_subdomain_orders;
    theta_subdomain_orders[10] = 1; // Increase the order by 1 on subdomain 10
    system.add_variable ("theta_x", libMesh::FIRST, nullptr, &theta_subdomain_orders);
    system.add_variable ("theta_y", libMesh::FIRST, nullptr, &theta_subdomain_orders);
    system.add_variable ("theta_z", libMesh::FIRST, nullptr, &theta_subdomain_orders);

    system.attach_assemble_function (assemble_matrix_and_rhs_variable_order);

    equation_systems.init ();

    system.solve ();

    // We set the solution to be 1 everywhere, so the final l1 norm of the
    // solution is equal to the number of dofs in the system.
    Real ref_l1_norm = static_cast<Real>(system.n_dofs());

    CPPUNIT_ASSERT_DOUBLES_EQUAL(system.solution->l1_norm(), ref_l1_norm, TOLERANCE*TOLERANCE);
  }

#ifdef LIBMESH_ENABLE_AMR
#ifdef LIBMESH_HAVE_METAPHYSICL
#ifdef LIBMESH_HAVE_PETSC
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
    gold_mat.init(n_new_dofs, n_old_dofs, n_new_dofs_local, n_old_dofs_local);

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
    gold_mat.init(n_new_dofs, n_old_dofs, n_new_dofs_local, n_old_dofs_local);

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
    gold_mat.init(n_new_dofs, n_old_dofs, n_new_dofs_local, n_old_dofs_local);

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
#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_HAVE_METAPHYSICL
#endif // LIBMESH_ENABLE_AMR


  void testProjectHierarchicEdge3() { testProjectLine(EDGE3); }
  void testProjectHierarchicQuad9() { testProjectSquare(QUAD9); }
  void testProjectHierarchicTri6()  { testProjectSquare(TRI6); }
  void testProjectHierarchicHex27() { testProjectCube(HEX27); }
  void testProjectMeshFunctionHex27() { testProjectCubeWithMeshFunction(HEX27); }

#ifdef LIBMESH_ENABLE_AMR
#ifdef LIBMESH_HAVE_METAPHYSICL
#ifdef LIBMESH_HAVE_PETSC
  // projection matrix tests
  void testProjectMatrixEdge2() { testProjectMatrix1D(EDGE2); }
  void testProjectMatrixQuad4() { testProjectMatrix2D(QUAD4); }
  void testProjectMatrixTri3() { testProjectMatrix2D(TRI3); }
  void testProjectMatrixHex8() { testProjectMatrix3D(HEX8); }
  void testProjectMatrixTet4() { testProjectMatrix3D(TET4); }
#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_HAVE_METAPHYSICL
#endif // LIBMESH_ENABLE_AMR

};

CPPUNIT_TEST_SUITE_REGISTRATION( SystemsTest );
