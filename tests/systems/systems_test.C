#include <libmesh/equation_systems.h>
#include <libmesh/int_range.h>
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
#include <libmesh/transient_system.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/node_elem.h>
#include <libmesh/edge_edge2.h>
#include <libmesh/dg_fem_context.h>
#include <libmesh/enum_solver_type.h>
#include <libmesh/enum_preconditioner_type.h>
#include <libmesh/linear_solver.h>
#include <libmesh/parallel.h>
#include <libmesh/face_quad4.h>
#include <libmesh/face_quad9.h>
#include <libmesh/face_quad8.h>
#include <libmesh/face_tri3.h>
#include <libmesh/face_tri6.h>
#include <libmesh/face_tri7.h>
#include <libmesh/cell_hex8.h>
#include <libmesh/cell_hex20.h>
#include <libmesh/cell_hex27.h>
#include <libmesh/cell_tet10.h>
#include <libmesh/cell_tet14.h>
#include <libmesh/boundary_info.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <string>

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

  virtual std::unique_ptr<GhostingFunctor> clone () const override
  {
    return std::make_unique<AugmentSparsityOnNodes>(_mesh);
  }

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
                coupled_elements.emplace(elem, null_mat);

                const Elem * neighbor = _mesh.elem_ptr(node_elem_id_2);
                if (neighbor->processor_id() != p)
                  coupled_elements.emplace(neighbor, null_mat);
              }
          }
        if (elem->id() == node_elem_id_2)
          {
            if (elem->processor_id() != p)
              {
                coupled_elements.emplace(elem, null_mat);

                const Elem * neighbor = _mesh.elem_ptr(node_elem_id_1);
                if (neighbor->processor_id() != p)
                  coupled_elements.emplace(neighbor, null_mat);
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
                             const std::string&)
{
  const MeshBase& mesh = es.get_mesh();
  LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("test");
  const DofMap& dof_map = system.get_dof_map();

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  std::vector<dof_id_type> dof_indices;

  SparseMatrix<Number> & matrix = system.get_system_matrix();

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

      matrix.add_matrix (Ke, dof_indices);
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

    matrix.add_matrix (Ke, dof_indices);
    system.rhs->add_vector    (Fe, dof_indices);
  }

  system.rhs->close();
  matrix.close();
}

// Assembly function that uses a DGFEMContext
void assembly_with_dg_fem_context(EquationSystems& es,
                                  const std::string& /*system_name*/)
{
  const MeshBase& mesh = es.get_mesh();
  LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("test");

  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;

  std::vector<dof_id_type> dof_indices;
  SparseMatrix<Number> & matrix = system.get_system_matrix();

  DGFEMContext context(system);
  {
    // For efficiency, we should prerequest all
    // the data we will need to build the
    // linear system before doing an element loop.
    FEBase* elem_fe = NULL;
    context.get_element_fe(0, elem_fe);
    elem_fe->get_JxW();
    elem_fe->get_phi();
    elem_fe->get_dphi();

    FEBase* side_fe = NULL;
    context.get_side_fe( 0, side_fe );
    side_fe->get_JxW();
    side_fe->get_phi();

    FEBase* neighbor_side_fe = NULL;
    context.get_neighbor_side_fe(0, neighbor_side_fe);
    neighbor_side_fe->get_phi();
  }

  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      context.pre_fe_reinit(system, elem);
      context.elem_fe_reinit();

      // Element interior assembly
      {
        FEBase* elem_fe = NULL;
        context.get_element_fe(0, elem_fe);

        const std::vector<Real> &JxW = elem_fe->get_JxW();
        const std::vector<std::vector<Real> >& phi = elem_fe->get_phi();
        const std::vector<std::vector<RealGradient> >& dphi = elem_fe->get_dphi();

        unsigned int n_dofs = context.get_dof_indices(0).size();
        unsigned int n_qpoints = context.get_element_qrule().n_points();

        for (unsigned int qp=0; qp != n_qpoints; qp++)
          for (unsigned int i=0; i != n_dofs; i++)
            for (unsigned int j=0; j != n_dofs; j++)
              context.get_elem_jacobian()(i,j) += JxW[qp] * dphi[i][qp]*dphi[j][qp];

        for (unsigned int qp=0; qp != n_qpoints; qp++)
          for (unsigned int i=0; i != n_dofs; i++)
            context.get_elem_residual()(i) += JxW[qp] * phi[i][qp];
      }

      matrix.add_matrix (context.get_elem_jacobian(), context.get_dof_indices());
      system.rhs->add_vector (context.get_elem_residual(), context.get_dof_indices());

      // Element side assembly
      for (context.side = 0; context.side != elem->n_sides(); ++context.side)
        {
          // If there is a neighbor, then we proceed with assembly
          // that involves elem and neighbor
          const Elem* neighbor = elem->neighbor_ptr(context.get_side());
          if(neighbor)
            {
              context.side_fe_reinit();
              context.set_neighbor(*neighbor);

              // This call initializes neighbor data, and also sets
              // context.dg_terms_are_active() to true
              context.neighbor_side_fe_reinit();

              FEBase* side_fe = NULL;
              context.get_side_fe(0, side_fe);

              const std::vector<Real> &JxW_face = side_fe->get_JxW();
              const std::vector<std::vector<Real> >& phi_face = side_fe->get_phi();

              FEBase* neighbor_side_fe = NULL;
              context.get_neighbor_side_fe(0, neighbor_side_fe);

              // These shape functions have been evaluated on the quadrature points
              // for elem->side on the neighbor element
              const std::vector<std::vector<Real> >& phi_neighbor_face =
                neighbor_side_fe->get_phi();

              const unsigned int n_dofs = context.get_dof_indices(0).size();
              const unsigned int n_neighbor_dofs = context.get_neighbor_dof_indices(0).size();
              unsigned int n_sidepoints = context.get_side_qrule().n_points();

              for (unsigned int qp=0; qp<n_sidepoints; qp++)
                {
                  for (unsigned int i=0; i<n_dofs; i++)
                    for (unsigned int j=0; j<n_dofs; j++)
                      {
                        context.get_elem_elem_jacobian()(i,j) +=
                          JxW_face[qp] * phi_face[i][qp] * phi_face[j][qp];
                      }

                  for (unsigned int i=0; i<n_dofs; i++)
                    for (unsigned int j=0; j<n_neighbor_dofs; j++)
                      {
                        context.get_elem_neighbor_jacobian()(i,j) +=
                          JxW_face[qp] * phi_face[i][qp] * phi_neighbor_face[j][qp];
                      }

                  for (unsigned int i=0; i<n_neighbor_dofs; i++)
                    for (unsigned int j=0; j<n_neighbor_dofs; j++)
                      {
                        context.get_neighbor_neighbor_jacobian()(i,j) +=
                          JxW_face[qp] * phi_neighbor_face[i][qp] * phi_neighbor_face[j][qp];
                      }

                  for (unsigned int i=0; i<n_neighbor_dofs; i++)
                    for (unsigned int j=0; j<n_dofs; j++)
                      {
                        context.get_neighbor_elem_jacobian()(i,j) +=
                          JxW_face[qp] * phi_neighbor_face[i][qp] * phi_face[j][qp];
                      }
                }

              matrix.add_matrix (context.get_elem_elem_jacobian(),
                                 context.get_dof_indices(),
                                 context.get_dof_indices());

              matrix.add_matrix (context.get_elem_neighbor_jacobian(),
                                 context.get_dof_indices(),
                                 context.get_neighbor_dof_indices());

              matrix.add_matrix (context.get_neighbor_elem_jacobian(),
                                 context.get_neighbor_dof_indices(),
                                 context.get_dof_indices());

              matrix.add_matrix (context.get_neighbor_neighbor_jacobian(),
                                 context.get_neighbor_dof_indices(),
                                 context.get_neighbor_dof_indices());
            }
        }
      }
}


Number cubic_test (const Point& p,
                   const Parameters&,
                   const std::string&,
                   const std::string&)
{
  const Real & x = p(0);
  const Real & y = LIBMESH_DIM > 1 ? p(1) : 0;
  const Real & z = LIBMESH_DIM > 2 ? p(2) : 0;

  return x*(1-x)*(1-x) + x*x*(1-y) + x*(1-y)*(1-z) + y*(1-y)*z + z*(1-z)*(1-z);
}


Number new_linear_test (const Point& p,
                    const Parameters&,
                    const std::string&,
                    const std::string&)
{
  const Real & x = p(0);
  const Real & y = LIBMESH_DIM > 1 ? p(1) : 0;
  const Real & z = LIBMESH_DIM > 2 ? p(2) : 0;

  return x + 2*y + 3*z - 1;
}


Number disc_thirds_test (const Point& p,
                         const Parameters&,
                         const std::string&,
                         const std::string&)
{
  const Real & x = p(0);
  const Real & y = LIBMESH_DIM > 1 ? p(1) : 0;
  const Real & z = LIBMESH_DIM > 2 ? p(2) : 0;

  return (3*x < 1) + (3*y < 2) + (3*z > 2);
}


struct TripleFunction : public FunctionBase<Number>
{
  TripleFunction(Number _offset = 0) : offset(_offset) {}

  virtual std::unique_ptr<FunctionBase<Number>> clone () const
  { return std::make_unique<TripleFunction>(offset); }

  // We only really need the vector-valued output for projections
  virtual Number operator() (const Point &,
                             const Real /*time*/ = 0.) override
  { libmesh_error(); }

  virtual void operator() (const Point & p,
                           const Real,
                           DenseVector<Number> & output) override
  {
    libmesh_assert_greater(output.size(), 0);
    Parameters params;
    output(0) = cubic_test(p, params, "", "") + offset;
    if (output.size() > 0)
      output(1) = new_linear_test(p, params, "", "") + offset;
    if (output.size() > 1)
      output(2) = disc_thirds_test(p, params, "", "") + offset;
  }

  Number component (unsigned int i,
                    const Point & p,
                    Real /* time */) override
  {
    Parameters params;
    switch (i) {
    case 0:
      return cubic_test(p, params, "", "") + offset;
    case 1:
      return new_linear_test(p, params, "", "") + offset;
    case 2:
      return disc_thirds_test(p, params, "", "") + offset;
    default:
      libmesh_error();
    }
    return 0;
  }

  Number offset;
};


class SystemsTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( SystemsTest );

  CPPUNIT_TEST( test100KVariables );

  CPPUNIT_TEST( testPostInitAddVector );
  CPPUNIT_TEST( testAddVectorProjChange );
  CPPUNIT_TEST( testAddVectorTypeChange );
  CPPUNIT_TEST( testPostInitAddVectorTypeChange );

  CPPUNIT_TEST( testProjectHierarchicEdge3 );
#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testProjectHierarchicQuad9 );
  CPPUNIT_TEST( testProjectHierarchicTri6 );
  CPPUNIT_TEST( testProjectHierarchicTri7 );
  CPPUNIT_TEST( test2DProjectVectorFETri3 );
  CPPUNIT_TEST( test2DProjectVectorFEQuad4 );
  CPPUNIT_TEST( test2DProjectVectorFETri6 );
  CPPUNIT_TEST( test2DProjectVectorFETri7 );
  CPPUNIT_TEST( test2DProjectVectorFEQuad8 );
  CPPUNIT_TEST( test2DProjectVectorFEQuad9 );
#ifdef LIBMESH_HAVE_SOLVER
  CPPUNIT_TEST( testBlockRestrictedVarNDofs );
#endif
#endif // LIBMESH_DIM > 1
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testProjectHierarchicHex27 );
  CPPUNIT_TEST( testProjectMeshFunctionHex27 );
  CPPUNIT_TEST( testBoundaryProjectCube );
  CPPUNIT_TEST( test3DProjectVectorFETet4 );
  CPPUNIT_TEST( test3DProjectVectorFEHex8 );
  CPPUNIT_TEST( test3DProjectVectorFETet10 );
  CPPUNIT_TEST( test3DProjectVectorFETet14 );
  CPPUNIT_TEST( test3DProjectVectorFEHex20 );
  CPPUNIT_TEST( test3DProjectVectorFEHex27 );
#ifdef LIBMESH_HAVE_SOLVER
  CPPUNIT_TEST( testSetSystemParameterOverEquationSystem);
  CPPUNIT_TEST( testAssemblyWithDgFemContext );
#endif
#endif // LIBMESH_DIM > 2
#ifdef LIBMESH_HAVE_SOLVER
  CPPUNIT_TEST( testDofCouplingWithVarGroups );
#endif

#ifdef LIBMESH_ENABLE_AMR
#ifdef LIBMESH_HAVE_METAPHYSICL
#ifdef LIBMESH_HAVE_PETSC
  CPPUNIT_TEST( testProjectMatrixEdge2 );
#if LIBMESH_DIM > 1
  CPPUNIT_TEST( testProjectMatrixQuad4 );
  CPPUNIT_TEST( testProjectMatrixTri3 );
#endif // LIBMESH_DIM > 1
#if LIBMESH_DIM > 2
  CPPUNIT_TEST( testProjectMatrixHex8 );
  CPPUNIT_TEST( testProjectMatrixTet4 );
#endif // LIBMESH_DIM > 2
#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_HAVE_METAPHYSICL
#endif // LIBMESH_ENABLE_AMR

  CPPUNIT_TEST_SUITE_END();

private:
  void tripleValueTest (const Point & p,
                        const TransientExplicitSystem & sys,
                        const PointLocatorBase & locator,
                        std::set<subdomain_id_type> & u_subdomains,
                        std::set<subdomain_id_type> & v_subdomains,
                        std::set<subdomain_id_type> & w_subdomains,
                        const Parameters & param)
  {
    const Elem * elem = locator(p);
    subdomain_id_type sbd_id = elem ? elem->subdomain_id() : 0;
    TestCommWorld->max(sbd_id);

    if (u_subdomains.count(sbd_id))
      {
        LIBMESH_ASSERT_NUMBERS_EQUAL(cubic_test(p,param,"",""),
                                sys.point_value(0,p),
                                TOLERANCE*TOLERANCE*10);
        LIBMESH_ASSERT_NUMBERS_EQUAL
          (cubic_test(p,param,"","") + Number(10),
           sys.point_value(0,p,sys.old_local_solution),
           TOLERANCE*TOLERANCE*100);
        LIBMESH_ASSERT_NUMBERS_EQUAL
          (cubic_test(p,param,"","") + Number(20),
           sys.point_value(0,p,sys.older_local_solution),
           TOLERANCE*TOLERANCE*100);
      }
    if (v_subdomains.count(sbd_id))
      {
        LIBMESH_ASSERT_NUMBERS_EQUAL
          (new_linear_test(p,param,"",""), sys.point_value(1,p),
           TOLERANCE*TOLERANCE*10);
        LIBMESH_ASSERT_NUMBERS_EQUAL
          (new_linear_test(p,param,"","") + Number(10),
           sys.point_value(1,p,sys.old_local_solution),
           TOLERANCE*TOLERANCE*100);
        LIBMESH_ASSERT_NUMBERS_EQUAL
          (new_linear_test(p,param,"","") + Number(20),
           sys.point_value(1,p,sys.older_local_solution),
           TOLERANCE*TOLERANCE*100);
      }
    if (w_subdomains.count(sbd_id))
      {
        LIBMESH_ASSERT_NUMBERS_EQUAL
          (disc_thirds_test(p,param,"",""), sys.point_value(2,p),
           TOLERANCE*TOLERANCE*10);
        LIBMESH_ASSERT_NUMBERS_EQUAL
          (disc_thirds_test(p,param,"","") + Number(10),
           sys.point_value(2,p,sys.old_local_solution),
           TOLERANCE*TOLERANCE*100);
        LIBMESH_ASSERT_NUMBERS_EQUAL
          (disc_thirds_test(p,param,"","") + Number(20),
           sys.point_value(2,p,sys.older_local_solution),
           TOLERANCE*TOLERANCE*100);
      }
  }

public:
  void setUp()
  {}

  void tearDown()
  {}


  ExplicitSystem & simpleSetup(UnstructuredMesh & mesh,
                               EquationSystems & es)
  {
    ExplicitSystem &sys =
      es.add_system<ExplicitSystem> ("Simple");

    sys.add_variable("u", FIRST);

    MeshTools::Generation::build_line (mesh,
                                       10,
                                       0., 1.,
                                       EDGE3);

    return sys;
  }


  void test100KVariables()
  {
    Mesh mesh(*TestCommWorld);
    EquationSystems es(mesh);
    ExplicitSystem &sys =
      es.add_system<ExplicitSystem> ("100KVars");

    dof_id_type n_dofs = 100000;

    // This takes 16 seconds in opt mode for me, 200 in dbg!?
    /*
    for (auto i : make_range(n_dofs))
      sys.add_variable(std::to_string(i), FIRST);
    */

    std::vector<std::string> var_names(n_dofs);
    for (auto i : make_range(n_dofs))
      var_names[i] = std::to_string(i);

    sys.add_variables(var_names, FIRST);

    MeshTools::Generation::build_line (mesh,
                                       4,
                                       0., 1.,
                                       EDGE3);

    es.init();

    CPPUNIT_ASSERT_EQUAL(sys.n_dofs(), n_dofs*5);
    for (const Node * node : mesh.node_ptr_range())
      CPPUNIT_ASSERT_EQUAL(dof_id_type(node->n_vars(0)), n_dofs);

    std::vector<dof_id_type> each = sys.get_dof_map().n_dofs_per_processor(888);
    CPPUNIT_ASSERT_EQUAL(std::accumulate(each.begin(), each.end(), dof_id_type(0)), dof_id_type(5));
    CPPUNIT_ASSERT_EQUAL(sys.get_dof_map().n_dofs(888), dof_id_type(5));
  }


  void testPostInitAddVector()
  {
    Mesh mesh(*TestCommWorld);
    EquationSystems es(mesh);
    ExplicitSystem & sys = simpleSetup(mesh, es);
    es.init();

    auto & late_vec = sys.add_vector("late");

    CPPUNIT_ASSERT_EQUAL(sys.n_dofs(), dof_id_type(11));

    // late_vec should be initialized
    CPPUNIT_ASSERT_EQUAL(late_vec.size(), dof_id_type(11));
    CPPUNIT_ASSERT_EQUAL(late_vec.local_size(), sys.solution->local_size());
  }


  void testAddVectorProjChange()
  {
    Mesh mesh(*TestCommWorld);
    EquationSystems es(mesh);
    ExplicitSystem & sys = simpleSetup(mesh, es);

    sys.add_vector("late", /* projections = */ false);
    CPPUNIT_ASSERT_EQUAL(sys.vector_preservation("late"), false);

    sys.add_vector("late");
    CPPUNIT_ASSERT_EQUAL(sys.vector_preservation("late"), true);
  }


  void testAddVectorTypeChange()
  {
    // Vector types are all pretty much equivalent in serial.
    if (TestCommWorld->size() == 1)
      return;

    Mesh mesh(*TestCommWorld);
    EquationSystems es(mesh);
    ExplicitSystem & sys = simpleSetup(mesh, es);

    auto & late_vec = sys.add_vector("late");
    CPPUNIT_ASSERT_EQUAL(late_vec.type(), PARALLEL);
    CPPUNIT_ASSERT_EQUAL(sys.vector_preservation("late"), true);

    // We should never downgrade projection settings, so "false"
    // should be safely ignored here
    sys.add_vector("late", false, GHOSTED);
    CPPUNIT_ASSERT_EQUAL(late_vec.type(), GHOSTED);
    CPPUNIT_ASSERT_EQUAL(sys.vector_preservation("late"), true);

    // We should never downgrade storage settings, so this should be a
    // no-op.
    sys.add_vector("late", true, PARALLEL);
    CPPUNIT_ASSERT_EQUAL(late_vec.type(), GHOSTED);
    CPPUNIT_ASSERT_EQUAL(sys.vector_preservation("late"), true);
  }


  void testPostInitAddVectorTypeChange()
  {
    // Vector types are all pretty much equivalent in serial.
    if (TestCommWorld->size() == 1)
      return;

    Mesh mesh(*TestCommWorld);
    EquationSystems es(mesh);
    ExplicitSystem & sys = simpleSetup(mesh, es);

    auto & late_vec = sys.add_vector("late");
    CPPUNIT_ASSERT_EQUAL(late_vec.type(), PARALLEL);
    CPPUNIT_ASSERT_EQUAL(sys.vector_preservation("late"), true);

    es.init();

    auto & dof_map = sys.get_dof_map();

    // Set some data to make sure it's preserved
    CPPUNIT_ASSERT_EQUAL(sys.n_dofs(), dof_id_type(11));
    for (auto i : make_range(dof_map.first_dof(),
                             dof_map.end_dof()))
      late_vec.set(i, 2.0*i);
    late_vec.close();

    sys.add_vector("late", false, GHOSTED);
    CPPUNIT_ASSERT_EQUAL(late_vec.type(), GHOSTED);

    std::vector<dof_id_type> dof_indices;
    for (auto & elem : mesh.active_local_element_ptr_range())
      {
        dof_map.dof_indices (elem, dof_indices);

        for (auto d : dof_indices)
          CPPUNIT_ASSERT_EQUAL(late_vec(d), Number(2.0*d));
      }
  }



  void testProjectLine(const ElemType elem_type)
  {
    Mesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    TransientExplicitSystem &sys =
      es.add_system<TransientExplicitSystem> ("SimpleSystem");

    std::set<subdomain_id_type> u_subdomains {0, 1, 4, 5},
                                v_subdomains {1, 2, 3, 4},
                                w_subdomains {0, 1, 2, 3, 4};

    sys.add_variable("u", THIRD,    HIERARCHIC, &u_subdomains);
    sys.add_variable("v", FIRST,    LAGRANGE,   &v_subdomains);
    sys.add_variable("w", CONSTANT, MONOMIAL,   &w_subdomains);

    MeshTools::Generation::build_line (mesh,
                                       6,
                                       0., 1.,
                                       elem_type);

    for (auto & elem : mesh.element_ptr_range())
      elem->subdomain_id() = elem->id();

    es.init();
    TripleFunction tfunc;
    sys.project_solution(&tfunc);
    tfunc.offset = 10;
    sys.project_vector(*sys.old_local_solution, &tfunc);
    tfunc.offset = 20;
    sys.project_vector(*sys.older_local_solution, &tfunc);

    std::unique_ptr<PointLocatorBase> locator = mesh.sub_point_locator();
    locator->enable_out_of_mesh_mode();
    for (Real x = 0.1; x < 1; x += 0.2)
      tripleValueTest(Point(x), sys, *locator,
                      u_subdomains, v_subdomains, w_subdomains,
                      es.parameters);

#ifdef LIBMESH_ENABLE_AMR
    for (auto & elem : mesh.element_ptr_range())
      if ((elem->id()/2)%2)
        elem->set_refinement_flag(Elem::REFINE);
    es.reinit();

    locator = mesh.sub_point_locator();
    locator->enable_out_of_mesh_mode();
    for (Real x = 0.1; x < 1; x += 0.2)
      tripleValueTest(Point(x), sys, *locator,
                      u_subdomains, v_subdomains, w_subdomains,
                      es.parameters);
#endif
  }

  void test2DProjectVectorFE(const ElemType elem_type)
  {
    Mesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    TransientExplicitSystem &sys =
      es.add_system<TransientExplicitSystem> ("SimpleSystem");

    auto u_var = sys.add_variable("u", Elem::type_to_default_order_map[elem_type], LAGRANGE_VEC);

    MeshTools::Generation::build_square (mesh,
                                         1, 1,
                                         0., 1., 0., 1.,
                                         elem_type);

    es.init();

    // Manually set-up the solution because I'm too lazy to set-up all the generic
    // function projection code right now
    for (const auto & node : mesh.local_node_ptr_range())
    {
      for (unsigned int i : make_range(Elem::type_to_dim_map[elem_type]))
      {
        auto dof_index = node->dof_number(sys.number(), u_var, i);
        sys.solution->set(dof_index, (*node)(i));
      }
    }

    // After setting values, we need to assemble
    sys.solution->close();


#ifdef LIBMESH_ENABLE_AMR
    for (auto & elem : mesh.element_ptr_range())
        elem->set_refinement_flag(Elem::REFINE);
    es.reinit();
#endif

    for (const auto & node : mesh.local_node_ptr_range())
    {
      // 2D element here
      for (unsigned int i : make_range(Elem::type_to_dim_map[elem_type]))
      {
        auto dof_index = node->dof_number(sys.number(), u_var, i);
        auto value = (*sys.solution)(dof_index);
        LIBMESH_ASSERT_NUMBERS_EQUAL(value, (*node)(i), TOLERANCE*TOLERANCE);
      }
    }
  }

  void test3DProjectVectorFE(const ElemType elem_type)
  {
    Mesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    TransientExplicitSystem &sys =
      es.add_system<TransientExplicitSystem> ("SimpleSystem");

    auto u_var = sys.add_variable
      ("u", Elem::type_to_default_order_map[elem_type], LAGRANGE_VEC);

    MeshTools::Generation::build_cube (mesh,
                                       1, 1, 1,
                                       0., 1., 0., 1., 0., 1.,
                                       elem_type);

    es.init();

    // Manually set-up the solution because I'm too lazy to set-up all the generic
    // function projection code right now
    for (const auto & node : mesh.local_node_ptr_range())
    {
      for (unsigned int i : make_range(Elem::type_to_dim_map[elem_type]))
      {
        auto dof_index = node->dof_number(sys.number(), u_var, i);
        sys.solution->set(dof_index, (*node)(i));
      }
    }

    // After setting values, we need to assemble
    sys.solution->close();


#ifdef LIBMESH_ENABLE_AMR
    for (auto & elem : mesh.element_ptr_range())
        elem->set_refinement_flag(Elem::REFINE);
    es.reinit();
#endif

    for (const auto & node : mesh.local_node_ptr_range())
    {
      for (unsigned int i : make_range(Elem::type_to_dim_map[elem_type]))
      {
        auto dof_index = node->dof_number(sys.number(), u_var, i);
        auto value = (*sys.solution)(dof_index);
        LIBMESH_ASSERT_NUMBERS_EQUAL(value, (*node)(i), TOLERANCE*TOLERANCE);
      }
    }
  }

  void testProjectSquare(const ElemType elem_type)
  {
    Mesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    TransientExplicitSystem &sys =
      es.add_system<TransientExplicitSystem> ("SimpleSystem");

    std::set<subdomain_id_type> u_subdomains {0, 1, 4, 5},
                                v_subdomains {1, 2, 3, 4},
                                w_subdomains {0, 1, 2, 3, 4};

    sys.add_variable("u", THIRD,    HIERARCHIC, &u_subdomains);
    sys.add_variable("v", FIRST,    LAGRANGE,   &v_subdomains);
    sys.add_variable("w", CONSTANT, MONOMIAL,   &w_subdomains);

    MeshTools::Generation::build_square (mesh,
                                         3, 3,
                                         0., 1., 0., 1.,
                                         elem_type);

    for (auto & elem : mesh.element_ptr_range())
      elem->subdomain_id() = elem->id()/2;

    es.init();
    TripleFunction tfunc;
    sys.project_solution(&tfunc);
    tfunc.offset = 10;
    sys.project_vector(*sys.old_local_solution, &tfunc);
    tfunc.offset = 20;
    sys.project_vector(*sys.older_local_solution, &tfunc);

    std::unique_ptr<PointLocatorBase> locator = mesh.sub_point_locator();
    locator->enable_out_of_mesh_mode();
    for (Real x = 0.1; x < 1; x += 0.2)
      for (Real y = 0.1; y < 1; y += 0.2)
        tripleValueTest(Point(x,y), sys, *locator,
                        u_subdomains, v_subdomains, w_subdomains,
                        es.parameters);

#ifdef LIBMESH_ENABLE_AMR
    for (auto & elem : mesh.element_ptr_range())
      if ((elem->id()/2)%2)
        elem->set_refinement_flag(Elem::REFINE);
    es.reinit();

    locator = mesh.sub_point_locator();
    locator->enable_out_of_mesh_mode();
    for (Real x = 0.1; x < 1; x += 0.2)
      for (Real y = 0.1; y < 1; y += 0.2)
        tripleValueTest(Point(x,y), sys, *locator,
                        u_subdomains, v_subdomains, w_subdomains,
                        es.parameters);
#endif
  }

  void testProjectCube(const ElemType elem_type)
  {
    Mesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    TransientExplicitSystem &sys =
      es.add_system<TransientExplicitSystem> ("SimpleSystem");

    std::set<subdomain_id_type> u_subdomains {0, 1, 4, 5},
                                v_subdomains {1, 2, 3, 4},
                                w_subdomains {0, 1, 2, 3, 4};

    sys.add_variable("u", THIRD,    HIERARCHIC, &u_subdomains);
    sys.add_variable("v", FIRST,    LAGRANGE,   &v_subdomains);
    sys.add_variable("w", CONSTANT, MONOMIAL,   &w_subdomains);

    MeshTools::Generation::build_cube (mesh,
                                       3, 3, 3,
                                       0., 1., 0., 1., 0., 1.,
                                       elem_type);

    for (auto & elem : mesh.element_ptr_range())
      elem->subdomain_id() = elem->id()/6;

    es.init();
    TripleFunction tfunc;
    sys.project_solution(&tfunc);
    tfunc.offset = 10;
    sys.project_vector(*sys.old_local_solution, &tfunc);
    tfunc.offset = 20;
    sys.project_vector(*sys.older_local_solution, &tfunc);

    std::unique_ptr<PointLocatorBase> locator = mesh.sub_point_locator();
    locator->enable_out_of_mesh_mode();
    for (Real x = 0.1; x < 1; x += 0.2)
      for (Real y = 0.1; y < 1; y += 0.2)
        for (Real z = 0.1; z < 1; z += 0.2)
          tripleValueTest(Point(x,y,z), sys, *locator,
                          u_subdomains, v_subdomains, w_subdomains,
                          es.parameters);

  #ifdef LIBMESH_ENABLE_AMR
    for (auto & elem : mesh.element_ptr_range())
      if ((elem->id()/2)%2)
        elem->set_refinement_flag(Elem::REFINE);
    es.reinit();

    locator = mesh.sub_point_locator();
    locator->enable_out_of_mesh_mode();
    for (Real x = 0.1; x < 1; x += 0.2)
      for (Real y = 0.1; y < 1; y += 0.2)
        for (Real z = 0.1; z < 1; z += 0.2)
          tripleValueTest(Point(x,y,z), sys, *locator,
                          u_subdomains, v_subdomains, w_subdomains,
                          es.parameters);
  #endif
  }

  void testProjectCubeWithMeshFunction(const ElemType elem_type)
  {
    // The source mesh needs to exist everywhere it's queried, so we
    // use a ReplicatedMesh
    ReplicatedMesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    System &sys = es.add_system<System> ("SimpleSystem");
    sys.add_variable("u", THIRD, MONOMIAL);

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

    // use 3rd order again so we can expect exact results
    proj_sys.add_variable("u", THIRD, HIERARCHIC);

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
            LIBMESH_ASSERT_NUMBERS_EQUAL
              (cubic_test(p,es.parameters,"",""),
               proj_sys.point_value(0,p), TOLERANCE*TOLERANCE);
          }
  }

  void testBoundaryProjectCube()
  {
    LOG_UNIT_TEST;

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

    LIBMESH_ASSERT_FP_EQUAL(sys.solution->l1_norm(), ref_l1_norm, TOLERANCE*TOLERANCE);
  }

  void testDofCouplingWithVarGroups()
  {
    LOG_UNIT_TEST;

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
    auto new_edge_elem = mesh.add_elem(Elem::build(EDGE2));
    new_edge_elem->set_node(0, new_node_a);
    new_edge_elem->set_node(1, new_node_b);

    mesh.elem_ref(0).subdomain_id() = 10;
    mesh.elem_ref(1).subdomain_id() = 10;

    // Add NodeElems for coupling purposes
    auto node_elem_1 = mesh.add_elem(Elem::build(NODEELEM));
    node_elem_1->set_node(0, mesh.elem_ref(0).node_ptr(1));
    auto node_elem_2 = mesh.add_elem(Elem::build(NODEELEM));
    node_elem_2->set_node(0, new_node_a);

    mesh.prepare_for_use();

    // Create an equation systems object.
    EquationSystems equation_systems (mesh);
    LinearImplicitSystem & system =
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

    // LASPACK GMRES + ILU defaults don't like this problem, but it's
    // small enough to just use a simpler iteration.
    system.get_linear_solver()->set_solver_type(JACOBI);
    system.get_linear_solver()->set_preconditioner_type(IDENTITY_PRECOND);

    equation_systems.init ();

    system.solve();

    // We set the solution to be 1 everywhere, so the final l1 norm of the
    // solution is the product of the number of variables and number of nodes.
    Real ref_l1_norm = static_cast<Real>(mesh.n_nodes() * system.n_vars());

    LIBMESH_ASSERT_FP_EQUAL(system.solution->l1_norm(), ref_l1_norm, TOLERANCE*TOLERANCE);
  }


  void testSetSystemParameterOverEquationSystem()
  {
    LOG_UNIT_TEST;

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
    auto new_edge_elem = mesh.add_elem(Elem::build(EDGE2));
    new_edge_elem->set_node(0, new_node_a);
    new_edge_elem->set_node(1, new_node_b);

    mesh.elem_ref(0).subdomain_id() = 10;
    mesh.elem_ref(1).subdomain_id() = 10;

    mesh.prepare_for_use();

    // Create an equation systems object.
    EquationSystems equation_systems (mesh);

    // Set some parameters to the equation system that would cause a failed test
    equation_systems.parameters.set<unsigned int>("linear solver maximum iterations") = 0;

    // Setup Linear Implicit system
    LinearImplicitSystem & li_system =
      equation_systems.add_system<LinearImplicitSystem> ("test");

    // We must use a discontinuous variable type in this test or
    // else the sparsity pattern will not be correct
    li_system.add_variable("u", FIRST, L2_LAGRANGE);

    MeshTools::Generation::build_cube (mesh,
                                       5, 5, 5,
                                       0., 1., 0., 1., 0., 1.,
                                       HEX8);

    li_system.attach_assemble_function (assembly_with_dg_fem_context);
    li_system.get_linear_solver()->set_solver_type(GMRES);
    // Need 5 iterations, dont overdo the preconditioning
    li_system.get_linear_solver()->set_preconditioner_type(IDENTITY_PRECOND);

    // Set some parameters to the system that work for the solve
    li_system.parameters.set<unsigned int>("linear solver maximum iterations") = 5;
    li_system.parameters.set<Real>("linear solver tolerance") = 1e-100;

    // Need to init before we can access the system matrix
    equation_systems.init ();

    // See the solve pass, indicating system parameters are used over equation system parameters
    li_system.solve();

    // Check that the number of iterations from the systems got obeyed
    CPPUNIT_ASSERT_EQUAL(li_system.n_linear_iterations(), 5u);
}

  void testAssemblyWithDgFemContext()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);

    EquationSystems es(mesh);
    System &sys = es.add_system<LinearImplicitSystem> ("test");

    // We must use a discontinuous variable type in this test or
    // else the sparsity pattern will not be correct
    sys.add_variable("u", FIRST, L2_LAGRANGE);

    MeshTools::Generation::build_cube (mesh,
                                       5, 5, 5,
                                       0., 1., 0., 1., 0., 1.,
                                       HEX8);

    es.init();
    sys.attach_assemble_function (assembly_with_dg_fem_context);
    sys.solve();

    // We don't actually assert anything in this test. We just want to check that
    // the assembly and solve do not encounter any errors.
  }

  void testBlockRestrictedVarNDofs()
  {
    LOG_UNIT_TEST;

    ReplicatedMesh mesh(*TestCommWorld);

    MeshTools::Generation::build_cube (mesh,
                                      4,
                                      4,
                                      0,
                                      0., 1.,
                                      0., 1.,
                                      0., 0.,
                                      QUAD4);

    for (const auto & elem : mesh.element_ptr_range())
      {
        Point c = elem->vertex_average();
        if (c(0) <= 0.5 && c(1) <= 0.5)
          elem->subdomain_id() = 0;
        else
          elem->subdomain_id() = 1;
      }

    mesh.prepare_for_use();

    // Create an equation systems object.
    EquationSystems equation_systems (mesh);
    ExplicitSystem& system =
      equation_systems.add_system<LinearImplicitSystem> ("test");

    std::set<subdomain_id_type> block0;
    std::set<subdomain_id_type> block1;
    block0.insert(0);
    block1.insert(1);
    auto u0 = system.add_variable ("u0", libMesh::FIRST, &block0);
    auto u1 = system.add_variable ("u1", libMesh::FIRST, &block1);
    equation_systems.init();

    std::vector<dof_id_type> u0_dofs;
    system.get_dof_map().local_variable_indices(u0_dofs, mesh, u0);
    std::vector<dof_id_type> u1_dofs;
    system.get_dof_map().local_variable_indices(u1_dofs, mesh, u1);

    std::set<dof_id_type> sys_u0_dofs;
    system.local_dof_indices(u0, sys_u0_dofs);
    std::set<dof_id_type> sys_u1_dofs;
    system.local_dof_indices(u1, sys_u1_dofs);

    // Get local indices from other processors too
    mesh.comm().allgather(u0_dofs);
    mesh.comm().allgather(u1_dofs);
    mesh.comm().set_union(sys_u0_dofs);
    mesh.comm().set_union(sys_u1_dofs);

    const std::size_t c9 = 9;
    const std::size_t c21 = 21;
    CPPUNIT_ASSERT_EQUAL(c9, u0_dofs.size());
    CPPUNIT_ASSERT_EQUAL(c21, u1_dofs.size());
    CPPUNIT_ASSERT_EQUAL(c9, sys_u0_dofs.size());
    CPPUNIT_ASSERT_EQUAL(c21, sys_u1_dofs.size());
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
    int n_new_dofs_local = sys.get_dof_map().n_local_dofs();
    int ndofs_old_first = sys.get_dof_map().first_old_dof();
    int ndofs_old_end   = sys.get_dof_map().end_old_dof();
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
    int n_new_dofs_local = sys.get_dof_map().n_local_dofs();
    int ndofs_old_first = sys.get_dof_map().first_old_dof();
    int ndofs_old_end   = sys.get_dof_map().end_old_dof();
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

        Elem * elem = mesh.add_elem(Elem::build_with_id(TET4, 0));
        elem->set_node(0, mesh.node_ptr(0));
        elem->set_node(1, mesh.node_ptr(1));
        elem->set_node(2, mesh.node_ptr(2));
        elem->set_node(3, mesh.node_ptr(3));

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


  void testProjectHierarchicEdge3() { LOG_UNIT_TEST; testProjectLine(EDGE3); }
  void testProjectHierarchicQuad9() { LOG_UNIT_TEST; testProjectSquare(QUAD9); }
  void testProjectHierarchicTri6()  { LOG_UNIT_TEST; testProjectSquare(TRI6); }
  void testProjectHierarchicTri7()  { LOG_UNIT_TEST; testProjectSquare(TRI7); }
  void testProjectHierarchicHex27() { LOG_UNIT_TEST; testProjectCube(HEX27); }
  void testProjectMeshFunctionHex27() { LOG_UNIT_TEST; testProjectCubeWithMeshFunction(HEX27); }
  void test2DProjectVectorFETri3() { LOG_UNIT_TEST; test2DProjectVectorFE(TRI3); }
  void test2DProjectVectorFEQuad4() { LOG_UNIT_TEST; test2DProjectVectorFE(QUAD4); }
  void test2DProjectVectorFETri6() { LOG_UNIT_TEST; test2DProjectVectorFE(TRI6); }
  void test2DProjectVectorFETri7() { LOG_UNIT_TEST; test2DProjectVectorFE(TRI7); }
  void test2DProjectVectorFEQuad8() { LOG_UNIT_TEST; test2DProjectVectorFE(QUAD8); }
  void test2DProjectVectorFEQuad9() { LOG_UNIT_TEST; test2DProjectVectorFE(QUAD9); }
  void test3DProjectVectorFETet4() { LOG_UNIT_TEST; test3DProjectVectorFE(TET4); }
  void test3DProjectVectorFEHex8() { LOG_UNIT_TEST; test3DProjectVectorFE(HEX8); }
  void test3DProjectVectorFETet10() { LOG_UNIT_TEST; test3DProjectVectorFE(TET10); }
  void test3DProjectVectorFETet14() { LOG_UNIT_TEST; test3DProjectVectorFE(TET14); }
  void test3DProjectVectorFEHex20() { LOG_UNIT_TEST; test3DProjectVectorFE(HEX20); }
  void test3DProjectVectorFEHex27() { LOG_UNIT_TEST; test3DProjectVectorFE(HEX27); }

#ifdef LIBMESH_ENABLE_AMR
#ifdef LIBMESH_HAVE_METAPHYSICL
#ifdef LIBMESH_HAVE_PETSC
  // projection matrix tests
  void testProjectMatrixEdge2() { LOG_UNIT_TEST; testProjectMatrix1D(EDGE2); }
  void testProjectMatrixQuad4() { LOG_UNIT_TEST; testProjectMatrix2D(QUAD4); }
  void testProjectMatrixTri3() { LOG_UNIT_TEST; testProjectMatrix2D(TRI3); }
  void testProjectMatrixHex8() { LOG_UNIT_TEST; testProjectMatrix3D(HEX8); }
  void testProjectMatrixTet4() { LOG_UNIT_TEST; testProjectMatrix3D(TET4); }
#endif // LIBMESH_HAVE_PETSC
#endif // LIBMESH_HAVE_METAPHYSICL
#endif // LIBMESH_ENABLE_AMR

};

CPPUNIT_TEST_SUITE_REGISTRATION( SystemsTest );
