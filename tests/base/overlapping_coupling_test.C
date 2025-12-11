#include <libmesh/equation_systems.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/numeric_vector.h>
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/face_quad4.h>
#include <libmesh/ghosting_functor.h>
#include <libmesh/coupling_matrix.h>
#include <libmesh/quadrature_gauss.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/mesh_modification.h>
#include <libmesh/partitioner.h>
#include <libmesh/sparse_matrix.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"

#include <memory>

using namespace libMesh;


// This coupling functor is attempting to mimic in a basic way
// a particular use case. The idea is that the mesh has overlapping domains
// (we assume only two in this test) and that the elements are not
// topologically connected between the overlapping domains, but when the user
// solves their problem, some dofs on are coupled between the overlapping
// elements. So we want to exercise the GhostingFunctor for this scenario,
// in particular both algebraic ghosting (dof ghosting) and coupling
// ghosting (sparsity pattern). In the use case, the coupling happens at
// quadrature points, so we look for the overlap based on the quadrature points.
//
// Additionally, the use case also has an issue that the coupling functor
// cannot be added until after the System is initialized because the overlap
// depends on the solution. So we will exercise that case as well.
class OverlappingCouplingFunctor : public GhostingFunctor
{
  public:

  OverlappingCouplingFunctor( System & system )
    : GhostingFunctor(),
      _system(system),
      _mesh(system.get_mesh()),
      _point_locator(system.get_mesh().sub_point_locator()),
      _coupling_matrix(nullptr)
  {
    // This is a highly specalized coupling functor where we are assuming
    // only two subdomains that are overlapping so make sure there's only two.
    CPPUNIT_ASSERT_EQUAL(2, static_cast<int>(_mesh.n_subdomains()));

    // We're assuming a specific set of subdomain ids
    std::set<subdomain_id_type> ids;
    _mesh.subdomain_ids(ids);
    CPPUNIT_ASSERT( ids.find(1) != ids.end() );
    CPPUNIT_ASSERT( ids.find(2) != ids.end() );

    _point_locator->enable_out_of_mesh_mode();
  }

  virtual ~OverlappingCouplingFunctor(){};

  virtual std::unique_ptr<GhostingFunctor> clone () const override
  {
    auto clone_functor =
      std::make_unique<OverlappingCouplingFunctor>(_system);

    auto coupling_matrix =
      std::make_unique<CouplingMatrix>(*_coupling_matrix);
    clone_functor->set_coupling_matrix(coupling_matrix);
    return clone_functor;
  }

  void set_coupling_matrix (std::unique_ptr<CouplingMatrix> & coupling_matrix)
  { _coupling_matrix = std::move(coupling_matrix); }

  virtual void operator()
  (const MeshBase::const_element_iterator & range_begin,
   const MeshBase::const_element_iterator & range_end,
   processor_id_type p,
   map_type & coupled_elements) override
  {
    std::unique_ptr<PointLocatorBase> sub_point_locator = _mesh.sub_point_locator();

    for( const auto & elem : as_range(range_begin,range_end) )
      {
        std::set<subdomain_id_type> allowed_subdomains;
        unsigned int var = 0;
        if( elem->subdomain_id() == 1 )
          {
            var = _system.variable_number("V");
            allowed_subdomains.insert(2);
          }
        else if( elem->subdomain_id() == 2 )
          {
            var = _system.variable_number("U");
            allowed_subdomains.insert(1);
          }
        else
          CPPUNIT_FAIL("Something bad happpend");

        unsigned int dim = 2;
        FEType fe_type = _system.variable_type(var);
        std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
        QGauss qrule (dim, fe_type.default_quadrature_order());
        fe->attach_quadrature_rule (&qrule);

        const std::vector<libMesh::Point> & qpoints = fe->get_xyz();

        fe->reinit(elem);

        for ( const auto & qp : qpoints )
          {
            const Elem * overlapping_elem = (*sub_point_locator)( qp, &allowed_subdomains );

            if( overlapping_elem && overlapping_elem->processor_id() != p )
              coupled_elements.emplace(overlapping_elem, _coupling_matrix.get());
          }
      }
  }

private:

  System & _system;

  const MeshBase & _mesh;

  std::unique_ptr<PointLocatorBase> _point_locator;

  std::unique_ptr<CouplingMatrix> _coupling_matrix;
};


// We use a custom partioner for this test to help ensure we're
// testing the cases we want. In particular, we're going to
// ensure that elements in subdomain one are on different processors
// from elements subdomain two so we know ghosting will be required
// for algebraic and coupling ghosting between the two subdomains.
class OverlappingTestPartitioner : public Partitioner
{
public:

  OverlappingTestPartitioner () = default;
  OverlappingTestPartitioner (const OverlappingTestPartitioner &) = default;
  OverlappingTestPartitioner (OverlappingTestPartitioner &&) = default;
  OverlappingTestPartitioner & operator= (const OverlappingTestPartitioner &) = default;
  OverlappingTestPartitioner & operator= (OverlappingTestPartitioner &&) = default;
  virtual ~OverlappingTestPartitioner() = default;

  /**
   * \returns A copy of this partitioner wrapped in a smart pointer.
   */
  virtual std::unique_ptr<Partitioner> clone () const override
  {
    return std::make_unique<OverlappingTestPartitioner>(*this);
  }

protected:

  /**
   * Partition the \p MeshBase into \p n subdomains.
   */
  virtual void _do_partition (MeshBase & mesh,
                              const unsigned int n) override
  {
    // If we're on one partition, then everyone gets to be on that partition
    if (n == 1)
      this->single_partition_range (mesh.active_elements_begin(), mesh.active_elements_end());
    else
      {
        libmesh_assert_greater (n, 0);

        // If there are more partitions than elements, then just
        // assign one element to each partition
        if (mesh.n_active_elem() <= n)
          {
            processor_id_type e = 0;
            for (auto & elem : mesh.active_element_ptr_range() )
              {
                elem->processor_id() = e;
                e++;
              }
          }
        // Otherwise, we'll split up the partitions into two groups
        // and then assign the elements of each subdomain into each of
        // the two groups
        else
          {
            unsigned int n_sub_two_elems = std::distance( mesh.active_subdomain_elements_begin(2),
                                                          mesh.active_subdomain_elements_end(2) );

            unsigned int n_parts_sub_two = 0;
            if( n_sub_two_elems < n/2 )
              n_parts_sub_two = n_sub_two_elems;
            else
                n_parts_sub_two = n/2;

            const unsigned int n_parts_sub_one = n - n_parts_sub_two;

            const dof_id_type sub_two_blk_size = cast_int<dof_id_type>
              (std::distance( mesh.active_subdomain_elements_begin(2),
                              mesh.active_subdomain_elements_end(2) )/n_parts_sub_two );

            this->assign_proc_id_subdomain( mesh, 2, sub_two_blk_size, n_parts_sub_two, 0 );


            const dof_id_type sub_one_blk_size = cast_int<dof_id_type>
              (std::distance( mesh.active_subdomain_elements_begin(1),
                              mesh.active_subdomain_elements_end(1) )/n_parts_sub_one );

            this->assign_proc_id_subdomain( mesh, 1, sub_one_blk_size, n_parts_sub_one, n_parts_sub_two );
          }
      }
  }

  void assign_proc_id_subdomain( MeshBase & mesh,
                                 const subdomain_id_type sid,
                                 const dof_id_type blksize,
                                 const unsigned int n_parts,
                                 const processor_id_type offset )
  {
    dof_id_type e = 0;
    for (auto & elem : mesh.active_subdomain_elements_ptr_range(sid))
      {
        if ((e/blksize) < n_parts)
          elem->processor_id() = offset + cast_int<processor_id_type>(e/blksize);
        else
          elem->processor_id() = offset;

        e++;
      }
  }

};


// Common functionality for all the different tests we'll be running
// We're creating two elements that will live in subdomain 1 and one
// element that is lives in subdomain 2 that overlaps both elements in
// subdomain 1, but is not topologically connected.
class OverlappingTestBase
{
protected:

  std::unique_ptr<MeshBase> _mesh;

  std::unique_ptr<EquationSystems> _es;

  void build_quad_mesh(unsigned int n_refinements = 0)
  {
    // We are making assumptions in various places about the presence
    // of the elements on the current processor so we're restricting to
    // ReplicatedMesh for now.
    _mesh = std::make_unique<ReplicatedMesh>(*TestCommWorld);

    _mesh->set_mesh_dimension(2);

    _mesh->add_point( Point(0.0,0.0),0 );
    _mesh->add_point( Point(1.0,0.0),1 );
    _mesh->add_point( Point(1.0,1.0),2 );
    _mesh->add_point( Point(0.0,1.0),3 );

    {
      Elem * elem = _mesh->add_elem(Elem::build_with_id(QUAD4, 0));
      elem->subdomain_id() = 1;

      for (unsigned int n=0; n<4; n++)
        elem->set_node(n, _mesh->node_ptr(n));
    }

    _mesh->add_point( Point(1.0,2.0),4 );
    _mesh->add_point( Point(0.0,2.0),5 );

    {
      Elem * elem = _mesh->add_elem(Elem::build_with_id(QUAD4, 1));
      elem->subdomain_id() = 1;

      elem->set_node(0, _mesh->node_ptr(3));
      elem->set_node(1, _mesh->node_ptr(2));
      elem->set_node(2, _mesh->node_ptr(4));
      elem->set_node(3, _mesh->node_ptr(5));
    }

    _mesh->add_point( Point(0.0,0.0),6 );
    _mesh->add_point( Point(1.0,0.0),7 );
    _mesh->add_point( Point(1.0,2.0),8 );
    _mesh->add_point( Point(0.0,2.0),9 );

    {
      Elem* elem = _mesh->add_elem(Elem::build_with_id(QUAD4, 2));
      elem->subdomain_id() = 2;

      elem->set_node(0, _mesh->node_ptr(6));
      elem->set_node(1, _mesh->node_ptr(7));
      elem->set_node(2, _mesh->node_ptr(8));
      elem->set_node(3, _mesh->node_ptr(9));
    }

    _mesh->partitioner() = std::make_unique<OverlappingTestPartitioner>();

    _mesh->prepare_for_use();

#ifdef LIBMESH_ENABLE_AMR
    if (n_refinements > 0)
      {
        MeshRefinement refine(*_mesh);
        refine.uniformly_refine(n_refinements);
      }
#else
    CPPUNIT_ASSERT_EQUAL(n_refinements, 0u);
#endif // LIBMESH_ENABLE_AMR
  }

  void init()
  {
    _es = std::make_unique<EquationSystems>(*_mesh);
    LinearImplicitSystem & sys = _es->add_system<LinearImplicitSystem> ("SimpleSystem");

    std::set<subdomain_id_type> sub_one;
    sub_one.insert(1);

    std::set<subdomain_id_type> sub_two;
    sub_two.insert(2);

    sys.add_variable("U", FIRST, LAGRANGE, &sub_two);
    sys.add_variable("L", FIRST, LAGRANGE, &sub_two);

    sys.add_variable("V", FIRST, LAGRANGE, &sub_one);
    sys.add_variable("p", FIRST, LAGRANGE, &sub_one);

    _es->init();
  }

  void clear()
  {
    _es.reset();
    _mesh.reset();
  }

  void setup_coupling_matrix( std::unique_ptr<CouplingMatrix> & coupling )
  {
    System & system = _es->get_system("SimpleSystem");

    coupling = std::make_unique<CouplingMatrix>(system.n_vars());

    const unsigned int u_var = system.variable_number("U");
    const unsigned int l_var = system.variable_number("L");
    const unsigned int v_var = system.variable_number("V");
    const unsigned int p_var = system.variable_number("p");

    // Only adding the overlapping couplings since the primary ones should
    // be there by default.
    (*coupling)(u_var,v_var) = true;
    (*coupling)(l_var,v_var) = true;
    (*coupling)(l_var,p_var) = true;
    (*coupling)(v_var,u_var) = true;
    (*coupling)(v_var,l_var) = true;
  }

};

// Mainly just to sanity check the mesh construction and
// the assumptions made in the OverlappingCouplingFunctor
// as well as the custom Partitioner.
class OverlappingFunctorTest : public CppUnit::TestCase,
                               public OverlappingTestBase
{
public:

  LIBMESH_CPPUNIT_TEST_SUITE( OverlappingFunctorTest );

#ifdef LIBMESH_HAVE_SOLVER
  CPPUNIT_TEST( checkCouplingFunctorQuad );
  CPPUNIT_TEST( checkCouplingFunctorTri );
  CPPUNIT_TEST( checkOverlappingPartitioner );
# ifdef LIBMESH_ENABLE_AMR
  CPPUNIT_TEST( checkCouplingFunctorQuadUnifRef );
  CPPUNIT_TEST( checkCouplingFunctorTriUnifRef );
  CPPUNIT_TEST( checkOverlappingPartitionerUnifRef );
# endif
#endif

  CPPUNIT_TEST_SUITE_END();

public:

  public:
  void setUp()
  {
    this->build_quad_mesh();
    this->init();
  }

  void tearDown()
  { this->clear(); }

  void checkCouplingFunctorQuad()
  { this->run_coupling_functor_test(0); }

  void checkCouplingFunctorQuadUnifRef()
  { this->run_coupling_functor_test(1); }

  void checkCouplingFunctorTri()
  {
    MeshTools::Modification::all_tri(*_mesh);
    _es->reinit();
    this->run_coupling_functor_test(0);
  }

  void checkCouplingFunctorTriUnifRef()
  {
    MeshTools::Modification::all_tri(*_mesh);
    _es->reinit();
    this->run_coupling_functor_test(1);
  }

  void checkOverlappingPartitioner()
  {
    this->run_partitioner_test(0);
  }

  void checkOverlappingPartitionerUnifRef()
  {
    this->run_partitioner_test(1);
  }

private:

  // This is basically to sanity check the coupling functor
  // with the supplied mesh to make sure all the assumptions
  // are kosher. This test requires AMR and some kind of a
  // linear solver, so let's only run it if we have PETSc.
  void run_coupling_functor_test(unsigned int n_refinements)
  {
#if defined(LIBMESH_ENABLE_AMR) && defined(LIBMESH_HAVE_SOLVER)
    if( n_refinements > 0 )
      {
        MeshRefinement refine(*_mesh);
        refine.uniformly_refine(n_refinements);
        _es->reinit();
      }
#else
    CPPUNIT_ASSERT_EQUAL(n_refinements, 0u);
#endif

    System & system = _es->get_system("SimpleSystem");

    OverlappingCouplingFunctor coupling_functor(system);

    GhostingFunctor::map_type subdomain_one_couplings;
    GhostingFunctor::map_type subdomain_two_couplings;

    coupling_functor( _mesh->active_subdomain_elements_begin(1),
                      _mesh->active_subdomain_elements_end(1),
                      DofObject::invalid_processor_id,
                      subdomain_one_couplings );

    coupling_functor( _mesh->active_subdomain_elements_begin(2),
                      _mesh->active_subdomain_elements_end(2),
                      DofObject::invalid_processor_id,
                      subdomain_two_couplings );

    dof_id_type n_elems_subdomain_two = std::distance( _mesh->active_subdomain_elements_begin(2),
                                                       _mesh->active_subdomain_elements_end(2) );

    CPPUNIT_ASSERT_EQUAL(n_elems_subdomain_two, static_cast<dof_id_type>(subdomain_one_couplings.size()));
    CPPUNIT_ASSERT_EQUAL(static_cast<dof_id_type>(2*n_elems_subdomain_two), static_cast<dof_id_type>(subdomain_two_couplings.size()));
  }

  void run_partitioner_test(unsigned int n_refinements)
  {
#ifdef LIBMESH_ENABLE_AMR
    if( n_refinements > 0 )
      {
        MeshRefinement refine(*_mesh);
        refine.uniformly_refine(n_refinements);
        _es->reinit();
      }
#else
    CPPUNIT_ASSERT_EQUAL(n_refinements, 0u);
#endif // LIBMESH_ENABLE_AMR

    System & system = _es->get_system("SimpleSystem");

    std::set<processor_id_type> sub_one_proc_ids, sub_two_proc_ids;
    for (auto & elem : _mesh->active_subdomain_elements_ptr_range(1))
      sub_one_proc_ids.insert(elem->processor_id());

    for (auto & elem : _mesh->active_subdomain_elements_ptr_range(2))
      sub_two_proc_ids.insert(elem->processor_id());


    // Everyone should be on the same processor if only 1 processor
    if( system.n_processors() == 1 )
      {
        CPPUNIT_ASSERT_EQUAL(1, static_cast<int>(sub_one_proc_ids.size()));
        CPPUNIT_ASSERT_EQUAL(1, static_cast<int>(sub_two_proc_ids.size()));
        CPPUNIT_ASSERT( sub_one_proc_ids == sub_two_proc_ids );
      }
    // Otherwise these sets should be disjoint
    else
      {
        // Make sure no subdomain one ids in the subdomain two ids
        for (auto & id : sub_one_proc_ids )
          CPPUNIT_ASSERT( sub_two_proc_ids.find(id) == sub_two_proc_ids.end() );

        // Vice-versa
        for (auto & id : sub_two_proc_ids )
          CPPUNIT_ASSERT( sub_one_proc_ids.find(id) == sub_one_proc_ids.end() );
      }
  }

};

// In this testing, we're relying on the presence of PETSc
#ifdef LIBMESH_HAVE_PETSC


// This testing class the algebraic ghosting of the
// OverlappingCouplingFunctor.
class OverlappingAlgebraicGhostingTest : public CppUnit::TestCase,
                                         public OverlappingTestBase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE( OverlappingAlgebraicGhostingTest );

  CPPUNIT_TEST( testGhostingCouplingMatrix );
  CPPUNIT_TEST( testGhostingNullCouplingMatrix );
#ifdef LIBMESH_ENABLE_AMR
  CPPUNIT_TEST( testGhostingNullCouplingMatrixUnifRef );
#endif

  CPPUNIT_TEST_SUITE_END();

public:

  void setUp()
  {}

  void tearDown()
  { this->clear(); }

  void testGhostingCouplingMatrix()
  {
    LOG_UNIT_TEST;
    this->run_ghosting_test(0, true);
  }

  void testGhostingNullCouplingMatrix()
  {
    LOG_UNIT_TEST;
    this->run_ghosting_test(0, false);
  }

  void testGhostingNullCouplingMatrixUnifRef()
  {
    LOG_UNIT_TEST;
    std::unique_ptr<CouplingMatrix> coupling_matrix;

    this->run_ghosting_test(2, false);
  }

private:

  void run_ghosting_test(const unsigned int n_refinements, bool build_coupling_matrix)
  {
    this->build_quad_mesh(n_refinements);
    this->init();

    std::unique_ptr<CouplingMatrix> coupling_matrix;
    if (build_coupling_matrix)
      this->setup_coupling_matrix(coupling_matrix);

    LinearImplicitSystem & system = _es->get_system<LinearImplicitSystem>("SimpleSystem");

    // If we don't add this coupling functor and properly recompute the
    // sparsity pattern, then PETSc will throw a malloc error when we
    // try to assemble into the global matrix
    OverlappingCouplingFunctor functor(system);
    functor.set_coupling_matrix(coupling_matrix);

    DofMap & dof_map = system.get_dof_map();
    dof_map.add_algebraic_ghosting_functor(functor);
    dof_map.reinit_send_list(system.get_mesh());

    // Update current local solution
    system.current_local_solution = libMesh::NumericVector<libMesh::Number>::build(system.comm());

    system.current_local_solution->init(system.n_dofs(), system.n_local_dofs(),
                                        dof_map.get_send_list(), false,
                                        libMesh::GHOSTED);

    system.solution->localize(*(system.current_local_solution),dof_map.get_send_list());

    std::unique_ptr<PointLocatorBase> point_locator = _mesh->sub_point_locator();

    const unsigned int u_var = system.variable_number("U");

    DenseMatrix<Number> K;

    FEMContext subdomain_one_context(system);
    FEMContext subdomain_two_context(system);

    // The use case on which this test is based, we only add terms to the residual
    // corresponding to the dofs in the second subdomain, but that have couplings
    // to dofs in the first subdomain.
    for (const auto & elem : _mesh->active_local_subdomain_elements_ptr_range(2))
      {
        // A little extra unit testing on the range iterator
        CPPUNIT_ASSERT_EQUAL(2, static_cast<int>(elem->subdomain_id()));

        subdomain_one_context.get_element_fe(u_var)->get_nothing(); // for this unit test
        const std::vector<libMesh::Point> & qpoints = subdomain_two_context.get_element_fe(u_var)->get_xyz();

        // Setup the context for the current element
        subdomain_two_context.pre_fe_reinit(system,elem);
        subdomain_two_context.elem_fe_reinit();

        std::set<subdomain_id_type> allowed_subdomains;
        allowed_subdomains.insert(1);

        // Now loop over the quadrature points and find the subdomain-one element that overlaps
        // with the current subdomain-two element and then initialize the FEMContext for the
        // subdomain-one element. If the algebraic ghosting has not been done properly,
        // this will error.
        for ( const auto & qp : qpoints )
          {
            const Elem * overlapping_elem = (*point_locator)( qp, &allowed_subdomains );
            CPPUNIT_ASSERT(overlapping_elem);

            // Setup the context for the overlapping element
            subdomain_one_context.pre_fe_reinit(system,overlapping_elem);
            subdomain_one_context.elem_fe_reinit();
          }
      }
  }

};


// This testing class now exercises testing the sparsity
// pattern augmented by the OverlappingCouplingFunctor
class OverlappingCouplingGhostingTest : public CppUnit::TestCase,
                                        public OverlappingTestBase
{
public:
  LIBMESH_CPPUNIT_TEST_SUITE( OverlappingCouplingGhostingTest );

  CPPUNIT_TEST( testSparsityCouplingMatrix );
  CPPUNIT_TEST( testSparsityNullCouplingMatrix );
#ifdef LIBMESH_ENABLE_AMR
  CPPUNIT_TEST( testSparsityNullCouplingMatrixUnifRef );
#endif

  CPPUNIT_TEST_SUITE_END();

public:

  void setUp()
  {}

  void tearDown()
  { this->clear(); }

  void testSparsityCouplingMatrix()
  {
    LOG_UNIT_TEST;
    this->run_sparsity_pattern_test(0, true);
  }

  void testSparsityNullCouplingMatrix()
  {
    LOG_UNIT_TEST;
    this->run_sparsity_pattern_test(0, false);
  }

  void testSparsityNullCouplingMatrixUnifRef()
  {
    LOG_UNIT_TEST;
    this->run_sparsity_pattern_test(1, false);
  }

private:

  void run_sparsity_pattern_test(const unsigned int n_refinements, bool build_coupling_matrix)
  {
    this->build_quad_mesh(n_refinements);
    this->init();

    std::unique_ptr<CouplingMatrix> coupling_matrix;
    if (build_coupling_matrix)
      this->setup_coupling_matrix(coupling_matrix);

    LinearImplicitSystem & system = _es->get_system<LinearImplicitSystem>("SimpleSystem");

    // If we don't add this coupling functor and properly recompute the
    // sparsity pattern, then PETSc will throw a malloc error when we
    // try to assemble into the global matrix
    OverlappingCouplingFunctor coupling_functor(system);
    coupling_functor.set_coupling_matrix(coupling_matrix);

    DofMap & dof_map = system.get_dof_map();
    dof_map.add_coupling_functor(coupling_functor);
    dof_map.reinit_send_list(system.get_mesh());

    // Update current local solution
    system.current_local_solution = libMesh::NumericVector<libMesh::Number>::build(system.comm());

    system.current_local_solution->init(system.n_dofs(), system.n_local_dofs(),
                                        dof_map.get_send_list(), false,
                                        libMesh::GHOSTED);

    system.solution->localize(*(system.current_local_solution),dof_map.get_send_list());

    // Now that we've added the coupling functor, we need
    // to recompute the sparsity
    dof_map.clear_sparsity();
    dof_map.compute_sparsity(system.get_mesh());

    // Now that we've recomputed the sparsity pattern, we need
    // to reinitialize the system matrix.
    SparseMatrix<Number> & matrix = system.get_system_matrix();
    libmesh_assert(dof_map.is_attached(matrix));
    matrix.init();

    std::unique_ptr<PointLocatorBase> point_locator = _mesh->sub_point_locator();

    const unsigned int u_var = system.variable_number("U");
    const unsigned int v_var = system.variable_number("V");

    DenseMatrix<Number> K12, K21;

    FEMContext subdomain_one_context(system);
    subdomain_one_context.get_element_fe(u_var)->get_nothing();
    subdomain_one_context.get_element_fe(v_var)->get_nothing();

    FEMContext subdomain_two_context(system);
    subdomain_two_context.get_element_fe(u_var)->get_xyz();
    subdomain_two_context.get_element_fe(v_var)->get_nothing();

    // Add normally coupled parts of the matrix
    for (const auto & elem : _mesh->active_local_subdomain_elements_ptr_range(1))
      {
        subdomain_one_context.pre_fe_reinit(system,elem);
        subdomain_one_context.elem_fe_reinit();

        std::vector<dof_id_type> & rows = subdomain_one_context.get_dof_indices();

        // Fill with ones in case PETSc ignores the zeros at some point
        std::fill( subdomain_one_context.get_elem_jacobian().get_values().begin(),
                   subdomain_one_context.get_elem_jacobian().get_values().end(),
                   1);

        // Insert the Jacobian for the dofs for this element
        matrix.add_matrix( subdomain_one_context.get_elem_jacobian(), rows );
      }

    for (const auto & elem : _mesh->active_local_subdomain_elements_ptr_range(2))
      {
        // A little extra unit testing on the range iterator
        CPPUNIT_ASSERT_EQUAL(2, static_cast<int>(elem->subdomain_id()));

        const std::vector<libMesh::Point> & qpoints = subdomain_two_context.get_element_fe(u_var)->get_xyz();

        // Setup the context for the current element
        subdomain_two_context.pre_fe_reinit(system,elem);
        subdomain_two_context.elem_fe_reinit();

        // We're only assembling rows for the dofs on subdomain 2 (U,L), so
        // the current element will have all those dof_indices.
        std::vector<dof_id_type> & rows = subdomain_two_context.get_dof_indices();

        std::fill( subdomain_two_context.get_elem_jacobian().get_values().begin(),
                   subdomain_two_context.get_elem_jacobian().get_values().end(),
                   1);

        // Insert the Jacobian for the normally coupled dofs for this element
        matrix.add_matrix( subdomain_two_context.get_elem_jacobian(), rows );

        std::set<subdomain_id_type> allowed_subdomains;
        allowed_subdomains.insert(1);

        // Now loop over the quadrature points and find the subdomain-one element that overlaps
        // with the current subdomain-two element and then add a local element matrix with
        // the coupling to the global matrix to try and trip any issues with sparsity pattern
        // construction
        for ( const auto & qp : qpoints )
          {
            const Elem * overlapping_elem = (*point_locator)( qp, &allowed_subdomains );
            CPPUNIT_ASSERT(overlapping_elem);

            // Setup the context for the overlapping element
            subdomain_one_context.pre_fe_reinit(system,overlapping_elem);
            subdomain_one_context.elem_fe_reinit();

            // We're only coupling to the "V" variable so only need those dof indices
            std::vector<dof_id_type> & v_indices = subdomain_one_context.get_dof_indices(v_var);
            std::vector<dof_id_type> columns(rows);
            columns.insert( columns.end(), v_indices.begin(), v_indices.end() );

            // This will also zero the matrix so we can just insert zeros for this test
            K21.resize( rows.size(), columns.size() );

            std::fill(K21.get_values().begin(), K21.get_values().end(), 1);

            // Now adding this local matrix to the global would trip a PETSc
            // malloc error if the sparsity pattern hasn't been correctly
            // built to include the overlapping coupling.
            matrix.add_matrix (K21, rows, columns);

            // Now add the other part of the overlapping coupling
            K12.resize(v_indices.size(), rows.size());
            std::fill(K12.get_values().begin(), K12.get_values().end(), 1);
            matrix.add_matrix(K12,v_indices,rows);
          }
      } // end element loop

    // We need to make sure to close the matrix for this test. There could still
    // be PETSc malloc errors tripped here if we didn't allocate the off-processor
    // part of the sparsity pattern correctly.
    matrix.close();
  }

};
#endif // LIBMESH_HAVE_PETSC


CPPUNIT_TEST_SUITE_REGISTRATION( OverlappingFunctorTest );

#ifdef LIBMESH_HAVE_PETSC
CPPUNIT_TEST_SUITE_REGISTRATION( OverlappingAlgebraicGhostingTest );
CPPUNIT_TEST_SUITE_REGISTRATION( OverlappingCouplingGhostingTest );
#endif
