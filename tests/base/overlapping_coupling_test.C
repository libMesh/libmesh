// Ignore unused parameter warnings coming from cppunit headers
#include <libmesh/ignore_warnings.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/TestCase.h>
#include <libmesh/restore_warnings.h>

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
    CPPUNIT_ASSERT_EQUAL(2, (int)_mesh.n_subdomains());

    // We're assuming a specific set of subdomain ids
    std::set<subdomain_id_type> ids;
    _mesh.subdomain_ids(ids);
    CPPUNIT_ASSERT( ids.find(1) != ids.end() );
    CPPUNIT_ASSERT( ids.find(2) != ids.end() );

    _point_locator->enable_out_of_mesh_mode();
  }

  virtual ~OverlappingCouplingFunctor(){};

  void set_coupling_matrix (std::unique_ptr<CouplingMatrix> & coupling_matrix)
  { _coupling_matrix = std::move(coupling_matrix); }

  virtual void operator()
  (const MeshBase::const_element_iterator & range_begin,
   const MeshBase::const_element_iterator & range_end,
   processor_id_type p,
   std::unordered_map<const Elem *,const CouplingMatrix*> & coupled_elements) override
  {

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
            const Elem * overlapping_elem = (*_point_locator)( qp, &allowed_subdomains );

            if( overlapping_elem && overlapping_elem->processor_id() != p )
              coupled_elements.insert( std::make_pair(overlapping_elem,_coupling_matrix.get()) );
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
    return libmesh_make_unique<OverlappingTestPartitioner>(*this);
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

  void assign_proc_id_subdomain( const MeshBase & mesh,
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
    _mesh.reset(new ReplicatedMesh(*TestCommWorld));

    _mesh->set_mesh_dimension(2);

    _mesh->add_point( Point(0.0,0.0),0 );
    _mesh->add_point( Point(1.0,0.0),1 );
    _mesh->add_point( Point(1.0,1.0),2 );
    _mesh->add_point( Point(0.0,1.0),3 );

    {
      Elem* elem = _mesh->add_elem( new Quad4 );
      elem->set_id(0);
      elem->subdomain_id() = 1;

      for (unsigned int n=0; n<4; n++)
        elem->set_node(n) = _mesh->node_ptr(n);
    }

    _mesh->add_point( Point(1.0,2.0),4 );
    _mesh->add_point( Point(0.0,2.0),5 );

    {
      Elem* elem = _mesh->add_elem( new Quad4 );
      elem->set_id(1);
      elem->subdomain_id() = 1;

      elem->set_node(0) = _mesh->node_ptr(3);
      elem->set_node(1) = _mesh->node_ptr(2);
      elem->set_node(2) = _mesh->node_ptr(4);
      elem->set_node(3) = _mesh->node_ptr(5);
    }

    _mesh->add_point( Point(0.0,0.0),6 );
    _mesh->add_point( Point(1.0,0.0),7 );
    _mesh->add_point( Point(1.0,2.0),8 );
    _mesh->add_point( Point(0.0,2.0),9 );

    {
      Elem* elem = _mesh->add_elem( new Quad4 );
      elem->set_id(2);
      elem->subdomain_id() = 2;

      elem->set_node(0) = _mesh->node_ptr(6);
      elem->set_node(1) = _mesh->node_ptr(7);
      elem->set_node(2) = _mesh->node_ptr(8);
      elem->set_node(3) = _mesh->node_ptr(9);
    }

    _mesh->partitioner() = std::unique_ptr<Partitioner>(new OverlappingTestPartitioner);

    _mesh->prepare_for_use();

    if (n_refinements > 0)
      {
        MeshRefinement refine(*_mesh);
        refine.uniformly_refine(n_refinements);
      }
  }

  void init(MeshBase & mesh)
  {
    _es.reset( new EquationSystems(*_mesh) );
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


    coupling.reset( new CouplingMatrix(system.n_vars()) );

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

  CPPUNIT_TEST_SUITE( OverlappingFunctorTest );

  CPPUNIT_TEST( checkCouplingFunctorQuad );
  CPPUNIT_TEST( checkCouplingFunctorQuadUnifRef );
  CPPUNIT_TEST( checkCouplingFunctorTri );
  CPPUNIT_TEST( checkCouplingFunctorTriUnifRef );
  CPPUNIT_TEST( checkOverlappingPartitioner );
  CPPUNIT_TEST( checkOverlappingPartitionerUnifRef );

  CPPUNIT_TEST_SUITE_END();

public:

  public:
  void setUp()
  {
    this->build_quad_mesh();
    this->init(*_mesh);
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
  // are kosher.
  void run_coupling_functor_test(unsigned int n_refinements)
  {
    if( n_refinements > 0 )
      {
        MeshRefinement refine(*_mesh);
        refine.uniformly_refine(n_refinements);
        _es->reinit();
      }

    System & system = _es->get_system("SimpleSystem");

    OverlappingCouplingFunctor coupling_functor(system);

    std::unordered_map<const Elem *,const CouplingMatrix*> subdomain_one_couplings;
    std::unordered_map<const Elem *,const CouplingMatrix*> subdomain_two_couplings;

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

    CPPUNIT_ASSERT_EQUAL( n_elems_subdomain_two, (dof_id_type) subdomain_one_couplings.size() );
    CPPUNIT_ASSERT_EQUAL( 2*n_elems_subdomain_two, (dof_id_type) subdomain_two_couplings.size() );
  }

  void run_partitioner_test(unsigned int n_refinements)
  {
    if( n_refinements > 0 )
      {
        MeshRefinement refine(*_mesh);
        refine.uniformly_refine(n_refinements);
        _es->reinit();
      }

    System & system = _es->get_system("SimpleSystem");

    std::set<processor_id_type> sub_one_proc_ids, sub_two_proc_ids;
    for (auto & elem : _mesh->active_subdomain_elements_ptr_range(1))
      sub_one_proc_ids.insert(elem->processor_id());

    for (auto & elem : _mesh->active_subdomain_elements_ptr_range(2))
      sub_two_proc_ids.insert(elem->processor_id());


    // Everyone should be on the same processor if only 1 processor
    if( system.n_processors() == 1 )
      {
        CPPUNIT_ASSERT_EQUAL( 1, (int)sub_one_proc_ids.size() );
        CPPUNIT_ASSERT_EQUAL( 1, (int)sub_two_proc_ids.size() );
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

CPPUNIT_TEST_SUITE_REGISTRATION( OverlappingFunctorTest );
