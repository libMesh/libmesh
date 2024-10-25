// libMesh includes
#include <libmesh/dof_map.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/linear_implicit_system.h>
#include <libmesh/non_manifold_coupling.h>
#include <libmesh/partitioner.h>
#include <libmesh/replicated_mesh.h>
#include <libmesh/sides_to_elem_map.h>

// Unit test includes
#include "test_comm.h"
#include "libmesh_cppunit.h"

// C++ includes
#include <memory>

using namespace libMesh;

// We use a custom partioner for this test to help ensure we're
// testing the cases we want. In particular, we're going to ensure
// that elements attached to the "shared" side are on different
// processors from each other so we know ghosting will be required for
// algebraic and coupling ghosting between the two subdomains.
class NonManifoldTestPartitioner : public Partitioner
{
public:

  NonManifoldTestPartitioner () = default;
  NonManifoldTestPartitioner (const NonManifoldTestPartitioner &) = default;
  NonManifoldTestPartitioner (NonManifoldTestPartitioner &&) = default;
  NonManifoldTestPartitioner & operator= (const NonManifoldTestPartitioner &) = default;
  NonManifoldTestPartitioner & operator= (NonManifoldTestPartitioner &&) = default;
  virtual ~NonManifoldTestPartitioner() = default;

  /**
   * \returns A copy of this partitioner wrapped in a smart pointer.
   */
  virtual std::unique_ptr<Partitioner> clone () const override
  {
    return std::make_unique<NonManifoldTestPartitioner>(*this);
  }

protected:

  /**
   * Partition the \p MeshBase onto \p n processors.
   */
  virtual void _do_partition (MeshBase & mesh,
                              const unsigned int n) override
  {
    // If we're on one partition, then everyone gets to be on that partition
    if (n == 1)
      this->single_partition_range (mesh.active_elements_begin(), mesh.active_elements_end());
    else
      {
        // Assign a round-robin processor id to each Elem attached to
        // the "non-manifold" edge in this Mesh. In this case, the
        // test Meshes were created so that they contain a single
        // non-manifold edge, so we can stop partitioning once we have
        // found that one edge.
        libmesh_assert_greater (n, 0);

        MeshTools::SidesToElemMap stem = MeshTools::SidesToElemMap::build(mesh);

        auto success = [&]() -> bool
        {
          for (const auto & elem : mesh.element_ptr_range())
            for (auto s : elem->side_index_range())
              {
                const auto [side_neighbors_begin, side_neighbors_end] =
                  stem.get_connected_elems(elem, s);

                if (std::distance(side_neighbors_begin, side_neighbors_end) == 5)
                  {
                    for (auto [e, it] = std::make_tuple(0u, side_neighbors_begin);
                         it != side_neighbors_end; ++e, ++it)
                      {
                        // If there are more Elems attached to the
                        // non-manifold edge than processors, just wrap
                        // around back to 0. Setting the processor id is a
                        // bit convoluted because the SidesToElemMap only
                        // has const pointers, but since we have a
                        // non-const reference to the mesh, we can work
                        // around that.
                        Elem * elem = mesh.elem_ptr((*it)->id());
                        elem->processor_id() = e % n;

                        // Debugging
                        // libMesh::out << "Partitioning Elem " << elem->id()
                        //              << " onto proc " << elem->processor_id()
                        //              << std::endl;
                      }

                    // Stop once we have partitioned elements attached to
                    // the non-manifold Side with 5 neighbors. This should
                    // effectively leave all other Elems assigned to
                    // processor 0.
                    return true;
                  }
              }

          // If we made it here, then we didn't find a Side shared by
          // 5 elements and we'll throw an error later
          return false;
        }();

        libmesh_error_msg_if(!success, "Did not find expected non-manifold edge.");
      }
  }
};


// Common functionality for all the different tests we'll be running
class NonManifoldCouplingTestBase
{
protected:

  std::unique_ptr<MeshBase> _mesh;

  std::unique_ptr<EquationSystems> _es;

  void read_mesh(const std::string & mesh_filename)
  {
    // We are making assumptions in various places about the presence
    // of the elements on the current processor so we're restricting to
    // ReplicatedMesh for now.
    _mesh = std::make_unique<ReplicatedMesh>(*TestCommWorld);
    _mesh->read(mesh_filename);

    // Use our custom partitioner that assigns non-manifold edge
    // neighbors to different processors.
    _mesh->partitioner() = std::make_unique<NonManifoldTestPartitioner>();

    _mesh->prepare_for_use();

    // Debugging: print processor ids determined by partitioner
    // for (const auto & elem : _mesh->element_ptr_range())
    //   libMesh::out << "Elem " << elem->id() << " is on processor " << elem->processor_id() << std::endl;
  }

  void init_es()
  {
    _es = std::make_unique<EquationSystems>(*_mesh);

    // Add System with a single linear Lagrange variable on it.
    LinearImplicitSystem & sys = _es->add_system<LinearImplicitSystem> ("SimpleSystem");
    sys.add_variable("u", FIRST, LAGRANGE);

    // Attach ghosting functor to the System. We use the std::shared_ptr interface
    // so that the DofMap takes ownership of the object. Note: if you comment out
    // this line, then the test _can_ fail when run in parallel, but it depends on
    // the number of processors used. For small numbers of procs, sometimes we can
    // get "lucky" and wind up with all the non-manifold neighbor DOFs either local
    // to or ghosted on all procs where they are needed. In practice, this test should
    // be run on 3 or more procs to make it fail when the NonManifoldGhostingFunctor
    // is not used.
    auto ghosting_functor = std::make_shared<NonManifoldGhostingFunctor>(*_mesh);
    sys.get_dof_map().add_coupling_functor(ghosting_functor, /*to_mesh=*/true);

    // Initialize the DofMap, etc. for all Systems
    _es->init();
  }
};

// This class defines the actual unit tests
class NonManifoldGhostingFunctorTest : public CppUnit::TestCase,
                                       public NonManifoldCouplingTestBase
{
public:

  LIBMESH_CPPUNIT_TEST_SUITE( NonManifoldGhostingFunctorTest );

  // These tests all require Exodus
#if defined(LIBMESH_HAVE_EXODUS_API)
  CPPUNIT_TEST( verifySendListEntries0 );
  CPPUNIT_TEST( verifySendListEntries1 );
  CPPUNIT_TEST( verifySendListEntries2 );
  CPPUNIT_TEST( verifySendListEntries3 );
#endif

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {
  }

  void tearDown()
  {
  }

  void verifySendListEntries0() { this->verify_send_list_entries_helper("meshes/non_manifold_junction0.exo"); }
  void verifySendListEntries1() { this->verify_send_list_entries_helper("meshes/non_manifold_junction1.exo"); }
  void verifySendListEntries2() { this->verify_send_list_entries_helper("meshes/non_manifold_junction2.exo"); }
  void verifySendListEntries3() { this->verify_send_list_entries_helper("meshes/non_manifold_junction3.exo"); }

  // We call this helper function (which can take an argument) from
  // the CPPUNIT_TEST macro functions, (which I don't think can take an argument?)
  void verify_send_list_entries_helper(const std::string & mesh_filename)
  {
    // Call base class implementations
    this->read_mesh(mesh_filename);
    this->init_es();

    // Get reference to the System, corresponding DofMap, and send_list for this processor
    System & system = _es->get_system("SimpleSystem");
    DofMap & dof_map = system.get_dof_map();
    const std::vector<dof_id_type> & send_list = dof_map.get_send_list();

    MeshTools::SidesToElemMap stem = MeshTools::SidesToElemMap::build(*_mesh);

    // Use the SidesToElemMap to get a range of Elem pointers attached
    // to the non-manifold edge that we intentionally partitioned onto
    // different processors.
    MeshTools::SidesToElemMap::ElemIter beg, end;
    auto side_neighbors_found = [&]() -> bool
    {
      for (const auto & elem : _mesh->element_ptr_range())
        for (auto s : elem->side_index_range())
          {
            auto range = stem.get_connected_elems(elem, s);
            if (std::distance(range.first, range.second) == 5)
              {
                beg = range.first;
                end = range.second;
                return true;
              }
          }

      // If we made it here, then we didn't find a Side shared by
      // 5 elements and we'll throw an error later
      return false;
    }();

    // Require that neighbors were found
    CPPUNIT_ASSERT(side_neighbors_found);

    // For every pair of elements (e,f) on this edge owned by processors (proc_e, proc_f) respectively, check that:
    // 1.) The DOFs of e are either local to, or in the send-list of, proc_f
    // 2.) The DOFs of f are either local to, or in the send-list of, proc_e
    for (auto it_e = beg; it_e != end; ++it_e)
      for (auto it_f = std::next(it_e); it_f != end; ++it_f)
        {
          // Lambda to be used for error checking
          auto check_dofs = [&](const Elem * elem)
          {
            std::vector<dof_id_type> dof_indices;
            dof_map.dof_indices(elem, dof_indices);

            for (const auto & dof : dof_indices)
            {
              bool is_local = (dof >= dof_map.first_dof()) && (dof < dof_map.end_dof());
              bool is_in_send_list = (Utility::binary_find(send_list.begin(), send_list.end(), dof) != send_list.end());
              CPPUNIT_ASSERT(is_local || is_in_send_list);
            }
          };

          const Elem * elem_e = *it_e;
          const Elem * elem_f = *it_f;

          if (_mesh->comm().rank() == elem_e->processor_id())
            check_dofs(elem_f);

          if (_mesh->comm().rank() == elem_f->processor_id())
            check_dofs(elem_e);
        }
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( NonManifoldGhostingFunctorTest );
