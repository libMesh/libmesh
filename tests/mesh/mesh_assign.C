#include <libmesh/mesh.h>
#include <libmesh/distributed_mesh.h>
#include <libmesh/mesh_generation.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/boundary_info.h>
#include <libmesh/elem.h>
#include <libmesh/equation_systems.h>
#include <libmesh/zero_function.h>
#include <libmesh/dirichlet_boundaries.h>
#include <libmesh/dof_map.h>
#include <libmesh/parallel.h>

#include "test_comm.h"
#include "libmesh_cppunit.h"


using namespace libMesh;

class MeshAssignTest : public CppUnit::TestCase {
  /**
   * This test verifies the operations of the move assignment operator for MeshBase and
   * derived classes. Future extensions can include copy assignment, move and copy constructors.
   */
public:
  CPPUNIT_TEST_SUITE( MeshAssignTest );

/* Tests need a 2d mesh and Dirichlet boundary conditions */
#if LIBMESH_DIM > 1
# ifdef LIBMESH_ENABLE_DIRICHLET
  CPPUNIT_TEST( testMeshMoveAssign );
# endif
#endif

  CPPUNIT_TEST_SUITE_END();

protected:

public:
  void setUp()
  {
  }

  void tearDown()
  {
  }

  void testMeshMoveAssign()
  {

    // We will do 6 tests:
    // {Mesh, ReplicatedMesh ,DistributedMesh} x {Mesh 2 from Memory, From File}
    std::array<std::string, 3> mesh_types = {"Mesh", "replicated", "distributed"};
    std::array<std::string, 2> mesh_two_cases = {"from_memory", "from_file"};

    // Create two mesh ptrs.
    std::shared_ptr<UnstructuredMesh> mesh_one;
    std::shared_ptr<UnstructuredMesh> mesh_two;

    for(auto& mesh_type : mesh_types)
    {
      libMesh::out << "Running test with " << mesh_type << " mesh." << std::endl << std::endl;

      for(auto mesh_two_case_type : mesh_two_cases )
      {
       libMesh::out << "Running test with mesh two " << mesh_two_case_type << "." << std::endl << std::endl;

       /**
        * Build a 2d 2x2 square mesh (mesh_one) covering [0.0, 1.0] x [0.0, 1.0]
        * with linear Quad elements. A new mesh will later be moved to
        * this same object using MeshBase::operator = (&&).
        */

       if (mesh_type.compare("Mesh") == 0)
       {
        mesh_one = std::make_shared<Mesh>(*TestCommWorld);
       }
       else if (mesh_type.compare("replicated") == 0)
       {
        mesh_one = std::make_shared<ReplicatedMesh>(*TestCommWorld);
       }
       else if (mesh_type.compare("distributed") == 0)
       {
        mesh_one = std::make_shared<DistributedMesh>(*TestCommWorld);
       }
       else
       {
        libmesh_error_msg("Error: specified mesh_type not understood");
       }

       MeshTools::Generation::build_square(*mesh_one,
                                          2, 2,
                                          0., 1.,
                                          0., 1.,
                                          QUAD9);

       libMesh::out << "Info from mesh one." << std::endl << std::endl;
       mesh_one->print_info();

       // build_square adds boundary_ids 0,1,2,3 for the bottom, right,
       // top, and left sides, respectively.

       /** Create an ES based on the mesh, with a single variable.
        * We will set non-homogenous bc u = 1.0 on the bottom and right
        * sides of the mesh. After moving a new mesh into mesh_one, we will
        * verify that u is projected properly onto the new mesh_one topology.
        */
       EquationSystems es(*mesh_one);
       System & system = es.add_system<System> ("SimpleSystem");
       system.add_variable("u", FIRST);

       // Now set the Dirichlet bcs: u = 1.0, when boundary id = 0 and 1.
       // The function describing the boundary value
       ConstFunction<Number> one(1.0);

       // The boundary ids of the boundaries where the value is to be set.
       const boundary_id_type all_bdry_ids[4] = {0, 1, 2, 3};
       const std::set<boundary_id_type> bottom_right_bdry_ids(all_bdry_ids, all_bdry_ids + 1);

       // The variables for which we are defining the bc.
       std::vector<unsigned int> variables = {0};
       variables.push_back(0);

       DirichletBoundary dirichlet_bc(bottom_right_bdry_ids, variables, &one);

       system.get_dof_map().add_dirichlet_boundary(dirichlet_bc);

       es.init();

       if (mesh_type.compare("Mesh") == 0)
       {
        mesh_two = std::make_shared<Mesh>(*TestCommWorld);
       }
       else if (mesh_type.compare("replicated") == 0)
       {
        mesh_two = std::make_shared<ReplicatedMesh>(*TestCommWorld);
       }
       else if (mesh_type.compare("distributed") == 0)
       {
        mesh_two = std::make_shared<DistributedMesh>(*TestCommWorld);
       }
       else
       {
        libmesh_error_msg("Error: specified mesh_type not understood");
       }

       if(mesh_two_case_type.compare("from_memory") == 0)
       {
         /**
          * Build another 2d 2x2 square mesh (mesh_two) covering [0.0, 1.0] x [0.0, 1.0]
          * with linear Quad elements. We will later refine this mesh and replace
          * mesh_one with it.
          */

          MeshTools::Generation::build_square(*mesh_two,
                                              2, 2,
                                              0., 1.,
                                              0., 1.,
                                              QUAD4);

          // Uniformly refine mesh_two
          std::unique_ptr<MeshRefinement> mesh_refinement(libmesh_make_unique<MeshRefinement>(*mesh_two));

          mesh_refinement->uniformly_refine(1);
          libMesh::out << "Info from mesh two, from memory." << std::endl << std::endl;
          mesh_two->print_info();

          // Move mesh_two into mesh_one
          system.get_mesh().clear();
          system.get_mesh().assign(std::move(*mesh_two));
          mesh_two->clear();
          libMesh::out << "Info from mesh one, after moving." << std::endl << std::endl;
          system.get_mesh().print_info();

          // Assert that the moved into mesh has the right number of elements.
          CPPUNIT_ASSERT_EQUAL(system.get_mesh().n_elem(), dof_id_type(20));

          // Reinit the dofs and other system related data structures
          // based on its mesh
          system.get_equation_systems().reinit_mesh();

          // Reinit
          system.get_equation_systems().reinit();

          // Write out the mesh xda, this currently fails with DistributedMesh,
          // running on more than 8 processors, when the mesh to be moved in is
          // read from a file.
          libMesh::out << "Write out the mesh, read from memory." << std::endl << std::endl;
          system.get_mesh().write("mesh.out.memory.moved.xda");

          system.get_mesh().clear();
        }
        else if (mesh_two_case_type.compare("from_file") == 0)
        {
          mesh_two->read("meshes/mesh_assign_test_mesh.xda");
          libMesh::out << "Info from mesh two, from file." << std::endl << std::endl;
          mesh_two->print_info();

          // Read in the mesh at this time instant and reinit to project solutions on to the new mesh
          system.get_mesh().clear();
          system.get_mesh().assign(std::move(*mesh_two));

          mesh_two->clear();
          libMesh::out << "Info from mesh one, after moving." << std::endl << std::endl;
          system.get_mesh().print_info();

          // Assert that the moved into mesh has the right number of elements.
          CPPUNIT_ASSERT_EQUAL(system.get_mesh().n_elem(), dof_id_type(42));

          system.get_equation_systems().reinit_mesh();
          system.get_equation_systems().reinit();

          libMesh::out << "Write out the mesh, read from file." << std::endl << std::endl;
          system.get_mesh().write("mesh.out.file.moved.xda");

          system.get_mesh().clear();
        }
         else
        {
         libmesh_error_msg("Error: invalid mesh two case type.");
        }

      }// End loop over mesh two memory vs file case

    }// End loop over type of mesh to be tested

  } // End definition of MeshAssignTest

}; // End definition of class MeshAssignTest

CPPUNIT_TEST_SUITE_REGISTRATION( MeshAssignTest );
