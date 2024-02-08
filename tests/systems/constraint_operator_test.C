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

namespace {

// A simple u"+u=f
void assemble_matrix_and_rhs(EquationSystems& es,
                             const std::string&)
{
  const MeshBase& mesh = es.get_mesh();
  LinearImplicitSystem& system = es.get_system<LinearImplicitSystem>("test");
  const DofMap& dof_map = system.get_dof_map();
  SparseMatrix<Number> & matrix = system.get_system_matrix();

  FEMContext c(system);

  const unsigned short dim = mesh.mesh_dimension();
  FEBase * fe = nullptr;
  c.get_element_fe(0, fe, dim);

  const std::vector<Real> & JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> & phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

  std::vector<dof_id_type> & dof_indices = c.get_dof_indices();

  // The subvectors and submatrices we need to fill:
  DenseMatrix<Number> & K = c.get_elem_jacobian();
  DenseVector<Number> & F = c.get_elem_residual();

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      const Elem* elem = *el;

      if(elem->type() == NODEELEM)
        continue;

      c.pre_fe_reinit(system, elem);
      c.elem_fe_reinit();

      unsigned int n_dofs = dof_indices.size();
      unsigned int n_qpoints = c.get_element_qrule().n_points();

      for (unsigned int qp=0; qp != n_qpoints; qp++)
        {
          const Gradient grad_u = c.interior_gradient(0, qp);
          const Number u = c.interior_value(0, qp);

          Number forcing = 1;

          for (unsigned int i=0; i != n_dofs; i++)
            {
              F(i) += JxW[qp] * (grad_u * dphi[i][qp] + (forcing-u) * phi[i][qp]);
              for (unsigned int j=0; j != n_dofs; ++j)
                K(i,j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp] - phi[i][qp]*phi[j][qp]);
            }
        }

      dof_map.constrain_element_matrix_and_vector (K, F, dof_indices);

      matrix.add_matrix (K, dof_indices);
      system.rhs->add_vector (F, dof_indices);
    }

  system.rhs->close();
  matrix.close();
}
}



class ConstraintOperatorTest : public CppUnit::TestCase {
public:
  LIBMESH_CPPUNIT_TEST_SUITE( ConstraintOperatorTest );

  CPPUNIT_TEST( test1DCoarseningOperator );

  CPPUNIT_TEST_SUITE_END();

public:
  void setUp()
  {}

  void tearDown()
  {}


  void test1DCoarseningOperator ()
  {
    LOG_UNIT_TEST;

    Mesh mesh(*TestCommWorld);
    EquationSystems es(mesh);

    ExplicitSystem &sys =
      es.add_system<LinearImplicitSystem> ("test");

    sys.add_variable("u", FIRST);
    sys.attach_assemble_function (assemble_matrix_and_rhs);

    // Make sure numbering matches our constraint operator matrix
    mesh.allow_renumbering(false);

    MeshTools::Generation::build_line (mesh,
                                       10,
                                       0., 1.,
                                       EDGE2);

    auto matrix = SparseMatrix<Real>::build (mesh.comm());

    // Create a projection matrix from 10 elements to 5.  We might as
    // well just make it on processor 0.
    processor_id_type my_rank = mesh.comm().rank();
    matrix->init(11, 6,
                 (my_rank ? 0 : 11),
                 (my_rank ? 0 : 6),
                 (my_rank ? 0 : 6),
                 0);

    if (!my_rank)
      {
        for (auto i : make_range(5))
          {
            matrix->set(i*2, i, 1);
            matrix->set(i*2+1, i, 0.5);
            matrix->set(i*2+1, i+1, 0.5);
          }
        matrix->set(10, 5, 1);
      }
    matrix->close();

    mesh.copy_constraint_rows(*matrix);

    es.init ();

    sys.solve();

    const DofMap & dof_map = sys.get_dof_map();

    CPPUNIT_ASSERT_EQUAL(dof_id_type(5), dof_map.n_constrained_dofs());
  }
};

CPPUNIT_TEST_SUITE_REGISTRATION( ConstraintOperatorTest );
