// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

#include "libmesh/coupling_matrix.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/gmv_io.h"
#include "libmesh/libmesh.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"


using namespace libMesh;


void assemble(EquationSystems & es,
              const std::string & system_name);





#ifdef LIBMESH_ENABLE_AMR
int main (int argc, char ** argv)
{
  LibMeshInit init(argc, argv);

  if (argc < 4)
    libMesh::out << "Usage: ./prog -d DIM filename" << std::endl;

  // Variables to get us started
  const unsigned int dim = atoi(argv[2]);

  std::string meshname  (argv[3]);

  // declare a mesh...
  Mesh mesh(init.comm(), dim);

  // Read a mesh
  mesh.read(meshname);

  GMVIO(mesh).write ("out_0.gmv");

  mesh.elem_ref(0).set_refinement_flag (Elem::REFINE);

  MeshRefinement mesh_refinement (mesh);

  mesh_refinement.refine_and_coarsen_elements ();
  mesh_refinement.uniformly_refine (2);

  mesh.print_info();


  // Set up the equation system(s)
  EquationSystems es (mesh);

  LinearImplicitSystem & primary =
    es.add_system<LinearImplicitSystem>("primary");

  primary.add_variable ("U", FIRST);
  primary.add_variable ("V", FIRST);

  primary.get_dof_map()._dof_coupling->resize(2);
  (*primary.get_dof_map()._dof_coupling)(0,0) = 1;
  (*primary.get_dof_map()._dof_coupling)(1,1) = 1;

  primary.attach_assemble_function(assemble);

  es.init ();

  es.print_info ();
  primary.get_dof_map().print_dof_constraints ();

  // call the solver.
  primary.solve ();

  GMVIO(mesh).write_equation_systems ("out_1.gmv",
                                      es);



  // Refine uniformly
  mesh_refinement.uniformly_refine (1);
  es.reinit ();

  // Write out the projected solution
  GMVIO(mesh).write_equation_systems ("out_2.gmv",
                                      es);

  // Solve again. Output the refined solution
  primary.solve ();
  GMVIO(mesh).write_equation_systems ("out_3.gmv",
                                      es);

  return 0;
}
#else
int main (int, char **)
{
  return 1;
}
#endif // ENABLE_AMR






void assemble(EquationSystems & es,
              const std::string & libmesh_dbg_var(system_name))
{
  libmesh_assert_equal_to (system_name, "primary");

  const MeshBase & mesh   = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  // Also use a 3x3x3 quadrature rule (3D).  Then tell the FE
  // about the geometry of the problem and the quadrature rule
  FEType fe_type (FIRST);

  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, FIFTH);

  fe->attach_quadrature_rule (&qrule);

  UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_type));
  QGauss qface(dim-1, FIFTH);

  fe_face->attach_quadrature_rule(&qface);

  LinearImplicitSystem & system =
    es.get_system<LinearImplicitSystem>("primary");


  // These are references to cell-specific data
  const std::vector<Real> & JxW_face                   = fe_face->get_JxW();
  const std::vector<Real> & JxW                        = fe->get_JxW();
  const std::vector<Point> & q_point                   = fe->get_xyz();
  const std::vector<std::vector<Real> > & phi          = fe->get_phi();
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  std::vector<unsigned int> dof_indices_U;
  std::vector<unsigned int> dof_indices_V;
  const DofMap & dof_map = system.get_dof_map();

  DenseMatrix<Number> Kuu;
  DenseMatrix<Number> Kvv;
  DenseVector<Number> Fu;
  DenseVector<Number> Fv;

  Real vol=0., area=0.;

  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();

  for (; el != end_el; ++el)
    {
      const Elem * elem = *el;

      // recompute the element-specific data for the current element
      fe->reinit (elem);


      //fe->print_info();

      dof_map.dof_indices(elem, dof_indices_U, 0);
      dof_map.dof_indices(elem, dof_indices_V, 1);

      // zero the element matrix and vector
      Kuu.resize (phi.size(),
                  phi.size());

      Kvv.resize (phi.size(),
                  phi.size());

      Fu.resize (phi.size());
      Fv.resize (phi.size());

      // standard stuff...  like in code 1.
      for (unsigned int gp=0; gp<qrule.n_points(); gp++)
        {
          for (std::size_t i=0; i<phi.size(); ++i)
            {
              // this is tricky.  ig is the _global_ dof index corresponding
              // to the _global_ vertex number elem->node_id(i).  Note that
              // in general these numbers will not be the same (except for
              // the case of one unknown per node on one subdomain) so
              // we need to go through the dof_map

              const Real f = q_point[gp]*q_point[gp];
              //    const Real f = (q_point[gp](0) +
              //    q_point[gp](1) +
              //    q_point[gp](2));

              // add jac*weight*f*phi to the RHS in position ig

              Fu(i) += JxW[gp]*f*phi[i][gp];
              Fv(i) += JxW[gp]*f*phi[i][gp];

              for (std::size_t j=0; j<phi.size(); ++j)
                {

                  Kuu(i,j) += JxW[gp]*((phi[i][gp])*(phi[j][gp]));

                  Kvv(i,j) += JxW[gp]*((phi[i][gp])*(phi[j][gp]) +
                                       1.*((dphi[i][gp])*(dphi[j][gp])));
                };
            };
          vol += JxW[gp];
        };


      // You can't compute "area" (perimeter) if you are in 2D
      if (dim == 3)
        {
          for (unsigned int side=0; side<elem->n_sides(); side++)
            if (elem->neighbor_ptr(side) == libmesh_nullptr)
              {
                fe_face->reinit (elem, side);

                //fe_face->print_info();

                for (std::size_t gp=0; gp<JxW_face.size(); gp++)
                  area += JxW_face[gp];
              }
        }

      // Constrain the DOF indices.
      dof_map.constrain_element_matrix_and_vector(Kuu, Fu, dof_indices_U);
      dof_map.constrain_element_matrix_and_vector(Kvv, Fv, dof_indices_V);


      system.rhs->add_vector(Fu, dof_indices_U);
      system.rhs->add_vector(Fv, dof_indices_V);

      system.matrix->add_matrix(Kuu, dof_indices_U);
      system.matrix->add_matrix(Kvv, dof_indices_V);
    }

  libMesh::out << "Vol="  << vol << std::endl;

  if (dim == 3)
    libMesh::out << "Area=" << area << std::endl;
}
