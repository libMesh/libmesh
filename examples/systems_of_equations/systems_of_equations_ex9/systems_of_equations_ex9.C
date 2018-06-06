// The libMesh Finite Element Library.
// Copyright (C) 2002-2018 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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



// <h1> Systems Example 9 - Linear elasticity problem with periodic constraints </h1>
//
// In this example we illustrate periodic constraints with a linear elasticity example.
// We consider a sector of a circular domain, and hence we must impose an "azimuthal"
// periodic condition, which entails a dof-transformation between the two periodic
// boundaries to account for the rotation of the coordinate system.
//
// The baseline code is from systems_of_equations_ex6, and we edit to use a different mesh
// and to impose the periodic boundary condition.


// C++ include files that we need
#include <iostream>
#include <algorithm>
#include <math.h>

// libMesh includes
#include "libmesh/libmesh_config.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dof_map.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dense_subvector.h"
#include "libmesh/perf_log.h"
#include "libmesh/elem.h"
#include "libmesh/boundary_info.h"
#include "libmesh/zero_function.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/getpot.h"
#include "libmesh/solver_configuration.h"
#include "libmesh/petsc_linear_solver.h"
#include "libmesh/petsc_macro.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/enum_solver_package.h"

// Bring in everything from the libMesh namespace
using namespace libMesh;

// Here we define the azimuthal periodic boundary condition class.
// This class assumes that (u,v,w) are variables 0, 1, and 2 in
// the System, but this could be made more general if needed.
class AzimuthalPeriodicBoundary : public PeriodicBoundaryBase
{
public:
  /**
   * Constructor
   */
  AzimuthalPeriodicBoundary(
    Point center,
    Point axis,
    Real angle)
  :
    PeriodicBoundaryBase(),
    _center(center),
    _axis(axis),
    _theta(angle)
  {
    set_variable(0);
    set_variable(1);
    set_variable(2);
    set_transformation_matrix(get_rotation_matrix());
  }

  /**
   * Copy constructor, with option for the copy to represent an inverse transformation.
   */
  AzimuthalPeriodicBoundary(const AzimuthalPeriodicBoundary & o, TransformationType t = FORWARD)
    :
    PeriodicBoundaryBase(o),
    _center(o._center),
    _axis(o._axis),
    _theta(o._theta)
  {
    if (t == INVERSE)
      {
        std::swap(myboundary, pairedboundary);
        _theta *= -1.0;
      }

    set_variable(0);
    set_variable(1);
    set_variable(2);
    set_transformation_matrix(get_rotation_matrix());
  }

  /**
   * Destructor
   */
  virtual ~AzimuthalPeriodicBoundary() {}

  /**
   * Get the rotation matrix for this transformation.
   */
  DenseMatrix<Real> get_rotation_matrix() const
  {
    // Formula for rotation matrix about an axis is given on wikipedia:
    // en.wikipedia.org/wiki/Rotation_matrix
    // We rotate by angle theta about the axis defined by u, which is a
    // unit vector in the direction of _axis.
    Point u = _axis.unit();
    Real u_x = u(0);
    Real u_y = u(1);
    Real u_z = u(2);
    DenseMatrix<Real> R(3,3);
    R(0,0) = cos(_theta) + u_x*u_x*(1.0 - cos(_theta));
    R(0,1) = u_x*u_y*(1.0 - cos(_theta)) - u_z*sin(_theta);
    R(0,2) = u_x*u_z*(1.0 - cos(_theta)) + u_y*sin(_theta);
    R(1,0) = u_y*u_x*(1.0 - cos(_theta)) + u_z*sin(_theta);
    R(1,1) = cos(_theta) + u_y*u_y*(1.0 - cos(_theta));
    R(1,2) = u_y*u_z*(1.0 - cos(_theta)) - u_x*sin(_theta);
    R(2,0) = u_z*u_x*(1.0 - cos(_theta)) - u_y*sin(_theta);
    R(2,1) = u_z*u_y*(1.0 - cos(_theta)) + u_x*sin(_theta);
    R(2,2) = cos(_theta) + u_z*u_z*(1.0 - cos(_theta));

    return R;
  }

  /**
   * This function should be overridden by derived classes to
   * define how one finds corresponding nodes on the periodic
   * boundary pair.
   */
  virtual Point get_corresponding_pos(const Point & pt) const override
  {
    DenseVector<Real> translated_pt(3);
    for(unsigned int i=0; i<3; i++)
    {
      translated_pt(i) = pt(i) - _center(i);
    }

    // Note that since _theta defines the angle from "paired boundary" to
    // "my boundary", and we want the inverse of that here, we must use
    // vector_mult_transpose below.
    DenseVector<Real> rotated_pt;
    get_transformation_matrix().vector_mult_transpose(rotated_pt, translated_pt);

    Point corresponding_pos;
    for(unsigned int i=0; i<3; i++)
    {
      corresponding_pos(i) = rotated_pt(i) + _center(i);
    }
    return corresponding_pos;
  }

  /**
   * If we want the DofMap to be able to make copies of references and
   * store them in the underlying map, this class must be clone'able,
   * i.e. have a kind of virtual construction mechanism.
   */
  virtual std::unique_ptr<PeriodicBoundaryBase> clone(TransformationType t = FORWARD) const override
  {
    return libmesh_make_unique<AzimuthalPeriodicBoundary>(*this, t);
  }

private:

  // Define the properties needed for the azimuthal periodic boundary.
  // Note that _theta specifies the angle from "paired boundary" to
  // "my boundary".
  Point _center;
  Point _axis;
  Real _theta;
};

class LinearElasticity : public System::Assembly
{
private:
  EquationSystems & es;

public:

  LinearElasticity (EquationSystems & es_in) :
    es(es_in)
  {}

  /**
   * Kronecker delta function.
   */
  Real kronecker_delta(unsigned int i,
                       unsigned int j)
  {
    return i == j ? 1. : 0.;
  }

  /**
   * Evaluate the fourth order tensor (C_ijkl) that relates stress to strain.
   */
  Real elasticity_tensor(unsigned int i,
                         unsigned int j,
                         unsigned int k,
                         unsigned int l)
  {
    // Hard code material parameters for the sake of simplicity
    const Real poisson_ratio = 0.3;
    const Real young_modulus = 1.;

    // Define the Lame constants
    const Real lambda_1 = (young_modulus*poisson_ratio)/((1.+poisson_ratio)*(1.-2.*poisson_ratio));
    const Real lambda_2 = young_modulus/(2.*(1.+poisson_ratio));

    return lambda_1 * kronecker_delta(i, j) * kronecker_delta(k, l) +
      lambda_2 * (kronecker_delta(i, k) * kronecker_delta(j, l) + kronecker_delta(i, l) * kronecker_delta(j, k));
  }

  /**
   * Assemble the system matrix and right-hand side vector.
   */
  void assemble()
  {
    const MeshBase & mesh = es.get_mesh();

    const unsigned int dim = mesh.mesh_dimension();

    LinearImplicitSystem & system = es.get_system<LinearImplicitSystem>("Elasticity");

    const unsigned int u_var = system.variable_number ("u");

    const DofMap & dof_map = system.get_dof_map();
    FEType fe_type = dof_map.variable_type(u_var);
    std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
    QGauss qrule (dim, fe_type.default_quadrature_order());
    fe->attach_quadrature_rule (&qrule);

    std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type));
    QGauss qface(dim-1, fe_type.default_quadrature_order());
    fe_face->attach_quadrature_rule (&qface);

    const std::vector<Real> & JxW = fe->get_JxW();
    const std::vector<std::vector<Real>> & phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();

    DenseMatrix<Number> Ke;
    DenseSubMatrix<Number> Ke_var[3][3] =
      {
        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
        {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}
      };

    DenseVector<Number> Fe;

    DenseSubVector<Number> Fe_var[3] =
      {DenseSubVector<Number>(Fe),
       DenseSubVector<Number>(Fe),
       DenseSubVector<Number>(Fe)};

    std::vector<dof_id_type> dof_indices;
    std::vector<std::vector<dof_id_type>> dof_indices_var(3);

    for (const auto & elem : mesh.active_local_element_ptr_range())
      {
        dof_map.dof_indices (elem, dof_indices);
        for (unsigned int var=0; var<3; var++)
          dof_map.dof_indices (elem, dof_indices_var[var], var);

        const unsigned int n_dofs   = dof_indices.size();
        const unsigned int n_var_dofs = dof_indices_var[0].size();

        fe->reinit (elem);

        Ke.resize (n_dofs, n_dofs);
        for (unsigned int var_i=0; var_i<3; var_i++)
          for (unsigned int var_j=0; var_j<3; var_j++)
            Ke_var[var_i][var_j].reposition (var_i*n_var_dofs, var_j*n_var_dofs, n_var_dofs, n_var_dofs);

        Fe.resize (n_dofs);
        for (unsigned int var=0; var<3; var++)
          Fe_var[var].reposition (var*n_var_dofs, n_var_dofs);

        for (unsigned int qp=0; qp<qrule.n_points(); qp++)
          {
            // assemble \int_Omega C_ijkl u_k,l v_i,j \dx
            for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
              for (unsigned int dof_j=0; dof_j<n_var_dofs; dof_j++)
                for (unsigned int i=0; i<3; i++)
                  for (unsigned int j=0; j<3; j++)
                    for (unsigned int k=0; k<3; k++)
                      for (unsigned int l=0; l<3; l++)
                        Ke_var[i][k](dof_i,dof_j) +=
                          JxW[qp] * elasticity_tensor(i,j,k,l) * dphi[dof_j][qp](l) * dphi[dof_i][qp](j);

            // assemble \int_Omega f_i v_i \dx
            DenseVector<Number> f_vec(3);
            if(elem->subdomain_id() == 101)
              {
                f_vec(0) = 1.;
                f_vec(1) = 1.;
                f_vec(2) = 0.;
              }
            else if(elem->subdomain_id() == 1)
              {
                f_vec(0) = 0.36603;
                f_vec(1) = 1.36603;
                f_vec(2) = 0.;
              }
            for (unsigned int dof_i=0; dof_i<n_var_dofs; dof_i++)
              for (unsigned int i=0; i<3; i++)
                Fe_var[i](dof_i) += JxW[qp] * (f_vec(i) * phi[dof_i][qp]);
          }

        dof_map.constrain_element_matrix_and_vector (Ke, Fe, dof_indices);

        system.matrix->add_matrix (Ke, dof_indices);
        system.rhs->add_vector    (Fe, dof_indices);
      }
  }
};


// Begin the main program.
int main (int argc, char ** argv)
{
  // Initialize libMesh and any dependent libraries
  LibMeshInit init (argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // Initialize the cantilever mesh
  const unsigned int dim = 3;

  // Make sure libMesh was compiled for 3D
  libmesh_example_requires(dim == LIBMESH_DIM, "3D support");

  // Create a 3D mesh distributed across the default MPI communicator.
  Mesh mesh(init.comm());

  // Create an equation systems object.
  EquationSystems equation_systems (mesh);

  // Declare the system and its variables.
  // Create a system named "Elasticity"
  LinearImplicitSystem & system =
    equation_systems.add_system<LinearImplicitSystem> ("Elasticity");

  // Add two azimuthal periodic boundaries on two adjacent domains.
  // We do this to show that the periodic boundary condition that
  // we impose leads to a continuous solution across adjacent domains.
  //
  // We add the periodic boundaries *before* reading the Mesh, so
  // that periodic neighbors will be retained when a DistributedMesh
  // is distributed.
  //
  // The angle specified below defines the mapping
  // from "pairedboundary" to "myboundary".
  {
    Point center(0., 0., 0.);
    Point axis(0., 0., 1.);
    Real angle = 2*libMesh::pi/12.0;
    AzimuthalPeriodicBoundary periodic_bc(center, axis, angle);
    periodic_bc.myboundary = 301;
    periodic_bc.pairedboundary = 302;
    system.get_dof_map().add_periodic_boundary(periodic_bc);
  }
  {
    Point center(0., 0., 0.);
    Point axis(0., 0., 1.);
    Real angle = 2*libMesh::pi/12.0;
    AzimuthalPeriodicBoundary periodic_bc(center, axis, angle);
    periodic_bc.myboundary = 401;
    periodic_bc.pairedboundary = 402;
    system.get_dof_map().add_periodic_boundary(periodic_bc);
  }

  mesh.read("systems_of_equations_ex9.exo");

  // Print information about the mesh to the screen.
  mesh.print_info();

  // Add three displacement variables, u and v, to the system
  unsigned int u_var = system.add_variable("u", FIRST, LAGRANGE);
  unsigned int v_var = system.add_variable("v", FIRST, LAGRANGE);
  unsigned int w_var = system.add_variable("w", FIRST, LAGRANGE);

  LinearElasticity le(equation_systems);
  system.attach_assemble_object(le);

  std::vector<unsigned int> variables;
  variables.push_back(u_var);
  variables.push_back(v_var);
  variables.push_back(w_var);
  ZeroFunction<> zf;
  std::set<boundary_id_type> clamped_boundary_ids;
  clamped_boundary_ids.insert(300);
  clamped_boundary_ids.insert(400);
  DirichletBoundary clamped_bc(clamped_boundary_ids, variables, zf);
  system.get_dof_map().add_dirichlet_boundary(clamped_bc);

  // Initialize the data structures for the equation system.
  equation_systems.init();

  // Print information about the system to the screen.
  equation_systems.print_info();

  // Solve the system
  system.solve();

  // Plot the solution
#ifdef LIBMESH_HAVE_EXODUS_API

  ExodusII_IO (mesh).write_equation_systems("solution.exo",
                                            equation_systems);

#endif // #ifdef LIBMESH_HAVE_EXODUS_API

  // All done.
  return 0;
}
