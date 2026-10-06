// The libMesh Finite Element Library.
// Copyright (C) 2002-2026 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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

// <h1>Finite Volume Example 1 - Cell-Centered Advection-Diffusion on Polygonal Meshes</h1>
// \author Alexander Lindsay
// \date 2026
//
// This example solves the steady advection-diffusion equation
//
//   -div(D grad u) + div(v u) = f
//
// on the unit square with a cell-centered finite volume method on a
// mesh of C0Polygon elements. The solution is a CONSTANT MONOMIAL
// variable, so each cell carries a single degree of freedom,
// interpreted as the solution value at the cell centroid.
//
// Integrating the equation over a cell and applying the divergence
// theorem gives a balance of fluxes through the cell faces. Cell
// gradients, reconstructed by least squares from the face neighbors
// and the boundary data, enter these fluxes in two places:
//
// * The diffusive flux uses a two-point difference along the vector
//   joining adjacent cell centroids, plus a correction for the angle
//   between that vector and the face normal. On this mesh the two
//   are not parallel across the slanted hexagon faces, and without
//   the correction the scheme would not converge to the exact
//   solution.
//
// * The advective flux is upwinded, with the upwind cell value
//   linearly extrapolated to the face centroid, giving second-order
//   accuracy.
//
// A face flux therefore depends on the two cells sharing the face and
// on their face neighbors, so the matrix row of a cell couples it to
// two layers of side neighbors. The DofMap's DefaultCoupling functor
// is configured for two layers before the mesh is built, so that the
// required ghost elements are retained when a distributed mesh deletes
// remote elements.
//
// The example solves a manufactured-solution problem on a sequence of
// uniformly refined meshes and reports the centroid-sampled L2 error,
// which should decrease at second order.

// libMesh includes
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/dof_map.h"
#include "libmesh/default_coupling.h"
#include "libmesh/elem.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/getpot.h"
#include "libmesh/enum_solver_package.h"
#include "libmesh/int_range.h"
#include "libmesh/exodusII_io.h"

// C++ includes
#include <cmath>
#include <map>
#include <memory>
#include <vector>

using namespace libMesh;

namespace
{

// The manufactured solution, which also supplies the Dirichlet data
// on the whole boundary.
Real
exact_solution(const Point & p)
{
  return std::sin(libMesh::pi * p(0)) * std::cos(libMesh::pi * p(1));
}

// The source term -D lap(u) + v.grad(u) for the manufactured solution
// and a constant, hence divergence-free, velocity v.
Real
forcing(const Point & p, const Real diffusivity, const RealVectorValue & velocity)
{
  const Real pi = libMesh::pi;
  const Real x = p(0);
  const Real y = p(1);
  return 2 * pi * pi * diffusivity * exact_solution(p) +
         velocity(0) * pi * std::cos(pi * x) * std::cos(pi * y) -
         velocity(1) * pi * std::sin(pi * x) * std::sin(pi * y);
}

// An affine function of the cell solution values,
//   sum_i coefficients[i] * u_i + constant,
// whose value type T is Real for scalar quantities such as face fluxes
// and RealVectorValue for reconstructed gradients.
template <typename T>
struct AffineStencil
{
  std::map<dof_id_type, T> coefficients;
  T constant = T();
};

typedef AffineStencil<Real> ScalarStencil;
typedef AffineStencil<RealVectorValue> GradientStencil;

// Adds scale * (g . a) to out.
void
add_dot(ScalarStencil & out, const GradientStencil & g, const RealVectorValue & a, const Real scale)
{
  for (const auto & [dof, coef] : g.coefficients)
    out.coefficients[dof] += scale * (coef * a);
  out.constant += scale * (g.constant * a);
}

// Returns 0.5 * (a + b).
GradientStencil
average(const GradientStencil & a, const GradientStencil & b)
{
  GradientStencil avg;
  for (const auto & [dof, coef] : a.coefficients)
    avg.coefficients[dof] += 0.5 * coef;
  for (const auto & [dof, coef] : b.coefficients)
    avg.coefficients[dof] += 0.5 * coef;
  avg.constant = 0.5 * (a.constant + b.constant);
  return avg;
}

struct FaceGeometry
{
  Point centroid;
  Real area;
  Point normal;
};

// Computes the centroid, area, and outward unit normal of a planar
// element side. Polygons and polyhedra have no reference element, so
// the geometry is taken from the physical side directly.
FaceGeometry
face_geometry(const Elem & elem, const unsigned int side)
{
  const std::unique_ptr<const Elem> side_elem = elem.build_side_ptr(side);
  return {side_elem->true_centroid(), side_elem->volume(), elem.side_vertex_average_normal(side)};
}

// Returns the least-squares gradient of u in elem as an affine function
// of cell values. Each face contributes one sample: the neighbor's
// centroid value for interior faces and the Dirichlet value at the
// face centroid for boundary faces. The gradient g minimizes
//   sum_k ((x_k - x_elem) . g - (u_k - u_elem))^2,
// which reproduces linear functions exactly.
GradientStencil
gradient_stencil(const Elem & elem, const unsigned int sys_num)
{
  struct Sample
  {
    RealVectorValue offset;
    dof_id_type dof;
    Real value;
  };

  const unsigned int dim = elem.dim();
  const Point x_elem = elem.true_centroid();

  std::vector<Sample> samples;
  DenseMatrix<Real> normal_matrix(dim, dim);
  for (const auto s : elem.side_index_range())
    {
      Sample sample;
      if (const Elem * neighbor = elem.neighbor_ptr(s))
        sample = {neighbor->true_centroid() - x_elem, neighbor->dof_number(sys_num, 0, 0), 0};
      else
        {
          const FaceGeometry face = face_geometry(elem, s);
          sample = {face.centroid - x_elem, DofObject::invalid_id, exact_solution(face.centroid)};
        }

      for (const auto i : make_range(dim))
        for (const auto j : make_range(dim))
          normal_matrix(i, j) += sample.offset(i) * sample.offset(j);
      samples.push_back(sample);
    }

  // g = sum_k w_k (u_k - u_elem), with weights w_k solving
  // normal_matrix w_k = x_k - x_elem.
  GradientStencil g;
  const dof_id_type elem_dof = elem.dof_number(sys_num, 0, 0);
  for (const auto & sample : samples)
    {
      DenseVector<Real> rhs(dim), w;
      for (const auto i : make_range(dim))
        rhs(i) = sample.offset(i);
      normal_matrix.lu_solve(rhs, w);

      RealVectorValue weight;
      for (const auto i : make_range(dim))
        weight(i) = w(i);

      g.coefficients[elem_dof] -= weight;
      if (sample.dof == DofObject::invalid_id)
        g.constant += sample.value * weight;
      else
        g.coefficients[sample.dof] += weight;
    }

  return g;
}

} // anonymous namespace

// Assembles one matrix row per local cell from the balance
//   sum_faces (-D grad(u) + v u) . n |face| = f(x_elem) |elem|.
void
assemble_advection_diffusion(EquationSystems & es, const std::string & system_name)
{
  const MeshBase & mesh = es.get_mesh();
  auto & system = es.get_system<LinearImplicitSystem>(system_name);
  const unsigned int sys_num = system.number();

  const Real diffusivity = es.parameters.get<Real>("diffusivity");
  const RealVectorValue velocity = es.parameters.get<RealVectorValue>("velocity");

  for (const Elem * elem : mesh.active_local_element_ptr_range())
    {
      const dof_id_type elem_dof = elem->dof_number(sys_num, 0, 0);
      const Point x_elem = elem->true_centroid();
      const GradientStencil grad_elem = gradient_stencil(*elem, sys_num);

      // The net outward flux through the cell boundary
      ScalarStencil net_flux;

      for (const auto s : elem->side_index_range())
        {
          const FaceGeometry face = face_geometry(*elem, s);
          const Elem * neighbor = elem->neighbor_ptr(s);

          // On interior faces, the neighbor centroid and gradient; on
          // boundary faces, the face centroid, where the Dirichlet
          // value is known, and the gradient of this cell.
          const Point x_across = neighbor ? neighbor->true_centroid() : face.centroid;
          const GradientStencil grad_neighbor =
            neighbor ? gradient_stencil(*neighbor, sys_num) : GradientStencil();
          const GradientStencil grad_face = neighbor ? average(grad_elem, grad_neighbor) : grad_elem;

          // Diffusion: with d the offset from x_elem to x_across,
          //   n.grad(u) ~= (u_across - u_elem) / (d.n) + (n - d / (d.n)).grad(u)_face,
          // which is exact for linear u. The second term vanishes when d
          // is parallel to n.
          const RealVectorValue d = x_across - x_elem;
          const Real d_dot_n = d * face.normal;
          const Real diffusive_scale = -diffusivity * face.area;

          net_flux.coefficients[elem_dof] -= diffusive_scale / d_dot_n;
          if (neighbor)
            net_flux.coefficients[neighbor->dof_number(sys_num, 0, 0)] += diffusive_scale / d_dot_n;
          else
            net_flux.constant += diffusive_scale / d_dot_n * exact_solution(face.centroid);
          add_dot(net_flux, grad_face, face.normal - d / d_dot_n, diffusive_scale);

          // Advection: the upwind cell value, extrapolated to the face
          // centroid with that cell's gradient. Inflow boundary faces
          // take the Dirichlet value.
          const Real mass_flux = (velocity * face.normal) * face.area;
          if (mass_flux >= 0)
            {
              net_flux.coefficients[elem_dof] += mass_flux;
              add_dot(net_flux, grad_elem, face.centroid - x_elem, mass_flux);
            }
          else if (neighbor)
            {
              net_flux.coefficients[neighbor->dof_number(sys_num, 0, 0)] += mass_flux;
              add_dot(net_flux, grad_neighbor, face.centroid - x_across, mass_flux);
            }
          else
            net_flux.constant += mass_flux * exact_solution(face.centroid);
        }

      std::vector<dof_id_type> row_dofs = {elem_dof};
      std::vector<dof_id_type> column_dofs;
      DenseMatrix<Number> Ke(1, net_flux.coefficients.size());
      for (const auto & [dof, coef] : net_flux.coefficients)
        {
          Ke(0, column_dofs.size()) = coef;
          column_dofs.push_back(dof);
        }

      DenseVector<Number> Fe(1);
      Fe(0) = forcing(x_elem, diffusivity, velocity) * elem->volume() - net_flux.constant;

      system.matrix->add_matrix(Ke, row_dofs, column_dofs);
      system.rhs->add_vector(Fe, row_dofs);
    }
}

// Returns the L2 norm of the difference between the cell values and
// the exact solution sampled at cell centroids.
Real
centroid_l2_error(const System & system)
{
  const unsigned int sys_num = system.number();

  Real error_sq = 0;
  for (const Elem * elem : system.get_mesh().active_local_element_ptr_range())
    {
      const Real u_h = libmesh_real(system.current_solution(elem->dof_number(sys_num, 0, 0)));
      const Real diff = u_h - exact_solution(elem->true_centroid());
      error_sq += elem->volume() * diff * diff;
    }
  system.comm().sum(error_sq);

  return std::sqrt(error_sq);
}

int
main(int argc, char ** argv)
{
  // Initialize libMesh.
  LibMeshInit init(argc, argv);

  // This example requires a linear solver package.
  libmesh_example_requires(libMesh::default_solver_package() != INVALID_SOLVER_PACKAGE,
                           "--enable-petsc, --enable-trilinos, or --enable-eigen");

  // The mesh is two-dimensional.
  libmesh_example_requires(2 <= LIBMESH_DIM, "2D support");

  // Parse the input file, allowing the command line to override it.
  GetPot infile("finite_volume_ex1.in");
  infile.parse_command_line(argc, argv);

  const unsigned int coarse_grid_size = infile("coarse_grid_size", 4);
  const unsigned int n_refinements = infile("n_refinements", 3);
  const Real diffusivity = infile("diffusivity", 1.);
  const RealVectorValue velocity(infile("velocity_x", 1.), infile("velocity_y", 0.5));

  Real previous_error = 0;
  for (const auto level : make_range(n_refinements + 1))
    {
      const unsigned int grid_size = coarse_grid_size << level;

      Mesh mesh(init.comm());
      EquationSystems equation_systems(mesh);
      auto & system = equation_systems.add_system<LinearImplicitSystem>("AdvectionDiffusion");
      system.add_variable("u", CONSTANT, MONOMIAL);
      system.attach_assemble_function(assemble_advection_diffusion);

      // Couple each cell to two layers of side neighbors. The DofMap
      // registers its coupling functor with the mesh, so building the
      // mesh after this call also retains the corresponding ghost
      // elements on a distributed mesh.
      system.get_dof_map().default_coupling().set_n_levels(2);

      MeshTools::Generation::build_square(mesh, grid_size, grid_size, 0., 1., 0., 1., C0POLYGON);

      equation_systems.parameters.set<Real>("diffusivity") = diffusivity;
      equation_systems.parameters.set<RealVectorValue>("velocity") = velocity;
      equation_systems.init();

      system.solve();

      const Real error = centroid_l2_error(system);
      libMesh::out << "Grid size " << grid_size << ": " << mesh.n_active_elem()
                   << " cells, centroid L2 error " << error;
      if (level)
        libMesh::out << ", convergence rate " << std::log2(previous_error / error);
      libMesh::out << std::endl;
      previous_error = error;

#ifdef LIBMESH_HAVE_EXODUS_API
      if (level == n_refinements)
        ExodusII_IO(mesh).write_equation_systems("out.e", equation_systems);
#endif
    }

  return 0;
}
