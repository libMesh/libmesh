// The libMesh Finite Element Library.
// Copyright (C) 2002-2025 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

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


// <h1>Miscellaneous Example 14 - Hydrogen Atom Using Infinite Elements With Imaginary Frequency</h1>
// \author Hubert Weissmann
// \date 2017
//
// This example uses second order infinite elements to solve the Schroedinger equation to obtain
// the wave function of the Hydrogen atom.
// To speed up the calculation, the SlepcSolver is configured to use the shift-and-invert-transform for solving
// the resulting generalised eigensystem.
//
// The result is printed in exodus-format and along the x-axis in a free format.

#include <iostream>
#include <fstream>
// libMesh include files.
#include "libmesh/getpot.h" // for input-argument parsing
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_tetgen_interface.h"
#include "libmesh/elem.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/eigen_system.h"
#include "libmesh/equation_systems.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dof_map.h"
#include "libmesh/fe_interface.h"
#include "libmesh/explicit_system.h"
// for infinite elements:
#include "libmesh/inf_fe.h"
#include "libmesh/inf_elem_builder.h"
#include "libmesh/fe_interface.h"
#include "libmesh/fe_compute_data.h"
// for finding element for point
#include "libmesh/point_locator_tree.h"

// for the SlepcSolverConfiguration
#include "libmesh/petsc_solver_exception.h"
#include "libmesh/solver_configuration.h"
#include "libmesh/slepc_eigen_solver.h"
#include "libmesh/enum_eigen_solver_type.h"

#ifdef LIBMESH_HAVE_SLEPC
EXTERN_C_FOR_SLEPC_BEGIN
# include <slepceps.h>
EXTERN_C_FOR_SLEPC_END
#endif // LIBMESH_HAVE_SLEPC

// Bring in everything from the libMesh namespace
using namespace libMesh;

#ifdef LIBMESH_HAVE_SLEPC
/**
 * Defines an \p enum for spectral tronsformations
 * applied before solving the (generalised) eigenproblem
 */
enum SpectralTransform {SHIFT=0,
                        SINVERT,
                        CAYLEY,

                        INVALID_ST
};

/**
 * A class that interfaces the \p SolverConfiguration to add the SLEPC option SetST.
 * In case of dense eigenvalues, the application of SINVERT or CAYLEY is highly beneficial
 * see <a href='http://slepc.upv.es/documentation/'>the SLEPC-manual</a> for details.
 */
class SlepcSolverConfiguration : public libMesh::SolverConfiguration
{
public:

  /**
   * Constructor: get a reference to the \p SlepcEigenSolver variable to be able to manipulate it
   */
  SlepcSolverConfiguration( libMesh::SlepcEigenSolver<libMesh::Number> &slepc_eigen_solver):
    _slepc_solver(slepc_eigen_solver),
    _st(INVALID_ST)
  {}

  /**
   * empty destructor
   */
  ~SlepcSolverConfiguration() {}

  virtual void configure_solver() override;

  /**
   * This is the main functionality of this little class.
   */
  void SetST(SpectralTransform st)
  { _st=st;}

private:

  /**
   *The linear solver object that we are configuring
   */
  libMesh::SlepcEigenSolver<libMesh::Number>& _slepc_solver;
  SpectralTransform _st;

};
#endif // LIBMESH_HAVE_SLEPC

//prototypes of functions needed to set-up the system:
// The assemble-function is similar to miscellaneous_ex1
void assemble_SchroedingerEquation(EquationSystems & , const std::string &);

// print the solution along the x-axis.
void line_print(EquationSystems& es, std::string output, std::string SysName);


int main (int argc, char** argv)
{
  // Initialize libMesh and the dependent libraries.
  LibMeshInit init (argc, argv);

  // Check for proper usage first.
  libmesh_error_msg_if(argc < 2, "\nUsage: " << argv[0] << " <input-filename>");

  // Tell the user what we are doing.
  std::cout << "Running " << argv[0];

  for (int i=1; i<argc; i++)
    std::cout << " " << argv[i];
  std::cout << std::endl << std::endl;

  // Skip SLEPc examples on a non-SLEPc libMesh build
#ifndef LIBMESH_HAVE_SLEPC
  libmesh_example_requires(false, "--enable-slepc");
#else
#ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS
  libmesh_example_requires(false, "--enable-ifem");
#else

#ifndef LIBMESH_USE_COMPLEX_NUMBERS
  libmesh_example_requires(false, "--enable-complex");
#endif

#ifdef LIBMESH_DEFAULT_SINGLE_PRECISION
  // SLEPc currently gives us a nasty crash with Real==float
  libmesh_example_requires(false, "--disable-singleprecision");
#endif

  // parse the input file to getpot and obtain an object holding all options
  GetPot cl(argv[1]);

  // Skip this example if libMesh was compiled with <3 dimensions.
  // INFINITE ELEMENTS ARE IMPLEMENTED ONLY FOR 3 DIMENSIONS AT THE MOMENT.
  libmesh_example_requires(3 <= LIBMESH_DIM, "3D support");

  // read the options from the input-file.
  const unsigned int nev = cl("nev",1);
  Real E = cl("Energy", -.5);
  Real r=cl("radius", 4.);
  unsigned int maxiter=cl("maxiter", 700);

  int dim = 3;

  // creation of an empty mesh-object
  Mesh mesh(init.comm(), dim);

  // fill the meshes with a spherical grid of type HEX27 with radius r
  MeshTools::Generation::build_sphere (mesh, r, 2, HEX27, 2, true);

  // The infinite elements are attached to all elements that build the outer surface of the FEM-region.
  InfElemBuilder builder(mesh);
  builder.build_inf_elem(true);

  // Reassign subdomain_id() of all infinite elements.
  // and their neighbours. This makes finding the surface
  // between these elemests much easier.
  for (auto & elem : mesh.element_ptr_range())
    if (elem->infinite())
      {
        elem->subdomain_id() = 1;
        // the base elements are always the 0-th neighbor.
        elem->neighbor_ptr(0)->subdomain_id()=2;
      }

  // find the neighbours; for correct linking the two areas
  mesh.find_neighbors();

  // Create an equation systems object
  EquationSystems eq_sys (mesh);

  // This is the only system added here.
  EigenSystem & eig_sys = eq_sys.add_system<EigenSystem> ("EigenSE");

  //set the complete type of the variable
  FEType fe_type(SECOND, LAGRANGE, FOURTH, JACOBI_20_00, CARTESIAN);

  // Name the variable of interest 'phi' and approximate it as \p fe_type.
  eig_sys.add_variable("phi", fe_type);

  // assign the ground state energy. For infinite elements, the quality of this guess is crucial
  // otherwise the long-range limit will be wrong. If the 'current frequency' is too much off,
  // the solution will have no physical meaning at all.
  eq_sys.parameters.set<Number>("gsE")=E;

  //save the mesh-size in a parameter to scale the region for printing the solution later accordingly.
  eq_sys.parameters.set<Real>("radius")    = r;

  //set number of eigen values ( \p nev) and number of
  // basis vectors \p ncv for the solution.
  //Note that ncv >= nev must hold and ncv >= 2*nev is recommended.
  eq_sys.parameters.set<unsigned int>("eigenpairs")    = nev;

  eq_sys.parameters.set<unsigned int>("basis vectors") = nev*3+4;

  // Set the solver tolerance and the maximum number of iterations.
  eq_sys.parameters.set<Real> ("linear solver tolerance") = pow(TOLERANCE, 5./3.);
  eq_sys.parameters.set<unsigned int>("linear solver maximum iterations") = maxiter;

  // set numerical parameters for SLEPC on how to solve the system.
#if SLEPC_VERSION_LESS_THAN(2,3,2)
  eig_sys.get_eigen_solver().set_eigensolver_type(ARNOLDI);
#else
  eig_sys.get_eigen_solver().set_eigensolver_type(KRYLOVSCHUR);
#endif

  eig_sys.get_eigen_solver().set_position_of_spectrum(E);

  //fetch the solver-object used internally to be able to manipulate it using the self-written class
  // to set the transformation
  SlepcEigenSolver<Number> * solver =
    cast_ptr<SlepcEigenSolver<Number>* >( &eig_sys.get_eigen_solver() );

  // setup of our class @SlepcSolverConfiguration
  SlepcSolverConfiguration ConfigSolver(*solver);

  // set the spectral transformation: the default (SHIFT) will not converge
  // so we apply some more elaborate scheme.
  ConfigSolver.SetST(SINVERT);
  solver ->set_solver_configuration(ConfigSolver);

  // attach the name of the function that assembles the matrix equation:
  eig_sys.attach_assemble_function (assemble_SchroedingerEquation);

  // important to set the system to be generalised nonhermitian eigen problem.
  // By default it is HEP and so _matrix_B is not available.
  // In the infinite element scheme, the Hamiltonian is not hermitian because left and right
  // basis are not the same!
  eig_sys.set_eigenproblem_type(GNHEP);

  // Initialize the data structures for the equation system.
  eq_sys.init();

  // Solve system. This function calls the assemble-functions.
  eig_sys.solve();

  // get number of actually converged eigenpairs. It can be larger or smaller than \p nev.
  unsigned int nconv = eig_sys.get_n_converged();

  out << "Number of converged eigenpairs: " << nconv << "\n";

  // Write the eigen vector to file and the eigenvalues to libMesh::out.
  for(unsigned int i=0; i<nconv; i++)
    {
      std::pair<Real,Real> eigpair = eig_sys.get_eigenpair(i);

      // get the complete eigenvalue:
      Complex eigenvalue(eigpair.first, eigpair.second);
      out<<"        "<<eigenvalue<<std::endl;

      // Here, one needs to formally distinguish between bound states (negative energy)
      // and free states (positive energy) because bound states have an imaginary momentum.
      // It is important to ensure that the imaginary part is positive; otherwise we get an exp.
      // rise instead of exponential decay when computing the shape functions in the infinite element region.
      if (eigpair.first > 0.)
        // positive eigen energy
        eq_sys.parameters.set<Complex>("momentum")=sqrt(eigenvalue*2.);
      else
        {
          // for negative energies: ensure that imaginary contribution is positive!
          Complex k= sqrt(2.*eigenvalue);
          if (std::imag(k)>0)
            eq_sys.parameters.set<Complex>("momentum")=k;
          else
            eq_sys.parameters.set<Complex>("momentum")=std::conj(k);
        }

#ifdef LIBMESH_HAVE_EXODUS_API
      // set the name of the Exodus-file
      std::ostringstream eigenvector_output_name;
      eigenvector_output_name << "U" << "-" << i << "_inf.e";
      ExodusII_IO (mesh).write_equation_systems(eigenvector_output_name.str(), eq_sys);
#endif

      // set the base-name for the free-format file.
      // This file can be viewed e.g. with gnuplot using
      // p 're_infini_0.txt' w l
      std::ostringstream file;
      file << "infini_" << i << ".txt";
      // print the solution along the x-coordinate
      line_print(eq_sys, file.str(), "EigenSE");
    }

  // All done.
#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
#endif // LIBMESH_HAVE_SLEPC
  return 0;
}

/**
 * In this function, we assemble the actual system matrices. The initial equation is the eigen system \f$ H \Psi = \epsilon
 * \Psi \f$
 * where H is the Hamilton operator and Psi the eigen vector with corresponding eigen value \f$\epsilon\f$.
 * In the FEM-scheme, this becomes the generalised eigen problem  \f$ H \Psi = \epsilon M \Psi \f$ where M is the element mass matrix.
 */
void assemble_SchroedingerEquation(EquationSystems &es, const std::string &system_name)
{
#if defined(LIBMESH_ENABLE_INFINITE_ELEMENTS) && defined(LIBMESH_HAVE_SLEPC)
  // It is a good idea to make sure we are assembling
  // the proper system.
  libmesh_assert_equal_to (system_name, "EigenSE");

  // Get a constant reference to the mesh object.
  const MeshBase& mesh = es.get_mesh();

  // The dimension that we are running.
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to our system.
  EigenSystem & eigen_system = es.get_system<EigenSystem> (system_name);

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = eigen_system.get_dof_map().variable_type(0);

  // A reference to the system matrices
  SparseMatrix<Number>&  matrix_A = eigen_system.get_matrix_A();
  SparseMatrix<Number>&  matrix_B = eigen_system.get_matrix_B();

  // A  Gauss quadrature rule for numerical integration.
  // Use the default quadrature order.
  QGauss qrule (dim, fe_type.default_quadrature_order());

  // Build a Finite Element object of the specified type.  Since the
  // \p FEBase::build() member dynamically creates memory we will
  // store the object as an \p std::unique_ptr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
  std::unique_ptr<FEBase> inf_fe (FEBase::build_InfFE(dim, fe_type));

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (&qrule);
  inf_fe->attach_quadrature_rule (&qrule);

  Number E=es.parameters.get<Number>("gsE");

  /*
   * this is used as a parameter in the matrix assembly which is used to resemble the exponential
   * term in the infinite region.
   * For bound states, i*k is real (and negative), for free states i*k is imaginary.
   * The latter case is what is usually assumed for infinite elements.
   */
  Number ik=-sqrt(-2.*E);

  // set parameters used for the infinite elements:
  es.parameters.set<Real>("speed")=137.0359991;
  es.parameters.set<Number>("current frequency")=es.parameters.get<Real>("speed")*sqrt(2.*E)/(2*pi);

  // A reference to the \p DofMap object for this system.  The \p DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.
  const DofMap& dof_map = eigen_system.get_dof_map();

  // The element mass matrix M and Hamiltonian H.
  DenseMatrix<Number> M;
  DenseMatrix<Number> H;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  In case users
  // later modify this program to include refinement, we will
  // be safe and will only consider the active elements;
  // hence we use a variant of the \p active_elem_iterator.
  for (const auto & elem : mesh.active_local_element_ptr_range())
    {
      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      // unifyging finite and infinite elements
      FEBase * cfe = nullptr;

      if (elem->infinite())
        {
          // We have an infinite element.  Let \p cfe point
          // to our \p InfFE object.  This is handled through
          // an std::unique_ptr.  Through the \p std::unique_ptr::get() we "borrow"
          // the pointer, while the \p  std::unique_ptr \p inf_fe is
          // still in charge of memory management.
          cfe = inf_fe.get();
        }
      else
        {
          cfe = fe.get();
        }

      // The element Jacobian * quadrature weight at each integration point.
      const std::vector<Real> & JxW = cfe->get_JxWxdecay_sq();

      // The element shape functions evaluated at the quadrature points.
      const std::vector<std::vector<Real>> & phi = cfe->get_phi_over_decayxR();
      const std::vector<std::vector<RealGradient>> & dphi = cfe->get_dphi_over_decayxR();
      const std::vector<Point>& q_point = cfe->get_xyz();
      // get extra data needed for infinite elements
      const std::vector<RealGradient>& dphase = cfe->get_dphase();
      const std::vector<Real> &         weight  = cfe->get_Sobolev_weightxR_sq();
      const std::vector<RealGradient> & dweight = cfe->get_Sobolev_dweightxR_sq();

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape functions
      // (phi, dphi) for the current element.
      cfe->reinit (elem);

      // Zero the element matrices and rhs before
      // summing them.  We use the resize member here because
      // the number of degrees of freedom might have changed from
      // the last element.  Note that this will be the case if the
      // element type is different (i.e. the last element was a
      // triangle, now we are on a quadrilateral).
      M.resize (dof_indices.size(), dof_indices.size());
      H.resize (dof_indices.size(), dof_indices.size());

      // Now loop over the quadrature points.  This handles
      // the numeric integration.
      //For infinite elements, the number of quadrature points is asked and than looped over; works for finite elements as well.
      unsigned int max_qp = cfe->n_quadrature_points();
      for (unsigned int qp=0; qp<max_qp; qp++)
        {

          // compute the Coulomb potential
          Number potval=-1./(q_point[qp]).norm();

          // Now, get number of shape functions that are nonzero at this point::
          const unsigned int n_sf =
            FEInterface::n_dofs(cfe->get_fe_type(), elem);
          // loop over them:
          for (unsigned int i=0; i<n_sf; i++)
            {
              for (unsigned int j=0; j<n_sf; j++)
                {

                  // This is just the of derivatives of shape functions i and j.
                  // For finite elements, dweight==0 and dphase==0, thus
                  // temp=dphi[i][qp]*dphi[j][qp].
                  Number temp = (dweight[qp]*phi[i][qp]+weight[qp]*(dphi[i][qp]-ik*dphase[qp]*phi[i][qp]))
                                *(dphi[j][qp]+ik*dphase[qp]*phi[j][qp]);

                  // assemble the Hamiltonian: H=1/2 Nabla^2 + V
                  H(i,j) += JxW[qp]*(0.5*temp + potval*weight[qp]*phi[i][qp]*phi[j][qp]);

                  // assemble the mass matrix:
                  M(i,j) += JxW[qp]*weight[qp]*phi[i][qp]*phi[j][qp];
                }
            }
        }
      // On an unrefined mesh, constrain_element_matrix does
      // nothing.  If this assembly function is ever repurposed to
      // run on a refined mesh, getting the hanging node constraints
      // right will be important.  Note that, even with
      // asymmetric_constraint_rows = false, the constrained dof
      // diagonals still exist in the matrix, with diagonal entries
      // that are there to ensure non-singular matrices for linear
      // solves but which would generate positive non-physical
      // eigenvalues for eigensolves.
      dof_map.constrain_element_matrix(M, dof_indices, false);
      dof_map.constrain_element_matrix(H, dof_indices, false);

      // Finally, simply add the element contribution to the
      // overall matrices.
      matrix_A.add_matrix (H, dof_indices);
      matrix_B.add_matrix (M, dof_indices);

    } // end of element loop

   // After the setup of the main part, we add the contributions from the surface,
   // arising due to the partial integration.
   // Keeping these terms follows the suggestion of Bettess (P. Bettess "Infinite Elements", Penshaw Pres 1992)
   subdomain_id_type sd;
   // For easier syntax, we loop over infinite elements and their neighbours separately.
   // loop over INFINITE ELEMENTS
   sd=1;{
      QGauss qrule2 (dim-1, fe_type.default_quadrature_order());
      auto face_fe = FEBase::build_InfFE(dim, fe_type);
      face_fe->attach_quadrature_rule (&qrule2);

      for (const auto & elem : mesh.active_local_subdomain_elements_ptr_range(sd))
      {
         //dof_map.dof_indices (elem, dof_indices_lm,lm_num);
         dof_map.dof_indices (elem, dof_indices);

         /*
          * Here we use an inconsistent setup: the Jacobian is weighted while all other quantities are
          * not weighted.
          * But on the base face all weights are one; thus, it is only a formal problem.
          */

         // Having the correct face, we can start initializing all quantities:
         const std::vector<Real>& JxW = face_fe->get_JxWxdecay_sq();
         const std::vector<Point>& q_point = face_fe->get_xyz();
         const std::vector<std::vector<RealGradient> >& dphi = face_fe->get_dphi();
         const std::vector<std::vector<Real> >& phi = face_fe->get_phi();
         const std::vector<Point>& normal = face_fe->get_normals();
         const std::vector<RealGradient>& dphase = face_fe->get_dphase();
         const std::vector<Real>& weight = face_fe->get_Sobolev_weight(); // in publication called D

         // the base side always has number 0
         const unsigned int side=0;
         face_fe->reinit (elem, side);
         const unsigned int n_sf =
           FEInterface::n_dofs(face_fe->get_fe_type(), elem);

         H.resize (dof_indices.size(), dof_indices.size());

         for (unsigned int pts=0; pts<q_point.size(); pts++){
            for (unsigned int i=0; i<n_sf; i++){
               // A(i,j)=A(column, row) =A(bra,ket)
               for (unsigned int j=0; j<n_sf; j++){
                  //WARNING: The normal() of infinite elements points inwards! Thus, the sign is flipped.
                  libmesh_assert(normal[pts]*q_point[pts]/q_point[pts].norm() > 0);
                  H(i,j)-=.5*JxW[pts]*weight[pts]*phi[i][pts]*normal[pts]*(dphi[j][pts]+ik*dphase[pts]*phi[j][pts]);
               }
            }
         }
        matrix_A.add_matrix (H, dof_indices);
      } // end of element loop
   }
   // loop over NEIGHBOURS OF INFINITE ELEMENTS
   sd=2; {
      QGauss qrule2 (dim-1, fe_type.default_quadrature_order());
      auto face_fe = FEBase::build(dim, fe_type);
      face_fe->attach_quadrature_rule (&qrule2);

      for (const auto & elem : mesh.active_local_subdomain_elements_ptr_range(sd))
      {
         dof_map.dof_indices (elem, dof_indices);

#ifdef DEBUG
         unsigned int num_neighbors=0;

         for (unsigned int i=0; i<elem->n_neighbors(); i++){
            if (elem->neighbor_ptr(i)->infinite()){
               num_neighbors++;
            }
         }
         libmesh_assert_greater(num_neighbors,0);
         libmesh_assert_greater(elem->n_neighbors(), num_neighbors);
#endif

         // In contrast to the infinite elements, here it is not as easy to find the correct side.
         // Also, a single finite element can have multiple infinite neighbours.
         for(unsigned int prev_neighbor=0; prev_neighbor<elem->n_neighbors(); prev_neighbor++){

            // Having the correct face, we can start initializing all quantities:
            const std::vector<Real>& JxW =face_fe->get_JxW();
            const std::vector<Point>& q_point = face_fe->get_xyz();
            const std::vector<std::vector<RealGradient> >& dphi = face_fe->get_dphi();
            const std::vector<std::vector<Real> >& phi = face_fe->get_phi();
            const std::vector<Point>& normal = face_fe->get_normals();

            const Elem* relevant_neighbor=elem->neighbor_ptr(prev_neighbor);

            // Get the correct face to the finite element; for this, continue looping
            // until a neighbor is an infinite element.
            // If none is found, we leave this loop
            for (; prev_neighbor<elem->n_neighbors(); prev_neighbor++){
               if (elem->neighbor_ptr(prev_neighbor)->infinite()){
                  relevant_neighbor=elem->neighbor_ptr(prev_neighbor);
                  break;
               }
            }
            if ( prev_neighbor == elem->n_neighbors()){
               break;
            }

            // now, we need to find the side of this element which faces the infinite element neighbor.
            // This is used as a consistency-check in debug-mode.
            unsigned int side=0;
            unsigned int found_s=0;
            // Trying to get a reasonable bound; I would like to avoid that it breaks.
            Real tol=0.01*elem->hmin();
            while (found_s==0){
               for(unsigned int n=0; n< elem->n_sides(); n++){
                  if (relevant_neighbor->close_to_point(elem->side_ptr(n)->vertex_average(), tol)){
                     found_s++;
                     side=n;
#ifndef DEBUG
                     // We should never find more than 1 element...
                     break;
#endif
                  }
               }
               tol*=2.;

#ifdef DEBUG
               // In case of 1/3 elem->hmin(), we can expect the point to be near all neighboring elements.
               if (tol > 0.2*elem->hmin()){
                  err<<"Warning: tolerance reached "<<tol<<","<<std::endl;
                  err<<"       but element-height is "<<elem->hmin()<<std::endl;
                  err<<"       (element id: "<<elem->id()<<")"<<std::endl;
               }
#endif
            }

#ifdef DEBUG
            libmesh_assert_equal_to(found_s, 1);
#endif
            //out<<"neighbor id: "<<relevant_neighbor->id();
            //out<<"  this id: "<<elem->id()<<std::endl;
            face_fe->reinit (elem, side);
            const unsigned int n_sf =
              FEInterface::n_dofs(face_fe->get_fe_type(), elem);

            H.resize (dof_indices.size(), dof_indices.size());

            for (unsigned int pts=0; pts<q_point.size(); pts++){
               // i=bra, j=ket
               for (unsigned int j=0; j<n_sf; j++){
                  for (unsigned int i=0; i<n_sf; i++){
                     libmesh_assert(normal[pts]*q_point[pts]/q_point[pts].norm() > 0);
                     H(i,j)+=.5*JxW[pts]*normal[pts]*phi[i][pts]*dphi[j][pts];
                  }
               }
            }
            matrix_A.add_matrix (H, dof_indices);

         } // end infinite neighbor-loop
      } // end of element loop
   }

  /**
   * All done!
   */

#else
  // Avoid unused variable warnings when compiling without infinite
  // elements and/or SLEPc enabled.
  libmesh_ignore(es, system_name);
#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS && LIBMESH_HAVE_SLEPC
}

#ifdef LIBMESH_HAVE_SLEPC
//define the functions needed for the @ SlepcSolverConfiguration.
void SlepcSolverConfiguration::configure_solver()
{
  // if a spectral transformation was requested
  if (_st!=INVALID_ST)
    {
      // initialise the st with the default values and change the spectral transformation value.
      ST st;
      PetscErrorCode ierr = EPSGetST(_slepc_solver.eps(), &st);
      if (ierr)
        libmesh_error();

      // Set it to the desired type of spectral transformation.
      // The value of the respective shift is chosen to be the target
      // specified via \p set_position_of_spectrum().
      switch (_st)
        {
        case SHIFT:
          ierr = STSetType(st, STSHIFT);
          break;
        case SINVERT:
          // this method has been renamed in 3.1
#if SLEPC_VERSION_LESS_THAN(3,1,0)
          ierr = STSetType(st, STSINV);
#else
          ierr = STSetType(st, STSINVERT);
#endif
          break;
        case CAYLEY:
          ierr = STSetType(st, STCAYLEY);
          break;
        default:
          // print a warning but do nothing more.
          break;
        }
      if (ierr)
        libmesh_error();

      // since st is a reference to the particular object used by \p _slepc_solver,
      // we don't need to hand back the manipulated object. It will be applied before
      // solving the system automatically.
    }
}
#endif // LIBMESH_HAVE_SLEPC

// To see the solution, even in the infinite-element region.
// The analytic solution is c*exp(-abs(|r|)) where r is the distance from the origin
// and C is some normalization constant (here ~0.08).
// The gradient is given just to show that it works.
void line_print(EquationSystems& es, std::string output, std::string SysName)
{

  // this function does not work without infinite elements properly.
#ifndef LIBMESH_ENABLE_INFINITE_ELEMENTS
  libmesh_ignore(es, output, SysName);
#else

  //Since we don't need any functionality that is special for EigenSystem,
  // we can use the base class here, thus keeping this function more general.
  System & system = es.get_system<System> (SysName);

  Real r = 2*es.parameters.get<Real>("radius");

  // get the variable number to specify, at which quantity we want to look later:
  unsigned short int phi_var=system.variable_number("phi");


  // set the full name of output-files:
  std::ostringstream re_output;
  re_output<<"re_"<<output;
  std::ostringstream im_output;
  im_output<<"im_"<<output;
  std::ostringstream abs_output;
  abs_output<<"abs_"<<output;

  // create output-objects, writing in respective files:
  std::ofstream im_out(re_output.str());
  std::ofstream re_out(im_output.str());
  std::ofstream abs_out(abs_output.str());

  Real N = 100.;
  Point q_point;
  const Real start=-r;
  Number soln;
  Gradient grad;

  for (int pts=1;pts<=N;pts++)
    {

      //specify the point to look at; going from -r to +r.
      q_point = Point(start+ 2*pts*r/N, 0., 0.);

      soln=system.point_value(phi_var, q_point);
      grad=system.point_gradient(phi_var, q_point);

      // and print them to the output-file:
      re_out<<" "<<std::setw(12)<<q_point(0);
      im_out<<" "<<std::setw(12)<<q_point(0);
      abs_out<<" "<<std::setw(12)<<q_point(0);
      re_out<<"  "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<std::real(soln);
      im_out<<"  "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<std::imag(soln);
      abs_out<<"  "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<std::abs(soln);
      // print the probability-current as well:
      re_out<<"  "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<std::imag(grad(0)*soln)<<std::endl;
      im_out<<"  "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<std::imag(grad(0)*soln)<<std::endl;
      abs_out<<"  "<<std::setw(12)<<std::scientific<<std::setprecision(6)<<std::abs(std::imag(grad(0)*soln))<<std::endl;

      // check that the solution is close to the analytic answer.
      // Comparing abs(soln) because there is a degree of freedom in the global phase.
      libmesh_assert_less(std::abs(std::abs(soln)-0.08*exp(-q_point.norm())), .002);

    }
#endif //LIBMESH_ENABLE_INFINITE_ELEMENTS
}
