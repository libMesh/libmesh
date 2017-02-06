// Libmesh includes
#include "libmesh/mesh.h"
#include "libmesh/quadrature.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/fourth_error_estimators.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/periodic_boundaries.h"
#include "libmesh/periodic_boundary.h"
#include "libmesh/elem.h"

// Example includes
#include "biharmonic_jr.h"

using namespace libMesh;

Biharmonic::JR::JR(EquationSystems & eqSys,
                   const std::string & name_in,
                   const unsigned int number_in) :
  TransientNonlinearImplicitSystem(eqSys, name_in, number_in),
  _biharmonic(dynamic_cast<Biharmonic &>(eqSys))
{
  // Check that we can actually compute second derivatives
#ifndef LIBMESH_ENABLE_SECOND_DERIVATIVES
  libmesh_error_msg("Must have second derivatives enabled");
#endif

#ifdef LIBMESH_ENABLE_PERIODIC
  // Add periodicity to the mesh
  DofMap & dof_map = get_dof_map();
  PeriodicBoundary xbdry(RealVectorValue(1.0, 0.0, 0.0));
#if LIBMESH_DIM > 1
  PeriodicBoundary ybdry(RealVectorValue(0.0, 1.0, 0.0));
#endif
#if LIBMESH_DIM > 2
  PeriodicBoundary zbdry(RealVectorValue(0.0, 0.0, 1.0));
#endif

  switch(_biharmonic._dim)
    {
    case 1:
      xbdry.myboundary = 0;
      xbdry.pairedboundary = 1;
      dof_map.add_periodic_boundary(xbdry);
      break;
#if LIBMESH_DIM > 1
    case 2:
      xbdry.myboundary = 3;
      xbdry.pairedboundary = 1;
      dof_map.add_periodic_boundary(xbdry);
      ybdry.myboundary = 0;
      ybdry.pairedboundary = 2;
      dof_map.add_periodic_boundary(ybdry);
      break;
#endif
#if LIBMESH_DIM > 2
    case 3:
      xbdry.myboundary = 4;
      xbdry.pairedboundary = 2;
      dof_map.add_periodic_boundary(xbdry);
      ybdry.myboundary = 1;
      ybdry.pairedboundary = 3;
      dof_map.add_periodic_boundary(ybdry);
      zbdry.myboundary = 0;
      zbdry.pairedboundary = 5;
      dof_map.add_periodic_boundary(zbdry);
      break;
#endif
    default:
      libmesh_error_msg("Invalid dimension = " << _biharmonic._dim);
    }
#endif // LIBMESH_ENABLE_PERIODIC

  // Adds the variable "u" to the system.
  // u will be approximated using Hermite elements
  add_variable("u", THIRD, HERMITE);

  // Give the system an object to compute the initial state.
  attach_init_object(*this);

  // Attache the R & J calculation object
  nonlinear_solver->residual_and_jacobian_object = this;

  // Attach the bounds calculation object
  nonlinear_solver->bounds_object = this;
}





void Biharmonic::JR::initialize()
{
  if (_biharmonic._verbose)
    libMesh::out << ">>> Initializing Biharmonic::JR\n";

  Parameters parameters;
  parameters.set<Point>("center") = _biharmonic._initialCenter;
  parameters.set<Real>("width")   = _biharmonic._initialWidth;

  if (_biharmonic._initialState == Biharmonic::BALL)
    project_solution(Biharmonic::JR::InitialDensityBall, Biharmonic::JR::InitialGradientZero, parameters);

  if (_biharmonic._initialState == Biharmonic::ROD)
    project_solution(Biharmonic::JR::InitialDensityRod, Biharmonic::JR::InitialGradientZero, parameters);

  if (_biharmonic._initialState == Biharmonic::STRIP)
    project_solution(Biharmonic::JR::InitialDensityStrip, Biharmonic::JR::InitialGradientZero, parameters);

  // both states are equal
  *(old_local_solution) = *(current_local_solution);

  if (_biharmonic._verbose)
    libMesh::out << "<<< Initializing Biharmonic::JR\n";
}






Number Biharmonic::JR::InitialDensityBall(const Point & p,
                                          const Parameters & parameters,
                                          const std::string &,
                                          const std::string &)
{
  // Initialize with a ball in the middle, which is a segment in 1D, a disk in 2D and a ball in 3D.
  Point center = parameters.get<Point>("center");
  Real width = parameters.get<Real>("width");
  Point pc = p-center;
  Real r = pc.norm();
  return (r < width) ? 1.0 : -0.5;
}




Number Biharmonic::JR::InitialDensityRod(const Point & p,
                                         const Parameters & parameters,
                                         const std::string &,
                                         const std::string &)
{
  // Initialize with a rod in the middle so that we have a z-homogeneous system to model the 2D disk.
  Point center = parameters.get<Point>("center");
  Real width = parameters.get<Real>("width");
  Real r = sqrt((p(0)-center(0))*(p(0)-center(0)) + (p(1)-center(1))*(p(1)-center(1)));
  return (r < width) ? 1.0 : -0.5;
}





Number Biharmonic::JR::InitialDensityStrip(const Point & p,
                                           const Parameters & parameters,
                                           const std::string &,
                                           const std::string &)
{
  // Initialize with a wide strip in the middle so that we have a yz-homogeneous system to model the 1D.
  Point center = parameters.get<Point>("center");
  Real width = parameters.get<Real>("width");
  Real r = sqrt((p(0)-center(0))*(p(0)-center(0)));
  return (r < width) ? 1.0 : -0.5;
}




Gradient Biharmonic::JR::InitialGradientZero(const Point &,
                                             const Parameters &,
                                             const std::string &,
                                             const std::string &)
{
  return Gradient(0.0, 0.0, 0.0);
}




void Biharmonic::JR::residual_and_jacobian(const NumericVector<Number> & u,
                                           NumericVector<Number> * R,
                                           SparseMatrix<Number> * J,
                                           NonlinearImplicitSystem &)
{
#ifdef LIBMESH_ENABLE_SECOND_DERIVATIVES
  if (!R && !J)
    return;

  // Declare a performance log.  Give it a descriptive
  // string to identify what part of the code we are
  // logging, since there may be many PerfLogs in an
  // application.
  PerfLog perf_log ("Biharmonic Residual and Jacobian", false);

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the DofMap
  // in future examples.
  const DofMap & dof_map = get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> fe (FEBase::build(_biharmonic._dim, fe_type));

  // Quadrature rule for numerical integration.
  // With 2D triangles, the Clough quadrature rule puts a Gaussian
  // quadrature rule on each of the 3 subelements
  UniquePtr<QBase> qrule(fe_type.default_quadrature_rule(_biharmonic._dim));

  // Tell the finite element object to use our quadrature rule.
  fe->attach_quadrature_rule (qrule.get());

  // Here we define some references to element-specific data that
  // will be used to assemble the linear system.
  // We begin with the element Jacobian * quadrature weight at each
  // integration point.
  const std::vector<Real> & JxW = fe->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<Real> > & phi = fe->get_phi();

  // The element shape functions' derivatives evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

  // The element shape functions'  second derivatives evaluated at the quadrature points.
  const std::vector<std::vector<RealTensor> > & d2phi = fe->get_d2phi();

  // For efficiency we will compute shape function laplacians n times,
  // not n^2
  std::vector<Real> Laplacian_phi_qp;

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Je" and "Re". More detail is in example 3.
  DenseMatrix<Number> Je;
  DenseVector<Number> Re;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  // Old solution
  const NumericVector<Number> & u_old = *old_local_solution;

  // Now we will loop over all the elements in the mesh.  We will
  // compute the element matrix and right-hand-side contribution.  See
  // example 3 for a discussion of the element iterators.

  MeshBase::const_element_iterator       el     = _biharmonic._mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _biharmonic._mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      // Store a pointer to the element we are currently
      // working on.  This allows for nicer syntax later.
      const Elem * elem = *el;

      // Get the degree of freedom indices for the
      // current element.  These define where in the global
      // matrix and right-hand-side this element will
      // contribute to.
      dof_map.dof_indices (elem, dof_indices);

      // Compute the element-specific data for the current
      // element.  This involves computing the location of the
      // quadrature points (q_point) and the shape function
      // values/derivatives (phi, dphi, d2phi) for the current element.
      fe->reinit (elem);

      // Zero the element matrix, the right-hand side and the Laplacian matrix
      // before summing them.
      if (J)
        Je.resize(dof_indices.size(), dof_indices.size());

      if (R)
        Re.resize(dof_indices.size());

      Laplacian_phi_qp.resize(dof_indices.size());

      for (unsigned int qp=0; qp<qrule->n_points(); qp++)
        {
          // AUXILIARY QUANTITIES:
          // Residual and Jacobian share a few calculations:
          // at the very least, in the case of interfacial energy only with a constant mobility,
          // both calculations use Laplacian_phi_qp; more is shared the case of a concentration-dependent
          // mobility and bulk potentials.
          Number
            u_qp = 0.0,
            u_old_qp = 0.0,
            Laplacian_u_qp = 0.0,
            Laplacian_u_old_qp = 0.0;

          Gradient
            grad_u_qp(0.0, 0.0, 0.0),
            grad_u_old_qp(0.0, 0.0, 0.0);

          Number
            M_qp = 1.0,
            M_old_qp = 1.0,
            M_prime_qp = 0.0,
            M_prime_old_qp = 0.0;

          for (std::size_t i=0; i<phi.size(); i++)
            {
              Laplacian_phi_qp[i] = d2phi[i][qp](0, 0);
              grad_u_qp(0) += u(dof_indices[i])*dphi[i][qp](0);
              grad_u_old_qp(0) += u_old(dof_indices[i])*dphi[i][qp](0);

              if (_biharmonic._dim > 1)
                {
                  Laplacian_phi_qp[i] += d2phi[i][qp](1, 1);
                  grad_u_qp(1) += u(dof_indices[i])*dphi[i][qp](1);
                  grad_u_old_qp(1) += u_old(dof_indices[i])*dphi[i][qp](1);
                }
              if (_biharmonic._dim > 2)
                {
                  Laplacian_phi_qp[i] += d2phi[i][qp](2, 2);
                  grad_u_qp(2) += u(dof_indices[i])*dphi[i][qp](2);
                  grad_u_old_qp(2) += u_old(dof_indices[i])*dphi[i][qp](2);
                }
              u_qp     += phi[i][qp]*u(dof_indices[i]);
              u_old_qp += phi[i][qp]*u_old(dof_indices[i]);
              Laplacian_u_qp     += Laplacian_phi_qp[i]*u(dof_indices[i]);
              Laplacian_u_old_qp += Laplacian_phi_qp[i]*u_old(dof_indices[i]);
            } // for i

          if (_biharmonic._degenerate)
            {
              M_qp           = 1.0 - u_qp*u_qp;
              M_old_qp       = 1.0 - u_old_qp*u_old_qp;
              M_prime_qp     = -2.0*u_qp;
              M_prime_old_qp = -2.0*u_old_qp;
            }

          // ELEMENT RESIDUAL AND JACOBIAN
          for (std::size_t i=0; i<phi.size(); i++)
            {
              // RESIDUAL
              if (R)
                {
                  Number ri = 0.0, ri_old = 0.0;
                  ri     -= Laplacian_phi_qp[i]*M_qp*_biharmonic._kappa*Laplacian_u_qp;
                  ri_old -= Laplacian_phi_qp[i]*M_old_qp*_biharmonic._kappa*Laplacian_u_old_qp;

                  if (_biharmonic._degenerate)
                    {
                      ri       -= (dphi[i][qp]*grad_u_qp)*M_prime_qp*(_biharmonic._kappa*Laplacian_u_qp);
                      ri_old   -= (dphi[i][qp]*grad_u_old_qp)*M_prime_old_qp*(_biharmonic._kappa*Laplacian_u_old_qp);
                    }

                  if (_biharmonic._cahn_hillard)
                    {
                      if (_biharmonic._energy == DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_WELL)
                        {
                          ri += Laplacian_phi_qp[i]*M_qp*_biharmonic._theta_c*(u_qp*u_qp - 1.0)*u_qp;
                          ri_old += Laplacian_phi_qp[i]*M_old_qp*_biharmonic._theta_c*(u_old_qp*u_old_qp - 1.0)*u_old_qp;
                          if (_biharmonic._degenerate)
                            {
                              ri     += (dphi[i][qp]*grad_u_qp)*M_prime_qp*_biharmonic._theta_c*(u_qp*u_qp - 1.0)*u_qp;
                              ri_old += (dphi[i][qp]*grad_u_old_qp)*M_prime_old_qp*_biharmonic._theta_c*(u_old_qp*u_old_qp - 1.0)*u_old_qp;
                            }
                        }// if (_biharmonic._energy == DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_WELL)

                      if (_biharmonic._energy == DOUBLE_OBSTACLE || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
                        {
                          ri -= Laplacian_phi_qp[i]*M_qp*_biharmonic._theta_c*u_qp;
                          ri_old -= Laplacian_phi_qp[i]*M_old_qp*_biharmonic._theta_c*u_old_qp;
                          if (_biharmonic._degenerate)
                            {
                              ri     -= (dphi[i][qp]*grad_u_qp)*M_prime_qp*_biharmonic._theta_c*u_qp;
                              ri_old -= (dphi[i][qp]*grad_u_old_qp)*M_prime_old_qp*_biharmonic._theta_c*u_old_qp;
                            }
                        } // if (_biharmonic._energy == DOUBLE_OBSTACLE || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)

                      if (_biharmonic._energy == LOG_DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
                        {
                          switch(_biharmonic._log_truncation)
                            {
                            case 2:
                              break;
                            case 3:
                              break;
                            default:
                              break;
                            }// switch(_biharmonic._log_truncation)
                        }// if (_biharmonic._energy == LOG_DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
                    }// if (_biharmonic._cahn_hillard)
                  Re(i) += JxW[qp]*((u_qp-u_old_qp)*phi[i][qp]-_biharmonic._dt*0.5*((2.0-_biharmonic._cnWeight)*ri + _biharmonic._cnWeight*ri_old));
                } // if (R)

              // JACOBIAN
              if (J)
                {
                  Number M_prime_prime_qp = 0.0;
                  if (_biharmonic._degenerate) M_prime_prime_qp = -2.0;
                  for (std::size_t j=0; j<phi.size(); j++)
                    {
                      Number ri_j = 0.0;
                      ri_j -= Laplacian_phi_qp[i]*M_qp*_biharmonic._kappa*Laplacian_phi_qp[j];
                      if (_biharmonic._degenerate)
                        {
                          ri_j -=
                            Laplacian_phi_qp[i]*M_prime_qp*phi[j][qp]*_biharmonic._kappa*Laplacian_u_qp               +
                            (dphi[i][qp]*dphi[j][qp])*M_prime_qp*(_biharmonic._kappa*Laplacian_u_qp)                  +
                            (dphi[i][qp]*grad_u_qp)*(M_prime_prime_qp*phi[j][qp])*(_biharmonic._kappa*Laplacian_u_qp) +
                            (dphi[i][qp]*grad_u_qp)*(M_prime_qp)*(_biharmonic._kappa*Laplacian_phi_qp[j]);
                        }

                      if (_biharmonic._cahn_hillard)
                        {
                          if (_biharmonic._energy == DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_WELL)
                            {
                              ri_j +=
                                Laplacian_phi_qp[i]*M_prime_qp*phi[j][qp]*_biharmonic._theta_c*(u_qp*u_qp - 1.0)*u_qp +
                                Laplacian_phi_qp[i]*M_qp*_biharmonic._theta_c*(3.0*u_qp*u_qp - 1.0)*phi[j][qp]        +
                                (dphi[i][qp]*dphi[j][qp])*M_prime_qp*_biharmonic._theta_c*(u_qp*u_qp - 1.0)*u_qp      +
                                (dphi[i][qp]*grad_u_qp)*M_prime_prime_qp*_biharmonic._theta_c*(u_qp*u_qp - 1.0)*u_qp  +
                                (dphi[i][qp]*grad_u_qp)*M_prime_qp*_biharmonic._theta_c*(3.0*u_qp*u_qp - 1.0)*phi[j][qp];
                            }// if (_biharmonic._energy == DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_WELL)

                          if (_biharmonic._energy == DOUBLE_OBSTACLE || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
                            {
                              ri_j -=
                                Laplacian_phi_qp[i]*M_prime_qp*phi[j][qp]*_biharmonic._theta_c*u_qp                   +
                                Laplacian_phi_qp[i]*M_qp*_biharmonic._theta_c*phi[j][qp]                              +
                                (dphi[i][qp]*dphi[j][qp])*M_prime_qp*_biharmonic._theta_c*u_qp                        +
                                (dphi[i][qp]*grad_u_qp)*M_prime_prime_qp*_biharmonic._theta_c*u_qp                    +
                                (dphi[i][qp]*grad_u_qp)*M_prime_qp*_biharmonic._theta_c*phi[j][qp];
                            } // if (_biharmonic._energy == DOUBLE_OBSTACLE || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)

                          if (_biharmonic._energy == LOG_DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
                            {
                              switch(_biharmonic._log_truncation)
                                {
                                case 2:
                                  break;
                                case 3:
                                  break;
                                default:
                                  break;
                                }// switch(_biharmonic._log_truncation)
                            }// if (_biharmonic._energy == LOG_DOUBLE_WELL || _biharmonic._energy == LOG_DOUBLE_OBSTACLE)
                        }// if (_biharmonic._cahn_hillard)
                      Je(i,j) += JxW[qp]*(phi[i][qp]*phi[j][qp] - 0.5*_biharmonic._dt*(2.0-_biharmonic._cnWeight)*ri_j);
                    } // for j
                } // if (J)
            } // for i
        } // for qp

      // The element matrix and right-hand-side are now built
      // for this element.  Add them to the global matrix and
      // right-hand-side vector.  The SparseMatrix::add_matrix()
      // and NumericVector::add_vector() members do this for us.
      // Start logging the insertion of the local (element)
      // matrix and vector into the global matrix and vector
      if (R)
        {
          // If the mesh has hanging nodes (e.g., as a result of refinement), those need to be constrained.
          dof_map.constrain_element_vector(Re, dof_indices);
          R->add_vector(Re, dof_indices);
        }

      if (J)
        {
          // If the mesh has hanging nodes (e.g., as a result of refinement), those need to be constrained.
          dof_map.constrain_element_matrix(Je, dof_indices);
          J->add_matrix(Je, dof_indices);
        }
    } // for el
#endif // LIBMESH_ENABLE_SECOND_DERIVATIVES
}





void Biharmonic::JR::bounds(NumericVector<Number> & XL,
                            NumericVector<Number> & XU,
                            NonlinearImplicitSystem &)
{
  // sys is actually ignored, since it should be the same as *this.

  // Declare a performance log.  Give it a descriptive
  // string to identify what part of the code we are
  // logging, since there may be many PerfLogs in an
  // application.
  PerfLog perf_log ("Biharmonic bounds", false);

  // A reference to the DofMap object for this system.  The DofMap
  // object handles the index translation from node and element numbers
  // to degree of freedom numbers.  We will talk more about the DofMap
  // in future examples.
  const DofMap & dof_map = get_dof_map();

  // Get a constant reference to the Finite Element type
  // for the first (and only) variable in the system.
  FEType fe_type = dof_map.variable_type(0);

  // Build a Finite Element object of the specified type.  Since the
  // FEBase::build() member dynamically creates memory we will
  // store the object as a UniquePtr<FEBase>.  This can be thought
  // of as a pointer that will clean up after itself.
  UniquePtr<FEBase> fe (FEBase::build(_biharmonic._dim, fe_type));

  // Define data structures to contain the bound vectors contributions.
  DenseVector<Number> XLe, XUe;

  // These vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;

  MeshBase::const_element_iterator       el     = _biharmonic._mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = _biharmonic._mesh.active_local_elements_end();

  for ( ; el != end_el; ++el)
    {
      // Extract the shape function to be evaluated at the nodes
      const std::vector<std::vector<Real> > & phi = fe->get_phi();

      // Get the degree of freedom indices for the current element.
      // They are in 1-1 correspondence with shape functions phi
      // and define where in the global vector this element will.
      dof_map.dof_indices (*el, dof_indices);

      // Resize the local bounds vectors (zeroing them out in the process).
      XLe.resize(dof_indices.size());
      XUe.resize(dof_indices.size());

      // Extract the element node coordinates in the reference frame
      std::vector<Point> nodes;
      fe->get_refspace_nodes((*el)->type(), nodes);

      // Evaluate the shape functions at the nodes
      fe->reinit(*el, &nodes);

      // Construct the bounds based on the value of the i-th phi at the nodes.
      // Observe that this doesn't really work in general: we rely on the fact
      // that for Hermite elements each shape function is nonzero at most at a
      // single node.
      // More generally the bounds must be constructed by inspecting a "mass-like"
      // matrix (m_{ij}) of the shape functions (i) evaluated at their corresponding nodes (j).
      // The constraints imposed on the dofs (d_i) are then are -1 \leq \sum_i d_i m_{ij} \leq 1,
      // since \sum_i d_i m_{ij} is the value of the solution at the j-th node.
      // Auxiliary variables will need to be introduced to reduce this to a "box" constraint.
      // Additional complications will arise since m might be singular (as is the case for Hermite,
      // which, however, is easily handled by inspection).
      for (std::size_t i=0; i<phi.size(); ++i)
        {
          // FIXME: should be able to define INF and pass it to the solve
          Real infinity = 1.0e20;
          Real bound = infinity;
          for (std::size_t j = 0; j < nodes.size(); ++j)
            {
              if (phi[i][j])
                {
                  bound = 1.0/std::abs(phi[i][j]);
                  break;
                }
            }

          // The value of the solution at this node must be between 1.0 and -1.0.
          // Based on the value of phi(i)(i) the nodal coordinate must be between 1.0/phi(i)(i) and its negative.
          XLe(i) = -bound;
          XUe(i) = bound;
        }
      // The element bound vectors are now built for this element.
      // Insert them into the global vectors, potentially overwriting
      // the same dof contributions from other elements: no matter --
      // the bounds are always -1.0 and 1.0.
      XL.insert(XLe, dof_indices);
      XU.insert(XUe, dof_indices);
    }
}
