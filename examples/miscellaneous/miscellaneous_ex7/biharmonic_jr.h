#ifndef BIHARMONIC_JR_H
#define BIHARMONIC_JR_H

// LibMesh includes
#include "libmesh/transient_system.h"
#include "libmesh/nonlinear_solver.h"


// Example includes
#include "biharmonic.h"

// Bring in bits from the libMesh namespace.
// Just the bits we're using, since this is a header.
using libMesh::EquationSystems;
using libMesh::Gradient;
using libMesh::NonlinearImplicitSystem;
using libMesh::Number;
using libMesh::NumericVector;
using libMesh::Parameters;
using libMesh::Point;
using libMesh::SparseMatrix;
using libMesh::System;
using libMesh::TransientNonlinearImplicitSystem;


/**
 * Biharmonic's friend class definition
 */
class Biharmonic::JR : public TransientNonlinearImplicitSystem,
                       public NonlinearImplicitSystem::ComputeResidualandJacobian,
                       public NonlinearImplicitSystem::ComputeBounds,
                       public System::Initialization
{
public:
  /**
   * Constructor.
   */
  JR(EquationSystems & eqSys,
     const std::string & name,
     const unsigned int number);

  void initialize();

  /**
   * Static functions to be used for initialization
   */
  static Number InitialDensityBall(const Point & p,
                                   const Parameters & parameters,
                                   const std::string &,
                                   const std::string &);

  static Number InitialDensityRod(const Point & p,
                                  const Parameters & parameters,
                                  const std::string &,
                                  const std::string &);

  static Number InitialDensityStrip(const Point & p,
                                    const Parameters & parameters,
                                    const std::string &,
                                    const std::string &);

  static Gradient InitialGradientZero(const Point &,
                                      const Parameters &,
                                      const std::string &,
                                      const std::string &);

  /**
   * The residual and Jacobian assembly function for the Biharmonic system.
   */
  void residual_and_jacobian(const NumericVector<Number> & u,
                             NumericVector<Number> * R,
                             SparseMatrix<Number> * J,
                             NonlinearImplicitSystem &);

  /**
   * Function defining the bounds of the Biharmonic system.
   */
  void bounds(NumericVector<Number> & XL,
              NumericVector<Number> & XU,
              NonlinearImplicitSystem &);

private:
  Biharmonic & _biharmonic;
};


#endif // BIHARMONIC_JR_H
