#include <iostream>
#include <stdexcept>

#include "metaphysicl_config.h"

// If we have MASA we test ourselves against an MMS solution; if not
// we just test that this compiles.
#ifdef METAPHYSICL_HAVE_MASA
#  include <masa.h>
  using namespace MASA;
#else
#  define masa_get_param(arg) 1;
#endif // METAPHYSICL_HAVE_MASA

using namespace MetaPhysicL;

template <typename Vector>
static double evaluate_q (const Vector& xyz, const int);

using namespace MetaPhysicL;

int main(void)
{
  int N = 10; // mesh pts. in x and y
  double s2u,s2v,s2e,s2p;

#ifdef METAPHYSICL_HAVE_MASA
  int err = 0;
  double su,sv,sp,se;
  double pnorm, unorm, vnorm, enorm;
  double pnorm_max = 0., unorm_max = 0., vnorm_max = 0., enorm_max = 0.;
  double prnorm_max = 0., urnorm_max = 0., vrnorm_max = 0., ernorm_max = 0.;
#endif

  typedef VectorUnitVector<NDIM,0,RawScalar>::type XVector;
  XVector xvec = VectorUnitVector<NDIM,0,RawScalar>::value();

  typedef VectorUnitVector<NDIM,1,RawScalar>::type YVector;
  YVector yvec = VectorUnitVector<NDIM,1,RawScalar>::value();

  typedef DualNumber<RawScalar, XVector> XFirstDerivType;
  typedef DualNumber<RawScalar, YVector> YFirstDerivType;
//  typedef DualNumber<XFirstDerivType, XVector::rebind<XFirstDerivType>::other > XSecondDerivType;
//  typedef DualNumber<YFirstDerivType, YVector::rebind<YFirstDerivType>::other > YSecondDerivType;

  typedef XFirstDerivType XADType;
  typedef YFirstDerivType YADType;

  typedef VectorOf<NDIM, 0, XADType, 1, YADType>::type Vector;

#ifdef METAPHYSICL_HAVE_MASA
  // initialize the problem in MASA
  err += masa_init("euler-maple","euler_2d");
  
  // call the sanity check routine
  // (tests that all variables have been initialized)
  err += masa_sanity_check();
  //err += masa_printid<double>();
#endif // METAPHYSICL_HAVE_MASA

  Vector xy;
  xy.insert<0>() = XADType(1., xvec);
  xy.insert<1>() = YADType(1., yvec);

  // the input argument xyz is another NumberVector 
  // a vector just like Q_rho_u, a spatial location rather 
  // than a vector-valued forcing function.
  double h = 1.0/N;
  for (int i=0; i != N+1; ++i)
    {
      xy.get<0>() = XADType(i*h, xvec);

      for (int j=0; j != N+1; ++j)
	{
          xy.get<1>() = YADType(j*h, yvec);

	  // AD source terms
	  s2u = evaluate_q(xy,1);
	  s2v = evaluate_q(xy,2);
	  s2p = evaluate_q(xy,3);
	  s2e = evaluate_q(xy,4);

#ifdef METAPHYSICL_HAVE_MASA
	  // evaluate masa source terms
	  su  = masa_eval_source_rho_u<double>(i*h,j*h);
	  sv  = masa_eval_source_rho_v<double>(i*h,j*h);
	  sp  = masa_eval_source_rho  <double>(i*h,j*h);
	  se  = masa_eval_source_rho_e<double>(i*h,j*h);

	  unorm = fabs(su-s2u);	  
	  vnorm = fabs(sv-s2v);
	  pnorm = fabs(sp-s2p);	  
	  enorm = fabs(se-s2e);

	  double urnorm = fabs(su-s2u)/std::max(su,s2u);	  
	  double vrnorm = fabs(sv-s2v)/std::max(sv,s2v);
	  double prnorm = fabs(sp-s2p)/std::max(sp,s2p);	  
	  double ernorm = fabs(se-s2e)/std::max(se,s2e);

          unorm_max = std::max(unorm, unorm_max);
          vnorm_max = std::max(vnorm, vnorm_max);
          pnorm_max = std::max(pnorm, pnorm_max);
          enorm_max = std::max(enorm, enorm_max);

          urnorm_max = std::max(urnorm, urnorm_max);
          vrnorm_max = std::max(vrnorm, vrnorm_max);
          prnorm_max = std::max(prnorm, prnorm_max);
          ernorm_max = std::max(ernorm, ernorm_max);
#endif // METAPHYSICL_HAVE_MASA

	}
    }
 
#ifdef METAPHYSICL_HAVE_MASA
  std::cout << "max error in u      : " << unorm_max << std::endl;
  std::cout << "max error in v      : " << vnorm_max << std::endl;
  std::cout << "max error in density: " << pnorm_max << std::endl;
  std::cout << "max error in energy : " << enorm_max << std::endl;

  std::cout << "max relative error in u      : " << urnorm_max << std::endl;
  std::cout << "max relative error in v      : " << vrnorm_max << std::endl;
  std::cout << "max relative error in density: " << prnorm_max << std::endl;
  std::cout << "max relative error in energy : " << ernorm_max << std::endl;
#endif // METAPHYSICL_HAVE_MASA

  // steady as she goes...
  return 0;

}

// Note: ADScalar needs to be a FirstDerivType or better since we have
// first derivatives here.  Adding diffusion will require a
// SecondDerivType or better

template <typename Vector>
double evaluate_q (const Vector& xyz, const int ret)
{
  typedef typename Vector::value_type ADScalar;

  typedef typename RawType<ADScalar>::value_type Scalar;

  typedef typename Vector::template rebind<Scalar>::other RawVector;

  typedef typename Vector::template rebind<ADScalar>::other FullVector;

  const Scalar PI = std::acos(Scalar(-1));

  const Scalar u_0 = 200.23;
  const Scalar u_x = 1.1;
  const Scalar u_y = 1.08;
  const Scalar v_0 = 1.2;
  const Scalar v_x = 1.6;
  const Scalar v_y = .47;
  const Scalar rho_0 = 100.02;
  const Scalar rho_x = 2.22;
  const Scalar rho_y = 0.8;
  const Scalar p_0 = 150.2;
  const Scalar p_x = .91;
  const Scalar p_y = .623;
  const Scalar a_px = .165;
  const Scalar a_py = .612;
  const Scalar a_rhox = 1.0;
  const Scalar a_rhoy = 1.0;
  const Scalar a_ux = .1987;
  const Scalar a_uy = 1.189;
  const Scalar a_vx = 1.91;
  const Scalar a_vy = 1.0;
  const Scalar Gamma = 1.01;
  const Scalar L = 3.02;

  const typename Vector::template entry_type<0>::type& x =
    xyz.template get<0>();
  const typename Vector::template entry_type<1>::type& y =
    xyz.template get<1>();

  // Treat velocity as a vector
  FullVector U;

  // Arbitrary manufactured solution
  U.template insert<0>() = u_0 + u_x * std::sin(a_ux * PI * x / L) + u_y * std::cos(a_uy * PI * y / L);
  U.template insert<1>() = v_0 + v_x * std::cos(a_vx * PI * x / L) + v_y * std::sin(a_vy * PI * y / L);
  ADScalar RHO = rho_0 + rho_x * std::sin(a_rhox * PI * x / L) + rho_y * std::cos(a_rhoy * PI * y / L);
  ADScalar P = p_0 + p_x * std::cos(a_px * PI * x / L) + p_y * std::sin(a_py * PI * y / L);

  P/RHO;






  // Perfect gas energies
  ADScalar E = 1./(Gamma-1.)*P/RHO;
  ADScalar ET = E + .5 * U.dot(U);

  // Euler equation residuals
  Scalar Q_rho = raw_value(divergence(RHO*U));
  RawVector Q_rho_u = raw_value(divergence(RHO*U.outerproduct(U)) + P.derivatives());

  // energy equation
  Scalar Q_rho_e = raw_value(divergence((RHO*ET+P)*U));

  switch(ret)
    {

      // u
    case 1: 
      return Q_rho_u.template get<0>();
      break;

      // v
    case 2:
      return Q_rho_u.template get<1>();
      break;

      // rho
    case 3:
      return Q_rho;
      break;

      // energy
    case 4:
      return Q_rho_e;
      break;

    default:
      throw std::domain_error("Bad evaluate_q input request");
    }
  return 0.;
}
