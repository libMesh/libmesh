#include <iostream>
#include <typeinfo>

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

int main(void)
{
  int N   = 2; // mesh pts. in x and y
  double s2u,s2v,s2e,s2p;

#ifdef METAPHYSICL_HAVE_MASA
  int err = 0;
  double su,sv,sp,se;
  double pnorm, unorm, vnorm, enorm;
  double pnorm_max = 0., unorm_max = 0., vnorm_max = 0., enorm_max = 0.;
  double prnorm_max = 0., urnorm_max = 0., vrnorm_max = 0., ernorm_max = 0.;
#endif

  RawVector xvecinit(0), yvecinit(0);
  xvecinit.insert<0>() = 1.;
  yvecinit.insert<1>() = 1.;

  typedef VectorUnitVector<NDIM,0,RawScalar>::type XVector;
  XVector xvec = VectorUnitVector<NDIM,0,RawScalar>::value();

  typedef VectorUnitVector<NDIM,1,RawScalar>::type YVector;
  YVector yvec = VectorUnitVector<NDIM,1,RawScalar>::value();

  typedef DualNumber<RawScalar, XVector> XFirstDerivType;
  typedef DualNumber<RawScalar, YVector> YFirstDerivType;

  typedef DualNumber<XFirstDerivType, XVector::rebind<XFirstDerivType>::other > XSecondDerivType;
  typedef DualNumber<YFirstDerivType, YVector::rebind<YFirstDerivType>::other > YSecondDerivType;

  typedef XSecondDerivType XADType;
  typedef YSecondDerivType YADType;

  typedef VectorOf<NDIM, 0, XADType, 1, YADType>::type Vector;

#ifdef METAPHYSICL_HAVE_MASA
  // initialize the problem in MASA
  err += masa_init("ns-maple","navierstokes_2d_compressible");
  
  masa_set_param("u_0", 200.23);
  masa_set_param("u_y", 1.08);
  masa_set_param("v_0", 1.2);
  masa_set_param("v_y", .67);
  masa_set_param("rho_0", 100.02);
  masa_set_param("rho_x", 2.22);
  masa_set_param("rho_y", 0.8);
  masa_set_param("p_0", 150.2);
  masa_set_param("a_rhox", 1.0);
  masa_set_param("a_rhoy", 1.0);
  masa_set_param("a_vy", 1.0);
  masa_set_param("mu", 10.0);
  masa_set_param("k", 100.0);

  // call the sanity check routine
  // (tests that all variables have been initialized)
  err += masa_sanity_check();
  //err += masa_printid<double>();
#endif // METAPHYSICL_HAVE_MASA

  // we first set up the DualNumbers that correspond to independent
  // variables, spatial coordinates x and y.

  Vector xy;

  // When main() says "xy[0] = ADType(1., xvec);", that's saying "x = 1, and 
  // the gradient of f(x,y)=x is the constant vector xvec={1,0}"  
  // Likewise "xy[1] = ADType(1., yvec);" means "y = 1, and the gradient of f(x,y)=y 
  xy.insert<0>() = XADType(1., xvec);
  xy.insert<1>() = YADType(1., yvec);

  // For getting second derivatives, the way to set up a
  // twice-differentiable independent variable used to be more
  // complicated: first set up a once-differentiable variable, then
  // make sure *its* derivatives are also properly set.

  // However, if the new DualNumber constructors are working properly
  // then this is unnecessary.

  // xy[0] = ADType(FirstDerivType(1., xvec), xvec);
  // xy[1] = ADType(FirstDerivType(1., yvec), yvec);

  // the input argument xyz is another NumberVector 
  // a vector just like Q_rho_u, a spatial location rather 
  // than a vector-valued forcing function.
  double h = 1.0/N;
  for (int i=0; i != N+1; ++i)
    {
      //
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
#else
          // Avoid "set but not used" variable warnings;
          (void) s2u;
          (void) s2v;
          (void) s2p;
          (void) s2e;
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

  const Scalar R = masa_get_param("R");
  const Scalar u_0 = masa_get_param("u_0");
  const Scalar u_x = masa_get_param("u_x");
  const Scalar u_y = masa_get_param("u_y");
  const Scalar v_0 = masa_get_param("v_0");
  const Scalar v_x = masa_get_param("v_x");
  const Scalar v_y = masa_get_param("v_y");
  const Scalar rho_0 = masa_get_param("rho_0");
  const Scalar rho_x = masa_get_param("rho_x");
  const Scalar rho_y = masa_get_param("rho_y");
  const Scalar p_0 = masa_get_param("p_0");
  const Scalar p_x = masa_get_param("p_x");
  const Scalar p_y = masa_get_param("p_y");
  const Scalar a_px = masa_get_param("a_px");
  const Scalar a_py = masa_get_param("a_py");
  const Scalar a_rhox = masa_get_param("a_rhox");
  const Scalar a_rhoy = masa_get_param("a_rhoy");
  const Scalar a_ux = masa_get_param("a_ux");
  const Scalar a_uy = masa_get_param("a_uy");
  const Scalar a_vx = masa_get_param("a_vx");
  const Scalar a_vy = masa_get_param("a_vy");
  const Scalar Gamma = masa_get_param("Gamma");
  const Scalar L = masa_get_param("L");
  const Scalar mu = masa_get_param("mu");
  const Scalar k = masa_get_param("k");

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

  // Temperature
  ADScalar T = P / RHO / R;

  // Perfect gas energies
  ADScalar E = 1./(Gamma-1.)*P/RHO;
  ADScalar ET = E + .5 * U.dot(U);

  // The shear strain tensor
  typedef typename Vector::template rebind<typename ADScalar::derivatives_type>::other Tensor;
  Tensor GradU = gradient(U);

  // The identity tensor I
  // typedef typename Vector::template rebind<RawVector>::other RawTensor;

  // RawTensor Identity = RawVector::identity();

  // The shear stress tensor
  // Tensor Tau = mu * (GradU + transpose(GradU) - 2./3.*divergence(U)*Identity);
  Tensor Tau = mu * (GradU + transpose(GradU) - 2./3.*divergence(U)*RawVector::identity(2));

  // Temperature flux
  FullVector q = -k * T.derivatives();

  // Euler equation residuals
  Scalar Q_rho = raw_value(divergence(RHO*U));
  RawVector Q_rho_u = raw_value(divergence(RHO*U.outerproduct(U) - Tau) + P.derivatives());

  // energy equation
  Scalar Q_rho_e = raw_value(divergence((RHO*ET+P)*U + q - Tau.dot(U)));

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
    }

  std::cout << "something is wrong!\n";
  exit(1);

  return 0;
}
