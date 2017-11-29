#include <cstdlib> // rand()
#include <iostream>
#include <limits>

#include "metaphysicl_config.h"

#include "metaphysicl/dualdynamicsparsenumbervector.h"
#include "metaphysicl/dualnumbervector.h"
#include "metaphysicl/dualsparsenumbervector.h"

using namespace MetaPhysicL;

const unsigned int nx = 10, ny = 10, nz = 10;

template <typename Scalar>
bool test_error (Scalar computed,
                Scalar analytic)
{
  using std::max;
  using std::fabs;

  static const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;

  Scalar max_abs_error = fabs(computed - analytic);

  // Handle NaNs properly.  Testing max_abs_error for NaN is
  // impossible because IEEE sucks:
  // https://en.wikipedia.org/wiki/IEEE_754_revision#min_and_max
  if (max_abs_error > tol || computed != computed)
    {
      std::cerr << "Value = " << analytic <<
	           "\nComputed = " << computed <<
	           "\nTol   " << tol << std::endl;
      return true; // failure
    }

  return false; // success
}

template <typename Vector, typename VectorX, typename VectorY, typename VectorZ>
bool vectester (const Vector & zero_vec,
                const VectorX & unit_x,
                const VectorY & unit_y,
                const VectorZ & unit_z)
{
  typedef typename ValueType<Vector>::type ADScalar;
  typedef typename RawType<ADScalar>::value_type Scalar;

  bool return_failure = false;

  for (unsigned int i=0; i != nx; ++i)
    {
      const Scalar x1 = Scalar(i)/nx;
      const ADScalar x = ADScalar(x1, unit_x);
      for (unsigned int j=0; j != ny; ++j)
        {
          const Scalar y1 = Scalar(j)/ny;
          const ADScalar y = ADScalar(y1, unit_y);
          for (unsigned int k=0; k != nz; ++k)
            {
              const Scalar z1 = Scalar(k)/nz;
              const ADScalar z = ADScalar(z1, unit_z);

              Vector U;
              U[0] = x*x*x + y*y + z*z;
              U[1] = x*x + y*y*y + z*z;
              U[2] = x*x + y*y + z*z*z;
               
              Scalar divgrad0 = raw_value(divergence(gradient(U)))[0];
              Scalar divgrad1 = raw_value(divergence(gradient(U)))[1];
              Scalar divgrad2 = raw_value(divergence(gradient(U)))[2];
              Scalar true_divgrad0 = 6*raw_value(x)+4;
              Scalar true_divgrad1 = 6*raw_value(y)+4;
              Scalar true_divgrad2 = 6*raw_value(z)+4;
              return_failure = return_failure ||
                test_error(divgrad0, true_divgrad0);
              return_failure = return_failure ||
                test_error(divgrad1, true_divgrad1);
              return_failure = return_failure ||
                test_error(divgrad2, true_divgrad2);
            }
        }
    }

  return return_failure;
}

template <typename Scalar>
bool dense_tester()
{
  typedef DualNumber<Scalar, NumberVector<3, Scalar> >
    FirstDerivType;

  typedef DualNumber<FirstDerivType, NumberVector<3, FirstDerivType> >
    SecondDerivType;

  typedef NumberVector<3, Scalar>
    PlainVectorType;

  typedef NumberVector<3, SecondDerivType>
    SecondVectorType;

  const PlainVectorType unit_x = 
    NumberVectorUnitVector <3, 0, Scalar>::value();
  const PlainVectorType unit_y = 
    NumberVectorUnitVector <3, 1, Scalar>::value();
  const PlainVectorType unit_z = 
    NumberVectorUnitVector <3, 2, Scalar>::value();

  const SecondVectorType zero_vec = 0;

  return vectester(zero_vec, unit_x, unit_y, unit_z);
}

int main(void)
{
  int returnval = 0;

  returnval = returnval || dense_tester<float>();
  returnval = returnval || dense_tester<double>();
  returnval = returnval || dense_tester<long double>();


#if 0
  returnval = returnval || vectester(SparseNumberVectorOf
    <N, 0, DualNumber<double>, 1, DualNumber<double>,
        2, DualNumber<double>, 3, DualNumber<double> >::type());
  returnval = returnval || vectester(SparseNumberVectorOf
    <N, 0, DualNumber<long double>, 1, DualNumber<long double>,
        2, DualNumber<long double>, 3, DualNumber<long double> >::type());

  DynamicSparseNumberVector<DualNumber<float>, unsigned int> float_dsnv;
    float_dsnv.resize(4);
    float_dsnv.raw_index(1) = 1;
    float_dsnv.raw_index(2) = 2;
    float_dsnv.raw_index(3) = 3;
  returnval = returnval || vectester(float_dsnv);

  DynamicSparseNumberVector<DualNumber<double>, unsigned int> double_dsnv;
    double_dsnv.resize(4);
    double_dsnv.raw_index(1) = 1;
    double_dsnv.raw_index(2) = 2;
    double_dsnv.raw_index(3) = 3;
  returnval = returnval || vectester(double_dsnv);

  DynamicSparseNumberVector<DualNumber<long double>, unsigned int> long_double_dsnv;
    long_double_dsnv.resize(4);
    long_double_dsnv.raw_index(1) = 1;
    long_double_dsnv.raw_index(2) = 2;
    long_double_dsnv.raw_index(3) = 3;
  returnval = returnval || vectester(long_double_dsnv);
#endif

  return returnval;
}
