#include <cstdlib> // rand()
#include <iostream>
#include <limits>

#include "metaphysicl_config.h"

#include "metaphysicl/dynamicsparsenumberarray.h"
#include "metaphysicl/dynamicsparsenumbervector.h"
#include "metaphysicl/numberarray.h"
#include "metaphysicl/numbervector.h"
#include "metaphysicl/sparsenumberarray.h"
#include "metaphysicl/sparsenumbervector.h"

static const unsigned int N = 10; // test pts.

using namespace MetaPhysicL;

#define one_test(test_func) \
  error_vec = test_func; \
  { int new_returnval = test_error_vec(random_vec, error_vec); \
  if (new_returnval) \
    std::cerr << "Failed test: " << #test_func << std::endl; \
  returnval = returnval || new_returnval; }

template <typename Vector>
int test_error_vec (const Vector& random_vec,
                    const Vector& error_vec)
{
  using std::max;
  using std::fabs;

  typedef typename Vector::value_type Scalar;

  static const Scalar tol = std::numeric_limits<Scalar>::epsilon() * 10;

  Scalar max_abs_error = 0;

  for (unsigned int i=0; i != error_vec.size(); ++i)
    {
      max_abs_error = max(max_abs_error, fabs(error_vec[i]));

      // Handle NaNs properly.  Testing max_abs_error for NaN is
      // impossible because IEEE sucks:
      // https://en.wikipedia.org/wiki/IEEE_754_revision#min_and_max
      if (max_abs_error > tol || error_vec[i] != error_vec[i])
        {
	  std::cerr << "Value " << random_vec[i] <<
		       "\nError " << error_vec[i] <<
		       "\nTol   " << tol << std::endl;
	  return 1;
        }
    }

  return 0;
}

template <typename Vector>
int vectester (Vector zerovec)
{
  using std::abs;
  using std::acos;
  using std::asin;
  using std::atan;
  using std::floor;
  using std::pow;
  using std::sin;
  using std::sqrt;
  using std::tan;

  typedef typename Vector::value_type Scalar;

  Vector random_vec = zerovec;

  Vector error_vec = zerovec;

  std::srand(12345); // Fixed seed for reproduceability of failures

  for (unsigned int i=0; i != random_vec.size(); ++i)
    random_vec[i] = .25 + (static_cast<Scalar>(std::rand())/RAND_MAX);

  int returnval = 0;

  one_test(2*random_vec - random_vec - random_vec);

  one_test(3*random_vec - random_vec*3);

  one_test((random_vec + random_vec)/2 - random_vec);

  one_test(sqrt(random_vec) * sqrt(random_vec) - random_vec);
  one_test(random_vec*random_vec - pow(random_vec,2));
  one_test(sqrt(random_vec) - pow(random_vec,Scalar(.5)));

  one_test(random_vec - sin(asin(random_vec)));
  one_test(random_vec - tan(atan(random_vec)));

  one_test(floor(random_vec / 2));
  one_test(abs(random_vec) - random_vec);

  return returnval;
}



template <typename Vector>
int if_else_tester (Vector zerovec)
{
  using std::sin;

  typedef typename Vector::value_type Scalar;

  Vector random_vec = zerovec;

  Vector error_vec = zerovec;


  std::srand(12345); // Fixed seed for reproduceability of failures

  for (unsigned int i=0; i != random_vec.size(); ++i)
    random_vec[i] = .25 + (static_cast<Scalar>(std::rand())/RAND_MAX);

  int returnval = 0;

  one_test((random_vec > 0.75) * 1.0 - (2*random_vec > 1.5) * 1.0);

  one_test(if_else(random_vec > 0.75,  random_vec, sin(random_vec)) -
           if_else(2*random_vec > 1.5, random_vec, sin(random_vec)));

  return returnval;
}



int main(void)
{
  int returnval = 0;
  returnval = returnval || vectester(NumberArray<N, float>());
  returnval = returnval || vectester(NumberArray<N, double>());
  returnval = returnval || vectester(NumberArray<N, long double>());

  // We no longer treat vectors like arrays for built-in functions, so
  // most of the identities above make no sense.
  /*
  returnval = returnval || vectester(NumberVector<N, float>());
  returnval = returnval || vectester(NumberVector<N, double>());
  returnval = returnval || vectester(NumberVector<N, long double>());
  */

  returnval = returnval ||
              vectester(SparseNumberArrayOf
                          <4, 0, float, 1, float,
                              2, float, 3, float>::type());
  returnval = returnval ||
              vectester(SparseNumberArrayOf
                          <4, 0, double, 1, double,
                              2, double, 3, double>::type());
  returnval = returnval ||
              vectester(SparseNumberArrayOf
                          <4, 0, long double, 1, long double,
                              2, long double, 3, long double>::type());

  returnval = returnval ||
              vectester(SparseNumberVectorOf
                          <4, 0, float, 1, float,
                              2, float, 3, float>::type());
  returnval = returnval ||
              vectester(SparseNumberVectorOf
                          <4, 0, double, 1, double,
                              2, double, 3, double>::type());
  returnval = returnval ||
              vectester(SparseNumberVectorOf
                          <4, 0, long double, 1, long double,
                              2, long double, 3, long double>::type());

  DynamicSparseNumberArray<float, unsigned int> float_dsna;
    float_dsna.resize(4);
    float_dsna.raw_index(1) = 1;
    float_dsna.raw_index(2) = 2;
    float_dsna.raw_index(3) = 3;
  returnval = returnval || vectester(float_dsna);
  returnval = returnval || if_else_tester(float_dsna);

  DynamicSparseNumberArray<double, unsigned int> double_dsna;
    double_dsna.resize(4);
    double_dsna.raw_index(1) = 1;
    double_dsna.raw_index(2) = 2;
    double_dsna.raw_index(3) = 3;
  returnval = returnval || vectester(double_dsna);
  returnval = returnval || if_else_tester(double_dsna);

  DynamicSparseNumberArray<long double, unsigned int> long_double_dsna;
    long_double_dsna.resize(4);
    long_double_dsna.raw_index(1) = 1;
    long_double_dsna.raw_index(2) = 2;
    long_double_dsna.raw_index(3) = 3;
  returnval = returnval || vectester(long_double_dsna);
  returnval = returnval || if_else_tester(long_double_dsna);

  DynamicSparseNumberVector<float, unsigned int> float_dsnv;
    float_dsnv.resize(4);
    float_dsnv.raw_index(1) = 1;
    float_dsnv.raw_index(2) = 2;
    float_dsnv.raw_index(3) = 3;
  returnval = returnval || vectester(float_dsnv);
  returnval = returnval || if_else_tester(float_dsnv);

  DynamicSparseNumberVector<double, unsigned int> double_dsnv;
    double_dsnv.resize(4);
    double_dsnv.raw_index(1) = 1;
    double_dsnv.raw_index(2) = 2;
    double_dsnv.raw_index(3) = 3;
  returnval = returnval || vectester(double_dsnv);
  returnval = returnval || if_else_tester(double_dsnv);

  DynamicSparseNumberVector<long double, unsigned int> long_double_dsnv;
    long_double_dsnv.resize(4);
    long_double_dsnv.raw_index(1) = 1;
    long_double_dsnv.raw_index(2) = 2;
    long_double_dsnv.raw_index(3) = 3;
  returnval = returnval || vectester(long_double_dsnv);
  returnval = returnval || if_else_tester(long_double_dsnv);

  return returnval;
}
