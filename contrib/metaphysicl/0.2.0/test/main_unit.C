
#include <iostream>
#include "metaphysicl/dualshadowvector.h"
#include "metaphysicl/dualshadowsparsevector.h"
#include "metaphysicl/dualshadowsparsestruct.h"
#include "metaphysicl/metaphysicl_asserts.h"

using namespace MetaPhysicL;

struct print_functor {
  template <typename T>
  void operator()(const T& i) const { std::cout << i << std::endl; }
};

int main(void)
{
  const RawType<const ShadowNumber<double, long double> >::value_type testin = 0;
  double testout = testin;

  if (testout)
    metaphysicl_error();

  DualNumber<double, double> dndd = 1.;
  2. * dndd;

  typedef SetConstructor<IntType<1>, IntType<7>, IntType<9>, IntType<5>, IntType<7> >::type SType;

  typedef SetConstructor<IntType<9>, IntType<7> >::type SType2;

  typedef SetConstructor<IntType<10,SType>, IntType<11,SType2> >::type SSType;

  typedef SetConstructor<IntType<1> >::type DumbType1;
  typedef SetConstructor<IntType<1,DumbType1> >::type DumbType2;

  ctassert<TypesEqual<SetOfSetsUnion<SSType>::type, SType>::value>::apply();

  typedef SetConstructor<UnsignedIntType<2,double>, UnsignedIntType<9,double> >::type SType3;

  typedef SetConstructor<IntType<3,double> >::type SType4;

  typedef SetConstructor<UnsignedIntType<0>, UnsignedIntType<1> >::type IS2D;

  std::cout << "sizeof(NullContainer) = " << sizeof(NullContainer<IntType<0> >) << std::endl;
  std::cout << "sizeof(IntType) = " << sizeof(IntType<0>) << std::endl;
  std::cout << "sizeof(IntDbl)  = " << sizeof(IntType<0,double>) << std::endl;
  std::cout << "sizeof(SType)   = " << sizeof(SType)  << std::endl;
  std::cout << "sizeof(STypeDbl)= " << sizeof(SType::rebind<double>::other)  << std::endl;
  std::cout << "sizeof(SType2)  = " << sizeof(SType2)  << std::endl;
  std::cout << "sizeof(SType3)  = " << sizeof(SType3) << std::endl;
  std::cout << "sizeof(SType4)  = " << sizeof(SType4) << std::endl;
  std::cout << "sizeof(DumbType2)     = " << sizeof(DumbType2) << std::endl;
  std::cout << "sizeof(IS2D)    = " << sizeof(IS2D) << std::endl;

  SType4 s4;

  s4.element<IntType<3> >();

  SType2 s;

  s.element<IntType<9> >();

  SparseNumberStruct<SType3> a3 = 2, b3 = 3, c3 = 0;

  a3.get<9>() = 1;
  b3.get<9>() = 5;
  a3.get<2>() = 7;

  std::cout << "a3Xb3    = " << a3.outerproduct(b3) << std::endl;
  std::cout << "(a3Xb3)' = " << transpose(a3.outerproduct(b3)) << std::endl;
  std::cout << "(a3Xb3)'[2] = " << (transpose(a3.outerproduct(b3))).get<2>() << std::endl;
  std::cout << "(a3Xb3)'[2][2] = " << (transpose(a3.outerproduct(b3))).get<2>().get<2>() << std::endl;

  SType::RuntimeForEach()(print_functor());

  SparseNumberVector<double, SType2> a2(0), b2(0), c2(0);

  SparseNumberVector<double, SType> a(0), b(0), c(0);

  SparseNumberVector<SparseNumberVector<double, SType>, SType> T(0);

  std::cout << "Found 8 = " << SType::Contains<IntType<8> >::value << std::endl;

  a2.raw_at(1) = 2.;

  a = a2;

  T.raw_at(2).raw_at(2) = 5.;

  b[5] = 3.;

  c = 2*a;

  c = a / 3 + b;

  if (T)
    std::cout << "T is non-zero" << std::endl;

  std::cout << "c = " << c << std::endl;
  std::cout << "a = " << a << std::endl;
  std::cout << "T = " << T << std::endl;
  std::cout << "c.dot(a) = " << c.dot(a) << std::endl;
  std::cout << "T.dot(a) = " << T.dot(a) << std::endl;

  const unsigned int NDIM = 2;
  typedef double RawScalar;

  typedef SparseNumberStructUnitVector<NDIM,0,RawScalar>::type XVector;
  XVector xvec = SparseNumberStructUnitVector<NDIM,0,RawScalar>::value();
  typedef SparseNumberStructFullVector<NDIM,RawScalar>::type FullVector;

/*
  typedef SparseNumberVectorUnitVector<NDIM,0,RawScalar>::type XVector;
  XVector xvec = SparseNumberVectorUnitVector<NDIM,0,RawScalar>::value();
  typedef SparseNumberVectorFullVector<NDIM,RawScalar>::type FullVector;
*/

/*
  typedef NumberVectorUnitVector<NDIM,0,RawScalar>::type XVector;
  XVector xvec = NumberVectorUnitVector<NDIM,0,RawScalar>::value();
  typedef NumberVector<NDIM,RawScalar> FullVector;
*/

  typedef DualNumber<RawScalar, XVector>    XFirstDerivType;
  typedef DualNumber<RawScalar, FullVector> FullFirstDerivType;

  typedef DualNumber<XFirstDerivType, XVector::rebind<XFirstDerivType>::other> XSecondDerivType;
  typedef DualNumber<FullFirstDerivType, FullVector::rebind<FullFirstDerivType>::other> FullSecondDerivType;
  typedef DualNumberConstructor<FullFirstDerivType, FullVector::rebind<FullFirstDerivType>::other> FSDTConstructor;

  std::cout << "sizeof(FSDTC)   = " << sizeof(FSDTConstructor) << std::endl;

  FullFirstDerivType  ffdt = XFirstDerivType(1.,xvec);
  FullSecondDerivType fsdt = 0;
  XSecondDerivType xsdt = XSecondDerivType(1.,xvec);
  FullVector::rebind<FullFirstDerivType>::other first_deriv_vector;

  fsdt = xsdt;

  fsdt.derivatives() /= fsdt.value();
  fsdt.value()/(fsdt.value()*fsdt.value());
  ffdt = fsdt.value();
  first_deriv_vector = fsdt.derivatives();

//  ctprint<FullFirstDerivType,FullVector::rebind<FullFirstDerivType>::other>::apply();
  ffdt * first_deriv_vector;
//  ffdt * fsdt.derivatives();
//  fsdt.value()/(fsdt.value()*fsdt.value()) * fsdt.derivatives();
//  fsdt.derivatives() -= fsdt.value()/(fsdt.value()*fsdt.value()) * fsdt.derivatives();

  fsdt.value() /= fsdt.value();

  fsdt / fsdt;

  return 0;
}
