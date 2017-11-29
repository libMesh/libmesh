
// MetaPhysicL
#include "metaphysicl/metaphysicl_asserts.h"
#include "metaphysicl/sparsenumbervector.h"
#include "metaphysicl/namedindexarray.h"

#include "metaphysicl_config.h"

// VexCL
#ifdef METAPHYSICL_HAVE_VEXCL
#include "vexcl/vexcl.hpp"
#endif

// C++
#include <iostream>

using namespace MetaPhysicL;

int main(void)
{
  typedef
    NamedIndexArray
      <double,
       SparseNumberVector
         <long unsigned int,
          ULongSetConstructor<3>::type> >
    indexed_by_three;

  indexed_by_three test_val;
  test_val.raw_data() = 2;
  test_val.raw_sizes().get<3>() = 1;

  auto test_val_2 = test_val * test_val;
  metaphysicl_assert_equal_to(test_val_2.raw_data(), 4);

#ifdef METAPHYSICL_HAVE_VEXCL
  vex::Context ctx (vex::Filter::Env && vex::Filter::Count(1));
  std::cout << ctx << std::endl;
      

  typedef
    NamedIndexArray
      <vex::vector<double>,
       SparseNumberVector
         <long unsigned int,
          ULongSetConstructor<2>::type> >
    vex_indexed_by_two;

  typedef
    NamedIndexArray
      <vex::vector<double>,
       SparseNumberVector
         <long unsigned int,
          ULongSetConstructor<1>::type> >
    vex_indexed_by_one;

  vex_indexed_by_one test_one(vex::vector<double>(ctx, 5),0);
  test_one.raw_sizes().template get<1>() = 5;
  vex_indexed_by_two test_two(vex::vector<double>(ctx, 3),0);
  test_two.raw_sizes().template get<2>() = 3;

  test_one.raw_data()[2] = 7;
  test_two.raw_data()[1] = 2;

  auto test_three = test_one * test_two;

  if (test_three.raw_sizes().template get<1>() != 5)
    return 1;

  if (test_three.raw_sizes().template get<2>() != 3)
    return 1;

  vex::vector<double> test_output = test_three.raw_data();

  if (test_output[7] != 14)
    return 1;

#endif

  return 0;
}
