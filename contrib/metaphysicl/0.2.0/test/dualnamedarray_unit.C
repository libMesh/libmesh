
// MetaPhysicL
#include "metaphysicl_config.h"

// VexCL
#ifdef METAPHYSICL_HAVE_VEXCL
#include "vexcl/vexcl.hpp"
#endif

// MetaPhysicL
#include "metaphysicl/compare_types.h"

#ifdef METAPHYSICL_HAVE_VEXCL
namespace MetaPhysicL {

// vex::vector expressions don't nest anything except native types in
// this code, so we can treat them as native types.

  template <typename T>
  struct BuiltinTraits
    <T,typename std::enable_if<vex::is_vector_expression<T>::value>::type>
  {
    static const bool value = true;
  };

// vex::vector is the only expression in this code referring to
// underlying storage

  template <typename T>
  struct copy_or_reference <vex::vector<T>&>
  {
    typedef vex::vector<T>& type;

    static const bool copy = false;
  };

  template <typename T>
  struct copy_or_reference <const vex::vector<T>&>
  {
    typedef const vex::vector<T>& type;

    static const bool copy = false;
  };

} // namespace MetaPhysicL
#endif

#include "metaphysicl/sparsenumbervector.h"
//#include "metaphysicl/dualnamedarray.h"
#include "metaphysicl/namedindexarray.h"

#include "metaphysicl/metaphysicl_asserts.h"

// C++
#include <iostream>

using namespace MetaPhysicL;

#define test_assert_equal_to(expr1,expr2)  do { if (!(expr1 == expr2)) { std::cerr << "Assertion `" #expr1 " == " #expr2 "' failed.\n" #expr1 " = " << (expr1) << "\n" #expr2 " = " << (expr2) << std::endl; metaphysicl_error(); } } while(0)



int main(void)
{
  typedef
    NamedIndexArray
      <double,
       SparseNumberVector
         <long unsigned int,
          ULongSetConstructor<3>::type> >
    indexed_by_three;

  /*
  typedef
    DualExpression<indexed_by_three, indexed_by_three> dual_three;

  dual_three test_val;
  test_val.value().raw_data() = 0.5;
  test_val.derivatives().raw_data() = 1;
  test_val.value().raw_sizes().get<3>() = 1;
  test_val.derivatives().raw_sizes().get<3>() = 1;
  */

  indexed_by_three val;
  val.raw_data() = 0.5;
  val.raw_sizes().get<3>() = 1;

  /*
  auto test_eight = make_dual_expression_copy(val,val)*val;

  double test_eight_output =
    test_eight.derivatives().raw_data();

  metaphysicl_assert_equal_to(test_eight_output, 0.25);
  */
#ifdef METAPHYSICL_HAVE_VEXCL
  // Passes, as it should
  ctassert<BuiltinTraits<vex::vector<double> >::value>::apply();

  // Fails, as it should
  // ctassert<BuiltinTraits<DualExpression<double,double> >::value>::apply();

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
    DualExpression<vex_indexed_by_two, vex_indexed_by_two> dual_two;

  typedef
    NamedIndexArray
      <vex::vector<double>,
       SparseNumberVector
         <long unsigned int,
          ULongSetConstructor<1>::type> >
    vex_indexed_by_one;

  typedef
    DualExpression<vex_indexed_by_one, vex_indexed_by_one> dual_one;

  double zeros[5];
  for (unsigned int i = 0; i != 5; ++i)
    zeros[i] = 0;

  vex::vector<double> raw_vex{ctx, 5, zeros};

  vex_indexed_by_one test_one_val(vex::vector<double>(ctx, 5, zeros),0);
  test_one_val.raw_sizes().template get<1>() = 5;
  test_one_val.raw_data()[2] = 7;

  vex_indexed_by_one test_one_deriv(vex::vector<double>(ctx, 5, zeros),0);
  test_one_deriv.raw_sizes().template get<1>() = 5;

  dual_one test_one(test_one_val, test_one_deriv);

  vex_indexed_by_two test_two_val(vex::vector<double>(ctx, 3, zeros),0);
  test_two_val.raw_sizes().template get<2>() = 3;
  test_two_val.raw_data()[1] = 2;

  vex_indexed_by_two test_two_deriv(vex::vector<double>(ctx, 3, zeros),0);
  test_two_deriv.raw_sizes().template get<2>() = 3;

  dual_two test_two(test_two_val, test_two_deriv);

  test_assert_equal_to (test_one.value().raw_data()[2], 7);
  for (unsigned int i = 0; i != 5; ++i)
    if (i != 2)
      test_assert_equal_to (test_one.value().raw_data()[i], 0);

  test_assert_equal_to (test_two.value().raw_data()[1], 2);
  for (unsigned int i = 0; i != 3; ++i)
    if (i != 1)
      test_assert_equal_to (test_two.value().raw_data()[i], 0);

  auto test_make = make_dual_expression_reference(test_one_val,
                                                  test_one_val);
  auto test_multiply_raw = test_make * raw_vex;

  auto test_lmultiply_raw = raw_vex * test_make;

  auto test_make_deriv = (test_make.derivatives())/2.0;
  vex::vector<double> test_make_deriv_output =
    test_make_deriv.raw_data();

  auto test_make_direct = (test_make-1.0)/2.0;

  vex::vector<double> test_make_output =
    test_make_direct.derivatives().raw_data();


  auto test_three_val_val = test_one_val * test_two_val;
  auto test_three_deriv_val = test_one_deriv * test_two_val;

  std::cout << "test_one_val = " << test_one_val << std::endl;
  test_assert_equal_to (test_one_val.raw_data()[2], 7);
  for (unsigned int i = 0; i != 5; ++i)
    if (i != 2)
      test_assert_equal_to (test_one_val.raw_data()[i], 0);

  auto test_three_ag = (test_one_val * test_one_val);
  vex::vector<double> test_output_three_ag = test_three_ag.raw_data();
  std::cout << "test_output_three_ag = " << test_output_three_ag << std::endl;
  test_assert_equal_to (test_output_three_ag[2], 49);
  for (unsigned int i = 0; i != 5; ++i)
    if (i != 2)
      test_assert_equal_to (test_output_three_ag[i], 0);

  auto test_three_am = (test_one * test_one_val);
  vex::vector<double> test_output_three_am = test_three_am.value().raw_data();
  std::cout << "test_output_three_am = " << test_output_three_am << std::endl;
  test_assert_equal_to (test_output_three_am[2], 49);
  for (unsigned int i = 0; i != 5; ++i)
    if (i != 2)
      test_assert_equal_to (test_output_three_am[i], 0);

  auto test_three_a = test_one * test_two_val;
  test_assert_equal_to (test_three_a.value().raw_sizes().template get<1>(), 5);
  test_assert_equal_to (test_three_a.value().raw_sizes().template get<2>(), 3);
  vex::vector<double> test_output_three_a = test_three_a.value().raw_data();
  std::cout << "test_output_three_a = " << test_output_three_a << std::endl;
  test_assert_equal_to (test_output_three_a[7], 14);
  for (unsigned int i = 0; i != 15; ++i)
    if (i != 7)
      test_assert_equal_to (test_output_three_a[i], 0);

  auto test_three = test_one * test_two;
  test_assert_equal_to (test_three.value().raw_sizes().template get<1>(), 5);
  test_assert_equal_to (test_three.value().raw_sizes().template get<2>(), 3);
  vex::vector<double> test_output_three = test_three.value().raw_data();
  test_assert_equal_to (test_output_three[7], 14);
  for (unsigned int i = 0; i != 15; ++i)
    if (i != 7)
      test_assert_equal_to (test_output_three[i], 0);

  auto test_four = test_one / 2.0;
  vex::vector<double> test_output_four = test_four.value().raw_data();
  test_assert_equal_to (test_output_four[2], 3.5);
  for (unsigned int i = 0; i != 5; ++i)
    if (i != 2)
      test_assert_equal_to (test_output_four[i], 0);

  auto test_five = test_one * test_one;
  vex::vector<double> test_output_five = test_five.value().raw_data();
  test_assert_equal_to (test_output_five[2], 49);
  for (unsigned int i = 0; i != 5; ++i)
    if (i != 2)
      test_assert_equal_to (test_output_five[i], 0);

  auto test_six = test_five / 2.0;
  vex::vector<double> test_output_six = test_six.value().raw_data();
  test_assert_equal_to (test_output_six[2], 24.5);
  for (unsigned int i = 0; i != 5; ++i)
    if (i != 2)
      test_assert_equal_to (test_output_six[i], 0);

  auto test_seven = (test_one * test_one_val) / 2.0;
  vex::vector<double> test_output_seven = test_seven.value().raw_data();
  test_assert_equal_to (test_output_seven[2], 24.5);
  for (unsigned int i = 0; i != 5; ++i)
    if (i != 2)
      test_assert_equal_to (test_output_seven[i], 0);

  vex_indexed_by_one test_nine_val(vex::vector<double>(ctx, 1, zeros), 0);
  test_nine_val.raw_sizes().template get<1>() = 1;

  vex::vector<double> test_output_nine =
    (make_dual_expression_reference(test_nine_val, test_nine_val) +
     test_nine_val).derivatives().raw_data();

  vex_indexed_by_one test_ten_val(vex::vector<double>(ctx, 1, zeros), 0);

  auto test_ten =
    make_dual_expression_reference(test_ten_val, test_ten_val) + test_ten_val;

  vex::vector<double> test_output_ten =
    test_ten.derivatives().raw_data();

#endif

  return 0;
}
