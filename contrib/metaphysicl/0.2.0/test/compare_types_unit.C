
#include "metaphysicl/compare_types.h"

using namespace MetaPhysicL;

template <typename T1, typename T2>
struct SubInstantiator {
  typename CompareTypes<T1 ,       T2 >::supertype st1;
  typename CompareTypes<T1 , const T2 >::supertype st3;
  typename CompareTypes<T1 , const T2&>::supertype st4;
};

template <typename T1, typename T2>
struct Instantiator {
  SubInstantiator<      T1 , T2> si1;
  SubInstantiator<const T1 , T2> si2;
  SubInstantiator<const T1&, T2> si4;

  void foil_unused_variable_warning() {}
};


int main (void)
{
  Instantiator<float, float> i1;
  Instantiator<float, double> i2;
  Instantiator<double, double> i3;
  Instantiator<int, int> i4;

  i1.foil_unused_variable_warning();
  i2.foil_unused_variable_warning();
  i3.foil_unused_variable_warning();
  i4.foil_unused_variable_warning();

  return 0;
}
