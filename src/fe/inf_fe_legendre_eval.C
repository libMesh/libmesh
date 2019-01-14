// The libMesh Finite Element Library.
// Copyright (C) 2002-2019 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


// Local Includes
#include "libmesh/libmesh_config.h"
#ifdef LIBMESH_ENABLE_INFINITE_ELEMENTS
#include "libmesh/inf_fe.h"
#include "libmesh/inf_fe_macro.h"

namespace libMesh
{

// Anonymous namespace for local helper functions
namespace {

Real legendre_eval(Real v,  unsigned i)
{
  libmesh_assert (-1.-1.e-5 <= v && v < 1.);

  switch (i)
    {
    case 0:
      return 1.;

    case 1:
      return v+1.;

    case 2:
      return 1.5*v*v-1.5;

    case 3:
      return 1.+(-1.5+2.5*v*v)*v;

    case 4:
      return -.625+(-3.75+4.375*v*v)*v*v;

    case 5:
      return 1.+(1.875+(-8.75+7.875*v*v)*v*v)*v;

    case 6:
      return -1.3125+
        (6.5625+
         (-19.6875+14.4375*v*v)*v*v)*v*v;

    case 7:
      return 1.+
        (-2.1875+
         (19.6875+
          (-43.3125+26.8125*v*v)*v*v)*v*v)*v;

    case 8:
      return -.7265625+
        (-9.84375+
         (54.140625+
          (-93.84375+50.2734375*v*v)*v*v)*v*v)*v*v;

    case 9:
      return 1.+
        (2.4609375+
         (-36.09375+
          (140.765625+
           (-201.09375+94.9609375*v*v)*v*v)*v*v)*v*v)*v;

    case 10:
      return -1.24609375+
        (13.53515625+
         (-117.3046875+
          (351.9140625+
           (-427.32421875+180.42578125*v*v)*v*v)*v*v)*v*v)*v*v;

    case 11:
      return 1.+
        (-2.70703125+
         (58.65234375+
          (-351.9140625+
           (854.6484375+
            (-902.12890625+344.44921875*v*v)*v*v)*v*v)*v*v)*v*v)*v;

    case 12:
      return -.7744140625+
        (-17.595703125+
         (219.9462890625+
          (-997.08984375+
           (2029.7900390625+
            (-1894.470703125+660.1943359375*v*v)*v*v)*v*v)*v*v)*v*v)*v*v;

    case 13:
      return 1.+
        (2.9326171875+
         (-87.978515625+
          (747.8173828125+
           (-2706.38671875+
            (4736.1767578125+
             (-3961.166015625+1269.6044921875*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v;

    case 14:
      return -1.20947265625+
        (21.99462890625+
         (-373.90869140625+
          (2368.08837890625+
           (-7104.26513671875+
            (10893.20654296875+
             (-8252.42919921875+2448.52294921875*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v*v;

    case 15:
      return 1.+
        (-3.14208984375+
         (124.63623046875+
          (-1420.85302734375+
           (7104.26513671875+
            (-18155.34423828125+
             (24757.28759765625+
              (-17139.66064453125+4733.81103515625*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v;

    case 16:
      return -.803619384765625+
        (-26.707763671875+
         (592.0220947265625+
          (-4972.985595703125+
           (20424.76226806641+
            (-45388.36059570313+
             (55703.89709472656+
              (-35503.58276367188+9171.758880615234*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v*v;

    case 17:
      return 1.+
        (3.338470458984375+
         (-169.149169921875+
          (2486.492797851563+
           (-16339.80981445313+
            (56735.45074462891+
             (-111407.7941894531+
              (124262.5396728516+
               (-73374.07104492188+17804.00253295898*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v;

    case 18:
      return -1.185470581054688+
        (31.71546936035156+
         (-888.0331420898438+
          (9531.555725097656+
           (-51061.90567016602+
            (153185.717010498+
             (-269235.5026245117+
              (275152.766418457+
               (-151334.0215301514+34618.89381408691*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v*v;


    default:
      libmesh_error_msg("bad index i = " << i);
    }
} // legendre_eval()




Real legendre_eval_deriv(Real v, unsigned i)
{
  libmesh_assert (-1.-1.e-5 <= v && v < 1.);

  switch (i)
    {
    case 0:
      return 0.;

    case 1:
      return 1.;

    case 2:
      return 3.*v;

    case 3:
      return 7.5*v*v-1.5;

    case 4:
      return (-7.5+17.5*v*v)*v;

    case 5:
      return 1.875+(-26.25+39.375*v*v)*v*v;

    case 6:
      return
        (13.125+
         (-78.75+86.625*v*v)*v*v)*v;

    case 7:
      return -2.1875+
        (59.0625+
         (-216.5625+187.6875*v*v)*v*v)*v*v;

    case 8:
      return
        (-19.6875+
         (216.5625+
          (-563.0625+402.1875*v*v)*v*v)*v*v)*v;

    case 9:
      return 2.4609375+
        (-108.28125+
         (703.828125+
          (-1407.65625+854.6484375*v*v)*v*v)*v*v)*v*v;

    case 10:
      return
        (27.0703125+
         (-469.21875+
          (2111.484375+
           (-3418.59375+1804.2578125*v*v)*v*v)*v*v)*v*v)*v;

    case 11:
      return -2.70703125+
        (175.95703125+
         (-1759.5703125+
          (5982.5390625+
           (-8119.16015625+3788.94140625*v*v)*v*v)*v*v)*v*v)*v*v;

    case 12:
      return
        (-35.19140625+
         (879.78515625+
          (-5982.5390625+
           (16238.3203125+
            (-18944.70703125+7922.33203125*v*v)*v*v)*v*v)*v*v)*v*v)*v;

    case 13:
      return 2.9326171875+
        (-263.935546875+
         (3739.0869140625+
          (-18944.70703125+
           (42625.5908203125+
            (-43572.826171875+16504.8583984375*v*v)*v*v)*v*v)*v*v)*v*v)*v*v;

    case 14:
      return
        (43.9892578125+
         (-1495.634765625+
          (14208.5302734375+
           (-56834.12109375+
            (108932.0654296875+
             (-99029.150390625+34279.3212890625*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v;

    case 15:
      return -3.14208984375+
        (373.90869140625+
         (-7104.26513671875+
          (49729.85595703125+
           (-163398.0981445313+
            (272330.1635742188+
             (-222815.5883789063+71007.16552734375*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v*v;

    case 16:
      return
        (-53.41552734375+
         (2368.08837890625+
          (-29837.91357421875+
           (163398.0981445313+
            (-453883.6059570313+
             (668446.7651367188+
              (-497050.1586914063+146748.1420898438*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v;

    case 17:
      return 3.338470458984375+
        (-507.447509765625+
         (12432.46398925781+
          (-114378.6687011719+
           (510619.0567016602+
            (-1225485.736083984+
             (1615413.01574707+
              (-1100611.065673828+302668.0430603027*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v*v;

    case 18:
      return
        (63.43093872070313+
         (-3552.132568359375+
          (57189.33435058594+
           (-408495.2453613281+
            (1531857.17010498+
             (-3230826.031494141+
              (3852138.729858398+
               (-2421344.344482422+623140.0886535645*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v*v)*v;


    default:
      libmesh_error_msg("bad index i = " << i);
    }
} // legendre_eval_deriv()

} // anonymous namespace





  // Specialize the eval() function for 1, 2, and 3 dimensions and the CARTESIAN mapping type
  // to call the local helper function from the anonymous namespace.
template <> Real InfFE<1,LEGENDRE,CARTESIAN>::eval(Real v, Order, unsigned i) { return legendre_eval(v, i); }
template <> Real InfFE<2,LEGENDRE,CARTESIAN>::eval(Real v, Order, unsigned i) { return legendre_eval(v, i); }
template <> Real InfFE<3,LEGENDRE,CARTESIAN>::eval(Real v, Order, unsigned i) { return legendre_eval(v, i); }

// Specialize the eval_deriv() function for 1, 2, and 3 dimensions and the CARTESIAN mapping type
// to call the local helper function from the anonymous namespace.
template <> Real InfFE<1,LEGENDRE,CARTESIAN>::eval_deriv(Real v, Order, unsigned i) { return legendre_eval_deriv(v, i); }
template <> Real InfFE<2,LEGENDRE,CARTESIAN>::eval_deriv(Real v, Order, unsigned i) { return legendre_eval_deriv(v, i); }
template <> Real InfFE<3,LEGENDRE,CARTESIAN>::eval_deriv(Real v, Order, unsigned i) { return legendre_eval_deriv(v, i); }

} // namespace libMesh


#endif // LIBMESH_ENABLE_INFINITE_ELEMENTS
